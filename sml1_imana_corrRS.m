function varargout = sml1_imana_corrRS(what,varargin)
% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';          
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia'];
pcmDir          =[baseDir '/pcm_stats'];
repSupDir       =[baseDir '/repsup_stats'];
corrTestDir     =[baseDir,'/corrTest'];

% update glmDir when adding new glms
glmFoSExDir     ={[baseDir '/glmFoSEx/glmFoSEx1'],[baseDir '/glmFoSEx/glmFoSEx2'],[baseDir '/glmFoSEx/glmFoSEx3'],[baseDir '/glmFoSEx/glmFoSEx4']};    

% ------------------------- Experiment Info -------------------------------

% Stimuli - numbers given in SeqNumb
num_train = 1:6;
num_untrain = 7:12;   
num_seq = 1:12;
num_fing = 13:17;
num_seqtype = 2;
SeqType={'trained','untrained'};

% for group 1 - group 2 (1-6 and 7-12 reversed)

% per session
numruns_sess      = 10;  
numruns_task_sess = 8;
numruns_loc_sess  = 2;
          
% seqNumb - all sequences: 1-19
% seqType - types of sequences 
    % 1 - training 
    % 2 - untrained (other group)
    % 3 - finger mapping

% ------------------------- ROI things ------------------------------------
hem        = {'lh','rh'};                                                   % left & right hemi folder names/prefixes
hemName    = {'LeftHem','RightHem'};
regname         = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp','CaudateN' 'Pallidum', 'Putamen' 'Thalamus','CIV','CV','CVI'};
regname_cortex  = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp'};
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen', 'Thalamus'};
regname_cerebellum = {'LobIV','LobV','LobVI'};
numregions_surf = 8;
numregions_BG   = 4;
numregions_cerebellum = 3;
numregions = numregions_surf+numregions_BG+numregions_cerebellum;        
regSide=[ones(1,8) ones(1,8)*2]; % 1-left, 2-right
regType=[1:8  1:8]; % cortical areas: 1-8, BG: 8-12, cereb: 13-15


% ------------------------- Freesurfer things -----------------------------         
atlasA    = 'x';                                                            % freesurfer filename prefix
atlasname = 'fsaverage_sym';                                                % freesurfer average atlas
hemName   = {'LeftHem','RightHem'};                                         % freesurfer hemisphere folder names    

% ------------------------- Subject things --------------------------------

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22'};  
num_run = 8; % 8 functional runs
% Other random notes %

% -------------------------- For plotting ---------------------------------
gray=[120 120 120]/255;
lightgray=[170 170 170]/255;
silver=[230 230 230]/255;
black=[0 0 0]/255;
blue=[49,130,189]/255;
lightblue=[158,202,225]/255;
red=[222,45,38]/255;
lightred=[252,146,114]/255;
ms=12;
stySeq=style.custom({'red','blue'},'markersize',ms);
stySeqType_exe=style.custom({gray,lightgray},'markersize',ms);
styTrained_exe=style.custom({red,lightred},'markersize',ms);
styUntrained_exe=style.custom({blue,lightblue},'markersize',ms);
stySess=style.custom({gray,lightgray,silver,black},'markersize',ms);


switch(what)
    case 'PCM_data_repsupModels' 
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1,3]; % S1 - positive pcm, PMd - negative pcm corr
        sn=[1:9,11:22];
        sessN=1; % need to be two sessions at the time
        seqType='trained';
        models={'generic','specific','spec+genCorr'};
        M_type=1; % 1 - generic, 2 - specific
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','seqType','M_type'})
        
        switch seqType
            case 'trained'
                stIndx=1;
            case 'untrained'
                stIndx=2;
        end
        
        for ss = sessN
            AllReg=[];
            for r = reg
                for p=1:length(sn)
                    B=load(fullfile(regDir,sprintf('betas_FoSEx_sess%d.mat',ss)));
                    glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                    
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==r)}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                    end
                    
                    % conditions - 1-6 for first exe, 7-12 for second
                    cond = D.seqNumb;
                    if stIndx==1
                        % make condVec for second exe into 7-12
                        cond(D.FoSEx==2)=cond(D.FoSEx==2)+6;
                    else
                        % make condVec for first exe into 1-6
                        cond(D.FoSEx==1)=cond(D.FoSEx==1)-6;
                    end
                    condVec{p} = cond(D.seqType==stIndx); % conditions
                    partVec{p} = D.run(D.seqType==stIndx);
                    Data{p} = beta(:,D.seqType==stIndx)';  % Data is N x P (cond x voxels) - no intercept
                end; % subj
                
                % construct models
                switch M_type
                    case 1 % generic
                        %M = pcm_corrModel;
                        M = pcm_repsupModel_generic;
                    case 2 % specific
                        % M = pcm_corrModel_indSeq;
                        M = pcm_repsupModel_specific;
                    case 3
                        M=pcm_repsupModel_specific_genCorr;
                end
                T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                if M_type==3
                    C = pcm_correlation(Data,partVec,condVec,M{5},runEffect,M_type);
                else
                    C = pcm_correlation(Data,partVec,condVec,M{4},runEffect,M_type);
                end
                T.roi = ones(size(T.SN))*r;
                AllReg=addstruct(AllReg,C);
                AllReg=addstruct(AllReg,T);
            end % region
            %remove some fields
         %   a1='reg'; a2='theta'; a3='theta_hat'; a4='thetaCr';
          %  AllReg=rmfield(AllReg,a1); AllReg=rmfield(AllReg,a2); AllReg=rmfield(AllReg,a3); AllReg=rmfield(AllReg,a4);
            % save output
            save(fullfile(corrTestDir,sprintf('PCM_RS_corr_%s_%s_sess%d.mat',models{M_type},seqType,ss)),'-struct','AllReg');
        end % session
    case 'PCM_simulate'
        runEffect = 'fixed';
        algorithm='NR'; % minimize or NR
        models={'generic','specific'};
        M_type=1; % 1 - generic, 2 - specific
        corrType='positive';
        noise=1;
        scale=[0:0.1:1];
        vararginoptions(varargin,{'runEffect','theta','noise','scale','M_type','corrType'});
        %scale=[0:0.1:1];
        switch(corrType)
            case 'negative'
                if M_type==1
                    theta=[-0.0365 -0.0003 -0.0144 -0.0035 0.0060]';
                elseif M_type==2
                    theta=[-0.0321 0.0230 0.0003  -0.0474 -0.0156 -0.0245 -0.0282 -0.0267 -0.0001,...
                        -0.0254 0.0005 0.0002 -0.0414 0.0005 -0.0001 0.0085 0.0178 0.0170 0.0057 0.0094]';
                end
            case 'positive'
                if M_type==1
                    theta=[-0.0213 -0.0001 -0.0205 -0.0002 -0.0135]';
                elseif M_type==2
                    theta=[-0.0001 0.0209 -0.0166 -0.0249 -0.02 0.0037 -0.0351 -0.0182 0.0001 0.0052 -0.0001,...
                        -0.0001 0.0001 -0.0112 -0.0128 -0.0178 -0.0175 -0.0129 -0.0158 -0.0017]';
                end
        end
        part=8;
        cond=6;

        % model specifications
        D.numVox  = 1000;
        D.condVec = [kron(ones(part,1),[1:cond+cond]')];      
        %  D.partVec = [kron([1:part]',ones(cond,1)); kron([1:part]',ones(cond,1))];  % Partitions
        D.partVec = [kron([1:part]',ones(cond+cond,1))];
        
        sn = 25; % 25 subjects
        
        switch M_type
            case 1 % generic
                %M = pcm_corrModel;
                M = pcm_repsupModel_generic;
            case 2 % specific
                % M = pcm_corrModel_indSeq;
                M = pcm_repsupModel_specific;
        end
        
        SS=[]; TT=[];
        if M_type == 2
            % need thetas for each sequence - 6x per parameter
            var1       = (theta(3:8,:).^2)';
            var2       = (theta(9:14,:).^2+theta(15:20,:).^2)';
            cov12      = (theta(3:8,:).*theta(15:20,:))';
            trueCorr = mean(cov12)/sqrt(mean(var1)*mean(var2));
        else
            trueCorr = calcCorr_thetas(theta(3),theta(4),theta(5));
        end
        for s=1:length(scale)
            
            Data=pcm_generateData_eva(M{4},theta,D,sn,scale(s),noise);
            
            
            T = pcm_fitModels(Data,M,D.partVec,D.condVec,runEffect,algorithm);
            C = pcm_correlation(Data,D.partVec,D.condVec,M{4},runEffect,M_type);
            
            C.trueCorr=ones(size(C.r_model2))*trueCorr;
            C.signalLevel=ones(size(C.r_model2))*scale(s);
            
            S=eval_simulation(trueCorr,C);
            S.signalLevel=[scale(s);scale(s);scale(s)];
            
            SS=addstruct(SS,S);
            TT=addstruct(TT,C);
            TT=addstruct(TT,T);
            
        end
        
        keyboard;
        % save structures
        save(fullfile(corrTestDir,sprintf('PCM_simulation_RS_%s_%s.mat',models{M_type},corrType)),'-struct','TT');
        save(fullfile(corrTestDir,sprintf('PCM_simulation_RS_%s_%s_summaryStats.mat',models{M_type},corrType)),'-struct','SS');
    
    case 'PLOT_pcm_logBayes'
        reg = [1,3];
        sessN=[1:4];
        seqType='trained';
        modelType='generic'; % or specific
        vararginoptions(varargin,{'reg','sessN','seqType','modelType'});
        
        T=[]; 
        for t=sessN
            R=load(fullfile(corrTestDir,sprintf('PCM_RS_corr_%s_%s_sess%d.mat',modelType,seqType,t)));
            R.sessN=ones(size(R.SN))*t;
            T=addstruct(T,R);
        end
        
        % excluding first - repsup model
        TT=T;
        TT.bayesEst(:,[2:4])=bsxfun(@minus,TT.bayesEst(:,[2:4]),TT.bayesEst(:,2));
        TT2.modelInd=[ones(size(TT.SN));ones(size(TT.SN))*2;ones(size(TT.SN))*3;ones(size(TT.SN))*4];
        TT2.bayesEst=TT.bayesEst(:);
        TT2.roi=[TT.roi;TT.roi;TT.roi;TT.roi];
        TT2.sessN=[TT.sessN;TT.sessN;TT.sessN;TT.sessN];
        TT2.SN=[TT.SN;TT.SN;TT.SN;TT.SN];
        
        % all models - with repsup
        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3;ones(size(T.SN))*4];
        T2.bayesEst=T.bayesEst(:);
        T2.roi=[T.roi;T.roi;T.roi;T.roi];
        T2.sessN=[T.sessN;T.sessN;T.sessN;T.sessN];
        T2.SN=[TT.SN;TT.SN;TT.SN;TT.SN];
        % rearranging how the data structure is arranged
        for ft=1:2
            if ft==1
                t=TT2;
            else
                t=T2;
            end
            figure
            sub=1;
            for r=reg
                subplot(1,numel(reg),sub)
                plt.line([t.sessN>3 t.sessN],t.bayesEst,'subset',t.modelInd>1&t.roi==r,'split',t.modelInd,'leg',{'repsup','repsup+seq','repsup+seq+corr'},'leglocation','north');
                plt.match('y');
                drawline(0,'dir','horz');
                title(sprintf('%s',regname{r}));
                if r==1
                    ylabel(sprintf('Log-Bayes %s',seqType));
                else
                    ylabel('');
                end
                sub=sub+1;
            end
            
        end
        
        for r=reg
            for ss=sessN
                b=getrow(t,t.roi==r&t.sessN==ss);
                fprintf('%s sess-%d sequence model \t',regname_cortex{r},ss);
                ttestDirect(b.bayesEst,[b.SN],2,'onesample','subset',b.modelInd==3);
                fprintf('%s sess-%d seq+corr model \t',regname_cortex{r},ss);
                ttestDirect(b.bayesEst,[b.SN],2,'onesample','subset',b.modelInd==4);
            end
        end
    case 'SAVE_pcm_corr_allSess'
        reg = [1,3];
        sessN=[1:4];
        seqType='trained';
        modelType='generic';
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','fitM'});
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                R=load(fullfile(corrTestDir,sprintf('PCM_RS_corr_%s_%s_sess%d.mat',modelType,seqType,t)));
                R.sessN=ones(size(R.SN))*t;
            %    R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end  
              
        %remove some fields
        r{1}='reg';r{2}='theta'; r{3}='theta_hat'; r{4}='thetaCr'; r{5}='thetaCr2';
        for rmIndx=1:size(r,2)
            T=rmfield(T,r{rmIndx});
        end
        
        save(fullfile(corrTestDir,sprintf('PCM_RS_corr_%s_%s_allSess.mat',modelType,seqType)),'-struct','T');
    case 'PLOT_corrTypes_scatter'
        reg = [1,3];
        modelType = 'specific';
        seqType = 'trained';
        sessN=[1:4];
        
        vararginoptions(varargin,{'reg','sessN','modelType'});
        
        T=load(fullfile(corrTestDir,sprintf('PCM_RS_corr_%s_%s_allSess.mat',modelType,seqType)));
        
        for ss=sessN
            for r=reg
                Tt=getrow(T,T.roi==r&T.sessN==ss);
                figure
                subplot(1,3,1)
                plt.scatter(Tt.r_naive,Tt.r_crossval);
                drawline(mean(Tt.r_crossval),'dir','horz','color',[1 0 0]);
                drawline(mean(Tt.r_naive),'dir','vert','color',[1 0 0],'lim',[-1 1]);
                drawline([-1,0,1],'dir','horz');drawline(0,'dir','vert','lim',[-1 1]);
                xlabel('Naive'); ylabel('Crossvalidated');
                title(sprintf('%s - %d session',regname{r},ss));
                
                subplot(1,3,2)
                plt.scatter(Tt.r_naive,real(Tt.r_crossval_noMakePD));
                drawline(mean(real(Tt.r_crossval_noMakePD)),'dir','horz','color',[1 0 0]);
                drawline(mean(Tt.r_naive),'dir','vert','color',[1 0 0],'lim',[-1 1]);
                drawline([-1,0,1],'dir','horz');drawline(0,'dir','vert','lim',[-1 1]);
                xlabel('Naive'); ylabel('Crossvalidated-notPD');
                title(sprintf('%s - %d session',regname{r},ss));
                
                subplot(1,3,3)
                plt.scatter(Tt.r_naive,Tt.r_model2);
                drawline(mean(Tt.r_model2),'dir','horz','color',[1 0 0]);
                drawline(mean(Tt.r_naive),'dir','vert','color',[1 0 0],'lim',[-1 1]);
                drawline([-1,0,1],'dir','horz');drawline(0,'dir','vert','lim',[-1 1]);
                xlabel('Naive'); ylabel('PCM model');
                title(sprintf('%s - %d session',regname{r},ss));
                
                plt.match('y');
            end
        end
    case 'PLOT_corr_sess_reg'
        roi = [1,3];
        modelType = 'specific';
        seqType = 'trained';
        sessN=[1:4];
        
        vararginoptions(varargin,{'reg','sessN','modelType'});
        
        T=load(fullfile(corrTestDir,sprintf('PCM_RS_corr_%s_%s_allSess.mat',modelType,seqType)));
             
        sub=1;
        for r=roi
            figure(1)
            subplot(1,3,1)
            plt.line([T.sessN>3 T.sessN],T.r_naive,'subset',T.roi==r);
            drawline(0,'dir','horz');
            title(sprintf('%s-naive',regname{r}));
            xlabel('Session'); ylabel('Correlation');
            subplot(1,3,2)
            plt.line([T.sessN>3 T.sessN],T.r_crossval,'subset',T.roi==r);
            drawline(0,'dir','horz');
            title(sprintf('%s-crossvalidated',regname{r}));
            ylabel('');
            subplot(1,3,3)
            plt.line([T.sessN>3 T.sessN],T.r_model2,'subset',T.roi==r);
            drawline(0,'dir','horz');
            title(sprintf('%s-pcm-%s',regname{r},modelType));
            ylabel('');
            
            plt.match('y');
            
            
            figure(2)
            subplot(1,numel(roi),sub)
            plt.line([T.sessN>3 T.sessN],T.r_crossval,'subset',T.roi==r);
            drawline(0,'dir','horz');
            title(sprintf('%s-crossval',regname{r}));
            xlabel('Session'); ylabel('Correlation');
            
            figure(3)
            subplot(1,numel(roi),sub)
            plt.line([T.sessN>3 T.sessN],T.r_model2,'subset',T.roi==r);
            drawline(0,'dir','horz');
            title(sprintf('%s-crossval',regname{r}));
            xlabel('Session'); ylabel('Correlation');
            
            sub=sub+1;
        end
        
        
    case 'PLOT_simulations'
        modelType='specific';
        corrType='positive';
        vararginoptions(varargin,{'modelType','corrType'});
        
 
        TT=load(fullfile(corrTestDir,sprintf('PCM_simulation_RS_%s_%s.mat',modelType,corrType)));
        SS=load(fullfile(corrTestDir,sprintf('PCM_simulation_RS_%s_%s_summaryStats.mat',modelType,corrType)));
        
        TT=rmfield(TT,'theta');
        
        gray=[120 120 120]/255;
        styDots=style.custom({gray},'markersize',10);
        styMean=style.custom({'blue'},'markersize',10);
        
        % correlation estimates - naive, crossval, model
        N=tapply(TT,{'signalLevel'},{'r_naive','mean'});
        figure
        subplot(1,3,1)
        plt.scatter(TT.signalLevel,TT.r_naive,'style',styDots);
        hold on;
        plt.scatter(N.signalLevel,N.r_naive,'style',styMean);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title(sprintf('Naive correlation'));
        xlabel('Signal level'); ylabel('Corr values');
        
        C=tapply(TT,{'signalLevel'},{'r_crossval','mean'});
        subplot(1,3,2)
        plt.scatter(TT.signalLevel,TT.r_crossval,'style',styDots);
        hold on;
        plt.scatter(C.signalLevel,C.r_crossval,'style',styMean);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title('Crossval correlation');
        xlabel('Signal level'); 
        
        M=tapply(TT,{'signalLevel'},{'r_model2','mean'});
        subplot(1,3,3)
        plt.scatter(TT.signalLevel,TT.r_model2,'style',styDots);
        hold on;
        plt.scatter(M.signalLevel,M.r_model2,'style',styMean);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title(sprintf('PCM correlation - %s',modelType));
        xlabel('Signal level');
        

                
        % summary stats
        leg_labels={'naive','crossval',sprintf('model %s simple',modelType),sprintf('model %s with shared seq pattern',modelType)};

        A=getrow(SS,SS.corrType==3);
        figure
        subplot(1,3,1)
        plt.line(A.signalLevel,A.bias,'split',A.corrType,'leg','off');
        drawline(0,'dir','horz');
        ylabel('Bias'); xlabel('Signal level'); title(sprintf('True correlation %d',A.trueCorr(1)));
        subplot(1,3,2)
        plt.line(A.signalLevel,A.var,'split',A.corrType,'leg','off');
        drawline(0,'dir','horz');
        ylabel('Variance'); xlabel('Signal level');
        subplot(1,3,3)
        plt.line(A.signalLevel,A.mse,'split',A.corrType,'leg',leg_labels,'leglocation','northeast');
        drawline(0,'dir','horz');
        ylabel('MSE'); xlabel('Signal level'); 
        
        
    case 'PLOT_corrDist_exe'
        roi=[1,3];
        subtract_mean=0;
        vararginoptions(varargin,{'roi','subtract_mean'});
        
        S=load(fullfile(repSupDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        
        Exe={'1st exe','2nd exe'};

        for r=1:numel(roi)
            figure(r)
            for exe=1:2
            subplot(1,2,exe)
            plt.line([S.sessN>3 S.sessN],S.corrDist,'style',stySeq,'split',S.seqType,'subset',S.FoSEx==exe&S.roi==roi(r),'leg',{'trained','untrained'});
            plt.match('y');
            title(sprintf('%s-%s',regname{roi(r)},Exe{exe}));
            if exe==1
            xlabel('Session'); ylabel('Corr distance');
            else
                ylabel('');
            end
            end
        end    
    case 'PLOT_dist_exe'
       betaChoice = 'multiPW';
       roi=[1,3];
       exe=[1:2];
       vararginoptions(varargin,{'betaChoice','roi','exe'});
       
       exe_label={'1st','2nd'};
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
       
       for r=1:numel(roi)
           figure(r)
           for e=exe
               subplot(1,2,e);
               plt.line([[S.speed;S.speed] [S.sessN;S.sessN]],[S.dist_train;S.dist_untrain],'split',[ones(size(S.sn));ones(size(S.sn))*2],'subset',[S.FoSEx;S.FoSEx]==e & [S.roi;S.roi]==roi(r),'style',stySeq,'leg',{'trained','untrained'});
               plt.match('y');
               if e==1
                   ylabel('Mahalanobis distance');
                   xlabel('Session');
               else
                   ylabel('');
               end
               title(sprintf('%s-%s',regname{roi(r)},exe_label{e}))
           end
       end    
    case 'STATS_dist_exe'
       betaChoice = 'multiPW';
       roi=[1,3];
       sessN=[1:4]
       vararginoptions(varargin,{'betaChoice','roi','exe'});
              
       T = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
       
       for r=roi
           for ss=sessN
               for exe=1:2
                   fprintf('t-test on %s sess-%d execution %d \t',regname{r},ss,exe);
                   ttestDirect(T.dist_train,[T.sn],2,'onesample','subset',T.roi==r&T.sessN==ss&T.FoSEx==exe);
               end
           end
       end
        
    case 'PLOT_pcm_corr'
        reg = [1,3];
        sessN=[1:4];
        seqType='trained';
        modelType='generic';
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','fitM'});
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                % R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_sess%d_%s_%s.mat',t,seqType{st},modelType)));
                R=load(fullfile(corrTestDir,sprintf('PCM_RS_corr_%s_%s_sess%d.mat',modelType,seqType,t)));
                R.sessN=ones(size(R.SN))*t;
            %    R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
%         corrType={'naive','crossval','pcm'};
%         corrVar=[T.r_naive, T.r_crossval, T.r_model2];
%         for ct=1:3 % three types of correlation
%             figure(ct)
%             sub=1;
%             for r=reg
%                 subplot(1,numel(reg),sub)
%                 plt.line([T.sessN>3 T.sessN],corrVar(:,ct),'subset',T.roi==r,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
%                 plt.match('y');
%                 drawline(0,'dir','horz');
%                 title(sprintf('%s',regname{r}));
%                 if sub==1
%                     ylabel(sprintf('%s corr',corrType{ct}));
%                     xlabel('Session');
%                 else
%                     ylabel('');
%                 end
%                 sub=sub+1;
%             end
%         end
        
        for ss=sessN
            for r=reg
               % Tt=getrow(T,T.roi==r&T.sessN==ss);
                figure
                subplot(1,3,1)
                plt.scatter(T.r_naive,T.r_crossval,'subset',T.roi==r&T.sessN==ss);
                drawline(mean(T.r_crossval),'dir','horz','color',[1 0 0]);
                drawline(mean(T.r_naive),'dir','vert','color',[1 0 0],'lim',[-1 1]);
                drawline([-1,0,1],'dir','horz');drawline(0,'dir','vert','lim',[-1 1]);
                xlabel('Naive'); ylabel('Crossvalidated');
                title(sprintf('%s - %d session',regname{r},ss));
                
                subplot(1,3,2)
                plt.scatter(T.r_naive,real(T.r_crossval_noMakePD),'subset',T.roi==r&T.sessN==ss);
                drawline(mean(real(T.r_crossval_noMakePD(T.roi==r&T.sessN==ss))),'dir','horz','color',[1 0 0]);
                drawline(mean(T.r_naive),'dir','vert','color',[1 0 0],'lim',[-1 1]);
                drawline([-1,0,1],'dir','horz');drawline(0,'dir','vert','lim',[-1 1]);
                xlabel('Naive'); ylabel('Crossvalidated-notPD');
                title(sprintf('%s - %d session',regname{r},ss));

                subplot(1,3,3)
                plt.scatter(T.r_naive,T.r_model2,'subset',T.roi==r&T.sessN==ss);
                drawline(mean(T.r_model2(T.roi==r&T.sessN==ss)),'dir','horz','color',[1 0 0]);
                drawline(mean(T.r_naive),'dir','vert','color',[1 0 0],'lim',[-1 1]);
                drawline([-1,0,1],'dir','horz');drawline(0,'dir','vert','lim',[-1 1]);
                xlabel('Naive'); ylabel('PCM model');
                title(sprintf('%s - %d session',regname{r},ss));
                
                plt.match('y');
            end
        end
end
end

function M = pcm_repsupModel_generic

    % for sequence-specific modelling - one parameter for all seq           
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'RepSup';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 4;
    M{3}.name       = 'RepSup+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    M{3}.Ac(:,3:8,3)  = [A;zeros(6)];      % Unique exe1 sequence patterns
    M{3}.Ac(:,9:14,4)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    
    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 5;
    M{4}.name         = 'RepSup+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    M{4}.Ac(:,3:8,3)  = [A;zeros(6)];      % Unique exe1 sequence patterns
    M{4}.Ac(:,9:14,4) = [zeros(6);A];     % Unique exe2 sequence pattterns
    M{4}.Ac(:,3:8,5)  = [zeros(6);A];     % Correlation exe1-exe2
    
end
function M = pcm_repsupModel_specific
    
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'RepSup';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 14;
    M{3}.name       = 'RepSup+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,3:8,2+i)  = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{3}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    end;

    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 20;
    M{4}.name         = 'RepSup+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,3:8,2+i)   = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{4}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
        M{4}.Ac(:,3:8,14+i)  = [zeros(6);A];     % Correlation exe1-exe2
    end;

end
function M = pcm_repsupModel_specific_genCorr
    
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'RepSup';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 14;
    M{3}.name       = 'RepSup+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,3:8,2+i)  = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{3}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    end;

    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 20;
    M{4}.name         = 'RepSup+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,3:8,2+i)   = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{4}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
        M{4}.Ac(:,3:8,14+i)  = [zeros(6);A];     % Correlation exe1-exe2
    end;
  
    % Model 5: Repetition suppression + generic correlation
    M{5}.type       = 'feature';
    M{5}.numGparams = 3;
    M{5}.name       = 'RepSup+GenCorr';
    M{5}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{5}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    M{5}.Ac(:,1,3) = [zeros(6,1);ones(6,1)];
end


function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Crossvalidated model comparision:
    [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    [Tcross2,thetaCr2] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);

    T.cross_likelihood = Tcross.likelihood;
    T.cross2_likelihood = Tcross2.likelihood;
    T.bayesEst = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
    T.bayesEst2 = bsxfun(@minus,T.cross2_likelihood,T.cross_likelihood(:,1));
    T.theta_hat=theta_hat;
    for t=1:size(thetaCr,2)
        tC{t} = thetaCr{t}';
        tC2{t} = thetaCr2{t}';
    end
    T.thetaCr=tC;
    T.thetaCr2=tC2;
end
function C = pcm_correlation(Data,partVec,condVec,M,runEffect,M_type)
sn=1:size(Data,2);

% make condition / partition vectors into structures per subject if not
% given by default
if size(condVec,2)==1
    cV=condVec; clear condVec;
    for p=sn
        condVec{p}=cV;
    end
end
if size(partVec,2)==1
    pV=partVec; clear partVec;
    for p=sn
        partVec{p}=pV;
    end
end
% --------------------------------------
% 1a. Empirical correlation
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec{p});
    b = pinv(Z)*Data{p};           % Estimate mean activities
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>6)-1);
    b(1:6,:)  = bsxfun(@minus,b(1:6,:) ,mean(b(1:6,:))); % Subtract mean per condition - first exe
    b(7:12,:) = bsxfun(@minus,b(7:12,:),mean(b(7:12,:))); % second exe
    G=cov(b');
    C.r_naive(p,1) = calcCorr(G);
end;

% 1b. Empirical correlation - subtract the mean for each run and condition
for p=sn
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>6)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    b=R*Data{p};
    % make mean pattern
    Z=pcm_indicatorMatrix('identity',condVec{p});
    b2 = pinv(Z)*b;           % Estimate mean activities
    G=cov(b2');
    C.r_naive2(p,1) = calcCorr(G);
end;

% --------------------------------------
% 2a. Crossvalidated correlation - make PD
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec{p});
    % Subtract mean for each condition and run
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>6)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    C.r_crossval(p,1)=calcCorr(pcm_makePD(Gcv(:,:,p)));
    C.r_crossval_noMakePD(p,1)=calcCorr(Gcv(:,:,p));
end;

% --------------------------------------
% 3. Fit model 2  and infer correlations from the parameters
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);
C.theta=theta;
% Get the correlations
switch M_type
    case 1
        var1       = (theta{1}(3,:).^2)';
        var2       = (theta{1}(4,:).^2+theta{1}(5,:).^2)';
        cov12      = (theta{1}(3,:).*theta{1}(5,:))';
       % var1       = (theta{1}(1,:).^2)';
       % var2       = (theta{1}(2,:).^2+theta{1}(3,:).^2)';
       % cov12      = (theta{1}(1,:).*theta{1}(3,:))';
        C.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
    case 2
        var1       = (theta{1}(3:8,:).^2)';
        var2       = (theta{1}(9:14,:).^2+theta{1}(15:20,:).^2)';
        cov12      = (theta{1}(3:8,:).*theta{1}(15:20,:))';
        %var1       = (theta{1}(1:6,:).^2)';
        %var2       = (theta{1}(7:12,:).^2+theta{1}(13:18,:).^2)';
        %cov12      = (theta{1}(1:6,:).*theta{1}(13:18,:))';
        C.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
    case 3
        var1       = (theta{1}(1,:).^2)';
        var2       = (theta{1}(2,:).^2+theta{1}(3,:).^2)';
        cov12      = (theta{1}(1,:).*theta{1}(3,:))';
       % var1       = (theta{1}(1,:).^2)';
       % var2       = (theta{1}(2,:).^2+theta{1}(3,:).^2)';
       % cov12      = (theta{1}(1,:).*theta{1}(3,:))';
        C.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
end
% --------------------------------------
end

function r=calcCorr(G)
d0=diag(G);
v1 = d0(1:6)';    % Variances first exe
v2 = d0(7:12)';   % Variances 2nd exe
cv=diag(G,6);     % Covariance
r = mean(cv)/sqrt(mean(v1)*mean(v2));
end
function r=calcCorr_thetas(th1,th2,th3)

    v1=th1.^2;
    v2=th2.^2+th3.^2;
    cv=th1.*th3;
    r=cv./sqrt(v1.*v2);
end
function S=eval_simulation(trueCorr,C)
    % summary evaluation statistics - comparing true correlation with
    % estimations in C
    
    S.trueCorr=[trueCorr; trueCorr; trueCorr];
    S.bias=[trueCorr-mean(C.r_naive);trueCorr-mean(C.r_crossval);trueCorr-mean(C.r_model2)];
    S.var=[var(C.r_naive);var(C.r_crossval);var(C.r_model2)];
    S.mse=[sum((C.trueCorr-C.r_naive).^2);
        sum((C.trueCorr-C.r_crossval).^2);
        sum((C.trueCorr-C.r_model2).^2)];
    S.corrType=[1;2;3]; % 1-naive, 2-crossval, 3-model
    
end