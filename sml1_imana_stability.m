function varargout=sml1_imana_stability(what,varargin)

% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
behavDir        =[baseDir '/behavioral_data/data'];            
imagingDir      =[baseDir '/imaging_data'];                   
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia'];
suitDir         =[baseDir '/suit'];
pcmDir          =[baseDir '/pcm_stats'];
stabDir         =[baseDir '/stability_stats'];
qualContrDir    =[baseDir '/qual_control'];

% update glmDir when adding new glms
glmLocDir       ={[baseDir '/glmLoc/glmL1'],[baseDir '/glmLoc/glmL2'],[baseDir '/glmLoc/glmL3']};   % localiser glm
glmLocSessDir   ={[baseDir '/glmLocSess/glmLocSess1'],[baseDir '/glmLocSess/glmLocSess2'],[baseDir '/glmLocSess/glmLocSess3'],[baseDir '/glmLocSess/glmLocSess4']}; % one glm for loc run per session
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session

% ------------------------- Experiment Info -------------------------------

% Stimuli - numbers given in SeqNumb
num_train = 1:6;
num_untrain = 7:12;   
num_seq = 1:12;
num_fing = 13:17;
num_seqtype = 2;
seqType_label = {'trained','untrained'};
% for group 1 - group 2 (1-6 and 7-12 reversed)

% per session
numruns_sess      = 10;  
numruns_task_sess = 8;
numruns_loc_sess  = 2;

% total - per subject (the total in the end will always be 40)
numruns           = [40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40];
numruns_task      = 32;
numruns_loc       = 8;

sess = [repmat(1,1,10),repmat(2,1,10),repmat(3,1,10),repmat(4,1,10)];   % all sessions

sess_sn = [4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4];    % per subject

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
% The variables in this section must be updated for every new subject.
subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18'};  

% Other random notes %

% -------------------------- For plotting ---------------------------------
stySeq=style.custom({'red','blue'},'markersize',12);

% ------------------------------ Analysis Cases --------------------------------
switch(what)
    
% make sure ROI_beta, ROI_stats and psc complete in sml1_imana_dist
    case 'beta_consistPattern_witSess'                               
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 1 - RSS/TSS (residual sum of squares per partition / total)
        
        % This stat is useful for determining which GLM model yields least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.
        %
        sessN = [1:4];
        sn  = [1:9,11,12];
        roi = [1:8];
        betaChoice = 'uw';  % raw / uw / mw 
        removeMean = 'yes'; % are we removing pattern means for patternconsistency?
        vararginoptions(varargin,{'sn','glm','roi','betaChoice','removeMean','sessN'});
        
        if strcmp(removeMean,'yes')
             rm = 1; % we are removing the mean
        else rm = 0; % we are keeping the mean (yeilds higher consistencies but these are biased)
        end
  
        Rreturn=[];
        C=[];
        CAll=[];
        %========%
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_sess%d.mat',ss))); % loads in struct 'T'
            for r=roi
                Rall=[]; %prep output variable
                for s=sn
                    S = getrow(T,(T.SN==s & T.region==r));
                    runs = 1:numruns_task_sess; % 8 func runs
                    switch(betaChoice)
                        case 'raw'
                            betaW  = S.betaRAW{1};
                        case 'uw'
                            betaW  = S.betaUW{1};
                        case 'mw'
                            betaW  = S.betaW{1};
                    end
                    
                    % make vectors for pattern consistency func
                    conditionVec = kron(ones(numel(runs),1),[1:12]');
                    partition    = kron(runs',ones(12,1));
                    % calculate the pattern consistency
                    R2   = rsa_patternConsistency(betaW,partition,conditionVec,'removeMean',rm);
                    Rall = [Rall,R2];
                    C.sn=s;
                    C.roi=r;
                    C.consist=R2;
                    C.sessN=ss;
                    CAll=addstruct(CAll,C);
                end
            end
        end

        dircheck(stabDir);
        % save the mat file
        save(fullfile(stabDir,sprintf('Consist_witSess_removeMean%d_betas_%s.mat',rm,betaChoice)),'-struct','CAll');
        %_______________    
        case 'PLOT_consist_witSess'
        betaChoice = 'uw';  % raw / uw / mw 
        removeMean = 1; % are we removing pattern means for patternconsistency?
        roi=[1:8];
        
        vararginoptions(varargin,{'roi','betaChoice','removeMean'});
        
        T = load(fullfile(stabDir,sprintf('Consist_witSess_removeMean%d_betas_%s.mat',removeMean,betaChoice)));
        
        
        figure;
        for s=1:numel(roi)
            subplot(1,numel(roi),s)
            lineplot(T.sessN,T.consist,'style_thickline','subset',T.roi==roi(s));
            title(sprintf('%s %s betas',regname{s},betaChoice));
            if s==1
                xlabel('Session'); ylabel('Consistency');
            else
                xlabel(''); ylabel('');
            end
        end
        case 'PLOT_Consist_QC_motion'
        % relate the consistency of patterns to amount of motion
        betaChoice = 'raw';
        removeMean=0;
        reg=[1:8];
        sn=[1:9,11,12];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','removeMean','reg','sn','sessN'});
        
        
        C = load(fullfile(stabDir,sprintf('Consist_witSess_removeMean%d_betas_%s.mat',removeMean,betaChoice)));
        Q=load(fullfile(qualContrDir,'QC_motion.mat'));
        
            for r=reg
                Cc=getrow(C,C.roi==r);
                
                figure
                subplot(1,2,1)
                scatterplot(Q.rms,Cc.consist,'split',Q.sess,'leg',{'sess1','sess2','sess3','sess4'},'regression','linear');
                xlabel('overall motion'); ylabel('pattern consistency of betas (1-RSS/TSS)');
                title(sprintf('%s',regname{r}));
                subplot(1,2,2)
                scatterplot(Q.fwd,Cc.consist,'split',Q.sess,'leg',{'sess1','sess2','sess3','sess4'},'regression','linear');
                xlabel('framewise displacement'); 
                title(sprintf('%s',regname{r}));
            end
        case 'PLOT_Consist_QC_motion_sn'
       % relate the consistency of patterns to amount of motion - per sn
        betaChoice = 'raw';
        removeMean=0;
        reg=[1:8];
        sn=[1:9,11,12];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','removeMean','reg','sn','sessN'});
        
        
        C = load(fullfile(stabDir,sprintf('Consist_witSess_removeMean%d_betas_%s.mat',removeMean,betaChoice)));
        Q=load(fullfile(qualContrDir,'QC_motion.mat'));
        
            for r=reg
                Cc=getrow(C,C.roi==r);
                
                figure
                subplot(1,2,1)
                scatterplot(Q.rms,Cc.consist,'split',Q.sn,'leg',subj_name(sn));
                xlabel('overall motion'); ylabel('pattern consistency of betas (1-RSS/TSS)');
                title(sprintf('%s',regname{r}));
                subplot(1,2,2)
                scatterplot(Q.fwd,Cc.consist,'split',Q.sn,'leg',subj_name(sn));
                xlabel('framewise displacement'); 
                title(sprintf('%s',regname{r}));
            end

    case 'beta_crossvalCorr_witSess'
        % Splits data for each session into two partitions (even and odd runs).
        % Calculates correlation coefficients between splits for each
        % condition pair - within subject & session.

        reg = [1:8]; 
        sn  = [1:9,11,12];
        sessN = [1:4];
        remove_mean = 1; % subtract
        betaChoice = 'mw'; % uw / mw / raw

        vararginoptions(varargin,{'roi','sn','remove_mean','sessN','betaChoice'});
                
        partitions = [1:2:numruns_task_sess; 2:2:numruns_task_sess];
        numRuns    = [1:numruns_task_sess];
        numConds   = num_seq;
        conds   = repmat([numConds],1,length(numRuns))';
        runNums = kron([numRuns],ones(1,length(numConds)))';
        
        C=[];
        for ss = sessN
            D   = load(fullfile(regDir,sprintf('betas_sess%d.mat',ss)));
            splitcorrs = [];
            for roi = reg;
                T   = getrow(D,D.region==reg(roi));
                
                for s = 1:numel(sn) % for each subject
                    sbeta=[]; prepBetas=[];
                    t = getrow(T,T.SN==sn(s));
                    
                    prepBetas = [];
                    switch (betaChoice)
                        case 'uw'
                            beta = t.betaUW{1};
                        case 'mw'
                            beta = t.betaW{1};
                        case 'raw'
                            beta = t.betaRAW{1};
                    end
                    
                    if remove_mean
                        for r=numRuns
                            sbeta(runNums==r,:) = bsxfun(@minus,beta(runNums==r,:),mean(beta(runNums==r,:)));
                        end
                    else
                        sbeta=beta;
                    end
                    % prep betas (harvest and subtract partition mean)
                    for i = 1:size(partitions,1)
                        partitionIdx = logical(ismember(runNums,partitions(i,:)));
                        condIdx{i}   = conds(partitionIdx);
                        prepBetas{i} = sbeta(partitionIdx,:);
                    end
                    
                    % correlate patterns across partitions, both within and across
                    % conditions
                    for c1 = numConds % for each condition
                        % condition mean activity pattern for the odd run partition
                        oddCon   = condIdx{1}==c1;
                        oddBetas = mean(prepBetas{1}(oddCon,:));
                        % condition mean activity pattern for the even run partition
                        evenCon   = condIdx{2}==c1;
                        evenBetas = mean(prepBetas{2}(evenCon,:));
                        % correlate condition patterns across partitions
                        tmp = corrcoef(evenBetas,oddBetas);
                        splitcorrs(c1) = tmp(1,2);
                    end
                    
                    Corr.sn = sn(s);
                    Corr.reg = roi;
                    Corr.trainCorr = mean(splitcorrs(:,[1:6]));
                    Corr.untrainCorr = mean(splitcorrs(:,[7:12]));
                    Corr.overallCorr = mean(splitcorrs);
                    Corr.sessN = ss;
                    
                    C = addstruct(C,Corr);         
                end
            end
        end

        % save
        save(fullfile(stabDir,sprintf('CrossvalCorr_witSess_removeMean%d_betas_%s.mat',remove_mean,betaChoice)),'-struct','C');   
        case 'PLOT_crossvalCorr_witSess'
        
        betaChoice = 'mw';
        remove_mean=1;
        reg=[1:8];
        
        vararginoptions(varargin,{'betaChoice','remove_mean','reg'})
        
        C=load(fullfile(stabDir,sprintf('CrossvalCorr_witSess_removeMean%d_betas_%s.mat',remove_mean,betaChoice)));
        
        figure
        for r=reg
            subplot(1,numel(reg),r)
            barplot(C.sessN, [C.trainCorr C.untrainCorr],'subset',C.reg==reg(r),'leg',{'train','untrain'});
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Crossvalidated correlation of beta patterns'); xlabel('Session');
            else
                ylabel(''); xlabel('');
            end
        end
        
        figure
        for i=1:numel(reg)
            subplot(1,numel(reg),i)
            lineplot([[C.sessN>3; C.sessN>3] [C.sessN; C.sessN]], [C.trainCorr; C.untrainCorr],'split',[ones(size(C.trainCorr));ones(size(C.trainCorr))*2],'style_thickline','leg',{'train','untrain'},'subset',[C.reg;C.reg]==i);
            title(sprintf('%s',regname{i}));
            if r==1
                ylabel('Crossvalidated correlation of beta patterns'); xlabel('Session');
            else
                ylabel(''); xlabel('');
            end
        end
        case 'PLOT_crossvalCorr_QC_motion'
        % relate the consistency of patterns to amount of motion
        betaChoice = 'raw';
        removeMean=0;
        reg=[1:8];
        sn=[1:9,11,12];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','removeMean','reg','sn','sessN'});
        
        
        C=load(fullfile(stabDir,sprintf('CrossvalCorr_witSess_removeMean%d_betas_%s.mat',removeMean,betaChoice)));      
        Q=load(fullfile(qualContrDir,'QC_motion.mat'));
        
            for r=reg
                Cc=getrow(C,C.reg==r);
                
                figure
                subplot(1,2,1)
                scatterplot(Q.rms,Cc.overallCorr,'split',Q.sess,'leg',{'sess1','sess2','sess3','sess4'},'regression','linear');
                xlabel('overall motion'); ylabel('crossvalidated corr of betas');
                title(sprintf('%s',regname{r}));
                subplot(1,2,2)
                scatterplot(Q.fwd,Cc.overallCorr,'split',Q.sess,'leg',{'sess1','sess2','sess3','sess4'},'regression','linear');
                xlabel('framewise displacement'); ylabel('crossvalidated corr of betas');
                title(sprintf('%s',regname{r}));
            end
        case 'PLOT_crossvalCorr_QC_motion_sn'
        % plot consistency of regressors against motion - per subject 
        %
        betaChoice = 'raw';
        removeMean=0;
        reg=[1:8];
        sn=[1:9,11,12];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','removeMean','reg','sn','sessN'});
        
        
        C=load(fullfile(stabDir,sprintf('CrossvalCorr_witSess_removeMean%d_betas_%s.mat',removeMean,betaChoice)));      
        Q=load(fullfile(qualContrDir,'QC_motion.mat'));
        
            for r=reg
                Cc=getrow(C,C.reg==r);
                
                figure
                subplot(1,2,1)
                scatterplot(Q.rms,Cc.overallCorr,'split',Q.sn,'leg',subj_name(sn));
                xlabel('overall motion'); ylabel('crossvalidated corr of betas');
                title(sprintf('%s',regname{r}));
                subplot(1,2,2)
                scatterplot(Q.fwd,Cc.overallCorr,'split',Q.sn,'leg',subj_name(sn));
                xlabel('framewise displacement'); ylabel('crossvalidated corr of betas');
                title(sprintf('%s',regname{r}));
            end
        
    case 'ROI_beta_consist_betwSess'
        % evaluate consistency of measures (psc, beta, z-scores) across
        % sessions for finger mapping
        
        sn  = 1;
        reg = 1:7;
        removeMean = 'yes'; % are we removing pattern means for patternconsistency?
        betaChoice = 'multi'; % options: uni / multi / raw
        seq = 'untrained';
        vararginoptions(varargin,{'sn','reg','betaChoice','removeMean','seq'});
        
        if strcmp(removeMean,'yes')
             rm = 1; % we are removing the mean
        else rm = 0; % we are keeping the mean (yeilds higher consistencies but these are biased)
        end
        
       for  roi = reg;
        CS=[];  % separate per sequence
        PS=[];  % across all sequences
        for sessN = 1:4; % per session
            C=[];P=[];
            T = load(fullfile(regDir,sprintf('betas_sess%d.mat',sessN))); % loads region data (T)
        
            switch (betaChoice)
            case 'uni'
                beta = T.betaUW;
            case 'multi'
                beta = T.betaW;
            case 'raw'
                beta = T.betaRAW;
            end
            
            runs=1:numruns_task_sess;
            conditionVec = kron(ones(numel(runs),1),[1:12]');
            split_run=1:4;
            splitRunVec = kron(ones(numel(split_run),1),[ones(12,1); ones(12,1).*2]); % split even and odd runs
            
            switch(seq)
                case 'trained'
                    idx=1:6;
                case 'untrained'
                    idx=7:12;
            end

            %C.beta=beta{roi};   
            for d = 1:6 %sequences
                C.beta_seq(d,:)=mean(beta{roi}(conditionVec==idx(d),:),1);  % beta values for each digit (avrg across blocks)
            end
            
            %C.zscore_seq = bsxfun(@rdivide,C.beta_seq,sqrt(T.resMS{roi}));
           
            C.seq_ind=[1:6]';
            C.sessN=ones(6,1)*sessN;
            C.roi=ones(6,1)*roi;

            P.beta_mean=mean(C.beta_seq,1);   % mean pattern acros trained / untrained sequences in each session
            P.sessN=sessN;
            P.roi=roi;
            
            CS=addstruct(CS,C);
            PS=addstruct(PS,P);
        end
        
        ind = indicatorMatrix('allpairs',([1:4]));  % betwSess indicator
        for n=1:numel(unique(CS.seq_ind))
            T = getrow(CS,CS.seq_ind==n);
            for i=1:size(ind,1)
                [i1 i2] = find(ind(i,:)~=0);
                if rm == 1
                    AcrSess_b(i)=corr(T.beta_seq(i1(2),:)',T.beta_seq(i2(2),:)');
                elseif rm == 0
                    AcrSess_b(i)=corrN(T.beta_seq(i1(2),:)',T.beta_seq(i2(2),:)');
                end
            end
            AcrSess_beta(n)=mean(AcrSess_b);

        end
        Consist.beta_corr(roi,:) = AcrSess_beta;
        Consist.roi(roi,1) = roi;
        
       end
       
        figure(1)
        col=hsv(7);
        for i = reg
            a(i)=plot(Consist.beta_corr(i,:),'-o','Color',col(i,:));
            hold on;
        end
        title('Beta values')
        legend(a,regname(reg));
        xlabel('All across-session combinations');
        ylabel('Correlation of sequence-specific pattern');

        keyboard;  
    case 'ROI_beta_consist_betwSess_LOC'
        % evaluate consistency of measures (psc, beta, z-scores) across
        % sessions for finger mapping
        
        sn  = 1;
        sessN = 1;
        reg = 1:7;
        keepmean=0;
        betaChoice = 'multi'; % options: uni / multi / raw
        vararginoptions(varargin,{'sn','sessN','reg','betaChoice','keepmean'});

        
       for  roi = reg;
        CS=[];  % separate per digit
        PS=[];  % across all digits
        for sessN = 1:4; % per session
            C=[];P=[];
            T = load(fullfile(regDir,sprintf('betas_LOC_sess%d.mat',sessN))); % loads region data (T)
        
            switch (betaChoice)
            case 'uni'
                beta = T.betaUW;
            case 'multi'
                beta = T.betaW;
            case 'raw'
                beta = T.betaRAW;
            end
        
            %C.beta=beta{roi};   % 12 betas x voxels (1-5, 1-5; 2 intercept)
            for d = 1:5 %digits
                C.beta_digit(d,:)=mean(beta{roi}([d,d+6],:),1);  % beta values for each digit (avrg across two blocks)
               % C.psc_digit(d,:)=mean(median(max(beta{roi}(d,:)))/median(max(beta{roi}(6,:))),...
               %     median(max(beta{roi}(d+6,:)))/median(max(beta{roi}(12,:))));  % mean of psc of two blocks - median response / intercept
            end
            
            C.zscore_digit = bsxfun(@rdivide,C.beta_digit,sqrt(T.resMS{roi}));
            C.digit_ind=[1:5]';
            C.sessN=ones(5,1)*sessN;
            C.roi=ones(5,1)*roi;
            
            
            P.beta_mean=mean(C.beta_digit,1);   % mean pattern acros digit in each session
            P.zscore_mean=bsxfun(@rdivide,P.beta_mean,sqrt(T.resMS{roi}));
            P.sessN=sessN;
            P.roi=roi;
            
            CS=addstruct(CS,C);
            PS=addstruct(PS,P);
        end
        
        O.betas(roi,:) = mean(PS.beta_mean,2)';    % one value per session
        O.zscore(roi,:) = mean(PS.zscore_mean,2)';
        O.roi(roi,1) = roi;
        
        ind = indicatorMatrix('allpairs',([1:4]));  % betwSess indicator
        for i=1:size(ind,1)
            [i1 i2] = find(ind(i,:)~=0);
            if keepmean == 0
                Consist.beta(roi,i)=corr(PS.beta_mean(i1(2),:)',PS.beta_mean(i2(2),:)');
                Consist.zscore(roi,i)=corr(PS.zscore_mean(i1(2),:)',PS.zscore_mean(i2(2),:)');
            elseif keepmean == 1
                Consist.beta(roi,i)=corrN(PS.beta_mean(i1(2),:)',PS.beta_mean(i2(2),:)');
                Consist.zscore(roi,i)=corrN(PS.zscore_mean(i1(2),:)',PS.zscore_mean(i2(2),:)');
            end
        end
        
        Consist.beta_RSA(roi,1) = rsa_patternConsistency(CS.beta_digit,CS.sessN,CS.digit_ind,'removeMean',keepmean);
        Consist.zscore_RSA(roi,1) = rsa_patternConsistency(CS.zscore_digit,CS.sessN,CS.digit_ind,'removeMean',keepmean);
        Consist.roi(roi,1) = roi;
        
       end
       
        figure(1)
        col=hsv(7);
        for i = reg
            a(i)=plot(Consist.beta(i,:),'-o','Color',col(i,:));
            hold on;
            drawline(Consist.beta_RSA(i),'dir','horz','color',col(i,:));
        end
        title('Beta values')
        legend(a,regname(reg));
        xlabel('All across-session combinations');
        ylabel('Correlation / RSA consistency(line)')
        
        figure(2)
        for j=reg
            b(j)=plot(Consist.zscore(j,:),'-o','Color',col(j,:));
            hold on;
            drawline(Consist.zscore_RSA(j),'dir','horz','color',col(j,:));
        end
        title('Z scores')
        legend(b,regname(reg));
        xlabel('All across-session combinations');
        ylabel('Correlation / RSA consistency(line)');
        
        keyboard;
    case 'ROI_beta_consist_witSess_LOC'
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % This stat is useful for determining which GLM model yields least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.
        %
        % enter sn, region, glm #, beta: 0=betaW, 1=betaUW, 2=raw betas
        % (1) Set parameters
        sessN = 1;
        sn  = 1;
        roi = 2; % default LH primary motor cortex
        betaChoice = 'uw';  % raw / uw / mw 
        keepmean = 0; % are we removing pattern means for patternconsistency?
        vararginoptions(varargin,{'sn','glm','roi','betaChoice','keepmean','sessN'});

        Rreturn=[];
        %========%
        for s=sessN
            T = load(fullfile(regDir,sprintf('betas_LOC_sess%d.mat',s))); % loads in struct 'T'
            for r=roi
                Rall=[]; %prep output variable
                for s=sn
                    S = getrow(T,(T.SN==s & T.region==r));
                    runs = 1:numruns_loc_sess; % 2 func runs
                    switch(betaChoice)
                        case 'raw'
                            betaW  = S.betaW{1}; 
                        case 'uw'
                            betaW  = S.betaUW{1}; 
                        case 'mw'
                            betaW  = S.betaRAW{1}; 
                    end
                    
                    % make vectors for pattern consistency func
                    conditionVec = kron(ones(numel(runs),1),[1:5]');
                    partition    = kron(runs',ones(5,1));
                    % calculate the pattern consistency
                    R2   = rsa_patternConsistency(betaW,partition,conditionVec,'removeMean',keepmean);
                    Rall = [Rall,R2];
                end
                Rreturn = [Rreturn;Rall];
            end
        end
        varargout = {Rreturn};
        fprintf('The consistency for %s betas in region %s is',betaChoice,regname{roi});
        % output arranged such that each row is an roi, each col is subj
    
    case 'reliability_betweenSess' 
        % collects data of interest, submits it to the main case
        reg = [1:8];
        sn  = [1:9,11:18];
        sessN = [1:2];
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        regType='Brodmann';
        % options for regularisation 'reg'
        % - minvalue (if variance <0.001 -> 0.001
        % - zero (if variance negative, make it 0)
        
        vararginoptions(varargin,{'sn','reg','sessN','regcorrType','seqTypeCorr','regType'});
        CAll = [];
        SeqT = [];
        for s = sn;
            for roi = reg;
                for  ss = 1:numel(sessN)
                    D   = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',regType,sessN(ss))));
                    t   = getrow(D,D.region==roi & D.SN==s);
                    data{ss} = t.betaUW{1};
                    C=[];
                end
                % send data to another case
                
                switch(seqTypeCorr)
                    case 'seqSpec' % sequence-specific
                        C = sml1_imana_stability('reliability_seqcorr','data',data,'removeMean',1,'numruns',numruns_task_sess,'reg',regcorrType);    % subtract mean, num of runs
                    case 'seqType'  % sequence-type (average across 6 sequences)
                         C = sml1_imana_stability('reliability_seqtypecorr','data',data,'numruns',numruns_task_sess,'reg',regcorrType);
                end
                 C.roi = ones(size(C.w1)).*roi;
                 C.sn = ones(size(C.w1)).*s;
                 CAll=addstruct(CAll,C);
                 
            end
        end
        
        save(fullfile(stabDir,sprintf('Stability_%s_sess%d-sess%d_regular_%s',seqTypeCorr,sessN(1),sessN(2),regcorrType)),'-struct','CAll');
        case 'reliability_simulate'
        % create different data types for simulating correlations within /
        % between sessions
        type = 'corr+noise';
        run_type = '2runs';
        vararginoptions(varargin,{'type','run_type'});
        P = 100;    % number of voxels simulated
        cond = 12;  % number of conditions
        switch (run_type)
            case '2runs'
                run=2;
            case '8runs'
                run=8;
        end
        
        switch(type)
            case 'random'
                % both datasets completely random
                % no correlation within / across runs
                data{1} = randn(run*cond,P);
                data{2} = randn(run*cond,P);      
            case 'corr1' % perfect correlation
                data1 = randn(cond,P);
                data{1}=repmat(data1,[run,1]);
                data{2}=data{1};
            case 'corr1+noise' % add noise levels on top of perfect correlation
                data1 = randn(cond,P);
                noise_lev = 0.80;   % amount of noise level
                
                true_X = repmat(data1,[run,1]);
                noise_X  = randn(cond,P).*noise_lev;

                data{1} = true_X + noise_X;
                data{2} = true_X + noise_X;
            case 'corr+noise' % lower correlation + noise
                
                % decide on the weights
                th.a = 0.3;  % common pattern
                th.b1 = 0.2; % trained
                th.b2 = 0.2; % untrained
                th.c1 = 0.6; % session 1
                th.c2 = 0.6; % session 2
                th.d1 = 0.3; % run1
                th.d2 = 0.3;
                th.d3 = 0.3;
                th.d4 = 0.3;
                
                % underlying patterns
                X_common = randn(6,P)*th.a;
                X_trained = randn(6,P)*th.b1;
                X_untrained = randn(6,P)*th.b2;
                X_sess1 = randn(6,P)*th.c1;
                X_sess2 = randn(6,P)*th.c2;
                X_run1 = randn(6,P)*th.d1;
                X_run2 = randn(6,P)*th.d2;
                X_run3 = randn(6,P)*th.d3;
                X_run4 = randn(6,P)*th.d4;
                
                % data sess1
                data{1}(1:6,:) = X_common+X_trained+X_sess1+X_run1;
                data{1}(7:12,:) = X_common+X_untrained+X_sess1+X_run1;
                data{1}(13:18,:) = X_common+X_trained+X_sess1+X_run2;
                data{1}(19:24,:) = X_common+X_untrained+X_sess1+X_run2;
                
                % data sess2
                data{2}(1:6,:) = X_common+X_trained+X_sess2+X_run3;
                data{2}(7:12,:) = X_common+X_untrained+X_sess2+X_run3;
                data{2}(13:18,:) = X_common+X_trained+X_sess2+X_run4;
                data{2}(19:24,:) = X_common+X_untrained+X_sess2+X_run4;
                
                % send to different function - calculate var / cov
                S = analytic_cov(th);

        end
        CAll = [];
        
        for perm = 1:1000
            C = sml1_imana('reliability_seqcorr',data,1,run);
            C.perm = ones(size(C.w1)).*perm;
            CAll=addstruct(CAll,C);
        end
        CAll.dif = CAll.acr - CAll.geoMean; 
        % CAll.dif - positive when rel across sess bigger than geometric
        % mean of reliability of each session
        figure;
        histogram(CAll.dif,'Normalization','probability');
        figure;
        barplot(CAll.seqType,[CAll.w1 CAll.w2 CAll.acr CAll.geoMean]);
        keyboard;
        case 'reliability_seqcorr'
        
        vararginoptions(varargin,{'data','removeMean','numruns','reg'});
        
        sess = size(data,2);
        partitions = [1:2:numruns; 2:2:numruns];
        numRuns    = 1:numruns;
        numConds   = num_seq;
        conds   = repmat([numConds],1,length(numRuns))';
        runNums = kron([numRuns],ones(1,length(numConds)))';
        
        for ss = 1:sess
            
            % subtract mean per run
            if removeMean
                for r=numRuns
                    data{ss}(runNums==r,:) = bsxfun(@minus,data{ss}(runNums==r,:),mean(data{ss}(runNums==r,:)));
                end
            else
                data{ss}=data{ss};
            end
            
            % split datas per partition
            for i = 1:size(partitions,1)
                partitionIdx = logical(ismember(runNums,partitions(i,:)))';
                condIdx{i}   = conds(partitionIdx);
                prepBetas{ss}{i} = data{ss}(partitionIdx,:); % session / partition
            end
            
            % calculate average pattern per partition
            for c1 = numConds % for each condition
                % condition mean activity pattern for this run partition
                oddCon   = condIdx{1}==c1;
                oddBetas{ss}(c1,:) = mean(prepBetas{ss}{1}(oddCon,:));
                % condition mean activity pattern for the other run partition
                evenCon   = condIdx{2}==c1;
                evenBetas{ss}(c1,:) = mean(prepBetas{ss}{2}(evenCon,:));
            end
        end
        
        % correlate patterns across partitions, both within and across
        % conditions
        for c1 = numConds % for each condition
            % correlate within sessions
            W1 = corrcoef(oddBetas{1}(c1,:),evenBetas{1}(c1,:));
            W2 = corrcoef(oddBetas{2}(c1,:),evenBetas{2}(c1,:));
            Acr1 = corrcoef(oddBetas{1}(c1,:),oddBetas{2}(c1,:));
            Acr2 = corrcoef(evenBetas{1}(c1,:),evenBetas{2}(c1,:));
            Acr3 = corrcoef(oddBetas{1}(c1,:),evenBetas{2}(c1,:));
            Acr4 = corrcoef(oddBetas{2}(c1,:),evenBetas{1}(c1,:));
            % correlate condition patterns across partitions
            C.w1(c1,:) = W1(1,2);
            C.w2(c1,:) = W2(1,2);
            C.acr(c1,:) = mean([Acr1(1,2),Acr2(1,2),Acr3(1,2),Acr4(1,2)]);
        end
        
        for cond = numConds
            % if one of the correlations is negative, make geoMean 0
            
            switch(reg)
                case 'zero'
                    % if one of the correlations is negative, make all data 0
                    C.geoMean(cond,:)=sqrt(C.w1(cond).*C.w2(cond));
                    if (C.w1(cond)<0 | C.w2(cond)<0)
                        C.w1(cond,:)=0;
                        C.w2(cond,:)=0;
                        C.geoMean(cond,:)=0;
                        C.acr(cond,:)=0;
                    end
                case 'minvalue'
                    if (C.w1(cond)<0.001)
                        C.w1(cond,:)=0.001;
                        C.acr(cond,:)=0.001;
                    end
                    if (C.w2(cond)<0.001)
                        C.w2(cond,:)=0.001;
                        C.acr(cond,:)=0.001;
                    end
                    C.geoMean(cond,:)=sqrt(C.w1(cond).*C.w2(cond));
            end
            
        end
        C.acrSess_corr = C.acr./C.geoMean; % correlation across sessions
        C.acrSess_corr(isnan(C.acrSess_corr))=0;
        % corrected by geometric mean
        C.seqType(1:6,:)=1;
        C.seqType(7:12,:)=2;
        
      
        varargout{1}=C;
        case 'reliability_seqtypecorr'
        vararginoptions(varargin,{'data','numruns','reg'});
        
        sess = size(data,2);
        partitions = [1:2:numruns_task_sess; 2:2:numruns_task_sess];
        numRuns    = 1:numruns_task_sess;
        numConds   = num_seq; 
        conds   = repmat([numConds],1,length(numRuns));
        cond_type = repmat([ones(6,1); ones(6,1)*2]',1,length(numRuns));
        runNums = kron([numRuns],ones(1,length(numConds)));
        
        % for constructing average patterns
        numConds_new = 1:num_seqtype; 
        cond_new = repmat([numConds_new],1,length(numRuns));
        runNums_new = kron([numRuns],ones(1,length(numConds_new)));
        
        for ss = 1:sess
            
            % calculate average trained / untrained pattern
            for r=numRuns
                for s=numConds_new
                data2{ss}(runNums_new==r & cond_new==s,:) = mean(data{ss}(runNums==r & cond_type==s,:),1); % 6 sequences
                end
            end
 
            % split datas per partition
            for i = 1:size(partitions,1)
                partitionIdx = logical(ismember(runNums_new,partitions(i,:)))';
                condIdx{i}   = cond_new(partitionIdx);
                prepBetas{ss}{i} = data2{ss}(partitionIdx,:); % session / partition
            end
            
            % calculate average pattern per partition
            for c1 = numConds_new % for each condition
                % condition mean activity pattern for this run partition
                oddCon   = condIdx{1}==c1;
                oddBetas{ss}(c1,:) = mean(prepBetas{ss}{1}(oddCon,:));
                % condition mean activity pattern for the other run partition
                evenCon   = condIdx{2}==c1;
                evenBetas{ss}(c1,:) = mean(prepBetas{ss}{2}(evenCon,:));
            end
        end
        
        % correlate patterns across partitions, both within and across
        % conditions
        for c1 = numConds_new % for each seqType
            % correlate within sessions
            W1 = corrcoef(oddBetas{1}(c1,:),evenBetas{1}(c1,:));
            W2 = corrcoef(oddBetas{2}(c1,:),evenBetas{2}(c1,:));
            Acr1 = corrcoef(oddBetas{1}(c1,:),oddBetas{2}(c1,:));
            Acr2 = corrcoef(evenBetas{1}(c1,:),evenBetas{2}(c1,:));
            Acr3 = corrcoef(oddBetas{1}(c1,:),evenBetas{2}(c1,:));
            Acr4 = corrcoef(oddBetas{2}(c1,:),evenBetas{1}(c1,:));
            % correlate condition patterns across partitions
            C.w1(c1,:) = W1(1,2);
            C.w2(c1,:) = W2(1,2);
            C.acr(c1,:) = mean([Acr1(1,2),Acr2(1,2),Acr3(1,2),Acr4(1,2)]);
        end
        
        for cond = numConds_new
             switch(reg)
                case 'zero'
                    % if one of the correlations is negative, make all data 0
                    C.geoMean(cond,:)=sqrt(C.w1(cond).*C.w2(cond));
                    if (C.w1(cond)<0 | C.w2(cond)<0)
                        C.w1(cond,:)=0;
                        C.w2(cond,:)=0;
                        C.geoMean(cond,:)=0;
                        C.acr(cond,:)=0;
                    end
                    
                case 'minvalue'
                    if (C.w1(cond)<0.001)
                        C.w1(cond,:)=0.001;
                    end
                    if (C.w2(cond)<0.001)
                        C.w2(cond,:)=0.001;
                    end
                    C.geoMean(cond,:)=sqrt(C.w1(cond).*C.w2(cond));
            end

        end
        C.seqType=[1 2]';
        C.acrSess_corr = C.acr./C.geoMean; % correlation across sessions
        C.acrSess_corr(isnan(C.acrSess_corr))=0;
        varargout{1}=C;
      
        case 'PLOT_reliability'
        
        sessN= [1:2];
        reg = [1:8];
        regcorrType='zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        vararginoptions(varargin,{'sessN','reg','regcorrType','seqTypeCorr'});
        
        R=load(fullfile(stabDir,sprintf('Stability_%s_sess%d-sess%d_regular_%s',seqTypeCorr,sessN(1),sessN(2),regcorrType)));

        figure;
        for r = 1:numel(reg)
            subplot(1,numel(reg),r)
            barplot(R.seqType,[R.w1 R.w2 R.acr R.geoMean],'subset',R.roi==reg(r),'leg',{'within sess1','within sess2','across sess','geometric mean'});
            if r==1
                ylabel('correlation across sessions');
            else
                ylabel('');
            end
            title(sprintf(regname{r}));
        end
        
        
        figure;
        for r = 1:numel(reg)
            subplot(1,numel(reg),r)
            hold on
            barplot(R.seqType,R.acrSess_corr,'subset',R.roi==reg(r));
            
            set(gca,'XTick',[1 2],'XtickLabel',{'trained','untrained'});
            if r==1
                ylabel('correlation - corrected for geoMean');
            else
                ylabel('');
            end
            title(sprintf(regname{r}));
        end
        case 'PLOT_reliability_allSess'
        reg = [1:8];
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        sessTr=[1:3]; % which session transitions to include
        %1: sess1-2
        %2: sess2-3
        %3: sess3-4
        sessName={'1-2','2-3','3-4'};
        
        vararginoptions(varargin,{'reg','regcorrType','seqTypeCorr','sessTr'});
        T=[];
        for t=1:numel(sessTr)
            R=load(fullfile(stabDir,sprintf('Stability_%s_sess%d-sess%d_regular_%s',seqTypeCorr,t,t+1,regcorrType)));
            R.sessTr=ones(size(R.w1))*t;
            T=addstruct(T,R);
        end
        
        figure
        for r=1:numel(reg)
            subplot(1,numel(reg),r)
            if numel(sessTr)==3
                 %plt.line([T.sessTr>2 T.sessTr],T.acr,'split',T.seqType,'style',stySeq,'subset',T.roi==r,'leg',{'trained','untrained'});
                plt.line([T.sessTr>2 T.sessTr],T.acrSess_corr,'split',T.seqType,'style',stySeq,'subset',T.roi==r,'leg',{'trained','untrained'});
            else
                %plt.line(T.sessTr,T.acr,'split',T.seqType,'style',stySeq,'subset',T.roi==r,'leg',{'trained','untrained'});
               plt.line(T.sessTr,T.acrSess_corr,'split',T.seqType,'style',stySeq,'subset',T.roi==r,'leg',{'trained','untrained'});
            end
            
            plt.match('y');
            set(gca,'XTickLabel',sessName(sessTr));
            drawline(0,'dir','horz');
            if r==1
                xlabel('Sess transitions');
                ylabel(sprintf('Correlation of %s pattern',seqTypeCorr));
            else
                ylabel('');
            end
            title(sprintf('%s',regname{r}));
        end
        case 'STATS_reliability'
        roi = [1:8];
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        sessTr=[1:2]; % which session transitions to include
        %1: sess1-2
        %2: sess2-3
        %3: sess3-4
        
        vararginoptions(varargin,{'reg','regcorrType','seqTypeCorr','sessTr'});
        T=[];
        for t=sessTr
            R=load(fullfile(stabDir,sprintf('Stability_%s_sess%d-sess%d_regular_%s',seqTypeCorr,t,t+1,regcorrType)));
            R.sessTr=ones(size(R.w1))*t;
            T=addstruct(T,R);
        end
        
        % Anova - sessTransition x seqType - per roi
        for r=roi
            fprintf('\n sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.acrSess_corr,T.sn,'within',[T.sessTr T.seqType],{'sessTrans','seqType'},'subset',T.roi==r);
        end
        
        % Post-hoc t-tests for effect of seqType per sessTransition
        for r=roi
            for ss=sessTr
                fprintf('\n post-hoc t-test on the effect of seqType in sessTr %d in %s \n',ss,regname{r});
                ttestDirect(T.acrSess_corr,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessTr==ss);
            end
        end
        
         % Post-hoc t-tests for effect of sessTransition per seqType
        for r=roi
            for ss=1:numel(seqType_label)
                fprintf('\n post-hoc t-test on the effect of sessTr in seqType %s in %s \n',seqType_label{ss},regname{r});
                ttestDirect(T.acrSess_corr,[T.sessTr T.sn],2,'paired','subset',T.roi==r&T.seqType==ss);
            end
        end
        
        case 'PLOT_reliability_speed'
        reg = [1:8];
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        sessTr=3; % which session transitions to include
        %3: sess3-4
        
        vararginoptions(varargin,{'reg','regcorrType','seqTypeCorr'});
        T=load(fullfile(stabDir,sprintf('Stability_%s_sess%d-sess%d_regular_%s',seqTypeCorr,sessTr,sessTr+1,regcorrType)));
        
        figure
        plt.bar(T.roi,T.acrSess_corr,'split',T.seqType,'style',stySeq,'leg',{'trained','untrained'});
        set(gca,'XTick',[1.5 4.5 8 11.5 14 17.5 21 24],'XTickLabel',regname);
        drawline(0,'dir','horz');
        ylabel(sprintf('Correlation of %s pattern',seqTypeCorr));
        case 'STATS_reliability_speed'
        reg = [1:8];
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        sessTr=3; % which session transitions to include
        %3: sess3-4
        
        vararginoptions(varargin,{'reg','regcorrType','seqTypeCorr'});
        T=load(fullfile(stabDir,sprintf('Stability_%s_sess%d-sess%d_regular_%s',seqTypeCorr,sessTr,sessTr+1,regcorrType)));
        for r=reg
            fprintf('\n t-test on stability of sequence type across sess 3-4 in %s \n',regname{r})
            ttestDirect(T.acrSess_corr,[T.seqType T.sn],2,'paired','subset',T.roi==r);
        end
        
        keyboard;
   
    case 'stability_ceiling'
      % collects data of interest, submits it to the main case
        reg = [1:8];
        sn  = [1:9,11:18];
        sessN = [1:2];
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        regType='Brodmann';
        seqType='trained';
        % options for regularisation 'reg'
        % - minvalue (if variance <0.001 -> 0.001
        % - zero (if variance negative, make it 0)
        
        vararginoptions(varargin,{'sn','reg','sessN','seqTypeCorr','regType'});
        CAll = [];
        switch seqType
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        for s = sn;
            for roi = reg;
                for  ss = 1:numel(sessN)
                    D   = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',regType,sessN(ss))));
                    glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{s});
                    T   = load(fullfile(glmDirSubj,'SPM_info.mat'));
                    t   = getrow(D,D.region==roi & D.SN==s);
                    data{ss} = t.betaUW{1}(T.seqType==st,:);
                    C=[];
                end
                % send data to another case
                
                switch(seqTypeCorr)
                    case 'seqSpec' % sequence-specific
                        C = sml1_imana_stability('stability_seqcorr','data',data,'removeMean',1,'numruns',numruns_task_sess);    % subtract mean, num of runs
                    case 'seqType'  % sequence-type (average across 6 sequences)
                         C = sml1_imana_stability('stability_seqtypecorr','data',data,'numruns',numruns_task_sess);
                end
                 C.roi = ones(size(C.w1)).*roi;
                 C.sn = ones(size(C.w1)).*s;
                 CAll=addstruct(CAll,C);          
            end
        end
        
        save(fullfile(stabDir,sprintf('StabilityCeiling_%s_%s_sess%d-sess%d',seqType,seqTypeCorr,sessN(1),sessN(2))),'-struct','CAll');  
        case 'stability_seqcorr'
        removeMean=1;
        vararginoptions(varargin,{'data','removeMean','numruns'});
        
        sess = size(data,2);
        partitions = [1:2:numruns; 2:2:numruns];
        numRuns    = 1:numruns;
        numConds   = 1:6; % only trained or untrained
        conds   = repmat([numConds],1,length(numRuns))';
        runNums = kron([numRuns],ones(1,length(numConds)))';
        
        for ss = 1:sess 
            % subtract mean per run
            if removeMean
                for r=numRuns
                    data{ss}(runNums==r,:) = bsxfun(@minus,data{ss}(runNums==r,:),mean(data{ss}(runNums==r,:)));
                end
            else
                data{ss}=data{ss};
            end           
            % split datas per partition
            for i = 1:size(partitions,1)
                partitionIdx = logical(ismember(runNums,partitions(i,:)))';
                condIdx{i}   = conds(partitionIdx);
                prepBetas{ss}{i} = data{ss}(partitionIdx,:); % session / partition
            end            
            % calculate average pattern per partition
            for c1 = numConds % for each condition
                % condition mean activity pattern for this run partition
                oddCon   = condIdx{1}==c1;
                oddBetas{ss}(c1,:) = mean(prepBetas{ss}{1}(oddCon,:));
                % condition mean activity pattern for the other run partition
                evenCon   = condIdx{2}==c1;
                evenBetas{ss}(c1,:) = mean(prepBetas{ss}{2}(evenCon,:));
            end
        end
        
        % correlate patterns across partitions, both within and across
        % conditions
        C.w1=corr(oddBetas{1}(:),evenBetas{1}(:));
        C.w2=corr(oddBetas{2}(:),evenBetas{2}(:));
        Acr1 = corr(oddBetas{1}(:),oddBetas{2}(:));
        Acr2 = corr(evenBetas{1}(:),evenBetas{2}(:));
        Acr3 = corr(oddBetas{1}(:),evenBetas{2}(:));
        Acr4 = corr(oddBetas{2}(:),evenBetas{1}(:));
        C.acr=mean([Acr1 Acr2 Acr3 Acr4]);
        C.geoMean=sqrt(C.w1*C.w2);
        
        if C.w1<0 || C.w2<0
            C.geoMean=NaN;
        end

        C.acrSess_corr = C.acr./C.geoMean; % correlation across sessions
        C.acrSess_corr(isnan(C.acrSess_corr))=0;
        % corrected by geometric mean
      
        varargout{1}=C;
        case 'stability_seqtypecorr'
        vararginoptions(varargin,{'data','numruns'});
        
        sess = size(data,2);
        partitions = [1:2:numruns_task_sess; 2:2:numruns_task_sess];
        numRuns    = 1:numruns_task_sess;
        numConds   = 1:6; 
        runNums = kron([numRuns],ones(1,length(numConds)));
        
        % for constructing average patterns
        runNums_new = 1:8;
        
        for ss = 1:sess
            
            % calculate average trained / untrained pattern
            for r=numRuns
                data2{ss}(runNums_new==r,:) = mean(data{ss}(runNums==r,:),1); % 6 sequences  
            end
            % split datas per partition
            for i = 1:size(partitions,1)
                partitionIdx = logical(ismember(runNums_new,partitions(i,:)))';
                prepBetas{ss}{i} = data2{ss}(partitionIdx,:); % session / partition
            end
            % calculate average pattern per partition
            % condition mean activity pattern for this run partition
            oddBetas{ss}(1,:) = mean(prepBetas{ss}{1});
            % condition mean activity pattern for the other run partition
            evenBetas{ss}(1,:) = mean(prepBetas{ss}{2});  
        end
        
        % correlate patterns across partitions, both within and across
        % conditions
        C.w1=corr(oddBetas{1}(:),evenBetas{1}(:));
        C.w2=corr(oddBetas{2}(:),evenBetas{2}(:));
        Acr1 = corr(oddBetas{1}(:),oddBetas{2}(:));
        Acr2 = corr(evenBetas{1}(:),evenBetas{2}(:));
        Acr3 = corr(oddBetas{1}(:),evenBetas{2}(:));
        Acr4 = corr(oddBetas{2}(:),evenBetas{1}(:));
        C.acr=mean([Acr1 Acr2 Acr3 Acr4]);
        C.geoMean=sqrt(C.w1*C.w2);
        
        if C.w1<0 || C.w2<0
            C.geoMean=NaN;
        end

        C.acrSess_corr = C.acr./C.geoMean; % correlation across sessions
        C.acrSess_corr(isnan(C.acrSess_corr))=0;
        varargout{1}=C;
        case 'stability_PLOT'
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        sessTr = [1:3]; % how many transitions
        roi=[1:8];
        seqType={'trained','untrained'};
        sessName={'1-2','2-3','3-4'};
        
        vararginoptions(varargin,{'seqTypeCorr','sessTr','roi'});
       
        TT=[];
        for ss=sessTr
            for st=1:2 % both seqTypes
                T=load(fullfile(stabDir,sprintf('StabilityCeiling_%s_%s_sess%d-sess%d.mat',seqType{st},seqTypeCorr,ss,ss+1)));
                T.sessTr=ones(size(T.sn))*ss;
                T.seqType=ones(size(T.sn))*st;
                TT=addstruct(TT,T);
            end
        end
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.line(TT.sessTr,TT.acr,'subset',TT.roi==r,'split',TT.seqType,'leg',seqType,'leglocation','northeast','style',stySeq);
           % drawline(mean(mean([TT.w1(TT.roi==r) TT.w2(TT.roi==r)],1)),'dir','horz');
            title(regname{r});
            if indx==1
                ylabel(sprintf('Correlation of %s pattern',seqTypeCorr));
            else
                ylabel('');
            end
            indx=indx+1;
            plt.match('y');
            set(gca,'XTickLabel',sessName(sessTr));
        end
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.line(TT.sessTr,TT.acrSess_corr,'subset',TT.roi==r,'split',TT.seqType,'leg',seqType,'leglocation','northeast','style',stySeq);
            % drawline(mean(mean([TT.w1(TT.roi==r) TT.w2(TT.roi==r)],1)),'dir','horz');
            title(regname{r});
            if indx==1
                ylabel(sprintf('Correlation of %s pattern - corrected for witSess corr',seqTypeCorr));
            else
                ylabel('');
            end
            indx=indx+1;
            plt.match('y');
            set(gca,'XTickLabel',sessName(sessTr));
        end
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.bar(TT.sessTr,[TT.w1 TT.w2 TT.acr],'subset',TT.roi==r&TT.seqType==1,'leglocation','northeast','leg',{'with1','with2','across'});
            title(regname{r});
            if indx==1
                ylabel('Correlation within / across sessions trained');
            else
                ylabel('');
            end
            indx=indx+1;
        end
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.bar(TT.sessTr,[TT.w1 TT.w2 TT.acr],'subset',TT.roi==r&TT.seqType==2,'leglocation','northeast','leg',{'with1','with2','across'});
            title(regname{r});
            if indx==1
                ylabel('Correlation within / across sessions untrained');
            else
                ylabel('');
            end
            indx=indx+1;
        end
        
    case 'PCM_simulateReliability'
        runEffect  = 'fixed';
        noise=1;
        scale=0.01;
        theta = [0 0 0 0 -200 -200 -200 -200]';
        algorithm = 'NR';
        part=8;
        cond=12;
        %theta = zeros(8,1);
        vararginoptions(varargin,{'runEffect','theta','noise','scale','algorithm'});
        
        numSim=10; % number of simulations
        % model specifications
        D.numVox  = 1000;
        
        % partitions 1-8 for session 1, 9-16 for session 2
        % conditions given different labels for the two sessions
        D.partVec = [kron([1:part]',ones(cond,1)); kron([1+part:part+part]',ones(cond,1))];  % Partitions
        D.condVec = [kron(ones(part,1),[1:cond]');kron(ones(part,1),[1+cond:cond+cond]')];   % Conditions   
        sn = 8; % 8 subjects
        M = pcm_corrModel;
        TT=[];
        
        for i=1:numSim
            [Data] = pcm_generateData_eva(M{end},theta,D,sn,scale,noise);
            T = pcm_fitModels(Data,M,part,cond,runEffect,algorithm); 
            T.numSim = ones(size(T.bayesEst))*numSim;
            TT=addstruct(TT,T);
        end
        keyboard;
    case 'PCM_simulateReliability_seqType'
        runEffect  = 'fixed';
        mNum=2; % which model used for data generation - 1: ind, 2: flex, 3: same
        noise=1;
        scale=[0:0.1:1];
        theta=[0.1 0.2 0.2]';        %[0.1 0.2 0.5] - very high corr [0.1 0.2 0] - 0 corr
        part=8;
        cond=6;
        modelCorr={'ind','flex','perfect'};
        fitModel={'generic','specific'}; 
        thetaScale=[1 1 1 1 1 1]';
        trueModel=1; % generic - 1 param for all seq, or specific - 1 param for each seq
        algorithm = 'NR';
        sessType='between';
        sessScale = [1 1];

        vararginoptions(varargin,{'runEffect','noise','scale','algorithm','mNum','trueModel','theta','sessType','thetaScale','sessScale'});
        
        numSim=1; % number of simulations
        
        % model specifications
        D.numVox  = 1000;
        
        % partitions 1-8 for session 1, 9-16 for session 2
        % conditions given different labels for the two sessions
        D.partVec = [kron([1:part]',ones(cond,1)); kron([1+part:part+part]',ones(cond,1))];  % Partitions
        D.condVec = [kron(ones(part,1),[1:cond]');kron(ones(part,1),[1+cond:cond+cond]')];   % Conditions
        signalPartScale=[ones(48,1);ones(48,1)*5];
        switch(sessType)
            case 'within'
                D.partVec = [kron([1:part]',ones(cond,1)); kron([1:part]',ones(cond,1))];  % Partitions
            case 'between'
                D.partVec = [kron([1:part]',ones(cond,1)); kron([1+part:part+part]',ones(cond,1))];  % Partitions
        end 
        sn = 25; % 25 subjects
        % create both models
        M{1}=pcm_corrModel_seqType_fixed;
        M{2}=pcm_corrModel_seqType_fixed_indSeq;
        
        TT=[];
        
        for s=1:length(scale)
            for i=1:numSim
                switch trueModel
                    case 2
                        theta = [repmat(theta(1),6,1).*thetaScale; repmat(theta(2),6,1).*thetaScale; repmat(theta(3),6,1).*thetaScale];
                end
                Data = pcm_generateData_eva(M{trueModel}{mNum},theta,D,sn,scale(s),noise,'signalPartScale',signalPartScale);
                
                switch(sessType)
                    case 'between' % scale data for each session
                        for ss=1:sn
                            Data{ss}(D.partVec<9,:)=Data{ss}(D.partVec<9,:).*sessScale(1);
                            Data{ss}(D.partVec>8,:)=Data{ss}(D.partVec>8,:).*sessScale(2);
                        end
                end
                for f=1:length(fitModel)
                    C = pcm_correlation(Data,D.partVec,D.condVec,M{f}{2},runEffect,fitModel{f});
                    [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M{f},D.partVec,D.condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);;
                    T.bayesFactor=bsxfun(@minus,T.likelihood,T.likelihood(:,1));
                    C.signalLevel=ones(size(T.SN))*scale(s);
                    C.fitModel=ones(size(T.SN))*f;
                    TT=addstruct(TT,C);
                    TT=addstruct(TT,T);
                end
            end
        end
        
        
        trueCorr  = calcCorr_thetas(theta(1),theta(2),theta(3));
        
        figure
        for fm=1:2
            subplot(1,2,fm)
            scatterplot(TT.signalLevel,TT.bayesFactor(:,2),'subset',TT.fitModel==fm,'markercolor',[0 0 1]);
            hold on
            scatterplot(TT.signalLevel,TT.bayesFactor(:,3),'subset',TT.fitModel==fm,'markercolor',[1 0 0]);
            xlabel('Signal level'); ylabel('Bayes factor');
            title(sprintf('True model %s %s - corr %d - %s session - fit %s model',fitModel{trueModel},modelCorr{mNum},trueCorr,sessType,fitModel{fm}));
            legend({'flex','perfect'});
            drawline(0,'dir','horz');
        end
        
        keyboard;
        figure
        subplot(1,4,1)
        scatterplot(TT.signalLevel,TT.r_naive,'subset',TT.fitModel==trueModel);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(trueCorr,'dir','horz','color',[1 0 0]);
        title('Naive correlation');
        xlabel('Signal level'); ylabel('Corr values');
        
        subplot(1,4,2)
        scatterplot(TT.signalLevel,TT.r_crossval,'subset',TT.fitModel==trueModel);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(trueCorr,'dir','horz','color',[1 0 0]);
        title('Crossval correlation');
        xlabel('Signal level'); 
        
        subplot(1,4,3)
        scatterplot(TT.signalLevel,TT.r_model2,'subset',TT.fitModel==1);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(trueCorr,'dir','horz','color',[1 0 0]);
        title('PCM correlation - generic model');
        xlabel('Signal level');
        
        subplot(1,4,4)
        scatterplot(TT.signalLevel,TT.r_model2,'subset',TT.fitModel==2);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(trueCorr,'dir','horz','color',[1 0 0]);
        title('PCM correlation - specific model');
        xlabel('Signal level');      
    case 'PCM_simulateCorr_seqType'
        runEffect = 'fixed';
        mNum=2; % which model used for data generation - 1: ind, 2: flex, 3: same
        noise=1;
        scale=[0:0.1:1];
        corrTheta=0;
        thetaFix=[0.1 0.2];
        part=8;
        cond=6;
        fitModel={'generic','specific'}; 
        thetaScale=[1 1 1 1 1 1]';
        trueModel=2; % generic - 1 param for all seq, or specific - 1 param for each seq

       % N.noise_corrTheta=[0 1; 0 0.2; 0.5 0.1; 0.5 0.2; 1 0; 1 0.2; 2 0; 2 0.2]; 

        sessType = 'within'; % within or between sessions
        
        vararginoptions(varargin,{'runEffect','thetaFix','noise','scale','mNum','trueModel','corrTheta','sessType','thetaScale'});
        
        numSim=1; % number of simulations
        % model specifications
        D.numVox  = 1000;  
        D.condVec = [kron(ones(part,1),[1:cond]');kron(ones(part,1),[1+cond:cond+cond]')];   % Conditions
        switch(sessType)
            case 'within'
                D.partVec = [kron([1:part]',ones(cond,1)); kron([1:part]',ones(cond,1))];  % Partitions
            case 'between'
                D.partVec = [kron([1:part]',ones(cond,1)); kron([1+part:part+part]',ones(cond,1))];  % Partitions
        end
        sn = 25; % 25 subjects
        
        % create both models
        M1=pcm_corrModel_seqType_fixed;
        M2=pcm_corrModel_seqType_fixed_indSeq;
        M{1}=M1{mNum}; % independent model
        M{2}=M2{mNum}; % all sequences equal model
        
        TT=[];
        SS=[];
        for t=1:length(corrTheta)
            for s=1:length(scale)
                for i=1:numSim
                    switch trueModel
                        case 1
                            theta = [thetaFix(1) thetaFix(2) corrTheta(t)]';
                        case 2
                            theta = [repmat(thetaFix(1),6,1).*thetaScale; repmat(thetaFix(2),6,1).*thetaScale; repmat(corrTheta(t),6,1).*thetaScale];
                    end
                    Data = pcm_generateData_eva(M{trueModel},theta,D,sn,scale(s),noise);
                    
                    % fit both models
                    for f=1:length(fitModel)
                        C = pcm_correlation(Data,D.partVec,D.condVec,M{f},runEffect,fitModel{f});
                        
                        trueCorr  = calcCorr_thetas(thetaFix(1),thetaFix(2),corrTheta(t));
                        C.trueCorr=ones(size(C.r_model2))*trueCorr;
                        C.signalLevel=ones(size(C.r_model2))*scale(s);
                        C.fitModel=ones(size(C.r_model2))*f;
                        
                        % summary evaluation statistics
                        S.trueCorr=[trueCorr; trueCorr; trueCorr];
                        S.bias=[trueCorr-mean(C.r_naive);trueCorr-mean(C.r_crossval);trueCorr-mean(C.r_model2)];
                        S.var=[var(C.r_naive);var(C.r_crossval);var(C.r_model2)];
                        S.mse=[sum((C.trueCorr-C.r_naive).^2);
                            sum((C.trueCorr-C.r_crossval).^2);
                            sum((C.trueCorr-C.r_model2).^2)];
                        S.corrType=[1;2;3]; % 1-naive, 2-crossval, 3-model
                        S.signalLevel=[scale(s);scale(s);scale(s)];
                        S.fitModel=[f;f;f];
                        SS=addstruct(SS,S);
                        TT=addstruct(TT,C);
                    end
                end
            end
        end
        
        figure
        subplot(1,4,1)
        scatterplot(TT.signalLevel,TT.r_naive,'subset',TT.fitModel==trueModel);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title('Naive correlation');
        xlabel('Signal level'); ylabel('Corr values');
        
        subplot(1,4,2)
        scatterplot(TT.signalLevel,TT.r_crossval,'subset',TT.fitModel==trueModel);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title('Crossval correlation');
        xlabel('Signal level'); 
        
        subplot(1,4,3)
        scatterplot(TT.signalLevel,TT.r_model2,'subset',TT.fitModel==1);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title('PCM correlation - generic model');
        xlabel('Signal level');
        
        subplot(1,4,4)
        scatterplot(TT.signalLevel,TT.r_model2,'subset',TT.fitModel==2);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title('PCM correlation - specific model');
        xlabel('Signal level');
        
        
        leg_labels={'naive','crossval','model generic','model specific'};
        
        % summary stats
        A=getrow(SS,SS.fitModel==trueModel | SS.corrType==3);
        A.corrType(A.fitModel==2&A.corrType==3)=4;
        figure
        subplot(1,3,1)
        plt.line(A.signalLevel,A.bias,'split',A.corrType,'leg','off');
        ylabel('Bias'); xlabel('Signal level'); title(sprintf('True correlation %d',A.trueCorr(1)));
        subplot(1,3,2)
        plt.line(A.signalLevel,A.var,'split',A.corrType,'leg','off');
        ylabel('Variance'); xlabel('Signal level');
        subplot(1,3,3)
        plt.line(A.signalLevel,A.mse,'split',A.corrType,'leg',leg_labels,'leglocation','northeast');
        ylabel('MSE'); xlabel('Signal level');

        keyboard;
        
    case 'PCM_extremeSimulation'
        % patterns are identical or completely random
        numPart = 8;
        sn=8;
        numVox=100;
        numCond = 12;
        runEffect='random';
        algorithm='NR';
        typeModel='perfectCorr';    % options - perfectCorr or zeroCorr
        
        vararginoptions(varargin,{'numPart','sn','numVox','numCond','runEffect','algorithm','typeModel'});
        
        
        M = pcm_corrModel;
        
        for s=1:sn
            D{s} = [randn(numCond*numPart,numVox)];
            switch (typeModel)
                case 'perfectCorr'
                    D{s}=[D{s}; D{s}];
                case 'zeroCorr'
                    D{s}=[D{s};randn(numCond*numPart,numVox)];
            end
        end
        
        partVec      = [kron([1:numPart]',ones(numCond,1));kron([numPart+1:numPart+numPart]',ones(numCond,1))];            % Partitions
        condVec      = [kron(ones(numPart,1),[1:numCond]');kron(ones(numPart,1),[numCond+1:numCond+numCond]')];            % Conditions
        %call PCM
        T = pcm_fitModels(D,M,partVec,condVec,runEffect,algorithm);
        keyboard;
    case 'PCM_constructReliability'
        runEffect  = 'random';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1:8];
        sn=[1:9,11,12];
        sessN=[1:2]; % need to be two sessions at the time
        AllReg=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm'})
        
        for r = reg
            for p=1:length(sn)
                partVec{p}=[];
                condVec{p}=[];
                for ss = 1:numel(sessN)
                    
                    B=load(fullfile(regDir,sprintf('betas_sess%d.mat',sessN(ss))));
                    glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                    
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==r)}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                    end
                    
                    indx = ones(size(D.run));
                    if ss == 1
                        condVec{p} = D.seqNumb; % conditions
                        partVec{p} = D.run;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    else
                        condVec{p} = [condVec{p}; D.seqNumb + 12]; % treat 2nd session as additional conditions
                        partVec{p} = [partVec{p}; D.run + 8];  % runs/partitions of 2nd session as additional runs
                        Data{p} = [Data{p}; beta(:,indx==1)'];  % Data is N x P (cond x voxels) - no intercept
                    end;
                end; % session
            end; % subj
            
            % construct models
            %M = pcm_corrModel;
            %M = pcm_corrModel_seqSPEC_only;
            M=pcm_corrModel_seqTypestable_seqSpec;
            T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
            T.roi = ones(size(T.SN))*r;
            AllReg=addstruct(AllReg,T);
        end
        
        % save output
        %save(fullfile(stabDir,sprintf('PCM_reliability_sess%d-sess%d.mat',sessN(1),sessN(2))),'-struct','AllReg');
        %save(fullfile(stabDir,sprintf('PCM_reliability_seqSPEC_only_sess%d-sess%d.mat',sessN(1),sessN(2))),'-struct','AllReg');
        save(fullfile(stabDir,sprintf('PCM_reliability_seqSPECchange_seqTypestable_sess%d-sess%d.mat',sessN(1),sessN(2))),'-struct','AllReg');
    case 'PCM_constructReliability_seqType'
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1:8];
        sn=[1:18];
        sessN=[1:2]; % need to be two sessions at the time
        AllReg=[];
        seqType='trained';
        parcelType='Brodmann'; % Brodmann or 162tessels
        models={'generic','specific'}; 
        fitModel=1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','seqType','runEffect','fitModel','parcelType'})
        
        switch seqType
            case 'trained'
                stIndx=1;
            case 'untrained'
                stIndx=2;
        end
            for r = reg
                for p=1:length(sn)
                    partVec{p}=[];
                    condVec{p}=[];
                    for ss = 1:numel(sessN)
                        
                        B=load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
                        glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                        
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p)&B.region==r)}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                        end
                        
                        indx = ones(size(D.run));
                        if ss == 1
                            condVec{p} = D.seqNumb(D.seqType==stIndx); % conditions
                            partVec{p} = D.run(D.seqType==stIndx);
                            Data{p} = beta(:,indx==1&D.seqType==stIndx)';  % Data is N x P (cond x voxels) - no intercept
                            if stIndx==2
                                condVec{p}=condVec{p}-6;    % make conditions always from 1
                            end
                        else
                            if stIndx==1
                                condVec{p} = [condVec{p}; D.seqNumb(D.seqType==stIndx) + 6]; % treat 2nd session as additional conditions
                            elseif stIndx==2
                                condVec{p} = [condVec{p}; D.seqNumb(D.seqType==stIndx)]; % for untrained it already starts with 6
                            end
                            partVec{p} = [partVec{p}; D.run(D.seqType==stIndx) + 8];  % runs/partitions of 2nd session as additional runs
                            Data{p} = [Data{p}; beta(:,indx==1 & D.seqType==stIndx)'];  % Data is N x P (cond x voxels) - no intercept
                        end;
                    end; % session
                end; % subj
                
                % construct models
                %M = pcm_corrModel_seqType;
                switch(fitModel) % both now assume run effect = fixed
                    case 1 % generic
                        % M = pcm_corrModel_seqType_fixed;
                        M = pcm_stabilityModel_generic;
                    case 2 % specific
                        %M = pcm_corrModel_seqType_fixed_indSeq;
                        M = pcm_stabilityModel_specific;
                end
                T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                C = pcm_correlation(Data,partVec{1},condVec{1},M{4},runEffect,models{fitModel});
                %C = pcm_correlation(Data,partVec{1},condVec{1},M{2},runEffect,models{fitModel});
                T.roi = ones(size(T.SN))*r;
                AllReg=addstruct(AllReg,T);
                AllReg=addstruct(AllReg,C);
            end
            
            %remove some fields
            a3='theta_hat'; a4='thetaCr';
            AllReg=rmfield(AllReg,a3); AllReg=rmfield(AllReg,a4);
            % save output
            save(fullfile(stabDir,sprintf('PCM_reliability_NEW_%s_%s_%s_sess%d-sess%d.mat',models{fitModel},runEffect,seqType,sessN(1),sessN(2))),'-struct','AllReg');
    case 'PCM_constructReliability_hyperModel'
        runEffect  = 'random';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1:8];
        sn=[1:9,11,12];
        sessN=[1:2]; % need to be two sessions at the time
        AllReg=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm'})
        
        for r = reg
            for p=1:length(sn)
                partVec{p}=[];
                condVec{p}=[];
                for ss = 1:numel(sessN)
                    
                    B=load(fullfile(regDir,sprintf('betas_sess%d.mat',sessN(ss))));
                    glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                    
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==r)}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                    end
                    
                    indx = ones(size(D.run));
                    if ss == 1
                        condVec{p} = D.seqNumb; % conditions
                        partVec{p} = D.run;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    else
                        condVec{p} = [condVec{p}; D.seqNumb + 12]; % treat 2nd session as additional conditions
                        partVec{p} = [partVec{p}; D.run + 8];  % runs/partitions of 2nd session as additional runs
                        Data{p} = [Data{p}; beta(:,indx==1)'];  % Data is N x P (cond x voxels) - no intercept
                    end;
                end; % session
            end; % subj
            
            % construct models
            M = pcm_corrHyperModel;
            T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
            T.roi = ones(size(T.SN))*r;
            AllReg=addstruct(AllReg,T);
        end
        
        % save output
        save(fullfile(stabDir,sprintf('PCM_reliability_hyperModel_sess%d-sess%d.mat',sessN(1),sessN(2))),'-struct','AllReg');
       
    case 'PCM_plot_allSess'
        reg = [1:8];
        sessName={'Sess 1-2','2-3','3-4'};
        vararginoptions(varargin,{'reg','sessN'});
        
        R1=load(fullfile(stabDir,'PCM_reliability_sess1-sess2.mat'));
        %R1=load(fullfile(stabDir,'PCM_reliability_seqSPEC_sess1-sess2.mat'));
        R1.sessTr=ones(size(R1.SN));
        R2=load(fullfile(stabDir,'PCM_reliability_sess2-sess3.mat'));
        %R2=load(fullfile(stabDir,'PCM_reliability_seqSPEC_sess2-sess3.mat'));
        R2.sessTr=ones(size(R2.SN))*2;
        R3=load(fullfile(stabDir,'PCM_reliability_sess3-sess4.mat'));
        %R3=load(fullfile(stabDir,'PCM_reliability_seqSPEC_sess3-sess4.mat'));
        R3.sessTr=ones(size(R3.SN))*3;
        
        T=[];
        T=addstruct(T,R1); T=addstruct(T,R2); T=addstruct(T,R3);
        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3];
        T2.bayesEst=T.bayesEst(:);
        T2.roi=[T.roi;T.roi;T.roi];
        T2.sessTr=[T.sessTr;T.sessTr;T.sessTr];
        % rearranging how the data structure is arranged
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line(T2.sessTr,T2.bayesEst,'subset',T2.modelInd>1&T2.roi==r,'split',T2.modelInd,'leg',{'flex','perfect'},'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel('Log-Bayes');
            else
                ylabel('');
            end
        end
    case 'PCM_plot_seqType_allSess'
        reg = [1:8];
        sessName={'1-2','2-3','3-4'};
        seqType='trained';
        runEffect='fixed';
        modelType='specific'; % generic or specific
        numTrans=3; % number of session transitions - 2:sess1-2-3; 3:sess1-2-3-4
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','numTrans','modelType'});

        T=[];
        for tr=1:numTrans
            R=load(fullfile(stabDir,sprintf('PCM_reliability_%s_%s_%s_sess%d-sess%d.mat',modelType,runEffect,seqType,tr,tr+1)));
            R.sessTr=ones(size(R.SN))*tr;
            T=addstruct(T,R);
        end
        
        T2.bayesEst=T.bayesEst(:);
        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3]; 
        T2.roi=[T.roi;T.roi;T.roi]; 
        T2.sessTr=[T.sessTr;T.sessTr;T.sessTr];
        T2.cross_likelihood=T.cross_likelihood(:);

        % rearranging how the data structure is arranged
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line([T2.sessTr>2 T2.sessTr],T2.bayesEst,'subset',T2.modelInd>1&T2.roi==r,'split',T2.modelInd,'leg',{'flex','perfect'},'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel(sprintf('Log-Bayes %s - %s model',seqType,modelType));
                xlabel('Sess transitions');
            else
                ylabel('');
            end
        end
    case 'PCM_plot_logBayes_NEW'
         reg = [1:8];
         sessName={'1-2','2-3','3-4'};
         seqType='trained';
         runEffect='fixed';
         modelType='generic'; % or specific
         numTrans=3; % number of session transitions - 2:sess1-2-3; 3:sess1-2-3-4
         vararginoptions(varargin,{'reg','sessN','seqType','runEffect','numTrans','modelType'});

        T=[];
        for tr=1:numTrans
            R=load(fullfile(stabDir,sprintf('PCM_reliability_NEW_%s_%s_%s_sess%d-sess%d.mat',modelType,runEffect,seqType,tr,tr+1)));
            R.sessTr=ones(size(R.SN))*tr;
            T=addstruct(T,R);
        end
        
        r='reg';
        T=rmfield(T,r);
        TT=T;
        TT2.bayesEst(:,[2:5])=bsxfun(@minus,TT.bayesEst(:,[2:5]),TT.bayesEst(:,2));
        TT2.modelInd=[ones(size(TT.SN));ones(size(TT.SN))*2;ones(size(TT.SN))*3;ones(size(TT.SN))*4;ones(size(TT.SN))*5];
        TT2.bayesEst=TT.bayesEst(:);
        TT2.roi=[TT.roi;TT.roi;TT.roi;TT.roi;TT.roi];
        TT2.sessTr=[TT.sessTr;TT.sessTr;TT.sessTr;TT.sessTr;TT.sessTr];
        

        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3;ones(size(T.SN))*4;ones(size(T.SN))*5];
        T2.bayesEst=T.bayesEst(:);
        T2.roi=[T.roi;T.roi;T.roi;T.roi;T.roi];
        T2.sessTr=[T.sessTr;T.sessTr;T.sessTr;T.sessTr;T.sessTr];
        % rearranging how the data structure is arranged
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line([T2.sessTr>2 T2.sessTr],T2.bayesEst,'subset',T2.modelInd>1&T2.roi==r,'split',T2.modelInd,'leg',{'session','session+seq','session+seq+corr','session+seq+perfectCorr'},'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            ylim([-5 170]);
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel(sprintf('Log-Bayes %s',seqType));
            else
                ylabel('');
            end
        end
    case 'PCM_plot_logBayes_seqType'
         reg = [1:8];
         sessName={'1-2','2-3','3-4'};
         seqType={'trained','untrained'};
         runEffect='fixed';
         modelType='generic'; % or specific
         numTrans=3; % number of session transitions - 2:sess1-2-3; 3:sess1-2-3-4
         vararginoptions(varargin,{'reg','sessN','runEffect','numTrans','modelType'});

        T=[];
        for st=1:2
            for tr=1:numTrans
                R=load(fullfile(stabDir,sprintf('PCM_reliability_NEW_%s_%s_%s_sess%d-sess%d.mat',modelType,runEffect,seqType{st},tr,tr+1)));
                R.sessTr=ones(size(R.SN))*tr;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        r='reg';
        T=rmfield(T,r);
        TT=T;
        TT2.bayesEst(:,[2:5])=bsxfun(@minus,TT.bayesEst(:,[2:5]),TT.bayesEst(:,2));
        TT2.modelInd=[ones(size(TT.SN));ones(size(TT.SN))*2;ones(size(TT.SN))*3;ones(size(TT.SN))*4;ones(size(TT.SN))*5];
        TT2.bayesEst=TT.bayesEst(:);
        TT2.roi=[TT.roi;TT.roi;TT.roi;TT.roi;TT.roi];
        TT2.sessTr=[TT.sessTr;TT.sessTr;TT.sessTr;TT.sessTr;TT.sessTr];
        TT2.seqType=[TT.seqType;TT.seqType;TT.seqType;TT.seqType;TT.seqType];
        

        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3;ones(size(T.SN))*4;ones(size(T.SN))*5];
        T2.bayesEst=T.bayesEst(:);
        T2.roi=[T.roi;T.roi;T.roi;T.roi;T.roi];
        T2.sessTr=[T.sessTr;T.sessTr;T.sessTr;T.sessTr;T.sessTr];
        T2.seqType=[T.seqType;T.seqType;T.seqType;T.seqType;T.seqType];
        % rearranging how the data structure is arranged
        figure
        for r=reg
            subplot(1,numel(reg),r)
            % plot trained vs. untrained for the flexible model
            plt.line([T2.sessTr>2 T2.sessTr],T2.bayesEst,'subset',T2.modelInd==5&T2.roi==r,'split',T2.seqType,'leg',seqType,'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            ylim([-5 170]);
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Log-Bayes');
            else
                ylabel('');
            end
        end
    case 'PCM_plot_seqType_corr'
        reg = [1:8];
        sessName={'1-2','2-3','3-4'};
        seqType={'trained','untrained'};
        runEffect='fixed';
        modelType='specific'; % generic or specific
        Trans=[1 2]; % which transitions: 1: 1-2, 2: 2-3, 3: 3-4
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','Trans','modelType'});
        
        T=[];
        for st=1:length(seqType)
            for tr=Trans
                R=load(fullfile(stabDir,sprintf('PCM_reliability_NEW_%s_%s_%s_sess%d-sess%d.mat',modelType,runEffect,seqType{st},tr,tr+1)));
                %R=load(fullfile(stabDir,sprintf('PCM_reliability_%s_%s_%s_sess%d-sess%d.mat',modelType,runEffect,seqType{st},tr,tr+1)));
                R.sessTr=ones(size(R.SN))*tr;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        indx=1;
        for r=reg
            figure(1)
            subplot(1,numel(reg),indx)
            plt.line([T.sessTr>2 T.sessTr],T.r_naive,'subset',T.roi==r,'split',T.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if indx==1
                ylabel('Naive correlation'); xlabel('Sess transitions');
            else
                ylabel('');
            end
            
            figure(2)
            subplot(1,numel(reg),indx)
            plt.line([T.sessTr>2 T.sessTr],T.r_crossval,'subset',T.roi==r,'split',T.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel('Crossval correlation - seqSpec');xlabel('Sess transitions');
            else
                ylabel('');
            end
            
            figure(3)
            subplot(1,numel(reg),indx)
            plt.line([T.sessTr>2 T.sessTr],T.r_crossval_seqType,'subset',T.roi==r,'split',T.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel('Crossval correlation - seqType');xlabel('Sess transitions');
            else
                ylabel('');
            end
            
            figure(4)
            subplot(1,numel(reg),indx)
            plt.line([T.sessTr>2 T.sessTr],T.r_model2,'subset',T.roi==r,'split',T.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel(sprintf('Model correlation %s',modelType));xlabel('Sess transitions');
            else
                ylabel('');
            end
            
            % comparing sessions 3-4
            figure(5)
            subplot(1,numel(reg),indx)
            plt.bar([T.sessTr],T.r_crossval,'subset',T.roi==r,'split',T.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel(sprintf('Crossval seqSpec correlation %s',modelType));xlabel('Sess transitions');
            else
                ylabel('');
            end
            
            figure(6)
            subplot(1,numel(reg),indx)
            plt.bar([T.sessTr],T.r_crossval_seqType,'subset',T.roi==r,'split',T.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel(sprintf('Crossval seqType correlation %s',modelType));xlabel('Sess transitions');
            else
                ylabel('');
            end
            
            indx=indx+1;
        end
    case 'PCM_plot_hyperModel'
        % choose which models you would like to plot
        % trained - untrained
        modelName={'ind-ind','flex-flex','one-one','ind-flex','ind-one','flex-ind','flex-one','one-ind','one-flex'};
        reg = [1:8];
        sessName={'Sess 1-2','2-3','3-4'};
        vararginoptions(varargin,{'reg','sessN','modelType'});
        
        R1=load(fullfile(stabDir,'PCM_reliability_hyperModel_sess1-sess2.mat'));
        R1.sessTr=ones(size(R1.SN));
        R2=load(fullfile(stabDir,'PCM_reliability_hyperModel_sess2-sess3.mat'));
        R2.sessTr=ones(size(R2.SN))*2;
        R3=load(fullfile(stabDir,'PCM_reliability_hyperModel_sess3-sess4.mat'));
        R3.sessTr=ones(size(R3.SN))*3;
        
                
        T=[];
        T=addstruct(T,R1); T=addstruct(T,R2); T=addstruct(T,R3);
        
        numModel=size(modelType,2); % number of models
        
        mIndx=[];
        for n=1:numModel
            mIndx=[mIndx;ones(size(T.SN))*n];
        end
          
        T2.modelInd=mIndx;
        bayesEst=T.bayesEst(:,modelType);
        T2.bayesEst=bayesEst(:);
        T2.roi=repmat(T.roi,numModel,1);
        T2.sessTr=repmat(T.sessTr,numModel,1);
        % rearranging how the data structure is arranged
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line(T2.sessTr,T2.bayesEst,'subset',T2.roi==r,'split',T2.modelInd,'leg',modelName(modelType),'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            set(gca,'XTickLabel',sessName);
            if r==1
                ylabel('Log-Bayes');
            else
                ylabel('');
            end
        end
    case 'PCM_plot_speedCorr'
        reg = [1:5,7];
        regExcl=[6,8]; % V1/2
        seqType={'trained','untrained'};
        runEffect='fixed';
        modelType='generic'; % generic or specific
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType','regExcl'});
        
        T=[];
        for st=1:length(seqType)    
            R=load(fullfile(stabDir,sprintf('PCM_reliability_NEW_%s_%s_%s_sess3-sess4.mat',modelType,runEffect,seqType{st})));
            R.seqType=ones(size(R.SN))*st;
            T=addstruct(T,R);    
        end
        
        T=getrow(T,~ismember(T.roi,regExcl));
        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        figure
        plt.line(T.reg,T.r_model2,'split',T.seqType,'style',stySeq);
        
    case 'PCM_constructSimpleModel_oneSess'
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1:8];
        sn=[1:9,11,12];
        sessN=3; % one session
        AllReg=[];
        seqType='trained';
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','seqType','runEffect'})
        
        switch seqType
            case 'trained'
                stIndx=1;
            case 'untrained'
                stIndx=2;
        end
        
        for r = reg
            for p=1:length(sn)
                partVec{p}=[];
                condVec{p}=[];
                for ss = sessN
                    
                    B=load(fullfile(regDir,sprintf('betas_sess%d.mat',ss)));
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==r)}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                    end
                    
                    indx = ones(size(D.run));
                    condVec{p} = D.seqNumb(D.seqType==stIndx); % conditions
                    partVec{p} = D.run(D.seqType==stIndx);
                    Data{p} = beta(:,indx==1&D.seqType==stIndx)';  % Data is N x P (cond x voxels) - no intercept
                    if stIndx==2
                        condVec{p}=condVec{p}-6;    % make conditions always from 1
                    end
                end; % session
            end; % subj
            
            % construct models
            M = pcm_simpleModel_oneSess;
            T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
            T.roi = ones(size(T.SN))*r;
            AllReg=addstruct(AllReg,T);
        end
         % save output
        save(fullfile(stabDir,sprintf('PCM_simpleModel_%s_sess%d.mat',seqType,sessN)),'-struct','AllReg');
    case 'PCM_plot_simpleModel'
        seqType={'trained','untrained'};
        sessN=[1:4];
        roi=[1:8];
        vararginoptions(varargin,{'sessN','seqType','roi'});
        
        T=[];
        for s=1:numel(seqType)
            for ss=sessN
                R=load(fullfile(stabDir,sprintf('PCM_simpleModel_%s_sess%d.mat',seqType{s},ss)));
                R.modelBayes=R.bayesEst(:,2);
                R.sessN=ones(size(R.modelBayes))*ss;
                R.seqType=ones(size(R.modelBayes))*s;
                
                T=addstruct(T,R);
            end
        end
       
        figure
        for r=roi
            subplot(1,numel(roi),r)
            plt.line([T.sessN>3 T.sessN],T.modelBayes,'split',T.seqType,'subset',T.roi==r,'style',stySeq,'leg',{'trained','untrained'});
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            if r==1
                xlabel('Session'); ylabel('Bayes factor for sequence model');
            else
                ylabel('');
            end
            fprintf('\nANOVA: session x seqType for %s \n',regname{r});
            anovaMixed(T.modelBayes,T.SN,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r);
            for ss=sessN
                % post-hoc t-test on the effect of sequence type for each session
                fprintf('\npost-hoc t-test on the effect of seqType in session %d in %s \n',ss,regname{r});
                ttestDirect(T.modelBayes,[T.seqType T.SN],2,'paired','subset',T.roi==r & T.sessN==ss);
            end
        end
        
        figure
        plt.line([T.sessN>3 T.sessN],T.modelBayes,'split',T.seqType,'style',stySeq,'leg',{'trained','untrained'});
        drawline(0,'dir','horz');
        xlabel('Session'); ylabel('Bayes factor for sequence model');
        title('All regions together');
        
        fprintf('\nANOVA: session x seqType for all regions \n');
        anovaMixed(T.modelBayes,T.SN,'within',[T.sessN T.seqType],{'session','seqType'});
        for ss=sessN
            fprintf('\npost-hoc t-test on the effect of seqType in session %d in all regions \n',ss);
            ttestDirect(T.modelBayes,[T.seqType T.SN],2,'paired','subset',T.sessN==ss);
        end
    case 'PCM_simpleModel_thetas'
        
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1:8];
        sn=[1:9,11,12];
        sessN=3; % one session
        AllReg=[];
        seqType='trained';
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','seqType','runEffect'})
        
        switch seqType
            case 'trained'
                stIndx=1;
            case 'untrained'
                stIndx=2;
        end
        
        for r = reg
            for p=1:length(sn)
                partVec{p}=[];
                condVec{p}=[];
                for ss = sessN
                    
                    B=load(fullfile(regDir,sprintf('betas_sess%d.mat',ss)));
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==r)}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                    end
                    
                    indx = ones(size(D.run));
                    condVec{p} = D.seqNumb(D.seqType==stIndx); % conditions
                    partVec{p} = D.run(D.seqType==stIndx);
                    Data{p} = beta(:,indx==1&D.seqType==stIndx)';  % Data is N x P (cond x voxels) - no intercept
                    if stIndx==2
                        condVec{p}=condVec{p}-6;    % make conditions always from 1
                    end
                end; % session
            end; % subj
            
            % construct models
            M = pcm_simpleModel_oneSess;
            [D,theta,G_hat] = pcm_fitModelIndivid(Data,M{2},partVec,condVec,'runEffect',runEffect);
            D.theta=theta{1}';
            D.roi = ones(size(D.SN))*r;
            AllReg=addstruct(AllReg,D);
        end
         % save output
        save(fullfile(stabDir,sprintf('PCM_simpleModel_thetas_%s_sess%d.mat',seqType,sessN)),'-struct','AllReg');

    case 'PCM_plot_thetas'
        seqType={'trained','untrained'};
        sessN=[1:4];
        roi=[1:8];
        vararginoptions(varargin,{'sessN','seqType','roi'});
        
        TT=[];
        for s=1:numel(seqType)
            for ss=sessN
                    R=load(fullfile(stabDir,sprintf('PCM_simpleModel_thetas_%s_sess%d.mat',seqType{s},ss)));
                    T.theta=[R.theta(:,1);R.theta(:,2)];
                    T.thetaInd=[ones(size(R.theta,1),1);ones(size(R.theta,1),1)*2]; % 1 - noise, 2 - signal
                    T.roi=[R.roi;R.roi];
                    T.sessN=ones(size(T.theta))*ss;
                    T.seqType=ones(size(T.theta))*s;
                    TT=addstruct(TT,T);
            end
        end
       
        keyboard;
        
        for r=roi
            figure(1)
            subplot(1,numel(roi),r)
            plt.line(TT.sessN,exp(TT.theta),'split',TT.seqType,'subset',TT.thetaInd==1&TT.roi==r,'leg',{'trained','untrained'},'style',stySeq,'leglocation','north');
            title(sprintf('theta %s',regname{r}));
            if r==1
                xlabel('Session'); ylabel('Thetas noise');
            else
                ylabel('');
            end
            
            figure(2)
            subplot(1,numel(roi),r)
            plt.line(TT.sessN,exp(TT.theta),'split',TT.seqType,'subset',TT.thetaInd==2&TT.roi==r,'leg',{'trained','untrained'},'style',stySeq,'leglocation','north');
            title(sprintf('theta %s',regname{r}));
            if r==1
                xlabel('Session'); ylabel('Thetas signal');
            else
                ylabel('');
            end
            
            N=getrow(TT,TT.thetaInd==1);
            S=getrow(TT,TT.thetaInd==2);
            
            figure(3)
            subplot(1,numel(roi),r)
            plt.line(N.sessN,(exp(S.theta))./(exp(N.theta.^2)),'split',N.seqType,'subset',N.roi==r,'leg',{'trained','untrained'},'style',stySeq,'leglocation','north');
            title(sprintf('theta %s',regname{r}));
            if r==1
                xlabel('Session'); ylabel('Thetas signal/noise');
            else
                ylabel('');
            end

        end
        
    case 'STATS_PCM_corr'
        reg = [1:8];
        seqType={'trained','untrained'};
        runEffect='fixed';
        modelType='specific'; % generic or specific
        Trans=[1 2]; % which transitions: 1: 1-2, 2: 2-3, 3: 3-4
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','Trans','modelType'});
        
        T=[];
        for st=1:length(seqType)
            for tr=Trans
                R=load(fullfile(stabDir,sprintf('PCM_reliability_NEW_%s_%s_%s_sess%d-sess%d.mat',modelType,runEffect,seqType{st},tr,tr+1)));
                R.sessTr=ones(size(R.SN))*tr;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        for r=reg
            % naive correlation
            fprintf('\nNaive corr: sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.r_naive,T.SN,'within',[T.sessTr T.seqType],{'sessTrans','seqType'},'subset',T.roi==r);
            
            % crossval seq-specific
            fprintf('\nCrossval seqSPEC corr: sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.r_crossval,T.SN,'within',[T.sessTr T.seqType],{'sessTrans','seqType'},'subset',T.roi==r);
            
            % crossval seqType
           fprintf('\nCrossval seqType corr: sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.r_crossval_seqType,T.SN,'within',[T.sessTr T.seqType],{'sessTrans','seqType'},'subset',T.roi==r);
        end
        
        
         % Post-hoc t-tests for effect of seqType per sessTransition
        for r=reg
            for ss=Trans
                fprintf('\n post-hoc t-test on the crossval seqSPEC - effect of seqType in sessTr %d in %s \n',ss,regname{r});
                ttestDirect(T.r_crossval,[T.seqType T.SN],2,'paired','subset',T.roi==r&T.sessTr==ss);
                fprintf('\n post-hoc t-test on the crossval seqTYPE effect of seqType in sessTr %d in %s \n',ss,regname{r});
                ttestDirect(T.r_crossval_seqType,[T.seqType T.SN],2,'paired','subset',T.roi==r&T.sessTr==ss);
            end
        end
        
         % Post-hoc t-tests for effect of sessTransition per seqType
        for r=reg
            for ss=1:numel(seqType)
                fprintf('\n post-hoc t-test on the crossval seqSPEC - effect of sessTr in seqType %s in %s \n',seqType{ss},regname{r});
                ttestDirect(T.r_crossval,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.seqType==ss);
                fprintf('\n post-hoc t-test on the crossval seqTYPE - effect of sessTr in seqType %s in %s \n',seqType{ss},regname{r});
                ttestDirect(T.r_crossval_seqType,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.seqType==ss);
            end
        end       
    case 'STATS_PCM_corr_speed'
        reg = [1:8];
        seqType={'trained','untrained'};
        runEffect='fixed';
        modelType='generic'; % generic or specific
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType'});
        
        T=[];
        for st=1:length(seqType)    
            R=load(fullfile(stabDir,sprintf('PCM_reliability_NEW_%s_%s_%s_sess3-sess4.mat',modelType,runEffect,seqType{st})));
            R.seqType=ones(size(R.SN))*st;
            T=addstruct(T,R);    
        end
        for r=reg
            % naive correlation
            fprintf('\nNaive corr: sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.r_naive,T.SN,'within',[T.seqType],{'seqType'},'subset',T.roi==r);
            
            % crossval seq-specific
            fprintf('\nCrossval seqSPEC corr: sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.r_crossval,T.SN,'within',[T.seqType],{'seqType'},'subset',T.roi==r);
            
            % crossval seqType
           fprintf('\nCrossval seqType corr: sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.r_crossval_seqType,T.SN,'within',[T.seqType],{'seqType'},'subset',T.roi==r);
        end
        
        
         % Post-hoc t-tests for effect of seqType per sessTransition
        for r=reg
                fprintf('\n post-hoc t-test on the crossval seqSPEC - effect of seqType in %s \n',regname{r});
                ttestDirect(T.r_crossval,[T.seqType T.SN],2,'paired','subset',T.roi==r);
                fprintf('\n post-hoc t-test on the crossval seqTYPE effect of seqType in %s \n',regname{r});
                ttestDirect(T.r_crossval_seqType,[T.seqType T.SN],2,'paired','subset',T.roi==r);
                fprintf('\n post-hoc t-test on the model corr effect of seqType in %s \n',regname{r});
                ttestDirect(T.r_model2,[T.seqType T.SN],2,'paired','subset',T.roi==r);
        end
    case 'CORR_splitHalf'
        sn=[1:9,11:25];
        sessN=[1:4];
        betaChoice='multiPW';
        vararginoptions(varargin,{'sn','roi','sessN',betaChoice});
        
        AllCorr=[];
        % Make the components for estimation of loadings in G
        Gc{1} = [ones(6,6)];   % seqType
        Gc{2} = [eye(6)];      % sequence
        GX=[];
        for i=1:2
            GX=[GX Gc{i}(:)];
        end;
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_FoSEx_sess%d',ss)));
            for s=sn
                D = load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM_info'));
                for r=1:max(T.region)
                    for i=1:2 %repetition
                        for j=i:2
                            for st=1:2 %seqType
                                D1 = getrow(T,T.SN==s & T.region==r);
                                % get betas
                                switch(betaChoice)
                                    case 'multiPW'
                                        betas = D1.betaW{:};
                                    case 'uniPW'
                                        betas = D1.betaUW{:};
                                    case 'raw'
                                        betas = D1.betaRAW{:};
                                end
                                % split into first (i), second (j), seqtype st
                                betas1 = betas(D.FoSEx==i & D.seqType==st,:);
                                betas2 = betas(D.FoSEx==j & D.seqType==st,:);
                                % info
                                D1 = getrow(D,D.FoSEx==i&D.seqType==st);
                                D2 = getrow(D,D.FoSEx==j&D.seqType==st);
                                if st==2
                                    D1.seqNumb = D1.seqNumb-6;
                                    D2.seqNumb = D2.seqNumb-6;
                                end
                                if (i==j)       % Same exe - split-half in the exe
                                    condVec = D1.seqNumb;
                                    partVec = D1.run;
                                    [K,G] = splitHalfCorr(betas1,partVec,condVec,'withinSes');
                                else
                                    % corr across exe
                                    % define condVec, partVec
                                    partVec = [{D1.run} {D2.run}];
                                    condVec = [{D1.seqNumb} {D2.seqNumb}];
                                    data    = [{betas1} {betas2}];
                                    [K,G] = splitHalfCorr(data,partVec,condVec,'acrossSes');
                                end;
                                K.sn        = s;
                                K.sessN     = ss;
                                K.regNum    = r;
                                if r<9
                                    K.regType   = r;
                                    K.regSide   = 1;
                                else
                                    K.regType   = r-8;
                                    K.regSide   = 1;
                                end
                                K.exe1      = i;
                                K.exe2      = j;
                                K.seqType   = st;
                                K.G          = {G};
                                param        = (pinv(GX)*G(:)); % get the two parameters from regression on G
                                K.hcov       = param(1); % overall mean
                                K.dcov       = param(2); % digit-specific variance
                                AllCorr=addstruct(AllCorr,K);  
                            end; % seqtype
                        end; % exe j
                    end; % exe i
                end; % region
                fprintf('%d \t done sess%d: %s\n',s,ss,subj_name{s});
            end;
        end
        %  save the structure
        save(fullfile(repSupDir,'corr_splitHalf_exe'),'-struct','AllCorr');
        
        

    otherwise
        disp('there is no such case.')
end;    % switch(what)
end


%  % Local functions

function dircheck(dir)
% Checks existance of specified directory. Makes it if it does not exist.

if ~exist(dir,'dir');
    %warning('%s didn''t exist, so this directory was created.\n',dir);
    mkdir(dir);
end
end
function A = analytic_cov(th)

    A.var_tr_s1r1    = th.a^2 + th.b1^2 + th.c1^2 + th.d1^2;
    A.var_untr_s1r1  = th.a^2 + th.b2^2 + th.c1^2 + th.d1^2;
    A.var_tr_s1r2    = th.a^2 + th.b1^2 + th.c1^2 + th.d2^2;
    A.var_untr_s1r2  = th.a^2 + th.b2^2 + th.c1^2 + th.d2^2;
    A.var_tr_s2r1    = th.a^2 + th.b1^2 + th.c2^2 + th.d3^2;
    A.var_untr_s2r1  = th.a^2 + th.b2^2 + th.c2^2 + th.d3^2;
    A.var_tr_s2r2    = th.a^2 + th.b1^2 + th.c2^2 + th.d4^2;
    A.var_untr_s2r2  = th.a^2 + th.b2^2 + th.c2^2 + th.d4^2;
   
    A.seqType = [1;2];
    A.cov_with1(1,:) = th.a^2 + th.b1^2 + th.c1^2;
    A.cov_with1(2,:) = th.a^2 + th.b2^2 + th.c1^2;
    A.cov_with2(1,:) = th.a^2 + th.b1^2 + th.c2^2;
    A.cov_with2(2,:) = th.a^2 + th.b2^2 + th.c2^2;
    A.corr_with1(1,:) = A.cov_with1(1,:)/sqrt(A.var_tr_s1r1*A.var_tr_s1r2);
    A.corr_with1(2,:) = A.cov_with1(1,:)/sqrt(A.var_untr_s1r1*A.var_untr_s1r2);
    A.corr_with2(1,:) = A.cov_with2(1,:)/sqrt(A.var_tr_s2r1*A.var_tr_s2r2);
    A.corr_with2(2,:) = A.cov_with2(1,:)/sqrt(A.var_untr_s2r1*A.var_untr_s2r2);
    
    A.cov_acr(1,:) = th.a^2 + th.b1^2;
    A.cov_acr(2,:) = th.a^2 + th.b2^2;
    corr_train_acr = [A.cov_acr(1,:)/sqrt(A.var_tr_s1r1*A.var_tr_s2r1),...
                       A.cov_acr(1,:)/sqrt(A.var_tr_s1r1*A.var_tr_s2r2),...
                       A.cov_acr(1,:)/sqrt(A.var_tr_s1r2*A.var_tr_s2r1),...
                       A.cov_acr(1,:)/sqrt(A.var_tr_s1r2*A.var_tr_s2r2)];
    corr_untrain_acr = [A.cov_acr(1,:)/sqrt(A.var_untr_s1r1*A.var_untr_s2r1),...
                         A.cov_acr(1,:)/sqrt(A.var_untr_s1r1*A.var_untr_s2r2),...
                         A.cov_acr(1,:)/sqrt(A.var_untr_s1r2*A.var_untr_s2r1),...
                         A.cov_acr(1,:)/sqrt(A.var_untr_s1r2*A.var_untr_s2r2)];
    A.corr_acr(1,:) = mean(corr_train_acr);
    A.corr_acr(2,:) = mean(corr_untrain_acr);
    
end
function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Crossvalidated model comparision:
    [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Group fit with %s algorithm done.\n',algorithm);
 
    [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);

    T.cross_likelihood = Tcross.likelihood;
    T.bayesEst = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
    T.theta_hat=theta_hat;
    T.thetaCr = thetaCr;
   % T.modelNum = [1:length(T.cross_likelihood)];
end
function M = pcm_corrModel
% --------------------------------------
% Models of sequence specific changes across sessions
% same for trained + untrained (modelled together)
% seqType always changes in the same way as seq-spec patterns

        % Model1: Model with independent trained / untrained patterns
        % across sessions - 0 correlation
        
        % Model 1 - independent sequence patterns between sessions
        M{1}.type       = 'feature';
        M{1}.numGparams = 10;
        A=zeros(6);
        for i=1:6
            A(i,i)=1;
        end;
        M{1}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns      (theta_a)
        M{1}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   (theta_b)
        M{1}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns     (theta_c)
        M{1}.Ac(:,19:24,4) = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 pattterns   (theta_d)
   
        M{1}.Ac(:,25,5)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    (theta_e)
        M{1}.Ac(:,26,6)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  (theta_f)
        M{1}.Ac(:,27,7)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2    (theta_i)
        M{1}.Ac(:,28,8)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2  (theta_j)
        M{1}.Ac(:,29,9)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{1}.Ac(:,30,10)  = [zeros(12,1);ones(12,1)];  % Session 2

        M{1}.name       = 'ind-ind';
        
        % --------------------------------------
        % Model2: Model with a flexible across-sess correlation for sequences 
        M{2}.type       = 'feature';
        M{2}.numGparams = 14;
        
        M{2}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns      (theta_a)
        M{2}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   (theta_b)
        M{2}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns     (theta_c)
        M{2}.Ac(:,19:24,4) = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 pattterns   (theta_d)
        M{2}.Ac(:,1:6,5)   = [zeros(6);zeros(6);A;zeros(6)];     % Same trained1 patterns        (theta_e)
        M{2}.Ac(:,7:12,6)  = [zeros(6);zeros(6);zeros(6);A];     % Same untrained1 pattterns     (theta_f)
        
        M{2}.Ac(:,25,7)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1        (theta_g)
        M{2}.Ac(:,26,8)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1      (theta_h)
        M{2}.Ac(:,27,9)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2        (theta_i)
        M{2}.Ac(:,28,10) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2      (theta_j)
        M{2}.Ac(:,25,11) = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall SAME component trained1   (theta_k)
        M{2}.Ac(:,26,12) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall SAME component untrained1 (theta_l)
        M{2}.Ac(:,28,13)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{2}.Ac(:,29,14)  = [zeros(12,1);ones(12,1)];  % Session 2
        
        M{2}.name       = 'flex-flex';
        
        % --------------------------------------
        % Model3: Model with a fixed r=1 correlation (second session seq same as first)
        M{3}.type       = 'feature';
        M{3}.numGparams = 10;
        M{3}.Ac(:,1:6,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns      (theta_a)
        M{3}.Ac(:,7:12,2) = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   (theta_b)
        M{3}.Ac(:,1:6,3)  = [zeros(6);zeros(6);A;zeros(6)];     % Same trained2 pattterns       (theta_c)
        M{3}.Ac(:,7:12,4) = [zeros(6);zeros(6);zeros(6);A];     % Same untrained2 pattterns     (theta_d)
   
        M{3}.Ac(:,13,5)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1  (theta_e)
        M{3}.Ac(:,14,6)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  (theta_f)
        M{3}.Ac(:,13,7)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Same component as trained1  (theta_i)
        M{3}.Ac(:,14,8)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Same component as untrained1  (theta_j)
        M{3}.Ac(:,17,9)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{3}.Ac(:,18,10)  = [zeros(12,1);ones(12,1)];  % Session 2
        
        M{3}.name        = 'one-one';
            
end
function M = pcm_corrModel_seqType_fixed
% --------------------------------------
% Models separately trained / untrained sequences
% also session / sequence type factor 

        % Model1: Model with independent trained / untrained patterns
        % across sessions - 0 correlation
        
        % Model 1 - independent sequence patterns between sessions
        M{1}.type       = 'feature';
        M{1}.name       = 'ind';
        M{1}.numGparams = 2;
        A=zeros(6);
        for i=1:6
            A(i,i)=1;
        end;
        M{1}.Ac(:,1:6 ,1)  = [A;zeros(6)];     % Unique sess1 sequence patterns    (theta_a) (either trained or untrained)
        M{1}.Ac(:,7:12,2)  = [zeros(6);A];     % Unique sess2 sequence pattterns   (theta_b)
        
        % --------------------------------------
        % Model2: Model with a flexible across-sess correlation for sequences 
        M{2}.type       = 'feature';
        M{2}.numGparams = 3;
        M{2}.name       = 'flex';
        M{2}.Ac(:,1:6 ,1)  = [A;zeros(6)];     % Unique sess1 sequence patterns    (theta_a) (either trained or untrained)
        M{2}.Ac(:,7:12,2)  = [zeros(6);A];     % Unique sess2 sequence pattterns   (theta_b)
        M{2}.Ac(:,1:6,3)   = [zeros(6);A];     % Same sess1 patterns   (theta_c)
        
        % --------------------------------------
        % Model3: Model with a fixed r=1 correlation (second session seq same as first)
        M{3}.type       = 'feature';
        M{3}.numGparams = 2;
        M{3}.name        = 'one';
        M{3}.Ac(:,1:6,1)  = [A;zeros(6)];     % Unique sess1 patterns      (theta_a)
        M{3}.Ac(:,1:6,2)  = [zeros(6);A];     % Same sess2 pattterns       (theta_b)
        
            
end
function M = pcm_corrModel_seqType_random
% --------------------------------------
% Models separately trained / untrained sequences
% also session / sequence type factor 

        % Model1: Model with independent trained / untrained patterns
        % across sessions - 0 correlation
        
        % Model 1 - independent sequence patterns between sessions
        M{1}.type       = 'feature';
        M{1}.name       = 'ind';
        M{1}.numGparams = 4;
        A=zeros(6);
        for i=1:6
            A(i,i)=1;
        end;
        M{1}.Ac(:,1:6 ,1)  = [A;zeros(6)];     % Unique sess1 sequence patterns    (theta_a) (either trained or untrained)
        M{1}.Ac(:,7:12,2)  = [zeros(6);A];     % Unique sess2 sequence pattterns   (theta_b)
        M{1}.Ac(:,13,3)    = [ones(6,1);zeros(6,1)];    % Sess1
        M{1}.Ac(:,14,4)    = [zeros(6,1);ones(6,1)];    % Sess2
        
        % --------------------------------------
        % Model2: Model with a flexible across-sess correlation for sequences 
        M{2}.type       = 'feature';
        M{2}.numGparams = 5;
        M{2}.name       = 'flex';
        M{2}.Ac(:,1:6 ,1)  = [A;zeros(6)];     % Unique sess1 sequence patterns    (theta_a) (either trained or untrained)
        M{2}.Ac(:,7:12,2)  = [zeros(6);A];     % Unique sess2 sequence pattterns   (theta_b)
        M{2}.Ac(:,1:6,3)   = [zeros(6);A];     % Same sess1 patterns   (theta_c)
        M{2}.Ac(:,13,4)    = [ones(6,1);zeros(6,1)];    % Sess1
        M{2}.Ac(:,14,5)    = [zeros(6,1);ones(6,1)];    % Sess2
        
        % --------------------------------------
        % Model3: Model with a fixed r=1 correlation (second session seq same as first)
        M{3}.type       = 'feature';
        M{3}.numGparams = 4;
        M{3}.name        = 'one';
        M{3}.Ac(:,1:6,1)  = [A;zeros(6)];     % Unique sess1 patterns      (theta_a)
        M{3}.Ac(:,1:6,2)  = [zeros(6);A];     % Same sess2 pattterns       (theta_b)
        M{3}.Ac(:,7,3)    = [ones(6,1);zeros(6,1)];    % Sess1
        M{3}.Ac(:,8,4)    = [zeros(6,1);ones(6,1)];    % Sess2
            
end
function M = pcm_corrModel_seqType_fixed_indSeq
% --------------------------------------
% Models separately trained / untrained sequences
% also session / sequence type factor 

        % Model1: Model with independent trained / untrained patterns
        % across sessions - 0 correlation
        
        % Model 1 - independent sequence patterns between sessions
        M{1}.type       = 'feature';
        M{1}.name       = 'ind';
        M{1}.numGparams = 12;
        
        for i=1:6
            A=zeros(6);
            A(i,i)=1;
            M{1}.Ac(:,1:6 ,i)    = [A;zeros(6)];       % Unique seq1 patterns   (theta_a)
            M{1}.Ac(:,7:12,6+i) = [zeros(6);A];        % Unique seq2 pattterns   (theta_c)
        end;
        
        % --------------------------------------
        % Model2: Model with a flexible across-sess correlation for sequences 
        M{2}.type       = 'feature';
        M{2}.numGparams = 18;
        M{2}.name       = 'flex';
        
        for i=1:6
            A=zeros(6);
            A(i,i)=1;
            M{2}.Ac(:,1:6 ,i)    = [A;zeros(6)];       % Seq1 patterns   (theta_a)
            M{2}.Ac(:,7:12,6+i) = [zeros(6);A];       % Unique seq2 pattterns   (theta_b)
            M{2}.Ac(:,1:6 ,12+i)  = [zeros(6);A];       % Same seq2 patterns  (theta_c)
        end;
        
        % --------------------------------------
        % Model3: Model with a fixed r=1 correlation (second session seq same as first)
        M{3}.type       = 'feature';
        M{3}.numGparams = 12;
        M{3}.name        = 'one';
        
        for i=1:6
            A=zeros(6);
            A(i,i)=1;
            M{3}.Ac(:,1:6 ,i)    = [A;zeros(6)]; % Seq1 finger patterns   (theta_a)
            M{3}.Ac(:,1:6 ,6+i)  = [zeros(6);A]; % Same seq2 patterns  (theta_b)
        end;
end       
function M = pcm_simpleModel_oneSess
    % only for one session at a time - to estimate signal / noise


    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:6 ,1)  = zeros(6);
    
    % Model 2 - sequence patterns
    M{2}.type       = 'feature';
    M{2}.name       = 'sequence';
    M{2}.numGparams = 1;
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{2}.Ac(:,1:6 ,1)  = A;     % Unique sequence patterns

    % --------------------------------------

end
function M = pcm_stabilityModel_generic

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
    
    % Model 2: First vs. second session
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'Session';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Session + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 4;
    M{3}.name       = 'Session+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    M{3}.Ac(:,3:8,3)  = [A;zeros(6)];      % Unique sess1 sequence patterns
    M{3}.Ac(:,9:14,4)  = [zeros(6);A];     % Unique sess2 sequence pattterns
    
    % Model 4: Session + sequence specific + correlation in sequences
    M{4}.type         = 'feature';
    M{4}.numGparams   = 5;
    M{4}.name         = 'Session+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    M{4}.Ac(:,3:8,3)  = [A;zeros(6)];       % Unique sess1 sequence patterns
    M{4}.Ac(:,9:14,4) = [zeros(6);A];       % Unique sess2 sequence pattterns
    M{4}.Ac(:,3:8,5)  = [zeros(6);A];     % Correlation sess1-sess2
    
    % Model 5: Session + sequence specific + PERFECT correlation in sequences
    M{5}.type         = 'feature';
    M{5}.numGparams   = 4;
    M{5}.name         = 'Session+Seq+PerfectCorr';
    M{5}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{5}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    M{5}.Ac(:,3:8,3)  = [A;zeros(6)];       % Unique sess1 sequence patterns
    M{5}.Ac(:,3:8,4) = [zeros(6);A];        % Identical sess2 sequence pattterns

    
end
function M = pcm_stabilityModel_specific
    
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: First vs. second session
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'Session';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Session + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 14;
    M{3}.name       = 'Session+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,3:8,2+i)  = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{3}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    end;

    % Model 4: Session + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 20;
    M{4}.name         = 'Session+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per session
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,3:8,2+i)   = [A;zeros(6)];     % Unique sess1 sequence patterns
        M{4}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique sess2 sequence pattterns
        M{4}.Ac(:,3:8,14+i)  = [zeros(6);A];     % Correlation sess1-sess2
    end;

    % Model 5: Session + sequence specific + PERFECT correlation in session
    M{5}.type         = 'feature';
    M{5}.numGparams   = 14;
    M{5}.name         = 'Session+Seq+PerfectCorr';
    M{5}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{5}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{5}.Ac(:,3:8,2+i)   = [A;zeros(6)];       % Unique sess1 sequence patterns
        M{5}.Ac(:,3:8,8+i)   = [zeros(6);A];       % Same sess2 sequence pattterns
    end;
    
end



function M = pcm_corrHyperModel
% --------------------------------------
        % Model1: Model with independent trained / untrained patterns
        % across sessions - 0 correlation
        
        % Model 1 - independent sequence patterns between sessions
        M{1}.name       = 'ind-ind'; 
        M{1}.type       = 'feature';
        M{1}.numGparams = 10;
        A=zeros(6);
        for i=1:6
            A(i,i)=1;
        end;
        M{1}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns      (theta_a)
        M{1}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   (theta_b)
        M{1}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns     (theta_c)
        M{1}.Ac(:,19:24,4) = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 pattterns   (theta_d)
   
        M{1}.Ac(:,25,5)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    (theta_e)
        M{1}.Ac(:,26,6)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  (theta_f)
        M{1}.Ac(:,27,7)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2    (theta_i)
        M{1}.Ac(:,28,8)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2  (theta_j)
        M{1}.Ac(:,29,9)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{1}.Ac(:,30,10)  = [zeros(12,1);ones(12,1)];  % Session 2

        
        % --------------------------------------
        % Model2: Model with a flexible across-sess correlation for sequences 
        M{2}.name       = 'flex-flex';
        M{2}.type       = 'feature';
        M{2}.numGparams = 14;
        
        M{2}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns      (theta_a)
        M{2}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   (theta_b)
        M{2}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns     (theta_c)
        M{2}.Ac(:,19:24,4) = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 pattterns   (theta_d)
        M{2}.Ac(:,1:6,5)   = [zeros(6);zeros(6);A;zeros(6)];     % Same trained1 patterns        (theta_e)
        M{2}.Ac(:,7:12,6)  = [zeros(6);zeros(6);zeros(6);A];     % Same untrained1 pattterns     (theta_f)
        
        M{2}.Ac(:,25,7)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1        (theta_g)
        M{2}.Ac(:,26,8)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1      (theta_h)
        M{2}.Ac(:,27,9)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2        (theta_i)
        M{2}.Ac(:,28,10) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2      (theta_j)
        M{2}.Ac(:,25,11) = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall SAME component trained1   (theta_k)
        M{2}.Ac(:,26,12) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall SAME component untrained1 (theta_l)
        M{2}.Ac(:,28,13)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{2}.Ac(:,29,14)  = [zeros(12,1);ones(12,1)];  % Session 2
        
        
        % --------------------------------------
        % Model3: Model with a fixed r=1 correlation (second session seq same as first)
        M{3}.name       = 'one-one';
        M{3}.type       = 'feature';
        M{3}.numGparams = 10;
        M{3}.Ac(:,1:6,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns      (theta_a)
        M{3}.Ac(:,7:12,2) = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   (theta_b)
        M{3}.Ac(:,1:6,3)  = [zeros(6);zeros(6);A;zeros(6)];     % Same trained2 pattterns       (theta_c)
        M{3}.Ac(:,7:12,4) = [zeros(6);zeros(6);zeros(6);A];     % Same untrained2 pattterns     (theta_d)
   
        M{3}.Ac(:,13,5)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1  (theta_e)
        M{3}.Ac(:,14,6)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  (theta_f)
        M{3}.Ac(:,13,7)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Same component as trained1  (theta_i)
        M{3}.Ac(:,14,8)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Same component as untrained1  (theta_j)
        M{3}.Ac(:,17,9)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{3}.Ac(:,18,10)  = [zeros(12,1);ones(12,1)];  % Session 2

        
        
        % ------------ COMBINATION MODELS -----------------
        % Model 4: independent trained, flexible untrained
        M{4}.name       = 'ind-flex';
        M{4}.type       = 'feature';
        M{4}.numGparams = 12;
        M{4}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns     
        M{4}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   
        M{4}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns    
        M{4}.Ac(:,19:24,4) = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 pattterns   
        M{4}.Ac(:,7:12,5)  = [zeros(6);zeros(6);zeros(6);A];     % Same untrained1 pattterns     
   
        M{4}.Ac(:,25,6)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    
        M{4}.Ac(:,26,7)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  
        M{4}.Ac(:,27,8)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2    
        M{4}.Ac(:,28,9)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2  
        M{4}.Ac(:,26,10) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall SAME component untrained1 
        M{4}.Ac(:,29,11)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{4}.Ac(:,30,12)  = [zeros(12,1);ones(12,1)];  % Session 2
        
        
        % --------------------------------
        % Model 5: independent trained, perfect untrained
        M{5}.name       = 'ind-one';
        M{5}.type       = 'feature';
        M{5}.numGparams = 10;
        M{5}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns     
        M{5}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   
        M{5}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns     
        M{5}.Ac(:,7:12,4)  = [zeros(6);zeros(6);zeros(6);A];     % Same untrained1 pattterns     
   
        M{5}.Ac(:,19,5)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    
        M{5}.Ac(:,20,6)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  
        M{5}.Ac(:,21,7)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2    
        M{5}.Ac(:,20,8) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall SAME component untrained1 
        M{5}.Ac(:,22,9)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{5}.Ac(:,23,10)  = [zeros(12,1);ones(12,1)];  % Session 2
        
        % --------------------------------
        % Model 6: flexible trained, independent untrained
        M{6}.name       = 'flex-ind';
        M{6}.type       = 'feature';
        M{6}.numGparams = 12;
        M{6}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns     
        M{6}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   
        M{6}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns    
        M{6}.Ac(:,19:24,4) = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 pattterns   
        M{6}.Ac(:,1:6,5)  = [zeros(6);zeros(6);A;zeros(6)];     % Same trained1 pattterns     
   
        M{6}.Ac(:,25,6)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    
        M{6}.Ac(:,26,7)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  
        M{6}.Ac(:,27,8)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2    
        M{6}.Ac(:,28,9)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2  
        M{6}.Ac(:,25,10) = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall SAME component trained1 
        M{6}.Ac(:,29,11)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{6}.Ac(:,30,12)  = [zeros(12,1);ones(12,1)];  % Session 2
        
        % --------------------------------
        % Model 7: flexible trained, perfect untrained
        M{7}.name       = 'flex-one';
        M{7}.type       = 'feature';
        M{7}.numGparams = 12;
        M{7}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns     
        M{7}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   
        M{7}.Ac(:,13:18,3) = [zeros(6);zeros(6);A;zeros(6)];     % Unique trained2 pattterns
        M{7}.Ac(:,7:12,4)  = [zeros(6);zeros(6);zeros(6);A];     % Same untrained1 pattterns   
        M{7}.Ac(:,1:6,5)   = [zeros(6);zeros(6);A;zeros(6)];     % Same trained1 patterns  
   
        M{7}.Ac(:,19,6)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    
        M{7}.Ac(:,20,7)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  
        M{7}.Ac(:,21,8)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall component trained2    
        M{7}.Ac(:,20,9)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall SAME component untrained1 
        M{7}.Ac(:,19,10) = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall SAME component trained1  
        M{7}.Ac(:,22,11) = [ones(12,1);zeros(12,1)];  % Session 1
        M{7}.Ac(:,23,12) = [zeros(12,1);ones(12,1)];  % Session 2
        
        % --------------------------------
        % Model 8: perfect trained, independent untrained
        M{8}.name       = 'one-ind';
        M{8}.type       = 'feature';
        M{8}.numGparams = 10;
        M{8}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns     
        M{8}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   
        M{8}.Ac(:,1:6,3)   = [zeros(6);zeros(6);A;zeros(6)];     % Same trained1 pattterns     
        M{8}.Ac(:,13:18,4) = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 pattterns     
   
        M{8}.Ac(:,19,5)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    
        M{8}.Ac(:,20,6)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  
        M{8}.Ac(:,19,7)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall SAME component trained1    
        M{8}.Ac(:,21,8)  = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2
        M{8}.Ac(:,22,9)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{8}.Ac(:,23,10) = [zeros(12,1);ones(12,1)];  % Session 2
        
        
        % --------------------------------
        % Model 9: perfect trained, flexible untrained
        M{9}.name       = 'one-flex';
        M{9}.type       = 'feature';
        M{9}.numGparams = 12;
        M{9}.Ac(:,1:6 ,1)  = [A;zeros(6);zeros(6);zeros(6)];     % Unique trained1 patterns     
        M{9}.Ac(:,7:12,2)  = [zeros(6);A;zeros(6);zeros(6)];     % Unique untrained1 pattterns   
        M{9}.Ac(:,1:6,3) = [zeros(6);zeros(6);A;zeros(6)];       % Same trained1 pattterns
        M{9}.Ac(:,13:18,4)  = [zeros(6);zeros(6);zeros(6);A];     % Unique untrained2 patterns  
        M{9}.Ac(:,7:12,5)   = [zeros(6);zeros(6);zeros(6);A];     % Same untrained1 patterns  
   
        M{9}.Ac(:,19,6)  = [ones(6,1);zeros(6,1);zeros(6,1);zeros(6,1)];  % Overall component trained1    
        M{9}.Ac(:,20,7)  = [zeros(6,1);ones(6,1);zeros(6,1);zeros(6,1)];  % Overall component untrained1  
        M{9}.Ac(:,19,8)  = [zeros(6,1);zeros(6,1);ones(6,1);zeros(6,1)];  % Overall SAME component trained1   
        M{9}.Ac(:,21,9) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall component untrained2 
        M{9}.Ac(:,20,10) = [zeros(6,1);zeros(6,1);zeros(6,1);ones(6,1)];  % Overall SAME component untrained1  
        M{9}.Ac(:,22,11)  = [ones(12,1);zeros(12,1)];  % Session 1
        M{9}.Ac(:,23,12)  = [zeros(12,1);ones(12,1)];  % Session 2
        
end
function C = calcDist(D,betaW,G)
% calculates G and distances for all sequences, trained / untrained /
% between the two
% INPUT: D - structure with run / cond etc.
%        betaW - all betas
%        G - overall G structure (12x12 - all seq)
% OUTPUT: Do - containes distances, eigenvalues

    % calculate trained / untrained G
    [G_train, Sig_train] = pcm_estGCrossval(betaW(D.seqType==1,:),D.run(D.seqType==1),D.seqNumb(D.seqType==1));
    [G_untrain, Sig_untrain] = pcm_estGCrossval(betaW(D.seqType==2,:),D.run(D.seqType==2),D.seqNumb(D.seqType==2));
    % calculate average train / average untrain pattern - beta_seqType
    indx=1;
    indx_seqType=[];
    indx_run=[];
    for rr = 1:8 % number of runs
        for seq = 1:max(unique(D.seqType))
            beta_seqType(indx,:)=mean(betaW(D.seqType==seq & D.run==rr,:),1);
            indx=indx+1;
            indx_seqType=[indx_seqType;seq]; % 1 - train, 2 - untrain
            indx_run=[indx_run;rr]; % 1-8 func runs
        end
    end
    % calculate G between sequence types
    [G_seqType, Sig_seqType] = pcm_estGCrossval(beta_seqType,indx_run,indx_seqType);

    % indicator matrices
    ind6=indicatorMatrix('allpairs',[1:6]); % Trained / Untrained
    ind2=indicatorMatrix('allpairs',[1:2]); % SeqType
    ind12=indicatorMatrix('allpairs',[1:12]); % All

    % calculate distances
    dist_train = rsa.rdm.squareRDM(diag(ind6*G_train*ind6'));
    dist_untrain = rsa.rdm.squareRDM(diag(ind6*G_untrain*ind6'));
    dist_cross = diag(ind2*G_seqType*ind2');
    dist_all = rsa.rdm.squareRDM(diag(ind12*G*ind12'));

    C.dist_train = mean(dist_train(:));
    C.dist_untrain = mean(dist_untrain(:));
    C.dist_cross = mean(dist_cross);
    C.dist_all = mean(dist_all(:));
 
    % calculate eigenvalues
    H=eye(6)-ones(6,6)./6;  % centering matrix!
    G_trainCent = H*G_train*H;  % double centered G matrix - rows and columns
    C.eigTrain = sort(eig(G_trainCent)','descend');    % sorted eigenvalues - trainedSeq
    G_untrainCent = H*G_untrain*H;  % double centered G matrix - rows and columns
    C.eigUntrain = sort(eig(G_untrainCent)','descend'); % sorted eigenvalues - untrainedSeq
end
function PermDist = randomSeqDist(D,betaW,G)
% randomises labels for seqNumb and seqType, calculates distances
PermDist=[];
S=D;
numPerm=1000;
permIndx=randSeq(numPerm);

for k = 1:numPerm
   % re-label seqNumb / seqType
   S.seqNumb=repmat(permIndx(k,:)',max(S.run),1);
   S.seqType(S.seqNumb<7)=1;
   S.seqType(S.seqNumb>6)=2;
   
   %submit for distance calculation on shuffled set
   PD = calcDist(D,betaW,G);
   PD.numPerm = k;
   
   PermDist=addstruct(PermDist,PD);
end

end
function P = randSeq(numPerm)
% generate permutations of numbers 1-12
% used for shuffling trained / untrained labels for sequences
    set1  = [];
    fullset = 1:12;
    for k = 1:numPerm
       newset         = fullset;
       newset(set1) = [];
       set1         = newset(randperm(length(newset), 6)); % find first set (trained)
       set1Idx      = ismember(fullset,set1); % index of trained numbers
       set2         = fullset(~set1Idx); % second set of 6 numbers (untrained)

       P(k,:)=[set1 set2];
    end
end
function C = pcm_correlation(Data,partVec,condVec,M,runEffect,M_type)
sn=1:size(Data,2);
% --------------------------------------
% 1. Empirical correlation
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec);
    b = pinv(Z)*Data{p};           % Estimate mean activities
    b(1:6,:)  = bsxfun(@minus,b(1:6,:) ,mean(b(1:6,:))); % Subtract mean per condition - first exe
    b(7:12,:) = bsxfun(@minus,b(7:12,:),mean(b(7:12,:))); % second exe
    G=cov(b');
    C.r_naive(p,1) = calcCorr(G);
end;

% --------------------------------------
% 2. Crossvalidated correlation
 condSeqTypeVec=[ones(96/2,1);ones(96/2,1)*2];
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec);
    %Z_seqType=pcm_indicatorMatrix('identity',condSeqTypeVec);
    % Subtract mean for each condition and run
    X = pcm_indicatorMatrix('identity',partVec*2+(condVec>6)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec,condVec);
    Gcv_seqType(:,:,p)=pcm_estGCrossval(Data{p},partVec,condSeqTypeVec);
    C.r_crossval(p,1)=calcCorr(pcm_makePD(Gcv(:,:,p)));
    G_seqType=pcm_makePD(Gcv_seqType(:,:,p));
    C.r_crossval_seqType(p,1)=calcCorr_thetas(G_seqType(1,1),G_seqType(2,2),G_seqType(1,2));
end;

% --------------------------------------
% 3. Fit model 2  and infer correlations from the parameters
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);

% Get the correlations
switch M_type
    case 'specific'
        var1       = (theta{1}(3:8,:).^2)';
        var2       = (theta{1}(9:14,:).^2+theta{1}(15:20,:).^2)';
        cov12      = (theta{1}(3:8,:).*theta{1}(15:20,:))';
        
        %var1       = (theta{1}(1:6,:).^2)';
        %var2       = (theta{1}(7:12,:).^2+theta{1}(13:18,:).^2)';
        %cov12      = (theta{1}(1:6,:).*theta{1}(13:18,:))';
        C.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
    case 'generic'
        var1       = (theta{1}(3,:).^2)';
        var2       = (theta{1}(4,:).^2+theta{1}(5,:).^2)';
        cov12      = (theta{1}(3,:).*theta{1}(5,:))';
        %var1       = (theta{1}(1,:).^2)';
        %var2       = (theta{1}(2,:).^2+theta{1}(3,:).^2)';
        %cov12      = (theta{1}(1,:).*theta{1}(3,:))';
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

    v1=th1^2;
    v2=th2^2+th3^2;
    cv=th1*th3;
    r=cv./sqrt(v1*v2);
end


