function varargout=sml1_imana_simu(what,varargin)

% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';

behavDir        =[baseDir '/behavioral_data/data'];            
imagingDir      =[baseDir '/imaging_data'];              
imagingDirRaw   =[baseDir '/imaging_data_raw'];           
dicomDir        =[baseDir '/imaging_data_dicom'];         
anatomicalDir   =[baseDir '/anatomicals'];       
fieldmapDir     =[baseDir '/fieldmaps/'];
freesurferDir   =[baseDir '/surfaceFreesurfer'];          
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia'];
suitDir         =[baseDir '/suit'];
physioDir       =[baseDir '/physio'];
pcmDir          =[baseDir '/pcm_stats'];
QCDir           =[baseDir '/quality_control'];

% update glmDir when adding new glms
glmLocDir       ={[baseDir '/glmLoc/glmL1'],[baseDir '/glmLoc/glmL2'],[baseDir '/glmLoc/glmL3']};   % localiser glm
glmLocSessDir   ={[baseDir '/glmLocSess/glmLocSess1'],[baseDir '/glmLocSess/glmLocSess2'],[baseDir '/glmLocSess/glmLocSess3'],[baseDir '/glmLocSess/glmLocSess4']}; % one glm for loc run per session
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session
glmFoSExDir     ={[baseDir '/glmFoSEx/glmFoSEx1'],[baseDir '/glmFoSEx/glmFoSEx2'],[baseDir '/glmFoSEx/glmFoSEx3'],[baseDir '/glmFoSEx/glmFoSEx4']};    
glmTrialDir     ={[baseDir '/glmTrial/glmTrial1'],[baseDir '/glmTrial/glmTrial2'],[baseDir '/glmTrial/glmTrial3'],[baseDir '/glmTrial/glmTrial4']};

% ------------------------- Experiment Info -------------------------------
numDummys  = 4;        % per run
numTRs     = [440 440 440 160 440 440 440 160 440 440,...
              440 440 440 160 440 440 440 160 440 440,...
              440 440 440 160 440 440 440 160 440 440,...
              440 440 440 160 440 440 440 160 440 440];             
% per functional run (includes dummies)  
% 440 - task; 160 - localiser

% Stimuli - numbers given in SeqNumb
num_train = 1:6;
num_untrain = 7:12;   
num_seq = 1:12;
num_fing = 13:17;
num_seqtype = 2;
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
regname         = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp','CaudateN' 'Pallidum', 'Putamen' 'Thalamus'};
%regname         = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp','CaudateN' 'Pallidum', 'Putamen' 'Thalamus','CIV','CV','CVI'};
regname_cortex  = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp'};
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen', 'Thalamus'};
regname_cerebellum = {'LobIV','LobV','LobVI'};
numregions_surf = 8;
numregions_BG   = 4;
numregions_cerebellum = 3;
numregions = numregions_surf+numregions_BG;
%numregions = numregions_surf+numregions_BG+numregions_cerebellum; - after adding cerebellum       
regSide=[ones(1,8) ones(1,8)*2]; % 1-left, 2-right
regType=[1:8  1:8]; % cortical areas: 1-8, BG: 8-12, cereb: 13-15


% ------------------------- Freesurfer things -----------------------------         
atlasA    = 'x';                                                            % freesurfer filename prefix
atlasname = 'fsaverage_sym';                                                % freesurfer average atlas
hemName   = {'LeftHem','RightHem'};                                         % freesurfer hemisphere folder names    

% ------------------------- Subject things --------------------------------

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
              's21','s22','s23','s24','s25','s26','s27','s28'};  


% Other random notes %

% Dicom system changed with s05 (sessN 3,4), s06 (sessN 1-4) and s07 (sessN 2-4)


% ------------------------------ Analysis Cases --------------------------------
   switch(what)
    
    case 'PCM_simulate_repsupModels' %---------------------- PCM on CORRELATION 1st-2nd REPETITION ----------------
        runEffect  = 'fixed';
        mNum=2; % which model used for data generation - 1: ind, 2: flex, 3: same
        noise=1;
        scale=[0:0.05:0.2];
        theta=[-0.1 0.2 -0.2]';        %[0.1 0.2 0.5] - very high corr [0.1 0.2 0] - 0 corr
        part=8;
        cond=6;
        modelCorr={'ind','flex','perfect'};
        fitModel={'generic','specific'}; 
        thetaScale=[1 1 1 1 1 1]';
        trueModel=2; % generic - 1 param for all seq, or specific - 1 param for each seq
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
        switch(sessType)
            case 'within'
                D.partVec = [kron([1:part]',ones(cond,1)); kron([1:part]',ones(cond,1))];  % Partitions
            case 'between'
                D.partVec = [kron([1:part]',ones(cond,1)); kron([1+part:part+part]',ones(cond,1))];  % Partitions
        end 
        sn = 25; % 25 subjects
        % create both models
        M{1}=pcm_corrModel;
        M{2}=pcm_corrModel_indSeq;
        
        TT=[];
        
        for s=1:length(scale)
            for i=1:numSim
                switch trueModel
                    case 2
                        theta = [repmat(theta(1),6,1).*thetaScale; repmat(theta(2),6,1).*thetaScale; repmat(theta(3),6,1).*thetaScale];
                end
                Data = pcm_generateData_eva(M{trueModel}{mNum},theta,D,sn,scale(s),noise);
                
                switch(sessType)
                    case 'between' % scale data for each session
                        for ss=1:sn
                            Data{ss}(D.partVec<9,:)=Data{ss}(D.partVec<9,:).*sessScale(1);
                            Data{ss}(D.partVec>8,:)=Data{ss}(D.partVec>8,:).*sessScale(2);
                        end
                end
                trueCorr  = calcCorr_thetas(theta(1),theta(2),theta(3));
                for f=1:length(fitModel)
                    T = pcm_fitModels(Data,M{f},D.partVec,D.condVec,runEffect,algorithm);
                    C = pcm_correlation(Data,D.partVec,D.condVec,M{f}{2},runEffect,f);                   
                    C.signalLevel=ones(size(T.SN))*scale(s);
                    C.fitModel=ones(size(T.SN))*f;
                    
                    C.trueCorr=ones(size(T.SN))*trueCorr;
                    TT=addstruct(TT,C);
                    TT=addstruct(TT,T);
                end
            end
        end

        keyboard;
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
        % figure of correlation vs. logBayes
        figure
        scatterplot(TT.r_model2,TT.bayesEst(:,2),'split',TT.signalLevel,'leg',{'signal 0','signal 0.05','signal 0.1','signal 0.15','signal 0.2'})        
        drawline(0,'dir','vert'); drawline(0,'dir','horz');
        drawline(unique(TT.trueCorr),'dir','vert','color',[1 0 0]);
        xlabel('PCM correlation'); ylabel('PCM logBayes flex model');
        
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
    case 'PCM_data_repsupModels' % repsup models
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1:8];
        sn=[1:9,11:25];
        sessN=1; % need to be two sessions at the time
        seqType='trained';
        models={'generic','specific'};
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
                    
                    indx = ones(size(D.run));
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
                end
                T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                %  C = pcm_correlation(Data,partVec,condVec,M{2},runEffect,M_type);
                C = pcm_correlation(Data,partVec,condVec,M{4},runEffect,M_type);
                T.roi = ones(size(T.SN))*r;
                AllReg=addstruct(AllReg,C);
                AllReg=addstruct(AllReg,T);
            end % region
            
            %remove some fields
            a1='reg'; a2='theta'; a3='theta_hat'; a4='thetaCr';
            AllReg=rmfield(AllReg,a1); AllReg=rmfield(AllReg,a2); AllReg=rmfield(AllReg,a3); AllReg=rmfield(AllReg,a4);
            % save output
            save(fullfile(pcmDir,sprintf('PCM_repsup_reliability_NEW_%s_%s_sess%d_%s.mat',models{M_type},seqType,ss,runEffect)),'-struct','AllReg');
        end % session
    case 'PCM_simulateCorr'
        runEffect = 'fixed';
        mNum=2; % 1 - seq only, 2 - seq+meanSeq
        noise=1;
        scale=[0:0.05:0.5];
        %scale=[0:0.1:1];
        corrTheta=0.1;
        thetaFix=[-0.1 0.2 0.1];
        part=8;
        cond=6;
        modelType='generic'; 
        fitModel={'seq','seq+meanSeq'}; 
        voxSelect = 0;  % choose only voxels with top distances
        thetaScale=[1 1 1 1 1 1]';
        dataGen='pcm'; % or pcm
        
        sessType = 'within'; % within or between sessions
        
        vararginoptions(varargin,{'runEffect','thetaFix','noise','scale','mNum','corrTheta','sessType','thetaScale','modelType','voxSelect','dataGen'});
        
        numSim=1; % number of simulations
        % model specifications
        D.numVox  = 1000;  
       % D.condVec = [kron(ones(part,1),[1:cond]');kron(ones(part,1),[1+cond:cond+cond]')];   % Conditions
        D.condVec = [kron(ones(part,1),[1:cond+cond]')];
        switch(sessType)
            case 'within'
              %  D.partVec = [kron([1:part]',ones(cond,1)); kron([1:part]',ones(cond,1))];  % Partitions
                D.partVec = [kron([1:part]',ones(cond+cond,1))];
            case 'between'
                D.partVec = [kron([1:part]',ones(cond,1)); kron([1+part:part+part]',ones(cond,1))];  % Partitions
        end
        sn = 25; % 25 subjects
        
        switch modelType    % create both models - simple + meanSeq
            case 'generic'      
                M1=pcm_corrModel;
                M2=pcm_corrModel_sharedSeqType;
                mf=1;
            case 'specific'
                M1=pcm_corrModel_indSeq;
                M2=pcm_corrModel_indSeq_sharedSeqType;
                mf=2;
        end
        M{1}=M1{2}; % simple model
        M{2}=M2{2}; % added common pattern
        
        TT=[];
        SS=[];
        for t=1:length(corrTheta)
            for s=1:length(scale)
                for i=1:numSim
                    switch modelType
                        case 'generic'
                            theta = thetaFix';
                        case 'specific'
                            if mNum==1
                            theta = [repmat(thetaFix(1),6,1).*thetaScale; repmat(thetaFix(2),6,1).*thetaScale; repmat(thetaFix(3),6,1).*thetaScale];
                            elseif mNum==2
                            theta = [repmat(thetaFix(1),6,1).*thetaScale; repmat(thetaFix(2),6,1).*thetaScale; repmat(thetaFix(3),6,1).*thetaScale; thetaFix(4); thetaFix(5)];    
                            end
                    end
                    
                    switch dataGen
                        case 'pcm'
                            %Data = pcm_generateData_eva(M{1},theta,D,sn,scale(s),noise);
                             Data = pcm_generateData_eva(M{2},theta,D,sn,scale(s),noise);
                            trueCorr  = calcCorr_thetas(thetaFix(1),thetaFix(2),thetaFix(3));
                        case 'randn'
                            for p=1:sn
                                Data{p}=randn(96,D.numVox);
                            end
                            trueCorr  = 0;
                    end
                    % selection
                    if voxSelect==1
                        condVec=D.condVec;
                        condVec(condVec>6)=condVec(condVec>6)-6;
                        Data = voxelSelect(Data,D.partVec,condVec);
                        %Data = voxelSelect(Data,D.partVec,D.condVec);
                    end
                    
                    % fit both models
                    for f=1:length(fitModel)
                        C = pcm_correlation(Data,D.partVec,D.condVec,M{f},runEffect,mf);
                        
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

        % save structures
        save(fullfile(pcmDir,sprintf('PCM_sim_repsup_corr_%s_%s_voxSelect%d_%s.mat',modelType,fitModel{mNum},voxSelect,dataGen)),'-struct','TT');
        save(fullfile(pcmDir,sprintf('PCM_sim_repsup_corr_summaryStats_%s_%s_voxSelect%d_%s.mat',modelType,fitModel{mNum},voxSelect,dataGen)),'-struct','SS');   
    case 'PCM_dataCorr'
        % calculate correlation between 1st and 2nd execution
        % naive, crossval, pcm correlation
        % generic or specific model - parameter for all 6 seq / per seq
        runEffect  = 'fixed';
        beta_choice = 'mw';
        reg = [1:8];
        sn=[1:9,11:25];
        sessN=1; % need to be two sessions at the time
        seqType='trained';
        modelType='generic'; %generic or specific
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','seqType','modelType'})
        
        switch seqType
            case 'trained'
                stIndx=1;
            case 'untrained'
                stIndx=2;
        end
        CC=[];
        for ss=sessN
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
                        cond(D.FoSEx==1)=cond(D.FoSEx==1)-6;
                    end
                    condVec{p} = cond(D.seqType==stIndx); % conditions
                    
                    partVec{p} = D.run(D.seqType==stIndx);
                    Data{p} = beta(:,D.seqType==stIndx)';  % Data is N x P (cond x voxels) - no intercept
                    fprintf('Extracted data for %s...\n',subj_name{sn(p)});
                end; % subj
                
                % construct models
                switch modelType
                    case 'generic'
                        M = pcm_corrModel_sharedSeqType;
                        mf=1;
                    case 'specific'
                        M = pcm_corrModel_indSeq_sharedSeqType;
                        mf=2;
                end
               
                % with mean pattern
                C = pcm_correlation(Data,partVec,condVec,M{2},runEffect,mf);
                C.SN=[1:size(C.r_naive,1)]';
                C.roi=ones(size(C.r_naive))*r;
                CC=addstruct(CC,C);
                fprintf('Done calculating corr for %s\n\n\n',regname_cortex{reg(r)});
            end       
            CC=rmfield(CC,'theta');
            % save output
            save(fullfile(pcmDir,sprintf('PCM_repsup_corr_sess%d_%s_%s.mat',ss,seqType,modelType)),'-struct','CC');
            fprintf('Done all regions sess-%d\n\n\n\n\n\n',ss);
        end
    
    case 'CROSSVAL_corr'
        % looking both into simple model (generic / specific)
        % and model with added pattern across sequences (shared)
        beta_choice = 'mw';
        reg = [1:8];
        sn=[1:9,11:25];
        sessN=[1:4]; 
        seqType='trained';
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','seqType'})
        
        switch seqType
            case 'trained'
                stIndx=1;
            case 'untrained'
                stIndx=2;
        end
        CC=[];
        RR=[];
        for r = reg
            for p=1:length(sn)
                for ss = sessN
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
                    
                    % first exe
                    Data1 = beta(:,D.FoSEx==1&D.seqType==stIndx)';  
                    partVec1 = D.run(D.FoSEx==1&D.seqType==stIndx);
                    condVec1 = D.seqNumb(D.FoSEx==1&D.seqType==stIndx);

                    % second exe
                    Data2 = beta(:,D.FoSEx==2&D.seqType==stIndx)';  
                    partVec2 = D.run(D.FoSEx==2&D.seqType==stIndx);
                    condVec2 = D.seqNumb(D.FoSEx==2&D.seqType==stIndx);
                    
                    if stIndx==2 % for untrained seq make condition vec 1-6
                        condVec1=condVec1-6;
                        condVec2=condVec2-6;
                    end
                    %across
                    d1a     = Data1(rem(partVec1,2)==1,:);
                    part1a  = partVec1(rem(partVec1,2)==1,:);
                    d2a     = Data2(rem(partVec1,2)==0,:);
                    part2a  = partVec2(rem(partVec2,2)==0,:);
                    dA      = [d1a;d2a];
                    partA   = [part1a;part2a];
                    d1b     = Data1(rem(partVec1,2)==0,:);
                    part1b  = partVec1(rem(partVec1,2)==0,:);
                    d2b     = Data2(rem(partVec1,2)==1,:);
                    part2b  = partVec2(rem(partVec2,2)==1,:);
                    dB      = [d1b;d2b];
                    partB   = [part1b;part2b];
                          
                    %pattern consistency
                 %   R1=rsa_patternConsistency(Data1,partVec1,condVec1);
                 %   R2=rsa_patternConsistency(Data2,partVec2,condVec2);
                    [R2w1 Rw1]=rsa_patternConsistency_crossval(Data1,partVec1,condVec1,'removeMean',0);
                    [R2w2 Rw2]=rsa_patternConsistency_crossval(Data2,partVec2,condVec2,'removeMean',0);
                    % across
                    [R2Aacr RAacr]=rsa_patternConsistency_crossval(dA,partA,condVec1,'removeMean',0);
                    [R2Bacr RBacr]=rsa_patternConsistency_crossval(dB,partB,condVec1,'removeMean',0);
                    
                    R2acr=mean([R2Aacr,R2Bacr]);
                    Racr=mean([RAacr,RBacr]);
                    
                    C.sn    = [p;p];
                   % C.pConsist=[R1;R2];
                  %  C.pConsistCross=[R1c;R2c];
                    C.r2    = [R2w1; R2w2];
                    C.r     = [Rw1; Rw2];
                    C.FoSEx = [1;2];
                    C.roi   = [r;r];
                    C.sessN = [ss;ss];
                    
                    R.r2            = R2acr;
                    R.r             = Racr;
                    R.r2_correct    = R2acr/sqrt(R2w1*R2w2);
                    R.r_correct     = Racr/sqrt(Rw1*Rw2);
                    if ~isreal(R.r2_correct) | ~isreal(R.r_correct)
                        R.r_correct  = NaN;
                        R.r2_correct = NaN;
                    end
                    R.sn            = p;
                    R.roi           = r;
                    R.sessN         = ss;
                    
                    CC=addstruct(CC,C);
                    RR=addstruct(RR,R);
                end; % session
            end; % subj
        end

        keyboard;
        % save structure
        save(fullfile(repSupDir,sprintf('patternConsist_%s',seqType)),'-struct','CC');
        save(fullfile(repSupDir,sprintf('patternConsistCorr_%s',seqType)),'-struct','RR');
        
        %plot
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line([CC.sessN>3 CC.sessN],CC.r,'subset',CC.roi==r,'split',CC.FoSEx,'leg',{'1st','2nd'},'leglocation','north');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Pattern consistency');
                xlabel('Session');
            else
                ylabel('');
                xlabel('');
            end
        end
        
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line([RR.sessN>3 RR.sessN],RR.r_correct,'subset',RR.roi==r,'leg',{'1st','2nd'},'leglocation','north');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Pattern consistency');
                xlabel('Session');
            else
                ylabel('');
                xlabel('');
            end
        end
    case 'PLOT_crossval_corr'
        reg=[1:8];
        sessN=[1:4];
        metric = 'r'; % r2, r, r2_correct, r_correct
        vararginoptions(varargin,{'sessN','reg','metric'});
        
        seqType={'trained','untrained'};
        
        TT=[];
        for st=1:2
            T=load(fullfile(repSupDir,sprintf('patternConsistCorr_%s',seqType{st})));
            T.seqType=ones(size(T.sn))*st;
            TT=addstruct(TT,T);
        end
        
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line([TT.sessN>3 TT.sessN],TT.(metric),'subset',TT.roi==r,'split',TT.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Pattern consistency');
                xlabel('Session');
            else
                ylabel('');
                xlabel('');
            end
        end
        
    case 'PLOT_pcm_simulateCorr'
        
        modelType='specific';
        voxSelect=0;
        dataGen='pcm';
        fitModel={'seq','seq+meanSeq'}; 
        fitM=1;
        vararginoptions(varargin,{'modelType','voxSelect','dataGen','fitM'});
        
 
        TT=load(fullfile(pcmDir,sprintf('PCM_sim_repsup_corr_%s_%s_voxSelect%d_%s.mat',modelType,fitModel{fitM},voxSelect,dataGen)));
        SS=load(fullfile(pcmDir,sprintf('PCM_sim_repsup_corr_summaryStats_%s_%s_voxSelect%d_%s.mat',modelType,fitModel{fitM},voxSelect,dataGen)));
        
        
        % correlation estimates - naive, crossval, model
        N=tapply(TT,{'signalLevel','fitModel'},{'r_naive','mean'});
        figure
        subplot(1,4,1)
        scatterplot(TT.signalLevel,TT.r_naive,'subset',TT.fitModel==1);
        hold on;
        scatterplot(N.signalLevel,N.r_naive,'subset',N.fitModel==1,'markercolor',[1 0 0]);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title(sprintf('Naive correlation %s - %dvoxSelect',dataGen,voxSelect));
        xlabel('Signal level'); ylabel('Corr values');
        
        C=tapply(TT,{'signalLevel','fitModel'},{'r_crossval','mean'});
        subplot(1,4,2)
        scatterplot(TT.signalLevel,TT.r_crossval,'subset',TT.fitModel==1);
        hold on;
        scatterplot(C.signalLevel,C.r_crossval,'subset',C.fitModel==1,'markercolor',[1 0 0]);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title('Crossval correlation');
        xlabel('Signal level'); 
        
        M=tapply(TT,{'signalLevel','fitModel'},{'r_model2','mean'});
        subplot(1,4,3)
        scatterplot(TT.signalLevel,TT.r_model2,'subset',TT.fitModel==1);
        hold on;
        scatterplot(M.signalLevel,M.r_model2,'subset',M.fitModel==1,'markercolor',[1 0 0]);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title(sprintf('PCM correlation - %s simple model',modelType));
        xlabel('Signal level');
        
        subplot(1,4,4)
        scatterplot(TT.signalLevel,TT.r_model2,'subset',TT.fitModel==2);
        hold on;
        scatterplot(M.signalLevel,M.r_model2,'subset',M.fitModel==2,'markercolor',[1 0 0]);
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline(TT.trueCorr,'dir','horz','color',[1 0 0]);
        title(sprintf('PCM correlation - %s with common pattern across seq',modelType));
        xlabel('Signal level');
        
                
        % summary stats
        leg_labels={'naive','crossval',sprintf('model %s simple',modelType),sprintf('model %s with shared seq pattern',modelType)};

        A=getrow(SS,SS.fitModel==1 | SS.corrType==3);
        A.corrType(A.fitModel==2&A.corrType==3)=4;
        figure
        subplot(1,3,1)
        plt.line(A.signalLevel,A.bias,'split',A.corrType,'leg','off');
        drawline(0,'dir','horz');
        ylabel('Bias'); xlabel('Signal level'); title(sprintf('True correlation %d - %s %d-voxSelect',A.trueCorr(1),dataGen,voxSelect));
        subplot(1,3,2)
        plt.line(A.signalLevel,A.var,'split',A.corrType,'leg','off');
        drawline(0,'dir','horz');
        ylabel('Variance'); xlabel('Signal level');
        subplot(1,3,3)
        plt.line(A.signalLevel,A.mse,'split',A.corrType,'leg',leg_labels,'leglocation','northeast');
        drawline(0,'dir','horz');
        ylabel('MSE'); xlabel('Signal level'); 
    case 'PLOT_pcm_corrMetric_evals'
        
        sessN=1;
        seqType='trained';
        modelType='generic'; % or specific
        data='real'; % load real data or simulation
        vararginoptions(varargin,{'sessN','seqType','modelType','data'});
        
        switch data
            case 'real'
                T=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_sess%d_%s_%s.mat',sessN,seqType,modelType)));
            case 'simulation'
                T=load(fullfile(pcmDir,sprintf('PCM_simulation_repsup_corr_%s.mat',modelType)));
        end
        figure
        subplot(2,2,1)
        scatterplot(T.r_naive,T.r_model2,'split',T.fitModel,'leg',{'simple','added shared seq pattern'},'draworig');
        ylabel('PCM corr'); xlabel('Naive correlation - mean subtracted for 1st/2nd'); title(sprintf('%s',modelType));
        subplot(2,2,2)
        scatterplot(T.r_naive_meanSubtract,T.r_model2,'split',T.fitModel,'leg',{'simple','added shared seq pattern'},'draworig');
        ylabel('PCM corr'); xlabel('Naive correlation - mean overall subtracted'); title(sprintf('%s',modelType));
        subplot(2,2,3)
        scatterplot(T.r_naive,T.r_crossval,'split',T.fitModel,'leg',{'simple','added shared seq pattern'},'draworig');
        ylabel('Crossval corr'); xlabel('Naive correlation - mean subtracted for 1st/2nd');
        subplot(2,2,4)
        scatterplot(T.r_naive_meanSubtract,T.r_crossval,'split',T.fitModel,'leg',{'simple','added shared seq pattern'},'draworig');
        ylabel('Crossval corr'); xlabel('Naive correlation - mean overall subtracted');
    case 'PLOT_pcm_logBayes_allSess'
         reg = [1:8];
         sessN=[1:4];
         seqType='trained';
         runEffect='fixed';
         modelType='specific'; % or specific
         vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType'});

        T=[];
        for t=sessN
            R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_%s_%s_sess%d_%s.mat',modelType,seqType,t,runEffect)));
            R.sessN=ones(size(R.SN))*t;
            T=addstruct(T,R);
        end

        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3];
        T2.bayesEst=T.bayesEst(:);
        T2.roi=[T.roi;T.roi;T.roi];
        T2.sessN=[T.sessN;T.sessN;T.sessN];
        % rearranging how the data structure is arranged
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line([T2.sessN>3 T2.sessN],T2.bayesEst,'subset',T2.modelInd>1&T2.roi==r,'split',T2.modelInd,'leg',{'flex','perfect'},'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel(sprintf('Log-Bayes %s',seqType));
            else
                ylabel('');
            end
        end
    case 'PLOT_pcm_corr_allSess'
        reg = [1:8];
        sessN=[1:4];
        seqType={'trained','untrained'};
        modelType='generic';
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType','fitM'});
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_sess%d_%s_%s.mat',t,seqType{st},modelType)));
               % R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_%s_%s_sess%d_%s.mat',modelType,seqType{st},t,runEffect)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        corrType={'naive','crossval','pcm'};
        corrVar=[T.r_naive, T.r_crossval3, T.r_model2];
        for ct=1:3 % three types of correlation 
            figure(ct)
            sub=1;
            for r=reg
                subplot(1,numel(reg),sub)
                plt.line([T.sessN>3 T.sessN],corrVar(:,ct),'subset',T.roi==r,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','north');
                plt.match('y');
                drawline(0,'dir','horz');
                title(sprintf('%s',regname{r}));
                if sub==1
                    ylabel(sprintf('%s corr',corrType{ct}));
                    xlabel('Session');
                else
                    ylabel('');
                end
                sub=sub+1;
            end     
        end
    case 'PLOT_pcm_corr_logBayes'
        reg = [1:8];
        sessN=[1:4];
        seqType={'trained','untrained'};
        runEffect='fixed';
        modelType='specific';
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType'});
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_%s_%s_sess%d_%s.mat',modelType,seqType{st},t,runEffect)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
    figure
    scatterplot(T.r_model2,T.bayesEst(:,2));
    xlabel('PCM correlation');
    ylabel('logBayes flexible model');
    drawline(0,'dir','horz');
      
    figure
    scatterplot(T.r_crossval,T.bayesEst(:,2));
    xlabel('Crossval correlation');
    ylabel('logBayes flexible model');
    drawline(0,'dir','horz');
    
    figure
    scatterplot(T.r_crossval,T.bayesEst(:,2));
    xlabel('Naive correlation');
    ylabel('logBayes flexible model');
    drawline(0,'dir','horz');
    case 'PLOT_pcm_logBayes_NEW'
         reg = [1:8];
         sessN=[1:4];
         seqType='trained';
         runEffect='fixed';
         modelType='generic'; % or specific
         vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType'});
         
         T=[];
         for t=sessN
            R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_NEW_%s_%s_sess%d_%s.mat',modelType,seqType,t,runEffect)));
            R.sessN=ones(size(R.SN))*t;
            T=addstruct(T,R);
         end

        TT=T; 
        TT.bayesEst(:,[2:4])=bsxfun(@minus,TT.bayesEst(:,[2:4]),TT.bayesEst(:,2));
        TT2.modelInd=[ones(size(TT.SN));ones(size(TT.SN))*2;ones(size(TT.SN))*3;ones(size(TT.SN))*4];
        TT2.bayesEst=TT.bayesEst(:);
        TT2.roi=[TT.roi;TT.roi;TT.roi;TT.roi];
        TT2.sessN=[TT.sessN;TT.sessN;TT.sessN;TT.sessN];
        
        
         
        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3;ones(size(T.SN))*4];
        T2.bayesEst=T.bayesEst(:);
        T2.roi=[T.roi;T.roi;T.roi;T.roi];
        T2.sessN=[T.sessN;T.sessN;T.sessN;T.sessN];
        % rearranging how the data structure is arranged
        figure
        for r=reg
            subplot(1,numel(reg),r)
            plt.line([TT2.sessN>3 TT2.sessN],TT2.bayesEst,'subset',TT2.modelInd>1&TT2.roi==r,'split',TT2.modelInd,'leg',{'repsup','repsup+seq','repsup+seq+corr'},'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel(sprintf('Log-Bayes %s',seqType));
            else
                ylabel('');
            end
        end
    case 'PLOT_pcm_corr_allSess_NEW'
        reg = [1:8];
        sessN=[1:4];
        seqType={'trained','untrained'};
        runEffect='fixed';
        modelType='generic';
        fitM=1; % 1 - no added shared pattern; 2 - with added shared pattern
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType','fitM'});
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
               % R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_sess%d_%s_%s.mat',t,seqType{st},modelType)));
                R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_NEW_%s_%s_sess%d_%s.mat',modelType,seqType{st},t,runEffect)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        corrType={'naive','crossval','pcm'};
        corrVar=[T.r_naive, T.r_crossval, T.r_model2];
        for ct=1:3 % three types of correlation 
            figure(ct)
            sub=1;
            for r=reg
                subplot(1,numel(reg),sub)
                plt.line([T.sessN>3 T.sessN],corrVar(:,ct),'subset',T.roi==r,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
                plt.match('y');
                drawline(0,'dir','horz');
                title(sprintf('%s',regname{r}));
                if sub==1
                    ylabel(sprintf('%s corr',corrType{ct}));
                    xlabel('Session');
                else
                    ylabel('');
                end
                sub=sub+1;
            end     
        end
    case 'PLOT_pcm_corr_acrSess'
        sessN=[1:3];
        seqType='trained';
        regExcl=6;
        modelType='specific';
        runEffect='fixed';
        metric='r_model2'; % r_model2 or r_crossval or r_naive
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice','modelType','metric'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        
        % load structure, concatenate across sessions / seqType
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                % R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_sess%d_%s_%s.mat',t,seqType{st},modelType)));
                R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_NEW_%s_%s_sess%d_%s.mat',modelType,seqType,t,runEffect)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
                
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        switch(seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        figure
        plt.line(T.reg,T.(metric),'split',T.sessN,'subset',T.seqType==st,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        drawline(0,'dir','horz');
        ylabel('Distance repetition suppression');
        set(gca,'XTickLabel',regname_cortex(reg));
        title(sprintf('%s',seqType));
        xlabel('ROI');  
        
        
    case 'CORR_splithalf' % ------------------------------- NEW : SPLIT-HALF Correlation
        sn=[1:9,11:25];
        sessN=[1:4];
        betaChoice='multiPW';
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice'});
        
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
                        for j=1:2
                            for st=1:2 %seqType
                                T1 = getrow(T,T.SN==s & T.region==r);
                                % get betas
                                switch(betaChoice)
                                    case 'multiPW'
                                        betas = T1.betaW{:};
                                    case 'uniPW'
                                        betas = T1.betaUW{:};
                                    case 'raw'
                                        betas = T1.betaRAW{:};
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
                                    [K,G] = splitHalfCorr(betas1,D1.run,D1.seqNumb,'withinSes');
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
        save(fullfile(repSupDir,sprintf('corr_splitHalf_NEW_exe_%s',betaChoice)),'-struct','AllCorr');
    case 'PLOT_splithalf_corr'
        % plotting function for split-half correlation - within or across
        % repetition (for trained / untrained)
        roi=[1:8];
        betaChoice='multiPW';
        metric='corr'; % corr, corr_noMean, corr_mean
        vararginoptions(varargin,{'roi','metric','betaChoice'});
        T = load(fullfile(repSupDir,sprintf('corr_splitHalf_NEW_exe_%s',betaChoice)));
        
        for r=roi
            R = getrow(T,T.regNum==r);
            figure
            subplot(1,3,1)
            plt.line([R.sessN>3 R.sessN],R.(metric),'split',R.seqType,'subset',R.exe1~=R.exe2,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
            title(sprintf('%s - across exe',regname{r}));
            subplot(1,3,2)
            plt.line([R.sessN>3 R.sessN],R.(metric),'split',R.seqType,'subset',R.exe1==1 & R.exe2==1,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
            title(sprintf('%s - within exe1',regname{r}));
            subplot(1,3,3)
            plt.line([R.sessN>3 R.sessN],R.(metric),'split',R.seqType,'subset',R.exe1==2 & R.exe2==2,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
            title(sprintf('%s - within exe2',regname{r}));
            plt.match('y');

        end
        % split into within and across
    case 'CORR_splithalf_corrected'
        exe = [1 2];
        betaChoice='multiPW';
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice'});
        CC=[];
        T = load(fullfile(repSupDir,sprintf('corr_splitHalf_exe_%s',betaChoice)));
        sn=unique(T.sn);
        for ss=sessN
            for s=1:length(sn)
                for r=1:max(T.regNum)
                    D = getrow(T,T.sn==s & T.regNum==r & T.sessN==ss);
                    for st=1:2 % seqType
                        for i=1:length(exe)
                            for j=i:length(exe)
                                t1=getrow(D,D.seqType==st & D.exe1==i & D.exe2==i);
                                t2=getrow(D,D.seqType==st & D.exe1==j & D.exe2==j);
                                tcross=getrow(D,D.seqType==st & D.exe1==i & D.exe2==j);
                                if size(tcross.corr,1)~=0 % if there is data for both sessions
                                    C.corr = tcross.corr_noMean/sqrt(t1.corr_noMean*t2.corr_noMean);
                                    if ~isreal(C.corr)
                                        C.corr=0; % or NaN
                                    end
                                    % other info
                                    C.exe1=i;
                                    C.exe2=j;
                                    C.sessN=ss;
                                    C.seqType=st;
                                    C.regNum=r;
                                    if r<9
                                        C.regSide=1;
                                        C.regType=r;
                                    else
                                        C.regSide=2;
                                        C.regType=r-8;
                                    end
                                    CC=addstruct(CC,C);
                                end
                            end
                        end
                    end; % seqtype
                end; % reg
                fprintf('%d \t done sess%d - %s\n',s,ss,subj_name{s});
            end; % sn
        end; % sessN
        save(fullfile(repSupDir,sprintf('Corr_splitHalf_exe_%s_corrected',betaChoice)),'-struct','CC');
    case 'PLOT_splithalf_corr_corrected'
        roi=[1:8];
        betaChoice='multiPW';
        vararginoptions(varargin,{'roi','metric','betaChoice'});
        T = load(fullfile(repSupDir,sprintf('Corr_splitHalf_exe_%s_corrected',betaChoice)));
        
        figure
        for r=roi
            subplot(1,length(roi),r);
            R = getrow(T,T.regNum==r);
            plt.line([R.sessN>3 R.sessN],R.corr,'split',R.seqType,'subset',R.exe1~=R.exe2,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
            title(sprintf('%s - across exe',regname{r}));  
        end  
    case 'CORR_splithalf_structure'
        betaChoice='multiPW';
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice'});
        CC=[];
        T = load(fullfile(repSupDir,sprintf('corr_splitHalf_exe_%s',betaChoice)));
        sn=unique(T.sn);
        
        for ss=sessN
            for s=sn'
                for r=1:max(T.regNum)
                    for st=1:2 % seqType
                        K=getrow(T,T.sn==s&T.sessN==ss&T.regNum==r&T.seqType==st);
                        R(1) = K.G(1);
                        R(2) = K.G(4);
                        C.corr      = rsa_calcCorrRDMs(R);
                        C.distCorr  = rsa_calcDistCorrRDMs(R);
                        C.regNum    = r;
                        C.seqType   = st;
                        C.sn        = s;
                        C.sessN     = ss;
                        C.regType   = K.regType(1);
                        C.regSide   = K.regSide(1);
                        CC = addstruct(CC,C);
                    end
                end
                fprintf('%d \t done sess%d - %s\n',s,ss,subj_name{s});
            end
        end
        save(fullfile(repSupDir,sprintf('Corr_Gstruct_exe_%s',betaChoice)),'-struct','CC');
    case 'PLOT_splithalf_structure'
        roi=[1:8];
        betaChoice='multiPW';
        metric='corr';
        vararginoptions(varargin,{'roi','metric','betaChoice'});
        T = load(fullfile(repSupDir,sprintf('Corr_Gstruct_exe_%s',betaChoice)));
        
        T=normData(T,metric);
        metricNew=sprintf('norm%s',metric);
        figure
        for r=roi
            subplot(1,length(roi),r);
            R = getrow(T,T.regNum==r);
            plt.line([R.sessN>3 R.sessN],R.(metricNew),'split',R.seqType,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
            if r==1
                ylabel('Correlation');
            else
                ylabel('');
            end
            plt.match('y');
            title(sprintf('%s - across exe',regname{r}));  
        end  
        
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
function M = pcm_corrModel_sharedSeqType
% --------------------------------------
% Models separately trained / untrained sequences
% stability of representation 1st -> 2nd execution

        % Model1: Model with independent 1st / 2nd patterns
        % across sessions - 0 correlation
        
        % Model 1 - independent sequence patterns between sessions
        M{1}.type          = 'feature';
        M{1}.name          = 'ind';
        M{1}.numGparams    = 4;
        A=zeros(6);
        for i=1:6
            A(i,i)=1;
        end;
        M{1}.Ac(:,1:6 ,1)  = [A;zeros(6)];     % Unique sess1 sequence patterns    (theta_a) (either trained or untrained)
        M{1}.Ac(:,7:12,2)  = [zeros(6);A];     % Unique sess2 sequence pattterns   (theta_b)
        M{1}.Ac(:,13,3)    = [ones(6,1); zeros(6,1)]; % same pattern across all first exe
        M{1}.Ac(:,14,4)    = [zeros(6,1); ones(6,1)]; % same pattern second exe
        
        % --------------------------------------
        % Model2: Model with a flexible across-sess correlation for sequences 
        M{2}.type          = 'feature';
        M{2}.numGparams    = 5;
        M{2}.name          = 'flex';
        M{2}.Ac(:,1:6 ,1)  = [A;zeros(6)];     % Unique sess1 sequence patterns    (theta_a) (either trained or untrained)
        M{2}.Ac(:,7:12,2)  = [zeros(6);A];     % Unique sess2 sequence pattterns   (theta_b)
        M{2}.Ac(:,1:6,3)   = [zeros(6);A];     % Same sess1 patterns   (theta_c)
        M{2}.Ac(:,13,4)    = [ones(6,1); zeros(6,1)]; % same pattern across all first exe
        M{2}.Ac(:,14,5)    = [zeros(6,1); ones(6,1)]; % same pattern second exe
        
        % --------------------------------------
        % Model3: Model with a fixed r=1 correlation (second session seq same as first)
        M{3}.type         = 'feature';
        M{3}.numGparams   = 4;
        M{3}.name         = 'one';
        M{3}.Ac(:,1:6,1)  = [A;zeros(6)];     % Unique sess1 patterns      (theta_a)
        M{3}.Ac(:,1:6,2)  = [zeros(6);A];     % Same sess2 pattterns       (theta_b)
        M{3}.Ac(:,7,3)    = [ones(6,1); zeros(6,1)]; % same pattern across all first exe
        M{3}.Ac(:,8,4)    = [zeros(6,1); ones(6,1)]; % same pattern second exe
        
            
end
function M = pcm_corrModel_indSeq_sharedSeqType
% --------------------------------------
% Models separately for each first / second sequences
% including also common pattern for first / second exe

    % Model1: Model with independent trained / untrained patterns
    % across sessions - 0 correlation

    % Model 1 - independent sequence patterns between sessions
    M{1}.type       = 'feature';
    M{1}.name       = 'ind';
    M{1}.numGparams = 14;

    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{1}.Ac(:,1:6 ,i)    = [A;zeros(6)];       % Unique seq1 patterns   (theta_a)
        M{1}.Ac(:,7:12,6+i) = [zeros(6);A];        % Unique seq2 pattterns   (theta_c)
    end;
    
    M{1}.Ac(:,13,13) = [ones(6,1);zeros(6,1)];
    M{1}.Ac(:,14,14) = [zeros(6,1);ones(6,1)];
    
    
    % --------------------------------------
    % Model2: Model with a flexible across-sess correlation for sequences
    M{2}.type       = 'feature';
    M{2}.numGparams = 20;
    M{2}.name       = 'flex';

    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{2}.Ac(:,1:6 ,i)    = [A;zeros(6)];       % Seq1 patterns   (theta_a)
        M{2}.Ac(:,7:12,6+i) = [zeros(6);A];       % Unique seq2 pattterns   (theta_b)
        M{2}.Ac(:,1:6 ,12+i)  = [zeros(6);A];       % Same seq2 patterns  (theta_c)
    end;
    M{2}.Ac(:,13,19) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,14,20) = [zeros(6,1);ones(6,1)];

    % --------------------------------------
    % Model3: Model with a fixed r=1 correlation (second session seq same as first)
    M{3}.type       = 'feature';
    M{3}.numGparams = 14;
    M{3}.name        = 'one';

    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,1:6 ,i)    = [A;zeros(6)]; % Seq1 finger patterns   (theta_a)
        M{3}.Ac(:,1:6 ,6+i)  = [zeros(6);A]; % Same seq2 patterns  (theta_b)
    end;
    M{3}.Ac(:,7,13) = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,8,14) = [zeros(6,1);ones(6,1)];
end
function M = pcm_repsupModel_generic

    % for sequence-specific modelling - one parameter for all seq           
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    
    % Model 1: No sequence pattern
    M{1}.type           = 'feature';
    M{1}.numGparams     = 1;
    M{1}.name           = 'null';
    M{1}.Ac(:,1:12,1)   = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type           = 'feature';
    M{2}.numGparams     = 2;
    M{2}.name           = 'RepSup';
    M{2}.Ac(:,1,1)      = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2)      = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type           = 'feature';
    M{3}.numGparams     = 4;
    M{3}.name           = 'RepSup+Seq';
    M{3}.Ac(:,1:6,1)    = [A;zeros(6)];             % Unique exe1 sequence patterns
    M{3}.Ac(:,7:12,2)   = [zeros(6);A];             % Unique exe2 sequence pattterns
    M{3}.Ac(:,13,3)     = [ones(6,1);zeros(6,1)];   % Exe1 mean
    M{3}.Ac(:,14,4)     = [zeros(6,1);ones(6,1)];   % Exe2 mean
    
    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 5;
    M{4}.name         = 'RepSup+Seq+Corr';
    M{4}.Ac(:,1:6,1)  = [A;zeros(6)];           % Unique exe1 sequence patterns
    M{4}.Ac(:,7:12,2) = [zeros(6);A];           % Unique exe2 sequence pattterns
    M{4}.Ac(:,1:6,3)  = [zeros(6);A];           % Correlation exe1-exe2
    M{4}.Ac(:,13,3)   = [ones(6,1);zeros(6,1)]; % Exe1 mean
    M{4}.Ac(:,14,4)   = [zeros(6,1);ones(6,1)]; % Exe2 mean
    
end
function M = pcm_repsupModel_specific
    
    % Model 1: No sequence pattern
    M{1}.type           = 'feature';
    M{1}.numGparams     = 1;
    M{1}.name           = 'null';
    M{1}.Ac(:,1:12,1)   = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'RepSup';
    M{2}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 14;
    M{3}.name       = 'RepSup+Seq';
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,1:6,i)     = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{3}.Ac(:,7:12,6+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    end;
    M{3}.Ac(:,13,13)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,14,14)  = [zeros(6,1);ones(6,1)];
    
    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 20;
    M{4}.name         = 'RepSup+Seq+Corr';
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,1:6,i)     = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{4}.Ac(:,7:12,6+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
        M{4}.Ac(:,1:6,12+i)  = [zeros(6);A];     % Correlation exe1-exe2
    end;
    M{4}.Ac(:,13,19)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,14,20)    = [zeros(6,1);ones(6,1)];

end

function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Crossvalidated model comparision:
    [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    [Tcross2,thetaCr2] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);

    T.cross_likelihood = Tcross.likelihood;
    T.bayesEst = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
    T.theta_hat=theta_hat;
    for t=1:size(thetaCr,2)
        tC{t} = thetaCr{t}';
    end
    T.thetaCr=tC;
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
% 1. Empirical correlation
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec{p});
    b = pinv(Z)*Data{p};           % Estimate mean activities
    G=cov(b');
    C.r_naive_wMean(p,1) = calcCorr(G);
    b(1:6,:)  = bsxfun(@minus,b(1:6,:) ,mean(b(1:6,:))); % Subtract mean per condition - first exe
    b(7:12,:) = bsxfun(@minus,b(7:12,:),mean(b(7:12,:))); % second exe
    G=cov(b');
    C.r_naive(p,1) = calcCorr(G);
end;
% --------------------------------------
% 2. Crossvalidated correlation - make PD
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec{p});
    % Subtract mean for each condition and run
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>6)-1);
    R = eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)         = pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    Gcv2(:,:,p)        = crossval_estG(Data{p},Z,condVec{p});
    C.r_crossval(p,1)  = calcCorr(pcm_makePD(Gcv(:,:,p)));
    C.r_crossval2(p,1) = calcCorr(pcm_makePD(Gcv2(:,:,p)));
    
    A=[];
    % new - different way of subtracting mean for each run - 1st / 2nd
    for i=1:numel(unique(partVec{p}))
        Xa = Z(partVec{p}==i,:);
        Ya = Data{p}(partVec{p}==i,:);
        A(:,:,i) = pinv(Xa)*Ya;
        %  subtract the mean across conditions - 1st / 2nd
        A(1:6,:,i) = bsxfun(@minus,A(1:6,:,i),mean(A(1:6,:,i)));
        A(7:12,:,i) = bsxfun(@minus,A(7:12,:,i),mean(A(7:12,:,i)));
    end;
    % reshape data back
    data=[];
    for i=1:numel(unique(partVec{p}))
        data = [data; A(:,:,i)];
    end
    Gcv3(:,:,p) = crossval_estG(data,Z,condVec{p});
    b = pinv(Z)*data;
    G2 = cov(b');
    C.r_naive2(p,1) = calcCorr(G2);
    C.r_crossval3(p,1) = calcCorr(pcm_makePD(Gcv3(:,:,p)));
end;

% --------------------------------------
% 3. Fit model 2  and infer correlations from the parameters
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);
C.theta=theta;
% Get the correlations
switch M_type
    case 1
        var1       = (theta{1}(1,:).^2)';
        var2       = (theta{1}(2,:).^2+theta{1}(3,:).^2)';
        cov12      = (theta{1}(1,:).*theta{1}(3,:))';
        C.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
    case 2
        var1       = (theta{1}(1:6,:).^2)';
        var2       = (theta{1}(7:12,:).^2+theta{1}(13:18,:).^2)';
        cov12      = (theta{1}(1:6,:).*theta{1}(13:18,:))';
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
function [R,G] = splitHalfCorr(data,partVec,condVec,type)
% function R = splitHalfCorr(data,partVec,condVec,type)
% performs split-half correlations
switch(type)
    case 'withinSes'
        % calculate within session split-half correlation
        X = indicatorMatrix('identity_p',condVec);
        G = crossval_estG(data,X,partVec);
        mPat=[];
        split=mod(partVec,2);
        % new
        for i=1:numel(unique(partVec))
            Xa = X(partVec==i,:);
            Ya = data(partVec==i,:);
            A(:,:,i) = pinv(Xa)*Ya;
        end;
      %  subtract the mean across conditions
        A=bsxfun(@minus,A,sum(A,1)/numel(unique(condVec)));
        % reshape data back
        data=[];
        for i=1:numel(unique(partVec))
            data = [data; A(:,:,i)];
        end
        % end of new
        mPat(:,:,1)     = pinv(X(split==0,:))*data(split==0,:);
        mPat(:,:,2)     = pinv(X(split==1,:))*data(split==1,:);
        COR             = corr(mPat(:,:,1)',mPat(:,:,2)');
        R.corr          = fisherinv(mean(fisherz(diag(COR))));
        meanPat         = mean(mPat,1);             % Mean activity pattern over sequences
        mPat            = bsxfun(@minus,mPat,meanPat);  % Subtract out mean
        COR             = corr(mPat(:,:,1)',mPat(:,:,2)');
        R.corr_noMean   = fisherinv(mean(fisherz(diag(COR))));
        R.corr_mean     = corr(meanPat(:,:,1)',meanPat(:,:,2)');
        
    case 'acrossSes'
        Y     = [data{1};data{2}];
        part  = [partVec{1};partVec{2}];
        cond  = [condVec{1};condVec{2}];
        ses   = [ones(size(partVec{1}));ones(size(partVec{2}))*2];
        X     = indicatorMatrix('identity_p',cond);
        G = crossval_estG(Y,X,ses);
        mPat=[];  
        % new
        for p=1:size(partVec,2)
            for i=1:numel(unique(partVec{p}))
                Xa = X(partVec{p}==i,:);
                Ya = data{p}(partVec{p}==i,:);
                A(:,:,i) = pinv(Xa)*Ya;
            end;
            %  subtract the mean across conditions
            A=bsxfun(@minus,A,sum(A,1)/numel(unique(condVec{p})));
            % reshape data back
            data{p}=[];
            for i=1:numel(unique(partVec{p}))
                data{p} = [data{p}; A(:,:,i)];
            end
        end
        % end of new
        % Do split half correlations
        split = mod(part,2);
        mPat(:,:,1) = pinv(X(ses==1 & split==0,:)) * data{1}(split(ses==1)==0,:);
        mPat(:,:,2) = pinv(X(ses==1 & split==1,:)) * data{1}(split(ses==1)==1,:);
        mPat(:,:,3) = pinv(X(ses==2 & split==0,:)) * data{2}(split(ses==2)==0,:);
        mPat(:,:,4) = pinv(X(ses==2 & split==1,:)) * data{2}(split(ses==2)==1,:);
        COR=[];
        COR(:,:,1)  = corr(mPat(:,:,1)',mPat(:,:,3)');
        COR(:,:,2)  = corr(mPat(:,:,1)',mPat(:,:,4)');
        COR(:,:,3)  = corr(mPat(:,:,2)',mPat(:,:,3)');
        COR(:,:,4)  = corr(mPat(:,:,2)',mPat(:,:,4)');
        R.corr      = fisherinv(mean(diag(mean(fisherz(COR),3))));
        meanPat     = mean(mPat,1);             % Mean activity pattern over fingers
        mPat        = bsxfun(@minus,mPat,meanPat);  % Subtract out mean
        COR=[];
        COR(:,:,1)  = corr(mPat(:,:,1)',mPat(:,:,3)');
        COR(:,:,2)  = corr(mPat(:,:,1)',mPat(:,:,4)');
        COR(:,:,3)  = corr(mPat(:,:,2)',mPat(:,:,3)');
        COR(:,:,4)  = corr(mPat(:,:,2)',mPat(:,:,4)');
        R.corr_noMean = fisherinv(mean(diag(mean(fisherz(COR),3))));
        COR=[];
        COR(1) = corr(meanPat(:,:,1)',meanPat(:,:,3)');
        COR(2) = corr(meanPat(:,:,1)',meanPat(:,:,4)');
        COR(3) = corr(meanPat(:,:,2)',meanPat(:,:,3)');
        COR(4) = corr(meanPat(:,:,2)',meanPat(:,:,4)');
        R.corr_mean = fisherinv(mean(fisherz(COR)));
end
end
