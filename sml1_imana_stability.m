function varargout=sml1_imana_stability(what,varargin)

% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
codeDir         ='/Users/eberlot/Documents/MATLAB/projects/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
behavDir        =[baseDir '/behavioral_data/analyze'];
imagingDir      =[baseDir '/imaging_data'];                   
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia'];
suitDir         =[baseDir '/suit'];
pcmDir          =[baseDir '/pcm_stats'];
stabDir         =[baseDir '/stability'];
qualContrDir    =[baseDir '/qual_control'];
distPscDir      =[baseDir '/dist_psc_stats'];

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
subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18',...
                's19','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'};  

% Other random notes %

% -------------------------- For plotting ---------------------------------
stySeq=style.custom({'red','blue'},'markersize',12);
stySeq2=style.custom({'red','blue'},'linestyle','--','markersize',12);

style.file(fullfile(codeDir,'sml_style.m'));
style.use('default');

% ------------------------------ Analysis Cases --------------------------------
switch(what)
 
    case 'beta_consistPattern_wSess'  
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % enter sn, region,  beta: 0=betaW, 1=betaU, 2=raw betas
        % (1) Set parameters
        sn  = [4:9,11:31];
        roi = 1:8;
        sessN=1:4;
        betaChoice='multiPW'; % multiPW, uniPW, raw
        parcelType='Brodmann'; 
        seqType = 'overall'; % perform calculations 'overall' or 'separate'
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice','parcelType','seqType'});
        
        numRun=numruns_task_sess;
        numCond=numel(num_seq);
        
        RR=[];
        %========%
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss))); % loads in struct 'T'
            for h=1:2
                for r=roi
                    for s=sn
                        S = getrow(T,(T.SN==s & T.regType==r & T.regSide==h));
                        switch(betaChoice)
                            case 'raw'
                                beta  = S.betaRAW{1};
                            case 'uniPW'
                                beta  = S.betaUW{1};
                            case 'multiPW'
                                beta  = S.betaW{1};
                        end
                         res = S.resMS{1};
                        % make vectors for pattern consistency func
                        partVec = kron([1:numRun]',ones(numCond,1));
                        condVec = kron(ones(numRun,1),[1:numCond]');
                        % here split by seqType
                        switch seqType
                            case 'overall'           
                                % calculate the pattern consistency
                                R.cnr                           = cnr_QC(beta,res,numCond,numRun);  % contrast to noise ratio
                                R.r2_rm                         = rsa_patternConsistency(beta,partVec,condVec,'removeMean',1);  % pattern consistency
                                R.r2                            = rsa_patternConsistency(beta,partVec,condVec,'removeMean',0);
                                [R.r2_cross_rm, R.r_cross_rm]   = bsp_patternConsistency_crossval(beta,partVec,condVec,'removeMean',1); % correlation of patterns
                                [R.r2_cross, R.r_cross]         = bsp_patternConsistency_crossval(beta,partVec,condVec,'removeMean',0);        
                                R.sn        = s;
                                R.sessN     = ss;
                                R.region    = (h-1)*max(roi)+r;
                                R.regType   = r;
                                R.regSide   = h;    
                                RR=addstruct(RR,R);
                            case 'separate'
                                for st=1:2
                                    idx = ismember(condVec,[(st-1)*6+1:(st-1)*6+6]);
                                    R.cnr                           = cnr_QC(beta(idx,:),res,numCond,numRun);  % contrast to noise ratio
                                    R.r2_rm                         = rsa_patternConsistency(beta(idx,:),partVec(idx,:),condVec(idx,:),'removeMean',1);  % pattern consistency
                                    R.r2                            = rsa_patternConsistency(beta(idx,:),partVec(idx,:),condVec(idx,:),'removeMean',0);
                                    [R.r2_cross_rm, R.r_cross_rm]   = bsp_patternConsistency_crossval(beta(idx,:),partVec(idx,:),condVec(idx,:),'removeMean',1); % correlation of patterns
                                    [R.r2_cross, R.r_cross]         = bsp_patternConsistency_crossval(beta(idx,:),partVec(idx,:),condVec(idx,:),'removeMean',0);
                                    R.sn        = s;
                                    R.sessN     = ss;
                                    R.region    = (h-1)*max(roi)+r;
                                    R.regType   = r;
                                    R.regSide   = h;
                                    R.seqType   = st;
                                    RR=addstruct(RR,R);
                                end
                        end
                    end
                    fprintf('Done sess-%d hemi-%d roi-%d\n',ss,h,r);
                end
            end
        end   
        %save the structure
        save(fullfile(stabDir,sprintf('consist_wSess_%s_%sBetas_%s',parcelType,betaChoice,seqType)),'-struct','RR');
        case 'PLOT_pattern_consist' 
       betaChoice = 'multiPW';
       parcelType='Brodmann';
       hemi=1;
       vararginoptions(varargin,{'var','betaChoice','roi','parcelType','hemi'});
       
       T = load(fullfile(stabDir,sprintf('consist_wSess_%s_%sBetas_overall',parcelType,betaChoice)));
       
       style.use('Sess');
       figure
       subplot(211)
       plt.bar(T.regType,T.r_cross,'subset',T.regSide==hemi,'split',T.sessN,'leg',{'sess1','sess2','sess3','sess4'},...
           'leglocation','northeast');
       title('crossval-correlation');
       ylabel('r');
       subplot(212)
       plt.bar(T.regType,T.r_cross_rm,'subset',T.regSide==hemi,'split',T.sessN,'leg',{'sess1','sess2','sess3','sess4'},...
           'leglocation','northeast');
       title('crossval-correlation mean removed');
       ylabel('r');    
        case 'PLOT_pattern_consist_seqType'
       betaChoice = 'multiPW';
       parcelType='Brodmann';
       hemi=1;
       roi=[1:8];
       vararginoptions(varargin,{'var','betaChoice','roi','parcelType','hemi'});
       
       T = load(fullfile(stabDir,sprintf('consist_wSess_%s_%sBetas_separate',parcelType,betaChoice)));
       style.use('Seq');
       figure
       for r=1:length(roi)
           subplot(2,numel(roi),r)
           plt.bar(T.sessN,T.r_cross,'subset',T.regSide==hemi & T.regType==roi(r),'split',T.seqType,'leg',{'trained','untrained'},...
               'leglocation','northeast');
           title(sprintf('cross-corr %s',regname_cortex{roi(r)}));
           ylabel('r');
           subplot(2,numel(roi),numel(roi)+r)
           plt.bar(T.sessN,T.r_cross_rm,'subset',T.regSide==hemi & T.regType==roi(r),'split',T.seqType,'leg',{'trained','untrained'},...
               'leglocation','northeast');
           title(sprintf('noMean %s',regname_cortex{roi(r)}));
           ylabel('r');
       end    
        case 'STATS_pattern_consist'
            betaChoice = 'multiPW';
            parcelType='Brodmann';
            hemi=1;
            roi=[1:8];
            vararginoptions(varargin,{'var','betaChoice','roi','parcelType','hemi'});
            
            T = load(fullfile(stabDir,sprintf('consist_wSess_%s_%sBetas_overall',parcelType,betaChoice)));
            for r=roi
                fprintf('\n session ANOVA for %s : overall pattern\n',regname_cortex{r});
                anovaMixed(T.r_cross,T.sn,'within',[T.sessN],{'session'},'subset',T.regType==r & T.regSide==hemi);
                fprintf('mean pattern removed\n')
                anovaMixed(T.r_cross_rm,T.sn,'within',[T.sessN],{'session'},'subset',T.regType==r & T.regSide==hemi);
            end
            keyboard;
            % posthoc t-tests
            ind = indicatorMatrix('allpairs',[1:4]);
            for r=roi
                for i=1:size(ind,1)
                    ind_r = find(ind(i,:));
                    fprintf('\n t-test %s : overall pattern, sess:%d-%d\n',regname_cortex{r},ind_r(1),ind_r(2));
                    ttestDirect(T.r_cross,[T.sessN T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & ismember(T.sessN,ind_r));
                    fprintf('mean pattern removed\n')
                    ttestDirect(T.r_cross_rm,[T.sessN T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & ismember(T.sessN,ind_r));
                end
            end
        case 'STATS_pattern_consist_seqType'
            betaChoice = 'multiPW';
            parcelType='Brodmann';
            hemi=1;
            roi=[1:8];
            vararginoptions(varargin,{'var','betaChoice','roi','parcelType','hemi'});
            
            T = load(fullfile(stabDir,sprintf('consist_wSess_%s_%sBetas_separate',parcelType,betaChoice)));
            for r=roi
                fprintf('\n session ANOVA for %s : overall pattern\n',regname_cortex{r});
                anovaMixed(T.r_cross,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi);
                fprintf('mean pattern removed\n')
                anovaMixed(T.r_cross_rm,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi);
            end
            keyboard;
            % posthoc t-tests
            for r=roi
                for ss=1:4
                    fprintf('\n t-test %s : overall pattern, sess:%d\n',regname_cortex{r},ss);
                    ttestDirect(T.r_cross,[T.seqType T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & T.sessN==ss);
                    fprintf('mean pattern removed\n')
                    ttestDirect(T.r_cross_rm,[T.seqType T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & T.sessN==ss);
                end
            end
    case 'beta_runwiseDrift'
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % enter sn, region,  beta: 0=betaW, 1=betaU, 2=raw betas
        % (1) Set parameters
        sn  = [4:9,11:31];
        roi = [1:8];
        sessN=[1:4];
        betaChoice='multiPW'; % multiPW, uniPW, raw
        parcelType='Brodmann'; 
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice','parcelType'});
        
        numRun=numruns_task_sess;
        numCond=numel(num_seq);
        % make vectors for pattern consistency func
        partVec = kron([1:numRun]',ones(numCond,1));
        %condVec = kron(ones(numRun,1),[1:numCond]');
        RR=[];
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss))); % loads in struct 'T'
            for h=1:2
                for r=roi
                    for s=sn
                        S = getrow(T,(T.SN==s & T.regType==r & T.regSide==h));
                        switch(betaChoice)
                            case 'raw'
                                beta  = S.betaRAW{1};
                            case 'uniPW'
                                beta  = S.betaUW{1};
                            case 'multiPW'
                                beta  = S.betaW{1};
                        end
                        for r1=1:8
                            for r2=r1:8
                                beta1       = beta(partVec==r1,:);
                                beta2       = beta(partVec==r2,:);
                                R.corr      = mean(diag(corr(beta1',beta2')));
                                beta1       = bsxfun(@minus,beta1,mean(beta1,1));
                                beta2       = bsxfun(@minus,beta2,mean(beta2,1));
                                R.corr_rm   = mean(diag(corr(beta1',beta2')));
                                R.run1      = r1;
                                R.run2      = r2;
                                R.sn        = s;
                                R.sessN     = ss;
                                R.region    = (h-1)*max(roi)+r;
                                R.regType   = r;
                                R.regSide   = h;
                                RR=addstruct(RR,R);
                            end
                        end
                    end
                    fprintf('Done sess-%d hemi-%d roi-%d\n',ss,h,r);
                end
            end
        end
        %save the structure
        save(fullfile(stabDir,sprintf('consist_wSess_runwise_%s_%sBetas',parcelType,betaChoice)),'-struct','RR');    
        case 'PLOT_runwiseDrift'
        betaChoice='multiPW'; % multiPW, uniPW, raw
        parcelType='Brodmann'; 
        roi=[1:8];
        hemi=1;
        var='corr'; % corr or corr_rm
        vararginoptions(varargin,{'roi','betaChoice','parcelType','hemi','var'});
        T = load(fullfile(stabDir,sprintf('consist_wSess_runwise_%s_%sBetas',parcelType,betaChoice)));
        
        figure
        for r=1:numel(roi)
            subplot(numel(roi),1,r)
            plt.line(T.run1,T.(var),'subset',T.regType==roi(r) & T.regSide==hemi & T.run2==T.run1+1,'split',T.sessN,...
                'style',stySess);
            title(regname_cortex{roi(r)});
        end
        
        figure
        for r=1:numel(roi)
            subplot(numel(roi),1,r)
            plt.line(T.run2,T.(var),'subset',T.regType==roi(r) & T.regSide==hemi & T.run1==1 & T.run2~=1,'split',T.sessN,...
                'style',stySess);
            title(regname_cortex{roi(r)});
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
        reg = [1:16]; % across hemispheres 
        sn  = [4:9,11:31];
        sessN = [1:2];
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        parcelType='Brodmann';
        % options for regularisation 'reg'
        % - minvalue (if variance <0.001 -> 0.001
        % - zero (if variance negative, make it 0)
        
        vararginoptions(varargin,{'sn','reg','sessN','regcorrType','seqTypeCorr','parcelType'});
        CAll = [];
        D={}; data={};
        for ss=1:numel(sessN) % pre-allocate the betas for the whole group
            D{ss} = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
        end
        for s = sn;
            fprintf('\nSubject %d/%d:',find(s==sn),numel(sn));
            for roi = reg;
                for  ss = 1:numel(sessN)
                    t   = getrow(D{ss},D{ss}.region==roi & D{ss}.SN==s); % get subject-specific betas
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
                C.region  = ones(size(C.w1)).*roi;
                C.regType = ones(size(C.w1)).*t.regType;
                C.regSide = ones(size(C.w1)).*t.regSide;
                C.sn = ones(size(C.w1)).*s;
                CAll=addstruct(CAll,C);
                fprintf('%d.',roi);
            end
        end
        save(fullfile(stabDir,sprintf('stability_acrSess_%s_sess%d-%d_%s_%sRegularised',seqTypeCorr,sessN(1),sessN(2),parcelType,regcorrType)),'-struct','CAll');
        case 'reliability_simulate'
        % create different data types for simulating correlations within /
        % between sessions
        % ensure that cross-session stability working properly
        % geometric mean > across-sess correlation
        type = 'corr+noise';
        regcorrType = 'zero';
        vararginoptions(varargin,{'type','run_type','regcorrType'});
        P = 100;    % number of voxels simulated
        cond = 12;  % number of conditions
        run=2;
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
                noise1  = randn(cond*run,P).*noise_lev;
                noise2  = randn(cond*run,P).*noise_lev;
                data{1} = true_X + noise1;
                data{2} = true_X + noise2;
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
            C = sml1_imana_stability('reliability_seqcorr','data',data,'removeMean',1,'numruns',run,'reg',regcorrType);
            C.perm = ones(size(C.w1)).*perm;
            CAll=addstruct(CAll,C);
        end
        CAll.dif = CAll.acr - CAll.geoMean;
        % CAll.dif - positive when rel across sess bigger than geometric
        % mean of reliability of each session
        figure;
        subplot(121)
        histogram(CAll.dif,'Normalization','probability');
        title('Difference: across-session - geoMean')
        subplot(122)
        barplot(CAll.seqType,[CAll.w1 CAll.w2 CAll.acr CAll.geoMean]);
        case 'reliability_seqcorr'
        
        vararginoptions(varargin,{'data','removeMean','numruns','reg'});
        
        sess = size(data,2);
        partitions = [1:2:numruns; 2:2:numruns];
        numRuns    = 1:numruns;
        numConds   = num_seq;
        condVec   = repmat([numConds],1,length(numRuns))';
        partVec = kron([numRuns],ones(1,length(numConds)))';
        
        for ss = 1:sess         
            % subtract mean per run
            if removeMean
                for r=numRuns
                    data{ss}(partVec==r,:) = bsxfun(@minus,data{ss}(partVec==r,:),mean(data{ss}(partVec==r,:)));
                end
            else
                data{ss}=data{ss};
            end
            
            % split datas per partition
            for i = 1:size(partitions,1)
                partitionIdx = logical(ismember(partVec,partitions(i,:)))';
                condIdx{i}   = condVec(partitionIdx);
                prepBetas{ss}{i} = data{ss}(partitionIdx,:); % session / partition
            end
            
            % calculate average pattern per partition
            for c1 = numConds % for each condition
                % condition mean seq-specific activity pattern for this partition
                oddCon   = condIdx{1}==c1;
                oddBetas{ss}(c1,:) = mean(prepBetas{ss}{1}(oddCon,:),1);
                % condition mean activity pattern for the other partition
                evenCon   = condIdx{2}==c1;
                evenBetas{ss}(c1,:) = mean(prepBetas{ss}{2}(evenCon,:),1);
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
        cond_type = repmat([ones(6,1); ones(6,1)*2]',1,length(numRuns));
        partVec = kron([numRuns],ones(1,length(numConds)));
        
        % for constructing average patterns
        numConds_new = 1:num_seqtype; 
        cond_new = repmat([numConds_new],1,length(numRuns));
        runNums_new = kron([numRuns],ones(1,length(numConds_new)));
        
        for ss = 1:sess
            
            % calculate average trained / untrained pattern
            for r=numRuns
                for s=numConds_new
                data2{ss}(runNums_new==r & cond_new==s,:) = mean(data{ss}(partVec==r & cond_type==s,:),1); % 6 sequences
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
        sessN = [1:2];
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
        hemi = 1;
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        parcelType = 'Brodmann';
        sessTr=[1:3]; % which session transitions to include
        %1: sess1-2
        %2: sess2-3
        %3: sess3-4
        sessName={'1-2','2-3','3-4'};
        
        vararginoptions(varargin,{'reg','regcorrType','seqTypeCorr','sessTr','parcelType','hemi'});
        T=[];
        for t=1:numel(sessTr)
            R=load(fullfile(stabDir,sprintf('stability_acrSess_%s_sess%d-%d_%s_%sRegularised',seqTypeCorr,t,t+1,parcelType,regcorrType)));
            R.sessTr=ones(size(R.w1))*t;
            T=addstruct(T,R);
        end
        
        figure
        for r=1:numel(reg)
            
            subplot(3,numel(reg),r)
            %plt.line([T.sessTr>2 T.sessTr],T.acr,'split',T.seqType,'style',stySeq,'subset',T.roi==r,'leg',{'trained','untrained'});
            plt.line([T.sessTr>2 T.sessTr],T.w1,'split',T.seqType,'style',stySeq,'subset',T.regType==r & T.regSide==hemi,'leg',{'trained','untrained'});
            if r==1
                xlabel('Sess transitions');
                ylabel('geometric mean');
            else
                ylabel('');
            end
            set(gca,'XTickLabel',sessName(sessTr));
            drawline(0,'dir','horz');              
            title(sprintf('%s',regname{r}));
            subplot(3,numel(reg),r+numel(reg))
            plt.line([T.sessTr>2 T.sessTr],T.acr,'split',T.seqType,'style',stySeq,'subset',T.regType==r & T.regSide==hemi,'leg',{'trained','untrained'});
            if r==1
                xlabel('Sess transitions');
                ylabel(sprintf('%s correlation',seqTypeCorr));
            else
                ylabel('');
            end
            set(gca,'XTickLabel',sessName(sessTr));
            drawline(0,'dir','horz');
            subplot(3,numel(reg),r+2*numel(reg))
            plt.line([T.sessTr>2 T.sessTr],T.acrSess_corr,'split',T.seqType,'style',stySeq,'subset',T.regType==r & T.regSide==hemi,'leg',{'trained','untrained'});
            if r==1
                xlabel('Sess transitions');
                ylabel('corrected');
            else
                ylabel('');
            end
            set(gca,'XTickLabel',sessName(sessTr));
            drawline(0,'dir','horz');
        end
        case 'STATS_reliability'
        roi = [1:8];
        hemi = 1;
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        parcelType = 'Brodmann';
        var = 'acr'; % variable to assess: geoMean, acr, acrSess_corr
        sessTr=[1:2]; % which session transitions to include
        %1: sess1-2
        %2: sess2-3
        %3: sess3-4
        
        vararginoptions(varargin,{'reg','regcorrType','seqTypeCorr','sessTr','hemi','var'});
        T=[];
        for t=sessTr
            R=load(fullfile(stabDir,sprintf('stability_acrSess_%s_sess%d-%d_%s_%sRegularised',seqTypeCorr,t,t+1,parcelType,regcorrType)));
            R.sessTr=ones(size(R.w1))*t;
            T=addstruct(T,R);
        end
        
        % Anova - sessTransition x seqType - per roi
        for r=roi
            fprintf('\n sessTransition x seqType ANOVA for %s \n',regname{r});
            anovaMixed(T.(var),T.sn,'within',[T.sessTr T.seqType],{'sessTrans','seqType'},'subset',T.regType==r & T.regSide==hemi);
        end
        
        % Post-hoc t-tests for effect of seqType per sessTransition
        for r=roi
            for ss=sessTr
                fprintf('\n post-hoc t-test on the effect of seqType in sessTr %d in %s \n',ss,regname{r});
                ttestDirect(T.(var),[T.seqType T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi&T.sessTr==ss);
            end
        end
        
         % Post-hoc t-tests for effect of sessTransition per seqType
        for r=roi
            for ss=1:numel(seqType_label)
                fprintf('\n post-hoc t-test on the effect of sessTr in seqType %s in %s \n',seqType_label{ss},regname{r});
                ttestDirect(T.(var),[T.sessTr T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi&T.seqType==ss);
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
        hemi=1;
        regcorrType = 'zero';
        seqTypeCorr = 'seqSpec'; % seqSpec or seqType
        parcelType = 'Brodmann';
        sessTr=3; % which session transitions to include
        %3: sess3-4
        
        vararginoptions(varargin,{'reg','regcorrType','seqTypeCorr','parcelType','hemi'});
        T=load(fullfile(stabDir,sprintf('stability_acrSess_%s_sess%d-%d_%s_%sRegularised',seqTypeCorr,sessTr,sessTr+1,parcelType,regcorrType)));
        for r=reg
            fprintf('\n t-test on stability of sequence type across sess 3-4 in %s \n',regname{r})
            ttestDirect(T.acrSess_corr,[T.seqType T.sn],2,'paired','subset',T.regType==r&T.regSide==hemi);
        end
           
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
        
    case 'RDM_structure_corr'
        sn=[4:9,11:31];
        sessN=[1:4];
        betaChoice='multiPW';
        parcelType = 'Brodmann';
        roi=[1:8];
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice'});
        
        AllCorr=[];
        for ss=1:numel(sessN)
            %T = load(fullfile(regDir,sprintf('betas_FoSEx_sess%d',ss)));
            T{ss} = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,sessN(ss))));
        end
        for s=sn
            D = load(fullfile(glmSessDir{ss},subj_name{s},'SPM_info'));
            for s1=1:numel(sessN)
                for s2=s1:numel(sessN)
                    for r=roi          
                        T1 = getrow(T{s1},T{s1}.SN==s & T{s1}.region==r);
                        T2 = getrow(T{s2},T{s2}.SN==s & T{s2}.region==r);
                        % get betas
                        switch(betaChoice)
                            case 'multiPW'
                                betas{1} = T1.betaW{:};
                                betas{2} = T2.betaW{:};
                            case 'uniPW'
                                betas{1} = T1.betaUW{:};
                                betas{2} = T2.betaUW{:};
                            case 'raw'
                                betas{1} = T1.betaRAW{:};
                                betas{2} = T2.betaRAW{:};
                        end
                        for st=1:2 %seqType
                            % split into first (i), second (j), seqtype st
                            betas1 = betas{1}(D.seqType==st,:);
                            betas2 = betas{2}(D.seqType==st,:);
                            % info
                            D1 = getrow(D,D.seqType==st);
                            if st==2
                                D1.seqNumb = D1.seqNumb-6;
                            end
                            if (s1==s2)   % Same session
                                K = splitHalfCorr(betas1,D1.run,D1.seqNumb,'withinSes');
                            else
                                % corr across exe
                                % define condVec, partVec
                                partVec = [{D1.run} {D1.run}];
                                condVec = [{D1.seqNumb} {D1.seqNumb}];
                                data    = [{betas1} {betas2}];
                                K = splitHalfCorr(data,partVec,condVec,'acrossSes');
                            end;
                            K.sn        = s;
                            K.regNum    = r;
                            if r<9
                                K.regType   = r;
                                K.regSide   = 1;
                            else
                                K.regType   = r-8;
                                K.regSide   = 1;
                            end
                            K.ses1      = s1;
                            K.ses2      = s2;
                            K.seqType   = st;
                            AllCorr=addstruct(AllCorr,K);
                        end; % seqtype
                    end; % region
                end; %session 1
            end; % session 2
        fprintf('Done %s\n',subj_name{s});
        end;
        
        %  save the structure
        save(fullfile(stabDir,sprintf('corr_splitHalf_%s_%s',parcelType,betaChoice)),'-struct','AllCorr');
    case 'PLOT_structure_corr'
        % plotting function for split-half correlation - within or across
        % repetition (for trained / untrained)
        roi=[1:8];
        betaChoice='multiPW';
        parcelType='Brodmann';
        metric='corr_vox'; % corr_vox; corr_RDM
        vararginoptions(varargin,{'roi','metric','betaChoice'});
       % T = load(fullfile(repSupDir,sprintf('corr_splitHalf_NEW_exe_%s',betaChoice)));
        T = load(fullfile(stabDir,sprintf('corr_splitHalf_%s_%s',parcelType,betaChoice)));
        
        if strcmp(metric,'corr_RDMAll')
            for r=roi
                R = getrow(T,T.regNum==r);
                figure(1)
                subplot(1,numel(roi),r)
                plt.line([R.ses2>3 R.ses2],R.(metric),'split',R.seqType,'subset',R.ses2==R.ses1+1,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
                drawline(0,'dir','horz');
                title(regname{r});
                plt.match('y');
                if r==1
                    ylabel('RDM correlation 1st-2nd');
                else
                    ylabel('');
                end
            end
        else
            for r=roi
                R = getrow(T,T.regNum==r);
                figure
                subplot(1,2,1)
                plt.line([R.ses2>3 R.ses2],R.(metric),'split',R.seqType,'subset',R.ses2==R.ses1+1,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
                title(sprintf('%s - across ses',regname{r}));
                subplot(1,2,2)
                plt.line([R.ses2>3 R.ses2],R.(metric),'split',R.seqType,'subset',R.ses2==R.ses1,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
                title(sprintf('%s - within ses',regname{r}));
                plt.match('y');
            end
        end
        % split into within and across
   
    case 'DIMENSION_acrossSess'
        roi=1:16;
        parcelType='Brodmann';
        metricType = 'individual'; % individual or cumulative dimensions
        dataType = 'beta'; % beta or G
        sessN=1:4;
        sn=[5:9,11:31];
        vararginoptions(varargin,{'roi','parcelType','dataType','metricType','sn'});
        
        partVec=kron((1:8)',ones(12,1));
        condVec=kron(ones(8,1),(1:12)');
        seqVec=condVec;
        seqVec(seqVec<7)=1;
        seqVec(seqVec>6)=2;
        pVec = partVec(seqVec==1);
        DD=[];
        for ss1=sessN
            T1 = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss1)));
            for ss2=(ss1+1):max(sessN)
                T2 = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss2)));
                for r=roi
                    for s=sn
                        t1=getrow(T1,T1.SN==s&T1.region==r);
                        t2=getrow(T2,T2.SN==s&T2.region==r);
                        for st=1:2 % sequence type
                            data1 = t1.betaW{1}(seqVec==st,:);
                            data2 = t2.betaW{1}(seqVec==st,:);
                            if strcmp(dataType,'G') 
                                G1 = pcm_estGCrossval(data1,partVec(seqVec==st),condVec(seqVec==st));
                                G2 = pcm_estGCrossval(data2,partVec(seqVec==st),condVec(seqVec==st));
                                % now double center both Gs
                                H=eye(6)-ones(6,6)./6;  % centering matrix!
                                U1 = H*G1*H';
                                U2 = H*G2*H';
                            elseif strcmp(dataType,'beta')
                                Z=pcm_indicatorMatrix('identity',condVec(seqVec==st));
                                A1=zeros(size(data1));
                                A2=zeros(size(data2));
                                % first remove the mean for each run
                                for i=1:numel(unique(partVec))
                                    X = Z(pVec==i,:);
                                    Y1 = data1(pVec==i,:);
                                    Y2 = data2(pVec==i,:);
                                    A1(pVec==i,:) = pinv(X)*Y1;
                                    A2(pVec==i,:) = pinv(X)*Y2;
                                    %  subtract the mean across conditions 
                                    A1(pVec==i,:) = bsxfun(@minus,A1(pVec==i,:),mean(A1(pVec==i,:)));
                                    A2(pVec==i,:) = bsxfun(@minus,A2(pVec==i,:),mean(A2(pVec==i,:)));
                                end;
                                % now calculate the mean pattern
                                U1=pinv(Z)*A1;
                                U2=pinv(Z)*A2;
                            end
                            [u,v,w]=svd(U1);
                            % now calculate correlations for different
                            % dimensions
                            for d=1:6
                                if strcmp(metricType,'cumulative')
                                    k=1:d;
                                elseif strcmp(metricType,'individual')
                                    k=d;
                                end
                                U1_recon = u(:,k)*v(k,k)*w(:,k)';
                                D.corrDim   = corr(U1_recon(:),U2(:));
                                D.dim       = d;
                                D.sn        = s;
                                D.roi       = r;
                                D.seqType   = st;
                                D.regType   = t1.regType;
                                D.regSide   = t1.regSide;
                                D.ses1      = ss1;
                                D.ses2      = ss2;
                                DD = addstruct(DD,D);
                            end
                        end
                    end
                end
                fprintf('Done all reg: sess%d-sess%d\n',ss1,ss2);
            end
        end
        save(fullfile(stabDir,sprintf('Dim_acrossSess_%s_%sDim_%s',parcelType,metricType,dataType)),'-struct','DD');
    case 'PLOT_dim_acrossSess'
        % plot the effect of dimension in each session
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        metricType = 'individual'; % individual or cumulative dimensions
        dataType = 'beta'; % beta or G
        vararginoptions(varargin,{'parcelType','hemi','metricType','dataType'});
        
        T=load(fullfile(stabDir,sprintf('Dim_acrossSess_%s_%sDim_%s',parcelType,metricType,dataType)));
        for r=roi
            t=getrow(T,T.roi==r&T.regSide==hemi);
            figure
            for ss=1:3
                subplot(2,3,ss)
                plt.line(t.dim,t.corrDim,'split',t.seqType,'subset',t.ses1==1&t.ses2==ss+1&t.dim<6,'style',stySeq,'leg',{'trained','untrained'});
                ylabel('correlation');
                title(sprintf('%s sess1-%d',regname{r},ss+1));
                subplot(2,3,ss+3)
                plt.line(t.dim,t.corrDim,'split',t.seqType,'subset',t.ses1==ss&t.ses2==ss+1&t.dim<6,'style',stySeq,'leg',{'trained','untrained'});
                xlabel('dimension');
                ylabel('correlation');
                title(sprintf('%s sess%d-%d',regname{r},ss,ss+1));
            end
        end
    case 'PLOT_sess_acrossDim'
        % plot the effect of session for each dimension
                roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        metricType = 'individual'; % individual or cumulative dimensions
        dataType = 'beta'; % beta or G
        vararginoptions(varargin,{'parcelType','hemi','metricType','dataType'});
        
        T=load(fullfile(stabDir,sprintf('Dim_acrossSess_%s_%sDim_%s',parcelType,metricType,dataType)));
        for r=roi
            t=getrow(T,T.roi==r&T.regSide==hemi);
            figure
            for dim=1:6
                subplot(2,6,dim)
                plt.line(t.ses2,t.corrDim,'split',t.seqType,'subset',t.ses1==1&t.dim==dim,'style',stySeq,'leg',{'trained','untrained'});
                ylabel('correlation');
                title(sprintf('%s sess1 to others - dim%d',regname{r},dim));
                subplot(2,6,dim+6)
                plt.line(t.ses2,t.corrDim,'split',t.seqType,'subset',t.ses1==t.ses2-1&t.dim==dim,'style',stySeq,'leg',{'trained','untrained'});
                xlabel('session');
                ylabel('correlation');
                title(sprintf('%s subsequent ses - dim%d',regname{r},dim));
            end
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
        
    case 'PCM_extremeSimulation' %old
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
    case 'PCM_constructReliability' % in use
        runEffect  = 'random';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = 1:8;
        parcelType='Brodmann';
        %sn=[4:9,11:31];
        sn=[5:9,11:28,30,31];
        modelType='specific';
        sessN=1:2; % need to be two sessions at the time
        AllReg=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','modelType'});
        for ss=1:numel(sessN)
            B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
        end
        for r = reg
            Data=[];
            for st=1:2 % per seqType
                for p=1:length(sn)
                    for ss = 1:numel(sessN)
                        glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        t = getrow(B{ss},B{ss}.SN==sn(p) & B{ss}.region==r);
                        switch (beta_choice)
                            case 'uw'
                                beta = t.betaUW{:};
                            case 'mw'
                                beta = t.betaW{:};
                            case 'raw'
                                beta = t.betaRAW{:};
                        end
                        indx = D.seqType==st;
                        if ss == 1
                            condVec{p} = D.seqNumb(indx==1,:); % conditions
                            partVec{p} = D.run(indx==1,:);
                            Data{p} = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                        else
                            condVec{p} = [condVec{p}; D.seqNumb(indx==1,:) + 12]; % treat 2nd session as additional conditions
                            partVec{p} = [partVec{p}; D.run(indx==1,:) + 8];  % runs/partitions of 2nd session as additional runs
                            Data{p} = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                        end;
                    end; % session
                end; % subj 
                % construct models
                switch modelType
                    case 'specific_sess'
                        M = pcm_stabilityModel_specific;
                        C = pcm_correlation(Data,partVec,condVec,M{4},runEffect,modelType);
                    case 'generic_sess'
                        M = pcm_stabilityModel_generic;
                        C = pcm_correlation(Data,partVec,condVec,M{4},runEffect,modelType);
                    case 'specific_noSess'
                        M = pcm_stabilityModel_specific_noSess;
                        C = pcm_correlation(Data,partVec,condVec,M{3},runEffect,'specific');
                    case 'generic_noSess'
                        M = pcm_stabilityModel_generic_noSess;
                        C = pcm_correlation(Data,partVec,condVec,M{3},runEffect,'generic');
                end
                T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);                
                T.roi = ones(size(T.SN))*r;
                T.regType = ones(size(T.SN))*regType(r);
                T.regSide = ones(size(T.SN))*regSide(r);
                T.seqType = ones(size(T.SN))*st;
                T=rmfield(T,{'theta_hat','thetaCr'});
                AllReg=addstruct(AllReg,T);
                AllReg=addstruct(AllReg,C);
            end
            fprintf('Done sess: %d-%d model-%s reg: %d/%d\n\n',sessN(1),sessN(2),modelType,r,max(reg))
        end
        fprintf('Done sess: %d-%d model-%s\n\n\n\n',sessN(1),sessN(2),modelType);
        % save output
        save(fullfile(stabDir,sprintf('PCM_stability_%s_%s_sess%d-sess%d.mat',parcelType,modelType,sessN(1),sessN(2))),'-struct','AllReg');
    case 'PCM_constructReliability_seqType' %old
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
    case 'PCM_constructReliability_hyperModel' %old
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
    case 'PCM_reliability_specificCorr'
        % run models specifying correlation across sessions (from 0:0.1:1)
        runEffect  = 'random';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = 1:8;
        parcelType='Brodmann';
        sn=[5:9,11:31];
        modelType='specific';
        regSelect='all';
        corrSpec = 0:0.1:1;
        sessN=1:2; % need to be two sessions at the time
        AllReg=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','modelType'});
        for ss=1:numel(sessN)
            B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
        end
        if strcmp(parcelType,'162tessels')
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
            regType=reg;
            regType(regType>158)=regType(regType>158)-158;
            regSide=ones(size(regType));
            regSide(reg>158)=2;
        end
        
        for r = reg
            for st=1:2 % per seqType
                for c = 1:length(corrSpec)
                    for p=1:length(sn)
                        for ss = 1:numel(sessN)
                            glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                            t = getrow(B{ss},B{ss}.SN==sn(p) & B{ss}.region==r);
                            switch (beta_choice)
                                case 'uw'
                                    beta = t.betaUW{:};
                                case 'mw'
                                    beta = t.betaW{:};
                                case 'raw'
                                    beta = t.betaRAW{:};
                            end
                            indx = D.seqType==st;
                            if ss == 1
                                condVec{p} = D.seqNumb(indx==1,:); % conditions
                                partVec{p} = D.run(indx==1,:);
                                Data{p} = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                            else
                                condVec{p} = [condVec{p}; D.seqNumb(indx==1,:) + 12]; % treat 2nd session as additional conditions
                                partVec{p} = [partVec{p}; D.run(indx==1,:) + 8];  % runs/partitions of 2nd session as additional runs
                                Data{p} = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                            end;
                        end; % session
                    end; % subj
                    % construct models
                    switch modelType
                        case 'specific'
                            M = pcm_stabilityModel_specificCorrelation(corrSpec(c));
                        case 'generic'
                            M = pcm_stabilityModel_genericCorrelation(corrSpec(c)); 
                    end
                    T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*r;
                    T.corr = ones(size(T.SN))*c;
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.seqType = ones(size(T.SN))*st;
                    T=rmfield(T,{'reg','theta_hat','thetaCr'});
                    AllReg=addstruct(AllReg,T);
                    fprintf('Done corr %d/%d sess: %d-%d reg: %d/%d\n\n',c,length(corrSpec),sessN(1),sessN(2),r,length(reg))
                end; % specific correlation
            end; % seqType
            fprintf('Done all corr sess: %d-%d reg: %d/%d\n\n',sessN(1),sessN(2),r,max(reg))
        end; % region
        fprintf('Done sess: %d-%d\n\n\n\n',sessN(1),sessN(2));
        % save output
        save(fullfile(stabDir,sprintf('PCM_stability_specCorr_%s_%s_sess%d-sess%d.mat',parcelType,modelType,sessN(1),sessN(2))),'-struct','AllReg');
    case 'PCM_reliability_allCorr'
        % assess all correlations at the same time
         % run models specifying correlation across sessions (from 0:0.1:1)
        runEffect  = 'random';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = 1:8;
        parcelType='Brodmann';
        sn=[5:9,11:28,30,31];
        modelType='generic';
        regSelect='all';
        checkCorr=1; % check if correlation exact and plot predicted Gs
        sessN=1:2; % need to be two sessions at the time
        AllReg=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','modelType'});
        for ss=1:numel(sessN)
            B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
        end
        if strcmp(parcelType,'162tessels')
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
            regType=reg;
            regType(regType>158)=regType(regType>158)-158;
            regSide=ones(size(regType));
            regSide(reg>158)=2;
        end       
        
        % construct omdels
        switch modelType
            case 'generic_noSess'
                M = pcm_stability_genericCorr_noSess;
            case 'specific_noSess'
                M = pcm_stability_specificCorr_noSess;
            case 'generic_sess'
                M = pcm_stability_genericCorr_sess;
            case 'specific_sess'
                M = pcm_stability_specificCorr_sess;
            case 'meanPattern'
                M = pcm_stability_meanPatternCorr;
        end
        
        for r = reg
            for st=1:2 % per seqType
                for p=1:length(sn)
                    for ss = 1:numel(sessN)
                        glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        t = getrow(B{ss},B{ss}.SN==sn(p) & B{ss}.region==r);
                        switch (beta_choice)
                            case 'uw'
                                beta = t.betaUW{:};
                            case 'mw'
                                beta = t.betaW{:};
                            case 'raw'
                                beta = t.betaRAW{:};
                        end
                        indx = D.seqType==st;
                        if ss == 1
                            condVec{p} = D.seqNumb(indx==1,:); % conditions
                            partVec{p} = D.run(indx==1,:);
                            Data{p} = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                        else
                            condVec{p} = [condVec{p}; D.seqNumb(indx==1,:) + 12]; % treat 2nd session as additional conditions
                            partVec{p} = [partVec{p}; D.run(indx==1,:) + 8];  % runs/partitions of 2nd session as additional runs
                            Data{p} = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                        end;
                    end; % session
                end; % subj
                % here just test if correlations correct
                if checkCorr
                    [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);  
                    figure
                    corrPred=zeros(11,1);
                    for c=1:11
                        subplot(3,4,c)
                        imagesc(G_pred{c});
                        corrPred(c)=calcCorr(G_pred{c});
                    end
                end
                % group fit
                T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                % individual fit
                I = pcm_fitModelIndivid(Data,M,partVec{p},condVec{p},'runEffect',runEffect);
                T.individ_likelihood = I.likelihood;
                T.individ_bayes     = bsxfun(@minus,T.individ_likelihood,T.individ_likelihood(:,1));
                T.individ_noise = I.noise;
                T.roi = ones(size(T.SN))*r;
                T.regType = ones(size(T.SN))*regType(r);
                T.regSide = ones(size(T.SN))*regSide(r);
                T.seqType = ones(size(T.SN))*st;
            %    T=rmfield(T,{'reg','theta_hat','thetaCr'});
                T=rmfield(T,{'theta_hat','thetaCr'});
                AllReg=addstruct(AllReg,T);
                fprintf('Done modelType: %s seqType %d sess: %d-%d reg: %d/%d\n\n',modelType,st,sessN(1),sessN(2),r,length(reg));
            end; % seqType
        end; % region
        fprintf('Done all:\tmodelType: %s \tsess: %d-%d\n\n\n\n',modelType,sessN(1),sessN(2));
        % save output
        save(fullfile(stabDir,sprintf('PCM_allCorr_%s_%s_sess%d-sess%d.mat',parcelType,modelType,sessN(1),sessN(2))),'-struct','AllReg');
    case 'PCM_noiseCeiling'
        % construct noise ceiling models here
        runEffect  = 'random';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = 1:8;
        parcelType='Brodmann';
        sn=[5:9,11:28,30,31];
        regSelect='all';
        checkCorr=0; % check if correlation exact and plot predicted Gs
        sessN=1:2; % need to be two sessions at the time
        AllReg=[];
        modelType='generic';
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','modelType'});
        for ss=1:numel(sessN)
            B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
        end
        if strcmp(parcelType,'162tessels')
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
            regType=reg;
            regType(regType>158)=regType(regType>158)-158;
            regSide=ones(size(regType));
            regSide(reg>158)=2;
        end
        for r = reg
            for st=1:2 % per seqType
                for p=1:length(sn)
                    for ss = 1:numel(sessN)
                        glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        t = getrow(B{ss},B{ss}.SN==sn(p) & B{ss}.region==r);
                        switch (beta_choice)
                            case 'uw'
                                beta = t.betaUW{:};
                            case 'mw'
                                beta = t.betaW{:};
                            case 'raw'
                                beta = t.betaRAW{:};
                        end
                        indx = D.seqType==st;
                        if ss == 1
                            condVec{p} = D.seqNumb(indx==1,:); % conditions
                            partVec{p} = D.run(indx==1,:);
                            Data{p} = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                        else
                            condVec{p} = [condVec{p}; D.seqNumb(indx==1,:) + 12]; % treat 2nd session as additional conditions
                            partVec{p} = [partVec{p}; D.run(indx==1,:) + 8];  % runs/partitions of 2nd session as additional runs
                            Data{p} = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                            G(p,:)    = rsa_vectorizeIPMfull(pcm_estGCrossval(Data{p},partVec{1},condVec{1}));
                        end;
                    end; % session
                end; % subj
                                
                % individual fit
                for p=1:length(sn)
                    % construct models here
                    % 1) Null model:
                    M{1}.type       = 'feature';
                    M{1}.numGparams = 1;
                    M{1}.name       = 'null';
                    M{1}.Ac(:,1:12,1)  = zeros(12);
                    % 2) Session effect:
                    M{2}.type       = 'feature';
                    M{2}.numGparams = 2;
                    M{2}.name       = 'zero_corr';
                    M{2}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
                    M{2}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
                    % 3) Zero correlation:
                    M{3}.type       = 'feature';
                    M{3}.numGparams = 14;
                    M{3}.name       = 'zero_corr';
                    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
                    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
                    % for sequence-specific modelling- one parameter per sequence
                    for i=1:6
                        A=zeros(6);
                        A(i,i)=1;
                        M{3}.Ac(:,3:8,2+i)      = [A;zeros(6)];     % Unique sess1 sequence patterns
                        M{3}.Ac(:,9:14,8+i)     = [zeros(6);A];     % Unique sess2 sequence pattterns
                    end;
                    % 4) Perfect correlation:
                    M{4}.type         = 'feature';
                    M{4}.numGparams   = 14;
                    M{4}.name         = 'perfect_corr';
                    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
                    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
                    % for sequence-specific modelling- one parameter per sequence
                    for i=1:6
                        A=zeros(6);
                        A(i,i)=1;
                        M{4}.Ac(:,3:8,2+i)   = [A;zeros(6)];       % Unique sess1 sequence patterns
                        M{4}.Ac(:,3:8,8+i)   = [zeros(6);A];       % Same sess2 sequence pattterns
                    end;
                    % 5) Noise ceiling - upper: (including all subjects in that
                    % group)
                    M{5}.type       = 'feature';
                    M{5}.numGparams = 1;
                    M{5}.name       = 'noiseCeiling_upp';
                    if rem(sn(p),2)==1
                        M{5}.Ac = rsa_squareIPMfull(nanmean(G(rem(sn,2)==1,:),1));% G for model
                    else
                        M{5}.Ac = rsa_squareIPMfull(nanmean(G(rem(sn,2)==0,:),1));% G for model
                    end
                    % 6) Noise ceiling - lower: (include all subjects but the
                    % one considered here)
                    M{6}.type       = 'feature';
                    M{6}.numGparams = 1;
                    M{6}.name       = 'noiseCeiling_lower';
                    if rem(sn(p),2)==1
                        idx = rem(sn,2)==1;
                        idx(p) = 0;
                        M{6}.Ac = rsa_squareIPMfull(nanmean(G(idx,:),1));% G for model
                    else
                        idx = rem(sn,2)==0;
                        idx(p) = 0;
                        M{6}.Ac = rsa_squareIPMfull(nanmean(G(idx,:),1));% G for model
                    end
                    % 7) Noise ceiling 
                    M{7}.type       = 'feature';
                    M{7}.numGparams = 1;
                    M{7}.name       = 'noiseCeiling';
                    M{7}.Ac         = rsa_squareIPMfull(G(p,:));
                    % 8) Noise ceiling 
                    M{8}.type = 'freedirect';
                    M{8}.numGparams = 0;
                    M{8}.theta0 = [];
                    M{8}.name = 'noice_ceiling';
                    
                    
                    I = pcm_fitModelIndivid(Data(p),M,partVec{p},condVec{p},'runEffect',runEffect);
                    T.SN = p;
                    T.likelihood = I.likelihood;
                    T.logBayes   = bsxfun(@minus,T.likelihood,T.likelihood(:,1));
                    T.roi       = r;
                    T.regType   = regType(r);
                    T.regSide   = regSide(r);
                    T.seqType   = st;
                    AllReg = addstruct(AllReg,T);
                end
                fprintf('Done modelType: noiseCeiling seqType %d sess: %d-%d reg: %d/%d\n\n',st,sessN(1),sessN(2),r,length(reg));
            end; % seqType
        end; % region
        fprintf('Done all:\tmodelType: noiseCeiling \tsess: %d-%d\n\n\n\n',sessN(1),sessN(2));
        % save output
        save(fullfile(stabDir,sprintf('PCM_corr_%s_noiseCeiling_sess%d-sess%d.mat',parcelType,sessN(1),sessN(2))),'-struct','AllReg');
        
        
    case 'PCM_plot_allSess' % in use
        reg = 1:8;
        sessName={'1-2','2-3','3-4'};
        seqType={'trained','untrained'};
        modelType='specific'; % generic or specific
        parcelType = 'Brodmann';
        numTrans=3; % number of session transitions - 2:sess1-2-3; 3:sess1-2-3-4
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','numTrans','modelType','parcelType'});
        
        if any(strcmp(modelType,{'specific','generic','specific_sess','generic_sess'}))
            modelLabel = {'sess','sess+seq','sess+seq+corr','sess+seq+perfect'};
            nModel = 4;
        elseif any(strcmp(modelType,{'specific_noSess','generic_noSess'}))
            modelLabel = {'seq','seq+corr','seq+perfect'};
            nModel = 3;
        end 
        
        T=[];
        for tr=1:numTrans
            R=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_sess%d-sess%d.mat',parcelType,modelType,tr,tr+1)));
            R.sessTr=ones(size(R.SN))*tr;
            T=addstruct(T,R);
        end
        TT=T;
        for st=1:2 % for both sequence types
            T=getrow(TT,TT.seqType==st);
            % models: null, session, session+seq, session+seq+flexCorr,
            % session+seq+perfectCorr
            % subtract the null model away
            bayes_subtr     = bsxfun(@minus,T.bayesEst,T.bayesEst(:,1));
            bayes_subtr2    = bayes_subtr(:,(2:end));
            T2.bayesEst     = bayes_subtr2(:); % don't include the null model
            modelInd=[];
            for n=1:nModel
                modelInd    = [modelInd; ones(size(bayes_subtr(:,1)))*n];
            end
            T2.modelInd     = modelInd;
            T2.roi          = repmat(T.roi,nModel,1);
            T2.sessTr       = repmat(T.sessTr,nModel,1);
            T2.SN           = repmat(T.SN,nModel,1);
            T2.seqType      = ones(size(T2.SN))*st;
            TS{st}=T2;
            % rearranging how the data structure is arranged
            figure
            for r=reg
                subplot(1,numel(reg),r)
                plt.bar([T2.sessTr>2 T2.sessTr],T2.bayesEst,'subset',T2.roi==r,'split',T2.modelInd,'leg',modelLabel,'leglocation','north');
                plt.match('y');
                drawline(0,'dir','horz');
                title(sprintf('%s',regname{r}));
                set(gca,'XTickLabel',sessName);
                if r==1
                    ylabel(sprintf('Log-Bayes %s - %s model',seqType{st},modelType));
                    xlabel('Sess transitions');
                else
                    ylabel('');
                end
            end
            %subtract also the session effect
            bayes_subtr3    = bsxfun(@minus,bayes_subtr2,bayes_subtr2(:,1));
            bayes_subtr3    = bayes_subtr3(:,(2:end));
            T3.bayesEst     = bayes_subtr3(:);
            T3.modelInd     = T2.modelInd(T2.modelInd~=1)-1;
            T3.roi          = repmat(T.roi,nModel-1,1);
            T3.sessTr       = repmat(T.sessTr,nModel-1,1);
            T3.seqType      = ones(size(T3.roi))*st;
            T3.SN           = repmat(T.SN,nModel-1,1);
            % save this for direct sequence type comparison
            S{st}=T3;
            figure
            for r=reg
                subplot(1,numel(reg),r)
                plt.bar([T3.sessTr>2 T3.sessTr],T3.bayesEst,'subset',T3.roi==r,'split',T3.modelInd,'leg',modelLabel(2:end),'leglocation','north');
                plt.match('y');
                drawline(0,'dir','horz');
                title(sprintf('%s',regname{r}));
                set(gca,'XTickLabel',sessName);
                if r==1
                    ylabel(sprintf('Log-Bayes %s - %s model',seqType{st},modelType));
                    xlabel('Sess transitions');
                else
                    ylabel('');
                end
            end
        end
        
        SS = [];
        SS = addstruct(SS,S{1});
        SS = addstruct(SS,S{2});
        MM = [];
        MM = addstruct(MM,TS{1});
        MM = addstruct(MM,TS{2});
        for r=reg
            figure
            for m=1:nModel-1 % model types
                subplot(1,nModel-1,m)
                plt.bar([SS.sessTr>2 SS.sessTr],SS.bayesEst,'subset',SS.roi==r&SS.modelInd==m,'split',SS.seqType,'leg',{'trained','untrained'},'style',stySeq);
                title(sprintf('%s model %s',regname{r},modelLabel{m+1}));
                if m==1
                    ylabel('Log-Bayes factor');
                else
                    ylabel('');
                end
            end
            figure(20)
            subplot(numel(reg),1,r)
            plt.bar([SS.sessTr>2 SS.sessTr],SS.bayesEst,'subset',SS.roi==r&SS.modelInd==3,'split',SS.seqType,'leg',{'trained','untrained'},'style',stySeq);
            title(sprintf('%s model %s',regname{r},modelLabel{4}));   

        end
        % save SS,MM
        save(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans',parcelType,modelType)),'-struct','SS'); 
        save(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans_sess',parcelType,modelType)),'-struct','MM'); % includes the session effect
    case 'PCM_stats_allSess_perfect'
        reg = 1:8;
        modelType='specific'; % generic or specific
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','numTrans','modelType','parcelType'});
        
        % here assessing the stats for figure 20 in the above case
        % fit of perfect model for trained vs. untrained
        N=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans',parcelType,modelType))); % with session effect
      
        for r=reg
            % first test for effect of sequence over session only
            for t=1:3 % all transitions - seq specificity > session only
                fprintf('%s transition %d - trained vs. untrained - perfect model\n',regname{r},t)
                n=getrow(N,N.modelInd==2&N.roi==r&N.sessTr==t);
                ttestDirect(N.bayesEst(N.modelInd==3&N.roi==r&N.sessTr==t)-N.bayesEst(N.modelInd==2&N.roi==r&N.sessTr==t),[n.seqType n.SN],2,'paired');
            end
        end
    case 'PCM_plot_sessSeq'
        % here plot comparison between models
        reg = 1:8;
        modelType='specific'; % generic or specific
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','numTrans','modelType','parcelType'});
        T=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans_sess',parcelType,modelType))); % with session effect
        N=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans',parcelType,modelType))); % no session effect
        for r=reg
            figure
            subplot(221)
            style.use('Trained');
            plt.bar([T.sessTr>2 T.sessTr],T.bayesEst,'subset',T.roi==r&T.seqType==1&T.modelInd<3,'split',T.modelInd,'leg',{'session','sess+seq'},'leglocation','northeast');
            ylabel('log-likelihood'); title(sprintf('%s - trained',regname{r}));
            subplot(222)
            style.use('Untrained');
            plt.bar([T.sessTr>2 T.sessTr],T.bayesEst,'subset',T.roi==r&T.seqType==2&T.modelInd<3,'split',T.modelInd,'leg',{'session','sess+seq'},'leglocation','northeast');
            ylabel('log-likelihood'); title(sprintf('%s - untrained',regname{r}));
            subplot(2,2,3:4)
            style.use('Seq')
           % plt.bar([T.sessTr>2 T.sessTr],(T.bayesEst(T.modelInd==2))-T.bayesEst(T.modelInd==1),'subset',T.roi==r&T.seqType==1&T.modelInd<3,'split',T.modelInd,'leg',{'session','sess+seq'},'leglocation','northeast')
            plt.bar([N.sessTr>2 N.sessTr],N.bayesEst,'subset',N.roi==r&N.modelInd==1,'split',N.seqType,'leg',{'trained','untrained'});
            ylabel('log-likelihood'); title('trained vs. untrained - seq specific');
            
            style.use('SeqSess');
            figure
            plt.bar([N.sessTr>2 N.sessTr],N.bayesEst,'subset',N.roi==r,'split',[N.modelInd N.seqType],'leg',{'trained','untrained'});
            ylabel('log-likelihood'); title(sprintf('%s - all models trained vs. untrained',regname{r}));
            
            figure(100) % here just the effect of session in 3-4 for trained / untrained
            subplot(numel(reg),2,(r-1)*2+1)
            plt.bar(T.sessTr,T.bayesEst,'subset',T.roi==r & T.modelInd==2,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','northeast','style',stySeq);
            ylabel('log-likelihood'); title(sprintf('%s session effect',regname{r}));
            subplot(numel(reg),2,r*2); % direct comparison trained - untrained
            T2 = getrow(T,T.seqType==1);
            style.use('gray');
            plt.bar(T2.sessTr,T.bayesEst(T.seqType==1)-T.bayesEst(T.seqType==2),'subset',T2.roi==r & T2.modelInd==2);
            ylabel(''); title(sprintf('%s trained - untrained',regname{r}));
            hold on; drawline(0,'dir','horz');
            
        end
    case 'PCM_stats_sessSeq'
        reg = 1:8;
        modelType='specific'; % generic or specific
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','numTrans','modelType','parcelType'});
        
        T=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans_sess',parcelType,modelType))); % no session effect
        N=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans',parcelType,modelType))); % with session effect
      
        for r=reg
            % first test for effect of sequence over session only
            for t=1:3 % all transitions - seq specificity > session only
                fprintf('%s - transition %d - trained\n',regname{r},t);
                ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==t & T.seqType==1 & T.modelInd<3);
                fprintf('%s - transition %d - untrained\n',regname{r},t);
                ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==t & T.seqType==2 & T.modelInd<3);
                fprintf('transition %d - trained vs. untrained\n',t)
                ttestDirect(N.bayesEst,[N.seqType N.SN],2,'paired','subset',N.roi==r & N.sessTr==t & N.modelInd==1);
                fprintf('%s session effect in trained vs. untrained - transition %d',regname{r},t)
                ttestDirect(T.bayesEst,[T.seqType T.SN],2,'paired','subset',T.roi==r & T.sessTr==t & T.modelInd==1);
            end
            for st=1:2 % here session transitions for session type
                fprintf('%s seqType %d - session 1-2\n',regname{r},st);
                ttestDirect(T.bayesEst,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.seqType==st & ismember(T.sessTr,[1,2]) & T.modelInd==1);
                fprintf('%s seqType %d - session 2-3\n',regname{r},st);
                ttestDirect(T.bayesEst,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.seqType==st & ismember(T.sessTr,[2,3]) & T.modelInd==1);
                fprintf('%s seqType %d - session 1-3\n',regname{r},st);
                ttestDirect(T.bayesEst,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.seqType==st & ismember(T.sessTr,[1,3]) & T.modelInd==1);
            end
            keyboard;
        end
        
    case 'PCM_seqType_stats'
        modelType = 'specific';
        parcelType = 'Brodmann';
        roi=[1:8];
        
        vararginoptions(varargin,{'modelType','parcelType','roi'});
        T = load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allTrans',parcelType,modelType)));
        
        % comparison against 0
        for r=roi
            fprintf('relative to 0: reg %d\n',r);
            for m=unique(T.modelInd)'
                for st=1:2 
                    fprintf('model %d - seqType %d\n',m,st);
                    ttestDirect(T.bayesEst,[T.SN],2,'onesample','subset',T.roi==r & T.sessTr==1 & T.seqType==st & T.modelInd==m);
                    ttestDirect(T.bayesEst,[T.SN],2,'onesample','subset',T.roi==r & T.sessTr==2 & T.seqType==st & T.modelInd==m);
                    ttestDirect(T.bayesEst,[T.SN],2,'onesample','subset',T.roi==r & T.sessTr==3 & T.seqType==st & T.modelInd==m);
                end

            end
        end
        % trained vs. untrained
         for r=roi
            fprintf('trained vs. untrained: reg %d\n',r);
            for m=unique(T.modelInd)'
                    fprintf('model %d\n',m);
                    ttestDirect(T.bayesEst,[T.seqType T.SN],2,'paired','subset',T.roi==r & T.sessTr==1 & T.modelInd==m);
                    ttestDirect(T.bayesEst,[T.seqType T.SN],2,'paired','subset',T.roi==r & T.sessTr==2 & T.modelInd==m);
                    ttestDirect(T.bayesEst,[T.seqType T.SN],2,'paired','subset',T.roi==r & T.sessTr==3 & T.modelInd==m);
            end
         end
        % model comparison
        for r=roi
            fprintf('trained vs. untrained: reg %d\n',r);
            for ss=unique(T.sessTr)'
                    fprintf('sessTr%d - seqType 1\n',ss);
                    ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==ss & ismember(T.modelInd,[1,2]) & T.seqType==1);
                    ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==ss & ismember(T.modelInd,[1,3]) & T.seqType==1);
                    ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==ss & ismember(T.modelInd,[2,3]) & T.seqType==1);
                    fprintf('sessTr%d - seqType 2\n',ss);
                    ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==ss & ismember(T.modelInd,[1,2]) & T.seqType==2);
                    ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==ss & ismember(T.modelInd,[1,3]) & T.seqType==2);
                    ttestDirect(T.bayesEst,[T.modelInd T.SN],2,'paired','subset',T.roi==r & T.sessTr==ss & ismember(T.modelInd,[2,3]) & T.seqType==2);
            end
         end
        
         for r=roi
             fprintf('reg %d\n',r);
             for m=unique(T.modelInd)'
                 for st=1:2 % seqType
                     fprintf('Model%d - seqType %d\n',m,st);
                     ttestDirect(T.bayesEst,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.modelInd==m & ismember(T.sessTr,[1,2]) & T.seqType==st);
                     ttestDirect(T.bayesEst,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.modelInd==m & ismember(T.sessTr,[2,3]) & T.seqType==st);
                     ttestDirect(T.bayesEst,[T.sessTr T.SN],2,'paired','subset',T.roi==r & T.modelInd==m & ismember(T.sessTr,[1,3]) & T.seqType==st);
                 end
             end
         end
    case 'PCM_plot_corr'
        reg = 1:8;
        hemi=1;
        sessName={'1-2','2-3','3-4'};
        modelType='specific'; % generic or specific
        parcelType = 'Brodmann';
        numTrans=3; % number of session transitions - 2:sess1-2-3; 3:sess1-2-3-4
        metric='r_model'; % options: r_naive, r_crossval, r_crossval_seqType, r_model
        vararginoptions(varargin,{'reg','sessN','metric','modelType','parcelType','hemi'});

        T=[];
        for tr=1:numTrans
            R=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_sess%d-sess%d.mat',parcelType,modelType,tr,tr+1)));
            R.sessTr=ones(size(R.SN))*tr;
            T=addstruct(T,R);
        end
        save(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allSess.mat',parcelType,modelType)),'-struct','T');
        TT=T;
       figure
       for r=reg
           subplot(1,max(reg),r)
           plt.line([TT.sessTr>2 TT.sessTr],TT.(metric),'split',TT.seqType,'subset',TT.regType==r&TT.regSide==hemi,'style',stySeq,'leg',{'trained','untrained'});
           xlabel('Sess transition');
           ylabel('Correlation');
           title(regname_cortex{r});
       end
    case 'PCM_plot_corrViolin'
        reg = 1:8;
        modelType='specific'; % generic or specific, specific_noSess
        parcelType = 'Brodmann';
        hemi=1;
        metric='r_model'; % options: r_naive, r_crossval, r_crossval_seqType, r_model
        vararginoptions(varargin,{'reg','sessN','metric','modelType','parcelType'});

        
        T = load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allSess',parcelType,modelType)));
        colRed = [252 146 114]./255;
        colBlue = [158 202 225]./255;
        for r=reg
            t = getrow(T,T.roi==r & T.regSide==hemi);
            figure
            subplot(211)
            distributionPlot(t.(metric)(t.sessTr==1 & t.seqType==1),'histOri','left','color',colRed,'widthDiv',[2 1],'xValues',1,'showMM',4);
            distributionPlot(t.(metric)(t.sessTr==1 & t.seqType==2),'histOri','right','color',colBlue,'widthDiv',[2 2],'xValues',1,'showMM',4);
            
            distributionPlot(t.(metric)(t.sessTr==2 & t.seqType==1),'histOri','left','color',colRed,'widthDiv',[2 1],'xValues',2,'showMM',4);
            distributionPlot(t.(metric)(t.sessTr==2 & t.seqType==2),'histOri','right','color',colBlue,'widthDiv',[2 2],'xValues',2,'showMM',4);
            
            distributionPlot(t.(metric)(t.sessTr==3 & t.seqType==1),'histOri','left','color',colRed,'widthDiv',[2 1],'xValues',3,'showMM',4);
            distributionPlot(t.(metric)(t.sessTr==3 & t.seqType==2),'histOri','right','color',colBlue,'widthDiv',[2 2],'xValues',3,'showMM',4);
            hold on;
            drawline([0,1],'dir','horz');
          %  xlabel('Session transition');
            ylabel('Correlation');
            title(regname{r});
            
            subplot(212)
            plt.dot(t.sessTr,t.(metric),'split',t.seqType,'style',stySeq,'leg',{'trained','untrained'});
            drawline([0,1],'dir','horz');
          %  title(regname{r});
            xlabel('Session transition');
            ylabel('Correlation');
            fprintf('%s: trained vs. untrained\n',regname{r});
            for s=1:3
                fprintf('%s: transition %d: trained vs. untrained\n',regname{r});
                ttestDirect(t.(metric),[t.seqType t.SN],2,'paired','subset',t.sessTr==s);
            end
            for st=1:2
                for s=1:2
                    fprintf('%s: seqType %d: transition:%d-%d\n',regname{r},st,s,s+1);
                    ttestDirect(t.(metric),[t.sessTr t.SN],2,'paired','subset',t.seqType==st & ismember(t.sessTr,[s,s+1]));
                end
            end
            
            keyboard;
        end
    case 'PCM_plot_corrModel'
      % plot all correlation models (0:0.1:1)
      reg = 1:8;
      hemi=1;
      sessName={'1-2','2-3','3-4'};
      modelType='specific'; % generic or specific
      parcelType = 'Brodmann';
      vararginoptions(varargin,{'parcelType','reg','hemi','sessName','modelType'});
      % load in all the data
      TT=[]; SS=[];
      for ss=1:length(sessName)
          T=load(fullfile(stabDir,sprintf('PCM_stability_specCorr_%s_%s_sess%d-sess%d.mat',parcelType,modelType,ss,ss+1)));
          % extract the relevant logBayes - seq+corr vs. seq only
          T.bayesCorr = bsxfun(@minus,T.bayesEst(:,4),T.bayesEst(:,3));
          T.sessTr = ones(size(T.roi))*ss;
          TT=addstruct(TT,T);
          S=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_sess%d-sess%d.mat',parcelType,modelType,ss,ss+1))); % extract true correlation
          S.sessTr=ones(size(S.SN))*ss;
          SS=addstruct(SS,S);
      end
      for r=reg
          figure
          [a,~,~]=pivottable([TT.corr],[TT.sessTr,TT.seqType],TT.bayesCorr,'mean','subset',TT.regType==r&TT.regSide==hemi);
          % create 0s for in-between sessions
          subplot(121)
          plot1=ribbon(a);
          alpha(plot1,.6);
          hold on;
          X=[0;6;6;0;0];
          Y=[0;0;10;10;0];
          Z=[5;5;5;5;5];
        %  plot3(X,Y,Z);
          h=fill3(X,Y,Z,'k');
          set(h,'FaceAlpha',0.2);
          title(regname_cortex{r});
          xlabel('Session transition');
          ylabel('Correlation model evaluated');
          zlabel('LogBayes relative to 0 corr');
          %plot3(7,0:0.1:10,0,'-o','Color',[0,0,0]);
          subplot(122)
          plt.bar([SS.sessTr>2 SS.sessTr],SS.r_model,'split',SS.seqType,'style',stySeq,'subset',SS.regType==r&SS.regSide==hemi,'leg',{'trained','untrained'});
          ylabel('estimated correlation from flexModel');
          title(regname_cortex{r});
      end
    case 'PCM_plot_noiseCeiling'
      reg = 1:8;
      hemi=1;
      sessName={'1-2','2-3','3-4'};
      parcelType = 'Brodmann';
      vararginoptions(varargin,{'parcelType','reg','hemi','sessName'});
      
      TT=[]; 
      for ss=1:length(sessName)
          T=load(fullfile(stabDir,sprintf('PCM_corr_%s_noiseCeiling_sess%d-sess%d.mat',parcelType,ss,ss+1)));
          % extract the relevant logBayes - seq+corr vs. seq only
          T.logBayes_sess = bsxfun(@minus,T.logBayes,T.logBayes(:,2)); % relative to the session factor
          T.logBayes_sess = T.logBayes_sess(:,2:end);
          T.sessTr = ones(size(T.roi))*ss;
          TT=addstruct(TT,T);
      end
      for r=reg
          t1=getrow(TT,TT.regType==r & TT.regSide==hemi & TT.seqType==1);
          t2=getrow(TT,TT.regType==r & TT.regSide==hemi & TT.seqType==2);
          figure
          subplot(241)
          barplot(t1.sessTr,t1.logBayes_sess(:,[2,3]),'subset',t1.sessTr==1);
          hold on;
          drawline(mean(t1.logBayes_sess(t1.sessTr==1,4)),'dir','horz','linewidth',10,'color',[0.8 0.8 0.8]);
          drawline(mean(t1.logBayes_sess(t1.sessTr==1,7)),'dir','horz','linewidth',10,'color',[0.5 0.5 0.5]);
          ylim([0 mean(t1.logBayes_sess(t1.sessTr==1,7))+5]);
          ylabel('logBayes'); set(gca,'XTickLabel',{'zero','perfect'});
          title(sprintf('trained-sessTr1: %s',regname{r}));
          subplot(245)
          barplot(t2.sessTr,t2.logBayes_sess(:,[2,3]),'subset',t2.sessTr==1);
          hold on;
          drawline(mean(t2.logBayes_sess(t2.sessTr==1,4)),'dir','horz','linewidth',10,'color',[0.8 0.8 0.8]);
          drawline(mean(t2.logBayes_sess(t2.sessTr==1,7)),'dir','horz','linewidth',10,'color',[0.5 0.5 0.5]);
          ylim([0 mean(t2.logBayes_sess(t2.sessTr==1,7))+5]);
          ylabel('logBayes'); set(gca,'XTickLabel',{'zero','perfect'});
          title(sprintf('untrained-sessTr1: %s',regname{r}));
          subplot(242)
          barplot(t1.sessTr,t1.logBayes_sess(:,[2,3]),'subset',t1.sessTr==2);
          hold on;
          drawline(mean(t1.logBayes_sess(t1.sessTr==2,4)),'dir','horz','linewidth',10,'color',[0.8 0.8 0.8]);
          drawline(mean(t1.logBayes_sess(t1.sessTr==2,7)),'dir','horz','linewidth',10,'color',[0.5 0.5 0.5]);
          ylim([0 mean(t1.logBayes_sess(t1.sessTr==2,7))+5]);
          ylabel('logBayes'); set(gca,'XTickLabel',{'zero','perfect'});
          title(sprintf('trained-sessTr2: %s',regname{r}));
          subplot(246)
          barplot(t2.sessTr,t2.logBayes_sess(:,[2,3]),'subset',t2.sessTr==2);
          hold on;
          drawline(mean(t2.logBayes_sess(t2.sessTr==2,4)),'dir','horz','linewidth',10,'color',[0.8 0.8 0.8]);
          drawline(mean(t2.logBayes_sess(t2.sessTr==2,7)),'dir','horz','linewidth',10,'color',[0.5 0.5 0.5]);
          ylim([0 mean(t2.logBayes_sess(t2.sessTr==2,7))+5]);
          ylabel('logBayes'); set(gca,'XTickLabel',{'zero','perfect'});
          title(sprintf('untrained-sessTr2: %s',regname{r}));

          % create new structure
          t=getrow(TT,TT.regType==r & TT.regSide==hemi);
          t.perfZero = bsxfun(@minus,t.logBayes_sess(:,3),t.logBayes_sess(:,2));
          subplot(2,4,[3,4,7,8])
          plt.bar(t.sessTr,t.perfZero,'split',t.seqType,'style',stySeq,'subset',t.sessTr~=3,'leg',{'trained','untrained'});
          title('perfect - zero corr');ylabel('');
          set(gca,'XTickLabel',{'1-2','1-2','2-3','2-3'});
          
          % now here transition 3-4
          figure
          subplot(221)
          barplot(t1.sessTr,t1.logBayes_sess(:,[2,3]),'subset',t1.sessTr==3);
          hold on;
          drawline(mean(t1.logBayes_sess(t1.sessTr==3,6)),'dir','horz','linewidth',10,'color',[0.8 0.8 0.8]);
          drawline(mean(t1.logBayes_sess(t1.sessTr==3,7)),'dir','horz','linewidth',10,'color',[0.5 0.5 0.5]);
          ylim([0 mean(t1.logBayes_sess(t1.sessTr==3,7))+5]);
          ylabel('logBayes'); set(gca,'XTickLabel',{'zero','perfect'});
          title(sprintf('trained-sessTr3: %s',regname{r}));
          subplot(223)
          barplot(t2.sessTr,t2.logBayes_sess(:,[2,3]),'subset',t2.sessTr==3);
          hold on;
          drawline(mean(t2.logBayes_sess(t2.sessTr==3,6)),'dir','horz','linewidth',10,'color',[0.8 0.8 0.8]);
          drawline(mean(t2.logBayes_sess(t2.sessTr==3,7)),'dir','horz','linewidth',10,'color',[0.5 0.5 0.5]);
          ylim([0 mean(t2.logBayes_sess(t2.sessTr==3,7))+5]); 
          ylabel('logBayes');set(gca,'XTickLabel',{'zero','perfect'});
          title(sprintf('untrained-sessTr3: %s',regname{r}));
          subplot(2,2,[2,4])
          plt.bar(t.sessTr,t.perfZero,'split',t.seqType,'style',stySeq,'subset',t.sessTr==3,'leg',{'trained','untrained'});
          title('perfect - zero corr');ylabel('');
          set(gca,'XTickLabel',{'3-4','3-4'});
          
      end
    case 'PCM_stats_noiseCeiling'
      reg = 1:8;
      hemi=1;
      sessName={'1-2','2-3','3-4'};
      parcelType = 'Brodmann';
      vararginoptions(varargin,{'parcelType','reg','hemi','sessName'});
      
      TT=[]; 
      for ss=1:length(sessName)
          T=load(fullfile(stabDir,sprintf('PCM_corr_%s_noiseCeiling_sess%d-sess%d.mat',parcelType,ss,ss+1)));
          % extract the relevant logBayes - seq+corr vs. seq only
          T.logBayes_sess = bsxfun(@minus,T.logBayes,T.logBayes(:,2)); % relative to the session factor
          T.logBayes_sess = T.logBayes_sess(:,2:end);
          T.logBayes_perfZero = bsxfun(@minus,T.logBayes_sess(:,3),T.logBayes_sess(:,2)); % perfect relative to zero
          T.sessTr = ones(size(T.roi))*ss;
          TT=addstruct(TT,T);
      end
      % for now consider models 2 and 3 (zero, perfect correlation)
      % plotting so far models 4 and 7 as noise ceiling (for now check 4)
      for r=reg
          for t=1:3 % all transitions
            T2 = getrow(TT,TT.regType==r&TT.regSide==hemi&TT.sessTr==t);
            tt.logBayes = [T2.logBayes_sess(:,2);T2.logBayes_sess(:,3)];% here just zero and perfect
            tt.seqType  = [T2.seqType;T2.seqType];
            tt.SN       = [T2.SN;T2.SN];
            tt.modelInd = [ones(size(T2.SN));ones(size(T2.SN))*2];
            fprintf('%s: zero vs. perfect transition-%d trained\n',regname{r},t)
            ttestDirect(tt.logBayes,[tt.modelInd tt.SN],2,'paired','subset',tt.seqType==1);
            fprintf('%s: zero vs. perfect transition-%d untrained\n',regname{r},t)
            ttestDirect(tt.logBayes,[tt.modelInd tt.SN],2,'paired','subset',tt.seqType==2);
            fprintf('%s: perfect-zero transition-%d trained vs. untrained\n',regname{r},t)
            ttestDirect(T2.logBayes_perfZero,[T2.seqType T2.SN],2,'paired');
          end
          keyboard;
      end
      
    case 'PCM_bootstrap_corr'
        % bootstrap the correlation values derived from flexible model
        % sample with replacement across subjects
        % useful to determine stability of estimates across subjects, 
        % confidence bounds
        modelType   = 'specific'; % generic or specific
        parcelType  = 'Brodmann';
        reg         = 1:8;
        hemi        = 1;
        nPerm       = 1000; % number of permutations
        perc        = 0.025; % what percentile to take (upp, low)
        quart       = 0.25; % for each sample what proportion to conside (uppAll, lowAll)
        vararginoptions(varargin,{'modelType','parcelType','reg','hemi','nPerm','perc'});
        
        T = load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allSess',parcelType,modelType)));
        sessTr=unique(T.sessTr)';
        nSubj=length(unique(T.SN));
        SS=[];
        for r=reg
            for tr=sessTr
            %    for st=seqType
                    t=T.r_model(T.sessTr==tr & T.regType==r & T.regSide==hemi & T.seqType==1)-T.r_model(T.sessTr==tr & T.regType==r & T.regSide==hemi & T.seqType==2);
                    %t=T.r_model(T.sessTr==tr & T.regType==r & T.regSide==hemi & T.seqType==st);
                    samp=sample_wr(t,nSubj,nPerm);
                    tstat = zeros(nPerm,1);
                    for i=1:nPerm
                        tstat(i)=ttest(samp(:,i),[],1,'onesample');
                    end
                    S.prob      = 1-tcdf(mean(abs(tstat)),25);
                    %S.prob      = length(find(tstat<0))/nPerm;
                    % estimate confidence intervals
                    S.trueMean  = mean(t); % true mean from flexible model fit
                    S.meanAll   = mean(samp,1); % store all mean values
                    S.meanCorr  = mean(S.meanAll);
                    sSamp       = sort(S.meanAll); % sort to determine percentiles
                    sAll        = sort(samp,1); % sort all of the data
                    S.upp       = sSamp(perc*nPerm);
                    S.low       = sSamp(end-perc*nPerm);
                    S.uppAll    = sAll(end-ceil(quart*size(sAll,1)),:);
                    S.lowAll    = sAll(floor(quart*size(sAll,1)),:);
                    S.sessTr    = tr;
                %    S.seqType   = st;
                    S.roi       = r;
                    S.hemi      = hemi;
                    SS=addstruct(SS,S);
                end
           % end
        end
        % here save the new structure
        save(fullfile(stabDir,sprintf('PCM_corr_bootstrap_stats_%s_%s',parcelType,modelType)),'-struct','SS');
    case 'PLOT_bootstrap_corr_old'
        % plotting distributions
        modelType = 'specific';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        vararginoptions(varargin,{'modelType','parcelType','reg','hemi'});
        
        T = load(fullfile(stabDir,sprintf('PCM_corr_bootstrap_%s_%s',parcelType,modelType)));
        sessTr = unique(T.sessTr)';
        for r=reg
            figure
            for ss=sessTr
                t = getrow(T,T.sessTr==ss & T.roi==r & T.hemi==hemi);
                subplot(1,length(sessTr),ss)
                histogram(t.meanAll(2,:)); % untrained in blue
                hold on;
                histogram(t.meanAll(1,:)); % trained in redish
                drawline([t.upp(1),t.low(1)],'dir','vert','color',[1 0 0],'linestyle','--'); % upper / lower bound
                drawline([t.upp(2),t.low(2)],'dir','vert','color',[0 0 1],'linestyle','--');
                drawline(t.meanCorr(1),'dir','vert','color',[1 0 0]);
                drawline(t.meanCorr(2),'dir','vert','color',[0 0 1]);
                title(sprintf('corr for %s - sess:%d-%d',regname{r},ss,ss+1));
                set(gca,'fontsize',14);                
            end
            figure
           % subplot(121)
            t = getrow(T,T.roi==r & T.hemi==hemi);
           % plt.line(t.sessTr,abs(mean(t.meanAll-t.uppAll,2)),'split',t.seqType,'style',stySeq,'leg',{'trained','untrained'});
            plt.line(t.sessTr,abs(t.meanCorr-t.low),'split',t.seqType,'style',stySeq,'leg',{'trained','untrained'});
        %    title(sprintf('%s lower bound',regname{r}));
            hold on;
          %  plt.line(t.sessTr,abs(mean(t.meanAll-t.lowAll,2)),'split',t.seqType,'style',stySeq,'leg',{'trained','untrained'});
            plt.line(t.sessTr,abs(t.meanCorr-t.upp),'split',t.seqType,'style',stySeq2,'leg',{'trained','untrained'})
         %   title(sprintf('%s upper bound',regname{r}));
            title(sprintf('%s - upper ceiling in dash, lower solid',regname{r}));
            xlabel('Session transition');
            ylabel('Deviation in correlation mean - lower/upper bound');
            plt.match('y');  
            
            figure
            subplot(121)
            t = getrow(T,T.seqType==1 & T.roi==r & T.hemi==hemi);
            histogram(t.meanAll(1,:));
            hold on;
            histogram(t.meanAll(2,:));
            histogram(t.meanAll(3,:));
            drawline(t.meanCorr,'dir','vert');
            title(sprintf('corr for trained across sess - %s',regname{r}));
            subplot(122)
            t = getrow(T,T.seqType==2 & T.roi==r & T.hemi==hemi);
            histogram(t.meanAll(1,:));
            hold on;
            histogram(t.meanAll(2,:));
            histogram(t.meanAll(3,:));
            drawline(t.meanCorr,'dir','vert');
            title(sprintf('corr for untrained across sess - %s',regname{r}));
        end
    case 'PLOT_bootstrap_corr_stats'
        % here plot statistics derived from boostrapping
        modelType = 'specific';
        parcelType = 'Brodmann';
        reg = [1:3,7,8];
        vararginoptions(varargin,{'modelType','parcelType','reg','hemi'});
        
        T = load(fullfile(stabDir,sprintf('PCM_corr_bootstrap_stats_%s_%s',parcelType,modelType)));
        figure
        for r=1:numel(reg)
            subplot(1,numel(reg),r)
            style.use('gray');
            plt.bar(T.sessTr,T.prob,'subset',T.roi==reg(r));
            hold on; drawline(0.05,'dir','horz');
            if r==1
                ylabel('probability');
            else
                ylabel('');
            end
            xlabel('session transition');
            title(regname{reg(r)});
        end
    case 'PCM_bootstrap_log'
        % boostrap the log factors from the correlation model
        % from PCM_allCorr model
        modelType   = 'specific'; % generic or specific
        parcelType  = 'Brodmann';
        reg         = 1:8;
        hemi        = 1;
        nPerm       = 1000; % number of permutations
        perc        = 0.025; % what percentile to take (upp, low)
        vararginoptions(varargin,{'modelType','parcelType','reg','hemi','nPerm','perc'});
        
        T = load(fullfile(stabDir,sprintf('PCM_allCorr_%s_%s_allSess',parcelType,modelType)));
        sessTr=unique(T.sessTr)';
        seqType=unique(T.seqType)';
        nSubj=length(unique(T.SN));
        SS=[];
        for r=reg
            for tr=sessTr
                for st=seqType
                    t=getrow(T,T.sessTr==tr & T.regType==r & T.regSide==hemi & T.seqType==st);
                    % permute to estimate the correlation with highest logBayes 
                    [logB,numCorr] = max(t.bayesEst,[],2); % here consider just the best model
                    permLogB     = sample_wr(logB,nSubj,nPerm);
                    permCorr     = sample_wr(numCorr,nSubj,nPerm)./10; % make into correlation
                    % here permute the estimate of perfect correlation
                    permMaxB    = sample_wr(t.bayesEst(:,end),nSubj,nPerm);
                    % extract info
                    S.logAll        = mean(permLogB,1);
                    S.corrAll       = mean(permCorr,1);
                    S.logMaxAll     = mean(permMaxB,1);
                    S.logMean       = mean(S.logAll);
                    S.corrMean      = mean(S.corrAll);
                    S.logMaxMean    = mean(S.logMaxAll);
                    sLog                = sort(S.logAll);
                    sCorr               = sort(S.corrAll);
                    sLogMax             = sort(S.logMaxAll);
                    S.lowLog            = sLog(perc*nPerm);
                    S.uppLog            = sLog(end-perc*nPerm);
                    S.lowCorr           = sCorr(perc*nPerm);
                    S.uppCorr           = sCorr(end-perc*nPerm);
                    S.lowLogMax         = sLogMax(perc*nPerm);
                    S.uppLogMax         = sLogMax(end-perc*nPerm);
                    S.sessTr            = tr;
                    S.seqType           = st;
                    S.roi               = r;
                    S.hemi              = hemi;
                    SS = addstruct(SS,S);
                end
            end
        end
        save(fullfile(stabDir,sprintf('PCM_log_bootstrap_%s_%s',parcelType,modelType)),'-struct','SS');
    case 'PLOT_bootstrap_log'
        modelType = 'specific';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        vararginoptions(varargin,{'modelType','parcelType','reg','hemi'});
        
        T = load(fullfile(stabDir,sprintf('PCM_log_bootstrap_%s_%s',parcelType,modelType)));
        sessTr = unique(T.sessTr)';
        for r=reg
            figure
            for ss=sessTr
                t = getrow(T,T.sessTr==ss & T.roi==r & T.hemi==hemi);
                subplot(1,length(sessTr),ss)
                histogram(t.logMaxAll(2,:)); % untrained in blue
                hold on;
                histogram(t.logMaxAll(1,:)); % trained in redish
                drawline([t.uppLogMax(1),t.lowLogMax(1)],'dir','vert','color',[1 0 0],'linestyle','--'); % upper / lower bound
                drawline([t.uppLogMax(2),t.lowLogMax(2)],'dir','vert','color',[0 0 1],'linestyle','--');
                drawline(t.logMaxMean(1),'dir','vert','color',[1 0 0]);
                drawline(t.logMaxMean(2),'dir','vert','color',[0 0 1]);
                title(sprintf('perfect corr model for %s - sess:%d-%d',regname{r},ss,ss+1));
                xlabel('logBayes');
                set(gca,'fontsize',14);
            end
            figure
            subplot(121)
            t = getrow(T,T.seqType==1 & T.roi==r & T.hemi==hemi);
            histogram(t.logMaxAll(1,:));
            hold on;
            histogram(t.logMaxAll(2,:));
            histogram(t.logMaxAll(3,:));
            drawline(t.logMaxMean,'dir','vert');
            title(sprintf('perfect corr for trained - %s',regname{r}));
            set(gca,'fontsize',14);
            subplot(122)
            t = getrow(T,T.seqType==2 & T.roi==r & T.hemi==hemi);
            histogram(t.logMaxMean(1,:));
            hold on;
            histogram(t.logMaxAll(2,:));
            histogram(t.logMaxAll(3,:));
            drawline(t.logMaxMean,'dir','vert');
            title(sprintf('perfect corr for untrained - %s',regname{r}));
            set(gca,'fontsize',14);
        end    
    case 'PCM_plot_allCorr'
      reg = 1:8;
      hemi=1;
      sessName={'1-2','2-3','3-4'};
      modelType='specific'; % generic or specific
      parcelType = 'Brodmann';
      metric = 'bayesEst'; % bayesEst or individ_bayes
      vararginoptions(varargin,{'modelType','metric','parcelType','reg'});
      TT=[];
      for ss=1:length(sessName)
          T=load(fullfile(stabDir,sprintf('PCM_allCorr_%s_%s_sess%d-sess%d.mat',parcelType,modelType,ss,ss+1)));
          T.sessTr = ones(size(T.roi))*ss;
          TT=addstruct(TT,T);
      end
      
      TT.newMetric = bsxfun(@minus,TT.likelihood,TT.likelihood(:,1));
      for r=reg
       %   figure
          t1 = getrow(TT,TT.regType==r&TT.regSide==hemi&TT.seqType==1);
          s1=mean(t1.(metric)(t1.sessTr==1,:))';
          s2=mean(t1.(metric)(t1.sessTr==2,:))';
          s3=mean(t1.(metric)(t1.sessTr==3,:))';
          [r1,tm1]=max(s1);
          st1=s1-r1;
          [r1,tm2]=max(s2);
          st2=s2-r1;
          [r1,tm3]=max(s3);
          st3=s3-r1;
       %   st1=s1;
       %   st2=s2;
       %   st3=s3;
          t2 = getrow(TT,TT.regType==r&TT.regSide==hemi&TT.seqType==2);
          s1=mean(t2.(metric)(t2.sessTr==1,:))';
          s2=mean(t2.(metric)(t2.sessTr==2,:))';
          s3=mean(t2.(metric)(t2.sessTr==3,:))';
          [r1,um1]=max(s1);
          su1=s1-r1;
          [r1,um2]=max(s2);
          su2=s2-r1;
          [r1,um3]=max(s3);
          su3=s3-r1;
       %   su1=s1;
       %   su2=s2;
       %   su3=s3;
%           figure
%           subplot(121)
%           ribbon([st1 st2 st3]);
%           title(sprintf('%s - trained',regname{r}));  
%           subplot(122)
%           ribbon([su1 su2 su3]);
%           title(sprintf('%s - untrained',regname{r}));  
          % determine maxima to plot over
          [~,mt1]=max(st1(11:-1:1));
          [~,mt2]=max(st2(11:-1:1));
          [~,mt3]=max(st3(11:-1:1));
          [~,mu1]=max(su1(11:-1:1));
          [~,mu2]=max(su2(11:-1:1));
          [~,mu3]=max(su3(11:-1:1));
          
%           figure
%           subplot(121)
%           imagesc([st1(11:-1:1) st2(11:-1:1) st3(11:-1:1)]);
%           hold on;
%           scatter(1:3,[mt1,mt2,mt3],'o','filled','MarkerFaceColor',[0 0 0]);
%           plot(1:3,[mt1,mt2,mt3],'-k','LineWidth',2);
%         %  caxis([-10 0]);
%           title(sprintf('%s - trained',regname{r}));
%           colormap('bone');colorbar;
%           subplot(122)
%           imagesc([su1(11:-1:1) su2(11:-1:1) su3(11:-1:1)]);
%           hold on;
%           scatter(1:3,[mu1,mu2,mu3],'o','filled','MarkerFaceColor',[0 0 0]);
%           plot(1:3,[mu1,mu2,mu3],'-k','LineWidth',2);
%         %  caxis([-10 0]);
%           title(sprintf('%s - untrained',regname{r}));
%           colormap('bone');colorbar;
%  

          figure
          subplot(121)
          lineplot([1:11 1:11]',[st1;st2],'split',[ones(11,1);ones(11,1)*2],...
              'style_symbols4*2','linewidth',3,'linecolor',{[0.1 0.1 0.1],[0.7 0.7 0.7]},'markertype',{'o','v'}); 
   
          hold on; drawline(0,'dir','horz','linestyle','--');
          title(sprintf('trained %s',regname{r})); ylabel(sprintf('%s relative to best correlation',metric)); xlabel('correlation'); legend('1-2','2-3');
          subplot(122)
          lineplot([1:11 1:11]',[su1;su2],'split',[ones(11,1);ones(11,1)*2],'style_symbols4*2','linewidth',3);
          hold on; drawline(0,'dir','horz','linestyle','--');
          title(sprintf('untrained %s',regname{r})); 
          % create a new structure
          t = getrow(TT,TT.regType==r & TT.regSide==hemi);
         % [mt,~] =  max(t.(metric),[],2);
         % t.(metric) = bsxfun(@minus,t.(metric),mt);
          n.(metric)=t.(metric)(:);
          n.corr=[];
          for i=1:size(t.(metric),2)
            n.corr = [n.corr; ones(size(t.SN,1),1)*i];
          end
          n.seqType = repmat(t.seqType,11,1);
          n.sessTr  = repmat(t.sessTr,11,1);
          
          % plot of sess 1-3
          figure
          subplot(221)
          plt.line(n.corr,n.(metric),'split',n.seqType,'subset',n.sessTr==1,'style',stySeq,'leg',{'trained','untrained'},'leglocation','southeast');
          hold on;
          ylabel(metric);
          drawline(tm1,'dir','vert','color',[1 0 0],'linestyle','--');
          drawline(um1,'dir','vert','color',[0 0 1],'linestyle','--');
          drawline(0,'dir','horz');
          ylabel(metric);
          title(sprintf('%s transition 1',regname{r}));
          subplot(222)
          plt.line(n.corr,n.(metric),'split',n.seqType,'subset',n.sessTr==2,'style',stySeq,'leg',{'trained','untrained'},'leglocation','southeast');
          hold on;
          ylabel(metric);
          drawline(tm2,'dir','vert','color',[1 0 0],'linestyle','--');
          drawline(um2,'dir','vert','color',[0 0 1],'linestyle','--');
          drawline(0,'dir','horz');
          ylabel('');
          title(sprintf('%s transition 2',regname{r}));
          style.use('Sess');
          subplot(223)
          plt.line(n.corr,n.(metric),'split',n.sessTr,'subset',n.sessTr<3 & n.seqType==1,'leg',{'1-2','2-3'},'leglocation','southeast');
          hold on;
          drawline(tm1,'dir','vert','linestyle','--');
          drawline(tm2,'dir','vert','color',[0.7 0.7 0.7],'linestyle','--');
          drawline(0,'dir','horz');
          title('trained transition');
          ylabel(metric);
          subplot(224)
          plt.line(n.corr,n.(metric),'split',n.sessTr,'subset',n.sessTr<3 & n.seqType==2,'leg',{'1-2','2-3'},'leglocation','southeast');
          hold on;
          drawline(um1,'dir','vert','linestyle','--');
          drawline(um2,'dir','vert','color',[0.7 0.7 0.7],'linestyle','--');
          drawline(0,'dir','horz');
          title('untrained transition'); ylabel('');
          
          figure
          plt.line(n.corr,n.(metric),'split',n.seqType,'subset',n.sessTr==3,'style',stySeq,'leg',{'trained','untrained'},'leglocation','southeast');
          hold on;
          ylabel(metric);
          drawline(tm3,'dir','vert','color',[1 0 0],'linestyle','--');
          drawline(um3,'dir','vert','color',[0 0 1],'linestyle','--');
          drawline(0,'dir','horz');
          ylabel(metric);
          title(sprintf('Session 3-4 %s',regname{r}));          
          keyboard;
      end
    case 'PLOT_allCorr_normalised'
      % here normalise by the top fit
      reg = 1:8;
      hemi=1;
      sessName={'1-2','2-3','3-4'};
      modelType='specific'; % generic or specific
      parcelType = 'Brodmann';
      metric = 'bayesEst'; % bayesEst or individ_bayes
      vararginoptions(varargin,{'modelType','metric','parcelType','reg'});
      TT=[];
      for ss=1:length(sessName)
          T=load(fullfile(stabDir,sprintf('PCM_allCorr_%s_%s_sess%d-sess%d.mat',parcelType,modelType,ss,ss+1)));
          T.sessTr = ones(size(T.roi))*ss;
          TT=addstruct(TT,T);
      end 
      for r=reg
          t = getrow(TT,TT.regType==r&TT.regSide==hemi);
          [~,maxCorr]=max(t.(metric),[],2);
          t.(metric) = bsxfun(@minus,t.(metric),max(t.(metric),[],2));
          % reshape
          C.(metric)    = t.(metric)(:);
          C.SN          = repmat(t.SN,11,1);
          C.seqType     = repmat(t.seqType,11,1);
          C.corr        = kron(1:11,ones(1,size(t.SN,1)))';
          C.sessTr      = repmat(t.sessTr,11,1);
%           figure
%           style.use('Seq');
%           subplot(221)
%           plt.line(C.corr,C.(metric),'split',C.seqType,'subset',C.sessTr==1,'leg',{'trained','untrained'},'leglocation','southeast');
%           hold on; drawline(0,'dir','horz','linestyle','--');
%           ylabel('logBayes normalised'); title(sprintf('%s transition 1',regname{r}));
%           subplot(222)
%           plt.line(C.corr,C.(metric),'split',C.seqType,'subset',C.sessTr==2,'leg',{'trained','untrained'},'leglocation','southeast');
%           hold on; drawline(0,'dir','horz','linestyle','--');
%           ylabel(''); title(sprintf('%s transition 2',regname{r}));
%           subplot(223)
%           style.use('Trained');
%           plt.line(C.corr,C.(metric),'split',C.sessTr,'subset',C.seqType==1&ismember(C.sessTr,[1,2]),'leg',{'1-2','2-3'},'leglocation','southeast');
%           hold on; drawline(0,'dir','horz','linestyle','--');
%           ylabel('logBayes factor normalised'); xlabel('Correlation model');
%           title(sprintf('%s trained transitions',regname{r}));
%           subplot(224)
%           style.use('Untrained');
%           plt.line(C.corr,C.(metric),'split',C.sessTr,'subset',C.seqType==2&ismember(C.sessTr,[1,2]),'leg',{'1-2','2-3'},'leglocation','southeast');
%           hold on; drawline(0,'dir','horz','linestyle','--');
%           ylabel(''); xlabel('Correlation model');
%           title(sprintf('%s untrained transitions',regname{r}));
%           
          figure
          style.use('Seq');
          plt.bar(C.sessTr,maxCorr,'split',C.seqType,'leg',{'Trained','Untrained'},'plotfcn','mean');
          title(sprintf('%s - correlation model with highest logBayes',regname{r}));
          ylabel('Correlation');
      end
    case 'PCM_stats_allCorr'
      reg = 1:8;
      hemi=1;
      sessName={'1-2','2-3','3-4'};
      modelType='specific_noSess'; % generic or specific
      parcelType = 'Brodmann';
      metric = 'bayesEst'; % bayesEst or individ_bayes
      vararginoptions(varargin,{'modelType','metric','parcelType','reg'});
      TT=[];
      for ss=1:length(sessName)
          T=load(fullfile(stabDir,sprintf('PCM_allCorr_%s_%s_sess%d-sess%d.mat',parcelType,modelType,ss,ss+1)));
          T.sessTr = ones(size(T.roi))*ss;
          TT=addstruct(TT,T);
      end
      
      for r=reg
          t = getrow(TT,TT.regType==r & TT.regSide==hemi);
          n.(metric)=t.(metric)(:);
          n.corr=[];
          for i=1:size(t.(metric),2)
            n.corr = [n.corr; ones(size(t.SN,1),1)*i];
          end
          n.seqType = repmat(t.seqType,11,1);
          n.sessTr  = repmat(t.sessTr,11,1);
          n.SN      = repmat(t.SN,11,1);
          
          for ss=1:3 % trained vs. untrained for each session transition
              fprintf('%s - trained vs. untrained transition: %d\n',regname{r},ss);
              ttestDirect(n.(metric),[n.seqType n.SN],2,'paired','subset',n.sessTr==ss,'split',n.corr);
          end
          fprintf('%s - trained trans 1-2\n',regname{r});
          ttestDirect(n.(metric),[n.sessTr n.SN],2,'paired','subset',n.seqType==1 & ismember(n.sessTr,[1,2]),'split',n.corr);
          fprintf('%s - untrained trans 1-2\n',regname{r});
          ttestDirect(n.(metric),[n.sessTr n.SN],2,'paired','subset',n.seqType==2 & ismember(n.sessTr,[1,2]),'split',n.corr);
          keyboard;
      end
      
    case 'interindivid'
        % explore the relationship between ID in learning and PCM
        % correlation
        reg = 1:8;
        hemi=1;
        modelType='specific'; % generic or specific
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','sessN','metric','modelType','parcelType','hemi'});

        TT=load(fullfile(stabDir,sprintf('PCM_stability_%s_%s_allSess.mat',parcelType,modelType)));
        B = load(fullfile(behavDir,'exponentialFits'));
        B=getrow(B,B.sn~=4); % exclude that subject (not used in PCM);
        
        for r=reg
            figure(1)
            t=getrow(TT,TT.regType==r & TT.regSide==hemi);
            subplot(2,max(reg),r)
            [i,j]=corr(t.r_model(t.sessTr==3&t.seqType==1),B.x2_mean);
            plt.scatter(t.r_model(t.sessTr==3&t.seqType==1),B.x2_mean);
            title(sprintf('%s T speed corr: %2.1f, p=%1.3f',regname{r},i,j));
            subplot(2,max(reg),r+max(reg))
            [i,j]=corr(t.r_model(t.sessTr==3&t.seqType==2),B.x2_mean);
            plt.scatter(t.r_model(t.sessTr==3&t.seqType==2),B.x2_mean);
            title(sprintf('%s U speed corr: %2.1f, p=%1.3f',regname{r},i,j));
            
            figure(2)
            subplot(1,max(reg),r)
            style.use('Trained');
            plt.scatter(t.r_model(t.sessTr==3&t.seqType==1),B.x2_mean);
            hold on;
            style.use('Untrained');
            plt.scatter(t.r_model(t.sessTr==3&t.seqType==2),B.x2_mean);
            title(sprintf('%s - speed',regname{r}));
            
            figure(3)
            subplot(2,max(reg),r)
            [i,j]=corr(t.r_model(t.sessTr==1&t.seqType==1),B.x0_mean);
            plt.scatter(t.r_model(t.sessTr==1&t.seqType==1),B.x0_mean);
            title(sprintf('%s T start corr: %2.1f, p=%1.3f',regname{r},i,j));
            subplot(2,max(reg),r+max(reg))
            [i,j]=corr(t.r_model(t.sessTr==1&t.seqType==2),B.x0_mean);
            plt.scatter(t.r_model(t.sessTr==1&t.seqType==2),B.x0_mean);
            title(sprintf('%s U start corr: %2.1f, p=%1.3f',regname{r},i,j));   
            
            figure(4)
            subplot(2,max(reg),r)
            [i,j]=corr(t.r_model(t.sessTr==1&t.seqType==1),B.x1_mean);
            plt.scatter(t.r_model(t.sessTr==1&t.seqType==1),B.x1_mean);
            title(sprintf('%s T curve corr: %2.1f, p=%1.3f',regname{r},i,j));
            subplot(2,max(reg),r+max(reg))
            [i,j]=corr(t.r_model(t.sessTr==1&t.seqType==2),B.x1_mean);
            plt.scatter(t.r_model(t.sessTr==1&t.seqType==2),B.x1_mean);
            title(sprintf('%s U curve corr: %2.1f, p=%1.3f',regname{r},i,j)); 
            
            figure(5)
            subplot(3,max(reg),r)
            [i,j]=corr(t.r_model(t.sessTr==1&t.seqType==1)-t.r_model(t.sessTr==1&t.seqType==2),B.x0_mean);
            plt.scatter(t.r_model(t.sessTr==1&t.seqType==1)-t.r_model(t.sessTr==1&t.seqType==2),B.x0_mean);
            title(sprintf('%s - diff T-U start: %2.1f, p=%1.3f',regname{r},i,j));
            subplot(3,max(reg),r+max(reg))
            [i,j]=corr(t.r_model(t.sessTr==1&t.seqType==1)-t.r_model(t.sessTr==1&t.seqType==2),B.x1_mean);
            plt.scatter(t.r_model(t.sessTr==1&t.seqType==1)-t.r_model(t.sessTr==1&t.seqType==2),B.x1_mean);
            title(sprintf('%s - diff T-U curve: %2.1f, p=%1.3f',regname{r},i,j));
            subplot(3,max(reg),r+2*max(reg))
            [i,j]=corr(t.r_model(t.sessTr==3&t.seqType==1)-t.r_model(t.sessTr==3&t.seqType==2),B.x2_mean);
            plt.scatter(t.r_model(t.sessTr==3&t.seqType==1)-t.r_model(t.sessTr==3&t.seqType==2),B.x2_mean);
            title(sprintf('%s - diff T-U curve: %2.1f, p=%1.3f',regname{r},i,j));
            
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
        
    case 'ratio_sess3-4'
        sn=[5:9,11:31];
        roi=1:8;
        hemi=1:2;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','parcelType','roi'});
        
        P=load(fullfile(distPscDir,sprintf('psc_%s_ROI',parcelType)));
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)));

        P=getrow(P,ismember(P.sn,sn));
        D=getrow(D,ismember(D.sn,sn));
        distType={'dist_train','dist_untrain'};
        NN=[];
        for h=hemi
            for r=roi
                for st=1:2
                    p3=getrow(P,P.regType==r & P.regSide==h & P.sessN==3 & P.seqType==st);
                    p4=getrow(P,P.regType==r & P.regSide==h & P.sessN==4 & P.seqType==st);
                    d3=getrow(D,D.regType==r & D.regSide==h & D.sessN==3);
                    d4=getrow(D,D.regType==r & D.regSide==h & D.sessN==4);
                    N.psc3          = p3.psc;
                    N.psc4          = p4.psc;
                    N.psc_ratio     = p4.psc./p3.psc; % session 4 / session 3 (>1 if increase)
                    N.dist3         = d3.(distType{st});
                    N.dist4         = d4.(distType{st});
                    N.predDist      = N.psc_ratio.*d3.(distType{st});
                    N.dist_ratio    = d4.(distType{st})./d3.(distType{st});
                    N.seqType       = ones(size(N.psc_ratio))*st;
                    N.roi           = d4.roi;
                    N.regType       = d4.regType;
                    N.regSide       = d4.regSide;
                    N.sn            = d4.sn;
                    NN=addstruct(NN,N);
                end
            end
        end
        save(fullfile(distPscDir,sprintf('session3-4_dist_psc_%s',parcelType)),'-struct','NN');
    case 'plot_ratio_sess3-4'
        % plot in bargraphs (appropriate for Brodmann) % for 162 project
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','parcelType'});
        
        T = load(fullfile(distPscDir,sprintf('session3-4_dist_psc_%s',parcelType)));
        for r=roi
            figure
            style.use('Sess');
            subplot(221)
            plt.bar(T.seqType,[T.psc3 T.psc4],'subset',T.regType==r&T.regSide==hemi,'leg',{'sess-3','sess-4'});
            xlabel('seqType'); ylabel('psc');
            title(sprintf('%s - percent signal change',regname{r}));
            
            subplot(222)
            plt.bar(T.seqType,[T.dist3 T.dist4],'subset',T.regType==r&T.regSide==hemi,'leg',{'sess-3','sess-4'});
            drawline(mean(T.predDist(T.regType==r&T.regSide==hemi&T.seqType==1)),'dir','horz','lim',[1.7 2.7],'linestyle','--');
            drawline(mean(T.predDist(T.regType==r&T.regSide==hemi&T.seqType==2)),'dir','horz','lim',[4.9 5.9],'linestyle','--');
            xlabel('seqType'); ylabel('dist');
            title(sprintf('%s - distance',regname{r}));
            
            subplot(223)
            plt.bar(T.seqType,T.psc_ratio,'subset',T.regType==r & T.regSide==hemi);
            drawline(1,'dir','horz');
            xlabel('seqType'); ylabel('ratio dist');
          
            subplot(224)
            plt.bar(T.seqType,T.dist_ratio,'subset',T.regType==r&T.regSide==hemi);
            hold on;
            drawline(1,'dir','horz');
              xlabel('seqType'); ylabel('ratio dist');
            title(sprintf('%s - distance',regname{r}));

        end
    case 'stats_ratio_sess3-4'
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','parcelType'});
        
        T = load(fullfile(distPscDir,sprintf('session3-4_dist_psc_%s',parcelType)));
        
        for r=roi
            t = getrow(T,T.regType==r & T.regSide==hemi);
            fprintf('\n\n%s t-test distance vs. predicted - trained\n',regname_cortex{r});
            ttestDirect([t.predDist; t.dist4],[[ones(size(t.predDist));ones(size(t.predDist))*2] [t.sn;t.sn]],2,'paired','subset',[t.seqType;t.seqType]==1);
            fprintf('%s t-test distance vs. predicted - untrained\n',regname_cortex{r});
            ttestDirect([t.predDist; t.dist4],[[ones(size(t.predDist));ones(size(t.predDist))*2] [t.sn;t.sn]],2,'paired','subset',[t.seqType;t.seqType]==2);
            fprintf('%s t-test trained vs. untrained dist ratio\n',regname_cortex{r});
            ttestDirect(t.dist_ratio,[t.seqType t.sn],2,'paired');
            fprintf('%s t-test trained vs. untrained psc ratio\n',regname_cortex{r});
            ttestDirect(t.psc_ratio,[t.seqType t.sn],2,'paired');
            fprintf('%s t-test trained vs. untrained pred vs. ratio\n',regname_cortex{r});
            ttestDirect(t.dist4-t.predDist,[t.seqType t.sn],2,'paired');
            fprintf('%s t-test trained vs. untrained distance 3\n',regname_cortex{r});
            ttestDirect(t.dist3,[t.seqType t.sn],2,'paired');
            fprintf('%s t-test trained vs. untrained distance 4\n',regname_cortex{r});
            ttestDirect(t.dist4,[t.seqType t.sn],2,'paired');
            
        end
    case 'sess4_beh_dist'
        sn=[5:9,11:31];
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','parcelType','roi'});
        
        P=load(fullfile(distPscDir,sprintf('psc_%s_ROI',parcelType)));
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)));
        B=load(fullfile(behavDir,'alldata'));
        B=getrow(B,ismember(B.SN,sn) & B.blockType==9);
        [MT,sn]=pivottable(B.SN,B.seqType,B.MT,'median');
        clear B;
        B.sn = [sn;sn];
        B.MT = [MT(:,1);MT(:,2)];
        B.seqType = [ones(size(sn));ones(size(sn))*2];
        P=getrow(P,ismember(P.sn,sn) & P.sessN==4);
        D=getrow(D,ismember(D.sn,sn) & D.sessN==4);
        DD.dist = [D.dist_train;D.dist_untrain];
        DD.seqType = [ones(size(D.dist_train));ones(size(D.dist_train))*2];
        DD.sn   = [D.sn;D.sn];
        DD.regType = [D.regType;D.regType];
        DD.regSide = [D.regSide; D.regSide];
        DD.roi = [D.roi; D.roi];
        D=DD;
        clear DD;
        for h=hemi
            for r=roi
                p=getrow(P,P.regType==r & P.regSide==h);
                d=getrow(D,D.regType==r & D.regSide==h);
                figure
                subplot(121)
                corrP=plt.scatter(B.MT(B.seqType==2)-B.MT(B.seqType==1),p.psc(p.seqType==1)-p.psc(p.seqType==2));
                drawline(0,'dir','horz');
                xlabel('behavioural advantage');
                ylabel('percent signal difference (trained - untrained)');
                title(sprintf('%s - psc: %1.3f corr',regname{r},corrP));
                subplot(122)
                corrD=plt.scatter(B.MT(B.seqType==2)-B.MT(B.seqType==1),d.dist(d.seqType==1)-d.dist(d.seqType==2));
                drawline(0,'dir','horz');
                xlabel('behavioural advantage');
                ylabel('distance difference (trained - untrained)');
                title(sprintf('%s - dist: %1.3f corr',regname{r},corrD));
            end
        end
    case 'run_job'
        for t=1:3
            sml1_imana_stability('PCM_noiseCeiling','sessN',[t,t+1],'parcelType','BG-striatum','reg',[1:2]);
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
    % Run group fit in 2 stages - first NR, finish with minimize (any
    % further improvement)
    % [T,theta_hat,~,~] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize');
     %fprintf('Group fit with minimize algorithm (round 1) done.\n');
     [~,theta_hat,~,~] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR');
     fprintf('Group fit with NR algorithm done.\n');
     [T,theta_hat,~,~] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize','theta0',theta_hat);
     fprintf('Group fit with minimize algorithm done.\n');
     % previously only one round of NR fit - minimize is slower, so here
     % only used to so that parameter estimates more likely converge
     % previous implementation:
   % [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
   % fprintf('Group fit with %s algorithm done.\n',algorithm);
 
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
function M = pcm_stabilityModel_generic_noSess

    % for sequence-specific modelling - one parameter for all seq           
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12,1) = zeros(12);

    % Model 2: Sequence specific
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'Seq';
    M{2}.Ac(:,1:6,1)  = [A;zeros(6)];      % Unique sess1 sequence patterns
    M{2}.Ac(:,7:12,2) = [zeros(6);A];     % Unique sess2 sequence pattterns
    
    % Model 3: Sequence specific + correlation in sequences
    M{3}.type         = 'feature';
    M{3}.numGparams   = 3;
    M{3}.name         = 'Seq+Corr';
    M{3}.Ac(:,1:6,1)  = [A;zeros(6)];       % Unique sess1 sequence patterns
    M{3}.Ac(:,7:12,2) = [zeros(6);A];       % Unique sess2 sequence pattterns
    M{3}.Ac(:,1:6,3)  = [zeros(6);A];     % Correlation sess1-sess2
    
    % Model 4: Sequence specific + PERFECT correlation in sequences
    M{4}.type         = 'feature';
    M{4}.numGparams   = 2;
    M{4}.name         = 'Seq+PerfectCorr';
    M{4}.Ac(:,1:6,1)  = [A;zeros(6)];       % Unique sess1 sequence patterns
    M{4}.Ac(:,1:6,2)  = [zeros(6);A];        % Identical sess2 sequence pattterns

    
end
function M = pcm_stabilityModel_specific_noSess 
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: Sequence specific
    M{2}.type       = 'feature';
    M{2}.numGparams = 12;
    M{2}.name       = 'Session+Seq';
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{2}.Ac(:,1:6,i)     = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{2}.Ac(:,7:12,6+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    end;

    % Model 3: Session + sequence specific + correlation in repetition
    M{3}.type         = 'feature';
    M{3}.numGparams   = 18;
    M{3}.name         = 'Seq+Corr';
    % for sequence-specific modelling- one parameter per session
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,1:6,i)     = [A;zeros(6)];     % Unique sess1 sequence patterns
        M{3}.Ac(:,7:12,6+i)  = [zeros(6);A];     % Unique sess2 sequence pattterns
        M{3}.Ac(:,1:6,12+i)  = [zeros(6);A];     % Correlation sess1-sess2
    end;

    % Model 4: Session + sequence specific + PERFECT correlation in session
    M{4}.type         = 'feature';
    M{4}.numGparams   = 12;
    M{4}.name         = 'Seq+PerfectCorr';
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,1:6,i)     = [A;zeros(6)];       % Unique sess1 sequence patterns
        M{4}.Ac(:,1:6,6+i)   = [zeros(6);A];       % Same sess2 sequence pattterns
    end;
    
end

function M = pcm_stabilityModel_specificCorrelation(r)
% create a model with a specific correlation r
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

    % Model 4: Session + sequence specific + specific correlation in repetition
    % determine needed theta
    [th2,th3] = det_theta_corr(r);
    M{4}.type         = 'feature';
    M{4}.numGparams   = 14;
    M{4}.name         = 'Session+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per session
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,3:8,2+i)   = [A;zeros(6)];      % Unique sess1 sequence patterns
        M{4}.Ac(:,9:14,8+i)  = [zeros(6);A*th2];  % Unique sess2 sequence pattterns
        M{4}.Ac(:,3:8,8+i)   = [zeros(6);A*th3];  % Correlation sess1-sess2
    end;
end
function M = pcm_stabilityModel_genericCorrelation(r)
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
    [th2,th3] = det_theta_corr(r);
    M{4}.type         = 'feature';
    M{4}.numGparams   = 4;
    M{4}.name         = 'Session+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    M{4}.Ac(:,3:8,3)  = [A;zeros(6)];       % Unique sess1 sequence patterns
    M{4}.Ac(:,9:14,4) = [zeros(6);A*th2];       % Unique sess2 sequence pattterns
    M{4}.Ac(:,3:8,4)  = [zeros(6);A*th3];       % Correlation sess1-sess2
end

function M = pcm_stability_specificCorr_noSess
% specific correlation models with no explicit modelling of session
corrS = 0:0.1:1;
for c=1:length(corrS)
    [th2,th3] = det_theta_corr(corrS(c));
    M{c}.type='feature';
    M{c}.numGparams = 12;
    M{c}.name = sprintf('SpecCorr_%1.1f',corrS(c));
    % for sequence-specific modelling- one parameter per session
     for i=1:6
         A=zeros(6);
         A(i,i)=1;
         M{c}.Ac(:,1:6,i)       = [A;zeros(6)];      % Unique sess1 sequence patterns
         M{c}.Ac(:,7:12,6+i)    = [zeros(6);A*th2];  % Unique sess2 sequence pattterns
         M{c}.Ac(:,1:6,6+i)     = [zeros(6);A*th3];  % Correlation sess1-sess2
     end;    
end
end
function M = pcm_stability_genericCorr_noSess
% specific correlation models with no explicit modelling of session
corrS = 0:0.1:1;
for c=1:length(corrS)
    [th2,th3] = det_theta_corr(corrS(c));
    M{c}.type='feature';
    M{c}.numGparams = 2;
    M{c}.name = sprintf('SpecCorr_%1.1f',corrS(c));
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end; 
 
    M{c}.Ac(:,1:6,1)  = [A;zeros(6)];           % Unique sess1 sequence patterns
    M{c}.Ac(:,7:12,2) = [zeros(6);A*th2];       % Unique sess2 sequence pattterns
    M{c}.Ac(:,1:6,2)  = [zeros(6);A*th3];       % Correlation sess1-sess2
end
end
function M = pcm_stability_specificCorr_sess
% specific correlation models with additional modelling of session
% Model 1: Null model
M{1}.type       = 'feature';
M{1}.numGparams = 1;
M{1}.name       = 'null';
M{1}.Ac(:,1:12,1)  = zeros(12);

% Model 2: First vs. second session
M{2}.type       = 'feature';
M{2}.numGparams = 2;
M{2}.name       = 'Session';
M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];

% Other models: specific correlation
corrS = 0:0.1:1;
for c=1:length(corrS)
    [th2,th3] = det_theta_corr(corrS(c));
    M{c+2}.type='feature';
    M{c+2}.numGparams = 14;
    M{c+2}.name = sprintf('SpecCorr_%1.1f',corrS(c));
    M{c+2}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{c+2}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per session
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{c+2}.Ac(:,3:8,i+2)     = [A;zeros(6)];      % Unique sess1 sequence patterns
        M{c+2}.Ac(:,9:14,8+i)    = [zeros(6);A*th2];  % Unique sess2 sequence pattterns
        M{c+2}.Ac(:,3:8,8+i)     = [zeros(6);A*th3];  % Correlation sess1-sess2
    end;
end
end
function M = pcm_stability_genericCorr_sess
% specific correlation models with additional modelling of session
% Model 1: Null model
M{1}.type       = 'feature';
M{1}.numGparams = 1;
M{1}.name       = 'null';
M{1}.Ac(:,1:12,1)  = zeros(12);

% Model 2: First vs. second session
M{2}.type       = 'feature';
M{2}.numGparams = 2;
M{2}.name       = 'Session';
M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];

% Other models: specific correlation
corrS = 0:0.1:1;
for c=1:length(corrS)
    [th2,th3] = det_theta_corr(corrS(c));
    M{c+2}.type='feature';
    M{c+2}.numGparams = 4;
    M{c+2}.name = sprintf('SpecCorr_%1.1f',corrS(c));
    M{c+2}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{c+2}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end; 
    M{c+2}.Ac(:,3:8,3)  = [A;zeros(6)];         % Unique sess1 sequence patterns
    M{c+2}.Ac(:,9:14,4) = [zeros(6);A*th2];     % Unique sess2 sequence pattterns
    M{c+2}.Ac(:,3:8,3)  = [zeros(6);A*th3];     % Correlation sess1-sess2
end
end
function M = pcm_stability_meanPatternCorr
% specific correlation models with additional modelling of session
% for mean pattern
corrS = 0:0.1:1;
for c=1:length(corrS)
    [th2,th3] = det_theta_corr(corrS(c));
    M{c}.type='feature';
    M{c}.numGparams = 2;
    M{c}.name = sprintf('SpecCorr_%1.1f',corrS(c));
    M{c}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];        % Unique mean sess1 pattern
    M{c}.Ac(:,2,2) = [zeros(6,1);ones(6,1)*th2];    % Unique mean sess2 pattern
    M{c}.Ac(:,1,2) = [zeros(6,1);ones(6,1)*th3];    % Shared mean pattern
end
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
    Z=pcm_indicatorMatrix('identity',condVec{p});
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
    Z=pcm_indicatorMatrix('identity',condVec{p});
    %Z_seqType=pcm_indicatorMatrix('identity',condSeqTypeVec);
    % Subtract mean for each condition and run
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>6)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    Gcv_seqType(:,:,p)=pcm_estGCrossval(Data{p},partVec{p},condSeqTypeVec);
    C.r_crossval(p,1)=calcCorr(pcm_makePD(Gcv(:,:,p)));
    G_seqType=pcm_makePD(Gcv_seqType(:,:,p));
    C.r_crossval_seqType(p,1)=calcCorr_thetas(G_seqType(1,1),G_seqType(2,2),G_seqType(1,2));
end;

% --------------------------------------
% 3. Fit model 2  and infer correlations from the parameters
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec{p},condVec{p},'runEffect',runEffect);

% Get the correlations
switch M_type
    case 'specific'
        var1       = (theta{1}(3:8,:).^2)';
        var2       = (theta{1}(9:14,:).^2+theta{1}(15:20,:).^2)';
        cov12      = (theta{1}(3:8,:).*theta{1}(15:20,:))';
        C.r_model =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
    case 'generic'
        var1       = (theta{1}(3,:).^2)';
        var2       = (theta{1}(4,:).^2+theta{1}(5,:).^2)';
        cov12      = (theta{1}(3,:).*theta{1}(5,:))';
        C.r_model =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
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
% calculate correlation given thetas
    v1=th1^2;
    v2=th2^2+th3^2;
    cv=th1*th3;
    r=cv./sqrt(v1*v2);
end
function [th2,th3]=det_theta_corr(r)
% determine theta3 such that correlation is of specific value
if r ==1
    th2=0;
    th3=1;
else
    th2=1;
    th3=sqrt((r^2)/(1-r^2));
end
end

function R = splitHalfCorr(data,partVec,condVec,type)
% function R = splitHalfCorr(data,partVec,condVec,type)
% performs split-half correlations
switch(type)
    case 'withinSes'
        % calculate within session split-half correlation
        X = indicatorMatrix('identity_p',condVec);
        sPat=[]; % sequence pattern
        mPat=[]; % mean pattern
        split=mod(partVec,2);
        partN=numel(unique(partVec));
        condN=numel(unique(condVec));    
        % subtract the run mean
        data_rm = zeros(size(data));
        for i=1:partN
            data_rm(partVec==i,:) = bsxfun(@minus,data(partVec==i,:),sum(data(partVec==i,:),1)/condN);
        end;
        sPat(:,:,1)     = pinv(X(split==0,:))*data_rm(split==0,:);
        sPat(:,:,2)     = pinv(X(split==1,:))*data_rm(split==1,:);
        COR             = corr(sPat(:,:,1)',sPat(:,:,2)');
        R.corr_vox      = fisherinv(mean(fisherz(diag(COR))));
        mPat(:,:,1)     = pinv(X(split==0,:))*data(split==0,:);
        mPat(:,:,2)     = pinv(X(split==0,:))*data(split==1,:);
        COR_mean        = corr(mPat(:,:,1)',mPat(:,:,2)');
        R.corr_voxMean  = fisherinv(mean(fisherz(diag(COR_mean))));
        % create RDMs
        RDM(1,:)        = rsa.distanceLDC(data_rm(split==0,:),partVec(split==0,:),condVec(split==0,:));
        RDM(2,:)        = rsa.distanceLDC(data_rm(split==1,:),partVec(split==1,:),condVec(split==1,:));
        % normalize
        RDM_norm        = normalizeX(RDM);
        corrRDM         = RDM_norm*RDM_norm';
        R.corr_RDM      = corrRDM(2);
        R.corr_RDMAll   = 1;
    case 'acrossSes'
        Y     = [data{1};data{2}];
        %part  = [partVec{1};partVec{2}];
        part  = [partVec{1};partVec{2}+max(partVec{1})];
        cond  = [condVec{1};condVec{2}];
        partN = numel(unique(part));
        condN = numel(unique(cond));
        ses   = [ones(size(partVec{1}));ones(size(partVec{2}))*2];
        X     = indicatorMatrix('identity_p',cond);
        sPat=[];  % sequence pattern
        % subtract the run mean
        data_rm = zeros(size(Y));
        for i=1:partN
            data_rm(part==i,:) = bsxfun(@minus,Y(part==i,:),sum(Y(part==i,:),1)/condN);
        end;
        % Do split half correlations
        split = mod(part,2);
        sPat(:,:,1) = pinv(X(ses==1 & split==0,:)) * data_rm(ses==1 & split==0,:);
        sPat(:,:,2) = pinv(X(ses==1 & split==1,:)) * data_rm(ses==1 & split==1,:);
        sPat(:,:,3) = pinv(X(ses==2 & split==0,:)) * data_rm(ses==2 & split==0,:);
        sPat(:,:,4) = pinv(X(ses==2 & split==1,:)) * data_rm(ses==2 & split==1,:);
        COR=[];
        COR(:,:,1)  = corr(sPat(:,:,1)',sPat(:,:,3)');
        COR(:,:,2)  = corr(sPat(:,:,1)',sPat(:,:,4)');
        COR(:,:,3)  = corr(sPat(:,:,2)',sPat(:,:,3)');
        COR(:,:,4)  = corr(sPat(:,:,2)',sPat(:,:,4)');
        R.corr_vox  = fisherinv(mean(diag(mean(fisherz(COR),3))));
        % with mean present   
        mPat(:,:,1) = pinv(X(ses==1 & split==0,:)) * Y(ses==1 & split==0,:);
        mPat(:,:,2) = pinv(X(ses==1 & split==1,:)) * Y(ses==1 & split==1,:);
        mPat(:,:,3) = pinv(X(ses==2 & split==0,:)) * Y(ses==2 & split==0,:);
        mPat(:,:,4) = pinv(X(ses==2 & split==1,:)) * Y(ses==2 & split==1,:);
        COR_mean=[];
        COR_mean(:,:,1)  = corr(mPat(:,:,1)',mPat(:,:,3)');
        COR_mean(:,:,2)  = corr(mPat(:,:,1)',mPat(:,:,4)');
        COR_mean(:,:,3)  = corr(mPat(:,:,2)',mPat(:,:,3)');
        COR_mean(:,:,4)  = corr(mPat(:,:,2)',mPat(:,:,4)');
        R.corr_voxMean   = fisherinv(mean(diag(mean(fisherz(COR_mean),3))));
        
        % create RDMs
        RDM(1,:)    = rsa.distanceLDC(data_rm(ses==1 & split==0,:),part(ses==1 & split==0,:),cond(ses==1 & split==0,:));
        RDM(2,:)    = rsa.distanceLDC(data_rm(ses==1 & split==1,:),part(ses==1 & split==1,:),cond(ses==1 & split==1,:));
        RDM(3,:)    = rsa.distanceLDC(data_rm(ses==2 & split==0,:),part(ses==1 & split==0,:),cond(ses==2 & split==0,:));
        RDM(4,:)    = rsa.distanceLDC(data_rm(ses==2 & split==1,:),part(ses==2 & split==1,:),cond(ses==2 & split==1,:));
        % normalize
        RDM_norm    = normalizeX(RDM);
        corrRDM     = RDM_norm*RDM_norm';
        R.corr_RDM  = mean([corrRDM(1,3),corrRDM(1,4),corrRDM(2,3),corrRDM(2,4)]);
        % from all data (not splithalf)
        RDM_all(1,:)    = rsa.distanceLDC(Y(ses==1,:),part(ses==1,:),cond(ses==1,:));
        RDM_all(2,:)    = rsa.distanceLDC(Y(ses==2,:),part(ses==2,:),cond(ses==2,:));
        RDM_normAll     = normalizeX(RDM_all);
        corrRDM         = RDM_normAll*RDM_normAll';
        R.corr_RDMAll   = corrRDM(2);
        
end
end

