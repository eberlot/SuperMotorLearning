function varargout = sml_plasticity(what,varargin)
% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
caretDir        =[baseDir '/surfaceCaret'];     
behDir          =[baseDir '/behavioral_data'];
anaDir          =[behDir '/analyze'];
stabDir         =[baseDir '/stability'];
distPscDir      =[baseDir '/dist_psc_stats'];
regDir          =[baseDir '/RegionOfInterest/']; 
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session
codeDir         ='/Users/eberlot/Documents/MATLAB/projects/SuperMotorLearning';
wbDir           =[baseDir '/surfaceWB'];
% other info
regname_cortex  = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp'};
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen', 'Thalamus'};
hem             = {'lh','rh'};     
hemName         = {'LeftHem','RightHem'};   
atlasname       = 'fsaverage_sym';       
atlasA          = 'x';                                                            % freesurfer filename prefix

% ------------------------- Subject things --------------------------------
% The variables in this section must be updated for every new subject.
subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18',...
                's19','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'};  
sn          = [5:9,11:31]; % final selection of subjects
sessN       = 1:4; % all sessions unless specified otherwise 
style.file(fullfile(codeDir,'sml_style.m'));
style.use('default');
hemI      = {'L','R'}; 
hemname   = {'CortexLeft','CortexRight'};
seqTypeName = {'trained','control'};

switch what
    case 'PSC:surf_all'
        % run all surface analysis on percent signal change
       % sml_plasticity('PSC:create');
       % sml_plasticity('PSC:surface');
       % sml_plasticity('PSC:group_make');
       % sml_plasticity('PSC:group_cSPM');
       % sml_plasticity('PSC:group_smooth');
        sml_plasticity('PSC:create_rgb');
        sml_plasticity('PSC:trained_rgb');
    case 'PSC:create'
    % calculate psc for trained and untrained sequences - based on betas 
    % per subject; this is later used for surface projection
    sessN = 1:4;
    vararginoptions(varargin,{'sn','sessN'});
    name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','TrainSeq','UntrainSeq'};
    for ss=sessN
        for s=sn
           % cd(fullfile(glmSessDir{ss}, subj_name{s}));
           cd(fullfile(glmSessDir{ss},subj_name{s}));
            load SPM;
            X=(SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
            h=median(max(X));               % Height of response;
            P={};
            numB=length(SPM.xX.iB);         % Partitions - runs
            for p=SPM.xX.iB
                P{end+1}=sprintf('beta_%4.4d.nii',p);       % get the intercepts and use them to calculate the baseline (mean images)
            end;
            for con=1:length(name)    % 14 contrasts
                P{numB+1}=sprintf('con_%s.nii',name{con});
                outname=sprintf('psc_sess%d_%s.nii',ss,name{con}); % ,subj_name{s}       
                formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);        
                spm_imcalc_ui(P,outname,formula,...
                    {0,[],spm_type(16),[]});        % Calculate percent signal change
            end;
            fprintf('Sess %d: %s- %3.3f...psc done.\n',ss,subj_name{s},h);
        end;
    end
    fprintf('\nDone processing case PSC:create.\n');

    case 'PSC:cross-section'
        % create cross-section plots
        % of trained / control across sessions
        % use Z-score psc maps
        reg = 'SMA';
        width = 15;
        sessN = 1:4;
        seqName = {'trained','untrained'};
        vararginoptions(varargin,{'width','reg','sessN'});
        switch reg
            case 'SMA'
                from = [-74 125];
                to = [-10 150];
                reg = 5;
            case 'PMd1'
                from = [-27 86];
                to = [6 82];
                reg = 3;
            case 'PMd2'
                from = [-6 52];
                to = [3 101];
                reg = 3;
            case 'PMv'
                from = [-34 24];
                to = [0 31];
                reg = 4;
            case 'M1'
                from = [9 39];
                to = [19 105];
                reg = 2;
            case 'S1'
                from = [16 44];
                to = [38 102];
                reg = 1;
            case 'SPLa'
                from = [50 76];
                to = [114 121];
                reg = 7;
        end
        cd(fullfile(wbDir,'FS_LR_164'));
        T = gifti('ROI.164k.L.label.gii'); % load in the region file
        for st=1:2
            for ss=sessN
               % [Y{st}{ss},~,~] = surf_cross_section('fs_LR.164k.L.flat.surf.gii',sprintf('L.Z_psc_%s.sess-%d.func.gii',seqName{st},ss),'from',from,'to',to,'width',width);
                [Y{st}{ss},~,~] = surf_cross_section('fs_LR.164k.L.flat.surf.gii',sprintf('sL.psc_%s.sess-%d.func.gii',seqName{st},ss),'from',from,'to',to,'width',width);
                Yz{st}{ss} = zscore(Y{st}{ss});
            end
        end
        figure
        idx=1;
        for s2=1:3
            subplot(3,4,idx)
            traceplot(1:100,Yz{1}{1}','errorfcn','nanstd','linecolor','b','linewidth',2,'patchcolor','b','transp',0.2);
          %  traceplot(1:100,Yz{1}{1}','errorfcn','nanstd','linecolor',[0 0 0],'linewidth',2,'patchcolor',[100 100 100]./255,'transp',0.2);
            hold on;
            traceplot(1:100,Yz{1}{s2+1}','errorfcn','nanstd','linecolor','r','linewidth',2,'patchcolor','r','transp',0.2);
            drawline(0,'dir','horz'); title(sprintf('%s - sess-1-%d trained',regname_cortex{reg},s2+1));
            subplot(3,4,idx+1)
            traceplot(1:100,Yz{2}{1}','errorfcn','nanstd','linecolor',[0 0 0],'linewidth',2,'patchcolor',[100 100 100]./255,'transp',0.2);
            hold on;
            traceplot(1:100,Yz{2}{s2+1}','errorfcn','nanstd','linecolor','b','linewidth',2,'patchcolor','b','transp',0.2);
            drawline(0,'dir','horz'); title(sprintf('%s - sess-1-%d control',regname_cortex{reg},s2+1));
            subplot(3,4,idx+2)
            traceplot(1:100,Yz{1}{1}','errorfcn','nanstd','linecolor',[0 0 0],'linewidth',2,'patchcolor',[100 100 100]./255,'transp',0.2);
            hold on;
            traceplot(1:100,Yz{2}{s2+1}','errorfcn','nanstd','linecolor','b','linewidth',2,'patchcolor','b','transp',0.2);
            traceplot(1:100,Yz{1}{s2+1}','errorfcn','nanstd','linecolor','r','linewidth',2,'patchcolor','r','transp',0.2);
            drawline(0,'dir','horz'); title(sprintf('%s - sess-1-%d trained vs. control',regname_cortex{reg},s2+1));
            subplot(3,4,idx+3)
            traceplot(1:100,Yz{2}{s2+1}','errorfcn','nanstd','linecolor','b','linewidth',2,'patchcolor','b','transp',0.2);
            hold on;
            traceplot(1:100,Yz{1}{s2+1}','errorfcn','nanstd','linecolor','r','linewidth',2,'patchcolor','r','transp',0.2);
            drawline(0,'dir','horz'); title(sprintf('%s - sess-%d trained vs. control',regname_cortex{reg},s2+1));
            idx=idx+4;
        end
    case 'PSC:cross_section_allReg'
        % create cross-section plots
        % of trained across sessions 1 and 3
        % use Z-score psc maps
        reg = [1:3,7];
        width = 15;
        seqName = {'trained','untrained'};
        vararginoptions(varargin,{'width','reg','sessN'});
        cd(fullfile(wbDir,'FS_LR_164'));

        figure
        for r=1:4
            switch r
                case 3 % PMd
                    from = [-27 86];
                    to = [6 82];
                    %case 'PMd2'
                    %    from = [-6 52];
                    %    to = [3 101];
                    %    reg = 3;
                case 2 %M1
                    from = [9 39];
                    to = [19 105];
                case 1 %S1
                    from = [16 44];
                    to = [38 102];
                case 4 % SPLa
                    from = [50 76];
                    to = [114 121];
            end
                for ss=[1,3]
                    % [Y{st}{ss},~,~] = surf_cross_section('fs_LR.164k.L.flat.surf.gii',sprintf('L.Z_psc_%s.sess-%d.func.gii',seqName{st},ss),'from',from,'to',to,'width',width);
                    [Y{ss},~,~] = surf_cross_section('fs_LR.164k.L.flat.surf.gii',sprintf('sL.psc_%s.sess-%d.func.gii',seqName{1},ss),'from',from,'to',to,'width',width);
                    Yz{ss} = zscore(Y{ss});
                end
            subplot(2,2,r)
            traceplot(1:100,Yz{1}','errorfcn','nanstd','linecolor','b','linewidth',2,'patchcolor','b','transp',0.2);
            %  traceplot(1:100,Yz{1}{1}','errorfcn','nanstd','linecolor',[0 0 0],'linewidth',2,'patchcolor',[100 100 100]./255,'transp',0.2);
            hold on;
            traceplot(1:100,Yz{3}','errorfcn','nanstd','linecolor','r','linewidth',2,'patchcolor','r','transp',0.2);
            drawline(0,'dir','horz'); title(sprintf('%s - sess-1-3 trained',regname_cortex{reg(r)}));
            
        end
        
    case 'PSC:save_roi'
        % extract the PSC estimates per roi (from beta / psc ROI extraction)
        roi = 1:8;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','roi','seq','sessN','parcelType'});

        Stats = [];
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            for s=sn
                for r=roi
                    for h=1:2
                        R = getrow(T,T.regType==r & T.SN==s & T.regSide==h);
                        S.psc(1,:)  = nanmean(R.psc_train{1});
                        S.psc(2,:)  = nanmean(R.psc_untrain{1});
                        S.seqType   = [1;2];
                        S.sn        = repmat(R.SN,2,1);
                        S.roi       = repmat(R.region,2,1);
                        S.regType   = repmat(R.regType,2,1);
                        S.regSide   = repmat(R.regSide,2,1);
                        S.sessN     = repmat(ss,2,1);
                        % save into a new structure
                        Stats=addstruct(Stats,S);
                    end
                end
            end
            fprintf('Done session %d.\n',ss);
        end
        % save structure
        save(fullfile(distPscDir,sprintf('psc_%s_ROI',parcelType)),'-struct','Stats');
    case 'PSC:plot_roi'
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';      
        vararginoptions(varargin,{'roi','sessN','hemi','parcelType'});     
        T=load(fullfile(distPscDir,sprintf('psc_%s_ROI.mat',parcelType)));   
        switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        style.use('Seq');
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([T.sessN>3 T.sessN],T.psc,'split',T.seqType,'subset',T.regType==r&T.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north');
            elseif numel(sessN)==3
                plt.line(T.sessN,T.psc,'split',T.seqType,'subset',T.regType==r&T.sessN<4&T.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north');
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regLab{r}));
            if r==1
                ylabel('Percent signal change');
                xlabel('Session');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'PSC:stats_roi'
        roi=1:8;
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        T=load(fullfile(distPscDir,sprintf('psc_%s_ROI.mat',parcelType)));
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n PSC in %s across sess:1-4\n',regLab{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi);
            fprintf('\n PSC in %s across sess:1-3\n',regLab{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi & T.sessN<4);
            fprintf('\n PSC in %s across sess:3-4\n',regLab{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi & T.sessN>2);
        end
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regLab{r});
                ttestDirect(T.psc,[T.seqType T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & T.sessN==ss);
            end
            for str=1:3 % session transition
                fprintf('\n post-hoc t-test on the effect of session transition %d-%d in %s \n',str,str+1,regLab{r});
                ttestDirect(T.psc,[T.sessN T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & ismember(T.sessN,[str,str+1]),'split',T.seqType);
            end
        end
    
    case 'SEARCH:group_make'                 % DELETE ALL OF THESE                               % STEP 4.5   :  Make group metric files by condensing subjec contrast metric files
        % Calculate group metric files from the searchlight results. 
        % Takes the 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross') across
        % subjects and makes a group level metric file that contains each
        % subject's data for that contrast type.
        name = 'dist';
        inputcol   = [1 2 3 4];
        replaceNaN = [1 1 1 1];
        OUTname    = {'dist_all','dist_trained','dist_untrained','dist_cross'};
        vararginoptions(varargin,{'sessN','sn','name','OUTname','inputcol','replaceNaN'});
        
        for ss=sessN
            % Loop over hemispheres.
            for h = 1:2
                % initialise some variables
                outfilenames    = cell(length(OUTname),1);
                infilenames     = cell(length(OUTname),length(sn));
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = [caretDir filesep atlasname filesep hemName{h}];
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'OUTname' and make a group metric file
                for j = 1:length(OUTname);
                    % Loop over subjects...
                    for i = 1:length(sn);
                        % ...and define the names of their metric files
                        infilenames{j}{i} = fullfile(caretDir,[atlasA subj_name{sn(i)}], hemName{h}, sprintf('%s_sess%d_%s.metric',subj_name{sn(i)},ss,name));
                    end;
                    outfilenames{j} = [surfaceGroupDir filesep hem{h} '.' OUTname{j} '_sess' num2str(ss) '.metric'];
                    % Finally, make the group metric file for this metric type/contrast
                    caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
                    fprintf('Session %d: %s - image: %i/%i group searchlight file created.\n',ss,hemName{h},j,length(OUTname));
                end; % filenames
            end; % hemisphere
        end; % session
        fprintf('\nDone processing case SEARCH:group_make.\n');
    case 'SEARCH:group_cSPM'                                                % STEP 4.6   :  Generate a statistical surface map (onesample_t test) from smoothed group metric files. Also avgs. distances across subjs.
        % Calculate group stats files from the group metric files. 
        % Takes 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross')  and 
        % calculates group level stats (one sample t test that the mean 
        % effect is bigger than zero).  
        % Although we calculate the t-score, corresponding p-value for that
        % t-score, and finally the z-score 
        % subject's data for that contrast type.    
        SPMname={'dist_all','dist_trained','dist_untrained','dist_cross'};
        sqrtTransform=[1,1,1,1]; % Should you take ssqrt before submitting? yes for distances
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name'});
        s=1:length(sn);
        for ss=sessN
            SummaryName = sprintf('.summary_dist_sess%d.metric',ss);
            hemi = [1 2];
            for h=hemi
                sfilenames = cell(length(SPMname),1); % initialise
                column_name = cell(1,length(SPMname));
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
                %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
                for i=1:length(SPMname);
                    %sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '.metric']; % smoothed
                    sfilenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '_sess' num2str(ss) '.metric']; % no smoothing
                    %----loop over the metric files and calculate the cSPM of each with the non-smoothed metrics
                    Data=caret_load(sfilenames{i});
                    if sqrtTransform(i)
                        Data.data=ssqrt(Data.data);
                    end;
                    cSPM=caret_getcSPM('onesample_t','data',Data.data(:,s),'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                    caret_savecSPM([surfaceGroupDir filesep hem{h} '.' SPMname{i} '_stats.metric'],cSPM);
                    save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');
                    % here for surface plots
                    if i==1 % initialise data
                        data = zeros(size(Data.data,1),length(SPMname)*2);
                    end
                    data(:,i)=cSPM.con(1).con; % mean
                    data(:,i+length(SPMname))=cSPM.con(1).Z; % T
                    column_name{i}=['mean_' SPMname{i} '_sess' num2str(ss)];
                    column_name{i+length(SPMname)}=['T_' SPMname{i} '_sess' num2str(ss)];
                end;
                C = caret_struct('metric','data',data,'column_name',column_name);
                caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);   
            end;
            fprintf('Done group stats: sess-%d\n',ss);
        end
        fprintf('\nDone processing case SEARCH:group_cSPM.\n');
    case 'SEARCH:group_smooth'
        vararginoptions(varargin,{'sessN'});
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
            cd(surfaceGroupDir)
            %----define name of coord and topology
            coordfile=[hem{h} '.WHITE.coord'];
            topofile=[hem{h} '.CLOSED.topo'];
            for ss=sessN
                %----get the full directory name of the metric files and the smoothed metric files that we create below
                for i=1:length(ss);
                    filename=[surfaceGroupDir filesep hem{h} '.summary_dist_sess' num2str(ss) '.metric']; % unsmoothed
                    sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
                end;
                fprintf('Sess %d: %s - smoothing done.\n',ss,hemName{h});
            end; % session
        end; % hemisphere
        fprintf('\nDone processing case SEARCH:group_smooth.\n');
    case 'SEARCH:create_rgb'
        smooth=0;
        lowS=1.0;
        maxS=8;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            for ss=sessN
                if smooth == 1
                    Cd=caret_load(['s' hem{h} sprintf('.summary_dist_sess%d.metric',ss)]);
                else
                    Cd=caret_load([hem{h} sprintf('.summary_dist_sess%d.metric',ss)]);
                end
                %------Dist
                RGBdata(:,1)=Cd.data(:,6); % Red: trained distances
                RGBdata(:,3)=Cd.data(:,7); % Blue: untrained distances
                RGBdata(:,4)=Cd.data(:,6)-Cd.data(:,7); % difference: trained - untrained
                RGBdata(:,6)=Cd.data(:,7)-Cd.data(:,6); % difference: untrained - trained
                RGBdata(RGBdata(:,1)<lowS,1)=0;
                RGBdata(RGBdata(:,3)<lowS,3)=0;
                RGBdata(RGBdata(:,4)<lowS,4)=0;
                RGBdata(RGBdata(:,6)<lowS,6)=0;
                % Scaling for Distances
                sc1=[lowS maxS;lowS maxS;lowS maxS];               
                sc2=[0.5 4;0.5 4;0.5 4];        
                
                name={sprintf('dist_sess%d',ss),sprintf('dist_difference_sess%d',ss)};                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1,sc2},'column_name',name);        
                if smooth == 1
                    caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.dist_sess%d.RGB_paint',ss)],C);
                else
                    caret_save([surfaceGroupDir filesep hem{h} sprintf('.dist_sess%d.RGB_paint',ss)],C);
                end
                fprintf('Sess %d: %s - rgb map done.\n',ss,hemName{h});
            end;
        end
        fprintf('\nDone processing case SEARCH: create_rgb.\n');
        
    case 'CORR:withinSess'
        % within session correlation of trained - untrained
        beta_choice = 'mw';
        reg         = 1:16;
        parcelType  = 'Brodmann'; 
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN'})
        
        AllReg=[];
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r = reg
                for p=1:length(sn)
                    % load betas, info
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    t = getrow(B,B.SN==sn(p) & B.region==r);
                    switch (beta_choice)
                        case 'uw'
                            beta = t.betaUW{:};
                        case 'mw'
                            beta = t.betaW{:};
                        case 'raw'
                            beta = t.betaRAW{:};
                    end
                    condVec  = D.seqType; % trained vs. untrained
                    Data     = beta(1:size(D.seqNumb,1),:);
                    Z        = pcm_indicatorMatrix('identity',condVec);
                    U        = pinv(Z)*Data; % mean activities
                  %  U        = bsxfun(@minus,U,mean(U,1)); % Subtract mean per condition - trained
                    % here calculate correlation
                    corrU    = corr(U');
                    %corrTr   = rsa_vectorizeRDM(corrU(1:6,1:6));
                    %corrC    = rsa_vectorizeRDM(corrU(7:end,7:end));
                    %corrAcr  = rsa_vectorizeRDM(corrU(1:6,7:end));
                    % fit PCM, save results
                    %T.corrTrained   = mean(corrTr);
                    %T.corrControl   = mean(corrC);
                    %T.corrAcross    = mean(corrAcr);
                    T.corrAcross    = corrU(2);
                    T.roi           = r;
                    T.regType       = t.regType;
                    T.regSide       = t.regSide;
                    T.sessN         = ss;
                    T.sn            = p;
                    AllReg          = addstruct(AllReg,T);
                end; % subject
                fprintf('Done: session:%d - reg:%d/%d\n\n',ss,r,numel(reg));
            end; % region
        end; % session
        save(fullfile(stabDir,sprintf('CORR_cross_withinSess_%s.mat',parcelType)),'-struct','AllReg');    
    case 'CORR:withinSess_permute'
        % here permute the order of assignment into trained / untrained
        beta_choice = 'mw';
        reg         = 1:16;
        parcelType  = 'Brodmann'; 
        nPerm       = 1000; % number of permutations
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','nPerm'})
        AllReg=[];
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r = reg
                for p=1:length(sn)
                    % load betas, info
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    t = getrow(B,B.SN==sn(p) & B.region==r);
                    switch (beta_choice)
                        case 'uw'
                            beta = t.betaUW{:};
                        case 'mw'
                            beta = t.betaW{:};
                        case 'raw'
                            beta = t.betaRAW{:};
                    end
                    Data     = beta(1:size(D.seqNumb,1),:);
                    corrU = zeros(nPerm,1);        
                    pr = 1;
                    while pr<nPerm+1
                        perm = randperm(12);
                        if ~(ismember(perm,1:12,'rows') || ismember(perm,[7:12 1:6],'rows')) % ensure the group assignment is not correct
                            sqType   = zeros(12,1);
                            sqType(perm<7) = 1;
                            sqType(perm>6) = 2;
                            condVec  = repmat(sqType,8,1); % new condition vector
                            Z        = pcm_indicatorMatrix('identity',condVec);
                            U        = pinv(Z)*Data; % mean activities
                            corrTmp  = corr(U');
                            corrU(pr) = corrTmp(2);
                            pr=pr+1;
                        end
                    end
                    % here save variables
                    T.corrPerm      = corrU';
                    T.corrMean      = mean(corrU);
                    T.corrStd       = std(corrU);
                    T.roi           = r;
                    T.regType       = t.regType;
                    T.regSide       = t.regSide;
                    T.sessN         = ss;
                    T.sn            = p;
                    AllReg          = addstruct(AllReg,T);
                end; % subject
                fprintf('Done: session:%d - reg:%d/%d\n\n',ss,r,numel(reg));
            end; % region
        end; % session
        save(fullfile(stabDir,sprintf('CORR_cross_permute_withinSess_%s.mat',parcelType)),'-struct','AllReg');   
    case 'CORR:plot_withinSess'
        reg = [1:5,7,8];
        hemi = 1;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','hemi','parcelType'});
        
        T = load(fullfile(stabDir,sprintf('CORR_cross_withinSess_%s.mat',parcelType)));
        P = load(fullfile(stabDir,sprintf('CORR_cross_permute_withinSess_%s.mat',parcelType)));
        figure
        style.use('Region');
        subplot(121)
        plt.line([T.sessN>3 T.sessN],T.corrAcross,'subset',T.regSide==hemi & ismember(T.regType,reg),'split',T.regType,'leg',regname_cortex(reg));
        ylabel('correlation'); xlabel('session');
        title('mean trained - untrained pattern correlation');
        subplot(122)
        plt.line([P.sessN>3 P.sessN],P.corrMean,'subset',P.regSide==hemi & ismember(P.regType,reg),'split',P.regType,'leg',regname_cortex(reg));
        ylabel(''); title('permuted assignment');
        plt.match('y');
        
        for r=reg
            figure
            style.use('Region_shade');
            plt.line(T.sessN,T.corrAcross,'subset',T.regSide==hemi & ismember(T.regType,r),'split',T.regType,'leg',regname_cortex(r));
            hold on;
            style.use('Region_ceiling');
            plt.line(P.sessN,P.corrMean,'subset',P.regSide==hemi & ismember(P.regType,r),'split',P.regType,'leg',regname_cortex(r));
            title(regname_cortex{r});xlabel('session'); ylabel('seqType correlation');
        end
    case 'CORR:stats'
        reg = [1:5,7,8];
        hemi = 1;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','hemi','parcelType'});
        
        T = load(fullfile(stabDir,sprintf('CORR_cross_withinSess_%s.mat',parcelType)));
        P = load(fullfile(stabDir,sprintf('CORR_cross_permute_withinSess_%s.mat',parcelType)));
        
        for r=reg
            for ss=sessN
                t=getrow(T,T.regType==r&T.sessN==ss&T.regSide==hemi);
                p=getrow(P,P.regType==r&P.sessN==ss&P.regSide==hemi);
                j=sort(p.corrPerm,2);
                pr=bsxfun(@minus,j,t.corrAcross);
                pos = pr>0;
                [~, ind] = max( pos~=0,[],2);
                prob = median(ind./1000); % average probability
                fprintf('Session %d - %s: %1.4f\n',ss,regname_cortex{r},prob);
            end
        end
    
    case 'CORR:acrossSess'
        % do comparisons of mean pattern (trained / untrained) to sess1
        reg         = 1:16;
        parcelType  = 'Brodmann'; 
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN'})
        
        AllReg=[];
        
        for ss=1:4
            B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
        end
        
        for r = reg
            for p=1:length(sn)
                % load betas, info
                glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                condVec  = D.seqType; % trained vs. untrained
                Z        = pcm_indicatorMatrix('identity',condVec);
                for ss=1:4
                    t{ss} = getrow(B{ss},B{ss}.SN==sn(p) & B{ss}.region==r);
                    beta{ss} = t{ss}.betaW{:};
                    Data{ss} = beta{ss}(1:size(D.seqNumb,1),:);
                    U{ss}    = pinv(Z)*Data{ss};
                end
                % here calculate correlation
                corrT = corr([U{1}(1,:)' U{2}(1,:)' U{3}(1,:)' U{4}(1,:)']);
                corrU = corr([U{1}(2,:)' U{2}(2,:)' U{3}(2,:)' U{4}(2,:)']);
                T.corrT     = rsa_vectorizeRDM(corrT);
                T.corrU     = rsa_vectorizeRDM(corrU);
                T.roi       = r;
                T.regType   = t{1}.regType;
                T.regSide   = t{1}.regSide;
                T.sn        = p;
                AllReg      = addstruct(AllReg,T);
            end; % subject
            fprintf('Done: reg:%d/%d\n\n',r,numel(reg));
        end; % region    
        save(fullfile(stabDir,sprintf('CORR_acrossSess_%s.mat',parcelType)),'-struct','AllReg');
    case 'CORR:plot_acrossSess'
        reg = 1:8;
        hemi=1;
        numTr=3; % 3 or 2
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','numTr'});
        T = load(fullfile(stabDir,sprintf('CORR_acrossSess_%s.mat',parcelType)));
        style.use('Seq');
        switch numTr
            case 2
                subS = [1,4];
            case 3
                subS = [1,4,6];
        end
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            cT = t.corrT(:,1:numTr); cT = cT(:);
            cU = t.corrU(:,1:numTr); cU = cU(:);
            c.all = [cT(:); cU(:)];
            cTS = t.corrT(:,subS); cTS = cTS(:);
            cUS = t.corrU(:,subS); cUS = cUS(:);
            c.allS = [cTS(:); cUS(:)];
            c.seqType = [ones(size(cT(:)));ones(size(cT(:)))*2];
            if numTr==3
                c.sessN   = repmat([ones(size(t.corrT,1),1);ones(size(t.corrT,1),1)*2;ones(size(t.corrT,1),1)*3],2,1);
            else
                c.sessN   = repmat([ones(size(t.corrT,1),1);ones(size(t.corrT,1),1)*2],2,1);
            end
            c.SN      = repmat([1:26]',numTr*2,1);
            c=normData(c,'all');
            c=normData(c,'allS');
            figure
            subplot(121)
            plt.line(c.sessN,c.normall,'split',c.seqType); title(sprintf('%s - relative to 1',regname_cortex{r}));
            xlabel('session transition'); ylabel('correlation'); set(gca,'XTickLabel',{'1-2','1-3','1-4'});
            subplot(122)
            plt.line(c.sessN,c.normallS,'split',c.seqType); title(sprintf('%s - subsequent',regname_cortex{r}));
            plt.match('y');
            xlabel('session transition'); ylabel(''); set(gca,'XTickLabel',{'1-2','2-3','3-4'});
            ttestDirect(c.all,[c.seqType c.SN],2,'paired','split',c.sessN)
            anovaMixed(c.all,c.SN,'within',[c.seqType c.sessN],{'seqType','session'})
           % anovaMixed(c.allS,c.SN,'within',[c.seqType c.sessN],{'seqType','session'})
            keyboard;
        end
        
    case 'PCM:withinSess'
        runEffect   = 'fixed';
        beta_choice = 'mw';
        algorithm   = 'NR'; % minimize or NR
        reg         = 1:16;
        parcelType  = 'Brodmann'; 
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','runEffect'})
        
        AllReg=[];
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r = reg
                for st = 1:2 % trained / untrained
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
                    for p=1:length(sn)
                        % load betas, info
                        glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        t = getrow(B,B.SN==sn(p) & B.region==r);
                        switch (beta_choice)
                            case 'uw'
                                beta = t.betaUW{:};
                            case 'mw'
                                beta = t.betaW{:};
                            case 'raw'
                                beta = t.betaRAW{:};
                        end  
                        condVec{p}  = D.seqNumb(D.seqType==st); % conditions
                        partVec{p}  = D.run(D.seqType==st);
                        Data{p}     = beta(D.seqType==st,:);  
                        if st==2
                            condVec{p}=condVec{p}-6;    % make conditions always from 1
                        end
                    end; % subject
                    % construct models
                    M = pcm_withinSess;
                    % fit PCM, save results
                    T           = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                    T.roi       = ones(size(T.SN))*r;
                    T.regType   = ones(size(T.SN))*t.regType;
                    T.regSide   = ones(size(T.SN))*t.regSide;
                    T.seqType   = ones(size(T.SN))*st;
                    T.sessN     = ones(size(T.SN))*ss;
                    AllReg      = addstruct(AllReg,T);
                    fprintf('Done: session:%d - reg:%d/%d\n\n',ss,r,numel(reg));
                end; % seqType
            end; % region    
        end; % session
        % save output
        save(fullfile(stabDir,sprintf('PCM_withinSess_%s.mat',parcelType)),'-struct','AllReg');    
    case 'PCM:stability_toyModels'
        % here construct PCM toy models (zero correlation, perfect
        % correlation, flexible correlation)
        % flexible correlation model estimates the correlation
        runEffect   = 'random';
        beta_choice = 'mw'; % multivariately prewhitened betas
        algorithm   = 'NR'; % minimize or NR
        reg         = 1:8;
        parcelType  = 'Brodmann';
        modelType   = 'wSess';
        sessTr      = 1:3; % transissions
        AllReg=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessTr','algorithm','parcelType','modelType'});
        for str = sessTr % session transition
            sessN = [str str+1];
            B = cell(length(sessN),1);
            for ss=1:numel(sessN)
                B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
            end
            for r = reg
                for st=1:2 % per seqType: trained / untrained
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
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
                                condVec{p} = D.seqNumb(indx==1,:);  % conditions
                                partVec{p} = D.run(indx==1,:);
                                Data{p} = beta(indx==1,:);          % Data is N x P (cond x voxels) - no intercept
                            else
                                condVec{p} = [condVec{p}; D.seqNumb(indx==1,:) + 12];   % treat 2nd session as additional conditions
                                partVec{p} = [partVec{p}; D.run(indx==1,:) + 8];        % runs/partitions of 2nd session as additional runs
                                Data{p} = [Data{p}; beta(indx==1,:)];                   % Data is N x P (cond x voxels) - no intercept
                            end;
                        end; % session
                    end; % subj
                    % construct models
                    switch modelType
                        case 'wSess'
                            M = pcm_toyModel_wSess;
                            C = pcm_correlation(Data,partVec,condVec,M{4},runEffect,modelType);
                        case 'noSess'
                            M = pcm_toyModel_noSess; 
                            C = pcm_correlation(Data,partVec,condVec,M{3},runEffect,modelType);
                    end
                    
                    T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*r;
                    T.regType = ones(size(T.SN))*t.regType;
                    T.regSide = ones(size(T.SN))*t.regSide;
                    T.seqType = ones(size(T.SN))*st;
                    T.sessTr  = ones(size(T.SN))*str;
                    T = rmfield(T,'reg');
                    AllReg=addstruct(AllReg,T);
                    AllReg=addstruct(AllReg,C);       
                end
                fprintf('Done sess: %d-%d model-%s reg: %d/%d\n\n',sessN(1),sessN(2),modelType,r,max(reg))
            end
            fprintf('Done sess: %d-%d model-%s\n\n\n\n',sessN(1),sessN(2),modelType);
        end
        % save output
        save(fullfile(stabDir,sprintf('PCM_toyModels_%s_%s.mat',parcelType,modelType)),'-struct','AllReg');
    case 'PCM:corrModels'
        % here assessing models of different correlation values
        runEffect   = 'random';
        beta_choice = 'mw';
        algorithm   = 'NR'; % minimize or NR
        reg         = 1:16;
       % reg          = [3,7];
        parcelType  = 'Brodmann';
        modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessType    = 'transitions'; % options: relativeto1 or transitions
        checkCorr   = 0;        % check if correlation exact and plot predicted Gs
        sessTr      = 1:3;      % transitions
        corrLim     = [0 1];    % bounds for lower / upper correlation of models
        nModel      = 30;       % number of correlation models (determines how fine grained the corr estimates are)
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessTr','algorithm','parcelType','modelType','nModel','sessType'});
        AllReg=[];
        for str = sessTr % session transition
            switch sessType
                case 'transitions' 
                    sessN = [str str+1];
                case 'relativeto1'
                    sessN = [1 str+1];
            end
            B = cell(length(sessN),1);
            for ss=1:numel(sessN)
                B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
            end
            % construct omdels
            switch modelType
                case 'noSess'
                    M = pcm_corrModel_noSess(corrLim,nModel);
                case 'wSess'
                    M = pcm_corrModel_wSess(corrLim,nModel);
                case 'meanPattern'
                    M = pcm_corrModel_meanPattern(corrLim,nModel);
            end  
            for r = reg
                for st=1:2 % per seqType
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
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
                                condVec{p}  = D.seqNumb(indx==1,:); % conditions
                                partVec{p}  = D.run(indx==1,:);
                                Data{p}     = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                            else
                                condVec{p}  = [condVec{p}; D.seqNumb(indx==1,:)+12]; % treat 2nd session as additional conditions
                                partVec{p}  = [partVec{p}; D.run(indx==1,:)+8];  % runs/partitions of 2nd session as additional runs
                                Data{p}     = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                            end;
                        end; % session
                    end; % subj
                    % here just test if correlations correct
                    if checkCorr % later on remove this
                        [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
                        figure
                        corrPred=zeros(11,1);
                        for c=1:length(corrSpec)
                            subplot(3,4,c)
                            imagesc(G_pred{c});
                            corrPred(c)=calcCorr(G_pred{c});
                        end
                    end
                    % group fit
                   % T = pcm_fitModels_conserv(Data,M,partVec,condVec,runEffect,algorithm);
                    T = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitAlgorithm',algorithm,'fitScale',0); % only group fit
                    T.bayesEst = bsxfun(@minus,T.likelihood,T.likelihood(:,1));
                    %T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                    %T = pcm_fitModels_corr(Data,M,partVec,condVec,runEffect,algorithm);
                    % individual fit
                    I = pcm_fitModelIndivid(Data,M,partVec{p},condVec{p},'runEffect',runEffect);
                    T.individ_likelihood    = I.likelihood;
                    T.bayesEst_individ      = bsxfun(@minus,T.individ_likelihood,T.individ_likelihood(:,1));
                    T.individ_noise         = I.noise;
                    T.roi                   = ones(size(T.SN))*r;
                    T.regType               = ones(size(T.SN))*t.regType;
                    T.regSide               = ones(size(T.SN))*t.regSide;
                    T.seqType               = ones(size(T.SN))*st;
                    T.sessTr                = ones(size(T.SN))*str;
                    T                       = rmfield(T,'reg');
                    AllReg = addstruct(AllReg,T);
                    fprintf('Done modelType: %s seqType %d sess: %d-%d reg: %d/%d\n\n',modelType,st,sessN(1),sessN(2),r,length(reg));
                end; % seqType
                fprintf('Done reg %d/%d:\tmodelType: %s \tsess: %d-%d\n\n',find(r==reg),numel(reg),modelType,sessN(1),sessN(2)); 
            end; % region
            fprintf('Done all:\tmodelType: %s \tsess: %d-%d\n\n\n\n',modelType,sessN(1),sessN(2)); 
        end; % session transition
        % save output
        save(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)),'-struct','AllReg');
    case 'PCM:corrModels_mean_G'
        % here assessing models of different correlation values
        % based on G matrix
        reg         = 1:16;
       % reg          = [3,7];
        parcelType  = 'Brodmann';
        sessType    = 'transitions'; % options: relativeto1 or transitions
        sessTr      = 1:3;      % transitions
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessTr','parcelType','sessType'});
        TT=[];
        for str = sessTr % session transition
            switch sessType
                case 'transitions' 
                    sessN = [str str+1];
                case 'relativeto1'
                    sessN = [1 str+1];
            end
            B = cell(length(sessN),1);
            for ss=1:numel(sessN)
                B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
            end
            for r = reg
                for st=1:2 % per seqType
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
                    for p=1:length(sn)
                        for ss = 1:numel(sessN)
                            glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                            t = getrow(B{ss},B{ss}.SN==sn(p) & B{ss}.region==r);
                            beta = t.betaUW{:};
                            indx = D.seqType==st;
                            if ss == 1
                                condVec{p}  = ones(size(find(D.seqType==st))); % conditions
                                partVec{p}  = D.run(D.seqType==st);
                                Data{p}     = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                            else
                                condVec{p}  = [condVec{p}; ones(size(find(D.seqType==st)))*2]; % treat 2nd session as additional conditions
                                partVec{p}  = [partVec{p}; D.run(indx==1,:)+8];  % runs/partitions of 2nd session as additional runs
                                Data{p}     = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                            end;
                        end; % session
                        G = pcm_estGCrossval(Data{p},partVec{p},condVec{p});
                        corrG = G(1,2)/sqrt(G(1,1)*G(2,2));
                        T.corr = corrG;
                        T.seqType = st;
                        T.sessTr = str;
                        T.sn = sn(p);
                        T.region = r;
                        TT = addstruct(TT,T);
                    end; % subj
                end; % seqType
                fprintf('Done reg %d/%d:\tsess: %d-%d\n\n',find(r==reg),numel(reg),sessN(1),sessN(2)); 
            end; % region
            fprintf('Done all:\t\tsess: %d-%d\n\n\n\n',sessN(1),sessN(2)); 
        end; % session transition
        % save output
        save(fullfile(stabDir,sprintf('Corr_acrossSess_%s_%s.mat',parcelType,sessType)),'-struct','TT');
    case 'PCM:corrModels_mean' 
        reg         = 1:16;
       % reg          = [3,7];
        parcelType  = 'Brodmann';
        sessType    = 'transitions'; % options: relativeto1 or transitions
        sessTr      = 1:3;      % transitions
        corrLim     = [0 1]; % bounds for correlation
        nModel      = 30; % number of models
        
        vararginoptions(varargin,{'corrLim','nModel','sn','reg','sessTr','parcelType','sessType'});
        TT=[];
        % build models
        corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
        for c=1:length(corrS)
                 M{c} = pcm_buildCorrModel('type','nonlinear','withinCov','iid','numCond',2,'numItems',6,'r',corrS(c),'condEffect',0);
       %     M{c} = pcm_buildCorrModel('type','nonlinear','withinCov','iid','numCond',2,'numItems',1,'r',corrS(c),'condEffect',0);
        end
        for str = sessTr % session transition
            switch sessType
                case 'transitions' 
                    sessN = [str str+1];
                case 'relativeto1'
                    sessN = [1 str+1];
            end
            B = cell(length(sessN),1);
            for ss=1:numel(sessN)
                B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
            end
            for r = reg
                for st=1:2 % per seqType
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
                    for p=1:length(sn)
                        for ss = 1:numel(sessN)
                            glmDirSubj=fullfile(glmSessDir{sessN(ss)}, subj_name{sn(p)});
                            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                            t = getrow(B{ss},B{ss}.SN==sn(p) & B{ss}.region==r);
                            beta = t.betaUW{:};
                            if ss == 1
                                Data{p} = beta(D.seqType==st,:);
                                partVec{p}  = D.run(D.seqType==st,:); 
                                %condVec{p}  = ones(sum(D.seqType==st),1); % conditions
                                condVec{p} = D.seqNumb(D.seqType==st);
                                %for rr = 1:8 % number of runs
                                %    Data{p}(rr,:)=mean(beta(D.seqType==st & D.run==rr,:),1);
                                % end
                                %  condVec{p}  = ones(8,1); % conditions
                                %  partVec{p}  = (1:8)';
                            else
                                %condVec{p}  = [condVec{p}; ones(sum(D.seqType==st),1)+1]; % treat 2nd session as additional conditions
                                condVec{p} = [condVec{p}; D.seqNumb(D.seqType==st)+6];
                                partVec{p}  = [partVec{p}; D.run(D.seqType==st,:)+8];  % runs/partitions of 2nd session as additional runs
                                Data{p}     = [Data{p}; beta(D.seqType==st,:)];  % Data is N x P (cond x voxels) - no intercept
                                % condVec{p}  = [condVec{p}; ones(8,1)*2]; % treat 2nd session as additional conditions
                                %  partVec{p}  = [partVec{p}; (1:8)'+8];  % runs/partitions of 2nd session as additional runs
                                %  for rr = 1:8 % number of runs
                                %     Data{p}=[Data{p};mean(beta(D.seqType==st & D.run==rr,:),1)];
                                %  end
                            end;
                        end; % session
                    end; % subj
                    T = sml_plasticity('PCM:corrModel_run','Data',Data,'model',M,'partV',partVec,'condV',condVec,'algorithm','NR','runEffect','none');
                    % other variables of interest
                    T.roi                   = ones(size(T.SN))*r;
                    T.regType               = ones(size(T.SN))*t.regType;
                    T.regSide               = ones(size(T.SN))*t.regSide;
                    T.seqType               = ones(size(T.SN))*st;
                    T.sessTr                = ones(size(T.SN))*str;
                    T = rmfield(T,'reg');
                    TT = addstruct(TT,T);
                end; % seqType
                fprintf('Done reg %d/%d:\tsess: %d-%d\n\n',find(r==reg),numel(reg),sessN(1),sessN(2)); 
            end; % region
            fprintf('Done all:\t\tsess: %d-%d\n\n\n\n',sessN(1),sessN(2)); 
        end; % session transition
        save(fullfile(stabDir,sprintf('PCM_corrModels_%s_nonlinear_meanPattern_v2_%s.mat',parcelType,sessType)),'-struct','TT');
    case 'PCM:corrModels_nonlinear'
        % here assessing models of different correlation values
        % across sessions
        runEffect   = 'random';
        beta_choice = 'mw';
        algorithm   = 'NR'; % minimize or NR
        reg         = 1:16;
        hemi        = 1;
        parcelType  = 'Brodmann';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessType    = 'transitions'; % options: relativeto1 or transitions
        sessTr      = 1:3;      % transitions
        corrLim     = [0 1];    % bounds for lower / upper correlation of models - change this to [-1 1]!!!
        nModel      = 30;       % number of correlation models (determines how fine grained the corr estimates are)
        withinCov   = 'individual'; % covariance type: individual or iid (or meanPattern)
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessTr','algorithm','parcelType','withinCov','nModel','sessType'});
        AllSess=[];
        if any(strcmp(parcelType,{'tesselsWB_162','tesselsWB_362'}))
            reg = sml_plasticity('TESSEL:select','sessN',1:3,'betaChoice','multi','hemi',hemi,'tesselType',parcelType)';
        end
        for str = sessTr % session transition
            AllReg = [];
            % construct models
            corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
            for c=1:length(corrS)
                if strcmp(withinCov,'meanPattern')
                    M{c} = pcm_buildCorrModel('type','nonlinear','withinCov','iid','numCond',2,'numItems',6,'r',corrS(c),'condEffect',0);
                    %M{c} = pcm_buildCorrModel('type','nonlinear','withinCov','mean','numCond',2,'numItems',6,'r',corrS(c),'condEffect',0);
                else
                    M{c} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',2,'numItems',6,'r',corrS(c));
                end
            end
            switch sessType
                case 'transitions'
                    sessN = [str str+1];
                case 'relativeto1'
                    sessN = [1 str+1];
            end
            B = cell(length(sessN),1);
            for ss=1:numel(sessN)
                B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
            end
            for r = reg
                for st=1:2 % per seqType
                    tstart = tic;
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
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
                                condVec{p}  = D.seqNumb(indx==1,:); % conditions
                                partVec{p}  = D.run(indx==1,:);
                                Data{p}     = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                            else
                                condVec{p}  = [condVec{p}; D.seqNumb(indx==1,:)+6]; % treat 2nd session as additional conditions
                                partVec{p}  = [partVec{p}; D.run(indx==1,:)+8];  % runs/partitions of 2nd session as additional runs
                                Data{p}     = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                            end;
                        end; % session
                    end; % subj
                    % group fit
                    % create a separate case - loads Data, M, partVec,
                    % condVec, runEffect,algorithm
                    T = sml_plasticity('PCM:corrModel_run','Data',Data,'model',M,'partV',partVec,'condV',condVec,'algorithm',algorithm,'runEffect',runEffect);
                    % other variables of interest
                    T.roi                   = ones(size(T.SN))*r;
                    T.regType               = ones(size(T.SN))*t.regType;
                    T.regSide               = ones(size(T.SN))*t.regSide;
                    T.seqType               = ones(size(T.SN))*st;
                    T.sessTr                = ones(size(T.SN))*str;
                    T = rmfield(T,'reg');
                    AllReg = addstruct(AllReg,T);
                    fprintf('Done modelType: %s seqType %d sess: %d-%d reg: %d/%d\n\n',withinCov,st,sessN(1),sessN(2),r,length(reg));
                    toc(tstart);
                end; % seqType
                fprintf('Done reg %d/%d:\tmodelType: %s \tsess: %d-%d\n\n',r,numel(reg),withinCov,sessN(1),sessN(2)); 
            end; % region
            save(fullfile(stabDir,sprintf('PCM_corrModels_%s_nonlinear_%s_%s_sessTrans-%d-%d.mat',parcelType,withinCov,sessType,sessN(1),sessN(2))),'-struct','AllReg');
            fprintf('Done all:\tmodelType: %s \tsess: %d-%d\n\n\n\n',withinCov,sessN(1),sessN(2)); 
            AllSess = addstruct(AllSess,AllReg);
        end; % session transition
        % save output
        save(fullfile(stabDir,sprintf('PCM_corrModels_%s_nonlinear_%s_%s.mat',parcelType,withinCov,sessType)),'-struct','AllSess');
    case 'PCM:corrModels_within'
        % nonlinear models within session (even - odd runs)
        runEffect   = 'random';
        beta_choice = 'mw';
        algorithm   = 'NR'; % minimize or NR
        reg         = 1:16;
        hemi        = 1;
        parcelType  = 'Brodmann';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        %sessType    = 'transitions'; % options: relativeto1 or transitions
        sessN       = 1:4;      % session
        corrLim     = [0 1];    % bounds for lower / upper correlation of models - change this to [-1 1]!!!
        nModel      = 30;       % number of correlation models (determines how fine grained the corr estimates are)
        withinCov   = 'individual'; % covariance type: individual or iid (or meanPattern)
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessTr','algorithm','parcelType','withinCov','nModel','sessType'});
        AllSess=[];
        if any(strcmp(parcelType,{'tesselsWB_162','tesselsWB_362'}))
            reg = sml_plasticity('TESSEL:select','sessN',1:4,'betaChoice','multi','hemi',hemi)';
        end
        % construct models
        corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
        for c=1:length(corrS)
            M{c} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',2,'numItems',6,'r',corrS(c));
        end
        for ss = sessN % session transition
            AllReg = [];
            B = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
            for r = reg
                for st=1:2 % per seqType
                    tstart = tic;
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        t = getrow(B,B.SN==sn(p) & B.region==r);
                        switch (beta_choice)
                            case 'uw'
                                beta = t.betaUW{:};
                            case 'mw'
                                beta = t.betaW{:};
                            case 'raw'
                                beta = t.betaRAW{:};
                        end
                        indx = D.seqType==st;
                        % here split into odd and even
                        Data{p} = [beta(mod(D.run,2)==1 & indx==1,:);beta(mod(D.run,2)==0 & indx==1,:)];
                        condVec{p} = [D.seqNumb(mod(D.run,2)==1 & indx==1);D.seqNumb(mod(D.run,2)==0 & indx==1)+6];
                        partVec{p} = kron(1:8,ones(1,6))';
                    end; % subj
                    % group fit
                    T = sml_plasticity('PCM:corrModel_run','Data',Data,'model',M,'partV',partVec,'condV',condVec,'algorithm',algorithm,'runEffect',runEffect);
                    % other variables of interest
                    T.roi                   = ones(size(T.SN))*r;
                    T.regType               = ones(size(T.SN))*t.regType;
                    T.regSide               = ones(size(T.SN))*t.regSide;
                    T.seqType               = ones(size(T.SN))*st;
                    T.sessN                 = ones(size(T.SN))*ss;
                    T = rmfield(T,'reg');
                    AllReg = addstruct(AllReg,T);
                    fprintf('Done modelType: %s seqType %d sess: %d reg: %d/%d\n\n',withinCov,st,ss,r,length(reg));
                    toc(tstart);
                end; % seqType
                fprintf('Done reg %d/%d:\tmodelType: %s \tsess: %d\n\n',r,numel(reg),withinCov,ss); 
            end; % region
            save(fullfile(stabDir,sprintf('PCM_corrModels_%s_nonlinear_%s_withinSess_sessN-%d.mat',parcelType,withinCov,ss)),'-struct','AllReg');
            fprintf('Done all:\tmodelType: %s \tsessN: %d\n\n\n\n',withinCov,ss); 
            AllSess = addstruct(AllSess,AllReg);
        end; % session 
        % save output
        save(fullfile(stabDir,sprintf('PCM_corrModels_%s_nonlinear_%s_withinSess.mat',parcelType,withinCov)),'-struct','AllSess');
    case 'PCM:corrModels_allSeq'
        % here assessing models of different correlation values
        % across sessions - use all sequences together (not trained vs.
        % control)
        % used for sessions 3-4
        runEffect   = 'random';
        beta_choice = 'mw';
        algorithm   = 'NR'; % minimize or NR
        reg         = 1:16;
        hemi        = 1;
        parcelType  = 'Brodmann';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessType    = 'transitions'; % options: relativeto1 or transitions
        sessTr      = 3;      % transitions
        corrLim     = [0 1];    % bounds for lower / upper correlation of models - change this to [-1 1]!!!
        nModel      = 30;       % number of correlation models (determines how fine grained the corr estimates are)
        withinCov   = 'individual'; % covariance type: individual or iid (or meanPattern)
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessTr','algorithm','parcelType','withinCov','nModel','sessType'});
        AllSess=[];
        if any(strcmp(parcelType,{'tesselsWB_162','tesselsWB_362'}))
            reg = sml_plasticity('TESSEL:select','sessN',1:4,'betaChoice','multi','hemi',hemi)';
        end
        % construct models
        corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
        for c=1:length(corrS)
            M{c} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',2,'numItems',12,'r',corrS(c));
        end

        sessN = [3 4];
        for str = sessTr % session transition
            AllReg = [];
            B = cell(length(sessN),1);
            for ss=1:numel(sessN)
                B{ss}=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,sessN(ss))));
            end
            for r = reg
                tstart = tic;
                condVec = cell(1,length(sn));
                partVec = condVec; Data = condVec; % initialise
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
                        indx = D.seqType>0;
                        if ss == 1
                            condVec{p}  = D.seqNumb(indx==1,:); % conditions
                            partVec{p}  = D.run(indx==1,:);
                            Data{p}     = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                        else
                            condVec{p}  = [condVec{p}; D.seqNumb(indx==1,:)+12]; % treat 2nd session as additional conditions
                            partVec{p}  = [partVec{p}; D.run(indx==1,:)+8];  % runs/partitions of 2nd session as additional runs
                            Data{p}     = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                        end;
                    end; % session
                end; % subj
                T = sml_plasticity('PCM:corrModel_run','Data',Data,'model',M,'partV',partVec,'condV',condVec,'algorithm',algorithm,'runEffect',runEffect);
                % other variables of interest
                T.roi                   = ones(size(T.SN))*r;
                T.regType               = ones(size(T.SN))*t.regType;
                T.regSide               = ones(size(T.SN))*t.regSide;
                T.sessTr                = ones(size(T.SN))*str;
                T = rmfield(T,'reg');
                AllReg = addstruct(AllReg,T);
                fprintf('Done modelType: %s sess: %d-%d reg: %d/%d\n\n',withinCov,sessN(1),sessN(2),r,length(reg));
                toc(tstart);
                fprintf('Done reg %d/%d:\tmodelType: %s \tsess: %d-%d\n\n',r,numel(reg),withinCov,sessN(1),sessN(2));
            end; % region
            fprintf('Done all:\tmodelType: %s \tsess: %d-%d\n\n\n\n',withinCov,sessN(1),sessN(2));
            AllSess = addstruct(AllSess,AllReg);
        end; % session transition
        % save output
        save(fullfile(stabDir,sprintf('PCM_corrModels_%s_allSeq_nonlinear_%s_%s.mat',parcelType,withinCov,sessType)),'-struct','AllReg');
    case 'PCM:corrModel_run'
        algorithm='NR';
        vararginoptions(varargin,{'Data','model','partV','condV','runEffect','algorithm'});
        T = pcm_fitModelGroup(Data,model,partV,condV,'runEffect',runEffect,'fitAlgorithm',algorithm,'fitScale',1); % only group fit
        T.bayesEst = bsxfun(@minus,T.likelihood,T.likelihood(:,1)); % now relative to model with 0 correlation
        T.mean_likelihood = mean(T.likelihood,2);
        %T.norm_likelihood = bsxfun(@minus,T.likelihood,T.mean_likelihood); % normalised to mean
        T.posterior = exp(T.bayesEst);
        T.posterior = bsxfun(@rdivide,T.posterior,sum(T.posterior,2));
        % individual fit
%         I = pcm_fitModelIndivid(Data,model,partV{1},condV{1},'runEffect',runEffect);
%         T.individ_likelihood    = I.likelihood;
%         T.bayesEst_individ      = bsxfun(@minus,T.individ_likelihood,T.individ_likelihood(:,1));
%         T.mean_like_individ     = mean(I.likelihood,2);
%         T.norm_like_individ     = bsxfun(@minus,I.likelihood,T.mean_like_individ); % normalised to mean
%         T.posterior_individ     = exp(T.bayesEst_individ);
%         T.posterior_individ     = bsxfun(@rdivide,T.posterior_individ,sum(T.posterior_individ,2));
        varargout{1}=T;
    case 'PCM:corr_plotflatmap'
        % project correlation from the winner correlation model
        hemi        = 1;
        parcelType  = 'tesselsWB_162';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessType    = 'transitions'; % options: relativeto1 or transitions
        sessTr      = 1:2;      % transitions
        modelType   = 'nonlinear_individual'; % covariance type: individual or iid (or meanPattern)
        vararginoptions(varargin,{'parcelType','sessType'});
        T = load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        surfaceGroupDir = fullfile(wbDir,'FS_LR_164');
        for h=hemi
            C=gifti(fullfile(surfaceGroupDir,sprintf('Icosahedron-162.164k.%s.label.gii',hemI{h}))); % freesurfer
            Data1 = zeros(size(C.cdata,1),2*numel(sessTr)); 
            T1 = getrow(T,T.regSide==h);
            for r=unique(T1.roi)'
                indx = C.cdata == r; % indicate all the surface nodes
                col1=1;
                for str=sessTr % all transitions
                    for seqType=1:2 % trained / untrained
                        t = getrow(T1,T1.seqType==seqType & T1.sessTr==str & T1.regType==r);
                        [ev1,maxCorr]=max(t.bayesEst,[],2);
                        [ev2,~]=min(t.bayesEst,[],2);
                        maxCorr = mean(maxCorr)/30;
                        eviCorr = mean(ev1-ev2);
                        Data1(indx,col1)=maxCorr;
                        Data2(indx,col1)=eviCorr;
                        col1 = col1+1;
                    end; % st
                end; %str
            end; % tessels
            % arrange border files
            borderName = {'CS','IPS','PoCS','SF'};
            for b=1:length(borderName);
                surface = fullfile(surfaceGroupDir,'fs_LR.164k.L.flat.surf.gii');
                border  = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.L.border-%s.border',borderName{b}));
                outname = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.L.border-%s.func.gii',borderName{b}));
                system(['wb_command -border-to-vertices ' surface ' ' border ' ' outname]);
                allBorder = gifti(outname);
                allB(:,b) = allBorder.cdata;
            end
            G = gifti(fullfile(surfaceGroupDir,'fs_LR.164k.L.flat.surf.gii'));
            figure      % rgb    
            for met=1:4
                data = [Data1(:,met) zeros(size(Data1,1),2)];
                borderFile = gifti(outname);
                %  transp = ones(size(data,1),1); transp(Data2(:,1)<2)=Data2(Data2(:,1)<2)./2;
                %transp = ones(size(data,1),1); transp(Data2(:,1)<2)=.5;
                data(Data2(:,1)<2)=0;
                subplot(2,2,met)
                caret_plotflatmap_rgb('coord',G.vertices,'topo',double(G.faces),'data',data,'xlims',[-90 200],'ylims',[-35 165],...
                    'alpha',ones(size(data,1),1),'dscale',[.3 .8],'border',allB,'underlay',double(G.faces));
              
             %   caret_plotflatmap_rgb('coord',G.vertices,'topo',double(G.faces),'data',data,'xlims',[-90 200],'ylims',[-35 165],...
             %       'alpha',transp,'dscale',[.3 .8],'border',allB,'underlay',double(G.faces));
                axis off;
            end
            figure      % not rgb  
            for met=1:4
                data = Data1(:,met);
                borderFile = gifti(outname);
                %  transp = ones(size(data,1),1); transp(Data2(:,1)<2)=Data2(Data2(:,1)<2)./2;
                data(Data2(:,1)<2)=0;
                subplot(2,2,met)
                caret_plotflatmap('coord',G.vertices,'topo',double(G.faces),'data',data,'xlims',[-90 200],'ylims',[-35 165],'border',allB)
                axis off;
            end
            
            statsFile = gifti(fullfile(surfaceGroupDir,sprintf('%s.PCMcorr_Tstats_crossval.%s.func.gii',hemI{h},sessType)));
            figure % difference rgb
            for m2=1:2
                data = [zeros(size(Data1,1),1) Data1(:,(m2-1)*2+2)-Data1(:,(m2-1)*2+1) zeros(size(Data1,1),1)];
                data(abs(statsFile.cdata(:,m2))<=1.3,:)=0;
                data(Data2(:,2)<2,:)=0;
                %data(data<0)=0;
                %data(Data2(:,m2)<2)=0; 
                transp = ones(size(data,1),1); %transp(Data2(:,1)<2)=.5; 
                %transp(abs(statsFile.cdata(:,m2))<=1.3)=.5;
                transp(Data2(:,1)<2)=.5; 
                subplot(1,2,m2)
                caret_plotflatmap_rgb('coord',G.vertices,'topo',double(G.faces),'data',data,'xlims',[-90 200],'ylims',[-35 165],'alpha',(transp),'dscale',[0 .35],'border',allB)
                axis off;
            end
            
        end; % hemisphere
    case 'PCM:corr_projectWinnerToSurface'
        % project correlation from the winner correlation model
        hemi        = 1;
        parcelType  = 'tesselsWB_162';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessType    = 'transitions'; % options: relativeto1 or transitions
        sessTr      = 1:2;      % transitions
        modelType   = 'nonlinear_individual'; % covariance type: individual or iid (or meanPattern)
        vararginoptions(varargin,{'parcelType','sessType'});
        T = load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        surfaceGroupDir = fullfile(wbDir,'FS_LR_164');
        for h=hemi
            C=gifti(fullfile(surfaceGroupDir,sprintf('Icosahedron-162.164k.%s.label.gii',hemI{h}))); % freesurfer
            Data1 = zeros(size(C.cdata,1),2*numel(sessTr)); 
            Data2 = zeros(size(C.cdata,1),numel(sessTr));  % difference
            T1 = getrow(T,T.regSide==h);
            for r=unique(T1.roi)'
                indx = C.cdata == r; % indicate all the surface nodes
                col1=1; col2=1;
                for str=sessTr % all transitions
                    for seqType=1:2 % trained / untrained
                        t = getrow(T1,T1.seqType==seqType & T1.sessTr==str & T1.regType==r);
                        [~,maxCorr]=max(t.bayesEst,[],2);
                        maxCorr = mean(maxCorr)/30;
                        Data1(indx,col1)=maxCorr;
                        column_name{col1} = sprintf('%s_sessTr-%d_%s',sessType,str,seqTypeName{seqType});
                        col1 = col1+1;
                    end; % st
                        Data2(indx,col2) = Data1(indx,col1-1) - Data1(indx,col1-2);
                        col2 = col2+1;
                        column_name2{col2} = sprintf('%s_sessTr-%d_difference',sessType,str);
                end; %str
            end; % tessels
            G = surf_makeFuncGifti(Data1,'anatomicalStruct',hemname{h},'columnNames',column_name);
            outFile = fullfile(surfaceGroupDir,sprintf('%s.PCMcorr.%s.func.gii',hemI{h},sessType));
            save(G,outFile);
            G = surf_makeFuncGifti(Data2,'anatomicalStruct',hemname{h},'columnNames',column_name2);
            outFile = fullfile(surfaceGroupDir,sprintf('%s.PCMcorr.%s_seqType_difference.func.gii',hemI{h},sessType));
            save(G,outFile);
        end; % hemisphere
    case 'PCM:corr_projectWinnerToSurface_maskEvidence'
        % project correlation from the winner correlation model
        % mask - only plot models that have evidence difference between st
        % model >2
        hemi        = 1;
        parcelType  = 'tesselsWB_162';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessType    = 'transitions'; % options: relativeto1 or transitions
        sessTr      = 1:2;      % transitions
        modelType   = 'nonlinear_individual'; % covariance type: individual or iid (or meanPattern)
        vararginoptions(varargin,{'parcelType','sessType'});
        T = load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        surfaceGroupDir = fullfile(wbDir,'FS_LR_164');
        for h=hemi
            C=gifti(fullfile(surfaceGroupDir,sprintf('Icosahedron-162.164k.%s.label.gii',hemI{h}))); % freesurfer
            Data1 = zeros(size(C.cdata,1),2*numel(sessTr)); 
            Data2 = zeros(size(C.cdata,1),numel(sessTr));  % difference
            T1 = getrow(T,T.regSide==h);
            for r=unique(T1.roi)'
                indx = C.cdata == r; % indicate all the surface nodes
                col1=1; col2=1;
                % here decide on transparency
                for str=2
                    t = getrow(T1,T1.sessTr==str & T1.regType==r);
                    [ev1,~]=max(t.bayesEst,[],2);
                    [ev2,~]=min(t.bayesEst,[],2);
                    eviCorr = mean(ev1-ev2); 
                end
                for str=sessTr % all transitions
                    for seqType=1:2 % trained / untrained
                        t = getrow(T1,T1.seqType==seqType & T1.sessTr==str & T1.regType==r);
                        [~,maxCorr]=max(t.bayesEst,[],2);
                        maxCorr = mean(maxCorr)/30;
                        if eviCorr<1
                            Data1(indx,col1)=0;
                        else
                            Data1(indx,col1)=maxCorr;
                        end
                        column_name{col1} = sprintf('%s_sessTr-%d_%s',sessType,str,seqTypeName{seqType});
                        col1 = col1+1;
                    end; % st
                    if eviCorr>1
                        Data2(indx,col2) = Data1(indx,col1-2) - Data1(indx,col1-1);
                    else
                        Data2(indx,col2) = 0;
                    end
                    col2 = col2+1;
                    column_name2{col2} = sprintf('%s_sessTr-%d_difference',sessType,str);
                end; %str
            end; % tessels
            G = surf_makeFuncGifti(Data1,'anatomicalStruct',hemname{h},'columnNames',column_name);
            outFile = fullfile(surfaceGroupDir,sprintf('%s.PCMcorr_maskedByEvidence.%s.func.gii',hemI{h},sessType));
            save(G,outFile);
            G = surf_makeFuncGifti(Data2,'anatomicalStruct',hemname{h},'columnNames',column_name2);
            outFile = fullfile(surfaceGroupDir,sprintf('%s.PCMcorr.%s_seqType_difference_maskedByEvidence.func.gii',hemI{h},sessType));
            save(G,outFile);
        end; % hemisphere
    case 'PCM:corr_projectStatsToSurface'
       % project t values
       % project correlation from the winner correlation model
       hemi        = 1;
       parcelType  = 'tesselsWB_162';
       %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
       sessType    = 'transitions'; % options: relativeto1 or transitions
       sessTr      = 1:2;      % transitions
       modelType   = 'nonlinear_individual'; % covariance type: individual or iid (or meanPattern)
       vararginoptions(varargin,{'parcelType','sessType'});
       T = load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
       surfaceGroupDir = fullfile(wbDir,'FS_LR_164');
       for h=hemi
           C=gifti(fullfile(surfaceGroupDir,sprintf('Icosahedron-162.164k.%s.label.gii',hemI{h}))); % freesurfer
           Data1 = zeros(size(C.cdata,1),numel(sessTr));
           T1 = getrow(T,T.regSide==h);
           for r=unique(T1.roi)'
               indx = C.cdata == r; % indicate all the surface nodes
               col1=1; 
               for str=sessTr % all transitions
                   t1 = getrow(T1,T1.regType==r & T1.seqType==1 & T1.sessTr==str);
                   t2 = getrow(T1,T1.regType==r & T1.seqType==2 & T1.sessTr==str);
                   [~,j1]=max(t1.bayesEst,[],2);
                   [~,j2]=max(t2.bayesEst,[],2);
                   maxS1 = round(mean(j1));
                   maxS2 = round(mean(j2));
                   [t,~]=ttest(t1.bayesEst(:,maxS1),t1.bayesEst(:,maxS2),2,'paired');
                   Data1(indx,col1)=t;
                   column_name{col1} = sprintf('%s_sessTr%d_Tstatistic',sessType,str);
                   col1 = col1+1;
               end; %str
           end; % tessels
           G = surf_makeFuncGifti(Data1,'anatomicalStruct',hemname{h},'columnNames',column_name);
           outFile = fullfile(surfaceGroupDir,sprintf('%s.PCMcorr_Tstats.%s.func.gii',hemI{h},sessType));
           save(G,outFile);
       end; % hemisphere
    case 'PCM:corr_projectCrossStatsToSurface'
        % project t values - crossvalidated
       % project correlation from the winner correlation model
       hemi        = 1;
       allSubj     = 1:26; 
       parcelType  = 'tesselsWB_162';
       %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
       sessType    = 'transitions'; % options: relativeto1 or transitions
       sessTr      = 1:2;      % transitions
       modelType   = 'nonlinear_individual'; % covariance type: individual or iid (or meanPattern)
       metric = 'bayesEst'; % bayesEst, or bayesEst_individ
       vararginoptions(varargin,{'parcelType','sessType'});
       T = load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
       surfaceGroupDir = fullfile(wbDir,'FS_LR_164');
       for h=hemi
           C=gifti(fullfile(surfaceGroupDir,sprintf('Icosahedron-162.164k.%s.label.gii',hemI{h}))); % freesurfer
           Data1 = zeros(size(C.cdata,1),numel(sessTr));
           T1 = getrow(T,T.regSide==h);
           for r=unique(T1.roi)'
               NN=[];
               indx = C.cdata == r; % indicate all the surface nodes
               col1=1;
               for str=sessTr % all transitions
                   t = getrow(T1,T1.regType==r & T1.sessTr==str);
                   t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
                   for s=1:26
                       testSubj = allSubj(allSubj~=s)';
                       [~,corT1]=max(t.metric2(ismember(t.SN,testSubj)&t.sessTr==str&t.seqType==1,:),[],2);
                       [~,corU1]=max(t.metric2(ismember(t.SN,testSubj)&t.sessTr==str&t.seqType==2,:),[],2);
                       % here get the new values
                       N.val(1,:) = t.metric2(t.SN==s&t.seqType==1,round(mean(corT1)));
                       N.val(2,:) = t.metric2(t.SN==s&t.seqType==1,round(mean(corU1)));
                       N.seqType = [1;2];
                       N.sn = [s;s];
                       NN=addstruct(NN,N);
                   end
                   t=ttestDirect(NN.val,[NN.seqType NN.sn],2,'paired');
                 %  [t,~]=ttest(t1.bayesEst(:,maxS1),t1.bayesEst(:,maxS2),2,'paired');
                   Data1(indx,col1)=t;
                   column_name{col1} = sprintf('%s_sessTr%d_Tstatistic',sessType,str);
                   col1 = col1+1;
               end; %str
           end; % tessels
           G = surf_makeFuncGifti(Data1,'anatomicalStruct',hemname{h},'columnNames',column_name);
           outFile = fullfile(surfaceGroupDir,sprintf('%s.PCMcorr_Tstats_crossval.%s.func.gii',hemI{h},sessType));
           save(G,outFile);
       end; % hemisphere
    case 'PCM:corr_maskDifference_byStats'
        % here get the difference functional gifti and mask by crossval T
        % values >1.7
        hemi        = 1;
        sessType    = 'transitions'; % options: relativeto1 or transitions
        surfaceGroupDir = fullfile(wbDir,'FS_LR_164');
        for h=hemi
            dF = gifti(fullfile(surfaceGroupDir,sprintf('%s.PCMcorr.%s_seqType_difference_maskedByEvidence.func.gii',hemI{h},sessType)));
            statF = gifti(fullfile(surfaceGroupDir,sprintf('%s.PCMcorr_Tstats_crossval.%s.func.gii',hemI{h},sessType)));
            for i=1:size(dF.cdata,2)
                data(:,i) = dF.cdata(:,i);
                data(abs(statF.cdata(:,i))<1.3,i)=0;
                colName{i} = sprintf('difference_trans-%d',i);
            end
            G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',colName);
            outFile = (fullfile(surfaceGroupDir,sprintf('%s.PCMcorr.%s_seqType_difference_maskedByEvidence_stats.func.gii',hemI{h},sessType)));
            save(G,outFile);
        end
    case 'PCM:simulate'
        modelType = 'noSess'; % noSess or wSess
        corrM = 0.6; % correlation simulated
        noise = 0.0001;
        partN = 8;
        voxN  = 100;
        nSubj = 10;
        nModel = 11; % determines the step size for corrModels
        runEffect = 'random';
        vararginoptions(varargin,{'modelType','corrM','noise'})
        
        signal=1;
        % first get the correct model
        switch modelType
            case 'noSess'
                M_gen = pcm_corrModel_noSess([corrM corrM],1);
                M_gen = M_gen{2};
                theta = ones(M_gen.numGparams,1)*0.1;
            case 'wSess'
                M_gen = pcm_corrModel_wSess([corrM corrM],1);
                M_gen = M_gen{3};
                theta = ones(M_gen.numGparams,1)*0.1;
        end
        % now generate data
        D.numPart = partN;
        D.numVox  = voxN;
        [Data,partVec,condVec]=pcm_generateData(M_gen,theta,D,nSubj,signal,noise);
                      
        % here fit models
        % 1) fit the toy models to estimate correlation
        % 2) fit the correlation models
        switch modelType
            case 'wSess'
                M_est1  = pcm_toyModel_wSess;
                C       = pcm_correlation(Data,{partVec},{condVec},M_est1{4},runEffect,modelType);
                M_est2  = pcm_corrModel_wSess([0 1],nModel);
            case 'noSess'
                M_est1  = pcm_toyModel_noSess;
                C       = pcm_correlation(Data,{partVec},{condVec},M_est1{3},runEffect,modelType);
                M_est2  = pcm_corrModel_noSess([0 1],nModel);
        end
        % fit correlation models
        T = pcm_fitModels(Data,M_est2,partVec,condVec,runEffect,'NR');
        % individual fit
        I = pcm_fitModelIndivid(Data,M_est2,partVec,condVec,'runEffect',runEffect);
        T.individ_likelihood    = I.likelihood;
        T.bayesEst_individ      = bsxfun(@minus,T.individ_likelihood,T.individ_likelihood(:,1));
        
        [~,~,G_pred,~] = pcm_fitModelGroup(Data,M_est2,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR');
        figure
        corrPred=zeros(11,1);
        for c=1:size(T.iterations,2)
            subplot(3,5,c)
            imagesc(G_pred{c});
            corrPred(c)=calcCorr(G_pred{c});
            title(sprintf('model - correlation %1.2f',corrPred(c)));
        end
        figure
        subplot(4,4,1:4)
        plt.dot(repmat(ones(nSubj,1),3,1),[C.r_naive;C.r_crossval;C.r_model],...
            'split',[ones(size(C.r_naive));ones(size(C.r_naive))*2;ones(size(C.r_naive))*3],'leg',{'naive','crossval','pcm'});
        ylabel('Correlation');
        title(sprintf('Estimated correlation for model with %1.2f corr',corrM));
        subplot(4,4,[5,9,13]);
        bar(mean(T.bayesEst(1,2:end),1)); hold on; [~,p]=max(T.bayesEst,[],2); drawline(p,'dir','vert');
        set(gca,'XTickLabel',0:0.1:1);
        ylabel('Bayes factor'); title('Group fit');
        subplot(4,4,[6,10,14]);
        bar(mean(T.bayesEst(1,2:9),1)); hold on; [~,p]=max(T.bayesEst,[],2); drawline(p,'dir','vert');
        set(gca,'XTickLabel',0:0.1:0.8); title('Group fit - subset');
        subplot(4,4,[7,11,15])
        bar(mean(T.bayesEst_individ(1,2:end),1)); hold on; [~,p]=max(T.bayesEst_individ,[],2); drawline(p,'dir','vert');
        set(gca,'XTickLabel',0:0.1:1); title('Individ fit');
        subplot(4,4,[8,12,16]);
        bar(mean(T.bayesEst_individ(1,2:9),1)); hold on; [~,p]=max(T.bayesEst_individ,[],2); drawline(p,'dir','vert');
        set(gca,'XTickLabel',0:0.1:0.8);
        set(gca,'XTickLabel',0:0.1:1); title('Individ fit - subset');
        keyboard;
    
    case 'PCM:plot_withinSess' % OLD
        reg         = 1:8;
        hemi        = 1;
        parcelType  = 'Brodmann'; 
        vararginoptions(varargin,{'parcelType','reg','hemi'});
        T = load(fullfile(stabDir,sprintf('PCM_withinSess_%s.mat',parcelType)));    
    
        for r=1:length(reg)
            for ss=sessN
                figure(r)
                subplot(1,4,ss)
                barplot(T.seqType,T.bayesEst_cross,'subset',T.roi==reg(r)&T.regSide==hemi&T.sessN==ss);
                plt.match('y'); ylabel('logBayes');
                title(sprintf('%s - sess%d',regname_cortex{reg(r)},ss));
            end
            figure(max(reg)+r)
            subplot(121)
            barplot(T.sessN,T.bayesEst_cross,'subset',T.roi==reg(r)&T.regSide==hemi&T.seqType==1);
            title(sprintf('%s - trained',regname_cortex{reg(r)})); ylabel('logBayes');
            subplot(122)
            barplot(T.sessN,T.bayesEst_cross,'subset',T.roi==reg(r)&T.regSide==hemi&T.seqType==2);
            title(sprintf('%s - untrained',regname_cortex{reg(r)})); ylabel('logBayes');
            plt.match('y');
            figure(99)
            subplot(1,numel(reg),r)
            style.use('Seq');
            plt.line([T.sessN>3 T.sessN],T.bayesEst_cross(:,4),'split',T.seqType,'subset',T.roi==reg(r) & T.regSide==hemi,'leg',{'trained','untrained'});
            if r==1
                ylabel('Sequence-specific log-Bayes');
            else
                ylabel('');
            end
            xlabel('Session'); title(sprintf('%s',regname_cortex{reg(r)}));
        end
    case 'PCM:stats_withinSess_old'
        roi=1:8;
        parcelType='Brodmann';
        hemi=1;
        modelInd = 3; % 2 or 3: 2 all distance equal, 3 flexible
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        T = load(fullfile(stabDir,sprintf('PCM_withinSess_%s.mat',parcelType)));    
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\ndistance logBayes in %s across sess:1-4\n',regLab{r});
            anovaMixed(T.bayesEst_cross(:,modelInd),T.SN,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi & T.regSide==hemi);
            fprintf('\ndistance logBayes in %s across sess:1-3\n',regLab{r});
            anovaMixed(T.bayesEst_cross(:,modelInd),T.SN,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi & T.regSide==hemi & T.sessN<4);
        end
        keyboard;
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regLab{r});
                ttestDirect(T.bayesEst_cross(:,modelInd),[T.seqType T.SN],2,'paired','subset',T.regType==r & T.regSide==hemi & T.sessN==ss);
            end
            for str = 1:3 % session transition
                fprintf('\n post-hoc t-test on the effect of session transition %d-%d TRAINED in %s \n',str,str+1,regLab{r});
                ttestDirect(T.bayesEst_cross(:,modelInd),[T.sessN T.SN],2,'paired','subset',T.regType==r & T.regSide==hemi & T.seqType==1 & ismember(T.sessN,[str str+1]));
                fprintf('\n post-hoc t-test on the effect of session transition %d-%d UNTRAINED in %s \n',str,str+1,regLab{r});
                ttestDirect(T.bayesEst_cross(:,modelInd),[T.sessN T.SN],2,'paired','subset',T.regType==r & T.regSide==hemi & T.seqType==2 & ismember(T.sessN,[str str+1]));
            end
        end
    case 'PCM:plot_toyModels'
        reg         = 1:8;
        hemi        = 1;
        parcelType  = 'Brodmann';
        modelType   = 'noSess';
       
        T = load(fullfile(stabDir,sprintf('PCM_toyModels_%s_%s.mat',parcelType,modelType)));
        T = rmfield(T,'reg');
        T.noSess        = bsxfun(@minus,T.bayesEst,T.bayesEst(:,2));
        T.noSess_cross  = bsxfun(@minus,T.bayesEst_cross,T.bayesEst_cross(:,2));
        for r=reg
            figure
            for str = 1:3 % transition
                subplot(4,3,str)
                barplot(T.seqType,T.bayesEst,'subset',T.roi==r&T.regSide==hemi&T.sessTr==str);
                ylabel('logBayes'); title(sprintf('%s-sess:%d-%d',regname_cortex{r},str,str+1));
                subplot(4,3,str+3)
                barplot(T.seqType,T.bayesEst_cross,'subset',T.roi==r&T.regSide==hemi&T.sessTr==str);
                ylabel('crossval-logBayes'); title(sprintf('%s-sess:%d-%d',regname_cortex{r},str,str+1));
                subplot(4,3,str+6)
                barplot(T.seqType,T.noSess(:,2:end-1),'subset',T.roi==r&T.regSide==hemi&T.sessTr==str);
                ylabel('crossval-logBayes'); title(sprintf('%s-sess:%d-%d',regname_cortex{r},str,str+1));
                subplot(4,3,str+9)
                barplot(T.seqType,T.noSess_cross(:,2:end-1),'subset',T.roi==r&T.regSide==hemi&T.sessTr==str);
                ylabel('crossval-logBayes'); title(sprintf('%s-sess:%d-%d',regname_cortex{r},str,str+1));
            end
        end
    case 'PCM:plot_flexCorrelation'
        reg         = 1:8;
        hemi        = 1;
        parcelType  = 'Brodmann';
        modelType   = 'wSess';
        metric      = 'r_model'; % r_model, r_crossval, r_cross_seqType, r_model_seqType
        vararginoptions(varargin,{'modelType','reg','hemi','metric'});
        
        T = load(fullfile(stabDir,sprintf('PCM_toyModels_%s_%s.mat',parcelType,modelType)));
        if isfield(T,'reg')
            T = rmfield(T,'reg');
        end
        style.use('Seq');
        for r=reg
            t=getrow(T,T.roi==r & T.regSide==hemi);
            figure
            %plt.box(t.sessTr,t.(metric),'split',t.seqType,'plotall',2,'leg',{'trained','control'});
            plt.line(t.sessTr,t.(metric),'split',t.seqType,'leg',{'trained','control'});
%             figure
%             plt.dot(t.sessTr,t.(metric),'split',t.seqType,'leg',{'trained','control'},'subset',t.sessTr==1);
%             xLim = get(gca,'XTick');
%             hold on
%             for i=unique(t.SN)'
%                 plot(xLim,[t.(metric)(t.sessTr==1 & t.SN==i & t.seqType==1) t.(metric)(t.sessTr==1 & t.SN==i & t.seqType==2)],'-k');
%             end
             title(sprintf('%s',regname_cortex{r})); ylabel('correlation');      
        end
    case 'PCM:stats_flexCorrelation'
        reg         = 1:8;
        hemi        = 1;
        parcelType  = 'Brodmann';
        modelType   = 'wSess';
        metric      = 'r_model';
        vararginoptions(varargin,{'modelType','reg','hemi','metric'});
        
        T = load(fullfile(stabDir,sprintf('PCM_toyModels_%s_%s.mat',parcelType,modelType)));
        if isfield(T,'reg')
            T = rmfield(T,'reg');
        end
        for r=reg
            t=getrow(T,T.roi==r & T.regSide==hemi);
            fprintf('t-test %s - per session\n',regname_cortex{r});
            ttestDirect(t.(metric),[t.seqType t.SN],2,'paired','split',t.sessTr);
            fprintf('t-test %s - sessTransition (1-2) - trained\n',regname_cortex{r});
            ttestDirect(t.(metric),[t.sessTr t.SN],2,'paired','subset',t.sessTr~=3 & t.seqType==1);
            fprintf('t-test %s - sessTransition (1-2) - untrained\n',regname_cortex{r});
            ttestDirect(t.(metric),[t.sessTr t.SN],2,'paired','subset',t.sessTr~=3 & t.seqType==2);
        end
        
    case 'PCM:plot_corrModels_withinSess'
        modelType = 'nonlinear_individual';
        % for nonlinear models: nonlinear_indiviudal / nonlinear_individual, nonlinear_meanPattern
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, bayesEst_cross or bayesEst_individ
        %sessType = 'withinSess'; % transitions or relativeto1
        vararginoptions(varargin,{'modelType','parcelType','reg','metric','sessType'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_withinSess.mat',parcelType,modelType)));

        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
            % reshape
            nModel      = size(t.(metric),2);
            D.(metric)  = t.(metric)(:);
            D.metric2   = t.metric2(:);
            D.SN        = repmat(t.SN,nModel,1);
            D.sessN    = repmat(t.sessN,nModel,1);
            D.seqType   = repmat(t.seqType,nModel,1);
            D.model     = kron((1:nModel)',ones(size(t.(metric),1),1));               
            figure
            for ss=1:4 % all sessions
                subplot(1,4,ss)
                style.use('SeqShade_small')
                plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessN==ss); 
                hold on; drawline(0,'dir','horz');
                title(sprintf('%s sess-%d',regname_cortex{r},ss));
                xlabel('Model correlation'); ylabel('Evidence'); 
                t = gca;
                corrTick = linspace(min(t.XTick),max(t.XTick),11);
                set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
            end
            figure(99)
            subplot(1,numel(reg),find(reg==r))
            style.use('SeqShade_small')
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessN],'subset',ismember(D.sessN,[1,3]),'leg',{'trained','control'},'leglocation','southeast');
            hold on; drawline(0,'dir','horz');
            style.use('SeqShade_dotted')
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessN],'subset',ismember(D.sessN,[1,3]),'leg',{'trained','control'},'leglocation','southeast');
        end
    case 'PCM:plot_corrModels'
        modelType = 'nonlinear_individual';
        % for linear models: wSess
        % for nonlinear models: nonlinear_indiviudal / nonlinear_individual, nonlinear_meanPattern
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, bayesEst_cross or bayesEst_individ
        sessType = 'transitions'; % transitions or relativeto1
        vararginoptions(varargin,{'modelType','parcelType','reg','metric','sessType'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));

        if strcmp(modelType,'wSess')
            if isfield(T,'reg')
                T = rmfield(T,'reg');
            end
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
        end
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
            t.metric3 = bsxfun(@minus,t.(metric),max(t.(metric),2));
            % reshape
            nModel      = size(t.(metric),2);
            D.(metric)  = t.(metric)(:);
            D.metric2   = t.metric2(:);
            D.metric3   = t.metric3(:);
            D.SN        = repmat(t.SN,nModel,1);
            D.sessTr    = repmat(t.sessTr,nModel,1);
            D.seqType   = repmat(t.seqType,nModel,1);
            D.model     = kron((1:nModel)',ones(size(t.(metric),1),1));               
            anovaMixed(D.metric2,D.SN,'within',[D.model D.seqType],{'model','seqType'},'subset',D.sessTr==1 & ismember(D.model,[1:15:30])) 
            figure(96)
            subplot(1,numel(reg),find(r==reg))
            style.use('SeqShade_small')
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessTr],'subset',D.sessTr==1,'leg',{'trained','control'},'leglocation','southeast');
            hold on; drawline(0,'dir','horz');
            style.use('SeqShade_dotted')
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessTr],'subset',D.sessTr==2,'leg',{'trained','control'},'leglocation','southeast');
            figure(95)
            subplot(4,numel(reg),find(r==reg))
            style.use('SeqShade_small')
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessTr],'subset',D.sessTr==1 & D.seqType==1); hold on;  drawline(0,'dir','horz');
            title(regname_cortex{r})
            style.use('SessUntrained_shade');
            subplot(4,numel(reg),find(r==reg)+numel(reg))
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessTr],'subset',D.sessTr==1 & D.seqType==2);
            subplot(4,numel(reg),find(r==reg)+2*numel(reg))
            style.use('SeqShade_dotted')
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessTr],'subset',D.sessTr==2 & D.seqType==1);
            subplot(4,numel(reg),find(r==reg)+3*numel(reg))
            style.use('SessUntrained_dotted')
            plt.line(D.model,D.metric2,'split',[D.seqType D.sessTr],'subset',D.sessTr==2 & D.seqType==2);

            figure(94)
            subplot(1,numel(reg),find(r==reg))
            style.use('SeqShade_4');
            plt.line(D.model,D.metric2,'split',[D.sessTr D.seqType],'subset',D.sessTr==1); 
            hold on;
            style.use('SeqShade_light');
            plt.line(D.model,D.metric2,'split',[D.sessTr D.seqType],'subset',D.sessTr==2); legend off;
            drawline(0,'dir','horz');
            xlabel('Model correlation'); ylabel('Evidence'); title(regname_cortex{r});
            t = gca;
            corrTick = linspace(min(t.XTick),max(t.XTick),11);
            set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
            
            figure(92)
            subplot(1,numel(reg),find(r==reg))
            style.use('SeqShade_small');
            plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessTr==3);
            hold on;
            drawline(0,'dir','horz');
            xlabel('Model correlation'); ylabel('Evidence'); title(regname_cortex{r});
            t = gca;
            corrTick = linspace(min(t.XTick),max(t.XTick),11);
            set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
        end
    case 'PCM:plot_corrModels_topCorr'
        modelType = 'nonlinear_individual';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        sessType = 'transitions'; % transitions or relativeto1
        metric = 'bayesEst'; % bayesEst or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric','sessType'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        if strcmp(modelType,'wSess')
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
        end
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            [t.corVal,t.corIdx]=max(t.(metric),[],2);
            t=normData(t,'corVal');
            t=normData(t,'corIdx');
            figure
            subplot(121)
            style.use('SeqShade');
            plt.bar([t.sessTr>2 t.sessTr],t.normcorIdx./30,'split',t.seqType); % divide to get the true 'corr' value
            xlabel('session transition'); ylabel('maximum correlation'); title(regname_cortex{r});
            subplot(122)
            plt.bar([t.sessTr>2 t.sessTr],t.normcorVal,'split',t.seqType); % divide to get the true 'corr' value
            ylabel('evidence for max correlation');
            figure(99)
            style.use('Seq');
            subplot(1,numel(reg),find(r==reg))
            plt.line(t.sessTr,t.normcorIdx./30,'split',t.seqType,'subset',t.sessTr<3); title(regname_cortex{r});
        end
        keyboard;
        [T.corVal,T.corIdx]=max(T.(metric),[],2);
        style.use('4Reg');
        figure
        T1 = getrow(T,ismember(T.regType,reg) & T.seqType==1);
        T2 = getrow(T,ismember(T.regType,reg) & T.seqType==2);
        T1 = normData(T1,'corIdx');
        T2 = normData(T2,'corIdx');
        T = normData(T,'corIdx');
        plt.line([T2.sessTr>2 T2.sessTr],(T1.corIdx-T2.corIdx)./30,'subset',T1.regSide==1 & ismember(T1.regType,[1:3,7]),'split',T1.regType);
        hold on; drawline(0,'dir','horz');
        xlabel('Session transition'); ylabel('Trained - control correlation');
        figure
        subplot(121)
        plt.line([T1.sessTr>2 T1.sessTr],(T1.normcorIdx)./30,'subset',T1.regSide==1& ismember(T1.regType,[1:3,7]),'split',T1.regType);
        hold on; drawline(0,'dir','horz'); title('trained'); xlabel('Session transition'); ylabel('correlation');
        subplot(122) 
        plt.line([T1.sessTr>2 T1.sessTr],(T2.normcorIdx)./30,'subset',T1.regSide==1& ismember(T1.regType,[1:3,7]),'split',T1.regType);
        hold on; drawline(0,'dir','horz'); title('control'); xlabel('Session transition'); ylabel('correlation');
        plt.match('y');
        figure
        subplot(121)
        plt.line(T.seqType,(T.normcorIdx)./30,'subset',T.regSide==1& ismember(T.regType,[1:3,7])&T.sessTr==1,'split',T.regType);
        hold on; drawline(0,'dir','horz'); title('week 1-2'); xlabel('Sequence type'); ylabel('correlation');
        subplot(122) 
        plt.line(T.seqType,(T.normcorIdx)./30,'subset',T.regSide==1& ismember(T.regType,[1:3,7])&T.sessTr==2,'split',T.regType);
        hold on; drawline(0,'dir','horz'); title('week 2-5'); xlabel('Sequence type'); ylabel('correlation');
        plt.match('y');
        figure
        subplot(121)
        plt.line(T.seqType,T.corVal,'subset',T.regSide==1& ismember(T.regType,[1:3,7])&T.sessTr==1,'split',T.regType);
        hold on; drawline(0,'dir','horz'); title('week 1-2'); xlabel('Sequence type'); ylabel('logBayes');
        subplot(122) 
        plt.line(T.seqType,T.corVal,'subset',T.regSide==1& ismember(T.regType,[1:3,7])&T.sessTr==2,'split',T.regType);
        hold on; drawline(0,'dir','horz'); title('week 2-5'); xlabel('Sequence type'); ylabel('logBayes');
        plt.match('y');
        
        % just transitions 1-2
        figure
        subplot(121)
        plt.line(T1.sessTr,(T1.normcorIdx)./30,'subset',T1.regSide==1& ismember(T1.regType,[1:3,7]) & T1.sessTr<3,'split',T1.regType);
        hold on; drawline(0,'dir','horz'); title('trained'); xlabel('Session transition'); ylabel('correlation');
        subplot(122)
        style.use('4Reg_dashed');
        plt.line(T1.sessTr,(T2.normcorIdx)./30,'subset',T1.regSide==1& ismember(T1.regType,[1:3,7]) & T1.sessTr<3,'split',T1.regType);
        hold on; drawline(0,'dir','horz'); title('control'); xlabel('Session transition'); ylabel('correlation');
        plt.match('y');
    case 'PCM:plot_corrModels_surface'
        % project the results to the surface
        parcelType = 'tesselsWB_162';
        withinCov = 'individual';
        sessType = 'transitions';
        nNodes = 162;
        % general
        h=1;
        G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nNodes,hemI{h})));
        nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
        for str=1:3
            T = load(fullfile(stabDir,sprintf('PCM_corrModels_%s_nonlinear_%s_%s_sessTrans-%d-%d.mat',parcelType,withinCov,sessType,str,str+1)));
            for r=1:nReg
                t = getrow(T,T.regType==r&T.regSide==h);
                idx = find(G.cdata(:,1)==r);
                if sum(isnan(t.dist_train))/length(t.dist_train)>0.7% ensure more than 70% have data
                    data(idx,(ss-1)*2+1) = NaN;
                    data(idx,(ss-1)*2+2) = NaN;
                else
                    data(idx,(ss-1)*2+1) = nanmean(t.dist_train);
                    data(idx,(ss-1)*2+2) = nanmean(t.dist_untrain);
                end
                column_name{(ss-1)*2+1} = sprintf('dist_trained-sess%d',ss);
                column_name{(ss-1)*2+2} = sprintf('dist_untrained-sess%d',ss);
            end
        end
        outfile         = fullfile(wbDir,'FS_LR_164',sprintf('%s.tesselsWB_%d-dist.func.gii',hemI{h},nNodes));
        G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
        % filter out any crazy values (on the border of mask)
        for c=1:size(G.cdata,2)
            G.cdata(G.cdata(:,c)>5,c)=0;
            G.cdata(G.cdata(:,c)<-5,c)=0;
            save(G, outfile);
            fprintf('Hemisphere: %d\n',h);
        end
    case 'PCM:stats_withinSess'
        modelType = 'nonlinear_individual';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_withinSess.mat',parcelType,modelType)));
        T.metric2 = bsxfun(@minus,T.(metric),mean(T.(metric),2));
        allSubj = 1:26; 
        for r=reg
            NN=[];
            for ss=sessN
                t = getrow(T,T.regType==r & T.regSide==hemi & T.sessN==ss);
                for s=1:26
                    testSubj = allSubj(allSubj~=s)';
                    [~,corT1]=max(t.metric2(ismember(t.SN,testSubj)&t.seqType==1,:),[],2);
                    [~,corU1]=max(t.metric2(ismember(t.SN,testSubj)&t.seqType==2,:),[],2);
                    % here get the new values
                    N.val(1,:) = t.metric2(t.SN==s&t.seqType==1,round(mean(corT1)));
                    N.val(2,:) = t.metric2(t.SN==s&t.seqType==1,round(mean(corU1)));
                    N.winner(1,:) = t.metric2(t.SN==s&t.seqType==1,end);
                    N.winner(2,:) = t.metric2(t.SN==s&t.seqType==2,end);
                    N.seqType = [1;2];
                    N.sess = [ss;ss];
                    N.sn = [s;s];
                    NN=addstruct(NN,N);
                end
            end
            fprintf('Region %s:\ns',regname_cortex{r});
            % crossvalidated type of t-test (for evidence for best model)
            fprintf('Crossvalidated t-test:\n');
            ttestDirect(NN.val,[NN.seqType NN.sn],2,'paired','split',NN.sess);
            % here test evidence for perfect correlation model
            fprintf('Perfect correlation t-test:\n');
            ttestDirect(NN.winner,[NN.seqType NN.sn],2,'paired','split',NN.sess);
            [~,NN.maxCorr]=max(T.metric2(T.regType==r & T.regSide==hemi,:),[],2);
            fprintf('Maximum likelihood correlation t-test:\n');
            ttestDirect(NN.maxCorr,[NN.seqType NN.sn],2,'paired','split',NN.sess);
        end
    case 'PCM:stats_shape_corrModels'
        % new stats - assessing the evidence for best trained / untrained
        % model
        modelType = 'nonlinear_individual';
        sessType = 'transitions';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        for r=reg
            for str=1:3 
                fprintf('%s: transition-%d\n',regname_cortex{r},str);
                t1 = getrow(T,T.regType==r & T.regSide==hemi & T.seqType==1 & T.sessTr==str);
                t2 = getrow(T,T.regType==r & T.regSide==hemi & T.seqType==2 & T.sessTr==str);
                [~,j1]=max(t1.(metric),[],2);
                [~,j2]=max(t2.(metric),[],2);
                maxS1 = round(mean(j1));
                maxS2 = round(mean(j2));
                ttest(t1.(metric)(:,maxS1),t1.(metric)(:,maxS2),2,'paired')
                ttest(t2.(metric)(:,maxS1),t2.(metric)(:,maxS2),2,'paired')
                % do crossvalidated
                maxS1 = round(mean(j1(1:10)));
                maxS2 = round(mean(j2(1:10)));
                 maxS1 = round(mean(j1(14:end)));
                maxS2 = round(mean(j2(14:end)));
                ttest(t1.(metric)(t1.SN>10,maxS1),t1.(metric)(t1.SN>10,maxS2),2,'paired')
                ttest(t2.(metric)(t2.SN>10,maxS1),t2.(metric)(t2.SN>10,maxS2),2,'paired')
                 ttest(t1.(metric)(t1.SN<14,maxS1),t1.(metric)(t1.SN<14,maxS2),2,'paired')
                ttest(t2.(metric)(t2.SN<14,maxS1),t2.(metric)(t2.SN<14,maxS2),2,'paired')
                % restructure for ANOVAs
                t1.metric2 = bsxfun(@minus,t1.(metric),mean(t1.(metric),2));
                t2.metric2 = bsxfun(@minus,t2.(metric),mean(t2.(metric),2));
                t.(metric) = [t1.metric2(:);t2.metric2(:)];
                t.sn       = repmat(t1.SN,30*2,1);
                t.model    = repmat(kron([1:30]',ones(26,1)),2,1); % 1-6
                t.seqType  = [ones(size(t1.(metric)(:)));ones(size(t1.(metric)(:)))*2];
                anovaMixed(t.(metric),t.sn,'within',[t.model t.seqType],{'model','seqType'},'subset',ismember(t.model,[1,10,20,30]))
            end
            keyboard;
        end
    case 'PCM:stats_winner_corrModels'
        % assess the winning model
        modelType = 'nonlinear_individual';
        sessType = 'transitions';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        if strcmp(modelType,'wSess')
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
            T=rmfield(T,'reg');
        end
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
            t.metric3 = bsxfun(@minus,t.(metric),max(t.(metric),2));
            [corVal,corIdx]=max(t.(metric),[],2);
            fprintf('%s: ANOVA: seqType x sessTr\n',regname_cortex{r});
            anovaMixed(corIdx,t.SN,'within',[t.sessTr t.seqType],{'session','seqType'},'subset',t.sessTr~=3);
            fprintf('%s: t-tests (trained vs. control): corModel\n',regname_cortex{r});
            ttestDirect(corIdx,[t.seqType t.SN],2,'paired','split',t.sessTr);
            fprintf('%s: t-tests (sessTrans 1 vs. 2)\n',regname_cortex{r});
            ttestDirect(corIdx,[t.sessTr t.SN],2,'paired','split',t.seqType,'subset',t.sessTr~=3);
            fprintf('%s: t-tests (trained vs. control): evidence for corModel\n',regname_cortex{r});
            ttestDirect(corVal,[t.seqType t.SN],2,'paired','split',t.sessTr);
        end
    case 'PCM:stats_speed'
        % compare fits across sessions 3 and 4; and across regions
        modelType = 'nonlinear_individual';
        sessType = 'transitions';
        parcelType = 'Brodmann';
        reg = [1:3,7];
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        if strcmp(modelType,'wSess')
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
            T=rmfield(T,'reg');
        end
        VV = [];
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi & T.sessTr==3);
            t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
            [V.corVal,V.corIdx]=max(t.(metric),[],2);
            V.seqType = t.seqType;
            V.regType = t.regType;
            V.SN = t.SN;
            VV = addstruct(VV,V);
        end
        keyboard;
        anovaMixed(VV.corIdx,VV.SN,'within',[VV.regType VV.seqType],{'region','seqType'})
    case 'PCM:stats_shape_anova'
        % here do a week x seqType anova per model
        modelType = 'nonlinear_individual';
        sessType = 'transitions';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        if strcmp(modelType,'wSess')
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
            T=rmfield(T,'reg');
        end
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
            t.metric3 = bsxfun(@minus,t.(metric),max(t.(metric),2));
            for m=1:size(t.likelihood,2)
           %    fprintf('%s: ANOVA - model %1.2f\n',regname_cortex{r},m/30);
             %   anovaMixed(t.metric2(:,m),t.SN,'within',[t.sessTr t.seqType],{'session','seqType'},'subset',t.sessTr~=3);
                fprintf('%s: t-test - model %1.2f - transition 1\n',regname_cortex{r},m/30);
                ttestDirect(t.metric2(:,m),[t.seqType t.SN],2,'paired','subset',t.sessTr==1);
                %fprintf('%s: ANOVA - model %1.2f - transition 2\n',regname_cortex{r},m/30);
                %anovaMixed(t.(metric)(:,m),t.SN,'within',[t.seqType],{'seqType'},'subset',t.sessTr==2);
                
            end
            fprintf('All %s\n\n',regname_cortex{r});
        end
    case 'PCM:stats_crossval'
            % here do a week x seqType anova per model
        modelType = 'nonlinear_individual';
        sessType = 'transitions';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        if strcmp(modelType,'wSess')
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
            T=rmfield(T,'reg');
        end
        allSubj = 1:26; 
        for r=reg
            NN=[];
            t = getrow(T,T.regType==r & T.regSide==hemi);
            t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
            t.metric3 = bsxfun(@minus,t.(metric),max(t.(metric),2));
            for s=1:26
                testSubj = allSubj(allSubj~=s)';
                [~,corT1]=max(t.metric2(ismember(t.SN,testSubj)&t.sessTr==1&t.seqType==1,:),[],2);
                [~,corU1]=max(t.metric2(ismember(t.SN,testSubj)&t.sessTr==1&t.seqType==2,:),[],2);
                [~,corT2]=max(t.metric2(ismember(t.SN,testSubj)&t.sessTr==2&t.seqType==1,:),[],2);
                [~,corU2]=max(t.metric2(ismember(t.SN,testSubj)&t.sessTr==2&t.seqType==2,:),[],2);
                % here get the new values
                N.val(1,:) = t.metric2(t.SN==s&t.seqType==1&t.sessTr==1,round(mean(corT1)));
                N.val(2,:) = t.metric2(t.SN==s&t.seqType==1&t.sessTr==1,round(mean(corU1)));
                N.val(3,:) = t.metric2(t.SN==s&t.seqType==1&t.sessTr==2,round(mean(corT2)));
                N.val(4,:) = t.metric2(t.SN==s&t.seqType==1&t.sessTr==2,round(mean(corU2)));
                N.val(5,:) = t.metric2(t.SN==s&t.seqType==1&t.sessTr==3,round(mean(corT2)));
                N.val(6,:) = t.metric2(t.SN==s&t.seqType==1&t.sessTr==3,round(mean(corU2)));
                N.seqType = [1;2;1;2;1;2];
                N.sess = [1;1;2;2;3;3];
                N.sn = [s;s;s;s;s;s];
                NN=addstruct(NN,N);
            end
            fprintf('Region %s:\ns',regname_cortex{r});
            ttestDirect(NN.val,[NN.seqType NN.sn],2,'paired','split',NN.sess);
        end
    case 'PCM:stats_acrSess_seq'
        modelType = 'nonlinear_individual';
        parcelType = 'Brodmann';
        sessType = 'transitions';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, bayesEst_cross or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        if strcmp(modelType,'wSess')
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
        end
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            t.meanLike = mean(t.(metric)(:,2:end),2);
            fprintf('\ndistance logBayes for sequence>session in %s across transitions:1-3\n',regname_cortex{r});
            anovaMixed(t.meanLike,t.SN,'within',[t.sessTr t.seqType],{'session','seqType'});
            fprintf('\ndistance logBayes for sequence>session in %s across transitions:1-2\n',regname_cortex{r});
            anovaMixed(t.meanLike,t.SN,'within',[t.sessTr t.seqType],{'session','seqType'},'subset',t.sessTr<3);
            fprintf('\nt-tests per session: trained vs. untrained:%s\n',regname_cortex{r});
            ttestDirect(t.meanLike,[t.seqType t.SN],2,'paired','split',t.sessTr);
        end
    case 'PCM:behaviourCorr'
        modelType = 'nonlinear_individual';
        parcelType = 'Brodmann';
        sessType = 'transitions';
        reg = [1:3,7];
        hemi = 1;
        metric = 'bayesEst'; % bayesEst, bayesEst_cross or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric','sessType'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s_%s.mat',parcelType,modelType,sessType)));
        B=load(fullfile(behDir,'analyze','testsRH.mat'));
        B.seqType(B.seqType>1) = 2;
        B=getrow(B,B.SN~=4);
        [MT,idx,~] = pivottable([B.testDay B.SN],[B.seqType],B.MT,'nanmedian','subset',B.isError==0 & B.FoSEx==1);
        K.MT = MT(:);
        K.seqType = [ones(size(MT,1),1);ones(size(MT,1),1)*2];
        K.day = [idx(:,1);idx(:,1)];
        K.sn  = [idx(:,2);idx(:,2)];
        behDiff = K.MT(K.seqType==1 & K.day==1) - K.MT(K.seqType==1 & K.day==2);
        behMax = K.MT(K.seqType==1&K.day==1) - K.MT(K.seqType==1&K.day==4);
        for r=reg
            t1 = getrow(T,T.regType==r & T.regSide==hemi & T.SN~=6 & T.sessTr == 1 & T.seqType==1);
            t2 = getrow(T,T.regType==r & T.regSide==hemi & T.SN~=6 & T.sessTr == 1 & T.seqType==2);
            t3 = getrow(T,T.regType==r & T.regSide==hemi & T.SN~=6 & T.sessTr == 2 & T.seqType==1);
            t4 = getrow(T,T.regType==r & T.regSide==hemi & T.SN~=6 & T.sessTr == 3 & T.seqType==1);
            [i1,j1]=max(t1.(metric),[],2);
            [i2,j2]=max(t2.(metric),[],2);
            [i3,j3]=max(t3.(metric),[],2);
            [i4,j4]=max(t4.(metric),[],2);
            [rho,p]=corr(t2.likelihood(:,1)-t1.likelihood(:,1),behMax) 
            for s=1:25 % all subjects
                l_t(s) = t1.(metric)(s,j1(s));
                l_u(s) = t2.(metric)(s,j2(s)); 
                l_tc(s) = t1.(metric)(s,j2(s)); 
            end
            l_t = t1.(metric)(:,j1); % trained
            l_u = t2.(metric)(:,j2); % control
            l_tc = t1.(metric)(:,j2); % trained comparison
        end
        
        
    case 'TESSEL:select'
        % choose all tessels where data present in all subj
        sessN=1;
        hemi=1;
        betaChoice='multi';
        tesselType = 'tesselsWB_362';
        vararginoptions(varargin,{'sessN','betaChoice','hemi','tesselType'});
        
        for ss=sessN       
            indx=[];
            T=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',tesselType,betaChoice,ss)));
            T = getrow(T,ismember(T.regSide,hemi));
            for r=unique(T.region)'
                T1 = getrow(T,T.region==r);
                if ~any(isnan(T1.dist_train)) % make sure all subjects have data there
                    indx = [indx;r];  
                end
            end
            if ss==sessN(1)
                choice=indx;
            else
                choice=unique([indx; choice]);
            end
        end
        varargout={choice};
        
    case 'DIST:plot_mahalanobis'
        roi=1:8;
        hemi=1;
        sessN=1:3;
        parcelType='Brodmann';
        vararginoptions(varargin,{'parcelType','roi','sessN','hemi'});
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)));
        
        T.dist      = [D.dist_train; D.dist_untrain];
        T.sn        = [D.sn; D.sn];
        T.regType   = [D.regType; D.regType];
        T.regSide   = [D.regSide; D.regSide];
        T.sessN     = [D.sessN; D.sessN];
        T.seqType   = [ones(size(D.dist_train));ones(size(D.dist_train))*2];
        T = normData(T,'dist');
        style.use('Seq');
        figure
        for r=roi
            subplot(1,numel(roi),find(roi==r))
            if numel(sessN)==3
                plt.line(T.sessN,ssqrt(T.normdist),'split',T.seqType,'subset',T.regType==r & T.regSide==hemi & T.sessN~=4,'leg',{'trained','control'});
            elseif numel(sessN)==2
                plt.line(T.sessN,ssqrt(T.normdist),'split',T.seqType,'subset',T.regType==r & T.regSide==hemi & T.sessN>2,'leg',{'trained','control'});
            else
                plt.line([T.sessN>3 T.sessN],ssqrt(T.normdist),'split',T.seqType,'subset',T.regType==r & T.regSide==hemi,'leg',{'trained','control'});
            end
            ylabel('Distance'); xlabel('Session'); title(regname_cortex{r});
        end
    case 'DIST:plot_cosine'
        % here plot cosine distances
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        sessN=1:3;
        vararginoptions(varargin,{'parcelType','roi','sessN'});
        D=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI',parcelType)));
        D=normData(D,'corrDist');
        style.use('Seq');
        figure      
        for r=roi
            subplot(1,numel(roi),find(r==roi))
            if numel(sessN)==3
                plt.line(D.sessN,D.normcorrDist,'split',D.seqType,'subset',D.regType==r & D.regSide==hemi & D.sessN~=4,'leg',{'trained','control'});
            else
                plt.line([D.sessN>3 D.sessN],D.normcorrDist,'split',D.seqType,'subset',D.regType==r & D.regSide==hemi,'leg',{'trained','control'});
            end
            ylabel('Distance'); xlabel('Session'); title(regname_cortex{r});
        end
    case 'DIST:stats_mahalanobis'
        roi=1:8;
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI.mat',parcelType)));
        D=getrow(D,D.regSide==hemi);
        T.dist=[D.dist_train; D.dist_untrain];
        T.seqType=[ones(size(D.dist_train));ones(size(D.dist_train))*2];
        T.sn=[D.sn;D.sn];
        T.regType=[D.regType;D.regType];
        T.regSide=[D.regSide;D.regSide];
        T.sessN=[D.sessN;D.sessN];
        
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname_cortex{r});
            anovaMixed(T.dist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r&T.regSide==hemi);
        end
        keyboard;
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname_cortex{r});
                ttestDirect(T.dist,[T.seqType T.sn],2,'paired','subset',T.regType==r&T.regSide==hemi&T.sessN==ss);
                fprintf('trained seq vs. 0\n');
                ttestDirect(T.dist,[T.sn],2,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==1);
                fprintf('untrained seq vs. 0\n');
                ttestDirect(T.dist,[T.sn],2,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==2);
            end
        end
        keyboard;
        for r=roi
            for str=1:3 % session transitions
                fprintf('\n post-hoc t-test on the effect of transition %d-%d trained - %s \n',str,str+1,regname_cortex{r});
                p=[T.dist(T.sessN==str&T.seqType==1&T.regType==r&T.regSide==hemi) T.dist(T.sessN==str+1&T.seqType==1&T.regType==r&T.regSide==hemi)];
                [a,b,c,d]=ttest2(p(:,1),p(:,2))
                %ttestDirect(T.dist,[T.sessN T.sn],1,'paired','subset',T.regType==r & T.seqType==1 & T.regSide==hemi & ismember(T.sessN,[str,str+1]));
                fprintf('\n post-hoc t-test on the effect of transition %d-%d untrained - %s \n',str,str+1,regname_cortex{r});
                %ttestDirect(T.dist,[T.sessN T.sn],1,'paired','subset',T.regType==r & T.seqType==2 & T.regSide==hemi & ismember(T.sessN,[str,str+1]));
                p=[T.dist(T.sessN==str&T.seqType==2&T.regType==r&T.regSide==hemi) T.dist(T.sessN==str+1&T.seqType==2&T.regType==r&T.regSide==hemi)];
                [a,b,c,d]=ttest2(p(:,1),p(:,2))
                
            end
        end
    case 'DIST:stats_cosine'
        roi=1:8;
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        T=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI.mat',parcelType)));

        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname_cortex{r});
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r&T.regSide==hemi);
            fprintf('\n Dist in %s - sess 1-3\n',regname_cortex{r});
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r&T.regSide==hemi &T.sessN~=4);
            fprintf('Only session 3 & 4:\n');
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r&T.regSide==hemi&T.sessN>2);
        end
        keyboard;
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname_cortex{r});
                ttestDirect(T.corrDist,[T.seqType T.sn],2,'paired','subset',T.regType==r&T.sessN==ss&T.regSide==hemi);
                fprintf('trained seq vs. 0\n');
                ttestDirect(T.corrDist,[T.sn],1,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==1);
                fprintf('untrained seq vs. 0\n');
                ttestDirect(T.corrDist,[T.sn],1,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==2);
            end
        end
        keyboard;
        for r=roi
            for str=1:3 % session transitions
                 fprintf('\n post-hoc t-test on the effect of transition %d-%d trained - %s \n',str,str+1,regname_cortex{r});
                 ttestDirect(T.corrDist,[T.sessN T.sn],2,'paired','subset',T.regType==r & T.seqType==1 & T.regSide==hemi & ismember(T.sessN,[str,str+1]));
               % p=[T.corrDist(T.sessN==str&T.seqType==1&T.regType==r&T.regSide==hemi) T.corrDist(T.sessN==str+1&T.seqType==1&T.regType==r&T.regSide==hemi)];
                %[a,b,c,d]=ttest2(p(:,1),p(:,2))
                fprintf('\n post-hoc t-test on the effect of transition %d-%d untrained - %s \n',str,str+1,regname_cortex{r});
                %p=[T.corrDist(T.sessN==str&T.seqType==2&T.regType==r&T.regSide==hemi) T.corrDist(T.sessN==str+1&T.seqType==2&T.regType==r&T.regSide==hemi)];
                %[a,b,c,d]=ttest2(p(:,1),p(:,2))
                ttestDirect(T.corrDist,[T.sessN T.sn],2,'paired','subset',T.regType==r & T.seqType==2 & T.regSide==hemi & ismember(T.sessN,[str,str+1]));
            end
        end
    
    case 'DIST:MDS_plot'
        roi=[1:8];
        sessN=[1:4];
        betaChoice='multi';
        parcelType = 'Brodmann';
        sn=[4:9,11:31];
        
        vararginoptions(varargin,{'roi','sessN','betaChoice','sn','parcelType'});
        
        seqNumb = kron(ones(2,1),[1:6]'); % 1-6
        seqType  = kron([1:2]',ones(6,1)); % 1: trained; 2: untrained
        color={'r','b'};
        
        for ss=sessN
            T=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
            
            for r=roi
                D = getrow(T,T.region==r & ismember(T.SN,sn));
                IPM = mean(D.IPM,1);
                Z = indicatorMatrix('identity_p',1:12);
                Z=bsxfun(@minus,Z,mean(Z,2));
                [Y3,~,~]=rsa_classicalMDS(IPM,'mode','IPM','contrast',Z);
                figure(r)
                subplot(1,4,ss);
                scatterplot3(Y3(:,1),Y3(:,2),Y3(:,3),'split',seqType,'label',seqNumb,'markercolor',color,'markerfill',color);
                
                for i=1:2
                    hold on;
                    indx=find(seqType==i);
                    indx(end+1)=indx(1);
                    line(Y3(indx,1),Y3(indx,2),Y3(indx,3),'color',color{i});
                end;
                plot3(0,0,0,'k+');
                title(sprintf('Session %d %s',ss,regname_cortex{r}),'FontSize',8);
                axis equal;
                
                hold off;
                set(gcf,'PaperPosition',[2 2 4 8],'Color',[1 1 1]);
            end
        end
    case 'DIST:MDS_sess_plot'
        roi=1:8;
        sessN=1:4;
        parcelType = 'Brodmann';
        sn=[5:9,11:31];
        
        vararginoptions(varargin,{'roi','sessN','betaChoice','sn','parcelType'});
        TT = []; NN = [];
        for ss=sessN
            T=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            T.sessN = ones(size(T.SN))*ss;
            TT = addstruct(TT,T);
        end
        
        for s=sn
            for r=roi
                t=getrow(TT,TT.SN==s&TT.region==r);
                if (size(t.betaW{1},2)~=size(t.betaW{2},2)) ||  (size(t.betaW{3},2)~=size(t.betaW{4},2))
                    keyboard;
                end
            end
        end
        
        for r=roi
            for s=sn
                D = getrow(TT,TT.region==r & TT.SN==s);
                if max(sessN)==4
                   beta = [D.betaW{1}(1:96,:);D.betaW{2}(1:96,:);D.betaW{3}(1:96,:);D.betaW{4}(1:96,:)];
                    condV = [kron(ones(8,1),(1:12)');kron(ones(8,1),(13:24)');kron(ones(8,1),(25:36)');kron(ones(8,1),(37:48)')];
                elseif max(sessN)==3
                   beta = [D.betaW{1}(1:96,:);D.betaW{2}(1:96,:);D.betaW{3}(1:96,:)];
                    condV = [kron(ones(8,1),(1:12)');kron(ones(8,1),(13:24)');kron(ones(8,1),(25:36)')];
                else
                    error('wrong number of sessions!');
                end
                partV = repmat(kron((1:8)',ones(12,1)),max(sessN),1);
              %  dist = rsa.distanceLDC(beta,partV,condV);
                G = pcm_estGCrossval(beta,partV,condV);
                dist=corr_crossval(G,'reg','minvalue');      % use cosine distance       
                N.dist = dist;
                N.sn = s;
                N.roi = r;
                NN = addstruct(NN,N);
            end
            RDM = mean(NN.dist,1);
            sessV = kron((sessN)',ones(12,1));
            Z1 = indicatorMatrix('identity_p',kron(ones(max(sessN),1),[ones(6,1);ones(6,1)*2]));
            Z2 = indicatorMatrix('identity_p',sessV);
            Z4 = indicatorMatrix('identity_p',kron(ones(max(sessN),1),[1:12]'));
            [Y3,~,~]=rsa_classicalMDS(RDM,'mode','RDM','contrast',Z4);
            figure
            scatterplot3(Y3(:,1),Y3(:,2),Y3(:,3),'split',[kron(ones(max(sessN),1),[ones(6,1);ones(6,1)*2])],'markercolor',{'r','b'},'markerfill',{'r','b'});
            colS = {'r','b'};
            lineS = {'-','--','-.',':'};
            for i=1:2
                for j=1:max(sessN)
                    hold on;
                    indx=Z1(:,i)==1 & Z2(:,j)==1;
                    line(Y3(indx,1),Y3(indx,2),Y3(indx,3),'LineStyle',lineS{j},'color',colS{i});
                end
            end;
            
            %   plot3(0,0,0,'k+');
            title(regname_cortex{r});
            axis equal;
        end
        keyboard;
    case 'DIST:plot_seqType_mahalanobis'
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'parcelType','roi'});
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)));
        
        D = normData(D,'dist_cross');
        style.use('gray');
        for r=roi
            figure
            plt.line([D.sessN>3 D.sessN],ssqrt(D.normdist_cross),'subset',D.regType==r & D.regSide==hemi,'leg',{'trained','control'});
            ylabel('Distance'); xlabel('Session'); title(regname_cortex{r});
        end
    case 'DIST:plot_seqType_cosine'
        roi=1:8;
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','parcelType','hemi'});
        
        D=load(fullfile(distPscDir,sprintf('corrDist_seqType_%s_ROI.mat',parcelType)));
        style.use('gray');
        for r=roi
            figure
            plt.line([D.sessN>3 D.sessN],D.corr_seqType,'subset',D.regType==r&D.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname_cortex{r}));
            ylabel('Correlation distance betweeen seqTypes');
            xlabel('Session');
        end
    case 'DIST:plot_seqType_mahalanobis_allReg'
        % plot all regions in the same plot
        roi = [1,2,3,7]; % S1, M1, PMd, SPLa
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'parcelType','roi'});
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)));
        style.use('4Reg');
        D = normData(D,'dist_cross');
        D = getrow(D,ismember(D.regType,roi) & D.regSide==1); % only contralateral
        figure
        plt.line([D.sessN>3 D.sessN],ssqrt(D.normdist_cross),'split',D.roi ,'leg',{'S1','M1','PMd','SPLa'});
        xlabel('Scan'); ylabel('Sequence type distance'); title('mahalanobis');
    case 'DIST:plot_seqType_cosine_allReg'
        % plot all regions in the same plot
        roi = [1,2,3,7]; % S1, M1, PMd, SPLa
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'parcelType','roi'});
        D=load(fullfile(distPscDir,sprintf('corrDist_seqType_%s_ROI.mat',parcelType)));
        style.use('4Reg');
        D = normData(D,'corr_seqType');
        D = getrow(D,ismember(D.regType,roi) & D.regSide==1); % only contralateral
        figure
        plt.line([D.sessN>3 D.sessN],ssqrt(D.normcorr_seqType),'split',D.roi ,'leg',{'S1','M1','PMd','SPLa'});
        xlabel('Scan'); ylabel('Sequence type distance'); title('cosine');
        figure
        plt.line(D.sessN,D.normcorr_seqType,'split',D.roi,'subset',D.sessN>2,'leg',{'S1','M1','PMd','SPLa'});
        xlabel('Scan'); ylabel('Sequence type distance'); title('cosine');
    case 'DIST:stats_seqType_mahalanobis'
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'parcelType','roi'});
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)));
        
        for r=roi
            t = getrow(D,D.regType==r & D.regSide==hemi);
            for str = 1:3 % sessions
                fprintf('%s - session transition: %d-%d\n',regname_cortex{r},str,str+1)
                ttestDirect(t.dist_cross,[t.sessN,t.sn],2,'paired','subset',ismember(t.sessN,[str,str+1]));
            end
        end
    case 'DIST:stats_seqType_cosine'
        roi=1:8;
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','parcelType','hemi'});
        
        D=load(fullfile(distPscDir,sprintf('corrDist_seqType_%s_ROI.mat',parcelType)));
        
        for r=roi
            t = getrow(D,D.regType==r & D.regSide==hemi);
            fprintf('%s sessions 1-3\n:',regname_cortex{r});
            anovaMixed(t.corr_seqType,t.sn,'within',t.sessN,{'session'},'subset',t.sessN~=4);
            fprintf('%s sessions 3-4\n:',regname_cortex{r});
            anovaMixed(t.corr_seqType,t.sn,'within',t.sessN,{'session'},'subset',t.sessN>2);
            for str = 1:3 % sessions
                fprintf('%s - session transition: %d-%d\n',regname_cortex{r},str,str+1)
                ttestDirect(t.corr_seqType,[t.sessN,t.sn],2,'paired','subset',ismember(t.sessN,[str,str+1]));
            end
        end
    case 'DIST:striatum_shift'
        % here quantify if there is a shift within the striatum across
        % sessions
        sessN=1:4;
        parcelType = 'BG-striatum';
        sn=[5:9,11:31];
        
        vararginoptions(varargin,{'roi','sessN','betaChoice','sn','parcelType'});
        TT = []; NN = [];
        condVec = repmat([ones(6,1);ones(6,1)*2],8,1);
        Z = pcm_indicatorMatrix('identity',condVec);
        for ss=sessN
            T=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            T.sessN = ones(size(T.SN))*ss;
            TT = addstruct(TT,T);
        end
        sessTr = indicatorMatrix('allpairs',sessN);
        for s=sn
            for str=1:size(sessTr,1)
                s1 = find(sessTr(str,:)==1);
                s2 = find(sessTr(str,:)==-1);
                t1 = getrow(TT,TT.SN==s & TT.sessN==s1);
                t2 = getrow(TT,TT.SN==s & TT.sessN==s2);
                Data1 = [t1.betaW{1} t1.betaW{2}];
                Data2 = [t2.betaW{1} t2.betaW{2}];
                dataM1 = pinv(Z)*Data1(1:length(Z),:);           % Estimate mean activities
                dataM2 = pinv(Z)*Data2(1:length(Z),:);           % Estimate mean activities
                N.sess1 = [s1;s1];
                N.sess2 = [s2;s2];
                N.sessTr = [str;str];
                N.seqType = [1;2];
                N.corrMean(1,:) = corr((dataM1(1,:))',(dataM2(1,:))');
                N.corrMean(2,:) = corr((dataM1(2,:))',(dataM2(2,:))');
                N.cosMean(1,:) = pdist([(dataM1(1,:));(dataM2(1,:))],'cos');
                N.cosMean(2,:) = pdist([(dataM1(2,:));(dataM2(2,:))],'cos');
                N.sn = [s;s];
                NN=addstruct(NN,N);
                fprintf('Done subject %s: transitions %d.\n',subj_name{s},str);
            end
        end
        keyboard;
        % significance test
        ttestDirect(NN.corrMean,[NN.seqType,NN.sn],2,'paired','split',NN.sessTr)
        ttestDirect(NN.cosMean,[NN.seqType,NN.sn],2,'paired','split',NN.sessTr)

    case 'DIST:seqType_chance_sess1' % not currently used
        % calculate what the chance distance is between all splits into 6
        % and 6 sequences
        roi = [1,2,3,7];
        parcelType = 'Brodmann';
        betaChoice = 'multi'; 
        vararginoptions(varargin,{'parcelType','betaChoice'});
        
        allSeq = 1:12;
        comb = nchoosek(allSeq,6);
        comb = comb(2:end-1,:);
        T = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess1.mat',parcelType,betaChoice))); % loads region data (D)
        for r = roi
            t = getrow(T,T.regType==r & T.regSide==1 & T.SN~=4);
            RDM = rsa_squareRDM(nanmean(t.RDM,1));
            cross_null = zeros(size(comb,1),1);
            for i = 1:size(comb,1); % all combinations for null distribution
                c1 = comb(i,:);
                idx = ~ismember(1:12,c1);
                c2 = allSeq(idx);
                cross_null(i) = nanmean(nanmean(RDM(c1,c2)));
            end
            figure
            histogram(cross_null); hold on;  
            %drawline(ssqrt(nanmean(t.dist_cross,1)),'dir','vert'); title(regname_cortex{r});
            drawline(nanmean(nanmean(RDM(1:6,7:12))),'dir','vert'); title(regname_cortex{r});
        end
    case 'DIST:modularity_mahalanobis'
        % determine if across-seqType distance always higher than within
        sessN = 1:4;
        roi = 1:8;
        hemi = 1;
        parcelType = 'Brodmann';
        betaChoice = 'multi';
        DD = [];
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
            T = getrow(T,T.SN~=4);
            D.RDM = T.RDM;
            D.regType = T.regType;
            D.regSide = T.regSide;
            D.session = ones(size(T.regSide))*ss;
            D.SN = T.SN;
            DD = addstruct(DD,D);
        end
        for r=roi
            figure
            t = getrow(DD,DD.regType==r&DD.regSide==hemi);
            for ss=1:4
                subplot(1,4,ss)
                t1 = getrow(t,t.session==ss);
                rdm = squareform(nanmean(t1.RDM,1));
%                 plot(rsa_vectorizeRDM(rdm(7:12,7:12)),rsa_vectorizeRDM(rdm(1:6,7:12)),'.','MarkerSize',15);
%                 hold on;
%                 plot(rsa_vectorizeRDM(rdm(1:6,1:6)),rsa_vectorizeRDM(rdm(1:6,7:12)),'r.','MarkerSize',15);
%                 % adjust the axes
%                 g = gca;
%                 m = max(g.XLim(2),g.YLim(2));
%                 xlim([0,m]);ylim([0,m]);
%                 refline(1,0);
                histogram(rsa_vectorizeRDM(rdm(7:12,7:12))',5); hold on;
                histogram(rsa_vectorizeRDM(rdm(1:6,1:6))',5); hold on;
                cross = rsa_vectorizeIPMfull(rdm(1:6,7:12)'); % because there are more, randomly subsample 15
                order = randperm(size(cross,2));
                histogram(cross(order(1:15)),5);
                title(sprintf('%s - session %d',regname_cortex{r},ss));
            end
        end
        
    case 'SPEED:psc'
        % here plot psc in sessions 3 and 4
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';      
        vararginoptions(varargin,{'roi','sessN','hemi','parcelType'});     
        T=load(fullfile(distPscDir,sprintf('psc_%s_ROI.mat',parcelType)));   
        T = getrow(T,T.sessN>2);
        switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        style.use('Seq');
        figure
        indx=1;
        T=normData(T,'psc');
        for r=roi
            subplot(1,numel(roi),indx)
            plt.line(T.sessN,T.normpsc,'split',T.seqType,'subset',T.regType==r & T.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north');
            drawline(0,'dir','horz');
            title(sprintf('%s',regLab{r}));
            if r==1
                ylabel('Percent signal change');
                xlabel('Session');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'SPEED:seqType'
        roi=1:8;
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','parcelType','hemi'});
        
        D=load(fullfile(distPscDir,sprintf('corrDist_seqType_%s_ROI.mat',parcelType)));
        D = getrow(D,D.sessN>2);
        style.use('gray');
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.bar(D.sessN,D.corr_seqType,'subset',D.regType==r&D.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north');
            drawline(0,'dir','horz');
            title(sprintf('%s',regname_cortex{r}));
            ylabel('Correlation distance betweeen seqTypes');
            xlabel('Session');
            indx=indx+1;
        end 
        figure
        style.use('4Reg');
        plt.bar(D.regType,D.corr_seqType,'subset',ismember(D.regType,roi)&D.regSide==hemi,'split',D.sessN,'leg',{'trained','untrained'},'leglocation','north');
    case 'SPEED:cos_dist'
        % here plot cosine distance in sessions 3 and 4
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'parcelType','roi','sessN'});
        D=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI',parcelType)));
        D=getrow(D,D.sessN>2);
        D=normData(D,'corrDist');
        style.use('Seq');
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.line(D.sessN,D.normcorrDist,'split',D.seqType,'subset',D.regType==r & D.regSide==hemi,'leg',{'trained','control'});
            ylabel('Distance'); xlabel('Session'); title(regname_cortex{r});
            indx=indx+1;
        end
    case 'SPEED:ratio'
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
                fprintf('%s - roi %d/%d.\n',hemName{h},find(r==roi),length(roi));
            end
        end
        save(fullfile(distPscDir,sprintf('SPEED_sess3-4_dist_psc_%s',parcelType)),'-struct','NN');
    case 'SPEED:plot_ratio'
        % plot in bargraphs (appropriate for Brodmann) % for 162 project
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','parcelType'});
        
        T = load(fullfile(distPscDir,sprintf('SPEED_sess3-4_dist_psc_%s',parcelType)));
        for r=roi
            figure
            style.use('Sess');
            subplot(221)
            plt.bar(T.seqType,[T.psc3 T.psc4],'subset',T.regType==r&T.regSide==hemi,'leg',{'sess-3','sess-4'});
            xlabel('seqType'); ylabel('psc');
            title(sprintf('%s - percent signal change',regname_cortex{r}));
            
            subplot(223)
            plt.bar(T.seqType,[T.dist3 T.dist4],'subset',T.regType==r&T.regSide==hemi,'leg',{'sess-3','sess-4'});
            drawline(mean(T.predDist(T.regType==r&T.regSide==hemi&T.seqType==1)),'dir','horz','lim',[1.7 2.7],'linestyle','--');
            drawline(mean(T.predDist(T.regType==r&T.regSide==hemi&T.seqType==2)),'dir','horz','lim',[4.9 5.9],'linestyle','--');
            xlabel('seqType'); ylabel('dist');
            title(sprintf('%s - distance',regname_cortex{r}));
            
            style.use('Seq');
            subplot(222)
            plt.scatter(T.psc3,T.psc4,'split',T.seqType,'subset',T.regType==r & T.regSide==hemi);
            xlabel('psc(sess3)'); ylabel('psc(sess4)'); title('psc - individual subjects');
            axis equal;hline = refline(1,0); hline.Color = 'k';    
            subplot(224)
            plt.scatter(T.dist3,T.dist4,'split',T.seqType,'subset',T.regType==r & T.regSide==hemi);
            xlabel('dist(sess3)'); ylabel('dist(sess4)'); title('dist - individual subjects');
            axis equal;hline = refline(1,0); hline.Color = 'k';
        end
    case 'SPEED:stats_ratio'
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','parcelType'});
        T = load(fullfile(distPscDir,sprintf('SPEED_sess3-4_dist_psc_%s',parcelType)));        
        for r=roi
            t = getrow(T,T.regType==r & T.regSide==hemi);
            fprintf('\n\n%s t-test distance vs. predicted - trained\n',regname_cortex{r});
            ttestDirect([t.predDist; t.dist4],[[ones(size(t.predDist));ones(size(t.predDist))*2] [t.sn;t.sn]],1,'paired','subset',[t.seqType;t.seqType]==1);
            fprintf('%s t-test distance vs. predicted - untrained\n',regname_cortex{r});
            ttestDirect([t.predDist; t.dist4],[[ones(size(t.predDist));ones(size(t.predDist))*2] [t.sn;t.sn]],1,'paired','subset',[t.seqType;t.seqType]==2);
            fprintf('%s t-test trained vs. untrained dist ratio\n',regname_cortex{r});
            ttestDirect(t.dist_ratio,[t.seqType t.sn],1,'paired');
            fprintf('%s t-test trained vs. untrained psc ratio\n',regname_cortex{r});
            ttestDirect(t.psc_ratio,[t.seqType t.sn],1,'paired');
            fprintf('%s t-test trained vs. untrained pred vs. ratio\n',regname_cortex{r});
            ttestDirect(t.dist4-t.predDist,[t.seqType t.sn],1,'paired');
            fprintf('%s t-test trained vs. untrained distance 3\n',regname_cortex{r});
            ttestDirect(t.dist3,[t.seqType t.sn],1,'paired');
            fprintf('%s t-test trained vs. untrained distance 4\n',regname_cortex{r});
            ttestDirect(t.dist4,[t.seqType t.sn],1,'paired');
        end
    case 'SPEED:distance_behavior'
        % here determine if there's higher distances for trained >
        % untrained based on behavioral performance
        sn = [5:9,11:31];
        D = load(fullfile(behDir,'analyze','alldata.mat'));
        D = getrow(D,D.isScan==1 & D.isMetronome==2 & ismember(D.SN,sn));
        TT = [];
        for s=sn
            for st=1:2
                d = getrow(D,D.SN==s & D.seqType==st);     
                if st==2
                    d.seqNumb = d.seqNumb-6;
                end
                dist = rsa.distanceLDC(d.MT,d.BN,d.seqNumb);
                T.distLDC = mean(dist);
                T.seqType = st;
                T.sn = s;
                TT = addstruct(TT,T);
            end
        end
        ttestDirect(TT.distLDC,[TT.seqType,TT.sn],2,'paired');
        ttestDirect(TT.distLDC,[TT.sn],2,'onesample','subset',TT.seqType==1);
        ttestDirect(TT.distLDC,[TT.sn],2,'onesample','subset',TT.seqType==2);
        figure
        plt.box(TT.seqType,TT.distLDC,'plotall',0); ylabel('Distance'); xlabel('Sequence');
    case 'Fig1:psc_roi'
        roi = [2,1,3,7];
        hemi=1;
        parcelType = 'Brodmann';
        T=load(fullfile(distPscDir,sprintf('psc_%s_ROI.mat',parcelType)));
        
        figure
        style.use('Seq');
        for r=1:length(roi)
            subplot(1,numel(roi),r)
            plt.line(T.sessN,T.psc,'split',T.seqType,'subset',T.regType==r & T.sessN<4 & T.regSide==hemi,'leg',{'trained','untrained'},'leglocation','northeast');
            set(gca,'XTickLabel',{1,2,5})
            ylim([0 1.8]);
            drawline(0,'dir','horz');
            title(sprintf('%s',regname_cortex{roi(r)}));
            if r==1
                ylabel('Percent signal change');
                xlabel('Week');
            else
                ylabel('');
            end
        end
    case 'Fig1:behaviour'
        
        D = load(fullfile(behDir,'analyze','alldata.mat'));
        B = getrow(D,ismember(D.SN,sn)&D.blockType==6);
        [movTime,subjDay,~] = pivottable([B.day B.SN],[B.seqType],B.MT,'median','subset',B.FoSEx==1); % only 1st execution - numbers seen
        S = getrow(D,ismember(D.SN,sn)&D.ScanSess>0 & D.seqType<3);
        S.day(S.ScanSess==1)=2;
        S.day(S.ScanSess==2)=7;
        S.day(S.ScanSess==3)=18;
        S.day(S.ScanSess==4)=19;
        
        SS = load(fullfile(anaDir,'scanner_summary'));
        
%         [movTimeSc,idx,~] = pivottable([S.day S.SN],[S.seqType],S.MT,'median','subset',S.isError==0 & S.FoSEx==1 & S.SN~=6);
%         SS.MT       = movTimeSc(:);
%         SS.seqType  = [ones(size(movTimeSc,1),1);ones(size(movTimeSc,1),1)*2];
%         SS.day      = [idx(:,1);idx(:,1)];
%        SS.sn       = [idx(:,2); idx(:,2)];
        figure
        style.use('SeqShade');
        subplot(121)
        plt.box(SS.sessN,SS.MT,'split',SS.seqType,'plotall',0,'leg',{'trained','control'});
        % first scanner behaviour
        subplot(122)
        plt.line([subjDay(:,1)>6 subjDay(:,1)],movTime);
        ylabel('Movement time (msec)'); xlabel('Days');
        
        B=load(fullfile(behDir,'analyze','testsRH.mat'));
        B.seqType(B.seqType>1) = 2;
        B=getrow(B,B.SN~=4);
        B.day(B.day==19)=20;
        [MT,idx,~] = pivottable([B.day B.SN],[B.seqType],B.MT,'nanmedian','subset',B.isError==0);
     %   [MT,idx,~] = pivottable([B.testDay B.SN],[B.seqType],B.MT,'nanmedian','subset',B.isError==0 & B.FoSEx==1);
        K.MT = MT(:);
        K.seqType = [ones(size(MT,1),1);ones(size(MT,1),1)*2];
        K.day = [idx(:,1);idx(:,1)];
        K.sn  = [idx(:,2);idx(:,2)];
        style.use('Seq');
        figure
        subplot(131)
         plt.line(K.day,K.MT,'split',K.seqType); title('behavior');
       % plt.box(K.day,K.MT,'split',K.seqType); title('behavior');
        subplot(132)
       % plt.box(SS.sessN,SS.MT,'split',SS.seqType,'plotall',0,'leg',{'trained','control'},'subset',SS.FoSEx==1);
        plt.line(SS.sessN,SS.MT,'split',SS.seqType,'leg',{'trained','control'},'subset',SS.FoSEx==1);
        subplot(133)
        plt.line([subjDay(:,1)>6 subjDay(:,1)],movTime);
        ylabel('Movement time (msec)'); xlabel('Days');
        figure
        %plt.dot(K.day,K.MT,'split',K.seqType,'leg',{'trained','random'});
        plot(K.day(K.seqType==1)*2,K.MT(K.seqType==1),'ro','MarkerFaceColor','r');
        hold on;
        plot(K.day(K.seqType==2)*2+1,K.MT(K.seqType==2),'bv','MarkerFaceColor','b');
        subplot(122)
        plt.line(B.testDay,B.MT,'split',B.seqType,'leg',{'trained','random'},'plotfcn','median');
        xlabel('Test day'); ylabel('Movement time');
        plt.match('y');
        
        
    case 'STATS:behaviour'
        S = load(fullfile(anaDir,'scanner_summary'));
        anovaMixed(S.MT,S.sn,'within',[S.sessN S.seqType],{'session','seqType'},'subset',S.sessN~=4);
        ttestDirect(S.MT,[S.seqType S.sn],2,'paired','split',S.sessN);
        ttestDirect(S.MT,[S.sessN S.sn],2,'paired','subset',S.sessN<2|S.sessN>3,'split',S.seqType);
        ttestDirect(S.MT,[S.sessN S.sn],2,'paired','subset',S.sessN<2|S.sessN>3);
        keyboard;
        B=load(fullfile(behDir,'analyze','testsRH.mat'));
        B.seqType(B.seqType>1) = 2;
        [MT,idx,~] = pivottable([B.testDay B.SN],[B.seqType],B.MT,'nanmedian','subset',B.isError==0 & B.FoSEx==1);
        K.MT = MT(:);
        K.seqType = [ones(size(MT,1),1);ones(size(MT,1),1)*2];
        K.day = [idx(:,1);idx(:,1)];
        K.sn  = [idx(:,2);idx(:,2)];
        anovaMixed(K.MT,K.sn,'within',[K.day K.seqType],{'testDay','seqType'});
        ttestDirect(K.MT,[K.seqType K.sn],2,'paired','split',K.day);
        for str=1:3
            ttestDirect(K.MT,[K.day K.sn],2,'paired','split',K.seqType,'subset',ismember(K.day,[str,str+1]));
        end
        
        
    otherwise
        fprintf('No such case!\n');
end
end
%% Local functions
function M = pcm_toyModel_wSess
% function M = pcm_toyModel_wSess
% constructs pcm models considering 'toy models' of 0, perfect,
% flexible correlation between sessions
% here session effect is explicitly modelled
    
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
        M{3}.Ac(:,3:8,2+i)   = [A;zeros(6)];      % Unique sess1 sequence patterns
        M{3}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique sess2 sequence pattterns
    end;

    % Model 4: Session + sequence specific + FLEXIBLE correlation in session
    M{4}.type         = 'feature';
    M{4}.numGparams   = 20;
    M{4}.name         = 'Session+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling - one parameter per session
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
    
    % Model 6: Noise ceiling
    M{6}.type       = 'freedirect';
    M{6}.numGparams = 0;
    M{6}.theta0     = [];
    M{6}.name       = 'noice_ceiling';
    
end
function M = pcm_toyModel_noSess 
% function M = pcm_toyModel_noSess
% constructs pcm models considering 'toy models' of 0, perfect,
% flexible correlation between sessions
% here session effect is NOT modelled
    
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

    % Model 3: Session + sequence specific + FLEXIBLE correlation in session
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
       % M{3}.Ac(:,1:6,6+i)  = [zeros(6);A];     % Correlation sess1-sess2
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
    
    % Model 5: Noise ceiling
    M{5}.type       = 'freedirect';
    M{5}.numGparams = 0;
    M{5}.theta0     = [];
    M{5}.name       = 'noice_ceiling';
    
end

function M = pcm_withinSess
% function M = pcm_withinSess
% constructs models only for one session at a time
% essentially estimating distances 

    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:6,1)  = zeros(6);
    
    % Model 2 - sequence patterns - all given the same parameter
    M{2}.type       = 'feature';
    M{2}.name       = 'sequence';
    M{2}.numGparams = 1;
    M{2}.Ac         = eye(6);

    % Model 3 - sequence patterns - specific
    M{3}.type       = 'feature';
    M{3}.name       = 'individ_sequence';
    M{3}.numGparams = 6;
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,1:6,i)  = A;     % Unique sequence patterns
    end;
    
    % Model 4 - noise ceiling
    M{4}.type       = 'freedirect';
    M{4}.numGparams = 0;
    M{4}.theta0     = [];
    M{4}.name       = 'noice_ceiling';
    % --------------------------------------

end

function M = pcm_corrModel_noSess(corrLim,nModel)
% specific correlation models with no explicit modelling of session
corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
M = cell(1,length(corrS));
% Model 1: Null model
M{1}.type           = 'feature';
M{1}.numGparams     = 1;
M{1}.name           = 'null';
M{1}.Ac(:,1:12,1)   = zeros(12);

% Correlation models
for c=1:length(corrS)
    [th2,th3] = determine_thetaCorr(corrS(c)); % determine correct thetas
    % Build the models
    M{c+1}.type = 'feature';
    M{c+1}.numGparams = 12;
    M{c+1}.name = sprintf('SpecCorr_%1.3f',corrS(c));
    % for sequence-specific modelling - one parameter per sequence
     for i=1:6
         A=zeros(6);
         A(i,i)=1;
         M{c+1}.Ac(:,1:6,i)       = [A;zeros(6)];      % Unique sess1 sequence patterns
         M{c+1}.Ac(:,7:12,6+i)    = [zeros(6);A*th2];  % Unique sess2 sequence pattterns
         M{c+1}.Ac(:,1:6,6+i)     = [zeros(6);A*th3];  % Correlation sess1-sess2
     end;    
end
end
function M = pcm_corrModel_wSess(corrLim,nModel)
% specific correlation models with additional modelling of session
corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
% Model 1: Null model
M{1}.type       = 'feature';
M{1}.numGparams = 1;
M{1}.name       = 'null';
M{1}.Ac(:,1:12,1)  = zeros(12);

% Model 2: First vs. second session
M{2}.type       = 'feature';
M{2}.numGparams = 2;
M{2}.name       = 'Session';
M{2}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
M{2}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];

% Other models: specific correlation
for c=1:length(corrS)
    [th2,th3]       = determine_thetaCorr(corrS(c));
    % Build the models
    M{c+2}.type         = 'feature';
    M{c+2}.numGparams   = 14;
    M{c+2}.name         = sprintf('SpecCorr_%1.3f',corrS(c));
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
function M = pcm_corrModel_meanPattern(corrLim,nModel)
% specific correlation models for mean pattern
corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
M = cell(1,length(corrS)+1);
% Model 1: Null model
M{1}.type           = 'feature';
M{1}.numGparams     = 1;
M{1}.name           = 'null';
M{1}.Ac(:,1:2,1)    = zeros(12,2);

% Correlation models
for c=1:length(corrS)
    [th2,th3] = determine_thetaCorr(corrS(c)); % determine the correct thetas
    % Build models
    M{c+1}.type         = 'feature';
    M{c+1}.numGparams   = 2;
    M{c+1}.name         = sprintf('SpecCorr_%1.3f',corrS(c));
    M{c+1}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];        % Unique mean sess1 pattern
    M{c+1}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)*th2];    % Unique mean sess2 pattern
    M{c+1}.Ac(:,1,2)    = [zeros(6,1);ones(6,1)*th3];    % Shared mean pattern
end
end

function [th2,th3]=determine_thetaCorr(r)
% determine thetas such that correlation is of specific value
if r==1
    th2=0;
    th3=1;
else
    th2=1;
    th3=sqrt((r^2)/(1-r^2)); % equation for determining thetas
end
end

function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Run group fit in 2 stages - first NR, finish with minimize (any further improvement)
     [T,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
  %   fprintf('Group fit with NR algorithm done.\n');
  %   [T,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize','theta0',theta_hat);
  %   fprintf('Group fit with minimize algorithm done.\n');
    % Run crossvalidated group fit
    [Tcross,~]          = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);
    T.cross_likelihood  = Tcross.likelihood;
    T.bayesEst          = bsxfun(@minus,T.likelihood,T.likelihood(:,1));
    T.bayesEst_cross    = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
end
function T = pcm_fitModels_conserv(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Run group fit in 2 stages - first minimize, then NR (get the correct
    % thetas)
     [~,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize');
     fprintf('Group fit with minimize algorithm done.\n');
     [~,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR','theta0',theta_hat);
     fprintf('Group fit with NR algorithm done.\n');
     [T,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize','theta0',theta_hat);
     fprintf('Group fit with minimize algorithm done.\n');
    % Run crossvalidated group fit
    [Tcross,theta_hat]       = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);
    %[Tcross,~]          = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm','minimize');
    %fprintf('Crossvalidated fit with minimize algorithm done.\n');
    T.cross_likelihood  = Tcross.likelihood;
    T.bayesEst          = bsxfun(@minus,T.likelihood,T.likelihood(:,1));
    T.bayesEst_cross    = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
end
function T = pcm_fitModels_corr(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    nModel = size(M,2);
    % Run group fit for each corr model separately - fit the estimated
    % theta of previous as starting
    for i=1:nModel
        if i<4
            [t,theta1,~,~]     = pcm_fitModelGroup(Data,M{i},partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
        else
            [t,theta1,~,~]     = pcm_fitModelGroup(Data,M{i},partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm,'theta0',theta_hat(i-1));
        end
            T.iterations(:,i)  = t.iterations;
     T.time(:,i)        = t.time;
     T.noise(:,i)       = t.noise;
     T.scale(:,i)       = t.scale;
     T.run(:,i)         = t.run;
     T.likelihood(:,i)  = t.likelihood;
     theta_hat(i)       = theta1;     
    end
    T.SN = t.SN;
    
     %   fprintf('Group fit with NR algorithm done.\n');
  %   [T,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize','theta0',theta_hat);
  %   fprintf('Group fit with minimize algorithm done.\n');
    % Run crossvalidated group fit
    [Tcross,~]          = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);
    T.cross_likelihood  = Tcross.likelihood;
    T.bayesEst          = bsxfun(@minus,T.likelihood,T.likelihood(:,1));
    T.bayesEst_cross    = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
end
function C = pcm_correlation(Data,partVec,condVec,M,runEffect,M_type)
%function C = pcm_correlation(Data,partVec,condVec,M,runEffect,M_type)
% estimates the correlation between sessions
% naive, crossvalidated, through PCM model
sn=1:size(Data,2);
% --------------------------------------
% 1. Empirical correlation
for p=sn
    Z               = pcm_indicatorMatrix('identity',condVec{1});
    b               = pinv(Z)*Data{p};           % Estimate mean activities
    b(1:6,:)        = bsxfun(@minus,b(1:6,:) ,mean(b(1:6,:))); % Subtract mean per condition - first exe
    b(7:12,:)       = bsxfun(@minus,b(7:12,:),mean(b(7:12,:))); % second exe
    G               = cov(b');
    C.r_naive(p,1)  = calcCorr(G);
end;
% --------------------------------------
% 2. Crossvalidated correlation
 condSeqTypeVec = [ones(size(condVec{1},1)/2,1);ones(size(condVec{1},1)/2,1)*2];
for p=sn
    % Subtract mean for each condition and run
    X                       = pcm_indicatorMatrix('identity',partVec{1}*2+(condVec{1}>6)-1);
    R                       = eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv_all                 = pcm_estGCrossval(Data{p},partVec{1},condVec{1});      % include all the data
    Gcv_seq                 = pcm_estGCrossval(R*Data{p},partVec{1},condVec{1});    % only seq-specific pattern (subtract the mean)
    Gcv_seqType             = pcm_estGCrossval(Data{p},partVec{1},condSeqTypeVec);  % only the mean pattern
    C.r_cross_all(p,1)      = calcCorr(pcm_makePD(Gcv_all));
    C.r_cross_seq(p,1)      = calcCorr(pcm_makePD(Gcv_seq));
    C.r_cross_seqType(p,1)  = calcCorr(pcm_makePD(Gcv_seqType));
end;
% --------------------------------------
% 3. Fit model 2  and infer correlations from the parameters
[~,theta,G_pred] = pcm_fitModelIndivid(Data,M,partVec{1},condVec{1},'runEffect',runEffect);
% set up sequence type model
M_st.Ac(:,:,1) = [1 0; 0 0]; M_st.Ac(:,:,2) = [0 0; 0 1]; M_st.Ac(:,:,3) = [0 1; 0 0];
M_st.type = 'feature';
M_st.numGparams = 3;
[~,theta_st,G_st] = pcm_fitModelIndivid(Data,M_st,partVec{1},condSeqTypeVec,'runEffect',runEffect);
% note: last two thetas are estimates of noise and run effects
% Get the correlations
switch M_type
    case 'wSess'
        var1       = (theta{1}(3:8,:).^2)';
        var2       = (theta{1}(9:14,:).^2+theta{1}(15:20,:).^2)';
        cov12      = (theta{1}(3:8,:).*theta{1}(15:20,:))';
        C.r_model =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
    case 'noSess'
        var1       = (theta{1}(1:6,:).^2)';
        var2       = (theta{1}(7:12,:).^2+theta{1}(13:18,:).^2)';
        cov12      = (theta{1}(1:6,:).*theta{1}(13:18,:))';
        C.r_model =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
end
% here get the sequence type model correlation
% var1 = theta_st{1}(1,:).^2;
% var2 = theta_st{1}(2,:).^2;
% cov12 = theta_st{1}(3,:).*theta_st{1}(1);
% C.r_model_seqType = (cov12./(sqrt(var1.*var2)))';

for p=sn
    g = G_st{:}(:,:,p);
    C.r_model_seqType(p,1) = calcCorr(g);
end
% --------------------------------------
end
function r = calcCorr(G)
% function r = calcCorr(G)
% calculates the correlation between first and second half of conditions
% given G as input
    nCond = size(G,1)/2;
    d0 = diag(G);
    v1 = d0(1:nCond)';          % Variances 1st session conditions
    v2 = d0(nCond+1:end)';      % Variances 2nd session conditions
    cv=diag(G,nCond);           % Covariance
    r = mean(cv)/sqrt(mean(v1)*mean(v2));

end
