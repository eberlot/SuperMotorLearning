function varargout = sml_plasticity(what,varargin)
% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
caretDir        =[baseDir '/surfaceCaret'];     
behDir          =[baseDir '/behavioral_data'];
stabDir         =[baseDir '/stability'];
distPscDir      =[baseDir '/dist_psc_stats'];
regDir          =[baseDir '/RegionOfInterest/']; 
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session
codeDir         ='/Users/eberlot/Documents/MATLAB/projects/SuperMotorLearning';

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
    case 'PSC:surface'
     % create surface maps of percent signal change 
     % trained and untrained sequences - per subject
        smooth = 0;   
        vararginoptions(varargin,{'sn','sessN','smooth'});

        name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','TrainSeq','UntrainSeq'};
        for ss=sessN
            fileList = cell(length(name),1);% column_name = fileList; images = fileList;
            for n = 1:length(name)
                fileList{n}=fullfile(['psc_sess' num2str(ss) '_' name{n} '.nii']);
                column_name{n} = fullfile(sprintf('Sess%d_%s.nii',ss,name{n}));
            end
            for s=sn
                for h=1:length(hem);
                    caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                    white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                    pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                    % load all files
                    C1=caret_load(white);
                    C2=caret_load(pial);
                    
                    for f=1:length(fileList)
                        images{f}=fullfile(glmSessDir{ss},subj_name{s},fileList{f});
                    end;
                    metric_out = fullfile(caretSDir,sprintf('%s_Contrasts_sess%d.metric',subj_name{s},ss));
                    M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                    M.column_name = column_name;
                    caret_save(metric_out,M);
                    fprintf('Sess %d: %s - %s done...surface psc done.\n',ss,subj_name{s},hemName{h});
                    
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        Out = caret_smooth(metric_out, 'coord', white, 'topo', closed);
                        char(Out);  % if smoothed adds an 's'
                    end;
                end; % hemi
            end; % sn
        end; % session
       fprintf('\nDone processing case PSC:surface.\n');
    case 'PSC:group_make'
        % Calculate group metric files from contrast / psc calculations. 
        % Takes 2 contrast results ('TrainSeq','UntrainSeq') across
        % subjects and makes a group level metric file that contains each
        % subject's data for that contrast type.
        vararginoptions(varargin,{'sessN','sn'});
        % Some presets
        name        = 'Contrasts';
        OUTname     = {'TrainSeq','UntrainSeq'};
        inputcol    = [13 14]; % only averages - all trained / untrained
        replaceNaN  = [1 1];
        for ss=sessN
            % Loop over hemispheres.
            for h = 1:2
                % initialise some variables
                outfilenames    = cell(length(OUTname),1);
                infilenames     = cell(length(OUTname),length(sn));
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = [caretDir filesep atlasname filesep hemName{h}];
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'INname' and make a group metric file
                for j = 1:length(OUTname);
                    % Loop over subjects...
                    for i = 1:length(sn);
                        % ...and define the names of their metric files
                        infilenames{j}{i} = fullfile(caretDir,[atlasA subj_name{sn(i)}], hemName{h}, sprintf('%s_%s_sess%d.metric',subj_name{sn(i)},name,ss));
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    outfilenames{j} = [surfaceGroupDir filesep hem{h} '.' OUTname{j} '_sess' num2str(ss) '.metric'];
                    % Finally, make the group metric file for this metric type/contrast
                    caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
                    fprintf('Session %d: %s - image: %i/%i group file created.\n',ss,hemName{h},j,length(OUTname));
                end;
            end;
        end;
        fprintf('\nDone processing case PSC:group_make.\n');
    case 'PSC:group_cSPM'
        % Calculate group stats files from the group metric files. 
        % Takes 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross')  and 
        % calculates group level stats (one sample t test that the mean 
        % effect is bigger than zero). 
        % 
        % Although we calculate the t-score, corresponding p-value for that
        % t-score, and finally the z-score 
        % subject's data for that contrast type.
        s=1:length(sn);
        vararginoptions(varargin,{'sessN','sn'});
        
        SPMname={'TrainSeq','UntrainSeq'};    
        sqrtTransform=[0,0,0,0]; % Should you take ssqrt before submitting? 
        % no for psc
        for ss=sessN
            SummaryName = sprintf('.summary_psc_sess%d.metric',ss);
            hemi = [1 2];
            for h=hemi
                % initialise
                filenames = cell(length(SPMname),1);
                sfilenames = filenames;
                column_name = cell(1,length(SPMname)*2);
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
                %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
                for i=1:length(SPMname);
                    filenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '_sess' num2str(ss) '.metric']; % no smoothing
                    sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '_sess' num2str(ss) '.metric']; % smoothing
                end;
                %----loop over the metric files and calculate the cSPM of each with the non-smoothed metrics
                for i=1:length(SPMname);
                    Data=caret_load(filenames{i});
                    if sqrtTransform(i)
                        Data.data=ssqrt(Data.data);
                    end;
                    cSPM=caret_getcSPM('onesample_t','data',Data.data(:,s),'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                    caret_savecSPM([surfaceGroupDir filesep hem{h} '.' SPMname{i} '_stats.metric'],cSPM);
                    save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');
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
        fprintf('\nDone processing case PSC:group_cSPM.\n');
    case 'PSC:group_smooth'
        vararginoptions(varargin,{'sessN'});
        for ss=sessN
            for h=1:2
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                cd(surfaceGroupDir)
                %----define name of coord and topology
                coordfile=[hem{h} '.WHITE.coord'];
                topofile=[hem{h} '.CLOSED.topo'];
                %----get the full directory name of the metric files and the smoothed metric files that we create below
                for i=1:length(ss);
                    filename = [surfaceGroupDir filesep hem{h} '.summary_psc_sess' num2str(ss) '.metric']; % unsmoothed
                    sfilename = caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
                end;  
                fprintf('Sess %d: %s - smoothing done.\n',ss,hemName{h});
            end;
        end 
        fprintf('\nDone processing case PSC:smooth.\n');
    case 'PSC:create_rgb'
        smooth=0;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for ss=sessN
            for h=1:2  
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                cd(surfaceGroupDir)
                if smooth == 1
                    Cp=caret_load(['s' hem{h} sprintf('.summary_psc_sess%d.metric',ss)]);
                else
                    Cp=caret_load([hem{h} sprintf('.summary_psc_sess%d.metric',ss)]);
                end  
                %----PSC
                low_1=0.5;
            %    RGBdata = zeros(size(Cp.data,1),3);
                RGBdata(:,1)=Cp.data(:,1); % Red: trained BOLD signal
                RGBdata(:,3)=Cp.data(:,2); % Blue: untrained BOLD signal
                RGBdata(RGBdata(:,1)<low_1,1)=0;
                RGBdata(RGBdata(:,3)<low_1,3)=0;

                sc1=[low_1 3;low_1 3;low_1 3];   % Scaling for PSC
                % prepare to save
                name={sprintf('psc_sess%d',ss)};                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1},'column_name',name);  
                if smooth == 1
                    caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.PSC_sess%d.RGB_paint',ss)],C);
                else
                    caret_save([surfaceGroupDir filesep hem{h} sprintf('.PSC_sess%d.RGB_paint',ss)],C);
                end
                fprintf('Sess %d: %s - rgb map done.\n',ss,hemName{h});
            end;
        end
        fprintf('\nDone processing case PSC:create_rgb.\n');
    case 'PSC:trained_rgb'
        sessN=1:3;
        smooth=1;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for ss=sessN
            for h=1:2
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
                cd(surfaceGroupDir)
                if smooth == 1
                    Cp=caret_load(['s' hem{h} sprintf('.summary_psc_sess%d.metric',ss)]);
                else
                    Cp=caret_load([hem{h} sprintf('.summary_psc_sess%d.metric',ss)]);
                end  
                %----PSC - only trained sequences
                %specific per session - sess1 yellow, sess2 orange, sess3 red
                switch ss
                    case 1
                        rgb_indx = [255 175 14]/256;
                    case 2
                        rgb_indx = [242 95 3]/256;
                    case 3
                        rgb_indx = [255 20 5]/256;
                end
                low_1=0.5;
                RGBdata = zeros(size(Cp.data,1),3);
                RGBdata(:,1)=Cp.data(:,1);
                RGBdata(:,2)=Cp.data(:,1);
                RGBdata(:,3)=Cp.data(:,1);
                RGBdata(RGBdata(:,1)<low_1,1)=0;
                RGBdata(RGBdata(:,2)<low_1,2)=0;
                RGBdata(RGBdata(:,3)<low_1,3)=0;
                RGBdata(RGBdata(:,1)>low_1,1)=rgb_indx(1);
                RGBdata(RGBdata(:,2)>low_1,2)=rgb_indx(2);
                RGBdata(RGBdata(:,3)>low_1,3)=rgb_indx(3);
                
                sc1=[low_1 1;low_1 1;low_1 1];                 % Scaling for PSC
                % variables for saving
                name={sprintf('psc_trained_sess%d',ss)};                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1},'column_name',name);
                if smooth == 1
                    caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.PSC_trained_sess%d.RGB_paint',ss)],C);
                else
                    caret_save([surfaceGroupDir filesep hem{h} sprintf('.PSC_trained_sess%d.RGB_paint',ss)],C);
                end
                fprintf('Sess %d: %s - traiend rgb map done.\n',ss,hemName{h});
            end;
        end
        fprintf('\nDone processing case PSC:trained_rgb.\n');
    
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
        end
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regLab{r});
                ttestDirect(T.psc,[T.seqType T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & T.sessN==ss);
            end
        end
    
    case 'SEARCH:group_make'                                                % STEP 4.5   :  Make group metric files by condensing subjec contrast metric files
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
       % reg         = 1:8;
        reg          = 3;
        parcelType  = 'Brodmann';
        modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        checkCorr   = 0;        % check if correlation exact and plot predicted Gs
        sessTr      = 1;      % transitions
        corrLim     = [0 1];    % bounds for lower / upper correlation of models
        nModel      = 30;       % number of correlation models (determines how fine grained the corr estimates are)
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessTr','algorithm','parcelType','modelType','nModel'});
        AllReg=[];
        for str = sessTr % session transition
            sessN = [str str+1];
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
                    T = pcm_fitModels_conserv(Data,M,partVec,condVec,runEffect,algorithm);
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
                    %T                       = rmfield(T,'reg');
                    AllReg = addstruct(AllReg,T);
                    fprintf('Done modelType: %s seqType %d sess: %d-%d reg: %d/%d\n\n',modelType,st,sessN(1),sessN(2),r,length(reg));
                end; % seqType
                fprintf('Done reg %d/%d:\tmodelType: %s \tsess: %d-%d\n\n',r,numel(reg),modelType,sessN(1),sessN(2)); 
            end; % region
            fprintf('Done all:\tmodelType: %s \tsess: %d-%d\n\n\n\n',modelType,sessN(1),sessN(2)); 
        end; % session transition
        % save output
       % save(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s.mat',parcelType,modelType)),'-struct','AllReg');
        save(fullfile(stabDir,'PCM_corrModels_TEST4.mat'),'-struct','AllReg');
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
    
    case 'PCM:plot_withinSess'
        reg         = 1:8;
        hemi        = 1;
        parcelType  = 'Brodmann'; 
        vararginoptions(varargin,{'parcelType','reg','hemi'});
        T = load(fullfile(stabDir,sprintf('PCM_withinSess_%s.mat',parcelType)));    
    
        for r=reg
            for ss=sessN
                figure(r)
                subplot(1,4,ss)
                barplot(T.seqType,T.bayesEst_cross,'subset',T.roi==r&T.regSide==hemi&T.sessN==ss);
                plt.match('y'); ylabel('logBayes');
                title(sprintf('%s - sess%d',regname_cortex{r},ss));
            end
            figure(max(reg)+r)
            subplot(121)
            barplot(T.sessN,T.bayesEst_cross,'subset',T.roi==r&T.regSide==hemi&T.seqType==1);
            title(sprintf('%s - trained',regname_cortex{r})); ylabel('logBayes');
            subplot(122)
            barplot(T.sessN,T.bayesEst_cross,'subset',T.roi==r&T.regSide==hemi&T.seqType==2);
            title(sprintf('%s - untrained',regname_cortex{r})); ylabel('logBayes');
            plt.match('y');
            figure(99)
            subplot(1,numel(reg),r)
            style.use('Seq');
            plt.line([T.sessN>3 T.sessN],T.bayesEst_cross(:,4),'split',T.seqType,'subset',T.roi==r & T.regSide==hemi,'leg',{'trained','untrained'});
            if r==1
                ylabel('Sequence-specific log-Bayes');
            else
                ylabel('');
            end
            xlabel('Session'); title(sprintf('%s',regname_cortex{r}));
        end
    case 'PCM:stats_withinSess'
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
        modelType   = 'noSess';
        metric      = 'r_model';
        vararginoptions(varargin,{'modelType','reg','hemi','metric'});
        
        T = load(fullfile(stabDir,sprintf('PCM_toyModels_%s_%s.mat',parcelType,modelType)));
        T = rmfield(T,'reg');
        style.use('Seq');
        for r=reg
            t=getrow(T,T.roi==r & T.regSide==hemi);
            figure
            plt.box(t.sessTr,t.(metric),'split',t.seqType,'plotall',2,'leg',{'trained','control'});
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
        modelType   = 'noSess';
        metric      = 'r_model';
        vararginoptions(varargin,{'modelType','reg','hemi','metric'});
        
        T = load(fullfile(stabDir,sprintf('PCM_toyModels_%s_%s.mat',parcelType,modelType)));
        T = rmfield(T,'reg');
        for r=reg
            t=getrow(T,T.roi==r & T.regSide==hemi);
            fprintf('t-test %s - per session\n',regname_cortex{r});
            ttestDirect(t.(metric),[t.seqType t.SN],2,'paired','split',t.sessTr);
            fprintf('t-test %s - sessTransition (1-2) - trained\n',regname_cortex{r});
            ttestDirect(t.(metric),[t.sessTr t.SN],2,'paired','subset',t.sessTr~=3 & t.seqType==1);
            fprintf('t-test %s - sessTransition (1-2) - untrained\n',regname_cortex{r});
            ttestDirect(t.(metric),[t.sessTr t.SN],2,'paired','subset',t.sessTr~=3 & t.seqType==2);
        end
    case 'PCM:plot_corrModels'
        modelType = 'wSess';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst_cross'; % bayesEst, bayesEst_cross or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s.mat',parcelType,modelType)));
       % T=load(fullfile(stabDir,'PCM_corrModels_TEST4.mat'));
        
    %    T = rmfield(T,'reg');
        if strcmp(modelType,'wSess')
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
            
            figure % sessions 1-3
            subplot(221)
            style.use('SeqShade_small');
            plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessTr==1 & D.model~=1,'leg',{'trained','control'});
            hold on; drawline(0,'dir','horz');
            ylabel('logBayes'); xlabel('Correlation'); title(sprintf('%s - transition 1-2',regname_cortex{r}));
            subplot(222)
            plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessTr==2 & D.model~=1,'leg',{'trained','control'});
            ylabel('logBayes'); xlabel('Correlation'); title(sprintf('%s - transition 2-3',regname_cortex{r}))
            hold on; drawline(0,'dir','horz');
            subplot(223)
            style.use('SessTrained_shade');
            plt.line(D.model,D.metric2,'split',D.sessTr,'subset',D.seqType==1 & D.sessTr<3 & D.model~=1,'leg',{'trans1','trans2'});
            hold on; drawline(0,'dir','horz');
            ylabel('logBayes'); xlabel('Correlation'); title(sprintf('%s - trained seq transition 1-2',regname_cortex{r}));
            subplot(224)
            style.use('SessUntrained_shade');
            plt.line(D.model,D.metric2,'split',D.sessTr,'subset',D.seqType==2 & D.sessTr<3 & D.model~=1,'leg',{'trans1','trans2'});
            hold on; drawline(0,'dir','horz');
            ylabel('logBayes'); xlabel('Correlation'); title(sprintf('%s - control seq transition 1-2',regname_cortex{r}));
            
            figure % sessions 3-4
            subplot(121)
            style.use('SeqShade_small');
            plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessTr==3 & D.model~=1,'leg',{'trained','control'});
            hold on; drawline(0,'dir','horz');
            ylabel('logBayes'); xlabel('Correlation'); title(sprintf('%s - transition 3',regname_cortex{r}));
            subplot(122)
            style.use('Seq');
            plt.box(D.sessTr,D.(metric),'split',D.seqType,'subset',D.sessTr==3 & D.model~=1,'leg',{'trained','control'});
            ylabel('logBayes of seqSpecific');
            
            d = tapply(t,{'seqType','sessTr'},{(metric),'mean'});
            d.metric = d.(metric)(:,2:end);
            [y,~]=min(d.metric,[],2);
            d.metric = bsxfun(@minus,d.metric,(y-0.001));
           % d.metric = bsxfun(@minus,d.metric,mean(d.metric,2));
            d.metric = bsxfun(@rdivide,d.metric,max(d.metric,[],2));
            
            meanEvi = mean(d.(metric)(:,2:end),2); % overall evidence  
            figure
            hold on;
            subplot(3,8,3:8)
            style.use('SeqShade_small');
            plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessTr==1 & D.model~=1,'leg',{'trained','control'});
            XLim1 = get(gca,'XLim'); ylabel('log-likelihood'); title(regname_cortex{r});
            hold on; drawline(0,'dir','horz');
            subplot(3,8,[17,18])
            hold on;
            bar(5,meanEvi(1),'FaceColor',[222,45,38]./255);
            bar(4,meanEvi(4),'FaceColor',[49,130,189]./255);
            bar(2,meanEvi(2),'FaceColor',[251,177,168]/255); 
            bar(1,meanEvi(5),'FaceColor',[158,202,225]/255); camroll(90);
            XLimD = get(gca,'XLim');
            subplot(3,8,19:24)
            hold on
            for i=1:length(d.metric)
                scatter(i+1,5,d.(metric)(1,i+1),[222,45,38]./255,'o','filled','MarkerFaceAlpha',d.metric(1,i));
                scatter(i+1,4,d.(metric)(2,i+1),[49,130,189]./255,'o','filled','MarkerFaceAlpha',d.metric(2,i));
                scatter(i+1,2,d.(metric)(3,i+1),[251,177,168]./255,'o','filled','MarkerFaceAlpha',d.metric(3,i));
                scatter(i+1,1,d.(metric)(4,i+1),[158,202,225]./255,'o','filled','MarkerFaceAlpha',d.metric(4,i));
            end
            set(gca,'YLim',XLimD);
            set(gca,'XLim',XLim1);
            subplot(3,8,11:16)
            style.use('SeqShade_light')
            plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessTr==2 & D.model~=1,'leg',{'trained','control'});
            hold on; drawline(0,'dir','horz'); ylabel('log-likelihood');

            figure(98)
            t.meanLike = mean(t.(metric)(:,2:end),2);
            tt = normData(t,'meanLike');
            subplot(1,numel(reg),r)
            style.use('Seq');
            plt.line([tt.sessTr>2 tt.sessTr],tt.normmeanLike,'split',tt.seqType,'leg',{'trained','control'});
            ylabel('log-Bayes'); title(sprintf('%s - seqSpecific',regname_cortex{r}));
            xlabel('Session');
        end
    case 'PCM:stats_corrModels'
        modelType = 'wSess';
        parcelType = 'Brodmann';
        reg = 1:8;
        hemi = 1;
        metric = 'bayesEst_cross'; % bayesEst, bayesEst_cross or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s.mat',parcelType,modelType)));
      %  T=load(fullfile(stabDir,'PCM_corrModels_TEST4.mat'));
        
    %    T = rmfield(T,'reg');
        if strcmp(modelType,'wSess')
            T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
            T.(metric) = T.(metric)(:,2:end);
        end
        for r=reg
            t = getrow(T,T.regType==r & T.regSide==hemi);
            t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
            t.metric3 = bsxfun(@minus,t.(metric),max(t.(metric),2));
            [corVal,corIdx]=max(t.(metric),[],2);
            fprintf('%s: t-tests (trained vs. control): corModel\n',regname_cortex{r});
            ttestDirect(corIdx,[t.seqType t.SN],2,'paired','split',t.sessTr);
            fprintf('%s: t-tests (trained vs. control): evidence for corModel\n',regname_cortex{r});
            ttestDirect(corVal,[t.seqType t.SN],2,'paired','split',t.sessTr);
        end
    case 'PCM:stats_acrSess_seq'
        modelType = 'wSess';
        parcelType = 'Brodmann';
        reg = [1:8];
        hemi = 1;
        metric = 'bayesEst_cross'; % bayesEst, bayesEst_cross or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric'});
        
        T=load(fullfile(stabDir,sprintf('PCM_corrModels_%s_%s.mat',parcelType,modelType)));
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
     
        
    case 'DIST:plot_mahalanobis'
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'parcelType','roi'});
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)));
        
        T.dist      = [D.dist_train; D.dist_untrain];
        T.sn        = [D.sn; D.sn];
        T.regType   = [D.regType; D.regType];
        T.regSide   = [D.regSide; D.regSide];
        T.sessN     = [D.sessN; D.sessN];
        T.seqType   = [ones(size(D.dist_train));ones(size(D.dist_train))*2];
        T = normData(T,'dist');
        style.use('Seq');
        for r=roi
            figure
            plt.bar([T.sessN>3 T.sessN],T.normdist,'split',T.seqType,'subset',T.regType==r & T.regSide==hemi,'leg',{'trained','control'});
            ylabel('Distance'); xlabel('Session'); title(regname_cortex{r});
        end
    case 'DIST:plot_cosine'
        % here plot cosine distances
        roi=1:8;
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'parcelType','roi'});
        D=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI',parcelType)));
        D=normData(D,'corrDist');
        style.use('Seq');
        for r=roi
            figure
            plt.line([D.sessN>3 D.sessN],D.normcorrDist,'split',D.seqType,'subset',D.regType==r & D.regSide==hemi,'leg',{'trained','control'});
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
                ttestDirect(T.dist,[T.seqType T.sn],1,'paired','subset',T.regType==r&T.regSide==hemi&T.sessN==ss);
                fprintf('trained seq vs. 0\n');
                ttestDirect(T.dist,[T.sn],1,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==1);
                fprintf('untrained seq vs. 0\n');
                ttestDirect(T.dist,[T.sn],1,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==2);
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
        end
        keyboard;
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname_cortex{r});
                ttestDirect(T.corrDist,[T.seqType T.sn],1,'paired','subset',T.regType==r&T.sessN==ss&T.regSide==hemi);
            end
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
        [movTime,subjDay,~] = pivottable([B.day B.SN],[B.seqType],B.MT,'median');
        S = getrow(D,ismember(D.SN,sn)&D.ScanSess>0 & D.seqType<3);
        S.day(S.ScanSess==1)=2;
        S.day(S.ScanSess==2)=7;
        S.day(S.ScanSess==3)=18;
        S.day(S.ScanSess==4)=19;
        [movTimeSc,idx,~] = pivottable([S.day S.SN],[S.seqType],S.MT,'median');
        SS.MT       = movTimeSc(:);
        SS.seqType  = [ones(size(movTimeSc,1),1);ones(size(movTimeSc,1),1)*2];
        SS.day      = [idx(:,1);idx(:,1)];
        SS.sn       = [idx(:,2); idx(:,2)];
        figure
        style.use('SeqShade');
        subplot(121)
        plt.box(SS.day,SS.MT,'split',SS.seqType,'plotall',0,'leg',{'trained','control'});
        % first scanner behaviour
        subplot(122)
        plt.line([subjDay(:,1)>6 subjDay(:,1)],movTime);
        ylabel('Movement time (msec)'); xlabel('Days');
        plt.match('y');
    
    case 'run_job'
     %   sml_plasticity('PCM:withinSess');
     %   fprintf('Done within session!\n\n\n');
      %  sml_plasticity('PCM:corrModels');
      %  fprintf('Done corr models 1 - with session!\n\n\n');
        sml_plasticity('PCM:corrModels','modelType','noSess');
        fprintf('Done corr models 2 - no session!\n\n\n');
        
    case 'TESSEL:make_distMask'
        % here make a mask for distances
        mask_thres = 1.71; % mask threshold (here so that t-map is p<.5)
        maxcsize = 200;
        smoothparam = 15;
        
        for h=1:2
            % first load the surface
            caretGroupDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
            coord = fullfile(caretGroupDir,sprintf('%s.WHITE.coord',hem{h}));
            topo  = fullfile(caretGroupDir,sprintf('%s.CLOSED.topo',hem{h}));
            
            surface = caret_getsurface(coord,topo);
            surface.Edges = surface.Tiles;
            for ss=1:4
                D = caret_load(fullfile(caretGroupDir,sprintf('s%s.summary_dist_sess%d.metric',hem{h},ss)));
                DD(:,ss)=D.data(:,5);
            end
            maxDist(:,1)    = max(DD,[],2); % here the maximum distance 
            idx             = maxDist>=mask_thres; % just binary
            cidx            = caret_clusters(surface,idx);
            cs              = unique(cidx);
            for i = 2:length(cs)
                curridx = cidx==cs(i);
                if sum(curridx)<maxcsize % if cluster is small
                    idx(curridx) = 0;
                end
            end          
            C = caret_struct('metric','data',maxDist,'column_name',{'maxDist'});
            caret_save(fullfile(caretGroupDir,sprintf('%s.distMask_maxDist.metric',hem{h})),C);
            M = caret_struct('metric','data',double(idx),'column_name',{'distance_mask'});
            fname = fullfile(caretGroupDir,sprintf('%s.distMask_binary.metric',hem{h}));
            caret_save(fname,M);
            % smooth now
            fname_smooth = caret_smooth(fname,'coord',coord,'topo',topo,'iterations',smoothparam);
            S = caret_load(fname_smooth{1});
            % save as .mat file
            mask(:,h) = double(S.data>0.1); 
        end
        save(fullfile(regDir,'surf_distmask'),'mask');
    case 'TESSEL:define_nodes'
        % HERE ADD INFO
        figOn=1;
        separation = 7;
        shownum = 0;
        vararginoptions(varargin,{'figOn','separation','shownum'});
        load(fullfile(regDir,'surf_distmask')); % loading in the mask data
        % Load flat map and mask
        for h=1:2 % each hemisphere
            if figOn;figure('units','centimeters','position',[5,5,15,15]);end;    
            groupDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
            cd(groupDir);      
            % -- flat map
            Flat = caret_load([hem{h} '.FLAT.coord']);
            xrange      = [floor(min(Flat.data(:,1))),ceil(max(Flat.data(:,1)))];
            yrange      = [floor(min(Flat.data(:,2))),ceil(max(Flat.data(:,2)))];         
            % -- border files
            Border = caret_load([hem{h} '.display.borderproj']);
            M=caret_load([hem{h} '.propatlas.metric']);
            P=caret_load(['lateral.paint']);
            % Display mask on flat map
            if figOn
                % flat map
                plot(Flat.data(:,1),Flat.data(:,2),'.','color',[0.8 0.8 0.8]);
                hold on; axis equal;
                % border dots
                switch (h)
                    case 1
                        borderidx   = cat(1,Border.Border.vertex);
                        borderX     = Flat.data(borderidx(:),1);
                        borderY     = Flat.data(borderidx(:),2);
                    case 2
                        borderX     = -borderX;
                end
                plot(borderX,borderY,'k.');            
                caxis([0.3,0.7]);
                set(gca,'xlim',xrange,'ylim',yrange);
                axis off
                title(hemName{h})
            end          
            % Generate hexagonal grid
            [nodeidx, edgeidx, sp, XV, YV] = clusterRDM_hexagon_flat('Flat', Flat, 'separation', separation,...
                'mask', mask(:,h));     
            if figOn
                plot(XV,YV,'r');
                for i=1:numel(nodeidx)
                    if shownum
                        text(Flat.data(nodeidx(i),1),Flat.data(nodeidx(i),2),num2str(i),...
                            'horizontalalignment','center','verticalalignment','middle',...
                            'fontsize',10);
                    end
                end
                set(gca,'xlim',xrange,'ylim',yrange);
                axis equal;
            end;     
            Nodes{h} = nodeidx;
            Edges{h} = edgeidx;
        end
        % save result
        save(fullfile(regDir, sprintf('surf_hexNodes_%dmmSeparation.mat',round(separation))),...
            'Nodes', 'Edges');
    case 'TESSEL:define_searchlight'
        % now define the ROIs on the surface
        overlap = 0.2; % 0~1 (1 was too large)
        separation = 7; % mm
        nVox = 120; % if not explicitly given, automatically calculated using separation and overlap parameters
        vararginoptions(varargin,{'overlap','separation','nVox','sn'});

        % Load ultra-reduced node indices (created in 'define_nodes')
        load(fullfile(regDir, sprintf('surf_hexNodes_%dmmSeparation.mat',round(separation)))); % load in Nodes, Edges
        % Calc radius and nVox
        radius  = separation*(0.5*overlap+0.5);
        if isempty(nVox)
            nVox    = 160 * (radius/12)^2; % 2D approximation
            % adjust digit
            nVox = 10*ceil(nVox/10);
        end
        % Loop over subjects
        fprintf('Defining searchlights in subjects:\n');
        for s=sn
            for h=1:2 % per hemisphere
                % load caret files
                caret_subjDIR   = fullfile(caretDir,[atlasA subj_name{s}],hemName{h});
              %  caret_subjDIR   = fullfile(caretDir,'fsaverage_sym',hemName{h});
                coord_pial      = caret_load(fullfile(caret_subjDIR, sprintf('%s.PIAL.coord',hem{h})));
                coord_white     = caret_load(fullfile(caret_subjDIR, sprintf('%s.WHITE.coord',hem{h})));
                topo            = caret_load(fullfile(caret_subjDIR, sprintf('%s.CLOSED.topo',hem{h})));
                % get the mask file
                ref = fullfile(regDir,sprintf('mask_%s.nii',subj_name{s}));
                refV            = spm_vol(ref);
                epiInfo.mat     = refV.mat;
                epiInfo.dim     = refV.dim;
                epiInfo.mask    = spm_read_vols(refV);
                centernode      = Nodes{h}; 
                num_nodes       = coord_white.num_nodes;
                %- run surfing_voxelselection
                [LI,voxmin,voxmax,vORr] = ...
                    surfing_voxelselection(coord_white.data',...
                    coord_pial.data',...
                    topo.data',...
                    [radius nVox],...
                    epiInfo,...
                    centernode,...
                    [5,0,1]);
                % check validity
             %   badidx  = false(size(vORr));%isnan(vORr); don't change SL here otherwise you'll pass unequal SL nodes to the searchlight function.
             %   LI      = LI(~badidx);
             %   voxmin  = voxmin(~badidx,:);
             %   voxmax  = voxmax(~badidx,:);
             %   vORr    = vORr(~badidx,:);
             %   nodeID  = centernode(~badidx)';
                nodeID = centernode';
                %- save searchlight def in .mat format - in ROI folder
                LI = LI';
                save(fullfile(regDir,sprintf('surf_hexNodes_%s_%s_%dmmSeparation_%dvox.mat',hem{h},subj_name{s},separation,nVox)),...
                    'LI','voxmin','voxmax','vORr','nodeID');
                
                %- save searchlight def in .metric format
                M                   = caret_struct('metric','data',vORr);
                M.data              = NaN(size(coord_white.data,1),1);
                M.data(centernode)  = vORr;
                M.column_name       = {'nodes'};
                M.num_rows          = length(M.data);
                SLname = sprintf('%s.hexNodes_%s_%dmmSeparation_%dvox.metric',hem{h},subj_name{s},separation,nVox);
                caret_save(fullfile(caret_subjDIR,SLname),M);
            end
            fprintf('%d.',s);
        end
    case 'TESSEL:select'
        % select tessels so that all the subjects have data in all nodes
        vararginoptions(varargin,{'sn'});        
        for h=1:2 % hemisphere
            idx = 1:1000; % just randomly set up all numbers (for finding overlap)
            for s=sn
                idx_s=[];
                t=load(fullfile(regDir,sprintf('surf_hexNodes_%s_%s_7mmSeparation_120vox',hem{h},subj_name{s})));
                for i=1:size(t.nodeID,1)
                    if ~any(isnan(t.voxmin(i,:)))
                        idx_s = [idx_s i];
                    end
                end
                idx = intersect(idx,idx_s);
            end
            select_subset{h} = idx;
        end
        save(fullfile(regDir,'surf_hexNodes_selection'),'select_subset'); % save the selection
    case 'TESSEL:getBetas'
        % extract betas from defined tessels
        separation=7;
        nVox=120;
        sessN=1:4;
        vararginoptions(varargin,{'nVox','separation','sessN'})
        
        load(fullfile(regDir,'surf_hexNodes_selection')); % here the selection
        for ss=sessN
            for s=sn
                % load data
                SS=[];
                load(fullfile(glmSessDir{ss},subj_name{s},'SPM'));
                V = SPM.xY.VY;
                volIn = spm_vol(V(1).fname); % use the first file for dim estimation
                for h=1:2 % for each hemisphere
                    fprintf('\nExtracting %s - %s:',subj_name{s},hemName{h});
                    T = load(fullfile(regDir,sprintf('surf_hexNodes_%s_%s_%dmmSeparation_%dvox.mat',hem{h},subj_name{s},separation,nVox)));
                    T = getrow(T,select_subset{h});
                    tstart=tic;
                    for t=1:size(T.voxmin,1);  % for each tessel
                        linVox  = unique(cat(2,T.LI{t})');  % find the linear voxel indices
                        [I,J,K] = ind2sub(volIn(1).dim,linVox);
                        % extract the data
                        Y = zeros(length(V),length(I));
                        for i=1:length(V)
                            Y(i,:)=spm_sample_vol(V(i),I,J,K,0);
                        end
                        % prewhiten
                        [betaW,resMS,~,beta]    = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                        S.betaW                 = {betaW};                             % multivariate pw
                        S.betaUW                = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                        S.betaRAW               = {beta};
                        S.resMS                 = {resMS};
                        % other info
                        S.sn                    = s;
                        S.regSide               = h;
                        S.regType               = select_subset{h}(t);
                        S.reg                   = t;
                        S.nodeID                = T.nodeID(t);
                        SS = addstruct(SS,S);
                        fprintf('%d.',t);
                    end; % tessel
                    toc(tstart)
                    keyboard;
                end; % hemi
            end; % subject
        end; % session
    
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