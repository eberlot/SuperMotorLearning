function varargout = sml_plasticity(what,varargin)
% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
caretDir        =[baseDir '/surfaceCaret'];     
behDir          =[baseDir '/behavioral_data'];
pcmDir          =[baseDir '/pcm_stats'];
stabDir         =[baseDir '/stability'];
distPscDir      =[baseDir '/dist_psc_stats'];
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
                    % Verbose display to user
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
            fprintf('\n PSC in %s \n',regLab{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r & T.regSide==hemi);
        end
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regLab{r});
                ttestDirect(T.psc,[T.seqType T.sn],2,'paired','subset',T.regType==r & T.regSide==hemi & T.sessN==ss);
            end
        end
        
        
    case 'Fig1'
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
    case 'PLOT:behaviour'
        
        D = load(fullfile(behDir,'analyze','alldata.mat'));
        B = getrow(D,ismember(D.SN,sn)&D.blockType==6);
        [movTime,subjDay,~] = pivottable([B.day B.SN],[B.seqType],B.MT,'median');
        S = getrow(D,ismember(D.SN,sn)&D.ScanSess>0);
        S.day(S.ScanSess==1)=2;
        S.day(S.ScanSess==2)=7;
        S.day(S.ScanSess==3)=18;
        S.day(S.ScanSess==4)=19;
        [movTimeSc,subjDaySc,~] = pivottable([S.day S.SN],[S.seqType],S.MT,'median');
        figure
        style.use('SeqShade');
        plt.bar(subjDaySc(:,1),movTimeSc,'split',)
        % first scanner behaviour
        hold on;
        plt.line([subjDay(:,1)>6 subjDay(:,1)],movTime);
        ylabel('Movement time (msec)'); xlabel('Days');
        
    otherwise
        fprintf('No such case!\n');
end
