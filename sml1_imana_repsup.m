function varargout=sml1_imana_repsup(what,varargin)

% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
behavDir        =[baseDir '/behavioral_data/data'];  
anaDir          =[baseDir '/behavioral_data/analyze'];          
imagingDir      =[baseDir '/imaging_data'];                      
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia_new'];
pcmDir          =[baseDir '/pcm_stats'];
repSupDir       =[baseDir '/repsup_stats'];

% update glmDir when adding new glms
glmFoSExDir     ={[baseDir '/glmFoSEx/glmFoSEx1'],[baseDir '/glmFoSEx/glmFoSEx2'],[baseDir '/glmFoSEx/glmFoSEx3'],[baseDir '/glmFoSEx/glmFoSEx4']};    
glmTrialDir     ={[baseDir '/glmTrial/glmTrial1'],[baseDir '/glmTrial/glmTrial2'],[baseDir '/glmTrial/glmTrial3'],[baseDir '/glmTrial/glmTrial4']};

 
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

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22',...
                's23','s24','s25','s26','s27','s28','s29','s30','s31'};  
num_run = 8; % 8 functional runs
% Other random notes %

% -------------------------- Sequence stuff -------------------------------
Seq = [1 5 3 5 3 4 3 1 5; ...
        1 5 2 1 4 5 3 4 2; ...
        3 4 2 1 4 5 1 5 2; ...
        3 1 5 3 4 2 1 4 5; ...
        5 3 4 5 4 2 1 5 3; ...
        5 4 2 3 4 2 3 1 5; ...
        5 1 3 1 3 2 3 5 1; ...
        5 1 4 5 2 1 3 2 4; ...
        3 2 4 5 2 1 5 1 4; ...
        3 5 1 3 2 4 5 2 1; ...
        1 3 2 1 2 4 5 1 3; ...
        1 2 4 3 2 4 3 5 1];
% for group 1 - group 2 (1-6 and 7-12 reversed)
SeqChunks = [1,2,3; ...
    4,5,6; ...
    6,5,4; ...
    3,6,5; ...
    2,7,1; ...
    7,6,3; ...
    8,9,10; ...
    11,12,13;...
    13,12,11;...
    10,13,12;...
    9,14,8;...
    14,13,10];

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
styRS = style.custom({gray,silver,lightgray});
styRSbeh = style.custom({red,lightred,blue,lightblue},'markersize',ms);

%% --------------------------------------- Analysis Cases -----------------------------------------------
switch(what)


    case 'BEH_scanner_MT' %----------------------- QUANTIFY BEHAVIOUR IN SCANNER -----------------------
        sn=[4:9,11:31];
        sessN=[1:4];
        
        vararginoptions(varargin,{'sn','sessN'});
        
        B_fun=[];
        
        for ss = sessN
            for s = sn
                D = load(fullfile(anaDir,['sml1_',subj_name{s}]));
                
                % functional runs in the scanner
                R = getrow(D,D.ScanSess==ss & (D.blockType==3 | D.blockType==9));
                % calculate / save variables
                % Func - MT, ER (for trained and untrained separately)
                BF.MT = R.MT;
                BF.ER = R.isError;
                BF.sessN = R.ScanSess;
                BF.seqType = R.seqType;
                BF.FoSEx = R.FoSEx;
                BF.sn = s*ones(size(R.ScanSess));
                BF.trialNum = R.TN;
                
                block_indx = unique(R.BN);
                BF.run=ones(size(BF.trialNum));
                for b=1:length(block_indx)
                    BF.run(R.BN==block_indx(b)) = b;
                end
                % Add all var into structures B_loc and B_fun
                
                B_fun = addstruct(B_fun,BF);
            end
        end

         % % save
        save(fullfile(repSupDir,'beh_scanner.mat'),'-struct','B_fun');   
    case 'PLOT_behScanner'
        
        B = load(fullfile(repSupDir,'beh_scanner.mat'));
        
        figure
        plt.bar(B.sessN,B.MT,'split',[B.seqType B.FoSEx],'leg',{'1st train','2nd train','1st untrain','2nd untrain'},'leglocation','northeast'); 
        xlabel('Session'); ylabel('Movement time');
        
        figure
        plt.bar(B.sessN,B.ER,'split',[B.seqType B.FoSEx],'leg',{'1st train','2nd train','1st untrain','2nd untrain'},'leglocation','northeast'); 
        xlabel('Session'); ylabel('Error rate');
    
        % reduced structure
        [fa ra ca] = pivottable([B.sessN B.sn],[B.FoSEx B.seqType],B.MT,'mean');
        N.MT = [fa(:,1);fa(:,2);fa(:,3);fa(:,4)];
        N.FoSEx = [ones(size(fa(:,1)));ones(size(fa(:,1)));ones(size(fa(:,1)))*2;ones(size(fa(:,1)))*2];
        N.seqType = [ones(size(fa(:,1)));ones(size(fa(:,1)))*2;ones(size(fa(:,1)));ones(size(fa(:,1)))*2];
        N.sn = [ra(:,2);ra(:,2);ra(:,2);ra(:,2)];
        N.sessN = [ra(:,1);ra(:,1);ra(:,1);ra(:,1)];
        
        figure
        plt.dot(N.sessN,N.MT,'split',[N.seqType N.FoSEx],'leg',{'1st train','2nd train','1st untrain','2nd untrain'},'leglocation','northeast','style',styRSbeh);
        %
        figure
        plt.box(N.sessN,N.MT,'split',[N.seqType N.FoSEx],'plotall',0,'leg',{'1st train','2nd train','1st untrain','2nd untrain'},'leglocation','northeast','style',styRSbeh);
        %
    
    case 'PSC_create_RepSup' %-------------------- CALCULATE PSC, PROJECT ON SURFACE -------------------
        % calculate psc for trained and untrained sequences (1st/2nd) - based on betas    
        vararginoptions(varargin,{'sn','sessN'});
        name={'TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
        for ss=sessN
            for s=sn
                cd(fullfile(glmFoSExDir{ss}, subj_name{s}));
                load SPM;
                T=load('SPM_info.mat');
                X=(SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
                h=median(max(X));               % Height of response;
                P={};
                numB=length(SPM.xX.iB);         % Partitions - runs
                for p=SPM.xX.iB
                    P{end+1}=sprintf('beta_%4.4d.nii',p);       % get the intercepts and use them to calculate the baseline (mean images)
                end;
                for con=1:length(name)    % 4 contrasts
                    P{numB+1}=sprintf('con_%s.nii',name{con});
                    outname=sprintf('psc_sess%d_%s.nii',ss,name{con}); % ,subj_name{s}
                    
                    formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);
                    
                    spm_imcalc_ui(P,outname,formula,...
                        {0,[],spm_type(16),[]});        % Calculate percent signal change
                end;
                fprintf('Subject %d sess %d: %3.3f\n',s,ss,h);
            end;
        end
    case 'PSC_surface_RepSup'
     % create surface maps of percent signal change 
     % trained and untrained sequences - 1st / 2nd execution
     sn=[4:9,11:31];
     sessN=[1:4];
     smooth = 0;
     vararginoptions(varargin,{'sn','sessN','smooth'});
     
     hemisphere=1:length(hem);
     fileList = [];
     column_name = [];
     name={'TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
     for ss=sessN
         for n = 1:length(name)
             fileList{n}=fullfile(['psc_sess' num2str(ss) '_' name{n} '.nii']);
             column_name{n} = fullfile(sprintf('Sess%d_RS_%s.nii',ss,name{n}));
         end
         for s=sn
             for h=hemisphere
                 caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                 specname=fullfile(caretSDir,['x',subj_name{s} '.' hem{h}   '.spec']);
                 white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                 pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                 
                 C1=caret_load(white);
                 C2=caret_load(pial);
                 
                 for f=1:length(fileList)
                     images{f}=fullfile(glmFoSExDir{ss},subj_name{s},fileList{f});
                 end;
                 metric_out = fullfile(caretSDir,sprintf('%s_Contrasts_RS_sess%d.metric',subj_name{s},ss));
                 M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                 M.column_name = column_name;
                 caret_save(metric_out,M);
                 fprintf('Subj %d, Sess%d, Hem %d\n',s,ss,h);
                 
                 if smooth == 1;
                     % Smooth output .metric file (optional)
                     % Load .topo file
                     closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                     Out = caret_smooth(metric_out, 'coord', white, 'topo', closed);%,...
                     %'algorithm','FWHM','fwhm',12);
                     char(Out);  % if smoothed adds an 's'
                 else
                 end;
                 
             end;
         end;
     end
    case 'PSC_group_make_RepSup'
        % Calculate group metric files from contrast / psc calculations. 
        % Takes 4 contrast results (1st/2nd exe for trained/untrained) across
        % subjects and makes a group level metric file that contains each
        % subject's data for that contrast type.
        sessN=[1:4];
        sn=[4:9,11:31];
        vararginoptions(varargin,{'sessN','sn'});
        % Some presets
        name    = 'Contrasts_RS';

        OUTname    = {'RS_TrainSeq_1st','RS_TrainSeq_2nd','RS_UntrainSeq_1st','RS_UntrainSeq_2nd'};
        inputcol   = [1 2 3 4]; % only averages - all trained / untrained
        replaceNaN = [1 1 1 1];
        
        for ss=sessN
            % Loop over hemispheres.
            for h = 1:2
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
                    fprintf('sess-%d: hem: %i  image: %i \n',ss,h,j);
                end;
            end;
        end;
    case 'PSC_group_cSPM_RepSup'
        % Calculate group stats files from the group metric files. 
        % Takes 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross')  and 
        % calculates group level stats (one sample t test that the mean 
        % effect is bigger than zero). 
        % 
        % Although we calculate the t-score, corresponding p-value for that
        % t-score, and finally the z-score 
        % subject's data for that contrast type.
        
        sessN=[1:4];
        sn=[4:9,11:31];
        vararginoptions(varargin,{'sessN','sn'});
        s=1:length(sn);

        SPMname={'TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
        
        sqrtTransform=[0,0,0,0]; % Should you take ssqrt before submitting? 
        % no for psc
        for ss=sessN
            SummaryName = sprintf('.summary_RS_psc_sess%d.metric',ss);
            hemi = [1 2];
            
            for h=hemi
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
                %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
                for i=1:length(SPMname);
                    filenames{i}=[surfaceGroupDir filesep hem{h} '.RS_' SPMname{i} '_sess' num2str(ss) '.metric']; % no smoothing
                    %sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '_RS_sess' num2str(sessN) '.metric']; % smoothing
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
            fprintf('Done sess-%d\n',ss);
        end
    case 'PSC_group_smooth_RepSup'
        sessN=[1:4];
        vararginoptions(varargin,{'sessN'});
        
        for ss=sessN
            for h=1:2
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                cd(surfaceGroupDir)
                C=caret_load([hem{h} sprintf('.summary_RS_psc_sess%d.metric',ss)]);
                %----define name of coord and topology
                coordfile=[hem{h} '.WHITE.coord'];
                topofile=[hem{h} '.CLOSED.topo'];
                %----get the full directory name of the metric files and the smoothed metric files that we create below
                filename=[surfaceGroupDir filesep hem{h} '.summary_RS_psc_sess' num2str(ss) '.metric']; % unsmoothed
                sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
            end;
        end
    case 'SURF_flatmap_RS_absolute'
        % plotting RS on surface as absolute
        % subtracting 2nd - 1st activation - ABSOLUTE RS
        sessN=1;
        smooth=1;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for ss=sessN
            for h=1:2
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                cd(surfaceGroupDir)
                if smooth == 1
                    C=caret_load(['s' hem{h} sprintf('.summary_RS_psc_sess%d.metric',ss)]);
                else
                    C=caret_load([hem{h} sprintf('.summary_RS_psc_sess%d.metric',ss)]);
                end
                flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
                topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
                
                %----PSC
                low_1=0.05;
                RGBdata(:,1)=C.data(:,1)-C.data(:,2); % Red: trained 1st-2nd BOLD signal
                RGBdata(:,3)=C.data(:,3)-C.data(:,4); % Blue: untrained 1st-2nd BOLD signal
                RGBdata(RGBdata(:,1)<low_1,1)=0;
                RGBdata(RGBdata(:,3)<low_1,3)=0;
                
                sc=[low_1 0.4;low_1 0.4;low_1 0.4];  % scaling
                
                name={sprintf('RS_absolute_psc_sess%d',ss)};
                %name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN),sprintf('dist_cross_sess%d',sessN)};
                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
                
                if smooth == 1
                    caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.RS_sess%d.RGB_paint',ss)],C);
                else
                    caret_save([surfaceGroupDir filesep hem{h} sprintf('.RS_sess%d.RGB_paint',ss)],C);
                end
            end;
        end
    case 'SURF_flatmap_RS_metric'
        % plotting RS on surface as absolute
        % subtracting 2nd - 1st activation - ABSOLUTE RS
        sessN=1;
        imageType='psc'; % psc or dist
        seqType={'trained','untrained'};
        vararginoptions(varargin,{'sessN','imageType'});
        
        for h=1:2
            
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            
            C=caret_load([hem{h} sprintf('.summary_RS_%s_sess%d.metric',imageType,sessN)]);
            
            flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
            topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
            
            for st = 1:2 % trained / untrained
                data(:,st) = C.data(:,(st-1)*2+2)-C.data(:,(st-1)*2+1); % 2nd-1st
                colIndx{st} = sprintf('RS_diff_%s_%s_sess%d',imageType,seqType{st},sessN);
            end
            
            C2=caret_struct('RGBpaint','data',data,'column_name',colIndx);
            C2.index=C.index;
            C2.column_color_mapping=C.column_color_mapping(1:2,:);
            caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_DIFF_sess%d.metric',imageType,sessN)],C2);
            clear data;
        end;
    case 'SURF_flatmap_RS_smooth'
        sessN=4;
        imageType='psc'; % psc or dist
        vararginoptions(varargin,{'sessN','imageType'});
        
        for h=1:2
            
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            C=caret_load([hem{h} sprintf('.%s_DIFF_sess%d.metric',imageType,sessN)]);
            %----define name of coord and topology
            coordfile=[hem{h} '.WHITE.coord'];
            topofile=[hem{h} '.CLOSED.topo'];
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for i=1:length(sessN);
                filename=fullfile(surfaceGroupDir, sprintf('%s.%s_DIFF_sess%d.metric',hem{h},imageType,sessN)); % unsmoothed
                sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
            end;

        end;
    case 'SURF_flatmap_RS_relative'
        % magnitude of RS in terms of percent - 2nd / 1st
        % RELATIVE RS
        % use the absolute version rather - easier to interpret surface plots
        
        sessN=1;
        smooth=1;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for h=1:2
            
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            if smooth == 1
                C=caret_load(['s' hem{h} sprintf('.summary_RS_psc_sess%d.metric',sessN)]);
            else
                C=caret_load([hem{h} sprintf('.summary_RS_psc_sess%d.metric',sessN)]);
            end
            flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
            topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
            
             %----PSC
            low_1=0.1; % RS of 5%
            high_1=0.5; % RS of 40%
            RGBdata(:,1)=1-C.data(:,2)./C.data(:,1); % Red: 1 - trained 2nd/1st BOLD signal  -> ratio of how much RS present - 0.2 - 20% suppression
            RGBdata(:,3)=1-C.data(:,4)./C.data(:,3); % Blue: 1 - untrained 2nd/1st BOLD signal 
            RGBdata(RGBdata(:,1)<low_1,1)=0;
            RGBdata(RGBdata(:,3)<low_1,3)=0;

            sc=[low_1 high_1;low_1 high_1;low_1 high_1]; 
            
            name={sprintf('RS_relative_psc_sess%d',sessN)};
            %name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN),sprintf('dist_cross_sess%d',sessN)};

            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
            
            if smooth == 1
                caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.RS_relative_sess%d.RGB_paint',sessN)],C);
            else
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.RS_relative_sess%d.RGB_paint',sessN)],C);
            end
        end;
    
    case 'BETA_get'                                                     % STEP 5.6   :  Harvest betas from rois (raw, univ, multiv prewhit)    
        sessN = [1:4];
        sn  = [4:9,11:31];    
        roi = 1:16; 
        roiDefine = 'all'; % determine all regions from region file R
        parcelType = 'Brodmann'; % 162tessels, Brodmann, cortex_buckner, BG-striatum, thalamus
        vararginoptions(varargin,{'sn','sessN','roi','parcelType','roiDefine'});
       
        for ss=sessN            
            % harvest
            for s=sn % for each subj
                T=[];
                fprintf('\nSubject: %d\n',s) % output to user
                % load files
                load(fullfile(glmFoSExDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)]));  % load subject's region parcellation (R)
                
                if strcmp(roiDefine,'all')==1
                    roi=1:size(R,2);
                end
                cd(fullfile(glmFoSExDir{ss},subj_name{s})); % maybe change back when remove glm
                
                P=SPM.Vbeta(SPM.xX.iC);
                
                % Add a few extra images
                %----task against rest
                O{1}=sprintf('psc_sess%d_TrainSeq_1st.nii',ss); %psc trained - 1st execution
                O{2}=sprintf('psc_sess%d_UntrainSeq_1st.nii',ss); %psc untrained - 1st execution
                O{3}=sprintf('psc_sess%d_TrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                O{4}=sprintf('psc_sess%d_UntrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                oP=spm_vol(char(O));
                
                V = SPM.xY.VY;
                
                for r = roi % for each region
                    % get raw data for voxels in region
                    % determine if any voxels for that parcel
                    if size(R{r}.data,1)==0 % no voxels
                        % make data into NaN
                        S.betaW                   = {NaN};
                        S.betaUW                  = {NaN};
                        S.betaRAW                 = {NaN};
                        S.resMS                   = {NaN};
                        S.psc_train_1st     = NaN;
                        S.psc_train_2nd     = NaN;
                        S.psc_untrain_1st   = NaN;
                        S.psc_untrain_2nd   = NaN;
                    else
                        Y = region_getdata(V,R{r});  % Data Y is N x P
                        data = region_getdata(oP,R{r}); % from added images
                        % exclude any missing data in voxels
                        idx = find(Y(1,:));
                        % estimate region betas
                        [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y(:,idx),SPM,'normmode','overall');
                        S.betaW                   = {betaW};                             % multivariate pw
                        S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                        S.betaRAW                 = {beta};
                        S.resMS                   = {resMS};
                        % info from maps for surface
                        S.psc_train_1st   = {data(1,:)};
                        S.psc_untrain_1st = {data(2,:)};
                        S.psc_train_2nd   = {data(3,:)};
                        S.psc_untrain_2nd = {data(4,:)};
                    end
                    
                    % voxel position
                    S.volcoord = {R{r}.data'};
                    S.SN                      = s;
                    S.region                  = r;
                    if r<(numel(roi)/2)+1
                        S.regSide = 1;
                        S.regType = S.region;
                    else
                        S.regSide = 2;
                        S.regType = S.region-(numel(roi)/2);
                    end
                    if any(strcmp(parcelType,{'BG-striatum','thalamus'}))
                        S.regName = {R{r}.name};
                    end
                    T = addstruct(T,S);
                    fprintf('%d.',r)
                    %fprintf('elapsed %d\n',telapsed);
                end
                dircheck(fullfile(betaDir,subj_name{s}));
                save(fullfile(betaDir,subj_name{s},sprintf('betas_FoSEx_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'-struct','T');
                fprintf('\nDone beta extraction for sess%d-%s\n',ss,subj_name{s});
            end            
        end
    case 'BETA_combineGroup'
        % combine individual subject beta structures into the whole
        % structure
        sessN=[1:4];
        sn=[4:9,11:31];
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type','parcelType'});
        
        for ss=sessN
            switch(type)
                case 'new'
                    T=[];
                case 'add'
                    T=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            end   
            fprintf('subjects added for sess-%d:\n',ss);
            for s=sn
              S=load(fullfile(betaDir,subj_name{s},sprintf('betas_FoSEx_%s_%s_sess%d',parcelType,subj_name{s},ss)));
              T=addstruct(T,S);

                fprintf('%d.',s);
            end
            save(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)),'-struct','T','-v7.3');
        end
    case 'ROI_betas_RepSup' %--------------------------- REPSUP BETAS ------------------------------------
       % make a beta structure with betas for each subject / ROI
        sessN = 1;
        sn  = [4:9,11:25];    
        roi = [1:16];
        roiDefine = 'all'; % determine all regions from region file R
        parcelType='cortex_buckner';  % 162tessels, Brodmann, cortex_buckner
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        vararginoptions(varargin,{'sn','sessN','roi','type','parcelType','roiDefine'});
        
        for ss=sessN
            switch(type)
                case 'new'
                    T=[];
                case 'add'
                    T=load(fullfile(regDir,sprintf('betas_FoSEx_sess%d.mat',ss)));
            end
            if exist(fullfile(regDir,sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)))
                fprintf('ROI betas already defined for sess-%d - skipping\n',ss);
            else
                % harvest
                for s=sn % for each subj
                    fprintf('\nSubject: %d\n',s) % output to user
                    
                    % load files
                    load(fullfile(glmFoSExDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                    load(fullfile(regDir,[subj_name{s} '_' parcelType '_regions.mat']));        % load subject's region parcellation (R)
                    V = SPM.xY.VY;
                    
                    if strcmp(roiDefine,'all')==1
                        roi=1:size(R,2);
                    end
                    
                    glmSubjDir = fullfile(glmFoSExDir{ss},subj_name{s});
                    cd(glmSubjDir);
                    
                    % Add a few extra images
                    %----task against rest
                    O{1}=sprintf('psc_sess%d_TrainSeq_1st.nii',ss); %psc trained - 1st execution
                    O{2}=sprintf('psc_sess%d_UntrainSeq_1st.nii',ss); %psc untrained - 1st execution
                    O{3}=sprintf('psc_sess%d_TrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                    O{4}=sprintf('psc_sess%d_UntrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                    oP=spm_vol(char(O));
                    
                    
                    for r = roi % for each region
                        if size(R{r}.data,1)==0 % no voxels
                            % make data into NaN
                            S.betaW                   = NaN;
                            S.betaUW                  = NaN;
                            S.betaRAW                 = NaN;
                            S.resMS                   = NaN;
                            S.psc_train_1st     = NaN;
                            S.psc_train_2nd     = NaN;
                            S.psc_untrain_1st   = NaN;
                            S.psc_untrain_2nd   = NaN;
                        else
                            % get raw data for voxels in region
                            Y = region_getdata(V,R{r});  % Data Y is N x P
                            % exclude any missing data in voxels
                            idx = find(Y(1,:));
                            data = region_getdata(oP,R{r}); % from added images
                            % estimate region betas
                            [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y(:,idx),SPM,'normmode','overall');
                            S.betaW                   = {betaW};                             % multivariate pw
                            S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                            S.betaRAW                 = {beta};
                            S.resMS                   = {resMS};
                            
                            % info from maps for surface
                            S.psc_train_1st   = {data(1,:)};
                            S.psc_untrain_1st = {data(2,:)};
                            S.psc_train_2nd   = {data(3,:)};
                            S.psc_untrain_2nd = {data(4,:)};
                        end
                        S.SN                      = s;
                        S.region                  = r;
                        if strcmp(parcelType,'162tessels')
                            if r<159
                                S.regSide = 1;
                                S.regType = S.region;
                            else
                                S.regSide = 2;
                                S.regType = S.region-158;
                            end
                            if strcmp(roiDefine,'all')==1
                                S.tesselName          = {R{r}.name};
                            end
                        elseif strcmp(parcelType,'Brodmann')
                            S.flatcoord = {R{r}.flatcoord'};
                            S.depth = {R{r}.depth'};
                            if r<9
                                S.regSide = 1;
                                S.regType = r;
                            else
                                S.regSide = 2;
                                S.regType = r-8;
                            end
                        elseif strcmp(parcelType,'cortex_buckner')
                            if r<8
                                S.regSide = 1;
                                S.regType = r;
                            else
                                S.regSide = 2;
                                S.regType = r-7;
                            end
                        end
                        T = addstruct(T,S);
                        fprintf('%d.',r)
                    end
                end
                % save T
                save(fullfile(regDir,sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)),'-struct','T','-v7.3');
                fprintf('\n');
            end
        end
    case 'ROI_stats_RepSup'
       % make a summary stats structure from betas for each subject / ROI
        sessN = 1;
        sn  = [4:9,11:31];
        roi = [1:14];
        roiDefine = 'all'; % determine from region file
        betaChoice = 'multi'; % uni, multi or raw
        type='new';
        parcelType = 'Brodmann'; % 162tessels, Brodmann, BG-striatum
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type','parcelType','roiDefine'});
            
        for ss=sessN    
            T = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            if strcmp(roiDefine,'all')==1
                roi=unique(T.region)';
            end
            % output structures
            switch(type)
                case 'new'
                    To=[];
                case 'add'
                    To=load(fullfile(regDir,sprintf('stats_FoSEx_%sPW_sess%d.mat',betaChoice,ss)));
            end
            
            % do stats
            for s = sn % for each subject
                D = load(fullfile(glmFoSExDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                fprintf('\nSubject: %d session: %d\n',s, ss)
                num_run = numruns_task_sess;
                
                for r = roi % for each region
                    S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                    fprintf('%d.',r)
                    
                    for exe = 1:2   % FoSEx
                        % check first if betas are only nan
                        if size(S.betaW{1},1)==1 & isnan(S.betaW{1})
                            fprintf('..no voxels-skipping..');
                        else
                            switch (betaChoice)
                                case 'uni'
                                    betaW  = S.betaUW{1}(D.FoSEx==exe,:);
                                case 'multi'
                                    betaW  = S.betaW{1}(D.FoSEx==exe,:);
                                case 'raw'
                                    betaW  = S.betaRAW{1}(D.FoSEx==exe,:);
                            end
                            
                            D_exe = getrow(D,D.FoSEx==exe);     
                            
                            % % To structure stats (all seqNumb - 12 conditions)
                            % crossval second moment matrix
                            [G,Sig]     = pcm_estGCrossval(betaW(1:(12*num_run),:),D_exe.run,D_exe.seqNumb);
                            So.IPM      = rsa_vectorizeIPM(G);
                            So.Sig      = rsa_vectorizeIPM(Sig);
                            So.IPMfull  = rsa_vectorizeIPMfull(G);
                            % squared distances
                            So.RDM_nocv = distance_euclidean(betaW',D_exe.seqNumb)';
                            So.RDM      = rsa.distanceLDC(betaW,D_exe.run,D_exe.seqNumb);
                            % trained seq
                            H=eye(6)-ones(6,6)./6;  % centering matrix!
                            [G_train, Sig_train] = pcm_estGCrossval(betaW(D_exe.seqType==1,:),D.run(D_exe.seqType==1),D.seqNumb(D_exe.seqType==1));
                            G_trainCent = H*G_train*H;  % double centered G matrix - rows and columns
                            So.eigTrain = sort(eig(G_trainCent)','descend');    % sorted eigenvalues
                            % untrained seq
                            [G_untrain, Sig_untrain] = pcm_estGCrossval(betaW(D_exe.seqType==2,:),D.run(D_exe.seqType==2),D.seqNumb(D_exe.seqType==2));
                            G_untrainCent = H*G_untrain*H;  % double centered G matrix - rows and columns
                            So.eigUntrain = sort(eig(G_untrainCent)','descend');
                            
                            % stats from additional images
                            if exe==1
                                So.psc_train = nanmean(S.psc_train_1st{:});
                                So.psc_untrain = nanmean(S.psc_untrain_1st{:});
                            else
                                So.psc_train = nanmean(S.psc_train_2nd{:});
                                So.psc_untrain= nanmean(S.psc_untrain_2nd{:});
                            end
                            
                            % crossvalidated 'activation' estimate
                            % length of vector
                            vecLength = diag(G);
                            So.act_vecLength_train = nanmean(vecLength(1:6));
                            So.act_vecLength_untrain = nanmean(vecLength(7:12));
                            % indexing fields
                            So.SN       = s;
                            So.region   = r;
                            if r<(numel(roi)/2)+1
                                So.regSide = 1;
                                So.regType = So.region;
                            else
                                So.regSide = 2;
                                So.regType = So.region-(numel(roi)/2);
                            end
                            So.FoSEx    = exe;
                            To          = addstruct(To,So);
                        end
                    end; % FoSEx
                end; % each region
            end; % each subject
            
            % % save
            save(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)),'-struct','To');
            
            fprintf('\nDone stats for session %d.\n',ss);
        end
    case 'ROI_stats_redIPM_exe' 
        % reduced summary structure - MEAN TRAINED vs. MEAN UNTRAINED
        % G / IPM calculated as a an average version of 2 sequences
        % focusing on 1st vs. 2nd for subset of 2 sequences
        % separately for trained and untrained
        % e.g. IPM: x y z 
        %           x: var of 1st seq, 
        %           y: var of 2nd seq,
        %           z: cov of 1st / 2nd seq
        % all possible combinations of 2 sequences, creating average IPM
        % resultant G (Gc)
        % A1,A1  B1,A1  A2,A1, B2,A1
        % organisation: seq A - B (first exe); seq A - B (second exe)
        sessN=[1:4];
        sn  = [4:9,11:31];
        roi = [1:16];
        betaChoice = 'multi'; % uni, multi or raw
        parcelType = 'Brodmann'; % Brodmann or BG-striatum
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','parcelType'});
        
        for ss=sessN
             T = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss))); % loads region data (T)   
         
            % output structures
            Ts = [];
            To = []; 
            % do stats
            for s = sn % for each subject
                D = load(fullfile(glmFoSExDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                fprintf('\nSubject: %d session: %d\n',s, ss)
                num_run = numruns_task_sess;
                
                for r = roi % for each region
                    S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                    fprintf('%d.',r)
                    
                    for seqType = 1:2   % sequence type - 1: trained, 2: untrained
                        switch (betaChoice)
                            case 'uni'
                                betaW  = S.betaUW{1}(D.seqType==seqType,:);
                            case 'multi'
                                betaW  = S.betaW{1}(D.seqType==seqType,:);
                            case 'raw'
                                betaW  = S.betaRAW{1}(D.seqType==seqType,:);
                        end
                        
                        D_seq = getrow(D,D.seqType==seqType);
                        if seqType==2
                            D_seq.seqNumb=D_seq.seqNumb-6;  % ensuring it's always 1-6
                        end
                        
                        % new index - taking into account sequence number +
                        % repetition
                        % 1-6 first execution, 7-12 second
                        % same for trained and untrained
                        D_seq.seqExe = D_seq.seqNumb;
                        D_seq.seqExe(D_seq.FoSEx==2)=D_seq.seqNumb(D_seq.FoSEx==2)+6;
                        
                        
                        % % To structure stats (all seqNumb - 12 conditions)
                        % crossval second moment matrix
                        [G,Sig]     = pcm_estGCrossval(betaW(1:(12*num_run),:),D_seq.run,D_seq.seqExe);
                        
                        % pick all possible combinations of two sequences
                        comb=nchoosek([1:6],2);
                        for c=1:size(comb,1)
                            % indx for sequence A and sequence B
                            c1=comb(c,1); % seqA
                            c2=comb(c,2); % seqB
                            % seqA second exe - c1+6
                            % seqB second exe = c2+6
                            % create a 4x4 G matrix - same sequence 1st/2nd
                            Gc=[G(c1,c1),G(c1,c2),G(c1,c1+6),G(c1,c2+6);
                                G(c2,c1),G(c2,c2),G(c2,c1+6),G(c2,c2+6);
                                G(c1+6,c1),G(c1+6,c2),G(c1+6,c1+6),G(c1+6,c2+6);
                                G(c2+6,c1),G(c2+6,c2),G(c2+6,c1+6),G(c2+6,c2+6)];
                            % vectorise it
                            IPM(c,:)=rsa_vectorizeIPM(Gc);
                            IPM_full(c,:)=rsa_vectorizeIPMfull(Gc);
                        end
                        
                        % generate average IPM matrix (use rsa_squareIPM to
                        % obtain G)
                        So.IPM_red = mean(IPM);
                        So.IPM_red_full = mean(IPM_full);
                        
                        % indexing fields
                        So.SN       = s;
                        So.region   = r;
                        So.regSide  = regSide(r);
                        So.regType  = regType(r);
                        So.seqType  = seqType;
                        To          = addstruct(To,So);
                        
                    end; % seqType
                end; % each region
            end; % each subject
            
            % % save
            save(fullfile(regDir,sprintf('reduced_IPM_FoSEx_%s_sess%d.mat',parcelType,ss)),'-struct','To');
            fprintf('\nDone for session %d.\n',ss);
        end
    
    case 'CALC_act_dist' %-------------------------- ACTIVATION, MAHALANOBIS DIST  ----------------------
        % summary structure for 1st / 2nd execution - act, mahalanobis dist
        sn = [4:9,11:31];
        roi = [1:16];
        sessN = 1:4;
        betaChoice = 'multiPW';
        parcelType='cortex_buckner';  % 162tessels, Brodmann, cortex_buckner
        subjfig = 0;
        vararginoptions(varargin,{'sn','roi','seq','sessN','betaChoice','fig','parcelType'});

        Stats = [];
        
        for ss = sessN % do per session number
            D = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%s_sess%d.mat',parcelType,betaChoice,ss))); % loads region data (D)
            T = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            
            runs=1:numruns_task_sess;
            conditionVec = kron(ones(numel(runs),1),[1:12]');

            indx_train = 1:6;
            indx_untrain = 7:12;
            
            for s=1:length(sn)
                SI = load(fullfile(glmFoSExDir{ss}, subj_name{sn(s)}, 'SPM_info.mat'));   % load subject's trial structure
                for r=roi
                    for exe=1:2;    % FoSEx
                        D_exe = getrow(D,D.FoSEx==exe & D.region==r & D.SN==sn(s));
                        % if empty - skip (Buckner network 5 for some subjects)
                        if size(D_exe.IPM,1)~=0
                            switch (betaChoice)
                                case 'uniPW'
                                    beta = T.betaUW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe,:);
                                case 'multiPW'
                                    beta = T.betaW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe,:);
                                case 'raw'
                                    beta = T.betaRAW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe,:);
                            end
                            
                            clear C;
                            for d = 1:6 %sequences
                                C.beta_seq_train(d,:)=mean(beta(conditionVec==indx_train(d),:),1);  % beta values for each digit (avrg across blocks)
                                C.beta_seq_untrain(d,:)=mean(beta(conditionVec==indx_untrain(d),:),1);
                            end
                            
                            AllDist = ssqrt(rsa_squareRDM(D_exe.RDM));
                            %AllDist = ssqrt(rsa_squareRDM(D_exe.RDM(D_exe.region==r & D_exe.SN==sn(s),:)));
                            SeqTrain = triu(AllDist(indx_train,indx_train));
                            SeqUntrain = triu(AllDist(indx_untrain,indx_untrain));
                            SeqCross = triu(AllDist(indx_train,indx_untrain));
                            SeqTrainAll = SeqTrain(SeqTrain~=0);
                            SeqUntrainAll = SeqUntrain(SeqUntrain~=0);
                            SeqCrossAll = SeqCross(SeqCross~=0);
                            
                            switch (subjfig)
                                case 1
                                    figure(r)
                                    subplot(2,max(sessN),(exe-1)*max(sessN)+ss)
                                    imagesc(AllDist);
                                    drawline(6.5,'dir','vert');
                                    drawline(6.5,'dir','horz');
                                    colorbar; caxis([-0.02 0.04]);
                                    if ss == 4; set(gcf,'PaperPosition',[2 2 20 8],'Color',[1 1 1]); wysiwyg; end;
                                    title(sprintf('Subject %d session %d region %s exe %d',sn(s),ss,regname{r},exe));
                            end
                            
                            S.sn=sn(s);
                            S.roi=r;
                            S.regSide=D_exe.regSide;
                            S.regType=D_exe.regType;
                            S.dist_train=mean(SeqTrainAll);
                            S.dist_untrain=mean(SeqUntrainAll);
                            S.dist_cross=mean(SeqCrossAll);
                            S.beta_train=mean(mean(C.beta_seq_train));
                            S.beta_untrain=mean(mean(C.beta_seq_untrain));
                            S.psc_train=D.psc_train(D.SN==sn(s)&D.region==r&D.FoSEx==exe); % CHECK
                            S.psc_untrain=D.psc_untrain(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            S.act_vecLength_train=D.act_vecLength_train(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            S.act_vecLength_untrain=D.act_vecLength_untrain(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            %if exe==1
                            %    S.psc_train=D.psc_train_1st(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            %    S.psc_untrain=D.psc_untrain_1st(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            %else
                            %    S.psc_train=D.psc_train_2nd(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            %    S.psc_untrain=D.psc_untrain_2nd(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            %end
                            S.sessN=ss;
                            S.FoSEx=exe;
                            if S.sessN<4
                                S.speed=1;
                            else
                                S.speed=2;
                            end
                            Stats=addstruct(Stats,S);
                        end; % if no data
                    end; % FoSEx
                end; % roi
                fprintf('Done sess%d %s\n',ss,subj_name{s});
            end; % sn
        end

        % % save 
        save(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)),'-struct','Stats');

    case 'PLOT_beta'
       % repetition suppression in betas - for trained / untrained

       betaChoice = 'multiPW';
       roi=[1:8];
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
 
       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.beta_train,'split',S.FoSEx,'subset',S.roi==f,'style',styTrained_exe,'leg',{'1st','2nd'});
           ylim([0 0.07])
           if f==1
               ylabel('Betas trained');
               xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname{f}));
       end
       
       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.beta_untrain,'split',S.FoSEx,'subset',S.roi==f,'style',styUntrained_exe,'leg',{'1st','2nd'});
           ylim([0 0.07])
           if f==1
               ylabel('Betas untrained');
               xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname{f}));
       end
    case 'PLOT_psc'
       % repetition suppression in psc - for trained / untrained
       parcelType='Brodmann';  % 162tessels, Brodmann, cortex_buckner
       betaChoice = 'multiPW';
       roi=[1:8];
       hemi=1;
       vararginoptions(varargin,{'betaChoice','roi','parcelType','hemi'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
       switch parcelType
           case 'cortex_buckner'
               lab = {'Network-1','Network-2','Network-3','Network-4',...
                   'Network-5','Network-6','Network-7'};
           case 'Brodmann'
               lab=regname_cortex;
           case 'BG-striatum'
               lab=regname_BG;
       end
       for st=1:2 % sequence type
           figure
           if st==1
               style=styTrained_exe;
               var=S.psc_train;
           else
               style=styUntrained_exe;
               var=S.psc_untrain;
           end
           indx=1;
           for f = roi
               subplot(1,numel(roi),indx);
               plt.line([S.speed S.sessN],var,'split',S.FoSEx,'subset',S.regType==f&S.regSide==hemi,'style',style,'leg',{'1st','2nd'});
               drawline(0,'dir','horz');
               plt.match('y');
               if f==1
                   ylabel(sprintf('Psc %s',SeqType{st}));
                   xlabel('Session');
               else
                   ylabel('');
               end
               title(sprintf('%s',lab{f}));
               indx=indx+1;
           end
       end
    case 'PLOT_psc_sess1'
        % plot altogether trained and untrained
       betaChoice = 'multiPW';
       roi=[1:5,7,8];
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
           
           SS=[];
           SS.psc=[S.psc_train;S.psc_untrain];
           SS.seqType=[ones(size(S.psc_train));ones(size(S.psc_untrain))*2];
           SS.roi=[S.roi;S.roi];
           SS.sessN=[S.sessN;S.sessN];
           SS.FoSEx=[S.FoSEx;S.FoSEx];
           SS.sn=[S.sn;S.sn];
           
           figure
           plt.bar(SS.roi,SS.psc,'split',SS.FoSEx,'subset',SS.sessN==1&ismember(SS.roi,roi),'style',stySeqType_exe,'leg',{'1st','2nd'});
           ylabel('Activation in psc');
           set(gca,'XTick',[2 5 8 11 14.5 18 21],'XTickLabel',regname_cortex(roi));
           xlabel('ROI');
           
           for r=roi
               fprintf('%s \t',regname_cortex{r});
               ttestDirect(SS.psc,[SS.FoSEx SS.sn],2,'paired','subset',SS.roi==r&SS.sessN==1);
           end
          
    case 'PLOT_dist'
       % plot Mahalanobis distances 1st vs. 2nd execution
       % trained, untrained and between seq type
       parcelType='Brodmann';  % 162tessels, Brodmann, cortex_buckner
       betaChoice = 'multiPW';
       roi=[1:8];
       hemi=1;
       vararginoptions(varargin,{'betaChoice','roi','parcelType','hemi'});
       
       if strcmp(parcelType,'cortex_buckner')
           S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
       else
           S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
       end
       figure % trained
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.dist_train,'split',S.FoSEx,'subset',S.regType==f&S.regSide==hemi,'style',styTrained_exe,'leg',{'1st','2nd'});
           plt.match('y');
           if f==1
               ylabel('Distances trained'); xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname{f}))
       end
      
       figure % untrained
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.dist_untrain,'split',S.FoSEx,'subset',S.regType==f&S.regSide==hemi,'style',styUntrained_exe,'leg',{'1st','2nd'});
           plt.match('y');
           if f==1
               ylabel('Distances untrained'); xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname{f}))
       end
       
       figure % seqType
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.dist_cross,'split',S.FoSEx,'subset',S.regType==f&S.regSide==hemi,'style',stySeqType_exe,'leg',{'1st','2nd'});
           plt.match('y');
           if f==1
               ylabel('Distance between seq sets'); xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname{f}))
       end
    case 'PLOT_repsup_roi'
        % both activation and distance per ROI
        betaChoice = 'multiPW';
        roi=1;
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        for r=1:numel(roi)
            figure
            subplot(3,4,[1 2])
            plt.line([S.speed S.sessN],S.psc_train,'split',S.FoSEx,'subset',S.roi==roi(r),'style',styTrained_exe,'leg',{'1st','2nd'}); ylabel(''); title(sprintf('%s trained act',regname{roi(r)}));
            subplot(3,4,[3 4])
            plt.line([S.speed S.sessN],S.psc_untrain,'split',S.FoSEx,'subset',S.roi==roi(r),'style',styUntrained_exe,'leg',{'1st','2nd'}); title('untrained act'); ylabel('');
            
            subplot(3,4,[5 6])
            plt.line([S.speed S.sessN],S.dist_train,'split',S.FoSEx,'subset',S.roi==roi(r),'style',styTrained_exe,'leg',{'1st','2nd'}); title('trained dist'); ylabel('');
            subplot(3,4,[7 8])
            plt.line([S.speed S.sessN],S.dist_untrain,'split',S.FoSEx,'subset',S.roi==roi(r),'style',styUntrained_exe,'leg',{'1st','2nd'});title('untrained dist'); ylabel('');
            subplot(3,4,[10 11])
            plt.line([S.speed S.sessN],S.dist_cross,'split',S.FoSEx,'subset',S.roi==roi(r),'style',stySeqType_exe,'leg',{'1st','2nd'}); title('cross dist'); ylabel('');
            
        end
    case 'PLOT_psc_exe'
        % plot per region trained and untrained - separately 1st / 2nd execution
        % comparison of 1st: trained vs. untrained and 2nd: T vs. UT
       betaChoice = 'multiPW';
       roi=[1:8];
       exe=[1:2];
       vararginoptions(varargin,{'betaChoice','roi','exe'});
       
       exe_label={'1st','2nd'};
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
       
       for r=1:numel(roi)
           figure(r)
           for e = exe
               subplot(1,2,e);
               plt.line([[S.speed;S.speed] [S.sessN;S.sessN]],[S.psc_train;S.psc_untrain],'split',[ones(size(S.sn));ones(size(S.sn))*2],'subset',[S.FoSEx;S.FoSEx]==e & [S.roi;S.roi]==roi(r),'style',stySeq,'leg',{'trained','untrained'});
               plt.match('y');
               if e==1
                   ylabel('Percent signal change');
                   xlabel('Session');
               else
                   ylabel('')
               end
               title(sprintf('%s-%s',regname{r},exe_label{e}))
           end
       end    
    case 'PLOT_dist_exe'
       betaChoice = 'multiPW';
       roi=[1:8];
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
               title(sprintf('%s-%s',regname{r},exe_label{e}))
           end
       end
    case 'STATS_psc'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        seqType={'trained','untrained'};
        
        D = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        
        T.psc = [D.psc_train; D.psc_untrain];
        T.seqType = [ones(size(D.psc_train));ones(size(D.psc_train))*2];
        T.sn = [D.sn; D.sn];
        T.roi = [D.roi; D.roi];
        T.sessN = [D.sessN; D.sessN];
        T.FoSEx = [D.FoSEx; D.FoSEx];
        
        
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n PSC in %s \n',regname{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.FoSEx T.seqType],{'session','repsup','seqType'},'subset',T.roi==r);
        end
    case 'STATS_dist'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        seqType={'trained','untrained'};
        
        D = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        
        T.dist = [D.dist_train; D.dist_untrain];
        T.seqType = [ones(size(D.psc_train));ones(size(D.psc_train))*2];
        T.sn = [D.sn; D.sn];
        T.roi = [D.roi; D.roi];
        T.sessN = [D.sessN; D.sessN];
        T.FoSEx = [D.FoSEx; D.FoSEx];
           % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.dist,T.sn,'within',[T.sessN T.FoSEx T.seqType],{'session','repsup','seqType'},'subset',T.roi==r);
        end
        for r=roi
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS ratio in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.dist,[T.seq T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end        
    case 'STATS_dist_seqType'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        
           % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.dist_cross,T.sn,'within',[T.sessN T.FoSEx],{'session','repsup'},'subset',T.roi==r);
        end
        
    case 'SEARCH_run_LDC' % ------------------------- run searchlight 1st / 2nd --------------------------
        sessN  = [1:4];
        sn     = [4:9,11:25];
        repInd = [1,2]; % do searchlight for first or second execution
        vararginoptions(varargin,{'sn','sessN','repInd'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)

        for s=sn
            % load subject searchlight definition
            L = load(fullfile(anatomicalDir,subj_name{s},sprintf('%s_searchlight_120.mat',subj_name{s})));
            for ss=sessN
                % go to subject's glm directory
                cd(fullfile(glmFoSExDir{ss},subj_name{s}));
                for exe = repInd
                    % check first if the searchlight map already exists
                    if exist(sprintf('%s_sess%d_exe%d_LDC.nii',subj_name{s},ss,exe));
                        fprintf('Searchlight map: %s sess%d exe%d already exists - skipping\n',subj_name{s},ss,exe);
                    else
                        D=load('SPM_info.mat');           
                        % initialise vectors
                        conditionVec  = zeros(size(D.run));    % 12 sequences
                        partition     = zeros(size(D.run));
                        
                        % only fill the ones of the right repetition
                        conditionVec(D.FoSEx==exe,:)  = D.seqNumb(D.FoSEx==exe);
                        partition(D.FoSEx==exe,:)     = D.run(D.FoSEx==exe);
                        
                        % load SPM file
                        load SPM;
                        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                        
                        name = sprintf('%s_sess%d_exe%d',subj_name{s},ss,exe);
                        % run the searchlight
                        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name,'idealBlock',block);
                        fprintf('Done: %s sess-%d exe-%d \n\n\n',subj_name{s},ss,exe);
                    end
                end; % repetition
            end % session
        end % subject
        cd(cwd);
    case 'SEARCH_check'
        % check if file exists       
        sn=[4:9,11:25];
        sessN=4;
        exe=[1,2];      
        vararginoptions(varargin,{'sessN'});
        
        for ss=sessN
            for s=sn
                for e=exe
                    LDC_file = fullfile(glmFoSExDir{ss},subj_name{s},sprintf('%s_sess%d_exe%d_LDC.nii',subj_name{s},ss,e)); % searchlight nifti
                    if exist(LDC_file)
                        fprintf('Exists - %s sess-%d - exe-%d\n',subj_name{s},ss,e);
                    else
                        fprintf('Missing - %s sess-%d - exe-%d\t!!!\n',subj_name{s},ss,e);
                    end
                end
            end
        end
    case 'SEARCH_map'
        sn  = 4;
        sessN = 4;
        exe = [1,2];
        vararginoptions(varargin,{'sn','sessN','exe'});

        cWD = cd;
        for ss=sessN
            for s = sn
                for e=exe
                    % Load subject surface searchlight results (1 vol per paired conds)
                    LDC_file            = fullfile(glmFoSExDir{ss},subj_name{s},sprintf('%s_sess%d_exe%d_LDC.nii',subj_name{s},ss,e)); % searchlight nifti
                    [subjDir,fname,ext] = fileparts(LDC_file);
                    cd(subjDir);
                    
                    vol     = spm_vol([fname ext]);
                    vdat    = spm_read_vols(vol); % is searchlight data
                    
                    % average across all paired dists
                    Y.LDC   = nanmean(vdat,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat;
                    Y.fname   = sprintf('%s_sess%d_exe%d_dist.nii',subj_name{s},ss,e);
                    
                    % separate trained / untrained dists
                    indMat = indicatorMatrix('allpairs',[1:12]);
                    trainInd = [];
                    untrainInd = [];
                    crossInd = [];
                    for i=1:length(indMat)
                        if sum(any(indMat(i,num_train),1))==2
                            trainInd=[trainInd;i];
                        elseif sum(any(indMat(i,num_untrain),1))==2
                            untrainInd=[untrainInd;i];
                        elseif sum(any(indMat(i,num_train),1))==1 & sum(any(indMat(i,num_untrain),1))==1
                            crossInd=[crossInd;i];
                        else
                        end
                    end
                    
                    trainDist = vdat(:,:,:,trainInd);
                    untrainDist = vdat(:,:,:,untrainInd);
                    crossDist = vdat(:,:,:,crossInd);
                    
                    T = Y;
                    U = Y;
                    Z = Y;
                    T.LDC = nanmean(trainDist,4);
                    U.LDC = nanmean(untrainDist,4);
                    Z.LDC = nanmean(crossDist,4);
                    T.fname = sprintf('%s_sess%d_exe%d_dist_trained.nii',subj_name{s},ss,e);
                    U.fname = sprintf('%s_sess%d_exe%d_dist_untrained.nii',subj_name{s},ss,e);
                    Z.fname = sprintf('%s_sess%d_exe%d_dist_cross.nii',subj_name{s},ss,e);
                    
                    % save outputs
                    spm_write_vol(Y,Y.LDC);
                    spm_write_vol(T,T.LDC);
                    spm_write_vol(U,U.LDC);
                    spm_write_vol(Z,Z.LDC);
                    fprintf('Done %s_sess%d-exe%d \n',subj_name{s},ss,e);
                    
                    clear vol vdat LDC Y T U Z
                end
            end
        end
        cd(cWD);  % return to working directory
    case 'SEARCH_dist_surf'      
        % call dist function for mapping subject contrasts to surfaces
        sn  = [4:9,11:15,18:23,25];
        sessN = 4;   
        fileList = {'exe1_dist','exe2_dist','exe1_dist_trained','exe2_dist_trained','exe1_dist_untrained','exe2_dist_untrained','exe1_dist_cross','exe2_dist_cross'}; % for dist or repsup
        glmDir = glmFoSExDir;
        outname = 'dist_repsup'; % dist or dist_repsup
    
        vararginoptions(varargin,{'sn','sessN','fileList','glmDir','outname'});
        
        sml1_imana_dist('SEARCH_dist_surf','sn',sn,'sessN',sessN,'fileList',fileList,'glmDir',glmDir,'outname',outname);
    case 'SEARCH_group_make'                                                % STEP 4.5   :  Make group metric files by condensing subjec contrast metric files
        % call dist function for making the group contrast       
        sessN   = 4;
        sn      = [4:9,11:15,18:23,25];
        name    = 'RS_dist';
        OUTname = {'exe1_dist','exe2_dist','exe1_dist_trained','exe2_dist_trained','exe1_dist_untrained','exe2_dist_untrained','exe1_dist_cross','exe2_dist_cross'};
        inputcol = [1:length(OUTname)];
        replaceNaN = ones(size(OUTname));
        vararginoptions(varargin,{'sessN','sn','name','OUTname'});
        sml1_imana_dist('SEARCH_group_make','sn',sn,'sessN',sessN,'name',name,'OUTname',OUTname,'inputcol',inputcol,'replaceNaN',replaceNaN);
    case 'SEARCH_group_cSPM'                                                % STEP 4.6   :  Generate a statistical surface map (onesample_t test) from smoothed group metric files. Also avgs. distances across subjs.
        % call dist function for calculating the overall group mean / T      
        sessN   = 4;
        sn      = [4:9,11:15,18:23,25];
        name    = 'RS_dist';
        SPMname = {'exe1_dist','exe2_dist','exe1_dist_trained','exe2_dist_trained','exe1_dist_untrained','exe2_dist_untrained','exe1_dist_cross','exe2_dist_cross'};
        sqrtTransform = ones(size(SPMname)); % Should you take ssqrt before submitting?
        % Yes, b/c used rsa.distanceLDC to
        % calculate distances. This function
        % returns squared cv mahalanobis distance.
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname'});
        sml1_imana_dist('SEARCH_group_cSPM','sn',sn,'sessN',sessN,'SPMname',SPMname,'sqrtTransform',sqrtTransform,'name',name);
    case 'SEARCH_group_smooth'
        sessN=1;
        vararginoptions(varargin,{'sessN'});
        
        for h=1:2        
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            %----define name of coord and topology
            coordfile=[hem{h} '.WHITE.coord'];
            topofile=[hem{h} '.CLOSED.topo'];
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for i=1:length(sessN);
                filename=[surfaceGroupDir filesep hem{h} '.summary_RS_dist_sess' num2str(sessN) '.metric']; % unsmoothed
                sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
            end;
        end;
    case 'SEARCH_RS_RGB'
        % plotting RS on surface as absolute 
        % subtracting 1st - 2nd distance
        sessN=4;
        smooth=1;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for h=1:2   
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            if smooth == 1
                C=caret_load(['s' hem{h} sprintf('.summary_RS_dist_sess%d.metric',sessN)]);
            else
                C=caret_load([hem{h} sprintf('.summary_RS_dist_sess%d.metric',sessN)]);
            end
            flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
            topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
            
            low_1=0.005;
            RGBdata(:,1)=C.data(:,1)-C.data(:,2); % Red: trained 1st-2nd BOLD signal
            RGBdata(:,3)=C.data(:,3)-C.data(:,4); % Blue: untrained 1st-2nd BOLD signal
            RGBdata(RGBdata(:,1)<low_1,1)=0;
            RGBdata(RGBdata(:,3)<low_1,3)=0;

            sc=[low_1 0.03;low_1 0.03;low_1 0.03];  % scaling
            
            name={sprintf('RS_absolute_dist_sess%d',sessN)};
            %name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN),sprintf('dist_cross_sess%d',sessN)};

            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
            
            if smooth == 1
                caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.RS_dist_sess%d.RGB_paint',sessN)],C);
            else
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.RS_dist_sess%d.RGB_paint',sessN)],C);
            end
        end;
        
        
    case 'corr_vertex'
        sessN=4;
        sn = [4:14,17:22,24];
        structure = 'cortex'; % cortex or striatum
        vararginoptions(varargin,{'sessN','sn','structure'});
        for ss=sessN
            for h=1:2 % hemisphere
                for st=1:2
                    if st==1
                        lab = {'TrainSeq','trained'};
                    else
                        lab = {'UntrainSeq','untrained'};
                    end
                    switch structure
                        case 'cortex'
                            surfaceGDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
                            p1 = caret_load(fullfile(surfaceGDir,sprintf('%s.RS_%s_1st_sess%d.metric',hem{h},lab{1},ss)));
                            p2 = caret_load(fullfile(surfaceGDir,sprintf('%s.RS_%s_2nd_sess%d.metric',hem{h},lab{1},ss)));
                            d1 = caret_load(fullfile(surfaceGDir,sprintf('%s.exe1_dist_%s_sess%d.metric',hem{h},lab{2},ss)));
                            d2 = caret_load(fullfile(surfaceGDir,sprintf('%s.exe2_dist_%s_sess%d.metric',hem{h},lab{2},ss)));
                        case 'striatum'
                            surfaceGDir = fullfile(BGDir,'surface','group',hemName{h});
                            p1 = caret_load(fullfile(surfaceGDir,sprintf('%s.psc_RS_%s_1st_sess%d.metric',hem{h},lab{1},ss)));
                            p2 = caret_load(fullfile(surfaceGDir,sprintf('%s.psc_RS_%s_2nd_sess%d.metric',hem{h},lab{1},ss)));
                            d1 = caret_load(fullfile(surfaceGDir,sprintf('%s.dist_RS_%s_exe1_sess%d.metric',hem{h},lab{2},ss)));
                            d2 = caret_load(fullfile(surfaceGDir,sprintf('%s.dist_RS_%s_exe2_sess%d.metric',hem{h},lab{2},ss)));
                    end
                    nVert = size(p1.data,1);
                    for v=1:nVert
                        % difference in psc (dP) and dist (dD) for exe 1-2 at each vertex
                        dP = p1.data(v,sn)-p2.data(v,sn);
                        dD = d1.data(v,:)-d2.data(v,:);
                        Corr(v,:) = corr(dP',dD');
                    end
                    colName = sprintf('corr_%s_sess-%d',lab{2},ss);
                    C = caret_struct('metric','data',Corr,'column_name',{colName});
                    C.index=p1.index;
                    caret_save(fullfile(surfaceGDir, sprintf('%s.RS_correlation_%s_sess%d.metric',hem{h},lab{2},ss)),C);
                end
            end
        end
    case 'corr_smooth'
        sessN=4;
        structure = 'cortex';
        nIter = 20; % 20 for cortex, 3 for striatum
        vararginoptions(varargin,{'sessN','structure','nIter'});
        seqName = {'trained','untrained'};
        for h=1:2       
            switch structure
                case 'cortex'
                    surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                    coordfile=[hem{h} '.WHITE.coord'];
                    topofile=[hem{h} '.CLOSED.topo'];
                case 'striatum'
                    surfaceGroupDir = fullfile(BGDir,'surface','group',hemName{h});
                    coordfile=[hem{h} '.striatum.coord.gii'];
                    topofile=[hem{h} '.striatum.topo.gii'];
            end
            cd(surfaceGroupDir)
            %----define name of coord and topology
           
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for st=1:2
                for i=1:length(sessN);
                    filename=[surfaceGroupDir filesep hem{h} '.RS_correlation_' seqName{st} '_sess' num2str(sessN(i)) '.metric']; % unsmoothed
                    sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations',nIter);
                end;
            end
        end;
        
    case 'CALC_corrDist' %--------------------------- CORRELATION DISTANCE -------------------------------
        reg = [1:14];
        sn  = [4:9,11:31];
        sessN = [1:4];
        parcelType='Brodmann';  % 162tessels, Brodmann, cortex_buckner
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        vararginoptions(varargin,{'sn','reg','sessN','subtract_mean','parcelType'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
            
            D   = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            
            fName = fieldnames(D);
            if any(strcmp('psc_train',fName));
                D = rmfield(D,{'psc_train','psc_untrain'});
            end
            for roi = reg;
                for s = sn;
                    SI = load(fullfile(glmFoSExDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                    
                    for exe = 1:2
                        t   = getrow(D,D.region==roi & D.SN==s);
                        SI_exe = getrow(SI,SI.FoSEx==exe);
                        % check if any data
                        if ~isempty(t.SN) & sum(sum(isnan(t.betaRAW{1})))==0
                            data = t.betaW{1}(SI.FoSEx==exe,:);
                            numRuns    = 1:numruns_task_sess;
                            numConds   = num_seq;
                            conds   = repmat([numConds],1,length(numRuns))';
                            runNums = kron([numRuns],ones(1,length(numConds)))';
                            
                            if subtract_mean
                                for r=numRuns
                                    data(runNums==r,:) = bsxfun(@minus,data(runNums==r,:),mean(data(runNums==r,:)));
                                end
                            else
                                data=data;
                            end
                            
                            % non-crossvalidated
                            Z=pcm_indicatorMatrix('identity',SI_exe.seqNumb);
                            b = pinv(Z)*data;           % Estimate mean activities
                            b(1:6,:)  = bsxfun(@minus,b(1:6,:) ,mean(b(1:6,:))); % Subtract mean per condition - trained
                            b(7:12,:) = bsxfun(@minus,b(7:12,:),mean(b(7:12,:))); % untrained
                            G=cov(b');
                            C_ncross=corr_crossval(G,'reg','minvalue');
                            C_ncross=rsa_squareRDM(C_ncross);
                            % average trained dist
                            S.corrDist_nonCrossval(1,:) = sum(sum(C_ncross(1:6,1:6)))/(6*5);
                            % average untrained dist
                            S.corrDist_nonCrossval(2,:) = sum(sum(C_ncross(7:12,7:12)))/(6*5);
                            
                            % crossval second moment matrix
                            [G,Sig]     = pcm_estGCrossval(data(1:(12*num_run),:),SI_exe.run,SI_exe.seqNumb);
                            C=corr_crossval(G,'reg','minvalue');
                            C=rsa_squareRDM(C);
                            dist = rsa_squareRDM(rsa.distanceLDC(data,SI_exe.run,SI_exe.seqNumb));
                            % average trained dist
                            S.corrDist(1,:) = sum(sum(C(1:6,1:6)))/(6*5);
                            % average untrained dist
                            S.corrDist(2,:) = sum(sum(C(7:12,7:12)))/(6*5);
                            S.dist(1,:)     = sum(sum(triu(dist(1:6,1:6))))/(6*5/2);
                            S.dist(2,:)     = sum(sum(triu(dist(7:12,7:12))))/(6*5/2);
                            S.seqType=[1;2]; % trained, untrained
                            S.sessN=[ss;ss];
                            S.roi=[roi;roi];
                            S.regType=[t.regType;t.regType];
                            S.regSide=[t.regSide;t.regSide];
                            S.sn=[s;s];
                            S.FoSEx=[exe;exe];
                            SAll=addstruct(SAll,S);
                            
                            
                            ST.corr_seqType = sum(sum(C(7:12,1:6)))/(6*5);
                            ST.sessN=ss;
                            ST.roi=roi;
                            ST.regType=t.regType;
                            ST.regSide=t.regSide;
                            ST.sn=s;
                            ST.FoSEx=exe;
                            STAll=addstruct(STAll,ST);
                        end
                    end;
                end;
            end;
             fprintf('Done ses%d\n',ss);
        end;
        save(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)),'-struct','SAll');
        save(fullfile(repSupDir,sprintf('corrDist_seqType_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)),'-struct','STAll');
        
    case 'PLOT_corrDist'
        roi=[1:8];
        subtract_mean=0;
        distType='crossval';
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','subtract_mean','distType','parcelType','hemi'});
        
        S=load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)));
        
        
        switch(distType)
            case 'crossval'
                var=S.corrDist;
            case 'non-crossval'
                var=S.corrDist_nonCrossval;
        end
        
        for seq=1:2
            indx=1;
            figure
            for r=roi
                subplot(1,numel(roi),indx)
                if seq==1
                    style=styTrained_exe;
                else
                    style=styUntrained_exe;
                end
                
               % plt.line([S.sessN>3 S.sessN],var,'style',style,'split',S.FoSEx,'subset',S.seqType==seq&S.roi==r,'leg',{'1st','2nd'});
                plt.line([S.sessN>3 S.sessN],S.dist,'style',style,'split',S.FoSEx,'subset',S.seqType==seq&S.regType==r&S.regSide==hemi,'leg',{'1st','2nd'});
                drawline(0,'dir','horz');
                %plt.match('y');
                title(sprintf('%s',regname{r}));
                if r==1
                    xlabel('Session'); ylabel(sprintf('Corr distance %s',SeqType{seq}));
                else
                    ylabel('');
                end
                indx=indx+1;
            end
        end
    case 'PLOT_corrDist_sess1'    
        roi=[1:5,7,8];
        subtract_mean=0;
        vararginoptions(varargin,{'roi','subtract_mean','distType'});
        
        SS=load(fullfile(repSupDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        
        figure
        plt.bar(SS.roi,SS.corrDist,'split',SS.FoSEx,'subset',SS.sessN==1&ismember(SS.roi,roi),'style',stySeqType_exe,'leg',{'1st','2nd'});
        ylabel('Cosine distance');
        set(gca,'XTick',[2 5 8 11 14.5 18 21],'XTickLabel',regname_cortex(roi));
        xlabel('ROI');
        
        for r=roi
            fprintf('%s \t',regname_cortex{r});
            ttestDirect(SS.corrDist,[SS.FoSEx SS.sn],2,'paired','subset',SS.roi==r&SS.sessN==1);
        end
     
    case 'PLOT_corrDist_exe'
        roi=[1:8];
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
            title(sprintf('%s-%s',regname{r},Exe{exe}));
            if exe==1
            xlabel('Session'); ylabel('Corr distance');
            else
                ylabel('');
            end
            end
        end
    case 'PLOT_corrDist_seqType'
        roi=[1:8];
        hemi=1;
        subtract_mean=0;
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','subtract_mean','hemi'});
        
        S=load(fullfile(repSupDir,sprintf('corrDist_seqType_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)));
       
        Exe={'1st exe','2nd exe'};
        figure;
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.line([S.sessN>3 S.sessN],S.corr_seqType,'style',stySeqType_exe,'split',S.FoSEx,'subset',S.regType==r&S.regSide==hemi,'leg',{'1st','2nd'});
            drawline(0,'dir','horz');
            %plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                xlabel('Session'); ylabel('Corr distance between seqType');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'STATS_corrDist'
        roi=[1:8];
        sessN=[1:4];
        subtract_mean=0;
        vararginoptions(varargin,{'betaChoice','roi','sessN','subtract_mean'});
        
        T=load(fullfile(repSupDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType T.FoSEx],{'session','seqType','repsup'},'subset',T.roi==r);
        end
        
        for ss=sessN
            for r=roi
                fprintf('\n Distance in %s session%d - larger than 0?',regname{r},ss);
                fprintf('\n trained \n');
                ttestDirect(T.corrDist,[T.sn],2,'onesample','subset',T.roi==r&T.seqType==1&T.sessN==ss);
                fprintf('\n untrained \n');
                ttestDirect(T.corrDist,[T.sn],2,'onesample','subset',T.roi==r&T.seqType==2&T.sessN==ss);
                fprintf('\n 1st>2nd \n')
                ttestDirect(T.corrDist,[T.FoSEx,T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);
            end
        end
        
    case 'STATS_corrDist_seqType'
        roi=[1:8];
        sessN=[1:4];
        subtract_mean=0;
        vararginoptions(varargin,{'betaChoice','roi','sessN','subtract_mean'});
        
        S=load(fullfile(repSupDir,sprintf('corrDist_seqType_RS_meanSubtract%d.mat',subtract_mean)));
        
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.corr_seqType,T.sn,'within',[T.sessN T.FoSEx],{'session','repsup'},'subset',T.roi==r);
        end
  
    case 'CALC_trainDur_seqSpec' %------------------- TRAINING DURATION (RATHER THAN SEQ-TYPE)-------------
        % label trials / sequences with respect to how many trials
        % performed before scanning
        % number of trials PER SEQUENCE
        betaChoice = 'multiPW';
        
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));

        % new structure - absolute values for 1st / 2nd
        T.sn=[S.sn;S.sn];
        T.roi=[S.roi;S.roi];
        T.seqType=[ones(size(S.roi));ones(size(S.roi))*2];
        T.dist=[S.dist_train;S.dist_untrain];
        T.psc=[S.psc_train;S.psc_untrain];
        T.beta=[S.beta_train;S.beta_untrain];
        T.FoSEx=[S.FoSEx;S.FoSEx];
        T.speed=[S.speed;S.speed];
        T.sessN=[S.sessN;S.sessN];
        % add training duration
        T.trainNum=zeros(size(T.sessN));
        T.trainNum(T.sessN==1)=0;
        T.trainNum(T.sessN==2&T.seqType==1)=162;
        T.trainNum(T.sessN==2&T.seqType==2)=48;
        T.trainNum(T.sessN==3&T.seqType==1)=560;
        T.trainNum(T.sessN==3&T.seqType==2)=96;
        T.trainNum(T.sessN==4&T.seqType==1)=608;
        T.trainNum(T.sessN==4&T.seqType==2)=144;
        
        
        % new structure - percent 2nd/1st
        A = getrow(T,T.FoSEx==1);
        
        P.sn=A.sn;
        P.roi=A.roi;
        P.seqType=A.seqType;
        P.psc_beta=[T.beta(T.FoSEx==2)./T.beta(T.FoSEx==1)];
        P.psc_psc =[T.psc(T.FoSEx==2)./T.psc(T.FoSEx==1)];
        P.subtr_dist=[T.dist(T.FoSEx==1)-T.dist(T.FoSEx==2)];
        P.sessN=A.sessN;
        P.speed=A.speed;
        P.trainNum=A.trainNum;
        
        % save structures
        save(fullfile(repSupDir,sprintf('RepSup_trainDur_SEQ_SPEC.mat')),'-struct','T');
        save(fullfile(repSupDir,sprintf('RepSup_trainDur_SEQ_SPEC_percent.mat')),'-struct','P');
    case 'CALC_trainDur_seqGeneral'
        % label trials / sequences with respect to how many trials
        % performed before scanning
        % number of trials IN GENERAL (not per sequence)
        betaChoice = 'multiPW';
        
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));

        % new structure - absolute values for 1st / 2nd
        T.sn=[S.sn;S.sn];
        T.roi=[S.roi;S.roi];
        T.seqType=[ones(size(S.roi));ones(size(S.roi))*2];
        T.dist=[S.dist_train;S.dist_untrain];
        T.psc=[S.psc_train;S.psc_untrain];
        T.beta=[S.beta_train;S.beta_untrain];
        T.FoSEx=[S.FoSEx;S.FoSEx];
        T.speed=[S.speed;S.speed];
        T.sessN=[S.sessN;S.sessN];
        % add training duration
        T.trainNum=zeros(size(T.sessN));
        T.trainNum(T.sessN==1)=528;
        T.trainNum(T.sessN==2)=2160;
        T.trainNum(T.sessN==3)=6288;
        T.trainNum(T.sessN==4)=6864;
        
        
        % new structure - percent 2nd/1st
        A = getrow(T,T.FoSEx==1);
        
        P.sn=A.sn;
        P.roi=A.roi;
        P.seqType=A.seqType;
        P.psc_beta=[T.beta(T.FoSEx==2)./T.beta(T.FoSEx==1)];
        P.psc_psc =[T.psc(T.FoSEx==2)./T.psc(T.FoSEx==1)];
        P.subtr_dist=[T.dist(T.FoSEx==1)-T.dist(T.FoSEx==2)];
        P.sessN=A.sessN;
        P.speed=A.speed;
        P.trainNum=A.trainNum;
        
        % save structures
        save(fullfile(repSupDir,sprintf('RepSup_trainDur_SEQ_GEN.mat')),'-struct','T');
        save(fullfile(repSupDir,sprintf('RepSup_trainDur_SEQ_GEN_percent.mat')),'-struct','P');
    case 'CALC_trainDur_seqType'
        % label trials / sequences with respect to how many trials
        % performed before scanning
        % number of trials IN GENERAL (not per sequence)
        betaChoice = 'multiPW';
        
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));

        % new structure - absolute values for 1st / 2nd
        T.sn=[S.sn;S.sn];
        T.roi=[S.roi;S.roi];
        T.seqType=[ones(size(S.roi));ones(size(S.roi))*2];
        T.dist=[S.dist_train;S.dist_untrain];
        T.psc=[S.psc_train;S.psc_untrain];
        T.beta=[S.beta_train;S.beta_untrain];
        T.FoSEx=[S.FoSEx;S.FoSEx];
        T.speed=[S.speed;S.speed];
        T.sessN=[S.sessN;S.sessN];
        % add training duration
        T.trainNum=zeros(size(T.sessN));
        T.trainNum(T.sessN==1)=0;
        T.trainNum(T.sessN==2&T.seqType==1)=1260;
        T.trainNum(T.sessN==2&T.seqType==2)=288;
        T.trainNum(T.sessN==3&T.seqType==1)=3648;
        T.trainNum(T.sessN==3&T.seqType==2)=576;
        T.trainNum(T.sessN==4&T.seqType==1)=3936;
        T.trainNum(T.sessN==4&T.seqType==2)=864;
        
        
        % new structure - percent 2nd/1st
        A = getrow(T,T.FoSEx==1);
        
        P.sn=A.sn;
        P.roi=A.roi;
        P.seqType=A.seqType;
        P.psc_beta=[T.beta(T.FoSEx==2)./T.beta(T.FoSEx==1)];
        P.psc_psc =[T.psc(T.FoSEx==2)./T.psc(T.FoSEx==1)];
        P.subtr_dist=[T.dist(T.FoSEx==1)-T.dist(T.FoSEx==2)];
        P.sessN=A.sessN;
        P.speed=A.speed;
        P.trainNum=A.trainNum;
        
        % save structures
        save(fullfile(repSupDir,sprintf('RepSup_trainDur_SEQ_TYPE.mat')),'-struct','T');
        save(fullfile(repSupDir,sprintf('RepSup_trainDur_SEQ_TYPE_percent.mat')),'-struct','P');
    case 'PLOT_trainDur_sequence'
        roi=1:8;
        type='seqSpec'; % seqType
        % training duration specific to each sequence or sequence type
        vararginoptions(varargin,{'roi','type'});
        
        switch type
            case 'seqSpec'
                P = load(fullfile(repSupDir,'RepSup_trainDur_SEQ_SPEC_percent.mat'));
            case 'seqType'
                P = load(fullfile(repSupDir,'RepSup_trainDur_SEQ_TYPE_percent.mat'));
        end
        
        for r=roi
            figure(1)
            subplot(1,max(numel(roi)),r)
            plt.line(P.trainNum,P.psc_psc,'style',stySeq,'split',P.seqType,'subset',P.sessN<4&P.roi==roi(r),'leg',{'trained','untrained'},'leglocation','north');
            title(sprintf('%s',regname{r}));
            hold on; drawline(0,'dir','horz'); drawline(1,'dir','horz');
            ylim([-0.2 1.4]);
            if r==1
                ylabel('Activation 2nd/1st');
            else
                ylabel('');
            end
            
            figure(2)
            subplot(1,max(numel(roi)),r)
            plt.line(P.trainNum,P.subtr_dist,'style',stySeq,'split',P.seqType,'subset',P.sessN<4&P.roi==roi(r));
            hold on; drawline(0,'dir','horz'); drawline(1,'dir','horz');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Distance 1st-2nd');
            else
                ylabel('');
            end
        end 
    case 'PLOT_trainDur_general'
        roi=1:8;
        % training duration in general - number of sequence trials
        % performed (trained / untrained / random)
        vararginoptions(varargin,{'roi'});
        
        P = load(fullfile(repSupDir,'RepSup_trainDur_SEQ_GEN_percent.mat'));

        for r=roi
            figure(1)
            subplot(1,max(numel(roi)),r)
            plt.line(P.trainNum,P.psc_psc,'style',stySeq,'split',P.seqType,'subset',P.sessN<4&P.roi==roi(r));
            title(sprintf('Activation %s',regname{r}));
            hold on; drawline(0,'dir','horz'); drawline(1,'dir','horz');
            ylim([-0.2 1.4]);
            if r==1
                ylabel('Activation 2nd/1st');
            else
                ylabel('');
            end
            
            figure(2)
            subplot(1,max(numel(roi)),r)
            plt.line(P.trainNum,P.subtr_dist,'style',stySeq,'split',P.seqType,'subset',P.sessN<4&P.roi==roi(r));
            hold on; drawline(0,'dir','horz'); drawline(1,'dir','horz');
            title(sprintf('Distance %s',regname{r}));
            if r==1
                ylabel('Distance 1st-2nd');
            else
                ylabel('');
            end
        end 
    case 'PLOT_speed'
        % 3rd vs. 4th session  
        roi=1:8;
        betaChoice='multiPW';
        
        vararginoptions(varargin,{'roi','betaChoice'});
        
        P = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s.mat',betaChoice)));
        M = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s.mat',betaChoice)));
        C = load(fullfile(repSupDir,'RepSup_subtract_corrDist.mat'));
        
        indx=1;
        for r=roi
            figure(1)
            subplot(1,length(roi),indx);
            plt.line(P.sessN,1-P.psc_psc,'style',stySeq,'split',P.seq,'subset',P.sessN>2&P.roi==r,'leg',{'trained','untrained'}); 
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('Activation ratio 1-2nd/1st');
                xlabel('Session');
            else
                ylabel('');
            end
            
            figure(2)
            subplot(1,length(roi),indx);
            plt.line(M.sessN,M.subtr_dist,'style',stySeq,'split',M.seq,'subset',M.sessN>2&M.roi==r,'leg',{'trained','untrained'});
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('Mahalanobis distance subtraction 1st-2nd');
                xlabel('Session');
            else
                ylabel('');
            end
            
            figure(3)
            subplot(1,length(roi),indx);
            plt.line(C.sessN,C.corrDist,'style',stySeq,'split',C.seqType,'subset',C.sessN>2&C.roi==r,'leg',{'trained','untrained'});
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('Correlation distance subtraction 1st-2nd');
                xlabel('Session');
            else
                ylabel('');
            end
            indx=indx+1;
        end
           
    case 'CALC_subtract_actDist' %--------------------------- REPSUP AS SUBTRACTION: 1st - 2nd ---------------------
       % Grafton version
       parcelType='Brodmann';
       betaChoice = 'multiPW';
       vararginoptions(varargin,{'betaChoice','roi','parcelType'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
 
       A = getrow(S,S.FoSEx==1);
       
       P.subtr_beta = [S.beta_train(S.FoSEx==1)-S.beta_train(S.FoSEx==2);S.beta_untrain(S.FoSEx==1)-S.beta_untrain(S.FoSEx==2)];
       P.subtr_psc = [S.psc_train(S.FoSEx==1)-S.psc_train(S.FoSEx==2);S.psc_untrain(S.FoSEx==1)-S.psc_untrain(S.FoSEx==2)];
       P.subtr_dist = [S.dist_train(S.FoSEx==1)-S.dist_train(S.FoSEx==2);S.dist_untrain(S.FoSEx==1)-S.dist_untrain(S.FoSEx==2)];
       P.sn = [A.sn; A.sn];
       P.roi = [A.roi; A.roi];
       P.regType = [A.regType; A.regType];
       P.regSide = [A.regSide; A.regSide];
       P.seq = [ones(size(A.sn)); ones(size(A.sn))*2];  % 1 for trained, 2 for untrained
       P.sessN = [A.sessN; A.sessN];
       P.speed = [A.speed; A.speed];
       
       % % save
        save(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)),'-struct','P');    
    case 'PLOT_subtract_beta'
        betaChoice = 'multiPW';
        parcelType = 'Brodmann';
        roi=[1:8];
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));

       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.subtr_beta,'split',S.seq,'subset',S.regType==f&S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
           plt.match('y');
           drawline(0,'dir','horz');
           
          % ylim([-0.1 1.1]); hold on; drawline(1,'dir','horz'); drawline(0,'dir','horz');
           if f==1
               ylabel('Subtract 1st - 2nd betas');
               xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname{f}));
       end
    case 'PLOT_subtract_psc'
        betaChoice = 'multiPW';
        roi=[1:8];
        parcelType = 'Brodmann';
        hemi=1;
        vararginoptions(varargin,{'betaChoice','roi','hemi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));

       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.subtr_psc,'split',S.seq,'subset',S.regType==f&S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
           drawline(0,'dir','horz');
           plt.match('y');
          % ylim([-0.1 1.1]); hold on; drawline(1,'dir','horz'); drawline(0,'dir','horz');
           if f==1
               ylabel('Subtract 1st - 2nd psc')
               xlabel('Session')
           else
               ylabel('')
           end
           title(sprintf('%s',regname{f}))
       end
    case 'PLOT_subtract_psc_sessions'
        betaChoice = 'multiPW';
        sessN=[1:3];
        seqType='trained';
        roi=[1:5,7,8];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s.mat',betaChoice))); 
        
        T=getrow(T,ismember(T.roi,roi)&ismember(T.sessN,sessN)); % exclude V12

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
        plt.line(T.reg,T.subtr_psc,'split',T.sessN,'subset',T.seq==st,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        drawline(0,'dir','horz');
        title(sprintf('%s',seqType));
        ylabel('Repetition suppression of activation - subtraction');
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI'); 
    case 'PLOT_subtract_psc_seqType'
        betaChoice = 'multiPW';
        sessN=[1:3];
        roi=[1:5,7,8];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});

        T=load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s.mat',betaChoice))); 
        
        T=getrow(T,ismember(T.roi,roi)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        for ss=sessN
            figure
            plt.bar(T.reg,T.subtr_psc,'split',T.seq,'subset',T.sessN==ss,'leg',{'trained','untrained'},'leglocation','northeast','style',stySeq);
            drawline(0,'dir','horz');
            title(sprintf('sess-%d',ss));
            ylabel('Repetition suppression of activation - subtraction');
          %  set(gca,'XTickLabel',regname_cortex(reg));
            set(gca,'XTick',[2 5 8 11 14 18 21],'XTickLabel',regname_cortex(reg));
        end
    case 'PLOT_subtract_dist'
        betaChoice = 'multiPW';
        roi=[1:8];
        parcelType = 'Brodmann';
        hemi=1;
        vararginoptions(varargin,{'betaChoice','roi','hemi'});
     
        S = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));

       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.subtr_dist,'split',S.seq,'subset',S.regType==f&S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
          % ylim([-0.1 1.1]); hold on; drawline(1,'dir','horz'); drawline(0,'dir','horz');
           if f==1
               ylabel('Subtract 1st - 2nd dist')
           else
               ylabel('')
           end
           title(sprintf('%s',regname{f}))
       end 
    case 'CALC_subtract_distCross'
       betaChoice = 'multiPW';
       roi=[1:8];
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
 
       A = getrow(S,S.FoSEx==1);
       
       C.sn=A.sn;
       C.roi=A.roi;
       C.psc_distcross = [S.dist_cross(S.FoSEx==2)-S.dist_cross(S.FoSEx==1)];
       C.sessN=A.sessN;
       C.speed=A.speed;
       
       % % save
       save(fullfile(repSupDir,sprintf('RepSup_subtract_crossseq_%s.mat',betaChoice)),'-struct','C');
    case 'PLOT_subtract_distCross'
        % distance between sets of trained and untrained sequences
        betaChoice = 'multiPW';
        roi=[1:8];
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_subtract_crossseq_%s.mat',betaChoice)));
        
        figure
        for f = 1:numel(roi)
            subplot(1,numel(roi),f);
            plt.line([S.speed S.sessN],S.psc_distcross,'subset',S.roi==f);
            hold on; drawline(1,'dir','horz'); drawline(0,'dir','horz');
            if f==1
                ylabel('Cross-dist 1st-2nd')
            else
                ylabel('')
            end
            title(sprintf('%s',regname{f}))
        end         
    case 'CALC_subtract_corrDist'
        
       roi=[1:8];
       vararginoptions(varargin,{'roi'});
       
       S = load(fullfile(repSupDir,'corrDist_RS.mat')); 
 
       A = getrow(S,S.FoSEx==1);
       
       C.sn=A.sn;
       C.roi=A.roi;
       C.corrDist = [S.corrDist(S.FoSEx==1)-S.corrDist(S.FoSEx==2)];
       C.sessN=A.sessN;
       C.seqType=A.seqType;
       
       % % save
       save(fullfile(repSupDir,'RepSup_subtract_corrDist.mat'),'-struct','C');
    case 'PLOT_subtract_corrDist_sessions'
        sessN=[1:3];
        seqType='trained';
        regExcl=6;
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(repSupDir,'RepSup_subtract_corrDist.mat')); 
        
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
        plt.line(T.reg,T.corrDist,'split',T.sessN,'subset',T.seqType==st,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        drawline(0,'dir','horz');
        ylabel('Distance repetition suppression');
        set(gca,'XTickLabel',regname_cortex(reg));
        title(sprintf('%s',seqType));
        xlabel('ROI');  
    case 'PLOT_subtract_corrDist_seqType'
        sessN=[1:3];
        regExcl=6;
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        
        T=load(fullfile(repSupDir,'RepSup_subtract_corrDist.mat')); 
        
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        for ss=sessN
            figure
            plt.bar(T.reg,T.corrDist,'split',T.seqType,'subset',T.sessN==ss,'leg',{'trained','untrained'},'leglocation','northeast','style',stySeq);
            drawline(0,'dir','horz');
            ylabel('Distance repetition suppression');
            set(gca,'XTick',[2 5 8 11 14 18 21],'XTickLabel',regname_cortex(reg));
            title(sprintf('sess-%d',ss));
            xlabel('ROI');
        end
    case 'STATS_subtract_psc'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s.mat',betaChoice)));
        
        % main ANOVA - session x sequence
        for r=roi
            fprintf('\n PSC in %s \n',regname{r});
            anovaMixed(T.subtr_psc,T.sn,'within',[T.sessN T.seq],{'session','seqType'},'subset',T.roi==r);
        end
        
        % t-tests comparing sequence type per session
        for r=roi
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS subtraction in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.subtr_psc,[T.seq T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end
        
        % t-tests comparing across sessions - separately for T / UT seq
        for st=1:2
            fprintf('\n %s sequence\n',seqLabel{st});
            for sessTr=1:3
                for r=roi
                    fprintf('\n post-hoc t-test on comparison of sessions %d-%d in %s \n',sessTr,sessTr+1,regname{r});
                    ttestDirect(T.subtr_psc,[T.sessN T.sn],2,'paired','subset',T.roi==r&ismember(T.sessN,[sessTr,sessTr+1])&T.seq==st);
                end               
            end
        end
    case 'STATS_subtract_dist'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s.mat',betaChoice)));
        
        % main ANOVA - session x sequence
        for r=roi
            fprintf('\n PSC in %s \n',regname{r});
            anovaMixed(T.subtr_dist,T.sn,'within',[T.sessN T.seq],{'session','seqType'},'subset',T.roi==r);
        end
        
        % t-tests comparing sequence type per session
        for r=roi
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS subtraction in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.subtr_dist,[T.seq T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end
        
        % t-tests comparing across sessions - separately for T / UT seq
        for st=1:2
            fprintf('\n %s sequence\n',seqLabel{st});
            for sessTr=1:3
                for r=roi
                    fprintf('\n post-hoc t-test on comparison of sessions %d-%d in %s \n',sessTr,sessTr+1,regname{r});
                    ttestDirect(T.subtr_dist,[T.sessN T.sn],2,'paired','subset',T.roi==r&ismember(T.sessN,[sessTr,sessTr+1])&T.seq==st);
                end               
            end
        end
    case 'STATS_subtract_corrDist'
        roi=[1:8];
        sessN=[1:4];
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'roi','sessN'});
        
        T = load(fullfile(repSupDir,'RepSup_subtract_corrDist.mat'));
        
        % main ANOVA - session x sequence
        for r=roi
            fprintf('\n Corr dist in %s \n',regname{r});
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r);
        end
        
        % t-tests comparing sequence type per session
        for r=roi
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS distance subtraction in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.corrDist,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end
        
        % t-tests comparing across sessions - separately for T / UT seq
        for st=1:2
            fprintf('\n %s sequence\n',seqLabel{st});
            for sessTr=1:3
                for r=roi
                    fprintf('\n post-hoc t-test on comparison of sessions %d-%d in %s \n',sessTr,sessTr+1,regname{r});
                    ttestDirect(T.corrDist,[T.sessN T.sn],2,'paired','subset',T.roi==r&ismember(T.sessN,[sessTr,sessTr+1])&T.seqType==st);
                end               
            end
        end        
         
    case 'CALC_ratio_actDist' %--------------------------- REPSUP AS RATIO: 2nd/1st --------------------------
               
       betaChoice = 'multiPW';
       roi=[1:8];
       parcelType='Brodmann';
       vararginoptions(varargin,{'betaChoice','roi','parcelType'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
 
       A = getrow(S,S.FoSEx==1);
       
       P.psc_beta = [S.beta_train(S.FoSEx==2)./S.beta_train(S.FoSEx==1);S.beta_untrain(S.FoSEx==2)./S.beta_untrain(S.FoSEx==1)];
       P.psc_psc  = [S.psc_train(S.FoSEx==2)./S.psc_train(S.FoSEx==1);S.psc_untrain(S.FoSEx==2)./S.psc_untrain(S.FoSEx==1)];
       P.psc_dist = [S.dist_train(S.FoSEx==2)./S.dist_train(S.FoSEx==1);S.dist_untrain(S.FoSEx==2)./S.dist_untrain(S.FoSEx==1)];

       P.sn = [A.sn; A.sn];
       P.roi = [A.roi; A.roi];
       P.regType = [A.regType; A.regType];
       P.regSide = [A.regSide; A.regSide];
       P.seq = [ones(size(A.sn)); ones(size(A.sn))*2];  % 1 for trained, 2 for untrained
       P.sessN = [A.sessN; A.sessN];
       P.speed = [A.speed; A.speed];
       
       % % save
        save(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice)),'-struct','P');
    case 'PLOT_ratio_beta'
        betaChoice = 'multiPW';        
        roi=[1:8];
        parcelType = 'Brodmann';
        hemi=1;
        vararginoptions(varargin,{'betaChoice','roi','parcelType','hemi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice)));

       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.psc_beta,'split',S.seq,'subset',S.regType==f&S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
           ylim([-0.1 1.1]); hold on; drawline(1,'dir','horz'); drawline(0,'dir','horz');
           if f==1
               ylabel('Psc 2nd/1st betas')
           else
               ylabel('')
           end
           title(sprintf('%s',regname{f}))
       end
    case 'PLOT_ratio_psc'
        betaChoice = 'multiPW';
        roi=[1:8];
        parcelType = 'Brodmann';
        hemi=1;
        vararginoptions(varargin,{'betaChoice','roi','parcelType','hemi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice)));

       figure
       indx=1;
       for f = roi
           subplot(1,numel(roi),indx);
           plt.line([S.speed S.sessN],S.psc_psc,'split',S.seq,'subset',S.regType==f&S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
           ylim([-0.1 1.1]); hold on; drawline(1,'dir','horz'); drawline(0,'dir','horz');
           if indx==1
               ylabel('Psc 2nd/1st psc activation');
               xlabel('Session')
           else
               ylabel('');
           end
           title(sprintf('%s',regname{f}));
           indx=indx+1;
       end
    case 'PLOT_ratio_psc_sessions'
        betaChoice = 'multiPW';
        sessN=[1:3];
        seqType='trained';
        roi=[1:5,7,8];
        parcelType = 'Brodmann';
        hemi=1;
        vararginoptions(varargin,{'betaChoice','roi','parcelType','hemi','sessN'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice)));
        
        T=getrow(T,ismember(T.regType,roi)&ismember(T.sessN,sessN)&T.regSide==hemi); % exclude V12

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
        plt.line(T.reg,1-T.psc_psc,'split',T.sessN,'subset',T.seq==st,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        drawline(0,'dir','horz');
        title(sprintf('%s',seqType));
        ylabel('Percent repetition suppression of activation');
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI');        
    case 'PLOT_ratio_psc_seqType'
        sessN=[1:3];
        regExcl=6;
        roi=[1:5,7,8];
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        
        T=load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s.mat',betaChoice))); 
        
        T=getrow(T,ismember(T.roi,roi)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        for ss=sessN
            figure
            plt.bar(T.reg,1-T.psc_psc,'split',T.seq,'subset',T.sessN==ss,'leg',{'trained','untrained'},'leglocation','northeast','style',stySeq);
            drawline(0,'dir','horz');
            ylabel('Psc repetition suppression');
            set(gca,'XTick',[2 5 8 11 14 18 21],'XTickLabel',regname_cortex(reg));
            title(sprintf('sess-%d',ss));
            xlabel('ROI');
        end
    case 'STATS_ratio_psc'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s.mat',betaChoice)));
        
        % main ANOVA - session x sequence
        for r=roi
            fprintf('\n PSC in %s \n',regname{r});
            anovaMixed(T.psc_psc,T.sn,'within',[T.sessN T.seq],{'session','seqType'},'subset',T.roi==r);
        end
        
        % t-tests comparing sequence type per session
        for r=roi
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS ratio in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.psc_psc,[T.seq T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end
        
        % t-tests comparing across sessions - separately for T / UT seq
        for st=1:2
            fprintf('\n %s sequence\n',seqLabel{st});
            for sessTr=1:3
                for r=roi
                    fprintf('\n post-hoc t-test on comparison of sessions %d-%d in %s \n',sessTr,sessTr+1,regname{r});
                    ttestDirect(T.psc_psc,[T.sessN T.sn],2,'paired','subset',T.roi==r&ismember(T.sessN,[sessTr,sessTr+1])&T.seq==st);
                end               
            end
        end
    
    case 'tmp'
        betaChoice = 'multiPW';
        sessN=[1:4];
        seqType='trained';
        roi=[1:8];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        
        R = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s.mat',betaChoice)));  
        S = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s.mat',betaChoice)));
        TT=[];
        for ss=sessN
            for r=roi
                for st=1:2
                    r1 = getrow(R,R.sessN==ss&R.roi==r&R.seq==st);
                    s1 = getrow(S,S.sessN==ss&S.roi==r&S.seq==st);
                    T.cor = corr(r1.psc_psc,s1.subtr_dist);
                    T.roi = r;
                    T.sessN = ss;
                    T.seq = st;
                    TT=addstruct(TT,T);
                end
            end
        end
        keyboard;
    case 'CALC_ratio_distCross'
        % calculate percent repetition suppression across the two sequence
        % sets
    
       betaChoice = 'multiPW';
       roi=[1:8];
       parcelType='Brodmann';
       vararginoptions(varargin,{'betaChoice','roi','parcelType'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
 
       A = getrow(S,S.FoSEx==1);
       
       C.sn=A.sn;
       C.roi=A.roi;
       C.regType=A.regType;
       C.regSide=A.regSide;
       C.psc_distcross = [S.dist_cross(S.FoSEx==2)./S.dist_cross(S.FoSEx==1)];
       C.sessN=A.sessN;
       C.speed=A.speed;
       
       % % save
       save(fullfile(repSupDir,sprintf('RepSup_percent_crossseq_%s.mat',betaChoice)),'-struct','C');
        % calculate percent repetition suppression across the two sequence
        % sets
    case 'PLOT_ratio_dist'
        betaChoice = 'multiPW';
        roi=[1:8];
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice)));

       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.psc_dist,'split',S.seq,'subset',S.regType==f&S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
           hold on; drawline(1,'dir','horz'); drawline(0,'dir','horz');
           if f==1
               ylabel('Dist 2nd/1st betas')
           else
               ylabel('')
           end
           title(sprintf('%s',regname{f}))
       end
    case 'PLOT_ratio_dist_sessions'
        betaChoice = 'multiPW';
        sessN=[1:3];
        seqType='trained';
        roi=[1:5,7,8];
        vararginoptions(varargin,{'sessN','roi','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s.mat',betaChoice))); 
        
        T=getrow(T,ismember(T.roi,roi)&ismember(T.sessN,sessN)); % exclude V12

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
        plt.line(T.reg,T.psc_dist,'split',T.sessN,'subset',T.seq==st,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        drawline(0,'dir','horz');
        ylabel('Distance repetition suppression');
        set(gca,'XTickLabel',regname_cortex(reg));
        title(sprintf('%s',seqType));
        xlabel('ROI');          
    case 'PLOT_ratio_crossseq_sessions'
        % distance between sets of trained and untrained sequences
        betaChoice = 'multiPW';
        reg=[1:5,7,8];
        sessN=[1:3];
        sessLegend={'sess1','sess2','sess3','sess4'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_percent_crossseq_%s.mat',betaChoice)));
        
        T=getrow(T,ismember(T.roi,reg)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        figure
        plt.line(T.reg,T.psc_distcross,'split',T.sessN,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        drawline(0,'dir','horz');
        ylabel('Distance seqType repetition suppression');
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI');          
        
    case 'subject_RS_corr' % TO DO
        sessN=4;
        sn  = [1:9,11:16];
        roi = [1:8];
        betaChoice = 'multiPW'; % uni, multi or raw
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_subtract_%s.mat',betaChoice)));
        B = load(fullfile(repSupDir,'beh_scanner.mat'));
        
        RSS=[];
        
        for ss=sessN
            for s=sn
                for r=roi
               %TI = load(fullfile(glmFoSExDir{ss}, subj_name{s},'SPM_info.mat'));
                    RS.psc=S.subtr_psc(S.sn==s&S.sessN==ss&S.roi==r);
                    RS.MT(1,:)=mean(B.MT(B.sn==s&B.sessN==ss&B.FoSEx==1&B.seqType==1)-B.MT(B.sn==s&B.sessN==ss&B.FoSEx==2&B.seqType==1));
                    RS.MT(2,:)=mean(B.MT(B.sn==s&B.sessN==ss&B.FoSEx==1&B.seqType==2)-B.MT(B.sn==s&B.sessN==ss&B.FoSEx==2&B.seqType==2));
                    RS.seqType=[1 2]';
                    RS.roi=[r;r];
                    RS.sn=[s;s];
                    RS.sessN=[ss;ss];
                    RSS=addstruct(RSS,RS);
                end
            end
        end
        
        keyboard;
        
        for r=roi
            figure
            scatterplot(RSS.MT,RSS.psc,'split',RSS.seqType,'subset',RSS.roi==r,'leg',{'trained','untrained'});
            drawline(0,'dir','horz'); drawline(0,'dir','vert');
            xlabel('Movement time RS'); ylabel('Percent signal RS');
            title(sprintf('RS %s',regname{r}));
        end
        
    case 'trialwise_act_ROI' %------------------------------- TRIALWISE ANALYSIS (CLEAN-UP)--------------------
        % to check the model quality of the glm
        sn=[4:9,11:31];
        sessN=[1:4];
        vararginoptions(varargin,{'sn','sessN'});

        pre=0;         % How many TRs before the trial onset
        post=8;        % How many TRs after the trial onset
        T=[];
        for ss=sessN
            for s=sn
                fprintf('Extracting the onsets and events for subject %s, and session %d\n',subj_name{s},ss);
                load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM.mat'));
                SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s}));
                load(fullfile(regDir,[subj_name{s} '_regions.mat']));      % This is made in case 'ROI_define'
                [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
                
                % Create a structure with trial onset and trial type (event)
                D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
                
                % D.event - conditions (seqNumb: 1-12 - first exe; 13-24 - second exe)
                % D.block - run
                % add other indicators - FoSEx, seqNumb, seqType
                D.FoSEx=D.event;
                D.FoSEx(find(D.FoSEx<13))=1;
                D.FoSEx(find(D.FoSEx>12))=2;
                D.seqNumb=D.event;
                D.seqNumb(find(D.seqNumb>12))=D.seqNumb(find(D.seqNumb>12))-12;
                D.seqType=D.event;
                D.seqType(D.seqNumb<7)=1;
                D.seqType(D.seqNumb>6)=2;
                
                for r=1:size(y_raw,2)   % regions
                    S.block=D.block;
                    for i=1:(size(S.block,1));
                        S.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    end;
                    S.event=D.event;
                    S.FoSEx=D.FoSEx;
                    S.seqNumb=D.seqNumb;
                    S.seqType=D.seqType;
                    S.sn=ones(length(S.event),1)*s;
                    S.region=ones(length(S.event),1)*r; % region
                    S.regType=regType(S.region)';
                    S.regSide=regSide(S.region)';
                    
                    T=addstruct(T,S);
                end;
            end;
            cd(repSupDir);
            save(sprintf('trialwise_ACT_sess%d.mat',ss),'-struct','T');
        end
    case 'trialwise_BEH'
        sn=[1:9,11,12];
        vararginoptions(varargin,{'sn','sessN'});
        
        RS2=[];
        for s=sn
            % load behavioural data
            B=dload(fullfile(behavDir,sprintf('sml1_%s.dat',subj_name{s})));
            % only that scan = functional runs
            BB = getrow(B,B.ScanSess==sessN & (B.blockType==9 | B.blockType==3));
            
            for r = 1:numruns_task_sess
                uniqrun=unique(BB.BN);
                R = getrow(BB,BB.BN==uniqrun(r)); % 1-8 func runs of the session
                
                for c = 1:numel(num_seq) % each sequence
                    idx = find(R.seqNumb==c & R.FoSEx==1);   % find indx of all FIRST trials in run - 1:6 trained; 7-12 untrained
                    % difference in MT - 1st-2nd execution
                    RS.MT_dif=R.MT(idx)-R.MT(idx+1);
                    
                    for i=1:numel(idx)
                        % determine if no response in one of two trials
                        if (R.MT(idx(i))==0 | R.MT(idx(i)+1)==0)
                            RS.noResponse(i,:)=1;
                        else
                            RS.noResponse(i,:)=0;
                        end
                        
                        % determine error type
                        % 1 - no error
                        % 2 - 1st correct, 2nd error
                        % 3 - 1st error, 2nd correct
                        % 4 - both errors
                        if R.isError(idx(i))==0 & R.isError(idx(i)+1)==0
                            RS.errorType(i,:)=1;
                        elseif R.isError(idx(i))==0 & R.isError(idx(i)+1)==1
                            RS.errorType(i,:)=2;
                        elseif R.isError(idx(i))==1 & R.isError(idx(i)+1)==0
                            RS.errorType(i,:)=3;
                        else
                            RS.errorType(i,:)=4;
                        end
                        
                    end
                    
                    % some info about sequence
                    RS.seqType=R.seqType(idx);
                    RS.seqNumb=R.seqNumb(idx);
                    
                    % other general info - sn, sessN
                    RS.sn=ones(size(RS.seqType))*s;
                    RS.sessN=ones(size(RS.seqType))*sessN;
                    RS.run=ones(size(RS.seqType))*r;
                    
                    RS2=addstruct(RS2,RS);
                    
                end
            end
        end
        
        save(fullfile(repSupDir, sprintf('trialwise_BEH_sess%d.mat',sessN)),'-struct','RS2');
    case 'trialwise_calc_RS'
        sn=[1:9,11,12];
        sessN=1;
        vararginoptions(varargin,{'sn','sessN'});
        
        D=load(fullfile(repSupDir,sprintf('trialwise_ACT_sess%d.mat',sessN)));
        %T=load(fullfile(repSupDir,sprintf('trialwise_BEH_sess%d.mat',sessN)));
        exe1=getrow(D,D.FoSEx==1);
        exe2=getrow(D,D.FoSEx==2);
        
        % difference in maximum activation in each trial (1st - 2nd)
        B.RS_maxact=max(exe1.y_adj(:,[3:8]),[],2)-max(exe2.y_adj(:,[3:8]),[],2); 
        % mean activation
        B.RS_meanact=mean(exe1.y_adj(:,[3:8]),2)-mean(exe2.y_adj(:,[3:8]),2);
        % overall activation difference (cumulative) within the trial
        B.RS_cumsum=sum(exe1.y_adj(:,[3:8])-exe2.y_adj(:,[3:8]),2);
        B.region=exe1.region;
        B.regType=exe1.regType;
        B.regSide=exe1.regSide;
        B.seqNumb=exe1.seqNumb;
        B.seqType=exe1.seqType;
        B.sn=exe1.sn;
        
       save(fullfile(repSupDir, sprintf('trialwise_RS_sess%d.mat',sessN)),'-struct','B');
    case 'trialwise_error_ACT_corr'
        sn=[1:9,11];
        sessN=1;
        roi=[1:8];
        vararginoptions(varargin,{'sn','sessN','roi'});
        
        RS=load(fullfile(repSupDir,sprintf('trialwise_RS_sess%d.mat',sessN)));
        D=load(fullfile(repSupDir,sprintf('trialwise_BEH_sess%d.mat',sessN)));
        

        for r=roi
            R=getrow(RS,RS.region==r);

            figure
            subplot(1,3,1)
            plt.line(D.errorType,R.RS_maxact,'subset',D.noResponse==0,'style',stySeq,'split',D.seqType,'leg',{'trained','untrained'});
            drawline(0,'dir','horz'); 
            xlabel('Error type'); ylabel('Difference in max act 1st-2nd');
            title(sprintf('RS for different error types session %d - %s',sessN,regname{r}));


            subplot(1,3,2)
            plt.line(D.errorType,R.RS_meanact,'subset',D.noResponse==0,'style',stySeq,'split',D.seqType,'leg',{'trained','untrained'});
            drawline(0,'dir','horz'); 

            
            subplot(1,3,3)
            plt.line(D.errorType,R.RS_cumsum,'subset',D.noResponse==0,'style',stySeq,'split',D.seqType,'leg',{'trained','untrained'});
            drawline(0,'dir','horz'); 
           
            
        end
    case 'trialwise_error_ACT_sessions'
        sn=[1:9,11];
        sessN=[1:4];
        roi=[1:8];
        vararginoptions(varargin,{'sn','sessN','roi'});
        
        for ss=sessN
            RS=load(fullfile(repSupDir,sprintf('trialwise_RS_sess%d.mat',ss)));
            D=load(fullfile(repSupDir,sprintf('trialwise_BEH_sess%d.mat',ss)));
            for r=roi
                R=getrow(RS,RS.region==r);
                
                figure(r)
                subplot(1,max(sessN),ss)
                lineplot(D.errorType,R.RS_meanact,'subset',D.noResponse==0,'style_thickline','split',D.seqType,'leg',{'trained','untrained'});
                drawline(0,'dir','horz');
                set(gca,'XTickLabel',{'1-1','1-0','0-1','0-0'});xlabel('Error type'); 
                if ss==1
                    ylabel('Difference in mean act 1st-2nd');
                    title(sprintf('RS %s - session %d',regname{r},ss));
                else
                    ylabel('');
                    title(sprintf('session %d',ss));
                end
                
            end
        end
    case 'trialwise_MT_ACT_corr'
        sn=[1:9,11,12];
        sessN=[1:4];
        roi=[1:8];
        vararginoptions(varargin,{'sn','sessN','roi'});
        
        for ss=sessN
            RS=load(fullfile(repSupDir,sprintf('trialwise_RS_sess%d.mat',ss)));
            D=load(fullfile(repSupDir,sprintf('trialwise_BEH_sess%d.mat',ss)));
            
            for r=roi
                R=getrow(RS,RS.region==r);
                
                figure(r)
                subplot(1,max(sessN),ss)
                scatterplot(D.MT_dif,R.RS_maxact,'subset',D.noResponse==0&D.errorType,'split',D.seqType,'leg',{'trained','untrained'});
                drawline(0,'dir','horz'); drawline(0,'dir','vert');
                xlabel('Movement time diff 1st-2nd'); ylabel('Act difference 1st-2nd');
                title(sprintf('Behaviour activation RS session %d - %s',ss,regname{r}));
            end
        end
     
    case 'ROI_getBetas_trialwise'
        sn=[4:9,11:31];
        sessN=[1:4];
        roi=[1:8];
        vararginoptions(varargin,{'sn','sessN','roi'});
        
        for ss=sessN
            SS=[];
            for s=sn
                % load files
                load(fullfile(glmTrialDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                T=load(fullfile(glmTrialDir{ss},subj_name{s},'SPM_info.mat'));
                load(fullfile(regDir,[subj_name{s} '_regions.mat']));        % load subject's region parcellation (R)
                V = SPM.xY.VY;
                
                for r=roi
                    % get raw data for voxels in region
                    Y = region_getdata(V,R{r});  % Data Y is N x P
                    
                    % estimate region betas
                    [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                    indTrial=repmat([1:72],1,8)'; % without the intercept
                    %indTrial = repmat([1:73],1,8)'; % indicator for trials - 73 is intercept - exclude later!
                    
                    S.betaW                   = {betaW(indTrial,:)};                             % multivariate pw
                    S.betaUW                  = {bsxfun(@rdivide,beta(indTrial),sqrt(resMS))}; % univariate pw
                    S.betaRAW                 = {beta(indTrial,:)};
                    
                    %S.betaW                   = {betaW(indTrial~=73,:)};                             % multivariate pw
                    %S.betaUW                  = {bsxfun(@rdivide,beta(indTrial~=73),sqrt(resMS))}; % univariate pw
                    %S.betaRAW                 = {beta(indTrial~=73,:)};
                    
                    %S.trialInd                = {indTrial(indTrial~=73)};
                    S.trialInd                = {indTrial(indTrial)};
                    S.run                     = {T.run};
                    S.seqNumb                 = {T.seqNumb};
                    S.seqType                 = {T.seqType};
                    S.MT                      = {T.MT};
                    S.isError                 = {T.isError};
                    S.region                  = r;
                    S.sn                      = s;
                    S.sessN                   = ss;
                    SS = addstruct(SS,S);
       
                end
            end
            save(fullfile(repSupDir,sprintf('betas_trialwiseGLM_sess%d.mat',ss)),'-struct','SS');
        end
    case 'Betas_MT_trialwise'
        sn=[1:9,11,12];
        sessN=[1:4];
        roi=[1:8];
        betaChoice='multi';
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice'});
        
        D=[];
        
        for ss=sessN
            T=load(fullfile(repSupDir,sprintf('betas_trialwiseGLM_sess%d.mat',ss)));
            B=load(fullfile(repSupDir,sprintf('trialwise_BEH_sess%d.mat',ss)));
            
            for r=roi
                for s=sn
                    t=getrow(T,T.region==r&T.sn==s);
                    b=getrow(B,B.sn==s);
                    switch(betaChoice)
                        case 'multi'
                            beta=t.betaW{1};
                        case 'uni'
                            beta=t.betaUW{1};
                        case 'raw'
                            beta=t.betaRAW{1};
                    end
                    beta1=beta(1:2:end,:);  % first execution
                    beta2=beta(2:2:end,:);  % second execution
                   
                    d.betaDiff=nanmean(beta1-beta2,2);
                    d.MT_dif=b.MT_dif;
                    d.errorType=b.errorType;
                    d.noResponse=b.noResponse;
                    d.seqNumb=b.seqNumb;
                    d.seqType=b.seqType;
                    d.reg=ones(size(d.seqType))*r;
                    d.sn=ones(size(d.seqType))*s;
                    d.sessN=ones(size(d.seqType))*ss;

                    D=addstruct(D,d);
                end
            end
        end
        % save
        save(fullfile(repSupDir,'ACT_BEH_trialwise.mat'),'-struct','D');
    case 'PLOT_trialwise_Beta_MT'
        roi=[1:8];
        sessN=[1:4];
        sn=[1:9,11,12];
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(repSupDir,'ACT_BEH_trialwise.mat'));

        for ss=sessN
            for r=roi
                figure(r)
                subplot(1,numel(sessN),ss)
                scatterplot(T.MT_dif,T.betaDiff,'subset',(T.errorType==1&T.reg==roi(r)&T.sessN==ss),'split',T.seqType);
            end
        end
        keyboard;
    case 'PLOT_trialwise_Beta_errorType'
        roi=[1:8];
        sessN=[1:4];
        sn=[1:9,11,12];
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(repSupDir,'ACT_BEH_trialwise.mat'));

        for ss=sessN
            for r=roi
                figure(r)
                subplot(1,numel(sessN),ss)
                plt.bar(T.errorType,T.betaDiff,'subset',T.sessN==ss&T.reg==roi(r),'split',T.seqType,'leg',{'Trained','Untrained'},'leglocation','north');
                plt.match('y');
                xlabel('Error type');
                if ss==1
                    ylabel('Beta difference 1st - 2nd');
                else
                    ylabel('');
                end
                title(sprintf('%s sess%d',regname{r},ss));
            end
        end
        keyboard;
    case 'STATS_trialwise_Beta_errorType'
        roi=[1:8];
        sessN=[1:4];
        sn=[1:9,11,12];
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(repSupDir,'ACT_BEH_trialwise.mat'));
        T.error2=T.errorType;
        T.error2(T.errorType>1)=2;
        for r=roi
            for ss=sessN
                fprintf('\n RS with respect to errors in %s session %d \n',regname{r},ss);
                anovaMixed(T.betaDiff,T.sn,'within',[T.error2 T.seqType],{'errorType','seqType'},'subset',T.reg==r&T.sessN==ss);
            end
        end
        
        keyboard;
        
    case 'MDS_seqType_fromRest' %------------------------------------ MDS -----------------------------------
        % MDS from activation for each sequence type from rest 
         roi=[1:8]; 
         sessN=[1:4];
         betaChoice='multi';
         sn=[4:9,11:31];
         regType = 'cortex';
         vararginoptions(varargin,{'roi','sessN','betaChoice','sn','regType'});
         
         seqNumb = kron(ones(2,1),[1:6]'); % 1-6
         seqType  = kron([1:2]',ones(6,1)); % 1: trained; 2: untrained
         color={'r','b'};
         
         for ss=sessN
             switch regType
                 case 'cortex'
                     T=load(fullfile(betaDir,'group',sprintf('stats_FoSEx_Brodmann_%sPW_sess%d.mat',betaChoice,ss)));
                     regLabel = regname_cortex;
                 case 'striatum'
                     T=load(fullfile(betaDir,'group',sprintf('stats_FoSEx_BG-striatum_%sPW_sess%d.mat',betaChoice,ss)));
                     regLabel = {'Caudate','Putamen'};
             end
             for r=roi
                 D = getrow(T,T.region==r & ismember(T.SN,sn));
                 for exe=1:2
                     IPM = mean(D.IPM(D.FoSEx==exe,:));

                     % Overall task activity with rest included
                     figure(r)
                     ax(exe)=subplot(2,4,(exe-1)*4+ss);
                     Z = indicatorMatrix('identity',seqType);
                     [Y1,l1,V1]=rsa_classicalMDS(IPM,'mode','IPM','contrast',Z);
                     scatterplot(Y1(:,1),Y1(:,2),'split',seqType,...
                         'markercolor',color,'xaxisIncl',0,'yaxisIncl',0);
                     if ss==1&exe==1
                         legend('trained','untrained','Location','SouthEast');
                         legend(ax(exe),'boxoff'); 
                     end
                     axis equal;
                     hold on;
                     mVector=pinv(Z)*Y1;   % Category vectors
                     line([0;mVector(1,1)],[0;mVector(1,2)],'color','r');
                     line([0;mVector(2,1)],[0;mVector(2,2)],'color','b');
                     plot(0,0,'k+');
                     hold off;
                     if exe==2
                         linkaxes(ax);
                         axis equal;
                     end
                     title(sprintf('Session %d %s execution %d',ss,regLabel{r},exe),'FontSize',8);
                     
                     set(gcf,'PaperPosition',[2 2 4 8],'Color',[1 1 1]);
                 end
             end
         end
    case 'MDS_seqType_exe'
        roi=[1:8]; 
         sessN=[1:4];
         betaChoice='multi';
         parcelType='Brodmann';
         sn=[4:9,11:25];
         vararginoptions(varargin,{'roi','sessN','betaChoice','sn','regType'});
         
         seqNumb = kron(ones(2,1),[1:6]'); % 1-6
         seqType  = kron([1:2]',ones(6,1)); % 1: trained; 2: untrained
         color{1}={'r','b'}; % first execution
         color{2}={'m','c'}; % second execution
         
         for ss=sessN
            T=load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
             switch parcelType
                 case 'Brodmann'
                     regLabel = regname_cortex;
                 case 'BG-striatum'
                     regLabel = {'Caudate','Putamen'};
             end
             for r=roi
                 D = getrow(T,T.region==r & ismember(T.SN,sn));
                 
                 for exe=1:2
                     IPM = mean(D.IPM(D.FoSEx==exe,:));

                     % Overall task activity with rest included
                     figure(r)
                     ax(exe)=subplot(1,4,ss);
                     Z = indicatorMatrix('identity',seqType);
                     [Y1,l1,V1]=rsa_classicalMDS(IPM,'mode','IPM','contrast',Z);
                     scatterplot(Y1(:,1),Y1(:,2),'split',seqType,...
                         'markercolor',color{exe},'xaxisIncl',0,'yaxisIncl',0);
                     if ss==1&exe==1
                         legend('trained','untrained','Location','SouthEast');
                         legend(ax(exe),'boxoff'); 
                     end
                     axis equal;
                     
                     hold on;
                     mVector=pinv(Z)*Y1;   % Category vectors
                     line([0;mVector(1,1)],[0;mVector(1,2)],'color',color{exe}{1});
                     line([0;mVector(2,1)],[0;mVector(2,2)],'color',color{exe}{2});
                     plot(0,0,'k+');
                     
                     axis equal;
                     title(sprintf('Session %d %s',ss,regLabel{r}),'FontSize',8);
                     
                 end
             end
         end
    case 'MDS_seqType_relative'
        % Relative MDS when average activation 0 - optimised for seqType
        % difference
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
            T=load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
             switch parcelType
                 case 'Brodmann'              
                     regLabel = regname_cortex;
                 case 'BG-striatum'
                     regLabel = {'Caudate','Putamen'};
             end
            for r=roi
                D = getrow(T,T.region==r & ismember(T.SN,sn));
                
                for exe=1:2
                    IPM = mean(D.IPM(D.FoSEx==exe,:));
                    
                    Z = indicatorMatrix('identity_p',[1:12]);
                    Z=bsxfun(@minus,Z,mean(Z,2));
                    [Y3,l3,V3]=rsa_classicalMDS(IPM,'mode','IPM','contrast',Z);
                    figure(r)
                    ax(exe)=subplot(2,4,(exe-1)*4+ss);
                    scatterplot3(Y3(:,1),Y3(:,2),Y3(:,3),'split',seqType,'label',seqNumb,'markercolor',color,'markerfill',color);
       
                    for i=1:2
                        hold on;
                        indx=find(seqType==i);
                        indx(end+1)=indx(1);
                        line(Y3(indx,1),Y3(indx,2),Y3(indx,3),'color',color{i});
                    end;
                    plot3(0,0,0,'k+');
                    title(sprintf('Session %d %s execution %d',ss,regLabel{r},exe),'FontSize',8);
                    axis equal;
                    
                    hold off;
                    set(gcf,'PaperPosition',[2 2 4 8],'Color',[1 1 1]);    
                end
            end
        end
    case 'MDS_reducedSeq'
        % 'schematic' MDS for two random sequences per trained / untrained
        % 4 points - 1st / 2nd trained / untrained
        roi=[1:8];
        sessN=[1:4];
        sn=[4:9,11:31];
        parcelType='Brodmann';
        
        color={'r','b'}; 
        seq={'trained','untrained'};
        
        %seqExe=kron([1:2]',ones(2,1));
        seqID=[1;2;1;2];
        vararginoptions(varargin,{'roi','sessN','sn','regType'});
        
        for ss=sessN
            T=load(fullfile(regDir,sprintf('reduced_IPM_FoSEx_%s_sess%d.mat',parcelType,ss)));
            
            switch parcelType
                case 'Brodmann'
                     regLabel = regname_cortex;
                 case 'BG-striatum'
                     regLabel = {'Caudate','Putamen'};
             end
            
            for r=roi
                D = getrow(T,T.region==r);
                
                for seqT=1:2

                    IPM = mean(D.IPM_red(D.seqType==seqT,:));
                    Z = indicatorMatrix('identity',seqID);
                    
                    figure(r)
                    subplot(2,4,(seqT-1)*4+ss)
                    [Y1,l1,V1]=rsa_classicalMDS(IPM,'mode','IPM','contrast',Z);
                    scatterplot(Y1(:,1),Y1(:,2),'split',seqID,...
                        'markercolor',color,'xaxisIncl',0,'yaxisIncl',0);
                    
                    hold on;
                    %mVector=pinv(Z)*Y1;   % Category vectors
                    % Seq A first
                    line([0,Y1(1,1)],[0,Y1(1,2)],'color','r');
                    % Seq B first
                    line([0,Y1(2,1)],[0,Y1(2,2)],'color','b');
                    
                    % Seq A second
                    line([0,Y1(3,1)],[0,Y1(3,2)],'color','r','linestyle','--');
                    % Seq B second
                    line([0,Y1(4,1)],[0,Y1(4,2)],'color','b','linestyle','--');

                    plot(0,0,'k+');
                    
                    title(sprintf('Session %d %s %s',ss,seq{seqT},regLabel{r}),'FontSize',8);
                    
                end
            end
            
        end
     
    case 'dist_corr_Eucl' % ------------------------------- VECTOR LENGTH, ANGLE --------------------------------
        % predict the distance for second exe if correlation / Euclidean
        % distance same across repetitions
        roi=[1:8];
        sessN=[1:4];
        sn=[1:9,11:25];
        
        seqID=[1;2;1;2];
        seqExe=[1;1;2;2];
        vararginoptions(varargin,{'roi','sessN','sn'});
        CC=[];
        
        for ss=sessN
            T=load(fullfile(regDir,sprintf('reduced_IPM_FoSEx_sess%d.mat',ss)));
            for s=sn
                for r=roi
                    for st=1:2  % seqType - trained / untrained
                        D = getrow(T,T.region==r&T.SN==s&T.seqType==st);
                        
                        G=rsa_squareIPMfull(D.IPM_red_full);
                        % calculate distances for first execution -
                        % Euclidean and correlation
                        distEuc = G(1,1)+G(2,2)-2*G(1,2);
                        distCorr = G(1,2)/sqrt(G(1,1)*G(2,2));

                        % use distances to calculate expected cov for 2nd exe
                        C.covEuc = (distEuc-G(3,3)-G(4,4))/2;
                        C.covCorr = distCorr*sqrt(G(3,3)*G(4,4));
                        C.covReal = G(3,4);
                        C.seqType=st;
                        C.sn=s;
                        C.roi=r;
                        C.sessN=ss;
                        CC=addstruct(CC,C);
                    end
                    
                end
            end
        end

        CC.covEuc=abs(CC.covEuc);

        % save
        save(fullfile(repSupDir,'Distance_stats.mat'),'-struct','CC');   
    case 'STATS_dist_corrEucl'
        sessN=[1:4];
        vararginoptions(varargin,{'sessN'});
        
        D=load(fullfile(repSupDir,'Distance_stats.mat'));
        
        covEuc_diff=D.covReal-D.covEuc;
        covCorr_diff=D.covReal-D.covCorr;
        C.cov_diff=[covEuc_diff;covCorr_diff];
        C.cov_type=[ones(size(covEuc_diff));ones(size(covEuc_diff))*2];
        C.sn=[D.sn;D.sn];
        C.roi=[D.roi;D.roi];
        C.sessN=[D.sessN;D.sessN];
        C.seqType=[D.seqType;D.seqType];
        
        seqType={'trained','untrained'};
        % do t-tests per region, per session, per seqType
        for ss=sessN
            for st=1:2
                fprintf('\n Session %d - %s sequence onesample t-tests of cov values calculated with Euclidian distance against real covariance values \n',ss,seqType{st});
                [t,p]=ttestDirect(C.cov_diff,C.sn,2,'onesample','subset',C.sessN==ss&C.seqType==st&C.cov_type==1,'split',C.roi);
                fprintf('\n Session %d - %s sequence onesample t-tests of cov values calculated with corr distance against real covariance values \n',ss,seqType{st});
                [t,p]=ttestDirect(C.cov_diff,C.sn,2,'onesample','subset',C.sessN==ss&C.seqType==st&C.cov_type==2,'split',C.roi);
            end
        end
        
        keyboard;
    case 'Vector_AngleLength_repsup'
        % calculated from IPM
        % representation of patterns in multivariate space 1st->2nd exe
        % vector length as variance (diagonal)
        % vector angle as covariance (later divided by variance to normalise)
        roi=[1:8];
        sessN=[1:4];
        sn=[1:9,11:25];
        betaChoice='multi';

        vararginoptions(varargin,{'roi','sessN','sn','betaChoice'});
        VV=[];
        for ss=sessN
            T=load(fullfile(regDir,sprintf('stats_FoSEx_%sPW_sess%d.mat',betaChoice,ss)));
            for s=sn
                for r=roi
                    for exe=1:2
                        ipm=T.IPMfull(T.region==r&T.SN==s&T.FoSEx==exe,:);
                        G=rsa_squareIPMfull(ipm);
                        G_train=G(1:6,1:6);
                        G_untrain=G(7:12,7:12);
                        V.var(1,:)=mean(diag(G_train)); % mean of diagonal
                        V.var(2,:)=mean(diag(G_untrain));
                        V.cov(1,:)=sum(sum(triu(G_train,1)))/(6*5/2); % mean of off-diagonal elements - variances between sequences
                        V.cov(2,:)=sum(sum(triu(G_untrain,1)))/(6*5/2);
                        V.dist(1,:)=2*V.var(1,:)-2*V.cov(1,:);
                        V.dist(2,:)=2*V.var(2,:)-2*V.cov(2,:);
                        V.seqType=[1;2];
                        V.exe=ones(size(V.seqType))*exe;
                        V.sn=ones(size(V.seqType))*s;
                        V.roi=ones(size(V.seqType))*r;
                        V.sessN=ones(size(V.seqType))*ss;
                        VV=addstruct(VV,V);
                    end  
                end
            end
        end
        % save
        save(fullfile(repSupDir,'RepSup_Var_Cov.mat'),'-struct','VV');
    case 'PLOT_multivar_vectorLength'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        seqType={'trained','untrained'};
        
        T=load(fullfile(repSupDir,'RepSup_Var_Cov.mat'));
        
        for r=roi
            figure(r)
            for st=1:2
                subplot(1,2,st)
                plt.line([T.sessN>3 T.sessN],T.var,'split',T.exe,'subset',T.seqType==st&T.roi==r,'leg',{'1st','2nd'},'style',stySeqType_exe);
                plt.match('y');
                xlabel('Session number');
                title(sprintf('%s %s',regname{r},seqType{st}));
                if st==1
                    ylabel('Vector length in multivariate patterns');
                else
                    ylabel('');
                end
            end
        end
    case 'PLOT_multiCov_vectorAngle'
        roi=[1:8];
        vararginoptions(varargin,{'roi','sessN'});
        seqType={'trained','untrained'};
        
        T=load(fullfile(repSupDir,'RepSup_Var_Cov.mat'));
        
        for r=roi
            figure(r)
            for st=1:2
                subplot(1,2,st)
                plt.line([T.sessN>3 T.sessN],T.dist./T.var,'split',T.exe,'subset',T.seqType==st&T.roi==r,'leg',{'1st','2nd'},'style',stySeqType_exe);
                plt.match('y');
                xlabel('Session number');
                title(sprintf('%s %s',regname{r},seqType{st}));
                if st==1
                    ylabel('Vector angle / cos dist of patterns');
                else
                    ylabel('');
                end
            end
        end
    case 'PLOT_multiCov_offDiag'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        seqType={'trained','untrained'};
        
        T=load(fullfile(repSupDir,'RepSup_Var_Cov.mat'));
        
        for r=roi
            figure(r)
            for st=1:2
                subplot(1,2,st)
                plt.line([T.sessN>3 T.sessN],T.cov,'split',T.exe,'subset',T.seqType==st&T.roi==r,'leg',{'1st','2nd'},'style',stySeqType_exe);
                plt.match('y');
                xlabel('Session number');
                title(sprintf('%s %s',regname{r},seqType{st}));
                if st==1
                    ylabel('Cov of patterns (off-diagonal)');
                else
                    ylabel('');
                end
            end
        end
    
    case 'STATS_vectorLength'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        
        seqType={'trained','untrained'};
        
        T=load(fullfile(repSupDir,'RepSup_Var_Cov.mat'));
          
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.var,T.sn,'within',[T.sessN T.seqType T.exe],{'session','seqType','repsup'},'subset',T.roi==r);
        end
        
        keyboard;
    case 'STATS_vectorAngle'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        
        seqType={'trained','untrained'};
        
        T=load(fullfile(repSupDir,'RepSup_Var_Cov.mat'));
          
        for r=roi
            fprintf('\n Angular dist in %s \n',regname{r});
            anovaMixed(T.dist./T.var,T.sn,'within',[T.sessN T.seqType T.exe],{'session','seqType','repsup'},'subset',T.roi==r);
        end
        
        keyboard;
    
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
    case 'PCM_dataCorr'         % USE THIS FUNCTION
        % calculate correlation between 1st and 2nd execution
        % naive, crossval, pcm correlation
        % generic or specific model - parameter for all 6 seq / per seq
        runEffect  = 'fixed';
        beta_choice = 'mw';
        reg = [1:8];
        parcelType='cortex'; % cortex or striatum
        sn=[4:9,11:25];
        sessN=[1:4]; % need to be two sessions at the time
        seqType={'trained','untrained'};
        modelType='generic'; %generic or specific
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','seqType','modelType','parcelType'})
        
        for stIndx = 1:2    % seqType
            for ss=sessN
                CC=[];
                for r = reg
                    for p=1:length(sn)
                        switch parcelType
                            case 'cortex'
                                B=load(fullfile(regDir,sprintf('betas_FoSEx_sess%d.mat',ss)));
                            case 'cortex_networks'
                                B=load(fullfile(regDir,sprintf('betas_FoSEx_cortex_buckner_sess%d.mat',ss)));
                            case 'striatum'
                                B=load(fullfile(regDir,sprintf('betas_BG-striatum_FoSEx_sess%d.mat',ss)));
                            case 'striatum_networks'
                                B=load(fullfile(regDir,sprintf('betas_striatum_buckner_hemi_FoSEx_sess%d.mat',ss)));
                        end
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.SN==sn(p)&B.region==r)};
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p)&B.region==r)}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                        end
                        %squeeze dimension of betaW if voxels with no data
                        if strcmp(parcelType,'striatum_networks')
                            indx=1;
                            tmp=[];
                            for v=1:size(beta,2)
                                if sum(isnan(beta(:,v)))==size(beta,1)% no data
                                else
                                    tmp(:,indx)=beta(:,v);
                                    indx=indx+1;
                                end
                            end
                            beta=[];
                            beta=tmp';
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
                        %Data{p} = beta(:,D.seqType==stIndx)';  % Data is N x P
                        %(cond x voxels) - no intercept used with mw
                        Data{p} = beta(D.seqType==stIndx,:);
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
                    C     = pcm_correlation(Data,partVec,condVec,M{2},runEffect,mf);
                    C.SN  = [1:size(C.r_naive,1)]';
                    C.roi = ones(size(C.r_naive))*r;
                    CC    = addstruct(CC,C);
                    fprintf('Done calculating corr for region-%d\n\n\n',r);
                end
                CC=rmfield(CC,'theta');
                % save output
                save(fullfile(pcmDir,sprintf('PCM_repsup_corr_%s_sess%d_%s_%s.mat',parcelType,ss,seqType{stIndx},modelType)),'-struct','CC');
               % save(fullfile(pcmDir,sprintf('PCM_cov_%s_sess%d_%s_%s.mat',parcelType,ss,seqType{stIndx},modelType)),'-struct','CC');
                fprintf('Done all regions sess-%d\n\n\n\n\n\n',ss);
            end
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
    
    case 'PCM_constructModelFamily'  
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        parcelType='Brodmann'; % Brodmann or 162tessels
        fingType='Count';
        regSelect='all'; % all or subset
        sn=[4:9,11:31]; % change
        sessN=[1:4];
        AllReg=[];
        KK=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType'})
        if strcmp(parcelType,'162tessels')
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            if ~strcmp(parcelType,'162tessels')
                reg=unique(B.region)';
                regSide=(reg>size(reg,2)/2)+1;
                regType=reg-((regSide-1)*size(reg,2)/2);
            else % only for 162tessels
                regType=reg;
                regType(regType>158)=regType(regType>158)-158;
                regSide=ones(size(regType));
                regSide(reg>158)=2;
            end
            for r = 1:length(reg)
                for rep=1:2
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p)&B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                        end
                        
                        idx=D.FoSEx==rep;
                        m = pcm_defineSequenceModels_fixed_new(Seq,SeqChunks,p);
                        [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        condVec{p} = repmat(Z,max(D.run),1);
                        partVec{p} = D.run(idx);  % runs/partitions
                        Data{p} = beta(:,(idx==1)')';  % Data is N x P (cond x voxels) - no intercept
                    end;         
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T.exe   = ones(size(T.SN))*rep;
                    T=rmfield(T,{'reg','thetaCr','theta_hat'});
                    AllReg=addstruct(AllReg,T);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.sessN = ones(length(K.indx),1)*ss;
                    K.exe   = ones(length(K.indx),1)*rep;   
                    KK=addstruct(KK,K);
                end
            end
        end
        
        % save variables;
        dircheck(fullfile(pcmDir));
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_FoSEx_%s_new.mat',parcelType)),'Comb');
        save(fullfile(pcmDir,sprintf('ModelFamily_FoSEx_Fit_%s_fing%s_new.mat',parcelType,fingType)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_FoSEx_Stats_%s_fing%s_new.mat',parcelType,fingType)),'-struct','KK');
        
        
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
        parcelType = 'cortex';
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','parcelType'});
        
        switch parcelType
            case 'cortex'
                lab = regname_cortex;
            case 'striatum'
                lab = {'Caudate','Putamen'};
            case 'striatum_networks'
                lab = {'Network-1','Network-2','Network-3','Network-4',...
                    'Network-5','Network-6','Network-7'};
            case 'cortex_networks'
                lab = {'Network-1','Network-2','Network-3','Network-4',...
                    'Network-5','Network-6','Network-7'};
        end
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_%s_sess%d_%s_%s.mat',parcelType,t,seqType{st},modelType)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        corrType={'naive','crossval','pcm','mean_naive','mean_crossval'};
        corrVar=[T.r_naive, T.r_crossval3, T.r_model2, T.r_naive_meanOnly, T.r_crossval_meanOnly_evenOdd];
        for ct=1:length(corrVar) % three types of correlation 
            figure(ct)
            sub=1;
            for r=reg
                subplot(1,numel(reg),sub)
                plt.line([T.sessN>3 T.sessN],corrVar(:,ct),'subset',T.roi==r,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','north');
                plt.match('y');
                drawline(0,'dir','horz');
                title(sprintf('%s',lab{r}));
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
    case 'PLOT_pcm_corr_allSess_NEW' %delete?
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
        
    case 'PLOT_rdmCorr'
        reg = [1:8];
        sessN=[1:4];
        seqType={'trained','untrained'};
        modelType='generic';
        parcelType = 'cortex';
        metric = 'RDM_corrRSA'; % RDM_corrRSA or RDM_corrDist
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','parcelType','metric'});
        
        switch parcelType
            case 'cortex'
                lab = regname_cortex;
            case 'striatum'
                lab = {'Caudate','Putamen'};
        end
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_%s_sess%d_%s_%s.mat',parcelType,t,seqType{st},modelType)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        figure
        sub=1;
        for r=reg
            subplot(1,numel(reg),sub)
            plt.line([T.sessN>3 T.sessN],T.(metric),'subset',T.roi==r,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','north');
            plt.match('y');
            drawline(0,'dir','horz');
            %title(sprintf('%s',lab{r}));
            if sub==1
                ylabel(metric);
                xlabel('Session');
            else
                ylabel('');
            end
            sub=sub+1;
        end
    case 'STATS_modelPCM'
        reg = [1:8];
        sessN=[1:4];
        runEffect='fixed';
        vararginoptions(varargin,{'reg','sessN','runEffect'});
        
        seqType={'trained','untrained'};
        model_label={'flexible','perfect'};
        
        T=[];
        for t=sessN
            for st=1:2
                R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_%s_sess%d_%s.mat',seqType{st},t,runEffect)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        T2.modelInd=[ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3];
        T2.bayesEst=T.bayesEst(:);
        T2.roi=[T.roi;T.roi;T.roi];
        T2.sessN=[T.sessN;T.sessN;T.sessN];
        T2.sn=[T.SN;T.SN;T.SN];
        T2.seqType=[T.seqType;T.seqType;T.seqType];
        
        for r=reg
            for mIndx=1:2 % models 2 (flex) and 3 (perfect)
                fprintf('\n %s Effect seqType and session on %s model\n',regname{r},model_label{mIndx});
                anovaMixed(T2.bayesEst,T2.sn,'within',[T2.sessN T2.seqType],{'session','seqType'},'subset',T2.roi==r&T2.modelInd==mIndx+1);
            end
        end
        
        for ss=sessN
            for r=reg
                for mIndx=1:2 % models 2 (flex) and 3 (perfect)
                    fprintf('\n %s T-test on seqType in sess %d - %s model\n',regname{r},ss,model_label{mIndx});
                    ttestDirect(T2.bayesEst,[T2.seqType T2.sn],2,'paired','subset',T2.roi==r&T2.modelInd==mIndx+1&T2.sessN==ss);
                end
            end
        end
        
        
        keyboard;
    case 'STATS_corrPCM'
        reg = [1:8];
        sessN=[1:4];
        seqType={'trained','untrained'};
        runEffect='fixed';
        modelType='generic';
        vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType'});
        
        T=[];
        for st=1:size(seqType,2)
            for t=sessN
                 R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_NEW_%s_%s_sess%d_%s.mat',modelType,seqType{st},t,runEffect)));
              %  R=load(fullfile(pcmDir,sprintf('PCM_repsup_reliability_%s_sess%d_%s.mat',seqType{st},t,runEffect)));
                R.sessN=ones(size(R.SN))*t;
                R.seqType=ones(size(R.SN))*st;
                T=addstruct(T,R);
            end
        end
        
        for r=reg
            fprintf('\n %s Effect of seqType and session on crossval correlation\n',regname{r});
            anovaMixed(T.r_crossval,T.SN,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r);
            fprintf('\n %s Effect of seqType and session on pcm correlation\n',regname{r});
            anovaMixed(T.r_model2,T.SN,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r);  
         end
        
         for ss=sessN
             for r=reg
                 fprintf('\n %s T-test on crossval corr in sess %d\n',regname{r},ss);
                 ttestDirect(T.r_crossval,[T.seqType T.SN],2,'paired','subset',T.roi==r&T.sessN==ss);
                 fprintf('\n %s T-test on model corr in sess %d\n',regname{r},ss);
                 ttestDirect(T.r_model2,[T.seqType T.SN],2,'paired','subset',T.roi==r&T.sessN==ss);
             end
         end
         
         keyboard;
    
    case 'CORR_splithalf' % ------------------------------- NEW : SPLIT-HALF Correlation
        sn=[1:9,11:25];
        sessN=[1:4];
        betaChoice='multiPW';
        parcelType = 'striatum_buckner_hemi';
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
            %T = load(fullfile(regDir,sprintf('betas_FoSEx_sess%d',ss)));
            T = load(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d',parcelType,ss)));
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
        %save(fullfile(repSupDir,sprintf('corr_splitHalf_NEW_exe_%s',betaChoice)),'-struct','AllCorr');
        save(fullfile(repSupDir,sprintf('corr_splitHalf_%s_exe_%s',parcelType,betaChoice)),'-struct','AllCorr');
    case 'PLOT_splithalf_corr'
        % plotting function for split-half correlation - within or across
        % repetition (for trained / untrained)
        roi=[1:8];
        betaChoice='multiPW';
        parcelType='striatum_buckner_hemi';
        metric='corr'; % corr, corr_noMean, corr_mean
        vararginoptions(varargin,{'roi','metric','betaChoice'});
       % T = load(fullfile(repSupDir,sprintf('corr_splitHalf_NEW_exe_%s',betaChoice)));
        T = load(fullfile(repSupDir,sprintf('corr_splitHalf_%s_exe_%s',parcelType,betaChoice)));
        
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
    
    case 'PCM_PREDICT_2nd'
        roi=[1:8];
        regType='cortex';
        sessN=[1:4];
        sn=[4:9,11:25];
        seqType=[1,2];
        betaChoice='multiPW';
        vararginoptions(varargin,{'roi','regType','sessN','sn'});
        
        TT=[];   
        for ss=sessN
            switch regType
                case 'cortex'
                    D = load(fullfile(regDir,sprintf('betas_FoSEx_sess%d',ss)));
                    regLabel=regname_cortex;
                case 'striatum'
                    D = load(fullfile(regDir,sprintf('betas_BG-striatum_FoSEx_sess%d',ss)));
                    regLabel={'Caudate','Putamen'};
            end
            for st=seqType
                for r=roi
                    for s=sn
                        D1 = getrow(D,D.SN==s & D.region==r);
                        % load SPM info, extract 1st vs. 2nd exe data
                        I = load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM_info'));
                        for exe=1:2
                            condVec{exe}=I.seqNumb(I.seqType==st & I.FoSEx==exe);
                            if st==2
                                condVec{exe}=condVec{exe}-6; % make both trained / untrained start from 1
                            end
                            partVec{exe}=I.run(I.seqType==st & I.FoSEx==exe);
                            switch betaChoice % which betas to extract
                                case 'multiPW'
                                    data{exe}=D1.betaW{1}(I.seqType==st & I.FoSEx==exe,:);
                                case 'uniPW'
                                    data{exe}=D1.betaW{1}(I.seqType==st & I.FoSEx==exe,:);
                            end
                        end
                        % go across all voxels
                        for P = 1:size(data{1},2);
                            % G from 1st exe -> pcm model
                            % splithalf
                            sh = mod(partVec{1},2);
                            G1a = pcm_estGCrossval(data{1}(sh==1,P),partVec{1}(sh==1),condVec{1}(sh==1));
                            G1b = pcm_estGCrossval(data{1}(sh==0,P),partVec{1}(sh==0),condVec{1}(sh==0));
                            G2a = pcm_estGCrossval(data{2}(sh==1,P),partVec{2}(sh==1),condVec{2}(sh==1));
                            G2b = pcm_estGCrossval(data{2}(sh==0,P),partVec{2}(sh==0),condVec{2}(sh==0));
                            tmp1a = triu(G1a);
                            tmp1b = triu(G1b);
                            tmp2a = triu(G2a);
                            tmp2b = triu(G2b);
                            tmp1a(tmp1a==0)=NaN;
                            tmp1b(tmp1b==0)=NaN;
                            tmp2a(tmp2a==0)=NaN;
                            tmp2b(tmp2b==0)=NaN;
                            t1a = tmp1a(~isnan(tmp1a));
                            t1b = tmp1b(~isnan(tmp1b));
                            t2a = tmp2a(~isnan(tmp2a));
                            t2b = tmp1b(~isnan(tmp2b));
                            % correlation of first - second
                            cor_exe1(P)  = corrN(t1a,t1b);
                            cor_exe2(P)  = corrN(t2a,t2b);
                            cor_cross(P) = mean([corr(t1a(:),t2b(:)),corr(t1b(:),t2a(:))]);
                            % r2
                            %       SS_tot = sum((t2-mean(t2)).^2);
                            %       SS_res = sum((t1-t2).^2);
                            %       r2(P)  = 1 - (SS_res/SS_tot);
                        end
                        % make a structure
                        T.cor_exe1 = mean(cor_exe1);
                        T.cor_exe2 = mean(cor_exe2);
                        T.cor_cross = mean(cor_cross);
                        T.sn = s;
                        T.roi = r;
                        T.sessN = ss;
                        T.seqType = st;
                        TT = addstruct(TT,T);
                        fprintf('Done sess-%d: %s-%s\n',ss,regLabel{r},subj_name{s});
                    end; %sn
                end; %roi
            end; %seqType
        end; % sessN
        % save variable
        save(fullfile(repSupDir,sprintf('correlation_exe_%s',regType)),'-struct','TT');
    case 'PCM_PREDICT_2nd_plot'
        regType = 'cortex';
        roi = [1:8];
        vararginoptions(varargin,{'regType','roi'});
        T = load(fullfile(repSupDir,sprintf('correlation_exe_%s',regType)));
        
        switch regType
            case 'cortex'
                regLabel = regname_cortex;
            case 'striatum'
                regLabel = {'Caudate','Putamen'};
        end    
        figure
        for r=roi
            subplot(1,numel(roi),r)
            plt.line(T.sessN,T.cor_cross./ssqrt(T.cor_exe1.*T.cor_exe2),'subset',T.roi==r,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','northeast');
            title(sprintf('%s',regLabel{r}));
            if r==1
                xlabel('Session');
                ylabel('Correlation 1st-2nd');
            else
                ylabel('');
            end
            plt.match('y');
            drawline(0,'dir','horz');
        end
        
    case 'SAVE_psc_cortex_BG'   % ------------------------------  for CCN - CHECK ALL OF THAT + RE-DO
       betaChoice = 'multiPW';
       roi_cortex=[2,3];
       roi_BG=[1,2];
       parcelBG = 'BG-striatum';
       vararginoptions(varargin,{'betaChoice','roi_cortex','roi_BG'}); 
       
       % load cortex and BG
       C = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
       B = load(fullfile(BGDir,'stats',sprintf('RepSup_%s_stats_%s.mat',parcelBG,betaChoice))); 
       
       C=getrow(C,ismember(C.roi,roi_cortex));
       B=getrow(B,ismember(B.roi,roi_BG));
       
       C.reg=C.roi;
       C.reg(C.reg==2)=1;
       C.reg(C.reg==3)=2;
       B.reg=B.roi;
       B.reg(B.reg==2)=4;
       B.reg(B.reg==1)=3;
       
       S=[];
       S.psc=[C.psc_train;C.psc_untrain;B.psc_train;B.psc_untrain];
       S.seqType=[ones(size(C.psc_train));ones(size(C.psc_untrain))*2;ones(size(B.psc_untrain));ones(size(B.psc_untrain))*2];
       S.roi=[C.reg;C.reg;B.reg;B.reg];
       S.sessN=[C.sessN;C.sessN;B.sessN;B.sessN];
       S.FoSEx=[C.FoSEx;C.FoSEx;B.FoSEx;B.FoSEx];
       S.sn=[C.sn;C.sn;B.sn;B.sn];
       
       save(fullfile(repSupDir,'RepSup_cortex_BG_psc.mat'),'-struct','S');
    case 'SAVE_psc_cortex_BG_diff'
       roi=[1:4]; % 1 - M1, 2 - PMd, 3 - Caudate, 4 - Putamen
       sn=[4:9,11:25];
       sessN=[1:4];
       T=load(fullfile(repSupDir,'RepSup_cortex_BG_psc.mat'));
       
       TT=[];
         for ss=sessN
            for s=sn
                for st=1:2 % sequence type
                    for r=roi
                        t1=getrow(T,T.FoSEx==1&T.seqType==st&T.sn==s&T.sessN==ss&T.roi==r);
                        t2=getrow(T,T.FoSEx==2&T.seqType==st&T.sn==s&T.sessN==ss&T.roi==r);
                        diff=t1.psc-t2.psc;
                        
                        T2.psc=[t1.psc;t2.psc;diff];
                        T2.exeType=[1;2;3]; % first, second, difference
                        T2.seqType=[st;st;st];
                        T2.reg=[r;r;r];
                        T2.sn=[s;s;s];
                        T2.sessN=[ss;ss;ss];
                        
                        TT=addstruct(TT,T2);
                    end
                end
            end
         end
        
        save(fullfile(repSupDir,'RepSup_cortex_diff_BG_psc.mat'),'-struct','TT');
    case 'SAVE_psc_cortex_BG_ratio'
       roi=[1:4]; % 1 - M1, 2 - PMd, 3 - Caudate, 4 - Putamen
       sn=[4:9,11:25];
       sessN=[1:4];
       T=load(fullfile(repSupDir,'RepSup_cortex_BG_psc.mat'));
       
       TT=[];
       
         for ss=sessN
            for s=sn
                for st=1:2 % sequence type
                    for r=roi
                        t1=getrow(T,T.FoSEx==1&T.seqType==st&T.sn==s&T.sessN==ss&T.roi==r);
                        t2=getrow(T,T.FoSEx==2&T.seqType==st&T.sn==s&T.sessN==ss&T.roi==r);
                        diff=t2.psc/t1.psc;
                        
                        T2.psc=[t1.psc;t2.psc;diff];
                        T2.exeType=[1;2;3]; % first, second, ratio
                        T2.seqType=[st;st;st];
                        T2.reg=[r;r;r];
                        T2.sn=[s;s;s];
                        T2.sessN=[ss;ss;ss];
                        
                        TT=addstruct(TT,T2);
                    end
                end
            end
         end
        
        save(fullfile(repSupDir,'RepSup_cortex_ratio_BG_psc.mat'),'-struct','TT'); 
    case 'PLOT_psc_cortex_BG_sess1'
        roi=[1:4];
        reg_label={'M1','PMd','Caudate','Putamen'};
        
        T=load(fullfile(repSupDir,'RepSup_cortex_diff_BG_psc.mat'));
        
        figure
        for r=roi
            subplot(1,numel(roi),r)
            % exeType: 1-1st exe, 2-2nd exe, 3 - difference
            plt.bar([T.exeType>2 T.exeType],T.psc,'subset',T.sessN==1 & T.reg==r);
            plt.match('y');
            title(sprintf('%s',reg_label{r}));
            if r==1
                ylabel('Percent signal change');
            else
                ylabel('');
            end
        end
    case 'PLOT_psc_cortex_BG_sess4'
        roi=[1:4];
        reg_label={'M1','PMd','Caudate','Putamen'};
        
        T=load(fullfile(repSupDir,'RepSup_cortex_diff_BG_psc.mat'));
        
        T=getrow(T,(T.sessN==1 | T.sessN==4) & T.exeType==3);
        
        T.seqIndx=T.seqType;
        
        T.seqIndx(T.sessN==4)=T.seqIndx(T.sessN==4)+1;
        T.seqIndx(T.sessN==1)=1;

        figure
        for r=roi
            subplot(1,numel(roi),r)
            plt.bar([T.seqIndx>1 T.seqIndx],T.psc,'subset',T.reg==r);
            plt.match('y');
            title(sprintf('%s',reg_label{r}));
            if r==1
                ylabel('Percent signal change');
            else
                ylabel('');
            end
        end   
    case 'CALC_expectDist_2nd_dist'
        roi_cortex=[2,3];
        roi_BG=[1,2];
        subtract_mean=0;
        sessN=[1:4];
        sn=[4:9,11:25];
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'roi','subtract_mean'});
        
        Dc = load(fullfile(repSupDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        Pc = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s.mat',betaChoice)));
        
        Dbg = load(fullfile(BGDir,'stats',sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        Pbg = load(fullfile(BGDir,'stats',sprintf('RepSup_percent_%s.mat',betaChoice)));
        
        P = load(fullfile(repSupDir,'RepSup_cortex_ratio_BG_psc.mat'));
       
        CC=[]; BB=[];
        for ss=sessN
            for s=sn
                for st=1:2 % sequence type
                    for r=roi_cortex
                        % extract what needed for cortex
                       % Pc2=getrow(Pc,Pc.roi==r & Pc.sessN==ss & Pc.sn==s & Pc.seq==st);
                        if r==2
                            Pc2=getrow(P,P.reg==1 & P.sessN==ss & P.sn==s & P.seqType==st & P.exeType==3);
                        elseif r==3
                            Pc2=getrow(P,P.reg==2 & P.sessN==ss & P.sn==s & P.seqType==st & P.exeType==3);
                        end
                        Dc2=getrow(Dc,Dc.roi==r & Dc.sessN==ss & Dc.sn==s & Dc.seqType==st);
                        
                        C.dist1=Dc2.corrDist(Dc2.FoSEx==1);
                        C.dist2=Dc2.corrDist(Dc2.FoSEx==2);
                     %   C.expectDist=Pc2.psc_psc.*C.dist1;
                        C.expectDist=Pc2.psc.*C.dist1;
                        C.diffFromExpect=C.dist2-C.expectDist;
                        C.roi=r;

                        C.sn=s;
                        C.sessN=ss;
                        C.seqType=st;
                        CC=addstruct(CC,C);
                    end
                    
                    for r=roi_BG
                        % extract what needed for BG
                       % Pbg2=getrow(Pbg,Pbg.roi==r & Pbg.sessN==ss & Pbg.sn==s & Pbg.seq==st);
                       if r==1
                            Pbg2=getrow(P,P.reg==3 & P.sessN==ss & P.sn==s & P.seqType==st & P.exeType==3);
                        elseif r==2
                            Pbg2=getrow(P,P.reg==4 & P.sessN==ss & P.sn==s & P.seqType==st & P.exeType==3);
                        end
                        Dbg2=getrow(Dbg,Dbg.roi==r & Dbg.sessN==ss & Dbg.sn==s & Dbg.seqType==st);
                        
                        B.dist1=Dbg2.corrDist(Dbg2.FoSEx==1);
                        B.dist2=Dbg2.corrDist(Dbg2.FoSEx==2);
                        %B.expectDist=Pbg2.psc_psc.*B.dist1;
                        B.expectDist=Pbg2.psc.*B.dist1;
                        B.diffFromExpect=B.dist2-B.expectDist;
                        B.roi=r;
                        B.sn=s;
                        B.sessN=ss;
                        B.seqType=st;
                        BB=addstruct(BB,B);
                    end
                    
                end
            end
        end
          
        CC.reg=CC.roi;
        CC.reg(CC.reg==2)=1;
        CC.reg(CC.reg==3)=2;
        BB.reg=BB.roi;
        BB.reg(BB.reg==2)=4;
        BB.reg(BB.reg==1)=3;
        
        TT=[];TT=addstruct(TT,CC);TT=addstruct(TT,BB);
        
        save(fullfile(repSupDir,'RepSup_expectDist_cortex_BG.mat'),'-struct','TT');
    case 'PLOT_expectDist_2nd_sess1'
        roi=[1:4];
        reg_label={'M1','PMd','Caudate','Putamen'};
        sessN=1;
        vararginoptions(varargin,{'roi','sessN'});
        
        D=load(fullfile(repSupDir,'RepSup_expectDist_cortex_BG.mat'));
      
        Dd.dist=[D.dist1;D.dist2;D.diffFromExpect];
        Dd.exeType=[ones(size(D.dist1));ones(size(D.dist2))*2;ones(size(D.dist2))*3];
        % 1- first exe, 2 - second exe, 3 - difference from expected
        Dd.reg=[D.reg;D.reg;D.reg];
        Dd.sessN=[D.sessN;D.sessN;D.sessN];
        Dd.sn=[D.sn;D.sn;D.sn];
        Dd.seqType=[D.seqType;D.seqType;D.seqType];
        
        for ss=sessN
            figure
            sub=1;
            for r=roi
                subplot(1,numel(unique(roi)),sub)
                if r==4
                    plt.bar([Dd.exeType>2 Dd.exeType],Dd.dist,'split',Dd.exeType,'subset',Dd.reg==r&Dd.sessN==ss&Dd.seqType==2,...
                        'leg',{'1st distance','2nd distance'}, 'leglocation','northeast','style',stySess);
                else
                    plt.bar([Dd.exeType>2 Dd.exeType],Dd.dist,'split',Dd.exeType,'subset',Dd.reg==r&Dd.sessN==ss,...
                        'leg',{'1st distance','2nd distance'}, 'leglocation','northeast','style',stySess);
                end
                meanLine=mean(D.expectDist(D.reg==r&D.sessN==ss));
                drawline(meanLine,'dir','horz','linestyle','--');
                drawline(0,'dir','horz');
                xlabel('Execution'); title(sprintf('%s',reg_label{r}));
                if sub==1
                    ylabel('Distance');
                else
                    ylabel('');
                end
                sub=sub+1;
                plt.match('y');
            end
        end
    case 'PLOT_expectDist_2nd_seqType'
        roi=[1:4];
        reg_label={'M1','PMd','Caudate','Putamen'};
        sessN=1;
        vararginoptions(varargin,{'roi','sessN'});
        
        D=load(fullfile(repSupDir,'RepSup_expectDist_cortex_BG.mat'));
      
        Dd.dist=[D.dist1;D.dist2];
        Dd.exeType=[ones(size(D.dist1));ones(size(D.dist2))*2];
        % 1- first exe, 2 - second exe, 3 - difference from expected
        Dd.reg=[D.reg;D.reg];
        Dd.sessN=[D.sessN;D.sessN];
        Dd.sn=[D.sn;D.sn];
        Dd.seqType=[D.seqType;D.seqType];
        
        for ss=sessN
            figure
            sub=1;
            for r=roi
                subplot(1,numel(unique(roi)),sub)
                if r==4
                    plt.bar(Dd.seqType,Dd.dist,'split',Dd.exeType,'subset',Dd.reg==r&Dd.sessN==ss,...
                        'leg',{'1st distance','2nd distance'}, 'leglocation','northeast','style',stySess);
                else
                    plt.bar(Dd.seqType,Dd.dist,'split',Dd.exeType,'subset',Dd.reg==r&Dd.sessN==ss,...
                        'leg',{'1st distance','2nd distance'}, 'leglocation','northeast','style',stySess);
                end
                meanLine=mean(D.expectDist(D.reg==r&D.sessN==ss));
                drawline(meanLine,'dir','horz','linestyle','--');
                drawline(0,'dir','horz');
                xlabel('Sequence type'); title(sprintf('%s',reg_label{r}));
                if sub==1
                    ylabel('Distance');
                else
                    ylabel('');
                end
                sub=sub+1;
               % plt.match('y');
            end
        end
    case 'PLOT_expectedDist_2nd_sess4'
        roi=[1:4];
        reg_label={'M1','PMd','Caudate','Putamen'};
        sessN=[1:3];
        vararginoptions(varargin,{'roi','sessN'});
        
        D=load(fullfile(repSupDir,'RepSup_expectDist_cortex_BG.mat'));
        
        T=getrow(D,(D.sessN==1 | D.sessN==4));
        
        T.seqIndx=T.seqType;
        
        T.seqIndx(T.sessN==4)=T.seqIndx(T.sessN==4)+1;
        T.seqIndx(T.sessN==1)=1;
        
        figure
        for r=roi
            subplot(1,numel(unique(roi)),r)
            plt.bar([T.seqIndx>1 T.seqIndx],T.diffFromExpect,'subset',T.reg==r);
            drawline(0,'dir','horz');
            
            xlabel('SeqType'); title(sprintf('%s',reg_label{r}));
            if r==1
                ylabel('Distance');
            else
                ylabel('');
            end    
        end    
    case 'PLOT_expectDist_percent'
        roi=[1:5,7,8];
        vararginoptions(varargin,{'roi'});
        
        D=load(fullfile(repSupDir,'expectedDist_RS.mat'));
        
        D=getrow(D,D.indx>1);
        
        expDist=D.dist(D.indx==2);
        realDist=D.dist(D.indx==3);
        
        P.ratio=realDist-expDist;
        P.roi=D.roi(D.indx==2);
        P.sessN=D.sessN(D.indx==2);
        P.seqType=D.seqType(D.indx==2);
        P.sn=D.sn(D.indx==2);
        
        figure
        sub=1;
        for r=roi
            subplot(1,numel(unique(roi)),sub)
            plt.line(P.sessN,P.ratio,'split',P.seqType,'subset',P.roi==r,'leg',{'trained','untrained'},...
                'leglocation','northeast','style',stySeq);
            xlabel('Session'); title(sprintf('%s',regname{r}));
            if sub==1
                ylabel('Distance'); 
            else
                ylabel('');
            end
            sub=sub+1;
        end     
             
    case 'RDM_connect'
        sn=[4:9,11:25];
        sessN=[1:4];
        reg_cortex=[1:16];
        reg_BG=[1:4];
        parcelBG='BG-striatum';
        betaChoice='multiPW';
        vararginoptions(varargin,{'sn','sessN','reg_cortex','reg_BG','parcelBG','betaChoice'});
        
        regAll  = [reg_cortex reg_BG];
        numReg  = length(regAll);     
        regType = [ones(size(reg_cortex)) ones(size(reg_BG))*2];
        KK=[];
        
        for ss=sessN
            C = load(fullfile(regDir,sprintf('stats_FoSEx_%s_sess%d',betaChoice,ss)));
            B = load(fullfile(BGDir,'stats',sprintf('stats_FoSEx_%s_%s_sess%d',parcelBG,betaChoice,ss)));
            for st=1:2 % seqType
                for exe=1:2 % execution
                    for s=sn
                        % get rdms from all regions
                        for r=1:numReg
                            if regType(r) == 1
                                T = getrow(C,C.SN==s & C.FoSEx==exe & C.region==regAll(r));
                            else
                                T = getrow(B,B.SN==s & B.FoSEx==exe & B.region==regAll(r));
                            end
                            tmp = rsa_squareRDM(T.RDM);
                            % get the right entries for trained / untrained
                            if st==1
                                R{r} = tmp(1:6,1:6);
                            else
                                R{r} = tmp(7:12,7:12);
                            end
                        end
                        % perform calculations
                        K.rdmCorr = rsa_calcCorrRDMs(R);
                        K.distCorr = rsa_calcDistCorrRDMs(R);
                        K.sn = s;
                        K.seqType = st;
                        K.exe = exe;
                        K.sessN = ss;
                        KK=addstruct(KK,K);
                        fprintf('Done sess-%d exe-%d: %s\n',ss,exe,subj_name{s});
                    end; % sn
                end; % exe
            end; % seqType
        end; % session
        save(fullfile(repSupDir,sprintf('connect_RDM_repsup_%s',betaChoice)),'-struct','KK');
    case 'CONN_perHemi'
        connType = 'RDM_repsup_multiPW';
        var      = {'rdmCorr','distCorr'}; 
        % var for RDM_multiPW / RDM_uniPW: {'rdmCorr','distCorr'}
        vararginoptions(varargin,{'connType','var'});
        T = load(fullfile(repSupDir,sprintf('connect_%s',connType)));
        
        regType=[1:8,17,18;9:16,19,20]; % regType per hemi
        I = indicatorMatrix('allpairs',[1:20]);
        for h=1:2
            % get the right pairs of reg per hemisphere
            ind=[];
            for j=1:size(I,1)
                [r,c] = find(I(j,:));
                if sum(ismember(c,regType(h,:)))==2 % both regions in I
                    ind = [ind;j];
                end
            end
            H = T;
            for v=1:size(var,2)
                H.(var{v}) = T.(var{v})(:,ind);
            end
            % save the new structure
            save(fullfile(repSupDir,sprintf('connect_%s_%s',connType,hem{h})),'-struct','H');
        end
    case 'RDM_connect_plot'
        regName = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp','Caudate','Putamen'};
        seqType = 1;
        sessN=[1:4];
        betaChoice='multiPW';
        hemi = 'lh';
        seqLabel = {'trained','untrained'};
        vararginoptions(varargin,{'seqType','sessN','betaChoice','hemi'});
         T = load(fullfile(repSupDir,sprintf('connect_RDM_repsup_%s_%s',betaChoice,hemi)));
        
         for st=seqType
             for exe=1:2
                 figure
                 indx=1;
                 for dt=1:2
                     if dt==1
                         metric='rdmCorr';
                     else
                         metric='distCorr';
                     end
                     for ss=sessN
                         K = getrow(T,T.seqType==st & T.sessN==ss & T.exe==exe);
                         subplot(2,numel(sessN),indx);
                         imagesc(rsa_squareRDM(mean(K.(metric),1)));
                         caxis([-0.05 0.8]);
                         set(gca,'XTickLabel',regName,'YTickLabel',regName);
                         title(sprintf('sess-%d %s exe-%d %s',ss,metric,exe,seqLabel{st}));
                         indx=indx+1;
                     end
                 end
             end
         end
    case 'RDM_connect_diff'
        regName = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp','Caudate','Putamen'};
        seqType = 1;
        sessN=[1:4];
        betaChoice='multiPW';
        seqLabel = {'trained','untrained'};
        vararginoptions(varargin,{'seqType','sessN','betaChoice'});
         T = load(fullfile(repSupDir,sprintf('RDM_connect_%s',betaChoice)));
        
         for st=seqType
             indx=1;
             figure
             for dt=1:2
                 if dt==1
                     metric='rdmCorr';
                 else
                     metric='distCorr';
                 end
                 for ss=sessN
                     K1 = getrow(T,T.seqType==st & T.sessN==ss & T.exe==1);
                     K2 = getrow(T,T.seqType==st & T.sessN==ss & T.exe==2);
                     diff = K1.(metric)-K2.(metric);
                     subplot(2,numel(sessN),indx);
                     imagesc(rsa_squareRDM(mean(diff,1)));
                     caxis([-0.15 0.2]);
                     set(gca,'XTickLabel',regName,'YTickLabel',regName);
                     title(sprintf('sess-%d %s %s difference',ss,metric,seqLabel{st}));
                     indx=indx+1;
                 end
             end
         end
    case 'RDM_lineplot'
        seqType = 1;
        sessN=[1:4];
        betaChoice='multiPW';
        vararginoptions(varargin,{'seqType','sessN','betaChoice'});
        T = load(fullfile(repSupDir,sprintf('RDM_connect_%s',betaChoice)));
        
        
        for st=seqType
            figure
            indx=1;
            for dt=1:2
                if dt==1
                    metric='rdmCorr';
                else
                    metric='distCorr';
                end
                
                for ss=sessN
                    K1 = getrow(T,T.seqType==st & T.sessN==ss & T.exe==1);
                    K2 = getrow(T,T.seqType==st & T.sessN==ss & T.exe==2);
                    subplot(2,numel(sessN),indx);
                    plot(mean(K1.(metric),1),'k-','LineWidth',2);
                    hold on;
                    plot(mean(K1.(metric),1),'ko');
                    plot(mean(K2.(metric),1),'r-','LineWidth',2);
                    plot(mean(K2.(metric),1),'ro');
                    title(sprintf('sess-%d %s',ss,metric));
                    indx=indx+1;        
                end
            end
        end
        
    case 'ROI_timeseries'
        % to check the model quality of the glm
        glm         = 2;
        parcelType  = 'cortex';
        sessN       = 4; 
        sn          = [4:9,11:25];
        reg         = 2;
        vararginoptions(varargin,{'sn','sessN','glm','parcelType','sessN','reg'});
        
        pre     = 4;          % How many TRs before the trial onset
        post    = 20;        % How many TRs after the trial onset
        for ss=sessN
            for s=sn
                T=[];
                load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM.mat'));
                
                SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s}));
                switch parcelType
                    case 'cortex'
                        load(fullfile(regDir,sprintf('%s_Brodmann_regions.mat',subj_name{s})));      % This is made in case 'ROI_define'
                    case 'striatum'
                        load(fullfile(regDir,sprintf('%s_BG-striatum_regions.mat',subj_name{s})));      % This is made in case 'ROI_define'
                    case 'striatum_networks'
                        load(fullfile(regDir,sprintf('%s_striatum_buckner_hemi_regions.mat',subj_name{s})));   
                    case 'cortex_networks'
                        load(fullfile(regDir,sprintf('%s_cortex_buckner_regions.mat',subj_name{s})));   
                end
                Y = region_getdata(SPM.xY.VY,R);
                % Create a structure with trial onset and trial type (event)
                D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
                % extract only the first execution
                DD = getrow(D,D.event<13); % 1-6 trained, 7-12 untrained
                
                % per region extract data, get timeseries
                for r=1:numel(reg)
                    % get timeseries from region for all voxels
                    y_raw=Y{reg(r)};
                    % get the design matrix
                    reg_interest=[SPM.xX.iH SPM.xX.iC];
                    % filter and get adjusted timeseries
                    if (~isfield(SPM.xX,'pKX')) % SPM not run -
                        y_filt = spm_filter(SPM.xX.K,y_raw);
                        SPM.xX.xKXs.X = spm_filter(SPM.xX.K,SPM.xX.X);
                        SPM.xX.pKX = pinv(SPM.xX.xKXs.X);
                        B = SPM.xX.pKX*y_filt;                              %-Parameter estimates
                        y_res   = y_filt - SPM.xX.xKXs.X*B;             %-Residuals
                        y_hat = SPM.xX.xKXs.X(:,reg_interest)*B(reg_interest,:); %- predicted values
                    else
                        y_filt = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
                        B = SPM.xX.pKX*y_filt;                              %-Parameter estimates
                        y_res   = spm_sp('r',SPM.xX.xKXs,y_filt);             %-Residuals
                        y_hat = SPM.xX.xKXs.X(:,reg_interest)*B(reg_interest,:); %- predicted values
                    end;
                    y_adj  = y_hat + y_res;
                    
                    % prewhiten
                    res = ssqrt(sum(y_res.^2,1));
                    y_adj_pw = bsxfun(@rdivide,y_adj,res);
                    y_hat_pw = bsxfun(@rdivide,y_hat,res);
                    % extract timeseries per voxel
                    nVox = size(y_adj,2);
                    for i=1:(size(DD.block,1));
                        for v=1:nVox
                            % voxel x timepoint saved in structure 
                            % currently can only save one region per
                            % subject / structure because of 3D nature
                            % i - per each run / seqType etc.
                            Y_adj(i,v,:)    = cut(y_adj(:,v),pre,round(DD.ons(i)),post,'padding','nan')';
                            Y_hat(i,v,:)    = cut(y_hat(:,v),pre,round(DD.ons(i)),post,'padding','nan')';
                            Y_res(i,v,:)    = cut(y_res(:,v),pre,round(DD.ons(i)),post,'padding','nan')';
                            Y_raw(i,v,:)    = cut(y_raw(:,v),pre,round(DD.ons(i)),post,'padding','nan')';
                            Y_adj_pw(i,v,:) = cut(y_adj_pw(:,v),pre,round(DD.ons(i)),post,'padding','nan')';
                            Y_hat_pw(i,v,:) = cut(y_hat_pw(:,v),pre,round(DD.ons(i)),post,'padding','nan')';
                        end;
                    end
                    % make the structure per timepoint
                    S.y_adj=Y_adj;
                    S.y_hat=Y_hat;
                    S.y_res=Y_res;
                    S.y_raw=Y_raw;
                    S.y_adj_pw=Y_adj_pw;
                    S.y_hat_pw=Y_hat_pw;
                    S.run=DD.block;
                    S.seqNum=DD.event;
                    % indicate trained / untrained
                    S.seqType=zeros(size(S.run));
                    S.seqType(S.seqNum<7)=1;
                    S.seqType(S.seqNum>6)=2;
                    S.sn=ones(size(S.run))*s;
                    S.region=ones(size(S.run))*reg(r); % region
                    S.regType=regType(S.region)';
                    S.regSide=regSide(S.region)';
                    clear Y_adj Y_hat Y_res Y_raw Y_adj_pw Y_hat_pw;
                end
                fprintf('Done timeseries extraction: \t%s-region-%d\n',subj_name{s},r);
                cd(repSupDir);
                save(sprintf('hrf_%s_%s_sess%d_FoSEx.mat',parcelType,subj_name{s},ss),'-struct','S');
            end
            
        end;
    case 'ROI_plot_timeseries_sn'                
        sn    = 4; 
        sessN = 4;
        reg   = 2;
        regS  = 1; % 1 - LH, 2 - RH
        parcelType = 'cortex'; % cortex or striatum

        vararginoptions(varargin,{'sn','sessN','glm','reg','regS','parcelType'});
        
        for s=sn
            T=load(fullfile(repSupDir,sprintf('hrf_%s_%s_sess%d_FoSEx.mat',parcelType,subj_name{s},sessN)));
            
            figure
            traceplot([-4:20],nanmean(T.y_hat,2),'errorfcn','stderr','subset',T.regSide==regS & T.regType==reg,'split',T.seqType,'leg','auto');
            hold on;
           % traceplot([-4:20],nanmean(T.y_hat,2),'subset',T.regSide==regS & T.regType==reg,'linestyle',':','split',T.seqType);
           % hold off;
            xlabel('TR');
            ylabel('activation');
            drawline(0,'dir','horz'); drawline(0,'dir','vert');
        end
    case 'ROI_plot_timeseries_group'     
        sn    = [4:9,11:25]; 
        sessN = 4;
        reg   = 2;
        regS  = 1; % 1 - LH, 2 - RH
        parcelType = 'cortex'; % cortex or striatum

        vararginoptions(varargin,{'sn','sessN','glm','reg','regS','parcelType'});
        SS=[];
        for s=sn
            T = load(fullfile(repSupDir,sprintf('hrf_%s_%s_sess%d_FoSEx.mat',parcelType,subj_name{s},sessN)));
            S.data = squeeze(mean(T.y_hat_pw,2));
            S.regSide = T.regSide;
            S.regType = T.regType;
            S.seqType = T.seqType;
            S.sn      = ones(size(T.seqType))*s;
            SS=addstruct(SS,S);
        end
        SS.seqType=SS.seqType+1;
        SS.seqType(SS.seqType==3)=1;
        figure
        traceplot([-4:20],SS.data,'errorfcn','stderr','split',SS.seqType,'leg',{'untrained','trained'});
        hold on;
        xlabel('TR');
        ylabel('activation');
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline([5,10,15],'dir','vert','linestyle','--');
        % plot across subjects
    case 'ROI_ts_distances'
        sn    = [4:9,11:25]; 
        sessN = 4;
        reg   = 2;
        regS  = 1; % 1 - LH, 2 - RH
        parcelType = 'cortex'; % cortex or striatum

        vararginoptions(varargin,{'sn','sessN','glm','reg','regS','parcelType'});
        SS=[];
  
        for s=sn
            T = load(fullfile(repSupDir,sprintf('hrf_%s_%s_sess%d_FoSEx.mat',parcelType,subj_name{s},sessN)));
            % per timepoint
            for t=1:size(T.y_adj,3)
                rowNan = find(any(isnan(T.y_hat_pw(:,:,1)),2));
                if ~isempty(rowNan)
                    T = getrow(T,[1:rowNan-1 rowNan+1:size(T.sn,1)]); % exlude nan trials
                end
                dist = rsa_squareRDM(rsa.distanceLDC(T.y_adj_pw(:,:,t),T.run,T.seqNum));
                dist_train = triu(dist(1:6,1:6));
                dist_untrain = triu(dist(7:12,7:12));
                dist_cross = triu(dist(1:6,7:12));
                dist_train(dist_train==0)=NaN;
                dist_untrain(dist_untrain==0)=NaN;
               % dist_cross(dist_cross==0)=NaN;
                S.dist(1,:) = nanmean(dist_train(:));
                S.dist(2,:) = nanmean(dist_untrain(:));
                S.seqType   = [1;2];
              %  S.distCross = nanmean(dist_cross(:));
                S.timepoint = [t-5;t-5];
                S.sn        = [s;s];
                SS=addstruct(SS,S);               
            end
        end

        figure
        lineplot(SS.timepoint,SS.dist,'style_shade','split',SS.seqType,'leg',{'trained','untrained'});
        xlabel('TR');
        ylabel('distances');
        drawline(0,'dir','horz'); drawline(0,'dir','vert');
        drawline([5,10,15],'dir','vert','linestyle','--');
      
    case 'CALC_scaling'
        parcelType = 'cortex_buckner';
        roi=[1:7];
        sn=[4:9,11:25];
        sessN=4;
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'sn','roi','sessN','parcelType','betaChoice'});
        
        TT=[];
        for ss=sessN
            switch parcelType
                case 'cortex_buckner'
                    D = load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract0.mat',parcelType))); % also dist - Mahalanobis
                    P = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
                case 'striatum_buckner_hemi'
                    P = load(fullfile(BGDir,'stats',sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice))); 
                    D = P;
            end
            % calculate per subject, roi, seqType
            for s=sn
                for r=roi
                    D2=getrow(D,D.sn==s & D.roi==r & D.sessN==ss);
                    P2=getrow(P,P.sn==s & P.roi==r & P.sessN==ss);
                    for st=1:2
                        % calculate scaling of percent
                        if st==1
                            var1 = 'psc_train';
                            var2 = 'search_train';
                        else
                            var1 = 'psc_untrain';
                            var2 = 'search_untrain';
                        end
                      
                        if ~size(P2.FoSEx,1)==0
                            T.psc_scale  = P2.(var1)(P2.FoSEx==2)/P2.(var1)(P2.FoSEx==1);
                            switch parcelType
                                case 'cortex_buckner'
                                    tmp = getrow(D2,D2.seqType==st);
                                    T.dist_scale = tmp.dist(tmp.FoSEx==2)/tmp.dist(tmp.FoSEx==1);
                                case 'striatum_buckner_hemi'
                                    T.dist_scale = P2.(var2)(P2.FoSEx==2)/P2.(var2)(P2.FoSEx==1);         
                            end                
                            T.roi        = r;
                            T.sn         = s;
                            T.sessN      = ss;
                            T.seqType    = st;
                            TT=addstruct(TT,T);
                        end
                    end
                end
            end
        end
        figure
        for r=roi
            subplot(numel(roi),2,(r-1)*2+1)
            plt.bar(TT.seqType,TT.dist_scale,'subset',TT.roi==r);
            title(sprintf('network-%d dist scale',r));
            set(gca,'XTickLabel',{'trained','untrained'});
            ylabel('ratio exe2/exe1');
            subplot(numel(roi),2,(r-1)*2+2)
            plt.bar(TT.seqType,TT.psc_scale,'subset',TT.roi==r);
            title(sprintf('network-%d psc scale',r));
            set(gca,'XTickLabel',{'trained','untrained'});
            ylabel('');
           % hold on;
           % drawline(mean(TT.psc_scale(TT.seqType==1 & TT.roi==r)),'lim',[0,2],'dir','horz');
           % drawline(mean(TT.psc_scale(TT.seqType==2 & TT.roi==r)),'lim',[2,4],'dir','horz');
        end

        figure
        for r=roi
            subplot(1,numel(roi),r)
            plt.bar(TT.seqType,TT.dist_scale,'subset',TT.roi==r);
            title(sprintf('network-%d dist scale',r));
            set(gca,'XTickLabel',{'trained','untrained'});
            ylabel('ratio exe2/exe1');
            hold on;
            drawline(mean(TT.psc_scale(TT.seqType==1 & TT.roi==r)),'lim',[0,2],'dir','horz');
            drawline(mean(TT.psc_scale(TT.seqType==2 & TT.roi==r)),'lim',[2,4],'dir','horz');
        end
        
        keyboard;
        figure
        subplot(2,1,1)
        plt.bar(TT.roi,TT.dist_scale,'split',TT.seqType);
        subplot(2,1,2)
        plt.bar(TT.roi,TT.psc_scale,'split',TT.seqType);
        keyboard;
    case 'SCALE_striatum'
        parcelType = 'striatum_buckner_hemi';
        roi=[1:7];
        sessN=4;
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'sn','roi','sessN','parcelType','betaChoice'});
        
        P = load(fullfile(BGDir,'stats',sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice))); 
        
        for ss=sessN
            figure
            for r=1:numel(roi)
                P2 = getrow(P,P.sessN==ss&P.roi==roi(r));                  
             %   scaleT = mean(P2.psc_train(P2.FoSEx==2)-P2.psc_train(P2.FoSEx==1))/mean(P2.psc_train(P2.FoSEx==1));
             %   scaleU = mean(P2.psc_untrain(P2.FoSEx==2)-P2.psc_untrain(P2.FoSEx==1))/mean(P2.psc_untrain(P2.FoSEx==1));
                scaleT = mean(P2.act_vecLength_train(P2.FoSEx==2)-P2.act_vecLength_train(P2.FoSEx==1))/mean(P2.act_vecLength_train(P2.FoSEx==1));
                scaleU = mean(P2.act_vecLength_untrain(P2.FoSEx==2)-P2.act_vecLength_untrain(P2.FoSEx==1))/mean(P2.act_vecLength_untrain(P2.FoSEx==1));
                % calculate expected second distance
                expDistT = P2.search_train(P2.FoSEx==1)-abs(scaleT)*P2.search_train(P2.FoSEx==1);
                expDistU = P2.search_untrain(P2.FoSEx==1)-abs(scaleU)*P2.search_untrain(P2.FoSEx==1);
                
                % do the stats
                tmp = getrow(P2,P2.FoSEx==2);
                tmp.distT = tmp.search_train-expDistT;
                tmp.distU = tmp.search_untrain-expDistU;
                T.dist = [tmp.distT;tmp.distU];
                T.sn = [tmp.sn;tmp.sn];
                T.seqType = [ones(size(tmp.sn));ones(size(tmp.sn))*2];
                anovaMixed(T.dist,T.sn,'within',T.seqType,{'seqType'});
                ttestDirect(T.dist,[T.sn],2,'onesample','subset',T.seqType==1);
        
                subplot(1,numel(roi),r);
                %plt.bar(P2.FoSEx,P2.search_train);
                %plt.bar([P2.FoSEx;P2.FoSEx],[P2.search_train;P2.search_untrain],'split',[ones(size(P2.search_train));ones(size(P2.search_train))*2]);
                plt.bar([ones(size(P2.search_train));ones(size(P2.search_train))*2],[P2.search_train;P2.search_untrain],'split',[P2.FoSEx;P2.FoSEx],...
                    'style',stySeqType_exe,'leg',{'exe1','exe2'});
                drawline(mean(expDistT),'dir','horz','lim',[1.6,2.9],'linewidth',3,'color',[0.3 0.3 0.3]);
                drawline(mean(expDistU),'dir','horz','lim',[4.8,6],'linewidth',3,'color',[0.3 0.3 0.3]);
                set(gca,'XTick',[1.5,5],'XTickLabel',{'trained','untrained'});
                title(sprintf('Network-%d',r));
                if r==1
                    ylabel('Distance striatum');
                else
                    ylabel('');
                end
            end
        end
    case 'SCALE_cortex'
        parcelType = 'cortex_buckner';
        roi=[1:7];
        sessN=4;
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'sn','roi','sessN','parcelType','betaChoice'});
        
      %   D = load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract0.mat',parcelType))); % also dist - Mahalanobis
      %   P = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
         D = load(fullfile(repSupDir,sprintf('corrDist_RS_meanSubtract0.mat'))); % also dist - Mahalanobis
         P = load(fullfile(repSupDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        
        for ss=sessN
            figure
            for r=1:numel(roi)
                P2 = getrow(P,P.sessN==ss&P.roi==roi(r));        
                D2 = getrow(D,D.sessN==ss&D.roi==roi(r));     
            %    scaleT = mean(P2.psc_train(P2.FoSEx==2)-P2.psc_train(P2.FoSEx==1))/mean(P2.psc_train(P2.FoSEx==1));
            %    scaleU = mean(P2.psc_untrain(P2.FoSEx==2)-P2.psc_untrain(P2.FoSEx==1))/mean(P2.psc_untrain(P2.FoSEx==1));
                scaleT = mean(P2.act_vecLength_train(P2.FoSEx==2)-P2.act_vecLength_train(P2.FoSEx==1))/mean(P2.act_vecLength_train(P2.FoSEx==1));
                scaleU = mean(P2.act_vecLength_untrain(P2.FoSEx==2)-P2.act_vecLength_untrain(P2.FoSEx==1))/mean(P2.act_vecLength_untrain(P2.FoSEx==1));
                % calculate expected second distance
                %expDist = P2.dist_train(P2.FoSEx==1)-abs(scaleTr)*P2.dist_train(P2.FoSEx==1);
                expDistT = D2.dist(D2.FoSEx==1&D2.seqType==1)-abs(scaleT)*D2.dist(D2.FoSEx==1&D2.seqType==1);
                expDistU = D2.dist(D2.FoSEx==1&D2.seqType==2)-abs(scaleU)*D2.dist(D2.FoSEx==1&D2.seqType==2);
                         
                % do the stats
                tmp = getrow(D2,D2.FoSEx==2);
                distT = tmp.dist(tmp.seqType==1)-expDistT;
                distU = tmp.dist(tmp.seqType==2)-expDistU;
                T.dist = [distT;distU];
                T.sn = [tmp.sn(tmp.seqType==1);tmp.sn(tmp.seqType==2)];
                T.seqType = [ones(size(distT));ones(size(distT))*2];
                anovaMixed(T.dist,T.sn,'within',T.seqType,{'seqType'});
                ttestDirect(T.dist,[T.sn],2,'onesample','subset',T.seqType==1);
                
                
                subplot(1,numel(roi),r);
                %plt.bar(P2.FoSEx,P2.dist_train);
                %plt.bar(D2.FoSEx,D2.dist,'subset',D2.seqType==1);
                %plt.bar(D2.FoSEx,D2.dist,'split',D2.seqType);
                plt.bar(D2.seqType,D2.dist,'split',D2.FoSEx,'style',stySeqType_exe,'leg',{'exe1','exe2'});
                set(gca,'XTick',[1.5,5],'XTickLabel',{'trained','untrained'});
                drawline(mean(expDistT),'dir','horz','lim',[1.6,2.9],'linewidth',3,'color',[0.3 0.3 0.3]);
                drawline(mean(expDistU),'dir','horz','lim',[4.8,6],'linewidth',3,'color',[0.3 0.3 0.3]);
              %  title(sprintf('Network-%d',r));
                title(sprintf('%s',regname{r}));
                if r==1
                    ylabel('Distance cortex');
                else
                    ylabel('');
                end
            end
        end
    case 'SCALE_surface'
        sessN=4;
        sn = [4:14,17:22,24];
        structure = 'cortex'; % cortex or striatum
        vararginoptions(varargin,{'sessN','sn','structure'});
        for ss=sessN
            for h=1:2 % hemisphere
                for st=1:2
                    if st==1
                        lab = {'TrainSeq','trained'};
                    else
                        lab = {'UntrainSeq','untrained'};
                    end
                    switch structure
                        case 'cortex'
                            surfaceGDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
                            p1 = caret_load(fullfile(surfaceGDir,sprintf('%s.RS_%s_1st_sess%d.metric',hem{h},lab{1},ss)));
                            p2 = caret_load(fullfile(surfaceGDir,sprintf('%s.RS_%s_2nd_sess%d.metric',hem{h},lab{1},ss)));
                            d1 = caret_load(fullfile(surfaceGDir,sprintf('%s.exe1_dist_%s_sess%d.metric',hem{h},lab{2},ss)));
                            d2 = caret_load(fullfile(surfaceGDir,sprintf('%s.exe2_dist_%s_sess%d.metric',hem{h},lab{2},ss)));
                        case 'striatum'
                            surfaceGDir = fullfile(BGDir,'surface','group',hemName{h});
                            p1 = caret_load(fullfile(surfaceGDir,sprintf('%s.psc_RS_%s_1st_sess%d.metric',hem{h},lab{1},ss)));
                            p2 = caret_load(fullfile(surfaceGDir,sprintf('%s.psc_RS_%s_2nd_sess%d.metric',hem{h},lab{1},ss)));
                            d1 = caret_load(fullfile(surfaceGDir,sprintf('%s.dist_RS_%s_exe1_sess%d.metric',hem{h},lab{2},ss)));
                            d2 = caret_load(fullfile(surfaceGDir,sprintf('%s.dist_RS_%s_exe2_sess%d.metric',hem{h},lab{2},ss)));
                    end
                    nVert = size(p1.data,1);
                    for v=1:nVert
                        % difference in psc (dP) and dist (dD) for exe 1-2 at each vertex
                        scale = mean(p2.data(v,:)-p1.data(v,:))/mean(p2.data(v,:));
                        % calculate expected second distance
                        expDist = d1.data(v,:)-abs(scale)*d1.data(v,:);
                        % deviation from expected
                        % + if higher than expected, - if lower
                        % average across subjects
                        dev(v,:) = nanmean(d2.data(v,:)-expDist);
                    end
                    colName = sprintf('distDeviation_fromScaling_%s_sess-%d',lab{2},ss);
                    C = caret_struct('metric','data',dev,'column_name',{colName});
                    C.index=p1.index;
                    C.data=dev;
                    caret_save(fullfile(surfaceGDir, sprintf('%s.RS_distDeviation_%s_sess%d.metric',hem{h},lab{2},ss)),C);
                end
            end
        end
    case 'SCALE_surface_smooth'
        sessN=4;
        structure = 'cortex';
        nIter = 5; % 20 for cortex, 3 for striatum
        vararginoptions(varargin,{'sessN','structure','nIter'});
        seqName = {'trained','untrained'};
        for h=1:2       
            switch structure
                case 'cortex'
                    surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                    coordfile=[hem{h} '.WHITE.coord'];
                    topofile=[hem{h} '.CLOSED.topo'];
                case 'striatum'
                    surfaceGroupDir = fullfile(BGDir,'surface','group',hemName{h});
                    coordfile=[hem{h} '.striatum.coord.gii'];
                    topofile=[hem{h} '.striatum.topo.gii'];
            end
            cd(surfaceGroupDir)
            %----define name of coord and topology
           
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for st=1:2
                for i=1:length(sessN);
                    filename=[surfaceGroupDir filesep hem{h} '.RS_distDeviation_' seqName{st} '_sess' num2str(sessN(i)) '.metric']; % unsmoothed
                    sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations',nIter);
                end;
            end
        end;
    case 'SCALE_surface_sn'
        sessN=4;
        sn = [4:9,11:14,18:23];
        structure = 'cortex'; % cortex or striatum
        vararginoptions(varargin,{'sessN','sn','structure'});
        for ss=sessN
            for s=sn
                for h=1:2 % hemisphere
                    for st=1:2
                        if st==1
                            lab = {'TrainSeq','trained'};
                            colPsc = [1,2];
                            colDist = [3,4];
                        else
                            lab = {'UntrainSeq','untrained'};
                            colPsc = [3,4];
                            colDist = [5,6];
                        end
                        switch structure
                            case 'cortex'
                                surfaceDir = fullfile(caretDir,sprintf('x%s',subj_name{s}),hemName{h});
                                p = caret_load(fullfile(surfaceDir,sprintf('%s_Contrasts_RS_sess%d.metric',subj_name{s},sessN)));
                                d = caret_load(fullfile(surfaceDir,sprintf('%s_sess%d_dist_repsup.metric',subj_name{s},sessN)));
                            case 'striatum'
                                surfaceDir = fullfile(BGDir,'surface',sprintf('x%s',subj_name{s}),hemName{h});
                                p1 = caret_load(fullfile(surfaceDir,sprintf('%s.psc_RS_%s_1st_sess%d.metric',hem{h},lab{1},ss)));
                                p2 = caret_load(fullfile(surfaceDir,sprintf('%s.psc_RS_%s_2nd_sess%d.metric',hem{h},lab{1},ss)));
                                d1 = caret_load(fullfile(surfaceDir,sprintf('%s.dist_RS_%s_exe1_sess%d.metric',hem{h},lab{2},ss)));
                                d2 = caret_load(fullfile(surfaceDir,sprintf('%s.dist_RS_%s_exe2_sess%d.metric',hem{h},lab{2},ss)));
                        end
                        nVert = size(p.data,1);
                        for v=1:nVert
                            % difference in psc (dP) and dist (dD) for exe 1-2 at each vertex
                            scale = p.data(v,colPsc(1))-p.data(v,colPsc(2))/p.data(v,colPsc(1));
                            % calculate expected second distance
                            expDist = scale*d.data(v,colDist(1));
                            % deviation from expected
                            % + if higher than expected, - if lower
                            % average across subjects
                            dev(v,:) = d.data(v,colDist(2))-expDist;
                           % dev(v,:) = (d.data(v,colDist(2))-expDist)/expDist;
                           %ratioObs = d.data(v,colDist(1))/d.data(v,colDist(2));
                           %ratioExp = expDist/d.data(v,colDist(2));
                          % dev(v,:) = ratioObs-ratioExp;
                        end
                        colName = sprintf('distDeviation_fromScaling_%s_sess-%d',lab{2},ss);
                        C = caret_struct('metric','data',dev,'column_name',{colName});
                        C.index=p.index;
                        C.data=dev;
                        caret_save(fullfile(surfaceDir, sprintf('%s.RS_distDeviation_%s_sess%d.metric',hem{h},lab{2},ss)),C);
                    end; % seqType
                end; % hemi
            end; % sn
        end; % sessN
    case 'SCALE_surface_makeGroup'
        sessN       = 4;
        sn          = [4:9,11:14,18:23];
        OUTname     = {'trained','untrained'};
        vararginoptions(varargin,{'sessN','sn','imageType','OUTname','inputcol','replaceNaN'});

        inputcol   = 1:length(OUTname)*length(sessN);
        replaceNaN = ones(size(inputcol));
        
        % Loop over hemispheres.
        for h = 1:2
            % Go to the directory where the group surface atlas resides
            surfaceGroupDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
            cd(surfaceGroupDir);
            % Loop over each input metric file in 'OUTname' and make a group metric file
            indx=1; % re-set for each hemisphere
            for j = 1:length(OUTname);
                for  ss=sessN
                    % Loop over subjects...
                    for s = 1:length(sn);
                        % ...and define the names of their metric files
                        infilenames{indx}{s} = fullfile(caretDir,sprintf('x%s',subj_name{sn(s)}),hemName{h},sprintf('%s.RS_distDeviation_%s_sess%d.metric',hem{h},OUTname{j},ss));                     
                        outcolnames{indx}{s} = sprintf('%s_RS_distDeviation_%s_sess%d',subj_name{sn(s)},OUTname{j},ss);
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    outfilenames{indx} = [surfaceGroupDir filesep hem{h} '.' 'RS_distDeviation_group_' OUTname{j} '_sess' num2str(ss) '.metric'];                 
                    % Finally, make the group metric file for this metric type/contrast
                    caret_metricpermute(infilenames{indx},'outfilenames',outfilenames(indx),'inputcol',inputcol(indx),'replaceNaNs',replaceNaN(indx),'outcolnames',outcolnames{indx});
                    % Verbose display to user
                    fprintf('hem: %i  sess: %i \n',h,ss);
                end;
            end;
        end;
    case 'SCALE_surface_cSPM'
        sessN=4;
        sn=[4:9,11:14,18:23];
        SPMname={'trained','untrained'};
        imageType='RS_distDeviation';
        
        s=1:length(sn);
        for ss=sessN
            SummaryName = sprintf('.summary_%s_sess%d.metric',imageType,ss);
            
            for h=1:2
                surfaceGroupDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
                cd(surfaceGroupDir);
                %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
                for i=1:length(SPMname);                 
                    filenames{i}=[surfaceGroupDir filesep hem{h} '.' imageType '_group_' SPMname{i} '_sess' num2str(ss) '.metric'];     
                end;
                %----loop over the metric files and calculate the cSPM of each with the non-smoothed metrics
                for i=1:length(SPMname);
                    Data=caret_load(filenames{i});
                    cSPM=caret_getcSPM('onesample_t','data',Data.data(:,s),'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                    caret_savecSPM([surfaceGroupDir filesep hem{h} '.' imageType '_' SPMname{i} '_stats.metric'],cSPM);
                    save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');
                    data(:,i)=cSPM.con(1).con; % mean
                    data(:,i+length(SPMname))=cSPM.con(1).Z; % T
                    column_name{i}=['mean_' imageType '_' SPMname{i} '_sess' num2str(ss)];
                    column_name{i+length(SPMname)}=['T_' imageType '_' SPMname{i} '_sess' num2str(ss)];
                end;
                C = caret_struct('metric','data',data,'column_name',column_name);
                caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);
                clear Data C cSPM data;
            end;
            fprintf('Done sess-%d %s\n',ss,imageType);
        end
    case 'SCALE_surface_group_smooth'
        sessN=4;
        structure = 'cortex';
        nIter = 5; % 20 for cortex, 3 for striatum
        vararginoptions(varargin,{'sessN','structure','nIter'});
        seqName = {'trained','untrained'};
        for h=1:2       
            switch structure
                case 'cortex'
                    surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                    coordfile=[hem{h} '.WHITE.coord'];
                    topofile=[hem{h} '.CLOSED.topo'];
                case 'striatum'
                    surfaceGroupDir = fullfile(BGDir,'surface','group',hemName{h});
                    coordfile=[hem{h} '.striatum.coord.gii'];
                    topofile=[hem{h} '.striatum.topo.gii'];
            end
            cd(surfaceGroupDir)
            %----define name of coord and topology
           
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for st=1:2
                for i=1:length(sessN);
                    filename=[surfaceGroupDir filesep hem{h} '.summary_RS_distDeviation_sess' num2str(sessN(i)) '.metric']; % unsmoothed
                    sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations',nIter);
                end;
            end
        end;
        
    case 'PLOT_psc_CCN'
        parcelType = 'cortex_buckner'; % cortex or striatum
        sessN = 4;
        betaChoice = 'multiPW';
        roi=2;
        vararginoptions(varargin,{'betaChoice','roi','parcelType','sessN'});
        
        switch parcelType
            case 'cortex_buckner'
                S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
            case 'striatum_buckner_hemi'
                S = load(fullfile(BGDir,'stats',sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice)));
        end
        S = getrow(S,S.sessN==sessN & S.roi==roi);
        
        TT=[];
        figure
        for st=1:2 % sequence type
            subplot(1,2,st);
            if st==1
               var=S.psc_train;
            else
               var=S.psc_untrain;
            end  
            T.psc = var;
            T.FoSEx=S.FoSEx;
            T.seqType = ones(size(var))*st;
            T.sn = S.sn;
            TT=addstruct(TT,T);
        end
        
        figure
        plt.bar(TT.seqType,TT.psc,'split',TT.FoSEx,'style',stySeqType_exe,'leg',{'exe1','exe2'});
        set(gca,'XTick',[1.5,5],'XTickLabel',{'trained','untrained'});
        ylabel('percent signal change');
    case 'PLOT_pcm_corr_CCN'
        reg = 2;
        sessN=4;
        seqType={'trained'};
        modelType='specific';
        parcelType = {'cortex_networks','striatum_networks'}; % cortex_networks or striatum_networks
       % parcelType = {'striatum_networks'}; % cortex_networks or striatum_networks
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','parcelType'});
        
        T=[];
        for p=1:size(parcelType,2)
            for st=1:size(seqType,2)
                for t=sessN
                    R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_%s_sess%d_%s_%s.mat',parcelType{p},t,seqType{st},modelType)));
                    R=getrow(R,R.roi==reg);
                    R.sessN=ones(size(R.SN))*t;
                    R.parcel=ones(size(R.SN))*p;
                    R.seqType=ones(size(R.SN))*st;
                    T=addstruct(T,R);
                end
            end
        end
        corrType={'pcm','mean-naive','crossval','crossval mean'};
        corrVar=[T.r_model2, T.r_naive_meanOnly, T.r_crossval3, T.r_crossval_meanOnly_evenOdd];
        for ct=1:size(corrVar,2) % three types of correlation 
            figure
            plt.bar(T.parcel,corrVar(:,ct),'subset',T.roi==reg); 
            drawline(0,'dir','horz');
            title(sprintf('Correlation-%s',corrType{ct}));
            set(gca,'XTickLabel',{'cortex','striatum'});
            ylabel('pattern correlation exe1-2');
        end
    
    case 'job'
        % job over Xmas / NY
        sml1_imana_prep('GLM_contrast_RepSup','sn',[4:9,11:31],'sessN',[1:4]);
        sml1_imana_repsup('PSC_create_RepSup','sn',[4:9,11:31],'sessN',[1:4]);
        sml1_imana_repsup('BETA_get','sn',[4:9,11:31],'sessN',[1:4]);
        sml1_imana_repsup('BETA_get','sn',[4:9,11:31],'sessN',[1:4],'parcelType','162tessels');
        sml1_imana_repsup('BETA_get','sn',[4:9,11:31],'sessN',[1:4],'parcelType','BG-striatum');
        sml1_imana_repsup('BETA_combineGroup');
        sml1_imana_repsup('BETA_combineGroup','parcelType','162tessels');
        sml1_imana_repsup('BETA_combineGroup','parcelType','BG-striatum');
    case 'job2'
        sml1_imana_repsup('SEARCH_run_LDC','sn',[6:9,11:25],'sessN',3,'repInd',2);
    case 'job3'
        sml1_imana_dist('SEARCH_dist_runLDC','sn',[13:25],'sessN',3);
    
   
        
        
    otherwise
        disp('there is no such case.')
end;    % switch(what)
end


%% Local functions

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

function M = pcm_defineSequenceModels_fixed_new(Seq,Chunks,sn)
%  No fixed / run component added
% specific sequence + overall component modelled together
    % add the option of finger identity or natural statistics
    % --------------------------------------
       
    % Model1: First finger model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    nSeq=size(Seq,1);
    M{1}(:,:)=zeros(nSeq,5); % first finger
    for u = 1:nSeq    % 12 seq
        firstfing = Seq(:,1);
        M{1}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model2: All fingers
    M{2}(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{2}(u,j) = length(placenumb);
        end
    end
    %---------------------------------------
    % Model3: Finger transitions
    Trans = [nchoosek([1:5],2); fliplr(nchoosek([1:5],2))];
    M{3}(:,:)=zeros(nSeq,size(Trans,1)); % all transitions (no repetition)
    for k = 1:size(Seq,1)    % set of 12 seq
        for p = 1:size(Trans,1) % all possible transitions
            if any(ismember(Seq(k,1:8),Trans(p,1)))
                ind = find(ismember(Seq(k,1:8),Trans(p,1)));    % find matching digits to first finger in doublet
                if any(ismember(Seq(k,ind+1),Trans(p,2)))       % compare the second finger in doublet
                    M{3}(k,p) = M{3}(k,p) + sum(ismember(Seq(k,ind+1),Trans(p,2)));
                end
            end
        end
    end
    tmp=M{3};
    % remove columns with 0s only
    tmp(:,~any(tmp,1))=[];
    M{3}=tmp;
    %----------------------------------------
    % Model4: Chunk - trained model
    nChunks=max(max(Chunks));
    A=zeros(nSeq,nChunks);
    for i=1:nSeq
        idxChunk=Chunks(i,:);
        A(i,idxChunk)=1;
    end
    M{4} = [A(1:6,1:7);zeros(6,7)];
    %----------------------------------------
    % Model5: Chunk - untrained model
    M{5} = [zeros(6,7);A(7:12,8:end)];
    % ---------------------------------------
    % Model6: Model with trained & untrained labels - 
    M{6}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{6}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model7: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{7}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model8: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{8}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
end
function K = pcm_calc(likelihood,compI)
% calculates posterior probability, knockIN/OUT factors
% INPUT: likelihood: (numSubj x numModels) estimated marginal log-likelihood of each model
%        compI: (numModels x numComp) Set of indicator variables showing the presence / absence
%               of the components for each fitted (pcm_constructModelFamily)
% OUTPUT: K: structure with postProp, logBayes, knockIN / OUT factors 

    % calculate posterior probability, logBayes factor
    [postProp,logBayes]=pcm_componentPosterior(likelihood,compI);

    % knockIN, knockOUT
    [knockIN,knockOUT]=pcm_knockModels(likelihood,compI);

    K.postProp=postProp(:);
    K.logBayes=logBayes(:);
    K.knockIN=knockIN(:);
    K.knockOUT=knockOUT(:);
    numComp=size(knockIN,2);
    K.indx=[];
    for i=1:numComp
        K.indx=[K.indx; ones(length(postProp),1)*i];
    end
    K.sn = repmat([1:length(knockIN)]',numComp,1);
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
    C.r_naive_wMean(p,1)    = calcCorr(G);
    C.r_naive_meanOnly(p,1) = corr(mean(b(1:6,:))',mean(b(7:12,:))');
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
 %   Gcv(:,:,p)         = pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    %Gcv2(:,:,p)        = crossval_estG(Data{p},Z,condVec{p});
 %   C.r_crossval(p,1)  = calcCorr(pcm_makePD(Gcv(:,:,p)));
    %C.r_crossval2(p,1) = calcCorr(pcm_makePD(Gcv2(:,:,p)));
    A=[]; meanP=[];
    % new - different way of subtracting mean for each run - 1st / 2nd
    for i=1:numel(unique(partVec{p}))
        Xa = Z(partVec{p}==i,:);
        Ya = Data{p}(partVec{p}==i,:);
        A(:,:,i) = pinv(Xa)*Ya;
        meanP(1,:,i) = mean(A(1:6,:,i));
        meanP(2,:,i) = mean(A(7:12,:,i));
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
   % G2 = cov(b');
    C.r_crossval3(p,1) = calcCorr(pcm_makePD(Gcv3(:,:,p))); % this is the correct one
    % calculate correlation of mean pattern
    % make even and odd mean patterns for 1st / 2nd
    tmp11 = mean(meanP(1,:,[1:2:8]),3);
    tmp12 = mean(meanP(1,:,[2:2:8]),3);
    tmp21 = mean(meanP(2,:,[1:2:8]),3);
    tmp22 = mean(meanP(2,:,[2:2:8]),3);
    cor1 = corr(tmp11',tmp22');
    cor2 = corr(tmp12',tmp21');
    C.r_crossval_meanOnly_evenOdd(p,1) = mean([cor1 cor2]); 
    % calculate correlation of Gs for 1st / 2nd
    RDM = rsa_squareRDM(rsa.distanceLDC(data,partVec{p},condVec{p}));
    Exe{1}=RDM(1:6,1:6);
    Exe{2}=RDM(7:12,7:12);
    C.RDM_corrRSA(p,1)  = rsa_calcCorrRDMs(Exe);
    C.RDM_corrDist(p,1) = rsa_calcDistCorrRDMs(Exe);
    % calculate the covariance structure
    cov_tmp=pcm_estGCrossval(data,partVec{p},condVec{p});
    C.cov(p,:)=rsa_vectorizeIPMfull(cov_tmp);
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
function vox=voxelSelect(Data,partVec,condVec)
sn=size(Data,2);
for p=1:sn
    data=Data{p};
    numVox=size(data,2);
    % indicator matrix for all sequences
    indM=indicatorMatrix('allpairs',[1:6]);
    %indM=indicatorMatrix('allpairs',[1:12]);
    for n=1:numVox
        G1=pcm_estGCrossval(data(:,n),partVec,condVec);
        %G1=pcm_estGCrossval(data(condVec<7,n),partVec(condVec<7),condVec(condVec<7));
        %G2=pcm_estGCrossval(data(condVec>6,n),partVec(condVec>6),condVec(condVec>6));
        dist1=rsa.rdm.squareRDM(diag(indM*G1*indM'));
        dist1=triu(dist1);
        dist1(dist1==0)=NaN;
        distVox1(n)=nanmean(dist1(:));
        %dist2=rsa.rdm.squareRDM(diag(indM*G2*indM'));
        %dist2=triu(dist2);
        %dist2(dist2==0)=NaN;
        %distVox2(n)=nanmean(dist2(:));
    end
    % sort voxels
    [sortedX, sortedInds1]=sort(distVox1,'descend');
    %topVox1 = sortedInds1(1:numVox*0.15);
    topVox1 = sortedInds1(sortedX>0);
    %[sortedX, sortedInds2]=sort(distVox2,'descend');
    %topVox2 = sortedInds2(1:numVox*0.25);
    %topVox=intersect(topVox1,topVox2);
    %vox{p}=data(:,topVox);
    vox{p}=data(:,topVox1);
end
         
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
