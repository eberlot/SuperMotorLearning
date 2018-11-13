function varargout=sml1_imana_dist(what,varargin)

% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
behavDir        =[baseDir '/behavioral_data/data'];            
imagingDir      =[baseDir '/imaging_data'];                     
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
betaDir         =[baseDir '/betas'];
BGDir           =[baseDir '/basal_ganglia_new'];
suitDir         =[baseDir '/suit'];
physioDir       =[baseDir '/physio'];
pcmDir          =[baseDir '/pcm_stats'];
distPscDir      =[baseDir '/dist_psc_stats'];
QCDir           =[baseDir '/quality_control'];

% update glmDir when adding new glms
glmLocDir       ={[baseDir '/glmLoc/glmL1'],[baseDir '/glmLoc/glmL2'],[baseDir '/glmLoc/glmL3']};   % localiser glm
glmLocSessDir   ={[baseDir '/glmLocSeiess/glmLocSess1'],[baseDir '/glmLocSess/glmLocSess2'],[baseDir '/glmLocSess/glmLocSess3'],[baseDir '/glmLocSess/glmLocSess4']}; % one glm for loc run per session
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session


% ------------------------- Experiment Info -------------------------------

% Stimuli - numbers given in SeqNumb
num_train = 1:6;
num_untrain = 7:12;   
num_seq = 1:12;
num_fing = 13:17;
num_seqtype = 2;

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

% per session
numruns_sess      = 10;  
numruns_task_sess = 8;
numruns_loc_sess  = 2;

% total - per subject (the total in the end will always be 40)
numruns           = [40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 40 40 40 30 40];
numruns_task      = 32;
numruns_loc       = 8;
sess = [repmat(1,1,10),repmat(2,1,10),repmat(3,1,10),repmat(4,1,10)];   % all sessions

sess_sn = [4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,3,4];    % per subject

run_task   = [1:3 5:7 9:10;
              11:13 15:17 19:20;
              21:23 25:27 29:30;
              31:33 35:37 39:40];    % task
run_loc    = [4 8;
              14 18;
              24 28;
              34 38];             % localizer

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

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
              's11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
              's21','s22','s23','s24','s25','s26','s27','s28','s29','s30'};  


% -------------------------- For plotting ---------------------------------
gray=[120 120 120]/255;
lightgray=[170 170 170]/255;
silver=[230 230 230]/255;
black=[0 0 0]/255;
blue=[49,130,189]/255;
mediumblue=[128,231,207]/255;
lightblue=[158,202,225]/255;
red=[222,45,38]/255;
mediumred=[237,95,76]/255;
lightred=[252,146,114]/255;
ms=12;
stySeq=style.custom({'red','blue'},'markersize',ms);
stySeqType=style.custom(gray,'markersize',12);

stySess=style.custom({gray,lightgray,silver,black},'markersize',ms);
styTrained_sess=style.custom({red,mediumred,lightred},'markersize',ms);
styUntrained_sess=style.custom({blue,mediumblue,lightblue},'markersize',ms);

% ------------------------------ Analysis Cases --------------------------------
switch(what)
     
    case 'PSC_create'
    % calculate psc for trained and untrained sequences - based on betas    
    vararginoptions(varargin,{'sn','sessN'});
    name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','TrainSeq','UntrainSeq'};
    for ss=sessN
        for s=sn
           % cd(fullfile(glmSessDir{ss}, subj_name{s}));
           cd(fullfile(glmSessDir{ss},subj_name{s}));
            load SPM;
            T=load('SPM_info.mat');
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
            fprintf('Subject %d sess %d: %3.3f\n',s,ss,h);
        end;
    end
    case 'PSC_surface'
     % create surface maps of percent signal change 
     % trained and untrained sequences
     
        smooth = 0;   
        vararginoptions(varargin,{'sn','sessN','smooth'});

        hemisphere=1:length(hem);
        fileList = [];
        column_name = [];
        name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','TrainSeq','UntrainSeq'};
        
        for ss=sessN
            for n = 1:length(name)
                fileList{n}=fullfile(['psc_sess' num2str(ss) '_' name{n} '.nii']);
                column_name{n} = fullfile(sprintf('Sess%d_%s.nii',ss,name{n}));
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
                        images{f}=fullfile(glmSessDir{ss},subj_name{s},fileList{f});
                    end;
                    metric_out = fullfile(caretSDir,sprintf('%s_Contrasts_sess%d.metric',subj_name{s},ss));
                    M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                    M.column_name = column_name;
                    caret_save(metric_out,M);
                    fprintf('Sess %d, Subj %d, Hem %d\n',ss,s,h);
                    
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        Out = caret_smooth(metric_out, 'coord', white, 'topo', closed);%,...
                        %'algorithm','FWHM','fwhm',12);
                        char(Out);  % if smoothed adds an 's'
                    else
                    end;
                    
                end; % hemi
            end; % sn
        end; % session
    case 'PSC_create_loc'
        % calculate psc for all digits / ind digit vs. rest - based on betas    
        glm=2;
        vararginoptions(varargin,{'sn','glm'});
        name={'DigitAny','Digit1','Digit2','Digit3','Digit4','Digit5'};
        for s=sn
            cd(fullfile(glmLocDir{glm}, subj_name{s}));
            load SPM;
            T=load('SPM_info.mat');
            X=(SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
            h=median(max(X));               % Height of response;
            P={};
            numB=length(SPM.xX.iB);         % Partitions - runs
            for p=SPM.xX.iB
                P{end+1}=sprintf('beta_%4.4d.nii',p);       % get the intercepts and use them to calculate the baseline (mean images)
            end;
            for con=1:length(name)    % 6 contrasts
                P{numB+1}=sprintf('con_%s.nii',name{con});
                outname=sprintf('psc_digit_%s.nii',name{con}); % ,subj_name{s}
                
                formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);    % 8 runs overall
                  
                spm_imcalc_ui(P,outname,formula,...
                    {0,[],spm_type(16),[]});        % Calculate percent signal change
            end;
            fprintf('Subject %d: %3.3f\n',s,h);
        end;
    case 'PSC_surface_loc'
        % create surface maps of percent signal change 
     % trained and untrained sequences
     
        smooth=0;  
        glm=2;
        vararginoptions(varargin,{'sn','glm','smooth'});

        hemisphere=1:length(hem);
        fileList = [];
        column_name = [];
        name={'DigitAny','Digit1','Digit2','Digit3','Digit4','Digit5'};
        for n = 1:length(name)
            fileList{n}=fullfile(['psc_digit_' name{n} '.nii']);
            column_name{n} = fullfile(sprintf('%s.nii',name{n}));
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
                    images{f}=fullfile(glmLocDir{glm},subj_name{s},fileList{f});
                end;
                metric_out = fullfile(caretSDir,sprintf('%s_Digit_PSC.metric',subj_name{s}));
                M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                M.column_name = column_name;
                caret_save(metric_out,M);
                fprintf('Subj %d, Hem %d\n',s,h);
                
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
    case 'PSC_group_make'
        % Calculate group metric files from contrast / psc calculations. 
        % Takes 2 contrast results ('TrainSeq','UntrainSeq') across
        % subjects and makes a group level metric file that contains each
        % subject's data for that contrast type.
        sessN=1;
        sn=[1:18];
        vararginoptions(varargin,{'sessN','sn'});
        % Some presets
        name    = 'Contrasts';

        OUTname    = {'TrainSeq','UntrainSeq'};
        inputcol   = [13 14]; % only averages - all trained / untrained
        replaceNaN = [1 1];
        
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
                    infilenames{j}{i} = fullfile(caretDir,[atlasA subj_name{sn(i)}], hemName{h}, sprintf('%s_%s_sess%d.metric',subj_name{sn(i)},name,sessN));
                    % Name the output filename for this group metric file in average surface folder
                end;
                outfilenames{j} = [surfaceGroupDir filesep hem{h} '.' OUTname{j} '_sess' num2str(sessN) '.metric'];
                % Finally, make the group metric file for this metric type/contrast
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
                % Verbose display to user
                fprintf('hem: %i  image: %i \n', h,j);
            end;
        end;    
    case 'PSC_group_cSPM'
        % Calculate group stats files from the group metric files. 
        % Takes 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross')  and 
        % calculates group level stats (one sample t test that the mean 
        % effect is bigger than zero). 
        % 
        % Although we calculate the t-score, corresponding p-value for that
        % t-score, and finally the z-score 
        % subject's data for that contrast type.
        
        sessN=1;
        sn=[1:18];
        s=1:length(sn);
        vararginoptions(varargin,{'sessN','sn'});
        
        SPMname={'TrainSeq','UntrainSeq'};
        
        sqrtTransform=[0,0,0,0]; % Should you take ssqrt before submitting? 
                                % no for psc
        SummaryName = sprintf('.summary_psc_sess%d.metric',sessN);
        hemi = [1 2];
        
        for h=hemi
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
            %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
            for i=1:length(SPMname);
                filenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '_sess' num2str(sessN) '.metric']; % no smoothing
                sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '_sess' num2str(sessN) '.metric']; % smoothing
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
                column_name{i}=['mean_' SPMname{i} '_sess' num2str(sessN)];
                column_name{i+length(SPMname)}=['T_' SPMname{i} '_sess' num2str(sessN)];
            end;
            C = caret_struct('metric','data',data,'column_name',column_name);
            caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);

        end;
        fprintf('Done \n')
    case 'PSC_group_smooth'
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
                filename=[surfaceGroupDir filesep hem{h} '.summary_psc_sess' num2str(sessN) '.metric']; % unsmoothed
                sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
            end;

        end;
    case 'PSC_crossflat'
        sessN=1;
        vararginoptions(varargin,{'sessN'});
        
        LimDist = [-0.5 3.5];
        coord_start = [-33 8]; 
        coord_end = [40 10];
        cross_width = 20;
        
        for h=1
            surfGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            flatFile = fullfile(surfGroupDir,[hem{h} '.FLAT.coord']);
            flatcoord=caret_load(flatFile);            
          
            Y_trained=fullfile(surfGroupDir,sprintf('%s.TrainSeq_sess%d.metric',hem{h},sessN));
            Y_untrained=fullfile(surfGroupDir,sprintf('%s.UntrainSeq_sess%d.metric',hem{h},sessN));

            S=caret_crosssection_flat(flatcoord.data,fullfile(surfGroupDir,[hem{h} '.surface_shape']),coord_start(h,:),coord_end(h,:),cross_width);
            [Y.tr P coord]=caret_crosssection_flat(flatcoord.data,Y_trained,coord_start(h,:),coord_end(h,:),cross_width);
            Y.utr=caret_crosssection_flat(flatcoord.data,Y_untrained,coord_start(h,:),coord_end(h,:),cross_width);
            
            x=[1:size(Y.tr,1)];
            figure;
            title(sprintf('session %d',sessN))
            subplot(3,1,[1:2]);
            hold on;
            traceplot(x,Y.tr','errorfcn','stderr','linecolor','r','linewidth',2,'patchcolor','r');
            traceplot(x,Y.utr','errorfcn','stderr','linecolor','b','linewidth',2,'patchcolor','b');
            ylabel('Percent signal change');
            ylim(LimDist);
            set(gca,'XLim',[1 length(x)]);
            drawline(0,'dir','horz');
            drawline(1,'dir','horz');
            drawline(2,'dir','horz');
            drawline(3,'dir','horz');
            
            C=caret_struct('metric','data',P); % logical value - vertex included
            
            %caret_save(fullfile(surfaceGroupDir,'summary.RGB_paint'),C);
            
            caret_save([surfGroupDir filesep hem{h} sprintf('.crosssection_vertex.metric')],C);
        end;
        subplot(3,1,3);
        plot(x,S(:,2));
        set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.5 1.5]);
        set(gcf,'PaperPosition',[2 2 5 3]);

    case 'SEARCH_all'
        vararginoptions(varargin,{'sn','sessN'});
        if sessN == 1
        sml1_imana('SEARCH_define','sn',sn,'sessN',sessN);
        end
        sml1_imana('SEARCH_dist_runLDC','sn',sn,'sessN',sessN);
        sml1_imana('SEARCH_dist_map','sn',sn,'sessN',sessN);
        sml1_imana('SEARCH_dist_surf','sn',sn,'sessN',sessN);
    case 'SEARCH_define'                                                    % STEP 4.1a   :  Defines searchlights for 120 voxels in grey matter surface
        vararginoptions(varargin,{'sn','sessN'});
        
        for s=sn
            mask       = fullfile(regDir,sprintf('mask_%s.nii',subj_name{s}));
            Vmask      = spm_vol(mask);
            Vmask.data = spm_read_vols(Vmask);
            
            LcaretDir = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
            RcaretDir = fullfile(caretDir,['x',subj_name{s}],'RightHem');
            white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
            pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
            topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
            S         = rsa_readSurf(white,pial,topo);
            
            L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[15 120]);
            save(fullfile(anatomicalDir,subj_name{s},sprintf('%s_searchlight_120.mat',subj_name{s})),'-struct','L');
        end
    case 'SEARCH_dist_runLDC'                                               % STEP 4.2a   :  Runs LDC searchlight using defined searchlights (above)
        sn = [14:25];
        sessN = [1:4];
        vararginoptions(varargin,{'sn','sessN'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        runs = 1:numruns_task_sess;
        % make index vectors
        conditionVec  = kron(ones(numel(runs),1),[1:12]');      % 12 sequences
        partition     = kron(runs',ones(12,1));
        for s=sn
            for ss=sessN
                % go to subject's glm directory
                cd(fullfile(glmSessDir{sessN},subj_name{s}));
                if exist(sprintf('%s_sess%d_LDC.nii',subj_name{s},ss));
                    fprintf('Searchlight map: %s sess%d already exists - skipping\n',subj_name{s},ss);
                else
                    % load their searchlight definitions and SPM file
                    L = load(fullfile(anatomicalDir,subj_name{s},sprintf('%s_searchlight_120.mat',subj_name{s})));
                    load SPM;
                    SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                    
                    name = sprintf('%s_sess%d',subj_name{s},ss);
                    % run the searchlight
                    rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name,'idealBlock',block);
                    fprintf('Searchlight map: %s sess%d done\n\n\n',subj_name{s},ss);
                end
            end
        end
        cd(cwd);
    case 'SEARCH_dist_map'                                                  % STEP 4.3a   :  Averaged LDC values for specified contrasts
        sn  = 1;
        sessN = 1;
        vararginoptions(varargin,{'sn','sessN'});

        cWD = cd;
        for s = sn
            % Load subject surface searchlight results (1 vol per paired conds)
            LDC_file            = fullfile(glmSessDir{sessN},subj_name{s},sprintf('%s_sess%d_LDC.nii',subj_name{s},sessN)); % searchlight nifti
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
            Y.fname   = sprintf('%s_sess%d_dist.nii',subj_name{s},sessN);
            
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
            T.fname = sprintf('%s_sess%d_dist_trained.nii',subj_name{s},sessN);
            U.fname = sprintf('%s_sess%d_dist_untrained.nii',subj_name{s},sessN);
            Z.fname = sprintf('%s_sess%d_dist_cross.nii',subj_name{s},sessN);
            
            % save outputs
            spm_write_vol(Y,Y.LDC);
            spm_write_vol(T,T.LDC);
            spm_write_vol(U,U.LDC);
            spm_write_vol(Z,Z.LDC);
            fprintf('Done %s_sess%d \n',subj_name{s},sessN)
            
            clear vol vdat LDC Y T U Z
            
        end
        cd(cWD);  % return to working directory
    case 'SEARCH_dist_surf'                                                 % STEP 4.4a   :  Map searchlight results (.nii) onto surface (.metric)
        % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        sessN = 1;   
        fileList = {'dist','dist_trained','dist_untrained','dist_cross'}; % for dist or repsup
        glmDir = glmSessDir;
        outname = 'dist'; % dist or dist_repsup
        vararginoptions(varargin,{'sn','sessN','fileList','glmDir','outname'});

        for s = sn
            for h=1:2
                caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
                pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
                
                for f = 1:length(fileList)
                    images{f}    = fullfile(glmDir{sessN},subj_name{s},sprintf('%s_sess%d_%s.nii',subj_name{s},sessN,fileList{f}));
                    column_name{f} = fullfile(sprintf('Sess%d_%s.nii',sessN,fileList{f}));
                end;    % filename
                outfile   = sprintf('%s_sess%d_%s.metric',subj_name{s},sessN,outname);
                M         = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
                M.column_name = column_name;
                caret_save(fullfile(caretSDir,outfile),M);
                fprintf('Done subj %d sess %d hemi %d \n',s,sessN,h);
            end;    % hemi
        end;    % subj
    case 'SEARCH_group_make'                                                % STEP 4.5   :  Make group metric files by condensing subjec contrast metric files
        % Calculate group metric files from the searchlight results. 
        % Takes the 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross') across
        % subjects and makes a group level metric file that contains each
        % subject's data for that contrast type.
        sessN=1;
        sn=[1:7];
        name = 'dist';
        inputcol   = [1 2 3 4];
        replaceNaN = [1 1 1 1];
        OUTname    = {'dist_all','dist_trained','dist_untrained','dist_cross'};
        vararginoptions(varargin,{'sessN','sn','name','OUTname','inputcol','replaceNaN'});


        % Loop over hemispheres.
        for h = 1:2
            % Go to the directory where the group surface atlas resides
            surfaceGroupDir = [caretDir filesep atlasname filesep hemName{h}];
            cd(surfaceGroupDir);
            % Loop over each input metric file in 'OUTname' and make a group metric file
            for j = 1:length(OUTname);
                % Loop over subjects...
                for i = 1:length(sn);
                    % ...and define the names of their metric files
                    infilenames{j}{i} = fullfile(caretDir,[atlasA subj_name{sn(i)}], hemName{h}, sprintf('%s_sess%d_%s.metric',subj_name{sn(i)},sessN,name));
                    % Name the output filename for this group metric file in average surface folder
                end;
                outfilenames{j} = [surfaceGroupDir filesep hem{h} '.' OUTname{j} '_sess' num2str(sessN) '.metric'];
                % Finally, make the group metric file for this metric type/contrast
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
                % Verbose display to user
                fprintf('hem: %i  image: %i \n', h,j);
            end;
        end;
    case 'SEARCH_group_cSPM'                                                % STEP 4.6   :  Generate a statistical surface map (onesample_t test) from smoothed group metric files. Also avgs. distances across subjs.
        % Calculate group stats files from the group metric files. 
        % Takes 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross')  and 
        % calculates group level stats (one sample t test that the mean 
        % effect is bigger than zero). 
        % 
        % Although we calculate the t-score, corresponding p-value for that
        % t-score, and finally the z-score 
        % subject's data for that contrast type.    
        sessN=1;
        sn=[1:7];
        SPMname={'dist_all','dist_trained','dist_untrained','dist_cross'};
        name='dist';
        sqrtTransform=[1,1,1,1]; % Should you take ssqrt before submitting?
        % Yes, b/c used rsa.distanceLDC to
        % calculate distances. This function
        % returns squared cv mahalanobis distance.
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name'});
        s=1:length(sn);
        SummaryName = sprintf('.summary_%s_sess%d.metric',name,sessN);
        hemi = [1 2];
        
        for h=hemi
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
            %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
            for i=1:length(SPMname);
                %sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '.metric']; % smoothed
                sfilenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '_sess' num2str(sessN) '.metric']; % no smoothing
            end;
            %----loop over the metric files and calculate the cSPM of each with the non-smoothed metrics
            for i=1:length(SPMname);
                Data=caret_load(sfilenames{i});
                if sqrtTransform(i)
                    Data.data=ssqrt(Data.data);
                end;
                cSPM=caret_getcSPM('onesample_t','data',Data.data(:,s),'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                caret_savecSPM([surfaceGroupDir filesep hem{h} '.' SPMname{i} '_stats.metric'],cSPM);
                save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');
                data(:,i)=cSPM.con(1).con; % mean
                data(:,i+length(SPMname))=cSPM.con(1).Z; % T
                column_name{i}=['mean_' SPMname{i} '_sess' num2str(sessN)];
                column_name{i+length(SPMname)}=['T_' SPMname{i} '_sess' num2str(sessN)];
            end;
            C = caret_struct('metric','data',data,'column_name',column_name);
            caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);
        end;
        fprintf('Done \n')
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
                filename=[surfaceGroupDir filesep hem{h} '.summary_dist_sess' num2str(sessN) '.metric']; % unsmoothed
                sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
            end;

        end;
    case 'SEARCH_dist_crosssection'
        sessN=1;
        vararginoptions(varargin,{'sessN'});
        
        LimDist = [-0.001,0.002];
        
        for h=1
            surfGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            border=fullfile(surfGroupDir,'CS.borderproj');
            
            Y_trained=fullfile(surfGroupDir,sprintf('%s.dist_trained_sess%d.metric',hem{h},sessN));
            Y_untrained=fullfile(surfGroupDir,sprintf('%s.dist_untrained_sess%d.metric',hem{h},sessN));
            
            X=caret_crosssection(border,fullfile(surfGroupDir,[hem{h} '.INFLATED.coord']));
            S=caret_crosssection(border,fullfile(surfGroupDir,[hem{h} '.surface_shape']));
            Y.tr=caret_crosssection(border,Y_trained);
            Y.utr=caret_crosssection(border,Y_untrained);

            x=[1:size(Y.tr,1)];
            figure(sessN);
            subplot(3,1,[1:2]);
            hold on;
            title(sprintf('session %d',sessN))
            traceplot(x,Y.tr','errorfcn','stderr','linecolor','r','linewidth',2,'patchcolor','r');
            traceplot(x,Y.utr','errorfcn','stderr','linecolor','b','linewidth',2,'patchcolor','b');
            ylabel('Average distance');
            ylim(LimDist);
            set(gca,'XLim',[1 length(x)]);
            drawline(0,'dir','horz');
            legend({'','Trained','','Untrained'});
        end;
        subplot(3,1,3);
        plot(x,S(:,2));
        set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.5 1.5]);
        set(gcf,'PaperPosition',[2 2 5 3]);
    case 'SEARCH_dist_crossflat'
        sessN=1;
        vararginoptions(varargin,{'sessN'});
        
        LimDist = [-0.0005,0.0035];
        coord_start = [-33 8]; 
        coord_end = [40 10];
        cross_width = 20;
        
        for h=1
            surfGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            flatFile = fullfile(surfGroupDir,[hem{h} '.FLAT.coord']);
            flatcoord=caret_load(flatFile);            
          
            Y_trained=fullfile(surfGroupDir,sprintf('%s.dist_trained_sess%d.metric',hem{h},sessN));
            Y_untrained=fullfile(surfGroupDir,sprintf('%s.dist_untrained_sess%d.metric',hem{h},sessN));
            Y_cross=fullfile(surfGroupDir,sprintf('%s.dist_cross_sess%d.metric',hem{h},sessN));
            
            S=caret_crosssection_flat(flatcoord.data,fullfile(surfGroupDir,[hem{h} '.surface_shape']),coord_start(h,:),coord_end(h,:),cross_width);
            Y.tr=caret_crosssection_flat(flatcoord.data,Y_trained,coord_start(h,:),coord_end(h,:),cross_width);
            Y.utr=caret_crosssection_flat(flatcoord.data,Y_untrained,coord_start(h,:),coord_end(h,:),cross_width);
            Y.cross=caret_crosssection_flat(flatcoord.data,Y_cross,coord_start(h,:),coord_end(h,:),cross_width);
            
            x=[1:size(Y.tr,1)];
            keyboard;
            % trained / untrained
            figure;
            subplot(3,1,[1:2]);
            hold on;
            title(sprintf('session %d',sessN))
            traceplot(x,Y.tr','errorfcn','stderr','linecolor','r','linewidth',2,'patchcolor','r');
            traceplot(x,Y.utr','errorfcn','stderr','linecolor','b','linewidth',2,'patchcolor','b');
            ylabel('Average distance');
            ylim(LimDist);
            set(gca,'XLim',[1 length(x)]);
            drawline(0,'dir','horz');
        end;
        subplot(3,1,3);
        plot(x,S(:,2));
        set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.5 1.5]);
        set(gcf,'PaperPosition',[2 2 5 3]);
       
        
        % distances trained, untrained, across seq sets
        figure;
        subplot(3,1,[1:2]);
        hold on;
        title(sprintf('session %d',sessN))
        traceplot(x,Y.tr','errorfcn','stderr','linecolor','r','linewidth',2,'patchcolor','r');
        traceplot(x,Y.utr','errorfcn','stderr','linecolor','b','linewidth',2,'patchcolor','b');
        traceplot(x,Y.cross','errorfcn','stderr','linecolor','g','linewidth',2,'patchcolor','g');
        ylabel('Average distance');
        ylim(LimDist);
        set(gca,'XLim',[1 length(x)]);
        drawline(0,'dir','horz');
        subplot(3,1,3);
        plot(x,S(:,2));
        set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.5 1.5]);
        set(gcf,'PaperPosition',[2 2 5 3]);
   
    case 'SEARCH_dist_corr'
        % run correlation distance as searchlight
        sn = [4:9,11:25];
        sessN = [4];
        vararginoptions(varargin,{'sn','sessN'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        runs = 1:numruns_task_sess;
        nConds = length(num_seq);
        % make index vectors
        conditionVec  = kron(ones(numel(runs),1),[1:nConds]');      % 12 sequences
        partitionVec  = kron(runs',ones(nConds,1));
        for s=sn
            for ss=sessN
                % go to subject's glm directory
                cd(fullfile(glmSessDir{sessN},subj_name{s}));
                if exist(sprintf('%s_sess%d_corrDist.nii',subj_name{s},ss));
                    fprintf('Searchlight map: %s sess%d already exists - skipping\n',subj_name{s},ss);
                else
                    % load their searchlight definitions and SPM file
                    L = load(fullfile(anatomicalDir,subj_name{s},sprintf('%s_searchlight_120.mat',subj_name{s})));
                    load SPM;
                    SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                    
                    inFiles = SPM.xY.VY;
                    % define the outfiles
                    nfiles  = nConds * (nConds-1)/2;
                    name = sprintf('%s_sess%d_corrDist',subj_name{s},ss);
                    for i = 1:nfiles
                        outFiles{i} = fullfile(cd,sprintf('%s.nii,%d',name,i));
                    end;
                    % run the searchlight
                    rsa.runSearchlight(L,inFiles,outFiles,@calcCorrSearch,'optionalParams',{SPM,conditionVec,partitionVec},'idealBlock',block);
                    fprintf('Searchlight map corrDist: %s sess%d done\n\n\n',subj_name{s},ss);
                end
            end
        end
        cd(cwd);
        
    case 'SEARCH_fingmap_runLDC'                                            % STEP 4.1b   :  Runs LDC searchlight using defined searchlights (above)

        glm=2;
        vararginoptions(varargin,{'sn','glm'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        for s=sn
            runs = 1:numruns_loc;   % all 8 runs across 4 sessions
            % make index vectors           
            conditionVec  = kron(ones(numel(runs),1),[1:5]');      % 12 sequences
            partition     = kron(runs',ones(5,1));
            % go to subject's glm directory 
            cd(fullfile(glmLocDir{glm},subj_name{s}));
            % load their searchlight definitions and SPM file
            L = load(fullfile(anatomicalDir,subj_name{s},sprintf('%s_searchlight_120.mat',subj_name{s})));
            load SPM;
            SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));

            name = sprintf('%s_fingmap',subj_name{s});
            % run the searchlight
            rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name,'idealBlock',block);

        end
        cd(cwd);
    case 'SEARCH_fingmap_map'
        sn=1;
        glm=2;
        vararginoptions(varargin,{'sn','glm'});

        cWD = cd;
        for s = sn
            % Load subject surface searchlight results (1 vol per paired conds)
            LDC_file            = fullfile(glmLocDir{glm},subj_name{s},sprintf('%s_fingmap_LDC.nii',subj_name{s})); % searchlight nifti
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
            Y.fname   = sprintf('%s_fingdist_mean.nii',subj_name{s});
            
            % per each finger
            C = rsa.util.pairMatrix(5); % ?!?!??!

            % save outputs
            spm_write_vol(Y,Y.LDC);
            fprintf('Done %s_fingmap \n',subj_name{s})
            
            clear vol vdat LDC Y
            
        end
        cd(cWD);  % return to working directory
    case 'SEARCH_fingmap_surf'
         % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        glm = 2;
        
        vararginoptions(varargin,{'sn','glm'});
        fileList = {'fingdist_mean'};
        hemisphere = 1:2;
        
        for s = sn
            for h=hemisphere
                caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
                pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
                
                for f = 1:length(fileList)
                    images{f}    = fullfile(glmLocDir{glm},subj_name{s},sprintf('%s_%s.nii',subj_name{s},fileList{f}));
                    column_name{f} = fullfile(sprintf('%s.nii',fileList{f}));
                end;    % filename
                outfile   = sprintf('%s_fingdist.metric',subj_name{s});
                M         = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
                M.column_name = column_name;
                caret_save(fullfile(caretSDir,outfile),M);
                fprintf('Done subj %d hemi %d \n',s,h);
            end;    % hemi
        end;    % subj
    
    case 'BETA_get'                                                     % STEP 5.6   :  Harvest betas from rois (raw, univ, multiv prewhit)    
        sessN = 1;
        sn  = [23:25];    
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
                load(fullfile(glmSessDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)]));  % load subject's region parcellation (R)
                
                if strcmp(roiDefine,'all')==1
                    roi=1:size(R,2);
                end
                cd (fullfile(glmSessDir{ss},subj_name{s})); % maybe change back when remove glm
                
                P=SPM.Vbeta(SPM.xX.iC);
                
                % Add a few extra images
                %----task against rest
                O{1}=sprintf('psc_sess%d_TrainSeq.nii',ss); %psc trained
                O{2}=sprintf('psc_sess%d_UntrainSeq.nii',ss); %psc untrained
                %O{3}=sprintf('%s_sess%d_dist_trained.nii',subj_name{s},sessN); %dist trained
                %O{4}=sprintf('%s_sess%d_dist_untrained.nii',subj_name{s},sessN); %dist untrained
                %O{5}=sprintf('%s_sess%d_dist_cross.nii',subj_name{s},sessN); %dist cross
                oP=spm_vol(char(O));
                
                V = SPM.xY.VY;
                
                for r = roi % for each region
                    % get raw data for voxels in region
                    % determine if any voxels for that parcel
                    if size(R{r}.data,1)==0 % no voxels
                        % make data into NaN
                        S.betaW                   = NaN;
                        S.betaUW                  = NaN;
                        S.betaRAW                 = NaN;
                        S.resMS                   = NaN;
                        S.psc_train   = NaN;
                        S.psc_untrain = NaN;
                    else
                        Y = region_getdata(V,R{r});  % Data Y is N x P
                        data = region_getdata(oP,R{r}); % from added images
                        
                        % estimate region betas
                        [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                        S.betaW                   = {betaW};                             % multivariate pw
                        S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                        S.betaRAW                 = {beta};
                        S.resMS                   = {resMS};
                        % info from maps for surface
                        S.psc_train   = {data(1,:)};
                        S.psc_untrain = {data(2,:)};
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
-                save(fullfile(betaDir,subj_name{s},sprintf('betas_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'-struct','T');
                fprintf('\nDone beta extraction for sess%d-%s\n',ss,subj_name{s});
            end            
        end
    case 'BETA_combineGroup'
        % combine individual subject beta structures into the whole
        % structure
        sessN=[1:4];
        sn=[4:9,11:28];
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
                S=load(fullfile(betaDir,subj_name{s},sprintf('betas_%s_%s_sess%d',parcelType,subj_name{s},ss)));
                S.regName(1,:) = {sprintf('%s_Putamen_lh',subj_name{s})};
                S.regName(2,:) = {sprintf('%s_CaudateN_lh',subj_name{s})};
                S.regName(3,:) = {sprintf('%s_Putamen_rh',subj_name{s})};
                S.regName(4,:) = {sprintf('%s_CaudateN_rh',subj_name{s})};
                T=addstruct(T,S);
                fprintf('%d.',s);
            end
            save(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)),'-struct','T');
        end
    case 'BETA_splitHalf'
          % crossvalidated version - separately for even and odd runs
        sn=[4:9,11:28,30];
        sessN=[1:4];
        parcelType='162tessels'; % Brodmann or 162tessels
        betaChoice='multi';
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        
        TT=[];
        nRun  = 8;
        nCond = 12;
        pVec  = kron([1:nRun]',ones(nCond,1));
        idxA = mod(pVec,2)==1; % split even and odd runs
        idxB = mod(pVec,2)==0;
        cVec  = kron(ones(nRun/2,1),[1:nCond]');   % for each partition
        
        if strcmp(parcelType,'Brodmann')
            roi     = [1:16];
            regSide = [ones(1,8) ones(1,8)*2];
            regType = [1:8 1:8];
        elseif strcmp(parcelType,'BG-striatum')
            roi     = [1:4];
            regSide = [1 1 2 2];
            regType = [1 2 1 2];
        elseif strcmp(parcelType,'thalamus')
            roi     = [1,2];
            regSide = [1 2];
            regType = [1 1];
        elseif strcmp(parcelType,'162tessels')
            roi     = sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
            % choose only clusters where group dist >0
            regSide = ones(size(roi));
            regSide(roi>162)=2;
            regType = roi;
            regType(regSide==2)=regType(regSide==2)-162;
        end

        for ss=sessN
            fprintf('\n\nSession %d\n',ss);
            V=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss)));
            for s=sn
                fprintf('\nSubject: %d\n',s) % output to user
                for r=1:length(roi)
                    beta=[];
                    D=getrow(V,V.region==roi(r) & V.SN==s);
                    switch betaChoice
                        case 'multi'
                            beta = D.betaW{:};
                        case 'uni'
                            beta = D.betaUW{:};
                        case 'raw'
                            beta = D.betaRAW{:};
                    end
                    if isnan(beta)
                        T.partA_train     = ones(1,6*5/2)*NaN;
                        T.partA_untrain   = ones(1,6*5/2)*NaN;
                        T.partA_all       = ones(1,12*11/2)*NaN;
                        T.partB_train     = ones(1,6*5/2)*NaN;
                        T.partB_untrain   = ones(1,6*5/2)*NaN;
                        T.partB_all       = ones(1,12*11/2)*NaN;
                    else
                        distA = rsa_squareRDM(rsa.distanceLDC(beta(idxA,:),pVec(idxA),cVec));
                        distB = rsa_squareRDM(rsa.distanceLDC(beta(idxB,:),pVec(idxB),cVec));
                        T.partA_train     = rsa_vectorizeRDM(distA(1:6,1:6));
                        T.partA_untrain   = rsa_vectorizeRDM(distA(7:12,7:12));
                        T.partA_all       = rsa_vectorizeRDM(distA);
                        T.partB_train     = rsa_vectorizeRDM(distB(1:6,1:6));
                        T.partB_untrain   = rsa_vectorizeRDM(distB(7:12,7:12));
                        T.partB_all       = rsa_vectorizeRDM(distB);
                    end   
                    T.sn = s;
                    T.sessN = ss;
                    T.roi = roi(r);
                    T.regType = regType(r);
                    T.regSide = regSide(r);
                    fprintf('%d.',r);
                    TT=addstruct(TT,T);
                    fprintf('Done %d/%d\n',find(s==sn),numel(sn));
                end
            end
        end
        save(fullfile(betaDir,'group',sprintf('betas_partition_%s',parcelType)),'-struct','TT');
    case 'BETA_combine_splitHalf'
        % combine different split-half beta structures
        % Brodmann, BG-striatum, thalamus
        sessN=[1:4];
        parcelType = {'Brodmann','BG-striatum','thalamus'};
        
        TT=[];
        for ss=sessN
            idx = [0 0];
            for i=1:length(parcelType)
                T = load(fullfile(betaDir,'group',sprintf('betas_partition_%s.mat',parcelType{i})));
                T = getrow(T,T.sessN==ss);
                T.regStruct = ones(size(T.sn))*i;
                T.roi       = T.roi + idx(1);
                T.regType   = T.regType + idx(2);
                idx = [(max(T.roi)) (max(T.regType))];
                
                TT = addstruct(TT,T);
            end
        end
    save(fullfile(betaDir,'group','betas_partition_combined'),'-struct','TT');
    case 'BETA_get_LOC'
        % for localizer  
        sessN = 1;
        sn  = 1;    
        roi = [1:16];
        vararginoptions(varargin,{'sn','sessN','roi'});
        
        T=[];
            
        % harvest
        for s=sn % for each subj
            fprintf('\nSubject: %d  sess: %d\n',s,sessN) % output to user
            
            % load files
            load(fullfile(glmLocSessDir{sessN}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));          % load subject's region parcellation (R)

            V = SPM.xY.VY; 
            
            for r = roi % for each region
                % get raw data for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P 
                % estimate region betas
                [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                S.betaW                   = {betaW};                             % multivariate pw
                S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw 
                S.betaRAW                 = {beta};
                S.resMS                   = {resMS};
                S.SN                      = s;
                S.region                  = r;
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
        end
        % save T
        save(fullfile(regDir,sprintf('betas_LOC_sess%d.mat',sessN)),'-struct','T'); 
        fprintf('\n');  
    case 'BETA_stats'                                                        % STEP 5.8   :  Calculate stats/distances on activity patterns - train/untrain seq
        sessN = 4;
        sn=[4:9,11:28];
        roi = [1:16];
        roiDefine = 'all'; % determine from region file
        betaChoice = 'multi'; % uni, multi or raw
        checkG=0; % test across all regions / subjects 
        simulations=0;
        type = 'add'; % new or add - if creating from scratch (no subject or adding new ones only)
        parcelType = '162tessels'; % or Brodmann
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','checkG','simulations','type','parcelType','roiDefine'});
        
        for ss=sessN;
            T = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            if strcmp(roiDefine,'all')==1
                roi=unique(T.region)';
            end
            
            % output structures
            Po = [];
            switch(type)
                case 'new'
                    To=[];
                case 'add'
                    To=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
            end
            g_ind=1;
            % do stats
            for s = sn % for each subject
                D = load(fullfile(glmSessDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                fprintf('\nSubject: %d session: %d\n',s,ss)
                num_run = numruns_task_sess;
                
                for r = roi % for each region
                    S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                    fprintf('%d.',r)
                    
                    switch (betaChoice)
                        case 'uni'
                            betaW  = S.betaUW{1};
                        case 'multi'
                            betaW  = S.betaW{1};
                        case 'raw'
                            betaW  = S.betaRAW{1};
                    end
                    
                    % if no data - label distances, RDM as NaNs
                    if isnan(betaW)
                        So.IPM=ones(1,78)*NaN;
                        So.Sig=ones(1,78)*NaN;
                        So.RDM=ones(1,66)*NaN;
                        So.psc_train = NaN;
                        So.psc_untrain = NaN;
                        Do.dist_train=NaN;
                        Do.dist_untrain=NaN;
                        Do.dist_cross=NaN;
                        Do.dist_all=NaN;
                        Do.eigTrain=ones(1,6)*NaN;
                        Do.eigUntrain=ones(1,6)*NaN;
                    else % run distance / RDM analyses
                        % crossval second moment matrix
                        [G,Sig]     = pcm_estGCrossval(betaW(1:(12*num_run),:),D.run,D.seqNumb);
                        So.IPM      = rsa_vectorizeIPM(G);
                        So.Sig      = rsa_vectorizeIPM(Sig);
                        % squared distances
                        So.RDM      = rsa.distanceLDC(betaW,D.run,D.seqNumb);
                        % add all subjects, regions - G values
                        G_all(:,:,g_ind)=G;
                        g_ind=g_ind+1;
                        % calculate G and dist
                        Do = calcDist(D,betaW,G);
                        
                        % make permutations
                        if simulations==1
                            P=randomSeqDist(D,betaW,G);
                            P.SN       = ones(length(P.dist_train),1)*s;
                            P.region   = ones(length(P.dist_train),1)*r;
                            P.regSide  = ones(length(P.dist_train),1)*regSide(r);
                            P.regType  = ones(length(P.dist_train),1)*regType(r);
                            % permutation struct
                            Po          = addstruct(Po,P);
                        end
                        % stats from additional images - remains nan if no data
                        So.psc_train = nanmean(S.psc_train{:});
                        So.psc_untrain = nanmean(S.psc_untrain{:});
                        %    So.surfdist_train = nanmean(S.searchlightdist_train{:});
                        %    So.surfdist_untrain = nanmean(S.searchlightdist_untrain{:});
                    end         
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
                    % data structure
                    To          = addstruct(To,So); % indexing fields, other images
                    To          = addstruct(To,Do); % distances
                         
                end; % each region
            end; % each subject
            
            if checkG==1
                % check the G across regions / subjects
                sml1_imana_dist('ROI_G_checkup',G_all,To);
            end
            
            % % save - stats data and simulations
            save(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)),'-struct','To');
            if simulations==1
                save(fullfile(regDir,sprintf('stats_%s_%sPW_SIMULATIONS_sess%d.mat',regType,betaChoice,ss)),'-struct','Po');
            end
            fprintf('\nDone.\n')
        end
        case 'ROI_G_checkup'
        G_all=varargin{1};
        To=varargin{2};
        
        ind=indicatorMatrix('allpairs',[1:12]);
        G_average = nanmean(G_all,3);
        G_group1 = nanmean(G_all(:,:,ismember(To.SN,[1:2:max(To.SN)])),3);
        G_group2 = nanmean(G_all(:,:,ismember(To.SN,[2:2:max(To.SN)])),3);
        
        figure
        subplot(1,3,1)
        imagesc(rsa.rdm.squareRDM(diag(ind*G_average*ind')));
        title('average RDM across subjects / regions');
        drawline(6.5,'dir','horz');drawline(6.5,'dir','vert');
        subplot(1,3,2)
        imagesc(rsa.rdm.squareRDM(diag(ind*G_group1*ind')));
        title('average RDM for group 1');
        drawline(6.5,'dir','horz');drawline(6.5,'dir','vert');
        subplot(1,3,3)
        imagesc(rsa.rdm.squareRDM(diag(ind*G_group2*ind')));
        title('average RDM for group 2');
        drawline(6.5,'dir','horz');drawline(6.5,'dir','vert');
        
        reg_indx=unique(To.region);
        for r = 1:length(reg_indx)
            G_reg = nanmean(G_all(:,:,ismember(To.region,reg_indx(r))),3);
            G_reg_G1 = nanmean(G_all(:,:,ismember(To.region,reg_indx(r))&ismember(To.SN,[1:2:max(To.SN)])),3);
            G_reg_G2 = nanmean(G_all(:,:,ismember(To.region,reg_indx(r))&ismember(To.SN,[2:2:max(To.SN)])),3);
            figure
            subplot(1,3,1)
            imagesc(rsa.rdm.squareRDM(diag(ind*G_reg*ind')));
            title(sprintf('region %s',regname{r}));
            drawline(6.5,'dir','horz');drawline(6.5,'dir','vert');
            subplot(1,3,2)
            imagesc(rsa.rdm.squareRDM(diag(ind*G_group1*ind')));
            title(sprintf('average RDM region %s for group 1',regname{r}));
            drawline(6.5,'dir','horz');drawline(6.5,'dir','vert');
            subplot(1,3,3)
            imagesc(rsa.rdm.squareRDM(diag(ind*G_group2*ind')));
            title(sprintf('average RDM region %s for group 2',regname{r}));
            drawline(6.5,'dir','horz');drawline(6.5,'dir','vert');
        end
        case 'ROI_PLOT_simulations'
            sessN=1;
            betaChoice='multi';
            vararginoptions(varargin,{'sessN','roi','betaChoice'});
            
            % load data and simulations
            D = load(fullfile(regDir,sprintf('stats_%sPW_sess%d.mat',betaChoice,sessN)));
            S = load(fullfile(regDir,sprintf('stats_%sPW_SIMULATIONS_sess%d.mat',betaChoice,sessN)));
            
            reg_indx=unique(D.region);
            
            for r=1:length(reg_indx)
                figure(1);
                subplot(1,length(reg_indx),r);
                histplot(S.dist_cross,'subset',S.region==r); % simulations
                drawline(mean(D.dist_cross(D.region==r)),'dir','vert','color',[1 0 0]);
                title(sprintf('%s',regname{r}));
                if r==1
                    ylabel('Simulation count'); xlabel('Trained vs. Untrained dist');
                end
                
                figure(2);
                subplot(1,length(reg_indx),r);
                histplot(S.dist_train,'subset',S.region==r); % simulations
                drawline(mean(D.dist_train(D.region==r)),'dir','vert','color',[1 0 0]);
                title(sprintf('%s',regname{r}));
                if r==1
                    ylabel('Simulation count'); xlabel('Trained dist');
                end
                
                figure(3);
                subplot(1,length(reg_indx),r);
                histplot(S.dist_untrain,'subset',S.region==r); % simulations
                drawline(mean(D.dist_untrain(D.region==r)),'dir','vert','color',[1 0 0]);
                title(sprintf('%s',regname{r}));
                if r==1
                    ylabel('Simulation count'); xlabel('Untrained dist');
                end
                
            end
            keyboard;            
    case 'BETA_stats_LOC'
        sessN = 1;
        sn  = [1:5];
        roi = [1:16];
        betaChoice = 'multi';
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice'});
        
        T = load(fullfile(regDir,sprintf('betas_LOC_sess%d.mat',sessN))); % loads region data (T)
        
        % output structures
        Ts = [];
        To = [];
        
        % do stats
        for s = sn % for each subject
            D = load(fullfile(glmLocSessDir{sessN}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
            fprintf('\nSubject: %d session: %d\n',s, sessN)
            num_run = numruns_loc_sess;
            
            for r = roi % for each region
                S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                fprintf('%d.',r)
                
                switch (betaChoice)
                    case 'uni'
                        betaW  = S.betaUW{1};
                    case 'multi'
                        betaW  = S.betaW{1};
                    case 'raw'
                        betaW  = S.betaRAW{1};
                end
                
                % % TseqNumb structure stats (finger mapping - 5 conditions)
                % crossval second moment matrix
                [G,Sig]     = pcm_estGCrossval(betaW(1:(5*num_run),:),D.run,D.seqNumb);
                H=eye(5)-ones(5,5)./5;  % centering matrix!
                G_cent      = H*G*H;  % double centered G matrix - rows and columns
                So.eig      = sort(eig(G_cent)','descend');
                So.IPM      = rsa_vectorizeIPM(G);
                So.Sig      = rsa_vectorizeIPM(Sig);
                % squared distances
                So.RDM_nocv = distance_euclidean(betaW',D.seqNumb)';
                So.RDM      = rsa.distanceLDC(betaW,D.run,D.seqNumb);
                % indexing fields
                So.SN       = s;
                So.region   = r;
                So.regSide  = regSide(r);
                So.regType  = regType(r);
                To          = addstruct(To,So);
                
            end; % each region
        end; % each subject

        % % save
        save(fullfile(regDir,sprintf('stats_LOC_%sPW_sess%d.mat',betaChoice,sessN)),'-struct','To');
        fprintf('\nDone.\n')      
       
    case 'ROI_pattern_consist'  
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 

        % enter sn, region,  beta: 0=betaW, 1=betaU, 2=raw betas
        % (1) Set parameters
        sn  = [4:9,11:25];
        roi = [1:8];
        sessN=1;
        betaChoice='multiPW'; % multiPW, uniPW, raw
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice'});
        
        numRun=numruns_task_sess;
        numCond=numel(num_seq);
        

        RR=[];
        %========%
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_Brodmann_sess%d_glm3',ss))); % loads in struct 'T'
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
                        
                        % calculate the pattern consistency
                        R.cnr                           = cnr_QC(beta,res,numCond,numRun);  % contrast to noise ratio
                        R.r2_rm                         = rsa_patternConsistency(beta,partVec,condVec,'removeMean',1);  % pattern consistency
                        R.r2                            = rsa_patternConsistency(beta,partVec,condVec,'removeMean',0);
                        [R.r2_cross_rm, R.r_cross_rm]   = rsa_patternConsistency_crossval(beta,partVec,condVec,'removeMean',1); % correlation of patterns
                        [R.r2_cross, R.r_cross]         = rsa_patternConsistency_crossval(beta,partVec,condVec,'removeMean',0);
                        
                        R.sn =s;
                        R.sessN=ss;
                        R.regType=r;
                        R.hemi=h;
                        
                        RR=addstruct(RR,R);
                    end
                end
            end
        end
        
        %save the structure
        save(fullfile(QCDir,sprintf('QualityControl_%sBetas',betaChoice)),'-struct','RR');
    case 'ROI_PLOT_pattern_consist' 
       var = 'cnr'; % cnr, r2_rm, r2, r2_cross, r_cross, r2_cross_rm, r_cross_rm
       betaChoice = 'multiPW';
       roi=[1:8];
       vararginoptions(varargin,{'var','betaChoice','roi'});
       
       T = load(fullfile(QCDir,sprintf('QualityControl_%sBetas',betaChoice)));
       
       % split by hemisphere
       figure
       subplot(3,4,[1:4]);
       plt.bar(T.regType,T.(var),'split',T.hemi,'leg',{'Contra','Ipsi'},'leglocation','northeast');
       ylabel(sprintf('%s %s',var,betaChoice));
       xlabel('ROIs');
       
       for rr=roi
           subplot(3,4,rr+4)
           plt.hist(T.(var),'subset',T.regType==rr,'split',T.hemi,'leg',{'contra','ipsi'},'leglocation','north');
           title(sprintf('%s',regname_cortex{rr}));
       end
       
       % split by session
     %  figure
     %  subplot(3,4,[1:4]);
     %  plt.bar(T.regType,T.(var),'subset',T.hemi==1,'split',T.sessN,'leg',{'sess1','sess2','sess3','sess4'},'leglocation','northeast');
     %  ylabel(sprintf('%s %s',var,betaChoice));
     %  xlabel('ROIs');
       
       for rr=roi
           subplot(3,4,rr+4)
           plt.hist(T.(var),'subset',T.hemi==1,'split',T.sessN,'leg',{'sess1','sess2','sess3','sess4'},'leglocation','north');
           title(sprintf('%s',regname_cortex{rr}));
       end
    
    case 'SURF_flatmap_rgb'
        sessN=1;
        smooth=0;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for h=1:2
            
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            if smooth == 1
                Cp=caret_load(['s' hem{h} sprintf('.summary_psc_sess%d.metric',sessN)]);
                Cd=caret_load([hem{h} sprintf('.summary_dist_sess%d.metric',sessN)]);
            else
                Cp=caret_load([hem{h} sprintf('.summary_psc_sess%d.metric',sessN)]);
                Cd=caret_load([hem{h} sprintf('.summary_dist_sess%d.metric',sessN)]);
            end
            flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
            topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
            
             %----PSC
            low_1=1;
            RGBdata(:,1)=Cp.data(:,1); % Red: trained BOLD signal
            RGBdata(:,3)=Cp.data(:,2); % Blue: untrained BOLD signal
            RGBdata(RGBdata(:,1)<low_1,1)=0;
            RGBdata(RGBdata(:,3)<low_1,3)=0;
            
            %------Dist
            low_2=0.002;               
            RGBdata(:,4)=Cd.data(:,2); % Red: trained distances
            RGBdata(:,6)=Cd.data(:,3); % Blue: untrained distances
            RGBdata(RGBdata(:,4)<low_2,4)=0; 
            RGBdata(RGBdata(:,6)<low_2,6)=0; 

            sc1=[low_1 3;low_1 3;low_1 3];                 % Scaling for PSC
            sc2=[low_2 0.04;low_2 0.04;low_2 0.04];        % Scaling for Distances
            
            name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN)};
            %name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN),sprintf('dist_cross_sess%d',sessN)};

            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1,sc2},'column_name',name);
            
            if smooth == 1
                caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.sess%d.RGB_paint',sessN)],C);
            else
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.sess%d.RGB_paint',sessN)],C);
            end
        end;
    case 'SURF_flatmap_trainedPsc_rgb'
        sessN=1;
        smooth=1;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for h=1:2
            
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            if smooth == 1
                Cp=caret_load(['s' hem{h} sprintf('.summary_psc_sess%d.metric',sessN)]);
            else
                Cp=caret_load([hem{h} sprintf('.summary_psc_sess%d.metric',sessN)]);
            end
            flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
            topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
            
             %----PSC - only trained sequences 
             %specific per session - sess1 yellow, sess2 orange, sess3 red
             
             switch sessN
                 case 1
                     rgb_indx = [255 175 14]/256;
                 case 2
                     rgb_indx = [242 95 3]/256;
                 case 3 
                     rgb_indx = [255 20 5]/256;
                 end
            low_1=0.9;
            
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
            
            name={sprintf('psc_trained_sess%d',sessN)};
            %name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN),sprintf('dist_cross_sess%d',sessN)};

            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1},'column_name',name);
            
            if smooth == 1
                caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.PSC_trained_sess%d.RGB_paint',sessN)],C);
            else
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.PSC_trained_sess%d.RGB_paint',sessN)],C);
            end
        end;
    case 'SURF_mapTessels'
        sn=[5:9,11:18];
        sessN=[1:4];
        parcelType='162tessels';
        betaChoice='multi';
        thres=[0 0.002];
        vararginoptions(varargin,{'sn','sessN','parcelType','betaChoice','thres'});
        for ss=sessN
            T=load(fullfile(regDir,sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
            for s=sn
                load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},parcelType))); % load region file into R
                for h=1:2
                    caretSubjDir=fullfile(caretDir,sprintf('x%s',subj_name{s}),hemName{h});
                    S=getrow(T,T.SN==s & T.regSide==h); % data from that hemisphere
                    C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
                    % per hemisphere
                    tesselNums=[1:148,150:158]; % exclude medial wall (149)
                    RGBdata=zeros(size(C.data,1),3);
                    for r=tesselNums, % we don't want 'FreeSurfer_Defined_Medial_Wall' - tessel numbers 149 for each hemisphere
                        indx=C.index(C.data==r);
                        % red - trained, blue - untrained
                        RGBdata(indx,1)=S.dist_train(r);
                        RGBdata(indx,3)=S.dist_untrain(r);
                    end
                    
                    RGBdata(RGBdata(:,1)<thres(1))=0;
                    RGBdata(RGBdata(:,3)<thres(1))=0;
                    
                    scale=[thres;thres;thres];
                    name={sprintf('dist_sess%d',ss)};
                    R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                    caret_save(fullfile(caretSubjDir,sprintf('%s.Dist_%s_sess%d.RGB_paint',hem{h},parcelType,ss)),R);
                end
            end
        end
    case 'SURF_mapTessels_group'
        %sn=[5,7,8,9,11,13,14,15,16];
        sn=[4,5,7,8,9,11,13,14,15,16,17,18]; % 1-3: different resolution, 6, 12- excessive movement, 10- no session 4
        sessN=[1:4];
        regType='162tessels';
        thres=[-0.002 0.002];
        vararginoptions(varargin,{'sn','sessN','parcelType','betaChoice','thres'});
        
        for ss=sessN
            for h=1:2
                count=0;
                for s=sn
                    count=count+1; % count of subjects
                    caretSubjDir=fullfile(caretDir,sprintf('x%s',subj_name{s}),hemName{h});
                    
                    S_rgb=caret_load(fullfile(caretSubjDir,sprintf('%s.Dist_%s_sess%d.RGB_paint',hem{h},regType,ss)));
                    
                    Group_x(:,count)=S_rgb.data(:,1);
                    Group_y(:,count)=S_rgb.data(:,2);
                    Group_z(:,count)=S_rgb.data(:,3);                  

                end
                S_group=S_rgb;
                % difference between the two
                S_group.data(:,1)=nanmean(Group_x,2);
                S_group.data(:,2)=nanmean(Group_x,2)-nanmean(Group_z,2);
                S_group.data(:,3)=nanmean(Group_z,2);
                scale=[thres;thres;thres];
                name={sprintf('dist_sess%d',ss)};
                S=caret_struct('RGBpaint','data',S_group.data,'scales',{scale},'column_name',name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.Dist_%s_sess%d.RGB_paint',hem{h},regType,ss)),S);
            end
        end
    case 'SURF_mapTessels_group_sessDiff'
        sn=[4,5,7,8,9,11,13,14,15,16,17,18];
        sessTr=[1:3];
        regType='162tessels';
        thres=[-0.002 0.002];
        vararginoptions(varargin,{'sn','sessN','parcelType','betaChoice','thres'});
        
        for ss=sessTr
            for h=1:2
                count=0;
                for s=sn
                    count=count+1; % count of subjects
                    caretSubjDir=fullfile(caretDir,sprintf('x%s',subj_name{s}),hemName{h});
                    
                    S_rgb1=caret_load(fullfile(caretSubjDir,sprintf('%s.Dist_%s_sess%d.RGB_paint',hem{h},regType,ss)));
                    S_rgb2=caret_load(fullfile(caretSubjDir,sprintf('%s.Dist_%s_sess%d.RGB_paint',hem{h},regType,ss+1)));
                    Group_x(:,count)=S_rgb2.data(:,1)-S_rgb1.data(:,1);
                    Group_y(:,count)=S_rgb2.data(:,2)-S_rgb1.data(:,2);
                    Group_z(:,count)=S_rgb2.data(:,3)-S_rgb1.data(:,3);                  

                end
                S_group=S_rgb1;
                 % difference between the two
                S_group.data(:,1)=nanmean(Group_x,2);
                S_group.data(:,2)=nanmean(Group_x,2)-nanmean(Group_z,2);
                S_group.data(:,3)=nanmean(Group_z,2);
                scale=[thres;thres;thres];
                name={sprintf('dist_sessTR%d-%d',ss,ss+1)};
                S=caret_struct('RGBpaint','data',S_group.data,'scales',{scale},'column_name',name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.Dist_%s_sessTR%d-%d.RGB_paint',hem{h},regType,ss,ss+1)),S);
            end
        end
    case 'Tessel_selection' % DEPRECIATED
        sn=[5,7,8,9,11,13,14,15,16];
        sessN=[1:4];
        regType='162tessels';
        thres=0;
        betaChoice='multi';
        vararginoptions(varargin,{'sn','sessN','parcelType','betaChoice','thres'}); 
        
        DD=[];
        for ss=sessN
            T=load(fullfile(regDir,sprintf('stats_%s_%sPW_sess%d.mat',regType,betaChoice,ss)));
            for s=sn
                load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},regType))); % load region file into R
                for h=1:2
                    S=getrow(T,T.SN==s & T.regSide==h); % data from that hemisphere
                    tesselNums=[1:148,150:158]; % exclude medial wall (149)

                    for r=tesselNums, % we don't want 'FreeSurfer_Defined_Medial_Wall' - tessel numbers 149 for each hemisphere
                        % red - trained, blue - untrained
                        D.dist(1,:)=S.dist_train(r);
                        D.dist(2,:)=S.dist_untrain(r);
                        D.seqType=[1;2];
                        D.sn=ones(size(D.seqType))*s;
                        D.hemi=ones(size(D.seqType))*h;
                        D.sessN=ones(size(D.seqType))*ss;
                        D.tesselNum=ones(size(D.seqType))*r;
                        DD=addstruct(DD,D);
                    end
                end
            end
        end
        
        keyboard;
    case 'SURF_mapTessel_selection'
        tesselNum=[1,46,16,49,105,104,43,44,45,47,52,103,48,50,2,20];
        vararginoptions(varargin,{'tesselNum'});
        
        for h=1:2;    % only contralateral
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            RGBdata=zeros(size(C.data,1),3);
            for r=tesselNum, % we don't want 'FreeSurfer_Defined_Medial_Wall' - tessel numbers 149 for each hemisphere
                indx=C.index(C.data==r);
                % mark in green
                RGBdata(indx,2)=1;
            end
            scale=[0 1; 0 1; 0 1];
            name={'tessel-selection'};
            R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
            caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tesselSelection.RGB_paint',hem{h})),R);
        end
        
    case 'Eig_Dist_relation'
        roi=[1:8];
        sessN=1;
        
        vararginoptions(varargin,{'roi','sessN'});
        
        for ss=sessN
            A=load(fullfile(regDir,sprintf('stats_multiPW_sess%d.mat',sessN)));
            figure
            scatterplot(sum(A.eigTrain,2),A.dist_train);
            xlabel('Sum of eigenvalues - trained'); ylabel('Average distance - trained');

            figure
            scatterplot(sum(A.eigUntrain,2),A.dist_untrain);
            xlabel('Sum of eigenvalues - untrained'); ylabel('Average distance - untrained');
                
        end 
        
        R1=sum(A.eigUntrain,2)./A.dist_untrain;
        R2=sum(A.eigTrain,2)./A.dist_train;
        fprintf('The sum of eigenvalues is %d times as much as the average distance.\n',R1(1));
        keyboard;   
    case 'CALC_eig_dim'
        % estimating the dimensionality of patterns 
        % cumulative sum of eig of G
        sn=[1:9,11:18];
        reg=[1:8];
        sessN=[1:4];
        betaChoice = 'multiPW';
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','reg','sessN','betaChoice','parcelType'});
        AA=[];
        for ss = sessN;    % sessions
            for s=sn
                for r=reg;
                    To = load(fullfile(regDir,sprintf('stats_%s_%s_sess%d.mat',parcelType,betaChoice,ss)));
                    %T = getrow(To,To.region==r);
                    T = getrow(To,To.region==r & To.SN==s);
                    A.eig(1,:) = T.eigTrain;
                    A.eig(2,:) = T.eigUntrain;
                    A.eig_sum(1,:) = cumsum(T.eigTrain,2);
                    A.eig_sum(2,:) = cumsum(T.eigUntrain,2);

                    negEig_train = sum(T.eigTrain<0);
                    negEig_untrain = sum(T.eigUntrain<0);
                    A.dim(1,:)=6-negEig_train-negEig_train;
                    A.dim(2,:)=6-negEig_untrain-negEig_untrain;
                    A.seqType=[1;2];
                    A.sn=ones(size(A.seqType))*s;
                    A.roi=ones(size(A.seqType))*r;
                    A.sessN=ones(size(A.seqType))*ss;
                    AA=addstruct(AA,A);
                end
            end
        end
        % save eigenvalues
        save(fullfile(distPscDir,'eigenvalues.mat'),'-struct','AA');
    case 'PLOT_eig_dim'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        seqType={'trained','untrained'};
        T=load(fullfile(distPscDir,'eigenvalues.mat'));
        
        AA=[];
        for t=1:size(T.eig,1)
            eigV=T.eig(t,:);
            A.dim=[1:6]';
            A.eig=eigV(:);
            A.seqType=ones(size(eigV'))*T.seqType(t);
            A.roi=ones(size(eigV'))*T.roi(t);
            A.sessN=ones(size(eigV'))*T.sessN(t);
            A.sn=ones(size(eigV'))*T.sn(t);
            AA=addstruct(AA,A);
        end
        
        for r=roi
            figure(r)
            for st=1:2
                subplot(1,2,st)
                plt.line(AA.dim,AA.eig,'split',AA.sessN,'subset',AA.roi==r&AA.seqType==st);
                title(sprintf('%s-%s',regname{r},seqType{st}));
                plt.match('y');
                drawline(0,'dir','horz');
            end
            
        end
    case 'STATS_eig_dimensions'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        seqType={'trained','untrained'};
        T=load(fullfile(distPscDir,'eigenvalues.mat'));
        
        AA=[];
        for t=1:size(T.eig,1)
            eigV=T.eig(t,:);
            A.dim=[1:6]';
            A.eig=eigV(:);
            A.seqType=ones(size(eigV'))*T.seqType(t);
            A.roi=ones(size(eigV'))*T.roi(t);
            A.sessN=ones(size(eigV'))*T.sessN(t);
            A.sn=ones(size(eigV'))*T.sn(t);
            AA=addstruct(AA,A);
        end
        
        for ss=sessN
            for r=roi
                fprintf('\n %s\n', regname{r});
                anovaMixed(AA.eig,AA.sn,'within',[AA.dim AA.seqType AA.sessN],{'dimension','seqType','session'},'subset',AA.roi==r);
                
            end
        end
        
    case 'ROI_dim_LOC'
        % estimating the dimensionality of patterns 
        % cumulative sum of eig of G
        sn=1;
        reg=2;
        betaChoice='multiPW';
        vararginoptions(varargin,{'sn','reg','betaChoice'});
        for s = 1:4;    % all sessions
            To = load(fullfile(regDir,sprintf('stats_LOC_%s_sess%d.mat',betaChoice,s)));
            T = getrow(To,To.region==reg & To.SN==sn);
            eig_fing = T.eig;
            eigfing_sum = cumsum(eig_fing);
             
            figure(1)
            subplot(2,2,s)
            title(sprintf('Session %d',s));
            plot(eigfing_sum,'-o','Color','b');
            legend('finger');
        end
    
    case 'CALC_dimensions'
        roi=3;
        sn=[4:9,11:22];
        sessN=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','sessN','parcelType'});
       
        Dim=[];
        for ss=sessN
            T=load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r=roi
                for s=sn
                    Tt=getrow(T,T.region==r & T.SN==s);
                    % get SPM info
                    TI=load(fullfile(glmSessDir{ss},subj_name{s},'SPM_info'));
                    for st=1:2 % trained / untrained
                        % calculate identity matrix for 6 sequences
                        ind=repmat(indicatorMatrix('identity',num_train),8,1);
                        % extract betas for seqType in question (st)
                        tmp=Tt.betaW{:};
                        betas=tmp(TI.seqType==st,:);
                        % extract condition / partition vectors
                        condVec=TI.seqNumb(TI.seqType==st);
                        partVec=TI.run(TI.seqType==st);
                        [D,dD]=crossval_estDim(betas,ind,num_train,partVec);
                        [corrC,dC]=crossval_estDimCorr(betas,ind,num_train,partVec);
                        
                        N.dim=D(:);
                        N.dimD=dD(:);
                        N.corr=corrC(:);
                        N.corrD=dC(:);
                        
                        N.indDim=[1:6]';
                        N.sn=ones(size(N.corr))*s;
                        N.roi=ones(size(N.corr))*r;
                        N.sessN=ones(size(N.corr))*ss;
                        N.seqType=ones(size(N.corr))*st;
                        
                        
                        Dim=addstruct(Dim,N);
                    end
                    fprintf('Done dimensionality: sess-%d %s %s \n',ss,regname{r},subj_name{s});
                end
            end
        end

        % save dimensionality
        save(fullfile(distPscDir,'Dimensionality.mat'),'-struct','Dim');
    case 'PLOT_dimensions'
        roi=3;
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'roi'});
        
        D=load(fullfile(distPscDir,'Dimensionality.mat'));
        D=getrow(D,D.roi==roi);
        figure
        for st=1:2
            subplot(1,2,st)
            plt.line(D.indDim,D.dimD,'split',D.sessN,'subset',D.seqType==st);
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s-%s',regname{roi},seqLabel{st}));
            if st==1
                ylabel('Dimensionality - old');
                xlabel('Possible dimensions');
            else
                ylabel('');
            end
        end
            
        figure
        for st=1:2
            subplot(1,2,st)
            plt.line(D.indDim,D.corrD,'split',D.sessN,'subset',D.seqType==st);
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s-%s',regname{roi},seqLabel{st}));
            if st==1
                ylabel('Dimensionality - correlation');
                xlabel('Possible dimensions');
            else
                ylabel('');
            end
        end
        
        figure
        for st=1:2
            subplot(1,2,st)
            plt.line(D.indDim,D.corr,'split',D.sessN,'subset',D.seqType==st);
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s-%s',regname{roi},seqLabel{st}));
            if st==1
                ylabel('Correlations');
                xlabel('Possible dimensions');
            else
                ylabel('');
            end
        end
         
    case 'ROI_psc_surfdist'
       sn = [1:7];
        roi = [1:8];
        sessN = 1:4;
        betaChoice = 'multiPW';
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sn','roi','seq','sessN','betaChoice','fig','parcelType'});

        Stats = [];
        
        for ss = sessN % do per session number
            T = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            
            for s=1:numel(sn)
                for r=roi
                    
                    S.sn=s;
                    S.roi=r;
                    S.surfdist_train=nanmean(T.searchlightdist_train{(T.region==r & T.SN==sn(s)),:});
                    S.surfdist_untrain=nanmean(T.searchlightdist_untrain{(T.region==r & T.SN==sn(s)),:});
                    S.psc_train=nanmean(T.psc_train{(T.region==r & T.SN==sn(s)),:});
                    S.psc_untrain=nanmean(T.psc_untrain{(T.region==r & T.SN==sn(s)),:});
                    S.sessN=ss;
                    
                    Stats=addstruct(Stats,S);
                end
            end
            
        end
        
        
        Stats.sessIndx = zeros(size(Stats.sessN));
        Stats.sessIndx(Stats.sessN<4)=1;
        Stats.sessIndx(Stats.sessN==4)=2;
        
        a = [Stats.sessIndx; Stats.sessIndx];
        b = [Stats.sessN; Stats.sessN];
        
        keyboard;
        figure
        for f = 1:numel(roi)
            subplot(1,numel(roi),f);
            lineplot([a b],[Stats.psc_untrain;Stats.psc_train],'split',[ones(length(Stats.sessN),1);ones(length(Stats.sessN),1)*2],'style_thickline','subset',[Stats.roi;Stats.roi]==f,'leg',{'untrained','trained'});
            ylim([0 2])
            if f==1
                ylabel('Percent signal change')
            else
                ylabel('')
            end
            title(sprintf('%s',regname{f}))
        end
        
        
        figure
        for f = 1:numel(roi)
            subplot(1,numel(roi),f);
            lineplot([a b],[Stats.surfdist_untrain;Stats.surfdist_train],'split',[ones(length(Stats.sessN),1);ones(length(Stats.sessN),1)*2],'style_thickline','subset',[Stats.roi;Stats.roi]==f,'leg',{'untrained','trained'});
            ylim([-0.002 0.002])
            if f==1
                ylabel('Distances')
            else
                ylabel('')
            end
            title(sprintf('%s',regname{f}))
        end
      
        
        keyboard;
    
    case 'SAVE_dist'
        sn = [1:9,11:18];
        roi = [1:8];
        sessN = [1:4];
        betaChoice='multiPW';
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice','fig','parcelType'});
        
        Stats = [];
        for ss=sessN
            D = load(fullfile(regDir,sprintf('stats_%s_%s_sess%d.mat',parcelType,betaChoice,ss))); % loads region data (D)
            for s=sn
                for r=roi
                    S.dist_train=D.dist_train(D.SN==s & D.region==r);
                    S.dist_untrain=D.dist_untrain(D.SN==s & D.region==r);
                    S.dist_cross=D.dist_cross(D.SN==s & D.region==r);
                    S.sn=s;
                    S.roi=r;
                    S.sessN=ss;
                    Stats=addstruct(Stats,S); 
                end
            end
        end
        % save structure
        save(fullfile(distPscDir,'dist_ROI.mat'),'-struct','Stats');
    case 'PLOT_dist'
        roi=[1:8];
        sessN=[1:4];
        
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,'dist_ROI.mat'));
        
        % create structure with trained / untrained together, seqType
        D.dist=[T.dist_train;T.dist_untrain];
        D.seqType=[ones(size(T.dist_train));ones(size(T.dist_train))*2];
        D.sn=[T.sn;T.sn];
        D.roi=[T.roi;T.roi];
        D.sessN=[T.sessN;T.sessN];
        
        indx=1;
        figure
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.dist),'split',D.seqType,'subset',D.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.dist),'split',D.seqType,'subset',D.roi==r&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==2
                plt.line(D.sessN,ssqrt(D.dist),'split',D.seqType,'subset',D.roi==r&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Mahalanobis distance');
                xlabel('Session');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'PLOT_dist_sessions'
        sessN=[1:3];
        regExcl=6;
        seqType='trained';
        vararginoptions(varargin,{'sessN','regExcl','seqType'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        switch (seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        T=load(fullfile(distPscDir,'dist_ROI.mat'));
        
        D.dist=[T.dist_train;T.dist_untrain];
        D.seqType=[ones(size(T.dist_train));ones(size(T.dist_train))*2];
        D.sn=[T.sn;T.sn];
        D.roi=[T.roi;T.roi];
        D.sessN=[T.sessN;T.sessN];
        
        T=getrow(D,~ismember(D.roi,regExcl)&ismember(D.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        figure
        plt.line(T.reg,T.dist,'split',T.sessN,'subset',T.seqType==st,'leg',sessLegend(sessN),'leglocation','northeast');
        ylabel('Distances');
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI');
    case 'PLOT_dist_seqType'
        roi=[1:8];
        sessN=[1:4];
        
        vararginoptions(varargin,{'roi','sessN'});
        
        D=load(fullfile(distPscDir,'dist_ROI.mat'));
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.dist_cross),'subset',D.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.dist_cross),'subset',D.roi==r&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==2
                plt.line(D.sessN,ssqrt(D.dist_cross),'subset',D.roi==r&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Mahalanobis distance');
                xlabel('Session');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'PLOT_dist_seqType_sessions'
        sessN=[1:3];
        regExcl=6;
        vararginoptions(varargin,{'sessN','regExcl','seqType'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(distPscDir,'dist_ROI.mat'));
        
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        figure
        plt.line(T.reg,T.dist_cross,'split',T.sessN,'leg',sessLegend(sessN),'leglocation','northeast');
        ylabel('Distances');
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI');
        
    case 'STATS_dist'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        
        D=load(fullfile(distPscDir,'dist_ROI.mat'));
        T.dist=[D.dist_train; D.dist_untrain];
        T.seqType=[ones(size(D.dist_train));ones(size(D.dist_train))*2];
        T.sn=[D.sn;D.sn];
        T.roi=[D.roi;D.roi];
        T.sessN=[D.sessN;D.sessN];
        
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.dist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r);
        end
        keyboard;
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.dist,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);
            end
        end
    case 'STATS_dist_seqType'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,'dist_ROI.mat'));
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Distance seqType in %s \n',regname{r});
            if numel(sessN)==4
                anovaMixed(T.dist_cross,T.sn,'within',[T.sessN] ,{'session'},'subset',T.roi==r);
            elseif numel(sessN)==3
                anovaMixed(T.dist_cross,T.sn,'within',[T.sessN] ,{'session'},'subset',T.roi==r&T.sessN<4);
            elseif numel(sessN)==2
                anovaMixed(T.dist_cross,T.sn,'within',[T.sessN] ,{'session'},'subset',T.roi==r&(T.sessN==1|T.sessN==4));
            end
        end

        keyboard;   
        
        for r=roi
            for sTr=1:3
                fprintf('\n post-hoc t-test on the effect of seqType in transition sess%d-%d in %s \n',sTr,sTr+1,regname{r});
                ttestDirect(T.dist_cross,[T.sessN T.sn],2,'paired','subset',T.roi==r&T.sessN==sTr|T.sessN==sTr+1);
            end
        end    
        
    case 'CALC_corrdist'
        reg = [1:8];
        sn  = [1:9,11:18];
        sessN = [1:4];
        parcelType='Brodmann';
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        vararginoptions(varargin,{'sn','reg','sessN','subtract_mean','parcelType'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
            D   = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for roi = reg;
                for s = sn;
                    SI = load(fullfile(glmSessDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                    
                    t   = getrow(D,D.region==roi & D.SN==s);
                    
                    data = t.betaW{1}(1:size(SI.SN,1),:);
                    
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
                    
                    % crossval second moment matrix
                    
                    [G,Sig]     = pcm_estGCrossval(data,SI.run,SI.seqNumb);
                    C=corr_crossval(G,'reg','minvalue');
                    C=rsa_squareRDM(C);
                    
                    % average trained dist
                    S.corrDist(1,:) = sum(sum(triu(C(1:6,1:6))))/(6*5/2);
                    % average untrained dist
                    S.corrDist(2,:) = sum(sum(triu(C(7:12,7:12))))/(6*5/2);
                    S.seqType=[1;2]; % trained, untrained
                    S.sessN=[ss;ss];
                    S.roi=[roi;roi];
                    S.sn=[s;s];
                    SAll=addstruct(SAll,S);
                    
                    ST.corr_seqType = sum(sum(triu(C(7:12,1:6))))/(6*5);
                    ST.sessN=ss;
                    ST.roi=roi;
                    ST.sn=s;
                    STAll=addstruct(STAll,ST);
                end
            end
        end
        save(fullfile(distPscDir,'corrDist_ROI.mat'),'-struct','SAll');
        save(fullfile(distPscDir,'corrDist_seqType_ROI.mat'),'-struct','STAll');
    case 'PLOT_corrDist'
        roi=[1:8];
        sessN=[1:4];
        
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,'corrDist_ROI.mat'));
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([T.sessN>3 T.sessN],T.corrDist,'split',T.seqType,'subset',T.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(T.sessN,T.corrDist,'split',T.seqType,'subset',T.roi==r&T.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Correlation distance');
                xlabel('Session');
            else
                ylabel('');
            end
        indx=indx+1;    
        end
    case 'PLOT_corrDist_sessions'
        sessN=[1:3];
        regExcl=6;
        seqType='trained';
        vararginoptions(varargin,{'sessN','regExcl','seqType'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        switch (seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        T=load(fullfile(distPscDir,'corrDist_ROI.mat'));
        
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        figure
        plt.line(T.reg,T.corrDist,'split',T.sessN,'subset',T.seqType==st,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        drawline(0,'dir','horz');
        ylabel('Correlation distances');
        title(sprintf('%s',seqType))
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI');
    case 'PLOT_corrDist_sessions_both'
        sessN=[1:3];
        regExcl=6;
        seqType={'trained','untrained'};
        vararginoptions(varargin,{'sessN','regExcl','seqType'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(distPscDir,'corrDist_ROI.mat'));
        
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        for st=1:2
            if st==1
                sty=styTrained_sess;
            else
                sty=styUntrained_sess;
            end
            figure(1)
            plt.line(T.reg,T.corrDist,'split',T.sessN,'subset',T.seqType==st,'leg',sessLegend(sessN),'leglocation','northeast','style',sty);
            hold on;
            drawline(0,'dir','horz');
            ylabel('Correlation distances');
            title(sprintf('%s',seqType{st}))
            set(gca,'XTickLabel',regname_cortex(reg));
            xlabel('ROI');
            
        end
    case 'PLOT_corrDist_seqType'
      roi=[1:8];
        sessN=[1:4];
        
        vararginoptions(varargin,{'roi','sessN'});
        
        D=load(fullfile(distPscDir,'corrDist_seqType_ROI.mat'));
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],D.corr_seqType,'subset',D.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==3
                plt.line(D.sessN,D.corr_seqType,'subset',D.roi==r&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==2
                plt.line(D.sessN,D.corr_seqType,'subset',D.roi==r&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Correlation distance betweeen seqTypes');
                xlabel('Session');
            else
                ylabel('');
            end
        indx=indx+1;    
        end
    case 'PLOT_corrDist_seqType_sessions'
        sessN=[1:3];
        regExcl=6;
        vararginoptions(varargin,{'sessN','regExcl','seqType'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};

        T=load(fullfile(distPscDir,'corrDist_seqType_ROI.mat'));
        
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        figure
        plt.line(T.reg,T.corr_seqType,'split',T.sessN,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        ylabel('Correlation distances');
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI');
    case 'STATS_corrDist'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,'corrdist_ROI.mat'));

        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r&T.sessN<4);
            %anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r&(T.sessN==1|T.sessN==4));
        end
        
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.corrDist,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);
            end
        end
    case 'STATS_corrDist_seqType'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,'corrdist_seqType_ROI.mat'));
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Distance seqType in %s \n',regname{r});
            if numel(sessN)==4
                anovaMixed(T.corr_seqType,T.sn,'within',[T.sessN] ,{'session'},'subset',T.roi==r);
            elseif numel(sessN)==3
                anovaMixed(T.corr_seqType,T.sn,'within',[T.sessN] ,{'session'},'subset',T.roi==r&T.sessN<4);
            elseif numel(sessN)==2
                anovaMixed(T.corr_seqType,T.sn,'within',[T.sessN] ,{'session'},'subset',T.roi==r&(T.sessN==1|T.sessN==4));
            end
        end
  
        for r=roi
            for sTr=1:3
                fprintf('\n post-hoc t-test on the effect of seqType in transition sess%d-%d in %s \n',sTr,sTr+1,regname{r});
                ttestDirect(T.corr_seqType,[T.sessN T.sn],2,'paired','subset',T.roi==r&T.sessN==sTr|T.sessN==sTr+1);
            end
        end
    
    case 'RAND_GEN_corrDist'
        numPerm=1000; % number of permutations
        numVox=1000; % number of voxels
        SS=[];
        
        partVec=kron([1:8]',ones(12,1));
        condVec=kron(ones(8,1),[1:12]');
        
        for p=1:numPerm
            Data=randn(96,numVox);
              [G,Sig]     = pcm_estGCrossval(Data,partVec,condVec);
              C=corr_crossval(G,'reg','minvalue');
              C=rsa_squareRDM(C);
              
              % average trained dist
              S.corrDist1 = sum(sum(triu(C(1:6,1:6))))/(6*5/2);
              % average untrained dist
              S.corrDist2 = sum(sum(triu(C(7:12,7:12))))/(6*5/2);
              S.corr_seqType = sum(sum(triu(C(7:12,1:6))))/(6*5);
              S.numPerm=p;       
              SS=addstruct(SS,S);
        end
        
        keyboard;
        
    case 'SAVE_psc'
        sn = [1:9,11:18];
        roi = [1:8];
        sessN = [1:4];
        parcelType='Brodmann';

        vararginoptions(varargin,{'sn','roi','seq','sessN','parcelType'});

        Stats = [];
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            for s=sn
                for r=roi
                   
                    S.psc(1,:)=nanmean(T.psc_train{(T.region==r & T.SN==s),:});
                    S.psc(2,:)=nanmean(T.psc_untrain{(T.region==r & T.SN==s),:});
                    S.seqType=[1;2];
                    S.sn=ones(size(S.seqType))*s;
                    S.roi=ones(size(S.seqType))*r;
                    S.sessN=ones(size(S.seqType))*ss;
                    
                    Stats=addstruct(Stats,S);
                end
            end
        end
           
        % save structure
        save(fullfile(distPscDir,'psc_ROI.mat'),'-struct','Stats');
    case 'PLOT_psc'
        roi=[1:8];
        sessN=[1:4];
        
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,'psc_ROI.mat'));
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([T.sessN>3 T.sessN],T.psc,'split',T.seqType,'subset',T.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(T.sessN,T.psc,'split',T.seqType,'subset',T.roi==r&T.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            end
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if r==1
                ylabel('Percent signal change');
                xlabel('Session');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'PLOT_psc_sessions'
        sessN=[1:3];
        regExcl=6;
        seqType='trained';
        vararginoptions(varargin,{'sessN','regExcl','seqType'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        switch (seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        T=load(fullfile(distPscDir,'psc_ROI.mat'));
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN)); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.roi;
        for i=1:length(T.roi) 
            T.reg(i,:)=find(T.reg(i,:)==reg);
        end
        
        figure
        plt.line(T.reg,T.psc,'split',T.sessN,'subset',T.seqType==st,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
        ylabel('Percent signal change');
        set(gca,'XTickLabel',regname_cortex(reg));
        xlabel('ROI');
    case 'STATS_psc'
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,'psc_ROI.mat'));
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n PSC in %s \n',regname{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r);
        end
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.psc,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);
            end
        end
    
    case 'SURF_dist'
        % number of vertices above a certain threshold
        vararginoptions(varargin,{'sn'});
        thres = 0.001;
        T=[];
        for s=sn;
            for ss = 1:sess_sn(s)
                for hem = 1:2
                    
                    caretSubjDir = fullfile(caretDir,sprintf('x%s',subj_name{s}),hemName{hem});
                    cd(caretSubjDir);
                    surf=caret_load(sprintf('%s_sess%d_dist.metric',subj_name{s},ss));
                    C = caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},['ROI.paint']));
                    dist_train=surf.data(:,2);
                    dist_untrain=surf.data(:,3);
                    A.hem(1:2,:)=hem;
                    A.sn(1:2,:)=s;
                    A.sessN(1:2,:)=ss;
                    A.vert(1,:)=sum(dist_train>thres);
                    A.vert(2,:)=sum(dist_untrain>thres);
                    A.dist_pattern(1,:)=dist_train';
                    A.dist_pattern(2,:)=dist_untrain';
                    A.seqType=[1,2]';
                    T=addstruct(T,A);
                end
            end
            
            B_train=getrow(T,T.sn==s & T.hem==1 & T.seqType==1);
            B_untrain=getrow(T,T.sn==s & T.hem==1 & T.seqType==2);
            for ss1=1:4
                for ss2=1:4
                    Corr_train(ss1,ss2)=corr(B_train.dist_pattern(ss1,:)',B_train.dist_pattern(ss2,:)','rows','pairwise');
                    Corr_untrain(ss1,ss2)=corr(B_untrain.dist_pattern(ss1,:)',B_untrain.dist_pattern(ss2,:)','rows','pairwise');
                end
            end
            Corr.all((s-1)*6+1,:)=Corr_train(1,2); % sess1-2
            Corr.all((s-1)*6+2,:)=Corr_train(2,3); % sess2-3
            Corr.all((s-1)*6+3,:)=Corr_train(3,4); % sess3-4
            Corr.all((s-1)*6+4,:)=Corr_untrain(1,2);
            Corr.all((s-1)*6+5,:)=Corr_untrain(2,3);
            Corr.all((s-1)*6+6,:)=Corr_untrain(3,4);
            Corr.sess([(s-1)*6+1,(s-1)*6+4],:)=1;
            Corr.sess([(s-1)*6+2,(s-1)*6+5],:)=2;
            Corr.sess([(s-1)*6+3,(s-1)*6+6],:)=3;
            Corr.indx((s-1)*6+1:(s-1)*6+3,:)=1; %train
            Corr.indx((s-1)*6+4:(s-1)*6+6,:)=2; %untrain
            Corr.sn((s-1)*6+1:(s-1)*6+6,:)=s; %sn
        end
        keyboard;
        figure;
        lineplot(T.sessN,T.vert,'split',T.seqType,'subset',T.hem==1,'style_thickline','leg',{'train','untrain'});
        title('All subjects')
        
        figure;
        for s=1:numel(sn)
            subplot(1,numel(sn),s)
            lineplot(T.sessN,T.vert,'split',T.seqType,'subset',T.hem==1 & T.sn==s,'style_thickline','leg',{'train','untrain'});
            title(sprintf('subject %d',s));
        end
        
        figure;
        for s=1:numel(sn)
            subplot(1,numel(sn),s)
            barplot(Corr.sess,Corr.all,'split',Corr.indx,'subset',Corr.sn==s,'leg',{'train','untrain'});
            title(sprintf('subject %d',s));
        end
    
    case 'PLOT_psc_speed'
        % plot relative increase in psc for trained / untrained - sess3/4
        roi=[1:5,7:8];
        sessN=[3,4];
        operation='subtract'; % divide or subtract

        vararginoptions(varargin,{'roi','sessN','operation'});

        T=load(fullfile(distPscDir,'psc_ROI.mat'));

        T1=getrow(T,T.sessN==3&T.roi~=6&T.sn~=12);
        T2=getrow(T,T.sessN==4&T.roi~=6&T.sn~=12);
        switch operation
            case 'subtract'
                S.psc=T2.psc-T1.psc;
            case 'divide'
                 S.psc=T2.psc./T1.psc;
        end

        S.roi=T2.roi;
        S.sn=T2.sn;
        S.seqType=T2.seqType;

        figure
        plt.bar(S.roi,S.psc,'split',S.seqType,'style',stySeq,'leg',{'trained','untrained'});
        drawline(1,'dir','horz');
        xlabel('sequence type'); ylabel(sprintf('Psc 4th %s 3rd session',operation));
        set(gca,'XTick',[1.5 4.5 8 11 14.5 18 21],'XTickLabel',regname(roi));
        keyboard;
    case 'PLOT_corrDist_speed'
        % plot relative increase in corrDist for trained / untrained - sess3/4
        roi=[1:8];
        sessN=[3,4];
        operation='subtract'; % divide or subtract
        
        vararginoptions(varargin,{'roi','sessN','operation'});
        
        T=load(fullfile(distPscDir,'corrDist_ROI.mat')); 
        
        T1=getrow(T,T.sessN==3&T.roi~=6&T.sn~=12);
        T2=getrow(T,T.sessN==4&T.roi~=6&T.sn~=12);
        switch (operation)
            case 'subtract'
                S.corrDist=T2.corrDist-T1.corrDist;
            case 'divide'
                 S.corrDist=T2.corrDist./T1.corrDist;
        end     
        S.roi=T2.roi;
        S.sn=T2.sn;
        S.seqType=T2.seqType;
        
        figure
        plt.bar(S.roi,S.corrDist,'split',S.seqType,'style',stySeq,'leg',{'trained','untrained'});
        drawline([0 1],'dir','horz');
        xlabel('sequence type'); ylabel(sprintf('Corr dist 4th %s 3rd session',operation));
        set(gca,'XTick',[1.5 4.5 8 11 14.5 18 21],'XTickLabel',regname(roi));
        keyboard;   
    case 'PLOT_psc_dist_speed'
        % plot relative increase in corrDist relative to increase in pscfor trained / untrained - sess3/4
        roi=[1:8];
        sessN=[3,4];
        operation='divide';
        
        vararginoptions(varargin,{'roi','sessN','operation'});
        
        T=load(fullfile(distPscDir,'dist_ROI.mat')); 
        P=load(fullfile(distPscDir,'psc_ROI.mat'));
        D.dist=[T.dist_train;T.dist_untrain];
        D.sessN=[T.sessN;T.sessN];
        
        D1=getrow(D,D.sessN==3);
        D2=getrow(D,D.sessN==4);
        P1=getrow(P,P.sessN==3);
        P2=getrow(P,P.sessN==4);
        switch operation
            case 'divide'
                S.psc=P2.psc./P1.psc;
                S.dist=D2.dist./D1.dist;
            case 'subtract'
                S.psc=P2.psc-P1.psc;
                S.dist=D2.dist-D1.dist;
        end

        S.distPsc=S.dist./S.psc;
        S.roi=P2.roi;
        S.sn=P2.sn;
        S.seqType=P2.seqType;
        
        figure
        plt.bar(S.roi,S.distPsc,'split',S.seqType,'style',stySeq,'leg',{'trained','untrained'});
        drawline(1,'dir','horz');
        xlabel('sequence type'); ylabel('Relative increase in dist to activation 4th / 3rd session');
         set(gca,'XTick',[1.5 4.5 8 11 14.5 18 21],'XTickLabel',regname(roi));
        keyboard;
        
    case 'PCM_simulateModelFamily'
        runEffect  = 'random';
        noise=1;
        scale=0.01;
        theta = [-20; -20; 0];
        fig=1;

        vararginoptions(varargin,{'runEffect','theta','noise','scale','fig'});
        
        numSim=10; % number of simulations
        % model specifications

        numPart = 8;
        D.numVox  = 1000;
        sn = 8; % 8 subjects
        K=[];
        KK=[];
        
        % subset of models 
        m = pcm_defineSequenceModels(Seq,1);
        mm{1}=m{1};
        mm{2}=m{2};
        mm{3}=m{6};
        [M Z] = pcm_buildModelFromFeatures(mm,'style','encoding_style','type','component');
        [Mf,Comb] = pcm_constructModelFamily(M,'alwaysInclude',1,'fullModel',1);
        
        
        D.partVec=kron([1:numPart]',ones(size(Z,1),1));
        D.condVec=kron(ones(numPart,1),Z);
        
        for i=1:numSim
            Data = pcm_generateData_eva(Mf{end},theta,D,sn,scale,noise,'design',Z);
          %  [Data,part,cond] = pcm_generateData_eva(Mf{end},theta,D,sn,scale,noise,'design',Z);
            T = pcm_fitModels(Data,Mf,D.partVec,D.condVec,runEffect);
            
            % calculate posterior probability of each component
            [postProp,logBayes]=pcm_componentPosterior(T.cross_likelihood,Comb);
            
            % knockIN, knockOUT
            [knockIN,knockOUT]=pcm_knockModels(T.cross_likelihood,Comb);
            
            K.postProp=postProp(:);
            K.logBayes=logBayes(:);
            K.knockIN=knockIN(:);
            K.knockOUT=knockOUT(:);
            K.indx=[ones(length(postProp),1);ones(length(postProp),1)*2;ones(length(postProp),1)*3];
            K.numSim=ones(length(K.indx),1)*i;
            KK=addstruct(KK,K);
        end

        if fig==1
        figure
        subplot(2,2,1)
        histplot(KK.postProp,'split',KK.indx);
        legend({'feat1','feat2','feat3'});
        xlabel('Posterior probability value');
        ylabel('Proportion of simulations');
        subplot(2,2,2)
        histplot(KK.logBayes,'split',KK.indx);
        legend({'feat1','feat2','feat3'});
        xlabel('Log Bayes factor');
        ylabel('Proportion of simulations');
        subplot(2,2,3)
        histplot(KK.knockIN,'split',KK.indx);
        legend({'feat1','feat2','feat3'});
        xlabel('Knock IN value');
        ylabel('Proportion of simulations');
        subplot(2,2,4)
        histplot(KK.knockOUT,'split',KK.indx);
        legend({'feat1','feat2','feat3'});
        xlabel('Knock OUT value');
        ylabel('Proportion of simulations');
        end
        keyboard;
    case 'PCM_seq_ROI_sess' % WORK ON THIS!!!
        runEffect  = 'fixed';
        beta_choice = 'mw';
        reg = [1:8];
        sn=[1:8];
        sessN=[1:4];
        AllReg=[];
        parcelType='Brodmann';
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','parcelType'})
        
        for ss = sessN
            B=load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            
            for r = reg
                
                for p=1:length(sn)
                    
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==p&B.region==r)}';
                        case 'mw'
                            beta = B.betaW{(B.SN==p&B.region==r)}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==p&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                    end
          
                    partVec{p} = D.run;  % runs/partitions
                    condVec{p} = D.seqNumb;
                    indx = ones(size(D.run));
                    
                    Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept

                end;
                
                m = pcm_defineSequenceModels(Seq);
                % build all possible models
                M = pcm_buildFeatureModel('test','indicator',m{1},m{2},m{9});
                T = pcm_fitModels(Data,M,partVec,condVec,runEffect);
                
                T.roi = ones(8,1)*r;
                T.sessN = ones(8,1)*ss;
                AllReg=addstruct(AllReg,T);
                
            end
        end
        
        keyboard;
        figure
        s=style.custom({'orange','green'},'markersize',12);
        for sess=1:4
        subplot(1,4,sess)
        A=getrow(AllReg,AllReg.sessN==sess);
        plt.line([A.roi; A.roi],[A.bayesEst(:,2);A.bayesEst(:,3)],'split',[ones(size(A.roi));ones(size(A.roi))*2],'style',s,'leg',{'SeqType','SeqNumb'});
        plt.labels('ROIs','log likelihood',(sprintf('session %d',sess)),[],[]);
        set(gca,'XTickLabel',regname(1:8));
        %lineplot(AllReg.reg_ind,AllReg.bayesEst(:,2),'linecolor',[1 0 0],'markerfill');
        end 
        plt.match('y');
        keyboard;
    case 'PCM_constructModelFamily' % work on the SPECIFIC models
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        parcelType='162tessels'; % Brodmann or 162tessels
        reg = [1,2,16,20,43,44,45,46,47,48,49,50,52,103,104,105]; % for tessels, 1:8 for Brodmann
        hemi=[1:2];
        sn=[4,5,7,8,9,11,13,14,15,16]; % change
        sessN=[1:4];
        AllReg=[];
        K=[];
        KK=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi'})
        tic;
        for ss = sessN
            B=load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for h=hemi
                for r = reg
                    for p=1:length(sn)
                        
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
                        
                        partVec{p} = D.run;  % runs/partitions
                        switch(runEffect)
                            case 'fixed'
                                %m = pcm_defineSequenceModels_fixed(Seq,p);
                                m = pcm_defineSequenceModels_new(Seq,p);
                                [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                                [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                            case 'random'
                                m = pcm_defineSequenceModels(Seq,p);
                                [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                                [Mf,Comb] = pcm_constructModelFamily(M,'alwaysInclude',1,'fullModel',1);
                        end                       
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = ones(size(D.run));
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept                        
                    end;
                    
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*r;
                    T.hemi = ones(size(T.SN))*h;
                    T.sessN = ones(size(T.SN))*ss;
                    AllReg=addstruct(AllReg,T);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*r;
                    K.hemi = ones(length(K.indx),1)*h;
                    K.sessN = ones(length(K.indx),1)*ss;
                    KK=addstruct(KK,K);
                    
                end
            end
        end
        
        keyboard;
        % save variables;
        dircheck(fullfile(pcmDir));
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_NEW_%s.mat',parcelType)),'Comb');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_NEW_%s.mat',parcelType)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_NEW_%s.mat',parcelType)),'-struct','KK');
        name = {'ST','TS','US','FF','FA','FT'};
    case 'PCM_PLOT_modelFamily'
        %plots knock-IN / OUT / posterior probability
        parcelType='Brodmann';
        reg=1;
        vararginoptions(varargin,{'reg','parcelType'});
        KK=load(fullfile(pcmDir,sprintf('ModelFamily_%s_Stats.mat',parcelType)));
        name = {'ST','TS','US','FF','FA','FT'};
        fig_num=0;
        for r = reg
            % knock-in
            figure(fig_num+1);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.knockIN);
                
                if s==1
                    ylabel('Log Bayes - knock IN');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
                set(gca,'XTickLabel',name);
                drawline(0,'dir','horz');
                plt.match('y');
            end
            % knock-out
            figure(fig_num+2);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.knockOUT);
                if s==1
                    ylabel('Log Bayes - knock OUT');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
                set(gca,'XTickLabel',name);
                drawline(0,'dir','horz');
                plt.match('y');
            end
            % posterior probability
            figure(fig_num+3);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.postProp);             
                set(gca,'XTickLabel',name);
                drawline(0.5,'dir','horz');
                plt.match('y');
                if s==1
                    ylabel('Posterior probability');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
            % posterior Log
            figure(fig_num+4);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.logBayes);
                set(gca,'XTickLabel',name);
                drawline(0,'dir','horz');
                plt.match('y');
                if s==1
                    ylabel('Log Bayes - post prob');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
            fig_num=fig_num+5;  
        end    
    case 'PCM_construct_ModelStruct'  

        parcelType='Brodmann';
        vararginoptions(vararin,{'parcelType'});
        M=load(fullfile(pcmDir,sprintf('ModelFamily_%s_Fit.mat',parcelType)));
        load(fullfile(pcmDir,sprintf('ModelFamilyComb_%s.mat',parcelType))); % loads combination of models as Comb
        name = {'ST','TS','US','FF','FA','FT'};  
        
        TT=[];
        % construct a structure where every line is one model
        for i=1:size(M.bayesEst,1)
            for m=1:size(M.bayesEst,2)
                T.sn=M.SN(i);
                T.bayesEst=M.bayesEst(i,m);
                T.modelComb=Comb(m,:);
                T.modelIndx=m;
                T.roi=M.roi(i);
                T.sessN=M.sessN(i);
                TT=addstruct(TT,T);
            end
        end
        
        save(fullfile(pcmDir,sprintf('PCM_AllModelsFamily_%s.mat',parcelType)),'-struct','TT');
   
    case 'STATS_PCM_modelFamily'
        roi=[1:8];
        sessN=[1:4];
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','sessN','parcelType'});
       
        T=load(fullfile(pcmDir,sprintf('ModelFamily_%s_Stats.mat',parcelType)));
        
        for r=roi
            fprintf('\n logBayes factor in %s\n',regname{r});
            anovaMixed(T.logBayes,T.sn,'within',[T.sessN T.indx],{'session','model'},'subset',T.roi==r);
            fprintf('\n knockIN factor in %s\n',regname{r});
            anovaMixed(T.knockIN,T.sn,'within',[T.sessN T.indx],{'session','model'},'subset',T.roi==r);
            fprintf('\n knockOUT factor in %s\n',regname{r});
            anovaMixed(T.knockOUT,T.sn,'within',[T.sessN T.indx],{'session','model'},'subset',T.roi==r);
            fprintf('\n posterior probability in %s\n',regname{r});
            anovaMixed(T.postProp-0.5,T.sn,'within',[T.sessN T.indx],{'session','model'},'subset',T.roi==r);
        end
        
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-tests on the effect of model in sess %d in %s \n',ss,regname{r});
                fprintf('logBayes factor \n');
                ttestDirect(T.logBayes,[T.sn],2,'onesample','subset',T.roi==r&T.sessN==ss,'split',T.indx);
                fprintf('\n post-hoc t-tests on knockIN factors in sess %d in %s \n',ss,regname{r});
                fprintf('knockIN factor \n');
                ttestDirect(T.knockIN,[T.sn],2,'onesample','subset',T.roi==r&T.sessN==ss,'split',T.indx);
                fprintf('\n post-hoc t-tests on knockOUT factors in sess %d in %s \n',ss,regname{r});
                fprintf('knockOUT factor \n');
                ttestDirect(T.knockOUT,[T.sn],2,'onesample','subset',T.roi==r&T.sessN==ss,'split',T.indx);
                fprintf('\n post-hoc t-tests on posterior probability in sess %d in %s \n',ss,regname{r});
                fprintf('post prob \n');
                ttestDirect(T.postProp-0.5,[T.sn],2,'onesample','subset',T.roi==r&T.sessN==ss,'split',T.indx);
            end
        end    
    case 'PCM_constructModelFamily_seqType'
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        seqType='trained';  % trained or untrained
        parcelType='162tessels'; % Brodmann or 162tessels
        reg = [1,2,16,20,43,44,45,46,47,48,49,50,52,103,104,105]; % for tessels, 1:8 for Brodmann
        hemi=[1:2];
        sn=[4,5,7,8,9,11,13,14,15,16]; % change
        sessN=[1:4];
        AllReg=[];
        K=[];
        KK=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','seqType','parcelType','hemi'})
        tic;
        for ss = sessN
            B=load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for h = hemi
                for r = reg                    
                    for p=1:length(sn)
                        
                        glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                        
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.regType==r & B.regSide==h)}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                        end
                        
                        switch(seqType)
                            case 'trained'
                                st=1;
                            case 'untrained'
                                st=2;
                        end
                        
                        partVec{p} = D.run(D.seqType==st);  % runs/partitions
                        
                        m = pcm_defineSequenceModels_seqType(Seq,p,st);
                        [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = zeros(size(D.seqType));
                        indx(D.seqType==st)=1;
                        
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                        
                    end;
                    
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*r;
                    T.hemi = ones(size(T.SN))*h;
                    T.sessN = ones(size(T.SN))*ss;
                    AllReg=addstruct(AllReg,T);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*r;
                    K.hemi = ones(length(K.indx),1)*h;
                    K.sessN = ones(length(K.indx),1)*ss;
                    
                    KK=addstruct(KK,K);
                end
            end
        end
        
        keyboard;
        % save variables;
        dircheck(fullfile(pcmDir));
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_seqType_%s.mat',parcelType)),'Comb');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_%s_%s.mat',parcelType,seqType)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_%s.mat',parcelType,seqType)),'-struct','KK');
    
    case 'PCM_FingSeq_modelFamily'
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        parcelType='162tessels'; % Brodmann or 162tessels
        reg = [1,2,16,20,43,44,45,46,47,48,49,50,52,103,104,105]; % for tessels, 1:8 for Brodmann
        hemi=[1:2];
        sn=[4,5,7,8,9,11,13:18]; % change
        sessN=[1:4];
        AllReg=[];
        K=[];
        KK=[];
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi'})

        for ss = sessN
            B=load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for h=hemi
                for r = reg
                    for p=1:length(sn)
                        
                        glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});                        
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p)&B.regType==r&B.regSide==h)};
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p)&B.regType==r&B.regSide==h)};
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p)&B.regType==r)}; % no intercept - use T.betaRAWint otherwise
                        end
                        
                        partVec{p} = D.run;  % runs/partitions
                       % m = pcm_defineSequenceFingerModels(Seq,p);
                        m = pcm_defineSequenceFingerModels_new(Seq,p);
                        [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                       % [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1,'alwaysInclude',1);
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
      
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = ones(size(D.run));
                        if ~isnan(beta)
                            Data{p} = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                            fit(p)=1;
                        else
                            fit(p)=0;
                        end
                    end;
                    
                    if sum(fit)==length(sn)
                        % fit models
                        T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                        T.roi = ones(size(T.SN))*r;
                        T.hemi = ones(size(T.SN))*h;
                        T.sessN = ones(size(T.SN))*ss;
                        AllReg=addstruct(AllReg,T);
                        
                        % calculations - posterior probability, knockIN/OUT
                        K = pcm_calc(T.cross_likelihood,Comb);
                        K.roi = ones(length(K.indx),1)*r;
                        K.hemi = ones(length(K.indx),1)*h;
                        K.sessN = ones(length(K.indx),1)*ss;
                        KK=addstruct(KK,K);
                    end
                end
            end
        end
        
        % save variables;
        dircheck(fullfile(pcmDir));
        r1='reg'; r2='thetaCr'; AllReg=rmfield(AllReg,r1); AllReg=rmfield(AllReg,r2);
        %save(fullfile(pcmDir,sprintf('ModelFamilyComb_SeqFing_%s.mat',parcelType)),'Comb');
        %save(fullfile(pcmDir,sprintf('ModelFamily_Fit_SeqFing_%s.mat',parcelType)),'-struct','AllReg');
        %save(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_%s.mat',parcelType)),'-struct','KK');
        % ALL PARCELS! - added ALL
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_SeqFing_NEW_ALL_%s.mat',parcelType)),'Comb');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_SeqFing_NEW_ALL_%s.mat',parcelType)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_NEW_ALL_%s.mat',parcelType)),'-struct','KK');
       % name = {'ST','SS','FF','FA'};  
        name = {'TS','US','FF','FA'};   
    case 'PLOT_PCM_SeqFing_SURFACE'
        parcelType='162tessels';
        sessN=[1:4];
        thres=[0 1]; 

        vararginoptions(varargin,{'parcelType','sessN','thres','var'});
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_NEW_ALL_%s.mat',parcelType)));
        
        tesselNum=unique(A.roi);
        
        for h=1:2;    % only contralateral
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            for ss=sessN
                RGBdata=zeros(size(C.data,1),3);
                for r=tesselNum'
                    S=getrow(A,A.roi==r & A.hemi==h & A.sessN==ss);
                    % data from that hemi, tessel, session
                    
                    % Evidence for sequence-specific model
                   % E.seq=S.logBayes(S.indx==2);
                    E.seq=S.logBayes(S.indx==1);
                    % Evidence for finger model - take first or all fingers
                    if mean(S.logBayes(S.indx==3))-mean(S.logBayes(S.indx==4))<0
                        E.fing=S.logBayes(S.indx==4);
                    else
                        E.fing=S.logBayes(S.indx==3);
                    end
                     % indicate the right tessel
                    indx=C.index(C.data==r);   
                    
                    RGBdata(indx,1)=mean(E.seq-E.fing);
                    RGBdata(indx,2)=mean(E.fing-E.seq);

                end
                RGBdata(RGBdata(:,1)<thres(1),1)=0;
                RGBdata(RGBdata(:,2)<thres(1),2)=0;
                RGBdata(RGBdata(:,3)<thres(1),3)=0;
                
                scale=[thres;thres;thres];
                name={sprintf('sess%d_PCM_SeqOverFing',ss)};
                %name={sprintf('logBayes_%s_sess%d',featureName{f},ss)};
                R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.PCM_FingSeq_wholeBrain_sess%d.RGB_paint',hem{h},ss)),R);
                %caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.logBayes_%s_sess%d.RGB_paint',hem{h},featureName{f},ss)),R);
            end
        end
    case 'PLOT_PCM_SeqFing_ROI'
        parcelType='Brodmann';
        sessN=[1:4];
        roi=[1:5,7,8];
        
        vararginoptions(varargin,{'parcelType','sessN','roi'});
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_NEW_%s.mat',parcelType)));
        
        indx=1;
        for r = roi
            figure(1)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.logBayes,'subset',A.roi==r,'split',A.indx,'leg',{'Trained Seq','Untrained Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('logBayes factor');
            else
                ylabel('');
            end

            figure(2)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.knockIN,'subset',A.roi==r,'split',A.indx,'leg',{'Trained Seq','Untrained Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockIN factor');
            else
                ylabel('');
            end
            
            figure(3)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.knockOUT,'subset',A.roi==r,'split',A.indx,'leg',{'Trained Seq','Untrained Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockOUT factor');
            else
                ylabel('');
            end
            
            figure(4)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.postProp,'subset',A.roi==r,'split',A.indx,'leg',{'Trained Seq','Untrained Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('Posterior probability');
            else
                ylabel('');
            end
            
            
           indx=indx+1;
        end
    case 'PLOT_PCM_SeqFing_ROI_reduced'
        parcelType='Brodmann';
        sessN=[1:4];
        roi=[1:5,7,8];
        
        vararginoptions(varargin,{'parcelType','sessN','roi'});
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_NEW_%s.mat',parcelType)));
        
        FingSeqStyle=style.custom({'red','green'},'markersize',12);
        
        indx=1;
        for r = roi
            figure(1)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.logBayes,'subset',A.roi==r&(A.indx==1 | A.indx==4),'split',A.indx,'leg',{'Sequence','Fingers'},'leglocation','northeast','style',FingSeqStyle);
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('logBayes factor');
            else
                ylabel('');
            end

            figure(2)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.knockIN,'subset',A.roi==r&(A.indx==1 | A.indx==4),'split',A.indx,'leg',{'Sequence','Fingers'},'leglocation','northeast','style',FingSeqStyle);
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockIN factor');
            else
                ylabel('');
            end
            
            figure(3)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.knockOUT,'subset',A.roi==r&(A.indx==1 | A.indx==4),'split',A.indx,'leg',{'Sequence','Fingers'},'leglocation','northeast','style',FingSeqStyle);
            plt.match('y');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockOUT factor');
            else
                ylabel('');
            end
            
            figure(4)
            subplot(1,numel(roi),indx);
            plt.bar([A.sessN>3 A.sessN],A.postProp,'subset',A.roi==r&(A.indx==1 | A.indx==4),'split',A.indx,'leg',{'Sequence','Fingers'},'leglocation','northeast','style',FingSeqStyle);
            plt.match('y');
            drawline(0.5,'dir','horz');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('Posterior probability');
            else
                ylabel('');
            end
            
            
           indx=indx+1;
        end
    case 'STATS_PCM_FingSeq'
        parcelType='Brodmann';
        sessN=[1:4];
        roi=[1:5,7,8];
        modelName={'Seq','Fing'};
        modelIndx=[1,4];
        
        vararginoptions(varargin,{'parcelType','sessN','roi'});
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_NEW_%s.mat',parcelType)));
        
        T=getrow(A,A.indx==1|A.indx==4); % only sequence vs. allFing models
        
        for r=roi
            fprintf('\n logBayes factor in %s\n',regname{r});
            anovaMixed(T.logBayes,T.sn,'within',[T.sessN T.indx],{'session','model'},'subset',T.roi==r);
        end
        
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-tests on the fing vs. seq model in %s session %d \n',regname{r},ss);
                ttestDirect(T.logBayes,[T.indx T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);
            end
        end
        
        for r=roi
            for sTr=1:3
                for m=1:2
                fprintf('\n post-hoc t-tests on sess transition %d-%d on %s model in %s \n',sTr,sTr+1,modelName{m},regname{r});
                ttestDirect(T.logBayes,[T.sessN T.sn],2,'paired','subset',T.roi==r&T.indx==modelIndx(m)&T.sessN==sTr|T.sessN==sTr+1);
                end
            end
        end
        
        keyboard;
               
    case 'PLOT_PCM_allModelsBayes'
        reg=1;
        parcelType='Brodmann'; % Brodmann or 162tessels
        sessN=[1:4];
        vararginoptions(varargin,{'reg','parcelType'});
        load(fullfile(pcmDir,'ModelFamilyComb.mat')); % loads combination of models as Comb
        name = {'ST','TS','US','FF','FA','FT'};  
        
        TT=load(fullfile(pcmDir,sprintf('PCM_allModelsFamily_%s.mat',parcelType)));
        
        for r=reg;
            figure(r)
            for s=sessN
                b=getrow(TT,TT.roi==r&TT.sessN==s);
                subplot(4,1,s)
                plt.bar(b.modelIndx,b.bayesEst);
                drawline(0,'dir','horz');
                plt.match('y');
                if s==1
                    ylabel('Log Bayes - model');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
        end        
    case 'PLOT_PCM_modelFamily_seqType'
        %plots knock-IN / OUT / posterior probability
        reg=1;
        parcelType='Brodmann'; % Brodmann or 162tessels
        vararginoptions(varargin,{'reg','parcelType'});
        
        KK=[];
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_trained.mat',parcelType)));
        B=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_untrained.mat',parcelType)));
        KK=addstruct(KK,A); KK=addstruct(KK,B);
        KK.seqType=[ones(size(A.sn));ones(size(A.sn))*2];
        
        name = {'Seq','FF','FA','FT'};
        fig_num=0;
        for r = reg
            % knock-in
            figure(fig_num+1);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.knockIN,'split',a.seqType,'style',stySeq,'leg',{'trained','untrained'},'leglocation','north');
                
                if s==1
                    ylabel('Log Bayes - knock IN');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
                set(gca,'XTick',[1.5 4.5 8 11.5],'XTickLabel',name);
                drawline(0,'dir','horz');
                plt.match('y');
            end
            % knock-out
            figure(fig_num+2);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.knockOUT,'split',a.seqType,'style',stySeq,'leg',{'trained','untrained'},'leglocation','north');
                if s==1
                    ylabel('Log Bayes - knock OUT');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
                set(gca,'XTick',[1.5 4.5 8 11.5],'XTickLabel',name);
                drawline(0,'dir','horz');
                plt.match('y');
            end
            % posterior probability
            figure(fig_num+3);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.postProp,'split',a.seqType,'style',stySeq,'leg',{'trained','untrained'},'leglocation','north');             
                set(gca,'XTick',[1.5 4.5 8 11.5],'XTickLabel',name);
                drawline(0.5,'dir','horz');
                plt.match('y');
                if s==1
                    ylabel('Posterior probability');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
            % posterior Log
            figure(fig_num+4);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.logBayes,'split',a.seqType,'style',stySeq,'leg',{'trained','untrained'},'leglocation','north');
                set(gca,'XTick',[1.5 4.5 8 11.5],'XTickLabel',name);
                drawline(0,'dir','horz');
                plt.match('y');
                if s==1
                    ylabel('Log Bayes - post prob');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
        end    
        
        % save combined structure of trained / untrained
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_bothSeqType.mat',parcelType)),'-struct','KK');
    case 'PLOT_SURFACE_PCM_tessels_seqType'
        parcelType='162tessels';
        sessN=[1:4];
        thres=[0 1]; % [0 10] for knockIN [0 1 for knockOUT]
        featureSet=[1:4];
        featureName={'Seq','FirstFing','AllFing','FingTrans'};
        vararginoptions(varargin,{'parcelType','sessN','thres'});
        KK=[];
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_trained.mat',parcelType)));
        B=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_untrained.mat',parcelType)));
        KK=addstruct(KK,A); KK=addstruct(KK,B);
        KK.seqType=[ones(size(A.sn));ones(size(A.sn))*2]; 
        
        tesselNum=unique(KK.roi);

        for h=1:2;    % only contralateral
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            for ss=sessN
                for f=featureSet
                    RGBdata=zeros(size(C.data,1),3);
                    for r=tesselNum'
                        S=getrow(KK,KK.roi==r & KK.hemi==h & KK.sessN==ss & KK.indx==f);
                        % data from that hemi, tessel, session - per feature
                        
                        indx=C.index(C.data==r);    % indicate the right tessel
                        
                        RGBdata(indx,1)=mean(S.knockOUT(S.seqType==1)*(-1)); % trained
                        RGBdata(indx,3)=mean(S.knockOUT(S.seqType==2)*(-1)); % untrained
                        RGBdata(indx,2)=mean(S.knockOUT(S.seqType==1)*(-1))-mean(S.knockOUT(S.seqType==2)*(-1)); % difference between the two
                    end
                    RGBdata(RGBdata(:,1)<thres(1))=0;
                    RGBdata(RGBdata(:,3)<thres(1))=0;
                    
                    scale=[thres;thres;thres];
                    name={sprintf('knockOUT_%s_sess%d',featureName{f},ss)};
                    %name={sprintf('logBayes_%s_sess%d',featureName{f},ss)};
                    R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                    caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.knockOUT_%s_sess%d.RGB_paint',hem{h},featureName{f},ss)),R);
                    %caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.logBayes_%s_sess%d.RGB_paint',hem{h},featureName{f},ss)),R);
                end
            end
        end
    case 'PLOT_SURFACE_PCM_tessels'     
        parcelType='162tessels';
        sessN=[1:4];
        thres=[0 3]; % [0 10] for knockIN [0 1 for knockOUT]
        var='knockIN'; %knockIN / knockOUT / logBayes
        %featureSet=[1,2,4:6]; % skip untrained - taking trained / untrained together
        featureSet=[1,3,4,5];
        %featureName={'SeqType','SeqSpec','','FirstFing','AllFing','FingTrans'};
         featureName={'SeqSpec','','FirstFing','AllFing','FingTrans'};
        % colours: green - red - blue - red - blue - green
        vararginoptions(varargin,{'parcelType','sessN','thres','var'});
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_NEW_%s.mat',parcelType)));
        
        tesselNum=unique(A.roi);

        for h=1:2;    % only contralateral
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            for ss=sessN
                for f=featureSet
                    RGBdata=zeros(size(C.data,1),3);
                    for r=tesselNum'
                        S=getrow(A,A.roi==r & A.hemi==h & A.sessN==ss & A.indx==f);
                        % data from that hemi, tessel, session - per feature
                        
                        indx=C.index(C.data==r);    % indicate the right tessel
                        switch (var)
                            case 'knockIN'
                                I=S.knockIN;
                                if f==1
                                    T=getrow(A,A.roi==r & A.hemi==h & A.sessN==ss & A.indx==f+1);  % untrained
                                    TI=T.knockIN;
                                end
                            case 'knockOUT'
                                I=S.knockOUT*(-1);
                                if f==1
                                    T=getrow(A,A.roi==r & A.hemi==h & A.sessN==ss & A.indx==f+1);  % untrained
                                    TI=T.knockOUT*(-1);
                                end
                            case 'logBayes'
                                I=S.logBayes;
                                if f==1
                                    T=getrow(A,A.roi==r & A.hemi==h & A.sessN==ss & A.indx==f+1);  % untrained
                                    TI=T.logBayes;
                                end
                        end
                        
                        if f==1
                        %if f==2 % red
                            RGBdata(indx,1)=mean(I); % trained
                            RGBdata(indx,3)=mean(TI); % untrained
                        else % always green
                            RGBdata(indx,2)=mean(I);
                        end
                    end
                    RGBdata(RGBdata(:,1)<thres(1))=0;
                    RGBdata(RGBdata(:,2)<thres(1))=0;
                    RGBdata(RGBdata(:,3)<thres(1))=0;
                    
                    scale=[thres;thres;thres];
                    name={sprintf('%s_%s_sess%d_PCMfull',var,featureName{f},ss)};
                    %name={sprintf('logBayes_%s_sess%d',featureName{f},ss)};
                    R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                    caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.PCM_fullNEW_%s_%s_sess%d.RGB_paint',hem{h},var,featureName{f},ss)),R);
                    %caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.logBayes_%s_sess%d.RGB_paint',hem{h},featureName{f},ss)),R);
                end
            end
        end
        
    case 'STATS_PCM_modelFamily_seqType'
        roi=[1:8];
        sessN=[1:4];
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','sessN','parcelType'});
       
        T=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_bothSeqType.mat',parcelType)));
        
        modelName = {'Seq','FF','FA','FT'};
        
        for r=roi
            fprintf('\n logBayes factor in %s\n',regname{r});
            anovaMixed(T.logBayes,T.sn,'within',[T.sessN T.seqType T.indx],{'session','seqType','model'},'subset',T.roi==r);
            fprintf('\n knockIN factor in %s\n',regname{r});
            anovaMixed(T.knockIN,T.sn,'within',[T.sessN T.seqType T.indx],{'session','seqType','model'},'subset',T.roi==r);
            fprintf('\n knockOUT factor in %s\n',regname{r});
            anovaMixed(T.knockOUT,T.sn,'within',[T.sessN T.seqType T.indx],{'session','seqType','model'},'subset',T.roi==r);
            fprintf('\n posterior probability in %s\n',regname{r});
            anovaMixed(T.postProp-0.5,T.sn,'within',[T.sessN T.seqType T.indx],{'session','seqType','model'},'subset',T.roi==r);
        end
        
        for r=roi
            for ss=sessN
                for m=1:size(modelName,2)
                    fprintf('\n post-hoc t-tests on the effect of seqType in model %s sess %d in %s \n',modelName{m},ss,regname{r});
                    fprintf('logBayes factor \n');
                    ttestDirect(T.logBayes,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss&T.indx==m);
                    fprintf('\n post-hoc t-tests on knockIN factors in model %s sess %d in %s \n',modelName{m},ss,regname{r});
                    fprintf('knockIN factor \n');
                    ttestDirect(T.knockIN,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss&T.indx==m);
                    fprintf('\n post-hoc t-tests on knockOUT factors in model %s sess %d in %s \n',modelName{m},ss,regname{r});
                    fprintf('knockOUT factor \n');
                    ttestDirect(T.knockOUT,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss&T.indx==m);
                    fprintf('\n post-hoc t-tests on posterior probability in model %s sess %d in %s \n',modelName{m},ss,regname{r});
                    fprintf('post prob \n');
                    ttestDirect(T.postProp-0.5,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss&T.indx==m);
                end
            end
        end
   
    case 'CONN_timeseries'
        % univariate connectivity - of mean timeseries across regions
        sessN=[1:4];
        sn=[4:9,11:25];
        vararginoptions(varargin,{'reg_cortex','reg_BG','sessN','sn','hemi'});
        
        CC = [];
        for ss=sessN
            for s=sn
                % load SPM and all regions
                load(fullfile(glmSessDir{ss},subj_name{s},'SPM.mat'));
                SPM=spmj_move_rawdata(SPM,fullfile(imagingDir, subj_name{s}));
                
                load(fullfile(regDir,sprintf('%s_Brodmann_regions.mat',subj_name{s})));
                Rc = R;
                load(fullfile(regDir,sprintf('%s_BG-striatum_regions.mat',subj_name{s})));
                Rs = R;
                R       = [Rc Rs];              
                [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
                C.corr_adj = rsa_vectorizeRDM(corrcoef(y_adj));
                C.corr_hat = rsa_vectorizeRDM(corrcoef(y_hat));
                C.reg = reg;
                C.regType = regType;
                C.regSide = regSide;
                C.sn = s;
                C.sessN = ss;
                
                CC = addstruct(CC,C);
                fprintf('Done sess-%d %s\n',ss,subj_name{s});
            end
            fprintf('Done sess-%d all subjects.\n\n\n',ss);
        end
        % save structure
        save(fullfile(distPscDir,'connect_timeseries'),'-struct','CC');  
    case 'CONN_RDM'
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
            C = load(fullfile(regDir,sprintf('stats_Brodmann_%s_sess%d',betaChoice,ss)));
            B = load(fullfile(BGDir,'stats',sprintf('stats_%s_%s_sess%d',parcelBG,betaChoice,ss)));
            for st=1:2 % seqType
                for s=sn
                    % get rdms from all regions
                    for r=1:numReg
                        if regType(r) == 1
                            T = getrow(C,C.SN==s & C.region==regAll(r));
                        else
                            T = getrow(B,B.SN==s & B.region==regAll(r));
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
                    K.sessN = ss;
                    KK=addstruct(KK,K);
                    fprintf('Done %s sess-%d\n',subj_name{s},ss);
                end; % sn
            end; % seqType
        end; % session
        save(fullfile(distPscDir,sprintf('connect_RDM_%s',betaChoice)),'-struct','KK');
    case 'CONN_perHemi'
        connType = 'RDM_multiPW';
        var      = {'rdmCorr','distCorr'}; 
        % var for RDM_multiPW / RDM_uniPW: {'rdmCorr','distCorr'}
        % var for timeseries: {'corr_adj','corr_hat'}
        vararginoptions(varargin,{'connType','var'});
        T = load(fullfile(distPscDir,sprintf('connect_%s',connType)));
        
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
            save(fullfile(distPscDir,sprintf('connect_%s_%s',connType,hem{h})),'-struct','H');
        end
        
    case 'CONN_RDM_plot'
        regName = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp','V12','Caudate','Putamen'};
        seqType = [1 2];
        sessN = [1:4];
        hemi = 'lh';
        betaChoice='multiPW';
        seqLabel = {'trained','untrained'};
        vararginoptions(varargin,{'seqType','sessN','betaChoice','hemi'});
        T = load(fullfile(distPscDir,sprintf('connect_RDM_%s_%s',betaChoice,hemi)));
        
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
                    K = getrow(T,T.seqType==st & T.sessN==ss);
                    subplot(2,numel(sessN),indx);
                    imagesc(rsa_squareRDM(mean(K.(metric),1)));
                    caxis([0 0.8]);
                    set(gca,'XTickLabel',regName,'YTickLabel',regName);
                    title(sprintf('sess-%d %s %s',ss,seqLabel{st},metric));
                    indx=indx+1;
                end
            end
        end
    case 'CONN_timeseries_plot'
        regName = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp','V12','Caudate','Putamen'};
        sessN = [1:4];
        hemi = 'lh';
        vararginoptions(varargin,{'seqType','sessN','betaChoice','hemi'});
        T = load(fullfile(distPscDir,sprintf('connect_timeseries_%s',hemi)));
        
        figure
        indx=1;
        for dt=1:2
            if dt==1
                metric='corr_adj';
            else
                metric='corr_hat';
            end
            for ss=sessN
                K = getrow(T,T.sessN==ss);
                subplot(2,numel(sessN),indx);
                imagesc(rsa_squareRDM(mean(K.(metric),1)));
                caxis([0 1]);
                set(gca,'XTickLabel',regName,'YTickLabel',regName);
                title(sprintf('sess-%d %s',ss,metric));
                indx=indx+1;
            end
        end
    case 'CONN_RDM_vs_timeseries'
        regName = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp','Caudate','Putamen'};
        seqType = [1 2];
        sessN = [1:4];
        hemi = 'lh';
        betaChoice='uniPW';
        seqLabel = {'trained','untrained'};
        vararginoptions(varargin,{'seqType','sessN','betaChoice','hemi'});
        
        T1 = load(fullfile(distPscDir,sprintf('connect_RDM_%s_%s',betaChoice,hemi)));
        T2 = load(fullfile(distPscDir,sprintf('connect_timeseries_%s',hemi)));
        sType = [ones(size(T2.sn));ones(size(T2.sn))*2];
        T3 = [];
        T3 = addstruct(T3,T2);
        T3 = addstruct(T3,T2);
        T3.seqType = sType;
        keyboard;
        
        figure
        indx=1;
        for ss=sessN
            K1 = getrow(T1,T1.sessN==ss);
            K2 = getrow(T3,T3.sessN==ss);
            subplot(1,numel(sessN),indx);
            plt.scatter(K2.corr_hat(:),K1.rdmCorr(:),'split',K2.seqType);
            title(sprintf('sess-%d',ss));
            indx=indx+1;
            plt.match('y');
        end

        
    case 'CONN_RDM_subj'
        sn=[4:9,11:25];
        sessN=[1:4];
        reg_cortex=[1:5,7,8];
        reg_BG=[1,2];
        parcelBG='BG-striatum';
        betaChoice='multiPW';
        vararginoptions(varargin,{'sn','sessN','reg','parcelBG','betaChoice'});
        
        
        regAll  = [reg_cortex reg_BG];
        numReg  = length(regAll);     
        regType = [ones(size(reg_cortex)) ones(size(reg_BG))*2];
        
        KK=[];
        for ss=sessN
            C = load(fullfile(regDir,sprintf('stats_Brodmann_%s_sess%d',betaChoice,ss)));
            B = load(fullfile(BGDir,'stats',sprintf('stats_%s_%s_sess%d',parcelBG,betaChoice,ss)));
            for r=1:numReg
                if regType(r) == 1
                    R = getrow(C,C.region==regAll(r));
                else
                    R = getrow(B,B.region==regAll(r));
                end
                for st=1:2 % seqType
                    for g=1:2 % group
                        uniqueSN = sn(mod(sn,2)==g-1);
                        % get rdms from all subjects
                        for s=1:length(uniqueSN)
                            T = getrow(R,R.SN==uniqueSN(s));
                            tmp = rsa_squareRDM(T.RDM);
                            % get the right entries for trained / untrained
                            if st==1
                                RDM{s} = tmp(1:6,1:6);
                            else
                                RDM{s} = tmp(7:12,7:12);
                            end
                        end
                        % perform calculations
                        K.rdmCorr       = {rsa_calcCorrRDMs(RDM)};
                        K.distCorr      = {rsa_calcDistCorrRDMs(RDM)};
                        K.mean_RDMCorr  = mean(rsa_calcCorrRDMs(RDM));
                        K.mean_distCorr = mean(rsa_calcDistCorrRDMs(RDM));
                        K.seqType       = st;
                        K.reg           = r;
                        K.regLoc        = regType(r);
                        K.regType       = regAll(r);
                        K.sessN         = ss;
                        K.group         = g;
                        KK = addstruct(KK,K);
                        fprintf('Done sess-%d group-%d: %s\n',ss, g,regname{r});
                    end; % group
                end; % seqType
            end; % reg
        end; % session
        save(fullfile(distPscDir,sprintf('RDM_connect_acrSubj_%s',betaChoice)),'-struct','KK');   
    case 'CONN_RDM_subj_plot'
        regName = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp','CaudateN','Putamen'};
        reg=3;
        seqType = [1 2];
        sessN=[1:4];
        betaChoice='multiPW';
        metric='rdmCorr';
        vararginoptions(varargin,{'seqType','sessN','betaChoice','reg','metric'});
        T = load(fullfile(distPscDir,sprintf('RDM_connect_acrSubj_%s',betaChoice)));
        
        for r=reg
            for g=1:2 % group   
                figure
                indx=1;
                for st=seqType
                    for ss=sessN
                        K = getrow(T,T.seqType==st & T.sessN==ss & T.reg==r & T.group==g);
                        subplot(numel(seqType),numel(sessN),indx);
                        imagesc(rsa_squareRDM(K.(metric){:}));
                        caxis([0 0.5]);
                        title(sprintf('sess-%d %s, G%d',ss,regName{r},g));
                        indx=indx+1;
                    end
                end
            end
        end
    case 'CONN_RDM_subj_lineplot'
        regName = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp','CaudateN','Putamen'};
        reg=3;
        betaChoice='multiPW';
        vararginoptions(varargin,{'betaChoice','reg'});
        T = load(fullfile(distPscDir,sprintf('RDM_connect_acrSubj_%s',betaChoice)));
        
        figure
        indx=1;
        for r=reg
            subplot(1,length(reg),indx)
            plt.line([T.sessN>3 T.sessN],T.mean_RDMCorr,'split',T.seqType,'subset',T.reg==r,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
            drawline(0,'dir','horz');
            title(regName{r});
            if r==1
                xlabel('Session');
                ylabel('Across subj correlation');
            else
                ylabel('');
            end
            plt.match('y');
            indx=indx+1;
        end
        
    case 'SEARCH_check'
        sn=[4:9,11:25];
        sessN=[1:4];
        vararginoptions(varargin,{'sn','sessN'});
        count=0;
        for ss=sessN
            for s=sn
                dir = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_LDC.nii',subj_name{s},ss));
                if exist(dir)
                    fprintf('Searchlight exists: \t%s-sess%d\n',subj_name{s},ss);
                else
                    fprintf('Searchlight missing: \t%s-sess%d\t!!!\n',subj_name{s},ss);
                    count=count+1;
                end
            end
        end
        fprintf('\n\nAltogether missing %d searchlights\n',count);
    
    case 'CLUSTER_choose'
        sessN=1;
        betaChoice='multi';
        vararginoptions(varargin,{'sessN','betaChoice'});
        
        for ss=sessN       
            indx=[];
            T=load(fullfile(betaDir,'group',sprintf('stats_162tessels_%sPW_sess%d',betaChoice,ss)));
            for r=unique(T.region)'
                T1 = getrow(T,T.region==r);
                if ~any(isnan(T1.dist_train)) % make sure all subjects have data there
                    [t,p]=ttestDirect(T.dist_train,[T.SN],1,'onesample','subset',T.region==r);
                    if (p<0.01 & t>0) % include into selection
                        indx=[indx;r];
                    end
                end
            end
            if ss==sessN(1)
                choice=indx;
            else
                choice=unique([indx; choice]);
            end
        end
        varargout={choice}; 
    case 'CLUSTER_extractROI_search'
        sn=[4:9,11:25];
        sessN=[1:4];
        parcelType='Brodmann'; % Brodmann or 162tessels
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        TT=[];
        for ss=sessN
            fprintf('\n\nSession %d\n',ss);
            for s=sn
                fprintf('\nSubject: %d\n',s) % output to user
                load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)]));  % load subject's region parcellation (R)
                roi=1:size(R,2);
                cd(fullfile(glmSessDir{ss},subj_name{s}));
                file=sprintf('%s_sess%d_LDC.nii',subj_name{s},ss);
                V=spm_vol(file);
                for r=roi
                    data = region_getdata(V,R{r});
                    % calculate the mean of the region
                    Mdata           = nanmean(data,2);
                    RDM             = rsa_squareRDM(Mdata');
                    train           = RDM([1:6],[1:6]);
                    untrain         = RDM([7:12],[7:12]);
                    T.RDM_train     = rsa_vectorizeRDM(train);
                    T.RDM_untrain   = rsa_vectorizeRDM(untrain);
                    T.dist_mean_train   = nanmean(T.RDM_train);
                    T.dist_mean_untrain = nanmean(T.RDM_untrain);
                    T.sn = s;
                    T.sessN = ss;
                    T.roi = r;
                    if r<(numel(roi)/2+1)
                        T.regType = r;
                        T.regSide = 1;
                    else
                        T.regType = r/2;
                        T.regSide = 2;
                    end
                    fprintf('%d.',r);
                    TT=addstruct(TT,T);
                end
            end
        end
        keyboard;
        save(fullfile(distPscDir,sprintf('cluster_RDM_search_%s',parcelType)),'-struct','TT');
    case 'CLUSTER_extractROI_dist'
        sn=[4:9,11:28,30];
        sessN=[1:4];
        parcelType='Brodmann'; % Brodmann or 162tessels
        betaChoice='multi';
        sessType = 'within'; % within one session only / across sessions
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        TT=[];
        for ss=sessN
            fprintf('\n\nSession %d\n',ss);
            V=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,ss)));
            if strcmp(parcelType,'Brodmann')
                roi     = [1:16];
                regSide = [ones(1,8) ones(1,8)*2];
                regType = [1:8 1:8];
            elseif strcmp(parcelType,'162tessels')
                roi     = sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                regSide = ones(size(roi));
                regSide(roi>162)=2;
                regType = roi;
                regType(regSide==2)=regType(regSide==2)-162;
            end
            for s=sn
                fprintf('\nSubject: %d\n',s) % output to user
                load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)]));  % load subject's region parcellation (R)             
                for r=1:length(roi)
                    D=getrow(V,V.region==roi(r) & V.SN==s);
                    T.RDM_train     = D.RDM_train;
                    T.RDM_untrain   = D.RDM_untrain;
                    T.RDM_all       = D.RDM;
                    T.dist_mean_train   = nanmean(T.RDM_train);
                    T.dist_mean_untrain = nanmean(T.RDM_untrain);
                    T.dist_all          = nanmean(T.RDM_all);
                    T.sn = s;
                    T.sessN = ss;
                    T.roi = roi(r);
                    T.regType = regType(r);
                    T.regSide = regSide(r);
                    fprintf('%d.',r);
                    TT=addstruct(TT,T);
                end
            end
        end
        if strcmp(sessType,'within');
            save(fullfile(distPscDir,sprintf('cluster_RDM_dist_%sSess%d_%s',sessType,sessN,parcelType)),'-struct','TT');
        else
            save(fullfile(distPscDir,sprintf('cluster_RDM_dist_%sSess%d-%d_%s',sessType,sessN(1),sessN(end),parcelType)),'-struct','TT');
        end
    case 'CLUSTER_extractROI_crossvalBetas'
        % crossvalidated version - separately for even and odd runs
        sn=[4:9,11:28];
        sessN=[1:4];
        parcelType='162tessels'; % Brodmann or 162tessels
        betaChoice='multi';
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        
        TT=[];
        nRun  = 8;
        nCond = 12;
        pVec  = kron([1:nRun]',ones(nCond,1));
        idxA = mod(pVec,2)==1; % split even and odd runs
        idxB = mod(pVec,2)==0;
        cVec  = kron(ones(nRun/2,1),[1:nCond]');   % for each partition
        
        if strcmp(parcelType,'Brodmann')
            roi     = [1:16];
            regSide = [ones(1,8) ones(1,8)*2];
            regType = [1:8 1:8];
        elseif strcmp(parcelType,'162tessels')
            roi     = sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
            % choose only clusters where group dist >0
            regSide = ones(size(roi));
            regSide(roi>162)=2;
            regType = roi;
            regType(regSide==2)=regType(regSide==2)-162;
        end

        for ss=sessN
            fprintf('\n\nSession %d\n',ss);
            V=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss)));
            for s=sn
                fprintf('\nSubject: %d\n',s) % output to user
                for r=1:length(roi)
                    beta=[];
                    D=getrow(V,V.region==roi(r) & V.SN==s);
                    switch betaChoice
                        case 'multi'
                            beta = D.betaW{:};
                        case 'uni'
                            beta = D.betaUW{:};
                        case 'raw'
                            beta = D.betaRAW{:};
                    end
                    if isnan(beta)
                        T.partA_train     = ones(1,6*5/2)*NaN;
                        T.partA_untrain   = ones(1,6*5/2)*NaN;
                        T.partA_all       = ones(1,12*11/2)*NaN;
                        T.partB_train     = ones(1,6*5/2)*NaN;
                        T.partB_untrain   = ones(1,6*5/2)*NaN;
                        T.partB_all       = ones(1,12*11/2)*NaN;
                    else
                        distA = rsa_squareRDM(rsa.distanceLDC(beta(idxA,:),pVec(idxA),cVec));
                        distB = rsa_squareRDM(rsa.distanceLDC(beta(idxB,:),pVec(idxB),cVec));
                        T.partA_train     = rsa_vectorizeRDM(distA(1:6,1:6));
                        T.partA_untrain   = rsa_vectorizeRDM(distA(7:12,7:12));
                        T.partA_all       = rsa_vectorizeRDM(distA);
                        T.partB_train     = rsa_vectorizeRDM(distB(1:6,1:6));
                        T.partB_untrain   = rsa_vectorizeRDM(distB(7:12,7:12));
                        T.partB_all       = rsa_vectorizeRDM(distB);
                    end   
                    T.sn = s;
                    T.sessN = ss;
                    T.roi = roi(r);
                    T.regType = regType(r);
                    T.regSide = regSide(r);
                    fprintf('%d.',r);
                    TT=addstruct(TT,T);
                    fprintf('Done %d/%d\n',find(s==sn),numel(sn));
                end
            end
        end
        save(fullfile(betaDir,'group',sprintf('betas_partition_%s',parcelType)),'-struct','TT');
    case 'CLUSTER_regress'
        sn=[4:9,11:28];
        sessN=[1:4];
        parcelType='162tessels'; % Brodmann or 162tessels
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        
        T = load(fullfile(betaDir,'group',sprintf('betas_partition_%s',parcelType)));
        RR=[];
        roi=unique(T.roi);
        for s=1:numel(sn)
            for i=1:numel(roi)*numel(sessN)
                for j=i:numel(roi)*numel(sessN)
                    D = getrow(T,T.sn==sn(s));
                    D1 = getrow(D,i);
                    D2 = getrow(D,j);
                    
                    beta = pinv(D1.partA_train')*D2.partA_train';
                    res = D2.partB_train'-(D1.partB_train'*beta);
                    SSR = sum(res.^2);
                    SST = sum(D2.partB_train.^2);
                    R2 = 1 - SSR/SST;
                    R.R2 = R2;
                    r = corrcoef(D2.partB_train',D1.partB_train'*beta);
                    R.r = r(1,2);
                    R.beta = beta;
                    R.reg1 = i;
                    R.reg2 = j;
                    R.sn   = sn(s);
                    R.s_idx = s;
                    RR=addstruct(RR,R);
                    % save struct
                    % save matrix
                    M_R2(j,i,s)=R2;
                    M_r(j,i,s)=R.r;
                end
                  fprintf('Done %s\treg: %d/%d\n',subj_name{sn(s)},i,numel(roi)*numel(sessN));
            end
        end
        keyboard;
    case 'CLUSTER_acrossROI'
        sn=[4:9,11:28,30];
        sessN=[1:4];
        parcelType='162tessels'; % Brodmann or 162tessels
        distType='euclidean'; % euclidean, correlation or cosine
        seqType='all'; % all, trained, untrained
        roiType='dist'; % search or dist
        sessType = 'within';
        maxClust=10;
        vararginoptions(varargin,{'sn','sessN','distType','roiType','parcelType','maxClust'});
        
        switch sessType
            case 'within'
                T=load(fullfile(distPscDir,sprintf('cluster_RDM_%s_%sSess%d_%s',roiType,sessType,sessN,parcelType)));
            case 'between'
                T=load(fullfile(distPscDir,sprintf('cluster_RDM_%s_%sSess%d-%d_%s',roiType,sessType,sessN(1),sessN(end),parcelType)));
        end
        roi=unique(T.roi)';
        S=[];
        for s=sn
            indx=1;
            for ss=sessN
                for r=roi
                    D=getrow(T,T.sn==s&T.sessN==ss&T.roi==r);
                    switch seqType
                        case 'trained'
                            r_dist(indx,:)=D.RDM_train;
                        case 'untrained'
                            r_dist(indx,:)=D.RDM_train;
                        case 'all'
                            r_dist(indx,:)=D.RDM_all;
                    end
                    indx=indx+1;
                end
            end
            % D=rsa_calcCorrRDMs(R_t);
            s = rsa_calcdist(r_dist,distType);
            S = cat(3,S,full(s));
        end
        W=nanmean(S,3);
        %  figure
        %  imagesc(W);
    %    eval=evalclusters(W,'kmeans','CalinskiHarabasz','KList',[1:10]);
    %    Y = kmeans(W,kNum);
        [C L U]=SpectralClustering(W,maxClust,3);
        % create a structure
        M.simGraph  = W;
        M.cluster   = C;
        M.laplace   = L;
        M.eigenv    = U;
        reg = [];
        ses = [];
        for i=sessN
            reg = [reg;roi'];
            ses = [ses;ones(size(roi'))*i];
        end
        M.roi=reg;
        M.sessN=ses;
        switch parcelType
            case '162tessels'
                M.regSide   = ones(size(M.roi));
                M.regSide(M.roi>162)=2;
                M.regType   = M.roi;
                M.regType(M.regSide==2)=M.regType(M.regSide==2)-162;
        end
        
        if strcmp(sessType,'within');
            save(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d',parcelType,distType,maxClust,seqType,sessType,sessN)),'-struct','M');

        else
            save(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d-%d',parcelType,distType,maxClust,seqType,sessType,sessN(1),sessN(end))),'-struct','M');
        end
    case 'CLUSTER_color'
        distType = 'distcorr';
        parcelType = '162tessels';
        maxClust = 10;
        sessType = 'within';
        sessN=[1:2];
        seqType='all';
        vararginoptions(varargin,{'distType','parcelType','maxClust','sessN','sessType','seqType'});
        
        if strcmp(sessType,'within');
            T = load(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d',parcelType,distType,maxClust,seqType,sessType,sessN)));
            
        else
            T = load(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d-%d',parcelType,distType,maxClust,seqType,sessType,sessN(1),sessN(end))));
        end
              
        % Evaluate similarity graph (W)
        [Csort,idx]=sort(T.cluster);
        W = full(T.simGraph);
        Wsort=W(idx,:);
        Wsort=Wsort(:,idx);
        Ncluster = numel(unique(T.cluster));
        
        % Laplacian eigenvector
        Lsort = T.laplace(idx,:); Lsort=Lsort(:,idx);
        Usort = T.eigenv(idx,:);
        
        % Merge similarities according to clusters
        Wsort(logical(eye(size(Wsort)))) = NaN;
        for i=1:Ncluster
            for j=1:Ncluster
                idxi=Csort==i;
                idxj=Csort==j;
                rWsrot(i,j) = nanmean(vec(Wsort(idxi,idxj)));
                rLsort(i,j) = nanmean(vec(Lsort(idxi,idxj)));
            end
            rUsort(i,:) = nanmean(Usort(idxi,:),1);
        end
        
        % inspect sorted diagram
        figure
        imagesc(Wsort);
        hold on;
        % get borders
        border = [0; find(diff(Csort))+1; length(Wsort)];
        for b=2:length(border)
            drawline(border(b-1),'dir','vert','lim',[border(b-1) border(b)],'color',[1 1 1]);
            drawline(border(b),'dir','vert','lim',[border(b-1) border(b)],'color',[1 1 1]);
            drawline(border(b-1),'dir','horz','lim',[border(b-1) border(b)],'color',[1 1 1]);
            drawline(border(b),'dir','horz','lim',[border(b-1) border(b)],'color',[1 1 1]);
        end
        
        % get dendogram linkage
        Z = linkage(real(rUsort),'ward','euclidean');
        Col = colorDendrogram(Z,size(rUsort,1),'colorspace','rgb','order',[1,2,3],'fig',0,'weight',1);

        figure
        [h,t,per] = dendrogram(Z,'Orientation','right');
        hold on;
        xlim=get(gca,'xlim');
        range=diff(xlim);
        for i=1:numel(per)
            plot(xlim(1),i,'o','markersize',8,...
                'markeredgecolor','k',...
                'markerfacecolor',Col(per(i),:));
            hold on;
        end
        
        for i=1:numel(h)
            set(h(i),'color',[0 0 0],'linewidth',1);
        end
        cd(distPscDir);
        dlmwrite(sprintf('my%dColors.txt',maxClust),Col);
    case 'CLUSTER_extractRDM'
        sessType='within';
        sessN=4;
        distType='cosine';
        clustN = 10;
        parcelType='162tessels';
        seqType='all';
        sn=[4:9,11:28,30];
        vararginoptions(varargin,{'sessType','sessN','distType','clustN'});
        
         switch sessType
            case 'within'
                R=load(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d',parcelType,distType,clustN,seqType,sessType,sessN)));
                T=load(fullfile(distPscDir,sprintf('cluster_RDM_dist_%sSess%d_%s',sessType,sessN,parcelType)));
             case 'across'
                R=load(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d-%d',parcelType,distType,clustN,seqType,sessType,sessN(1),sessN(end))));
                T=load(fullfile(distPscDir,sprintf('cluster_RDM_dist_%sSess%d-%d_%s',sessType,sessN(1),sessN(end),parcelType)));
         end
         DD=[];
         for ss=sessN
             % get group betas
             B = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss)));
             for s=sn
                 for c=1:clustN % extract data for each cluster
                     C = getrow(R,R.cluster==c);
                     T1 = getrow(T,ismember(T.roi,C.roi) & T.sn==s);
                     B1 = getrow(B,ismember(B.region,C.roi) & B.SN==s);
                     % concatenate betas into one region
                     beta=[];
                     for i=1:length(unique(C.roi))
                         if sum(sum(isnan(B1.betaW{i})))==0
                            beta = [beta B1.betaW{i}];
                         end
                     end
                     D.beta = {beta};
                     D.dist = nanmean(T1.RDM_all,1);
                     D.cluster = c;
                     D.sn = s;
                     DD=addstruct(DD,D);
                 end
                 fprintf('Done subject %d/%d\n',find(s==sn),length(sn));
             end
         end
         switch sessType
             case 'within'
                save(fullfile(distPscDir,sprintf('clusterData_%sDist-%dcluster_%sSess%d',distType,clustN,sessType,sessN)),'-struct','DD');
             case 'between'
                 save(fullfile(distPscDir,sprintf('clusterData_%sDist-%dcluster_%sSess_%d-%d',distType,clustN,sessType,sessN(1),sessN(end))),'-struct','DD');
         end
    case 'CLUSTER_plotRDM'
        sessType='within';
        sessN=4;
        distType='cosine';
        clustN = 10;
        parcelType='162tessels';
        vararginoptions(varargin,{'sessType','sessN','distType','clustN'});
        
        switch sessType
             case 'within'
                T=load(fullfile(distPscDir,sprintf('clusterData_%sDist-%dcluster_%sSess%d',distType,clustN,sessType,sessN)));
             case 'between'
                T=load(fullfile(distPscDir,sprintf('clusterData_%sDist-%dcluster_%sSess_%d-%d',distType,clustN,sessType,sessN(1),sessN(end))));
        end
         
        keyboard;
        figure
        for i=1:clustN
            subplot(clustN/5,5,i)
            T1=getrow(T,T.cluster==i);
            imagesc(rsa_squareRDM(nanmean(T1.dist,1)));
            colorbar;
            title(sprintf('RDM-cluster%d',i));
        end
        figure % subtract the mean
        tmp = mean(T.dist,1);
        for i=1:clustN
            subplot(clustN/5,5,i)
            T1=getrow(T,T.cluster==i);
            imagesc(rsa_squareRDM(nanmean(T1.dist,1)-tmp));
            colorbar;
            title(sprintf('RDM-cluster%d',i));
        end
        
    case 'CLUSTER_consist_reg'
          % estimate the consistency of RDMs in each tessel
        % both within and across subjects
        sn=[4:9,11:28,30];
        sessN=[1:4];
        hemi=[1,2]; % 1- contra, 2 - ipsi, [1,2] - both
        parcelType='162tessels'; % Brodmann or 162tessels; or combined
        sessType = 'across';
        seqType='all'; % all, trained, untrained
        maxClust=10;
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','hemi','maxClust','seqType'});
        
        T = load(fullfile(betaDir,'group',sprintf('betas_partition_%s',parcelType)));
        % split participants by groups
        if strcmp(sessType,'within')
            T = getrow(T,T.sessN==sessN);
        end
        if size(hemi,2)==1
            T = getrow(T,T.regSide==hemi);
        end
        roi=unique(T.roi);   
        
        A_all=zeros(numel(roi)*numel(sessN),numel(roi)*numel(sessN),numel(sn));
        for s=1:numel(sn)
            fprintf('Subject %d/%d\n',s,numel(sn));
            T1 = getrow(T,T.sn==sn(s));
            switch seqType
                case 'all'
                    A_all(:,:,s) = corrN(T1.partA_all',T1.partB_all');
                    %Conf(:,:,s) = ssqrt((T1.partA_all*T1.partA_all') * (T1.partB_all*T1.partB_all'));
                case 'trained' 
                    A_all(:,:,s) = corrN(T1.partA_train',T1.partB_train');
                    %Conf(:,:,s) = ssqrt((T1.partA_train*T1.partA_train') * (T1.partB_train*T1.partB_train'));
                case 'untrained'
                    A_all(:,:,s) = corrN(T1.partA_untrain',T1.partB_untrain');
                    %Conf(:,:,s) = ssqrt((T1.partA_untrain*T1.partA_untrain') * (T1.partB_untrain*T1.partB_untrain'));
            end
        end
        A = nanmean(A_all,3);
        Conf = diag(A);
        A = 1-A; % make into correlation distance
        % now normalize so max = 1
        A = (A-min(min(A)))./max(max(A));
        A(1:size(A,1)+1:end)=ones(size(A,1),1); % normalize diagonals
        Clust = kmeans(A,maxClust);
        Community = community_louvain(A);
        Roi=T1.roi;
        RegType=T1.regType;
        RegSide=T1.regSide;
        Conf_all=Conf;
        Conf=nanmean(Conf_all,3);
        sessN=T1.sessN;
        save(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)),'A','Conf','Clust','Community','Roi','RegSide','RegType','sessN');
        save(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s_individSubj',seqType,parcelType)),'A_all','Conf_all','sessN');
       % W_thr = threshold_proportional(W, 0.75);
        % make into a similarity
%         thres = mean(quantile(W,0.05));
%         S = exp(-W.^2./(2*thres^2));
%         
%         t = rsa_vectorizeRDM(W);
%         W1 = rsa_squareRDM(t);
%         W1(1:(size(W,1)+1):end)=diag(W);
%         
 
    case 'CLUSTER_consist_subj'
        % estimate the consistency of RDMs in each tessel
        % both within and across subjects
        sn=[4:9,11:28,30];
        sessN=[1:4];
        parcelType='162tessels'; % Brodmann or 162tessels
        sessType = 'within';
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType'});
        
        T = load(fullfile(betaDir,'group',sprintf('betas_partition_%s',parcelType)));
        % split participants by groups
        S{1}=sn(mod(sn,2)==1);
        S{2}=sn(mod(sn,2)==0);
        if strcmp(sessType,'within')
            T = getrow(T,T.sessN==sessN);
        end
        RR=[];
        roi=unique(T.roi);
        
        if strcmp(parcelType,'Brodmann')
            regSide = [ones(1,8) ones(1,8)*2];
            regType = [1:8 1:8];
        elseif strcmp(parcelType,'162tessels')
            regSide = ones(size(roi));
            regSide(roi>162)=2;
            regType = roi;
            regType(regSide==2)=regType(regSide==2)-162;
        end
        
        for g=1:2 % groups 1 and 2
            for s1=1:size(S{g},2)
                for s2=1:size(S{g},2)
                    for i=1:numel(roi)
                        D1 = getrow(T,T.sn==S{g}(s1) & T.roi==roi(i));
                        D2 = getrow(T,T.sn==S{g}(s2) & T.roi==roi(i));
                        if s1==s2
                            corrDist = corr(D1.partA_train',D2.partB_train');
                        else
                            c1 = corr(D1.partA_train',D2.partB_train');
                            c2 = corr(D1.partB_train',D2.partA_train');
                            c3 = corr(D1.partA_train',D2.partA_train');
                            c4 = corr(D1.partB_train',D2.partB_train');
                            corrDist = mean([c1 c2 c3 c4]);
                        end
                        R.corrDist  = corrDist;
                        R.roi       = D1.roi;
                        R.regType   = D1.regType;
                        R.regSide   = D2.regSide;
                        R.sn1       = S{g}(s1);
                        R.sn2       = S{g}(s2);
                        RR=addstruct(RR,R);
                    end
                end
                fprintf('Done group%d: %d/%d\n',g,s1,size(S{g},2));
            end
        end
        save(fullfile(distPscDir,sprintf('consist_interSubj_%s_sess%d',parcelType,sessN)),'-struct','RR');
    case 'PLOT_consist'
        sessN=4; 
        parcelType='162tessels';
        type = {'within','between'};
        vararginoptions(varargin,{'sessN','parcelType'});
        
        T = load(fullfile(distPscDir,sprintf('consist_interSubj_%s_sess%d',parcelType,sessN)));
        for h=1:2;
            % per hemisphere
            caretSDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
            C=caret_load(fullfile(caretSDir,sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            data=zeros(size(C.data,1),2);
            for t=1:2 % within / between subjects
                roi = unique(T.roi(T.regSide==h));
                for r=roi'
                    T1 = getrow(T,T.roi==r);
                    % calculate within or between subject consistency
                    rS(1) = nanmean(T1.corrDist(T1.sn1==T1.sn2));
                    rS(2) = nanmean(T1.corrDist(T1.sn1~=T1.sn2));
                    indx = C.index(C.data==r);
                    % mark in green
                    data(indx,t)=rS(t);
                    column_name{t} = fullfile(sprintf('consistRDM_%sSubj_sess%d.nii',type{t},sessN));
                end
            end
              R=caret_struct('metric','data',data,'column_name',column_name);
              caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.consistRDM_sess%d.metric',hem{h},sessN)),R);
        end
        
     
    case 'CLUSTER_surface_crossval'
        distType='community'; % community or crosscorr
        parcelType='162tessels';
        sessN=[1:2];
        clustN=10;
        sessType='within';
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'distType','clustN','parcelType','sessN','sessType','seqType'});
        
        T=load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)));
      
        switch distType
            case 'community'
                clustN=length(unique(T.Community));
            case 'crosscorr'
                clustN=length(unique(T.Clust));
        end
        colMap = [155 0 95; 10 41 92; 26 169 166; 152 95 153; 230 57 70;...
                    17 86 95; 193 100 94; 44 98 99; 32 164 243; 249 189 205];
        colMap = colMap./255;
        
        for h=1:2
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            RGBdata=zeros(size(C.data,1),3);
            for ss=sessN
                for k=1:clustN
                    S=getrow(T,T.RegSide==h & T.Clust==k & T.sessN==ss);
                    % data from that hemi, tessel, per cluster
                    indx=C.index(ismember(C.data,S.RegType));    % indicate the right tessel
                    RGBdata(indx,1)=colMap(k,1);
                    RGBdata(indx,2)=colMap(k,2);
                    RGBdata(indx,3)=colMap(k,3);
                    scale=[1;1;1];
                end
                scale(1:clustN,1)=0;
                scale(1:clustN,2)=1;
                
                name={sprintf('clusters%d_%s_%sSeq_%sSess%d',clustN,distType,seqType,sessType,ss)};
                R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.clusters-%d_%s_%sSeq_%sSess%d.RGB_paint',hem{h},clustN,distType,seqType,sessType,ss)),R);
            end
        end
    case 'CLUSTER_consist'
        sessN=[1:4]; 
        parcelType='162tessels';
        seqType='all'; % all or trained
        vararginoptions(varargin,{'sessN','parcelType','seqType'});
        
        T=load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)));
        for ss=sessN
            for h=1:2
                % per hemisphere
                caretSDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
                C=caret_load(fullfile(caretSDir,sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
                data=zeros(size(C.data,1),1);
                T1=getrow(T,T.sessN==ss & T.RegSide==h);
                roi = unique(T1.RegType(T1.RegSide==h));
                for r=1:length(roi)
                    indx = C.index(C.data==roi(r));
                    data(indx,:)=T1.Conf(r);
                end
                column_name{1} = fullfile(sprintf('consistCorr_sess%d_%sSeq.nii',ss,seqType));
                
                R=caret_struct('metric','data',data,'column_name',column_name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.consistCorr_%sSeq_sess%d.metric',hem{h},seqType,ss)),R);
            end
        end
        
     
    case 'CLUSTER_surface'
        distType='cosine';
        parcelType='162tessels';
        sessN=[1:2];
        clustN=10;
        sessType='within';
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'distType','clustN','parcelType','sessN','sessType'});
        
        colMap = dlmread(sprintf('my%dColors.txt',clustN));
        switch sessType
            case 'within'
               T=load(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d',parcelType,distType,clustN,seqType,sessType,sessN)));
            case 'across'
                T=load(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d-%d',parcelType,distType,clustN,seqType,sessType,sessN(1),sessN(end))));
        end
        for h=1:2;
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            RGBdata=zeros(size(C.data,1),3);
            for ss=sessN
                for k=1:clustN
                    S=getrow(T,T.regSide==h & T.cluster==k & T.sessN==ss);
                    % data from that hemi, tessel, per cluster
                    indx=C.index(ismember(C.data,S.regType));    % indicate the right tessel
                    RGBdata(indx,1)=colMap(k,1);
                    RGBdata(indx,2)=colMap(k,2);
                    RGBdata(indx,3)=colMap(k,3);
                    scale=[1;1;1];
                end
                scale(1:clustN,1)=0;
                scale(1:clustN,2)=1;
                
                name={sprintf('clusters%d_%sDist_%sSeq_%sSess%d',clustN,distType,seqType,sessType,ss)};
                R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.cluster-%d_%sDist_%sSeq_%sSess%d.RGB_paint',hem{h},clustN,distType,seqType,sessType,ss)),R);
            end
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
function M = pcm_defineSequenceModels(Seq,sn)
 % --------------------------------------
    % Model1: Null model - just common component to all sequences
    M{1}(:,1)    = ones(12,1);  % Common component 
    % --------------------------------------
    % Model2: Model with trained & untrained labels - 
    M{2}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{2}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model3: Model for each specific TRAINED sequence  
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{3}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model4: Model for each specific UNTRAINED sequence 
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{4}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
    % ---------------------------------------
    % Model5: First finger model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    
    M{5}(:,:)=eye(5); % first finger
    for u = 1:size(Seq,1)    % 12 seq
        firstfing = Seq(:,1);
        M{5}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model6: All fingers
    M{6}(:,:)=eye(5); % all fingers
    for u = 1:size(Seq,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{6}(u,j) = length(placenumb);
        end
    end
    %---------------------------------------
    % Model7: Finger transitions
    Trans = [nchoosek([1:5],2); fliplr(nchoosek([1:5],2))];
    M{7}(:,:)=zeros(12,size(Trans,1)); % all transitions (no repetition)
    for k = 1:size(Seq,1)    % set of 12 seq
        for p = 1:size(Trans,1) % all possible transitions
            if any(ismember(Seq(k,1:8),Trans(p,1)))
                ind = find(ismember(Seq(k,1:8),Trans(p,1)));    % find matching digits to first finger in doublet
                if any(ismember(Seq(k,ind+1),Trans(p,2)))       % compare the second finger in doublet
                    M{7}(k,p) = M{7}(k,p) + sum(ismember(Seq(k,ind+1),Trans(p,2)));
                end
            end
        end  
    end
end
function M = pcm_defineSequenceModels_fixed(Seq,sn)
 %  No fixed / run component added

    % --------------------------------------
    % Model1: Model with trained & untrained labels - 
    M{1}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{1}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model2: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        %A=zeros(6);
        A(i,i)=1;
       % M{2}(:,1:6 ,i)    = [A;zeros(6)];
    end;
    M{2}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model3: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        %A=zeros(6);
        A(i,i)=1;
        %M{3}(:,1:6,i) = [zeros(6);A];
    end;
    M{3}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
    % ---------------------------------------
    % Model4: First finger model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    
    M{4}(:,:)=eye(5); % first finger
    for u = 1:size(Seq,1)    % 12 seq
        firstfing = Seq(:,1);
        M{4}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model5: All fingers
    M{5}(:,:)=eye(5); % all fingers
    for u = 1:size(Seq,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{5}(u,j) = length(placenumb);
        end
    end
    %---------------------------------------
    % Model6: Finger transitions
    Trans = [nchoosek([1:5],2); fliplr(nchoosek([1:5],2))];
    M{6}(:,:)=zeros(12,size(Trans,1)); % all transitions (no repetition)
    for k = 1:size(Seq,1)    % set of 12 seq
        for p = 1:size(Trans,1) % all possible transitions
            if any(ismember(Seq(k,1:8),Trans(p,1)))
                ind = find(ismember(Seq(k,1:8),Trans(p,1)));    % find matching digits to first finger in doublet
                if any(ismember(Seq(k,ind+1),Trans(p,2)))       % compare the second finger in doublet
                    M{6}(k,p) = M{6}(k,p) + sum(ismember(Seq(k,ind+1),Trans(p,2)));
                end
            end
        end  
    end
end
function M = pcm_defineSequenceModels_new(Seq,sn)
%  No fixed / run component added
% specific sequence + overall component modelled together

    % --------------------------------------
    % Model1: Model with trained & untrained labels - 
    M{1}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{1}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model2: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{2}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model3: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{3}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
    % ---------------------------------------
    % Model4: First finger model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    
    M{4}(:,:)=eye(5); % first finger
    for u = 1:size(Seq,1)    % 12 seq
        firstfing = Seq(:,1);
        M{4}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model5: All fingers
    M{5}(:,:)=eye(5); % all fingers
    for u = 1:size(Seq,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{5}(u,j) = length(placenumb);
        end
    end
end
function M = pcm_defineSequenceFingerModels(Seq,sn)
%  No fixed / run component added
% specific sequence + overall component modelled together

    % --------------------------------------
    % Model1: Model with trained & untrained labels - 
    M{1}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{1}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model2: Model for each specific TRAINED + UNTRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{2}(:,1:6)    = [A;A];       % Trained and Untrained sequence patterns
    % ---------------------------------------
    % Model3: First finger model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    M{3}(:,:)=zeros(12,5); % first finger
    firstfing = Seq(:,1);
    for u = 1:size(Seq,1)    % 12 seq
        M{3}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model4: All fingers
    M{4}(:,:)=zeros(12,5); % all fingers
    for u = 1:size(Seq,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{4}(u,j) = length(placenumb);
        end
    end
end
function M = pcm_defineSequenceFingerModels_new(Seq,sn)
%  No fixed / run component added
% specific sequence + overall component modelled together

    % --------------------------------------
    % --------------------------------------
    % Model1: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{1}(:,1:6)  = [A;zeros(6)];       % Trained sequence patterns
    M{1}(:,7)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    % --------------------------------------
    % Model2: Model for each specific UNTRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{2}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
    M{2}(:,7)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % ---------------------------------------
    % Model3: First finger model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    M{3}(:,:)=zeros(12,5); % first finger
    firstfing = Seq(:,1);
    for u = 1:size(Seq,1)    % 12 seq
        M{3}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model4: All fingers
    M{4}(:,:)=zeros(12,5); % all fingers
    for u = 1:size(Seq,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{4}(u,j) = length(placenumb);
        end
    end
end

function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Crossvalidated model comparision:

    [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Group fit with %s algorithm done.\n',algorithm);
    %[T1,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize');
    %fprintf('Group fit with minimize algorithm done.\n');
    %[T2,theta_hat2,G_pred2,theta02] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR');
    %fprintf('Group fit with NR algorithm done.\n');
    
  %  [T3,theta_hat3,G_pred3,theta03] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR','theta0',theta_hat);
  %  [T4,theta_hat4,G_pred4,theta04] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize','theta0',theta_hat2);
    
    [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);

    T.cross_likelihood = Tcross.likelihood;
    T.bayesEst = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
    T.thetaCr = thetaCr;
   % T.modelNum = [1:length(T.cross_likelihood)];
end
function M = pcm_defineSequenceModels_seqType(Seq,sn,seqType)
%  Model for one seqType - trained or untrained sequences
% Run separately for trained and untrained
% seqType is index whether trained (1) or untrained (2)

    % --------------------------------------
    % Model1: Model for each specific sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{1}(:,1:6)    = A;       % Trained sequence patterns
    % --------------------------------------
    % Model2: First finger model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    if seqType==1
        SeqSet=Seq([1:6],:);
    elseif seqType==2
        SeqSet=Seq([7:12],:);
    end
    M{2}(:,:)=eye(5); % first finger
    for u = 1:size(SeqSet,1)    % 12 seq
        firstfing = SeqSet(:,1);
        M{2}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model3: All fingers
    M{3}(:,:)=eye(5); % all fingers
    for u = 1:size(SeqSet,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(SeqSet(u,:)==j);
            M{3}(u,j) = length(placenumb);
        end
    end
    %---------------------------------------
    % Model4: Finger transitions
    Trans = [nchoosek([1:5],2); fliplr(nchoosek([1:5],2))];
    M{4}(:,:)=zeros(6,size(Trans,1)); % all transitions (no repetition)
    for k = 1:size(SeqSet,1)    % set of 12 seq
        for p = 1:size(Trans,1) % all possible transitions
            if any(ismember(SeqSet(k,1:8),Trans(p,1)))
                ind = find(ismember(SeqSet(k,1:8),Trans(p,1)));    % find matching digits to first finger in doublet
                if any(ismember(SeqSet(k,ind+1),Trans(p,2)))       % compare the second finger in doublet
                    M{4}(k,p) = M{4}(k,p) + sum(ismember(SeqSet(k,ind+1),Trans(p,2)));
                end
            end
        end  
    end
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
function C = calcDist(D,betaW,G)
% calculates G and distances for all sequences, trained / untrained /
% between the two
% INPUT: D - structure with run / cond etc.
%        betaW - all betas
%        G - overall G structure (12x12 - all seq)
% OUTPUT: Do - containes distances, eigenvalues

    % calculate trained / untrained G
    G_train = G([1:6],[1:6]);
    G_untrain = G([7:12],[7:12]);
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

    % calculate distances - equivalent to rsa.distanceLDC
    dist_train = rsa.rdm.squareRDM(diag(ind6*G_train*ind6'));
    dist_untrain = rsa.rdm.squareRDM(diag(ind6*G_untrain*ind6'));
    dist_cross = diag(ind2*G_seqType*ind2');
    dist_all = rsa.rdm.squareRDM(diag(ind12*G*ind12'));
    C.RDM_train     = rsa_vectorizeRDM(dist_train);
    C.RDM_untrain   = rsa_vectorizeRDM(dist_untrain);
    % check - new!!!!! Mar8th 2018
    dist_train = triu(dist_train); dist_train(dist_train==0)=NaN;
    dist_untrain = triu(dist_untrain); dist_untrain(dist_untrain==0)=NaN;
    dist_cross = triu(dist_cross); dist_cross(dist_cross==0)=NaN;
    dist_all = triu(dist_all); dist_all(dist_all==0)=NaN;
    
    C.dist_train    = nanmean(dist_train(:));
    C.dist_untrain  = nanmean(dist_untrain(:));
    C.dist_cross    = nanmean(dist_cross);
    C.dist_all      = nanmean(dist_all(:));
 
    % calculate eigenvalues
    H=eye(6)-ones(6,6)./6;  % centering matrix!
    G_trainCent = H*G_train*H';  % double centered G matrix - rows and columns
    C.eigTrain = sort(eig(G_trainCent)','descend');    % sorted eigenvalues - trainedSeq
    G_untrainCent = H*G_untrain*H';  % double centered G matrix - rows and columns
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

function output = calcCorrSearch(Y,SPM,conditionVec,partitionVec)
% Called from rsa_runSearchlight in 'SEARCH_calcCorr'
% calculate correlation distance
[~,res_MS,~,beta_hat] = rsa.spm.noiseNormalizeBeta(Y,SPM);  % get raw beta regressor weights
betaUW  = bsxfun(@rdivide,beta_hat,sqrt(res_MS));           % apply univariate whitening to beta regressors (divide by voxel's variation)
[G,Sig] = pcm_estGCrossval(betaUW,partitionVec,conditionVec);
C       = corr_crossval(G,'reg','minvalue');
output  = C';
end