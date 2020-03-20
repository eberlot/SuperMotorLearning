function varargout=sml1_imana_dist(what,varargin)

% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
%baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
behavDir        =[baseDir '/behavioral_data/analyze'];
connectDir      =[baseDir '/connectivity'];
clusterDir      =[connectDir '/cluster'];
transformDir    =[connectDir '/transform'];
imagingDir      =[baseDir '/imaging_data'];                     
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia_new'];
%suitDir         =[baseDir '/suit'];
%physioDir       =[baseDir '/physio'];
pcmDir          =[baseDir '/pcm_stats'];
distPscDir      =[baseDir '/dist_psc_stats'];
QCDir           =[baseDir '/quality_control'];
wbDir           =[baseDir '/surfaceWB'];
% update glmDir when adding new glms
glmLocDir       ={[baseDir '/glmLoc/glmL1'],[baseDir '/glmLoc/glmL2'],[baseDir '/glmLoc/glmL3']};   % localiser glm
glmLocSessDir   ={[baseDir '/glmLocSeiess/glmLocSess1'],[baseDir '/glmLocSess/glmLocSess2'],[baseDir '/glmLocSess/glmLocSess3'],[baseDir '/glmLocSess/glmLocSess4']}; % one glm for loc run per session
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session
glmErrorDir{1}  ={[baseDir '/glmSessError1/glmSessError1'],[baseDir '/glmSessError1/glmSessError2'],[baseDir '/glmSessError1/glmSessError3'],[baseDir '/glmSessError1/glmSessError4']}; % one glm per session
glmErrorDir{2}  ={[baseDir '/glmSessError2/glmSessError1'],[baseDir '/glmSessError2/glmSessError2'],[baseDir '/glmSessError2/glmSessError3'],[baseDir '/glmSessError2/glmSessError4']}; % one glm per session
glmErrorDir{3}  ={[baseDir '/glmSessError3/glmSessError1'],[baseDir '/glmSessError3/glmSessError2'],[baseDir '/glmSessError3/glmSessError3'],[baseDir '/glmSessError3/glmSessError4']}; % one glm per session

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
Chunks = [1 5 3; 5 3 4; 3 1 5; 1 5 2; 1 4 5; 3 4 2; 5 4 2;...
          5 1 3; 1 3 2; 3 5 1; 5 1 4; 5 2 1; 3 2 4; 1 2 4];

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

% per session
numruns_sess      = 10;  
numruns_task_sess = 8;
numruns_loc_sess  = 2;

% total - per subject (the total in the end will always be 40)
numruns           = [40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40];
numruns_task      = 32;
numruns_loc       = 8;
sess = [repmat(1,1,10),repmat(2,1,10),repmat(3,1,10),repmat(4,1,10)];   % all sessions

sess_sn = [4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4];    % per subject

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
hem       = {'lh','rh'};
hemI      = {'L','R'};                                                      % freesurfer WB: initials
hemname   = {'CortexLeft','CortexRight'};
% ------------------------- Subject things --------------------------------
% The variables in this section must be updated for every new subject.

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
              's11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
              's21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'};  


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

stySess=style.custom({black,gray,lightgray,blue},'markersize',ms);
styTrained_sess=style.custom({red,mediumred,lightred},'markersize',ms);
styUntrained_sess=style.custom({blue,mediumblue,lightblue},'markersize',ms);

% ------------------------------ Analysis Cases --------------------------------
switch(what)
     
    case 'PSC_create'
    % calculate psc for trained and untrained sequences - based on betas  
    sn = [5:9,11:31];
    sessN = 1:4;
    errorNum = 1;
    vararginoptions(varargin,{'sn','sessN','errorNum'});
    %name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','TrainSeq','UntrainSeq','AllSeq'};
    name={'TrainSeq','UntrainSeq','AllSeq'};
    for ss=sessN
        for s=sn
           cd(fullfile(glmErrorDir{errorNum}{ss},subj_name{s}));
           %cd(fullfile(glmSessDir{ss},subj_name{s}));
            load SPM;
            T=load('SPM_info.mat');
            X=(SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
            h=median(max(X));               % Height of response;
            P={};
            numB=length(SPM.xX.iB);         % Partitions - runs
            for p=SPM.xX.iB
                P{end+1}=sprintf('beta_%4.4d.nii',p);       % get the intercepts and use them to calculate the baseline (mean images)
            end;
            for con=1:length(name)    % all contrasts
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
        sessN=[1:4];
        sn=[5:9,11:31];
        smooth = 1;   
        vararginoptions(varargin,{'sn','sessN','smooth'});

        hemisphere=1:length(hem);
        fileList = [];
        column_name = [];
        name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','AllSeq','TrainSeq','UntrainSeq'};
        
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
                    metric_out = fullfile(caretSDir,sprintf('%s_psc_sess%d.metric',subj_name{s},ss));
                    M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                    M.column_name = column_name;
                    caret_save(metric_out,M);
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        cd(caretSDir);
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        Out = caret_smooth(metric_out, 'coord', white, 'topo', closed);%,...
                        %'algorithm','FWHM','fwhm',12);
                        char(Out);  % if smoothed adds an 's'
                    end;
                    fprintf('Sess %d, Subj %d, Hem %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session
    case 'PSC_create_loc' % depreciated
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
        sessN=[1:4];
        sn=[4:9,11:31];
        vararginoptions(varargin,{'sessN','sn'});
        % Some presets
        name    = 'Contrasts';
        
        OUTname    = {'TrainSeq','UntrainSeq'};
        inputcol   = [13 14]; % only averages - all trained / untrained
        replaceNaN = [1 1];
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
                    fprintf('hem: %i  image: %i \n', h,j);
                end;
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
        
        sessN=[1:4];
        sn=[4:9,11:31];
        s=1:length(sn);
        vararginoptions(varargin,{'sessN','sn'});
        
        SPMname={'TrainSeq','UntrainSeq'};
        
        sqrtTransform=[0,0,0,0]; % Should you take ssqrt before submitting? 
        % no for psc
        for ss=sessN
            SummaryName = sprintf('.summary_psc_sess%d.metric',ss);
            hemi = [1 2];
            
            for h=hemi
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
            fprintf('Done sess-%d\n',ss);
        end
    case 'PSC_group_smooth'
        sessN=[1:4];
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
                    filename=[surfaceGroupDir filesep hem{h} '.summary_psc_sess' num2str(ss) '.metric']; % unsmoothed
                    sfilename=caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', 15);
                end;  
            end;
        end
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

    case 'search_job'
        sml1_imana_dist('SURF_wb:non-permuted','sessN',3:4);
        sml1_imana_dist('SURF_wb:non-permuted','sessN',2:4,'metric','psc');
        sml1_imana_dist('SURF_wb:non-permuted','sessN',1,'metric','psc');
        sml1_imana_dist('SURF_wb:permute_stats','sessN',3:4);
        sml1_imana_dist('SURF_wb:permute_stats','sessN',2:4,'metric','psc');
        sml1_imana_dist('SURF_wb:permute_stats','sessN',1,'metric','psc');
    case 'SURF_psc'
     % create surface maps of percent signal change 
     % trained and untrained sequences
        sessN=[1:4];
        sn=[5:9,11:31];
        smooth = 1;   
        vararginoptions(varargin,{'sn','sessN','smooth'});

        hemisphere=1:length(hem);
        fileList = [];
        column_name = [];
        name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','AllSeq','TrainSeq','UntrainSeq'};
        
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
                    metric_out = fullfile(caretSDir,sprintf('%s_psc_sess%d.metric',subj_name{s},ss));
                    M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                    M.column_name = column_name;
                    caret_save(metric_out,M);
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        cd(caretSDir);
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        Out = caret_smooth(metric_out, 'coord', white, 'topo', closed);%,...
                        %'algorithm','FWHM','fwhm',12);
                        char(Out);  % if smoothed adds an 's'
                    end;
                    fprintf('Sess %d, Subj %d, Hem %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session
    case 'SURF_diff'                                                 % STEP 4.4a   :  Map searchlight results (.nii) onto surface (.metric)
        % map volume images to metric file and save them in individual surface folder
        sn      = [5:9,11:31];
        sessN   = [1:4];   
        fileList = {'dist','dist_trained','dist_untrained','dist_cross'}; % for dist or repsup
        glmDir = glmSessDir;
        outname = 'dist'; % dist or dist_repsup
        smooth = 1;
        vararginoptions(varargin,{'sn','sessN','fileList','glmDir','outname'});
        
        for ss = sessN
            for s = sn
                for h=1:2
                    caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                    white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
                    pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
                    
                    for f = 1:length(fileList)
                        images{f}    = fullfile(glmDir{ss},subj_name{s},sprintf('%s_sess%d_%s.nii',subj_name{s},ss,fileList{f}));
                        column_name{f} = fullfile(sprintf('Sess%d_%s.nii',ss,fileList{f}));
                    end;    % filename
                    outfile   = sprintf('%s_%s_sess%d.metric',subj_name{s},outname,ss);
                    M         = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
                    M.column_name = column_name;
                    caret_save(fullfile(caretSDir,outfile),M);
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        white = fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                        Out = caret_smooth(fullfile(caretSDir,outfile), 'coord',white,'topo',closed);%,...
                        %'algorithm','FWHM','fwhm',12);
                        char(Out);  % if smoothed adds an 's'
                    end;
                    fprintf('Done subj %d sess %d hemi %d \n',s,ss,h);
                end;    % hemi
            end;    % subj
        end; % sessN
    case 'SURF_seqDiff'
        % trained - untrained difference per session
        metric = 'dist';
        sn = [5:9,11:31];
        sessN = 1:4;
        smooth = 1;
        vararginoptions(varargin,{'dist','smooth','metric','sn'});
        switch metric
            case 'dist'
                  col = [2,3];
            case 'psc'
                  col = [14,15];
        end
        for s=sn
            for h=1:2 % hemi
                caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                for ss=sessN
                    M = caret_load(fullfile(caretSDir,sprintf('%s_%s_sess%d.metric',subj_name{s},metric,ss)));
                    data(:,1) = M.data(:,col(1))-M.data(:,col(2));
                    column_name{1} = sprintf('%s_sess%d',metric,ss);
                    C = caret_struct('metric','data',data,'column_name',column_name);
                    fileName = sprintf('s%s_%s_seqDiff_sess%d.metric',subj_name{s},metric,ss);
                    caret_save(fullfile(caretSDir,fileName),C);
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        Out = caret_smooth(fullfile(caretSDir,fileName), 'coord', white, 'topo', closed);%,...
                        %'algorithm','FWHM','fwhm',12);
                        char(Out);  % if smoothed adds an 's'
                    end;
                    clear data;
                end; % sess
            end; % hemi
            fprintf('Subj %d\n',s);
        end; % sn
    case 'SURF_sessDiff'
        % differences across sessions
        metric = 'dist';
        sn = [5:9,11:31];
        smooth = 1;
        vararginoptions(varargin,{'smooth','metric'});
        switch metric
            case 'dist'
                  col = [2,3];
            case 'psc'
                  col = [14,15];
        end
        colName = {'trained','untrained'};
        for s=sn
            for h=1:2 % hemi
                caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                for ss=1:3
                    M1 = caret_load(fullfile(caretSDir,sprintf('s%s_%s_sess%d.metric',subj_name{s},metric,ss)));
                    M2 = caret_load(fullfile(caretSDir,sprintf('s%s_%s_sess%d.metric',subj_name{s},metric,ss+1)));
                    for i=1:size(col,2)
                        data(:,i) = M2.data(:,col(i))-M1.data(:,col(i));
                        column_name{i} = sprintf('%s_%s_sessDiff',metric,colName{i});
                    end                    
                    C = caret_struct('metric','data',data,'column_name',column_name);
                    fileName = sprintf('%s_%s_sessDiff_sess%d.metric',subj_name{s},metric,ss);
                    caret_save(fullfile(caretSDir,fileName),C);
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        Out = caret_smooth(fullfile(caretSDir,fileName), 'coord', white, 'topo', closed);
                        char(Out);  % if smoothed adds an 's'
                    end;
                    clear data;
                end; % sess
            end; % hemi
            fprintf('Subj %d\n',s);
        end; % sn
    case 'SURF_groupMake'
        % Calculate group metric files 
        sessN=1:4;
        sn=[5:9,11:31];
        name = 'dist'; % psc / dist / psc_seqDiff / dist_seqDiff
        inputcol   = [1 2 3 4]; % [13,14,15] for psc
        replaceNaN = [1 1 1 1];
        smooth = 1;
        OUTname    = {'dist_all','dist_trained','dist_untrained','dist_cross'}; 
        %psc: {'psc_all','psc_trained','psc_untrained'}
        %psc_seqDiff: psc_seqDiff
        %dist_seqDiff: dist_seqDiff
        %psc_sessDif: {'psc_trained_sessDiff','psc_untrained_sessDiff'}
        vararginoptions(varargin,{'sessN','sn','name','OUTname','inputcol','replaceNaN'});
        
        % Loop over hemispheres.
        for h = 1:2
            % Go to the directory where the group surface atlas resides
            surfaceGroupDir = [caretDir filesep atlasname filesep hemName{h}];
            cd(surfaceGroupDir);
            for ss=sessN
                % Loop over each input metric file in 'OUTname' and make a group metric file
                for j = 1:length(OUTname);
                    % Loop over subjects...
                    for i = 1:length(sn);
                        % ...and define the names of their metric files
                        switch smooth
                            case 0
                                infilenames{j}{i} = fullfile(caretDir,[atlasA subj_name{sn(i)}], hemName{h}, sprintf('%s_%s_sess%d.metric',subj_name{sn(i)},name,ss));
                            case 1
                                infilenames{j}{i} = fullfile(caretDir,[atlasA subj_name{sn(i)}], hemName{h}, sprintf('s%s_%s_sess%d.metric',subj_name{sn(i)},name,ss));
                        end
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    switch smooth
                        case 0
                            outfilenames{j} = [surfaceGroupDir filesep hem{h} '.' OUTname{j} '_sess' num2str(ss) '.metric'];
                        case 1
                            outfilenames{j} = [surfaceGroupDir filesep 's' hem{h} '.' OUTname{j} '_sess' num2str(ss) '.metric'];
                    end
                    % Finally, make the group metric file for this metric type/contrast
                    caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
                    % Verbose display to user
                    fprintf('sess: %d  hem: %i  image: %i \n',ss,h,j);
                end;
            end
        end;
    case 'SURF_group_cSPM'
        % Calculate group stats files from the group metric files. 
        sessN=1:4;
        sn=[5:9,11:31];
        name={'dist_all','dist_trained','dist_untrained','dist_cross'};
        %{'psc_all','psc_trained','psc_untrained'} {'psc_seqDiff'},{'dist_seqDiff'},{'psc_sessDiff_trained','psc_sessDiff_untrained'}
        metric='dist'; % psc / dist / psc_seqDiff / dist_seqDiff
        sqrtTransform=[1,1,1,1]; % Should you take ssqrt before submitting? for distances
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric'});
        s=1:length(sn);
        hemi = [1 2];
        for ss=sessN
            SummaryName = sprintf('.summary_%s_sess%d.metric',metric,ss);
            for h=hemi
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
                %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
                for i=1:length(name);
                    %sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' name{i} '.metric']; % smoothed
                    sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' name{i} '_sess' num2str(ss) '.metric']; % smoothed maps
                    %sfilenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '_sess' num2str(ss) '.metric']; % no smoothing
                end;
                %----loop over the metric files and calculate the cSPM of each with the non-smoothed metrics
                for i=1:length(name);
                    Data=caret_load(sfilenames{i});
                    if sqrtTransform(i)
                        Data.data=ssqrt(Data.data);
                    end;
                    SPMname{i} = sprintf('%s_sess%d',name{i},ss);
                    cSPM=caret_getcSPM('onesample_t','data',Data.data(:,s),'maskthreshold',0.7); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                    caret_savecSPM([surfaceGroupDir filesep hem{h} '.' SPMname{i} '_stats.metric'],cSPM);
                    save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');
                    data(:,i)=cSPM.con(1).con; % mean
                    data(:,i+length(name))=cSPM.con(1).Z; % T
                    column_name{i}=['mean_' SPMname{i} '_sess' num2str(ss)];
                    column_name{i+length(name)}=['T_' SPMname{i} '_sess' num2str(ss)];
                end;
                C = caret_struct('metric','data',data,'column_name',column_name);
                caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);
            end;
        end
        fprintf('Done \n');
    case 'SURF_group_list'
         % TabDat=caret_list(Surface,cSPM,u,k,varargin);
        hemi        = 1;
        cSPMname    = 'cSPM_psc_seqDiff_sess2.mat';
        cSPMdir     = caretDir;
        th_hight    = 0.01; % p-val
        th_size     = 1;% cluster size mm^2 or corrected p-value
        vararginoptions(varargin,{'cSPMname','th_hight','th_size','cSPMdir'});
        for h=hemi
            if strcmp(cSPMdir,caretDir)
                surfaceGroupDir=[cSPMdir filesep 'fsaverage_sym' filesep hemName{h}];
            else
                surfaceGroupDir = cSPMdir; % for permutations
            end
            %---- define name of coord and topology
            coordfile   = [caretDir filesep 'fsaverage_sym' filesep hemName{h} filesep hem{h} '.WHITE.coord'];
            topofile    = [caretDir filesep 'fsaverage_sym' filesep hemName{h} filesep hem{h} '.CLOSED.topo'];        
            %---- get surface
            surface = caret_getsurface(coordfile,topofile);  
            %---- load reference .paint file
            P = caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},'ROI.paint'));      
            %---- load cSPM
            cSPMfile = fullfile(surfaceGroupDir,cSPMname);
            load(cSPMfile);         
            %---- get table
            table = caret_list(surface,cSPM,th_hight,th_size,'label',P,'mask',P.data>0,'sort_by','area');
            varargout{1} = table;
        end;
    case 'SURF_permute_stats'
        % implementing permutation tests for assessing statistical
        % significance of detected clusters
        sessN=1:4;
        sn=[5:9,11:31];
        nPerm = 1000;
        name={'dist_all','dist_trained','dist_untrained','dist_cross'};
        %{'psc_all','psc_trained','psc_untrained'} {'psc_seqDiff'},{'dist_seqDiff'},{'psc_sessDiff_trained','psc_sessDiff_untrained'}
        metric='dist'; % psc / dist / psc_seqDiff / dist_seqDiff
        sqrtTransform=[1,1,1,1]; % Should you take ssqrt before submitting? for distances
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','nPerm'});
        s=1:length(sn);
        hemi = 1;
        % load permutations
        load(fullfile(caretDir,'PermComb'));       
        perm = randi(size(PermComb,1),1,1000); % choose at random
        PermComb = PermComb(perm,:); % here subset of permutations - randomly picked
        fprintf('Running permutations for %s.',metric);
        for ss=sessN
            fprintf('\nSession-%d:\n');
            for h=hemi
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
                permDir = fullfile(surfaceGroupDir,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                dircheck(permDir);
                %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
                for i=1:length(name)
                    sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' name{i} '_sess' num2str(ss) '.metric']; % smoothed maps
                    Data(i)=caret_load(sfilenames{i});
                    if sqrtTransform(i)
                        Data(i).data=ssqrt(Data(i).data);
                    end;
                end
                for i=1:length(name);
                    for p=1:nPerm
                        tElapsed = tic;
                        %----loop over the metric files and calculate the cSPM for each permuted file
                        permData = bsxfun(@times,Data(i).data(:,s),PermComb(p,:));
                        SPMname = sprintf('%s_sess%d_perm-%d',name{i},ss,p);
                        cSPM=caret_getcSPM('onesample_t','data',permData,'maskthreshold',0.7); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                        save([permDir  filesep   'cSPM_' SPMname '.mat'],'cSPM');
                        table = sml1_imana_dist('SURF_group_list','cSPMname',sprintf('cSPM_%s.mat',SPMname),'cSPMdir',permDir);
                        save(fullfile(permDir,sprintf('table_permutation-%d.mat',p)),'-struct','table')
                        fprintf('%d.',p); toc(tElapsed);
                    end;
                end;
            end
        end
        fprintf('Done \n');
    case 'SURF_permute'
        % here calculate many possible permutations for flipping
        nSubj = 26;
        Comb = [];
        randRow = randi(nSubj-1,1,10);
        for i=1:randRow
            A=nchoosek([1:nSubj],i);
            n=size(A,1);
            X=zeros(n,nSubj);
            for j=1:n
                X(j,A(j,:))=1;
            end;
            Comb=[Comb;X];
        end;
        Comb(Comb==0)=-1;
        PermComb = Comb;
        save(fullfile(caretDir,'PermComb'),'-v7.3');
    case 'SURF_wb:permute_stats'
        sessN=1:4;
        nPerm = 1000;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='dist'; % psc / dist 
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','nPerm'});
        INname = sprintf('%sDifference',metric);
        hemi = 1;
        if strcmp(metric,'dist')
            kernel=1;
        else
            kernel=2;
        end
        % load permutations
        load(fullfile(caretDir,'PermComb'));       
        perm = randi(size(PermComb,1),1,1000); % choose at random
        PermComb = PermComb(perm,:); % here subset of permutations - randomly picked
        fprintf('Running permutations for %s.',metric);
        for ss=sessN
            fprintf('\nSession-%d:\n',ss);
            for h=hemi
                surfaceGroupDir = fullfile(wbDir,atlas);
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                dircheck(permDir);
                inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname,ss));
                G1 = gifti(inFile);
                for p=623:nPerm %p=1:nPerm 
                    tElapsed = tic;
                    %----loop over the metric files and calculate the cSPM for each permuted file
                    permData = bsxfun(@times,G1.cdata,PermComb(p,:));
                    cSPM = surf_getcSPM('onesample_t','data',permData,'maskthreshold',0.7);
                    SPMname = sprintf('%s_sess%d_perm-%d',INname,ss,p);
                    save(fullfile(permDir,sprintf('cSPM_%s.mat',SPMname)),'cSPM');
                    % here map
                    data(:,1)=cSPM.con(1).Z; % T
                    column_name(1)={'T_map'};
                    G2 = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
                    summaryName = fullfile(permDir,sprintf('map_perm-%d.func.gii',p));
                    save(G2,summaryName);
                    % here smooth
                    surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
                    surf_smooth(summaryName,'surf',surfFile,'kernel',kernel);
                    smoothName = fullfile(permDir,sprintf('smap_perm-%d.func.gii',p));
                    % here find clusters
%                     clustName = fullfile(permDir,sprintf('clusters_perm-%d_thres1.func.gii',p));
%                     comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.3,3,clustName);
%                     % 1.7 - T corresponding to p<0.05,one -sided
%                     [err,out]=system(comm);
%                     clustName = fullfile(permDir,sprintf('clusters_perm-%d_thres2.func.gii',p));
%                     comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.7,3,clustName);
%                     % 1.7 - T corresponding to p<0.05,two-sided
%                     [err,out]=system(comm);
                    clustName = fullfile(permDir,sprintf('clusters_perm-%d_thres3.func.gii',p));
                    comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,2.49,3,clustName);
                    % 2.49 - T corresponding to p<0.01
                    [err,out]=system(comm);
                    fprintf('%d.',p); toc(tElapsed);
                end;
            end
        end
        fprintf('Done \n');
    case 'SURF_wb:non-permuted'
        % same procedure but no permutation
        sessN=1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='dist'; % psc / dist 
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','nPerm'});
        INname = sprintf('%sDifference',metric);
        hemi = 1;
        if strcmp(metric,'dist')
            kernel=1;
        else
            kernel=2;
        end
        fprintf('Running permutations for %s.',metric);
        for ss=sessN
            fprintf('\nSession-%d:\n',ss);
            for h=hemi
                surfaceGroupDir = fullfile(wbDir,atlas);
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                dircheck(permDir);
                inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname,ss));
                G1 = gifti(inFile);
                
                permData = G1.cdata;
                cSPM = surf_getcSPM('onesample_t','data',permData,'maskthreshold',0.7);
                SPMname = sprintf('%s_sess%d_non-permuted',INname,ss);
                save(fullfile(permDir,sprintf('cSPM_%s.mat',SPMname)),'cSPM');
                % here map
                data(:,1)=cSPM.con(1).Z; % T
                column_name(1)={'T_map'};
                G2 = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
                summaryName = fullfile(permDir,sprintf('map_non-permuted.func.gii'));
                save(G2,summaryName);
                % here smooth
                surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
                surf_smooth(summaryName,'surf',surfFile,'kernel',kernel);
                smoothName = fullfile(permDir,sprintf('smap_non-permuted.func.gii'));
                % here find clusters
                clustName = fullfile(permDir,sprintf('clusters_non-permuted_thres1.func.gii'));
                comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.3,2,clustName);
                % 1.3 - T corresponding to p<0.05, one-sided
                [err,out]=system(comm);
                clustName = fullfile(permDir,sprintf('clusters_non-permuted_thres2.func.gii'));
                comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.7,2,clustName);
                % 1.3 - T corresponding to p<0.05, two-sided
                [err,out]=system(comm);
                clustName = fullfile(permDir,sprintf('clusters_non-permuted_thres3.func.gii'));
                comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,2.49,2,clustName);
                % 2.49 - T corresponding to p<0.01
                [err,out]=system(comm);
            end
        end
        fprintf('Done \n');
    case 'SURF_wb:non-permuted_clusters'
        sessN=1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='dist'; % psc / dist 
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','nPerm'});
        hemi = 1;
        for ss=sessN
            for h=hemi
                surfaceGroupDir = fullfile(wbDir,atlas);
                surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                for i=1:2
                    if i==1
                        smoothName = fullfile(permDir,sprintf('smap_non-permuted.func.gii'));
                        name = 'non-permuted';
                    else
                        S = gifti(fullfile(permDir,sprintf('smap_non-permuted.func.gii')));
                        column_name(1)={'T_map_reversed'};
                        G = surf_makeFuncGifti(S.cdata.*(-1),'anatomicalStruct',hemname{h},'columnNames',column_name);
                        smoothName = fullfile(permDir,sprintf('map_non-permuted_reversed.func.gii'));
                        save(G,smoothName);
                        name = 'non-permuted_reversed';
                    end
                    % here find clusters
                    clustName = fullfile(permDir,sprintf('clusters_%s_thres1.func.gii',name));
                    comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.3,2,clustName);
                    % 1.3 - T corresponding to p<0.05, one-sided
                    [err,out]=system(comm);
                    clustName = fullfile(permDir,sprintf('clusters_%s_thres2.func.gii',name));
                    comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.7,2,clustName);
                    % 1.3 - T corresponding to p<0.05, two-sided
                    [err,out]=system(comm);
                    clustName = fullfile(permDir,sprintf('clusters_%s_thres3.func.gii',name));
                    comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,2.49,2,clustName);
                    % 2.49 - T corresponding to p<0.01
                    [err,out]=system(comm);
                end
            end
            fprintf('Session-%d.\n',ss);
        end
    case 'SURF_wb:permuted_clusters'
        sessN=1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='dist'; % psc / dist 
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','nPerm'});
        hemi = 1;
        for ss=sessN
            fprintf('\nSession-%d.\n',ss);
            for h=hemi
                surfaceGroupDir = fullfile(wbDir,atlas);
                surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                for p=1:1000 % permutations
                    for i=1:2
                        if i==1
                            smoothName = fullfile(permDir,sprintf('smap_perm-%d.func.gii',p));
                            name = sprintf('perm-%d',p);
                        else
                            S = gifti(fullfile(permDir,sprintf('smap_perm-%d.func.gii',p)));
                            column_name(1)={'T_map_reversed'};
                            G = surf_makeFuncGifti(S.cdata.*(-1),'anatomicalStruct',hemname{h},'columnNames',column_name);
                            smoothName = fullfile(permDir,sprintf('smap_perm-%d_reversed.func.gii',p));
                            save(G,smoothName);
                            name = sprintf('perm-%d_reversed',p);
                        end
                        % here find clusters
                        if strcmp(metric,'dist')
                         %   clustName = fullfile(permDir,sprintf('clusters_%s_thres1.func.gii',name));
                         %   comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.3,.5,clustName);
                            % 1.3 - T corresponding to p<0.05, one-sided
                         %   [err,out]=system(comm);
                            clustName = fullfile(permDir,sprintf('clusters_%s_thres2.func.gii',name));
                            comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,1.7,2,clustName); % 2 / 0.5
                            % 1.7 - T corresponding to p<0.05, two-sided
                            [err,out]=system(comm);
                        end
                        clustName = fullfile(permDir,sprintf('clusters_%s_thres3.func.gii',name)); % for psc only the most stringent
                        comm = sprintf('wb_command -metric-find-clusters %s %s %d %d %s',surfFile,smoothName,2.49,2,clustName); % 2 / 0.5
                        % 2.49 - T corresponding to p<0.01
                        [err,out]=system(comm);
                    end
                    fprintf('%d.',p);
                end
            end
        end
    case 'SURF_wb:cluster_size_T'
        % here determine the size and peak T value for clusters
        sessN=1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='psc'; % psc / dist 
        thres = 3;
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','thres'});
        hemi = 1;
        TT = [];
        for ss=sessN
            for h=hemi
                G = gifti(fullfile(wbDir,atlas,'fs_LR.164k.coord.func.gii')); % for left hemisphere only
                S = gifti(fullfile(wbDir,atlas,'fs_LR.164k.surfArea.func.gii')); % for left hemisphere only - surface area
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                for i=1:2
                    if i==1
                        mapName = fullfile(permDir,sprintf('smap_non-permuted.func.gii'));
                        name = 'non-permuted';
                    else
                        mapName = fullfile(permDir,sprintf('map_non-permuted_reversed.func.gii'));
                        name = 'non-permuted_reversed';
                    end
                    fprintf('%s:\n\n',name)
                    M = gifti(mapName);
                    C = gifti(fullfile(permDir,sprintf('clusters_%s_thres%d.func.gii',name,thres)));
                    for c=1:max(unique(C.cdata))
                        fprintf('Cluster-%d.\n',c);
                        ind = C.cdata==c;
                        coord = double(G.cdata(ind,:));
                        shp = alphaShape(coord);
                        % here determine the shape
                        T.maxT = max(M.cdata(C.cdata==c));
                        T.numNodes = numel(M.cdata(C.cdata==c));
                        T.surfArea = sum(S.cdata(C.cdata==cc));
                        T.volume = volume(shp);
                        T.cluster = c;
                        T.positive_negative = i;
                        T.sessN = ss;
                        TT = addstruct(TT,T);
                    end
                end
            end
            fprintf('Session-%d done.\n',ss);
        end
        save(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterData_thres-%d',thres)),'-struct','TT');
    case 'SURF_wb:permutation_histogramT'
        % determine the distribution of Ts (+ / -) per session
        sessN=1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='psc'; % psc / dist 
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','thres'});
        hemi = 1;
        for ss=sessN
            fprintf('\nSession-%d:\n',ss);
            for h=hemi
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                T_hist = zeros(2,1000);
                for p=1:1000
                    G = gifti(fullfile(permDir,sprintf('smap_perm-%d.func.gii',p)));
                    T_hist(1,p) = max(G.cdata);
                    T_hist(2,p) = min(G.cdata);
                    fprintf('%d.',p);
                end
                fprintf('...done.\n');
            end
            save(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('T_perm_distribution_sess-%d',ss)),'T_hist');
        end
    case 'SURF_wb:permutation_histogram_clusterSize'
        % determine the distribution of cluster sizes (+ / -) per session
        sessN=1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='psc'; % psc / dist 
        thres = 3;
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','thres'});
        hemi = 1;
        for ss=sessN
            fprintf('\nSession-%d:\n',ss);
            for h=hemi
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                nNode = cell(2,1);
                for p=1:1000
                    for i=1:2 % positive and negative
                        if i==1
                            name = sprintf('perm-%d',p);
                        else
                            name = sprintf('perm-%d_reversed',p);
                        end
                        clustName = fullfile(permDir,sprintf('clusters_%s_thres%d.func.gii',name,thres));
                        G = gifti(clustName);
                        if numel(unique(G.cdata))~=1
                            for c=1:max(unique(G.cdata))
                                nNode{i} = [nNode{i},numel(find(G.cdata==c))];
                            end
                        end
                    end
                    fprintf('%d.',p);
                end
                fprintf('...done.\n');
            end
            save(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterSize_distribution_thres-%d_sess-%d',thres,ss)),'nNode');
        end
    case 'SURF_wb:permutation_histogram_clusterSize_v2'
        % determine the distribution of cluster sizes (+ / -) per session
        % only take the largest cluster
        sessN=2:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='psc'; % psc / dist 
        thres = 3;
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','thres'});
        hemi = 1;
        for ss=sessN
            fprintf('\nSession-%d:\n',ss);
            for h=hemi
                permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
                nNode = cell(2,1);
                for p=1:1000
                    for i=1:2 % positive and negative
                        if i==1
                            name = sprintf('perm-%d',p);
                        else
                            name = sprintf('perm-%d_reversed',p);
                        end
                        clustName = fullfile(permDir,sprintf('clusters_%s_thres%d.func.gii',name,thres));
                        G = gifti(clustName);
                        if numel(unique(G.cdata))~=1
                            nNode_perm = zeros(max(unique(G.cdata)),1);
                            for c=1:max(unique(G.cdata))
                                nNode_perm(c) = numel(find(G.cdata==c));
                            end
                        else
                            nNode_perm = 0;
                        end
                        nNode{i} = [nNode{i} max(nNode_perm)];
                    end
                    fprintf('%d.',p);
                end
                fprintf('...done.\n');
            end
            save(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterSize_distribution_thres-%d_sess-%d_v2',thres,ss)),'nNode');
        end
    case 'SURF_wb:determine_T_significance'
        % determine significance of clusters
        sessN=2:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='psc'; % psc / dist 
        thres = 3;
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','thres'});
        PP=[];
        T = load(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterData_thres-%d.mat',thres)));
        for ss=sessN
            load(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('T_perm_distribution_sess-%d.mat',ss))); % T_perm distribution
            for i=1:2 % positive_negative
                T1 = getrow(T,T.sessN==ss & T.positive_negative==i);
                for c=1:size(T1.cluster,1)
                    P = getrow(T1,T1.cluster==c);
                    P.T_p_unc = tpdf(P.maxT,25); % 25 degrees of freedom
                    P.T_p_corr = sum(numel(find(abs(T_hist(i,:))>P.maxT)))/1000;
                    PP = addstruct(PP,P);
                end
            end
        end
        save(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterData_thres-%d_T_p',thres)),'-struct','PP');
    case 'SURF_wb:determine_clusterSize_significance'
        % determine significance of size of clusters
        sessN=2:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='psc'; % psc / dist 
        thres = 3;
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','thres'});
        PP=[];
        T = load(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterData_thres-%d_T_p.mat',thres)));
        for ss=sessN
           % load(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterSize_distribution_thres-%d_sess-%d_v2.mat',thres,ss))); % T_perm distribution
            load(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterSize_distribution_thres-%d_sess-%d.mat',thres,ss))); % T_perm distribution
            for i=1:2 % positive_negative
                T1 = getrow(T,T.sessN==ss & T.positive_negative==i);
                for c=1:size(T1.cluster,1)
                    P = getrow(T1,T1.cluster==c);
                    P.clustSize_p_corr = sum(numel(find(nNode{i}>P.numNodes)))/size(nNode{1},2);
                    PP = addstruct(PP,P);
                end
            end
        end
        save(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterData_thres-%d_T_p_clusterSize_v2',thres)),'-struct','PP');
    case 'SURF_wb:add_surfaceArea'
        sessN=2:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric='psc'; % psc / dist
        thres = 3;
        hemi = 1;
        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','name','metric','thres'});
        
        T = load(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterData_thres-%d_T_p_clusterSize.mat',thres)));
        S = fullfile(wbDir,atlas,'fs_LR.164k.L.inflated.surf.gii');
        S2 = fullfile(wbDir,atlas,'fs_LR.164k.L.surfArea.func.gii');
        comm = sprintf('wb_command -surface-vertex-areas %s %s',S,S2);
        [err,out]=system(comm);
        PP = [];
        for ss=sessN
            permDir = fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('sess-%d',ss));
            for h=hemi
                S = gifti(S2); % for left hemisphere only - surface area
                for i=1:2 % positive_negative
                    T1 = getrow(T,T.sessN==ss & T.positive_negative==i);
                    if i==1
                        name = 'non-permuted';
                    else
                        name = 'non-permuted_reversed';
                    end
                    C = gifti(fullfile(permDir,sprintf('clusters_%s_thres%d.func.gii',name,thres)));
                    for c=1:size(T1.cluster,1)
                        P = getrow(T1,T1.cluster==c);
                        P.surfArea = sum(S.cdata(C.cdata==c));
                        PP = addstruct(PP,P);
                    end
                end
            end
        end
        save(fullfile(wbDir,atlas,sprintf('permutation_%s',metric),sprintf('clusterData_thres-%d_T_p_clusterSize.mat',thres)),'-struct','PP');

    case 'SEARCH_all'
        vararginoptions(varargin,{'sn','sessN'});
        if sessN == 1
        sml1_imana('SEARCH_define','sn',sn,'sessN',sessN);
        end
        sml1_imana('SEARCH_dist_runLDC','sn',sn,'sessN',sessN);
        sml1_imana('SEARCH_dist_map','sn',sn,'sessN',sessN);
        sml1_imana('SEARCH_dist_surf','sn',sn,'sessN',sessN);
    case 'SEARCH_define'                                                    % STEP 4.1a   :  Defines searchlights for 120 voxels in grey matter surface
        sn=4;
        sessN=1;
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
                cd(fullfile(glmSessDir{ss},subj_name{s}));
                if exist(sprintf('%s_sess%d_newTest_LDC.nii',subj_name{s},ss));
                    fprintf('Searchlight map: %s sess%d already exists - skipping\n',subj_name{s},ss);
                else
                    % load their searchlight definitions and SPM file
                    L = load(fullfile(anatomicalDir,subj_name{s},sprintf('%s_searchlight_120.mat',subj_name{s})));
                    load SPM;
                    SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                    
                    name = sprintf('%s_sess%d_newTest',subj_name{s},ss);
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
        for ss=sessN
            for s = sn
                % Load subject surface searchlight results (1 vol per paired conds)
                %LDC_file            =
                %fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_newTest_LDC.nii',subj_name{s},ss));
                %%NOTE! newTest is a newer version, done for session 2
                %LDC_file            = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_LDC.nii',subj_name{s},ss)); % searchlight nifti - use the recalculated version
                LDC_file            = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_newTest_LDC.nii',subj_name{s},ss)); % searchlight nifti - use the recalculated version
                %LDC_file            = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_LDC_recalc.nii',subj_name{s},ss)); % searchlight nifti - use the recalculated version
                %LDC_file            = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_LDC_recalc_test.nii',subj_name{s},ss)); % searchlight nifti - use the recalculated version
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
                Y.fname   = sprintf('%s_sess%d_dist.nii',subj_name{s},ss);
                
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
                T.fname = sprintf('%s_sess%d_dist_trained.nii',subj_name{s},ss);
                U.fname = sprintf('%s_sess%d_dist_untrained.nii',subj_name{s},ss);
                Z.fname = sprintf('%s_sess%d_dist_cross.nii',subj_name{s},ss);
                
                % save outputs
                spm_write_vol(Y,Y.LDC);
                spm_write_vol(T,T.LDC);
                spm_write_vol(U,U.LDC);
                spm_write_vol(Z,Z.LDC);
                fprintf('Done %s_sess%d \n',subj_name{s},ss)
                
                clear vol vdat LDC Y T U Z
                
            end
        end
        cd(cWD);  % return to working directory
    case 'SEARCH_dist_surf'                                                 % STEP 4.4a   :  Map searchlight results (.nii) onto surface (.metric)
        % map volume images to metric file and save them in individual surface folder
        sn      = [5:9,11:31];
        sessN   = [1:4];   
        fileList = {'dist','dist_trained','dist_untrained','dist_cross'}; % for dist or repsup
        glmDir = glmSessDir;
        outname = 'dist'; % dist or dist_repsup
        smooth = 1;
        vararginoptions(varargin,{'sn','sessN','fileList','glmDir','outname'});
        
        for ss = sessN
            for s = sn
                for h=1:2
                    caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                    white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
                    pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
                    
                    for f = 1:length(fileList)
                        images{f}    = fullfile(glmDir{ss},subj_name{s},sprintf('%s_sess%d_%s.nii',subj_name{s},ss,fileList{f}));
                        column_name{f} = fullfile(sprintf('Sess%d_%s.nii',ss,fileList{f}));
                    end;    % filename
                    outfile   = sprintf('%s_%s_sess%d.metric',subj_name{s},outname,ss);
                    M         = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
                    M.column_name = column_name;
                    caret_save(fullfile(caretSDir,outfile),M);
                    if smooth == 1;
                        % Smooth output .metric file (optional)
                        % Load .topo file
                        closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                        white = fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                        Out = caret_smooth(fullfile(caretSDir,outfile), 'coord',white,'topo',closed);%,...
                        %'algorithm','FWHM','fwhm',12);
                        char(Out);  % if smoothed adds an 's'
                    end;
                    fprintf('Done subj %d sess %d hemi %d \n',s,ss,h);
                end;    % hemi
            end;    % subj
        end; % sessN
    case 'SEARCH_group_make'                                                % depreciated STEP 4.5   :  Make group metric files by condensing subjec contrast metric files
        % Calculate group metric files from the searchlight results. 
        % Takes the 4 contrast results ('dist','dist_trained','dist_untrained','dist_cross') across
        % subjects and makes a group level metric file that contains each
        % subject's data for that contrast type.
        sessN=1;
        sn=[4:9,11:31];
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
        sn=[4:9,11:31];
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
            
           % X=caret_crosssection(border,fullfile(surfGroupDir,[hem{h} '.INFLATED.coord']));
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
    case 'SEARCH_recalc'
        % need to recalculate searchlights for s04-s28 (old prewhitening)
        sn=[4:9,11:31];
        vararginoptions(varargin,{'sn'});
         for s=sn
            for ss=1:4
                % load SPM
                glmDir  = fullfile(glmSessDir{ss},subj_name{s});
                cd(glmDir);
                load SPM;
                idx = mean(diag(SPM.xX.Bcov));
                %idx = mean(diag(SPM.xX.Bcov))*SPM.xX.trRV; % just tried,didn't work
                V = spm_vol(fullfile(sprintf('%s_sess%d_LDC.nii',subj_name{s},ss)));
                V2 = spm_read_vols(V);
             %   if s<29 % for new subjects no need
             for i=1:size(V2,4)
                 V3(:,:,:,i) = V2(:,:,:,i).*(k2./idx*SPM.xX.trRV);
             end
                    V3 = V2./idx;
             %   else
              %      V3 = V2;
              %  end
                newName = sprintf('%s_sess%d_LDC_recalc.nii',subj_name{s},ss);
                for i=1:size(V,1)
                    Z=NaN(V(1).dim);
                    %[d f t m]=spm_fileparts(outFiles{i});
                    outFiles{i}=fullfile(cd,sprintf('%s,%d',newName,i));
                    [d f t m]=spm_fileparts(outFiles{1});
                    Vo      = struct(...
                        'fname',    fullfile(d,[f t]),...
                        'dim',      V(1).dim,...
                        'dt',       [spm_type('float32') spm_platform('bigend')],...
                        'mat',      V(i).mat,...
                        'n',        [i 1],...
                        'descrip',  V(i).descrip);
                    Z=V3(:,:,:,i);            % Project back into voxel space
                    spm_write_vol(Vo,Z);
                end; % per distance pair
                fprintf('Done:\t %s   sess-%d\n',subj_name{s},ss);
            end; % session
        end; % subject
    case 'SEARCH_check'
        %check if searchlights exist
        sn = [5:9,11:31];
        sessN = 4;
        vararginoptions(varargin,{'sn','sessN'});
        for ss=sessN
            idx=0;
            for s=sn
                % go to subject's glm directory
                cd(fullfile(glmSessDir{ss},subj_name{s}));
                if exist(sprintf('%s_sess%d_newTest_LDC.nii',subj_name{s},ss));
                    fprintf('Searchlight map: %s sess%d already exists!\n',subj_name{s},ss);
                else
                    fprintf('Searchlight map MISSING: %s sess%d!!!!\n',subj_name{s},ss);
                    idx=idx+1;
                end
            end
            fprintf('\nMissing %d searchlights for sess-%d!!\n',idx,ss);
        end

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
    
    case 'SURF_wb:map_psc_individ'
        % projects individual percent signal change volume files to WorkBench surface
        sessN=1:4;
        sn=[5:9,11:31];
        name={'AllSeq','TrainSeq','UntrainSeq'};
        colName={'psc_all','psc_trained','psc_untrained'};
        vararginoptions(varargin,{'sn','sessN','name','colName'});

        for ss=sessN
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    white   = fullfile(subjDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},hemI{h}));
                    pial    = fullfile(subjDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},hemI{h}));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    for f=1:length(name)
                        images{f}=fullfile(glmSessDir{ss},subj_name{s},sprintf('psc_sess%d_%s.nii',ss,name{f}));
                        column_name{f} = fullfile(sprintf('Sess%d_%s.nii',ss,colName{f}));
                    end;
                    outfile         = fullfile(subjDir,sprintf('%s.%s.psc.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    G               = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names',column_name, ...
                                        'anatomicalStruct',hemname{h});
                    save(G, outfile);
                    fprintf('Sess: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session
    case 'SURF_wb:map_dist_individ'
        % projects individual distance files to WorkBench surface
        sessN=1:4;
        sn=[5:9,11:31];
        vararginoptions(varargin,{'sn','sessN'});
        name = {'dist','dist_trained','dist_untrained','dist_cross'};  
        OUTname = {'dist_all','dist_trained','dist_untrained','dist_cross'};  
        for ss=sessN
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    white   = fullfile(subjDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},hemI{h}));
                    pial    = fullfile(subjDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},hemI{h}));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    for f = 1:length(name)
                        images{f}    = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_%s.nii',subj_name{s},ss,name{f}));
                        column_name{f} = fullfile(sprintf('Sess%d_%s.nii',ss,OUTname{f}));
                    end;    
                    outfile         = fullfile(subjDir,sprintf('%s.%s.dist.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    G               = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names',column_name, ...
                                        'anatomicalStruct',hemname{h});
                    % filter out any crazy values (on the border of mask)
                    for c=1:size(G.cdata,2)
                        G.cdata(G.cdata(:,c)>5,c)=0;
                        G.cdata(G.cdata(:,c)<-5,c)=0;
                    end
                    save(G, outfile);
                    fprintf('Sess: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session
    case 'SURF_wb:map_seq_diff_individ'
        %calculate the difference in distances per individual (trained -
        %untrained)
        sessN=1:4;
        sn=[5:9,11:31];
        metric = 'dist'; % dist or psc
        name = {'dist_trained','dist_untrained'};
        vararginoptions(varargin,{'sn','sessN','name','metric'});
        for ss=sessN
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    white   = fullfile(subjDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},hemI{h}));
                    pial    = fullfile(subjDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},hemI{h}));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    switch metric
                        case 'dist'
                            INname1 = sprintf('%s_sess%d_%s.nii',subj_name{s},ss,name{1});
                            INname2 = sprintf('%s_sess%d_%s.nii',subj_name{s},ss,name{2});
                        case 'psc'
                            INname1 = sprintf('%s_sess%d_%s.nii',metric,ss,name{1});
                            INname2 = sprintf('%s_sess%d_%s.nii',metric,ss,name{2});
                    end
                    t       = fullfile(glmSessDir{ss},subj_name{s},INname1); % trained
                    ut      = fullfile(glmSessDir{ss},subj_name{s},INname2); % untrained
                    
                    outfile         = fullfile(subjDir,sprintf('%s.%s.%sDiff.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss));
                    G_trained       = surf_vol2surf(C1.vertices,C2.vertices,t,'column_names',{sprintf('%sDifference',metric)}, ...
                                        'anatomicalStruct',hemname{h});                     
                    G_untrained     = surf_vol2surf(C1.vertices,C2.vertices,ut,'column_names',{sprintf('%sDifference',metric)}, ...
                                        'anatomicalStruct',hemname{h});
                    G = G_trained;
                    G.cdata = G.cdata - G_untrained.cdata;
                    % filter out any crazy values (on the border of mask)
                    for c=1:size(G.cdata,2)
                        G.cdata(G.cdata(:,c)>5,c)=0;
                        G.cdata(G.cdata(:,c)<-5,c)=0;
                    end
                    save(G, outfile);
                    fprintf('Sess: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session
    case 'SURF_wb:map_dist_sess_diff_individ'
        %calculate the difference in distances across sessions
        %(1-2, 2-3, 3-4) for trained and control
        sessTr = 1:3;
        sn=[5:9,11:31];
        vararginoptions(varargin,{'sn','sessN'});
        name = {'dist','dist_trained','dist_untrained','dist_cross'};  
        OUTname = {'dist_all','dist_trained','dist_untrained','dist_cross'};  
        for ss=sessTr
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    white   = fullfile(subjDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},hemI{h}));
                    pial    = fullfile(subjDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},hemI{h}));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    for f = 1:length(name)
                        s1  = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_%s.nii',subj_name{s},ss,name{f}));
                        s2  = fullfile(glmSessDir{ss+1},subj_name{s},sprintf('%s_sess%d_%s.nii',subj_name{s},ss+1,name{f}));
                        column_name{f} = fullfile(sprintf('dist_sessTr_%d-%d_%s.nii',ss,ss+1,OUTname{f}));
                        G1  = surf_vol2surf(C1.vertices,C2.vertices,s1,'anatomicalStruct',hemname{h});
                        G2  = surf_vol2surf(C1.vertices,C2.vertices,s2,'anatomicalStruct',hemname{h});
                        data(:,f) = G2.cdata - G1.cdata;
                    end;      
                    outfile         = fullfile(subjDir,sprintf('%s.%s.dist_sessTr_%d-%d.func.gii',subj_name{s},hemI{h},ss,ss+1));
                    G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
                    % filter out any crazy values (on the border of mask)
                 %   for c=1:size(G.cdata,2)
                 %       G.cdata(G.cdata(:,c)>5,c)=0;
                  %      G.cdata(G.cdata(:,c)<-5,c)=0;
                  %  end
                    save(G, outfile);
                    fprintf('Sess transition: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session transition
    case 'SURF_wb:map_psc_sess_diff_individ'
         %calculate the difference in distances across sessions
        %(1-2, 2-3, 3-4) for trained and control
        sessTr = 1:3;
        sn=[5:9,11:31];
        vararginoptions(varargin,{'sn','sessTr'});
        name = {'psc_all','psc_trained','psc_untrained'};  
        nameSeq = {'AllSeq','TrainSeq','UntrainSeq'};
        for ss=sessTr
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    white   = fullfile(subjDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},hemI{h}));
                    pial    = fullfile(subjDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},hemI{h}));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    for f = 1:length(name)
                        s1  = fullfile(glmSessDir{ss},subj_name{s},sprintf('psc_sess%d_%s.nii',ss,nameSeq{f}));
                        s2  = fullfile(glmSessDir{ss+1},subj_name{s},sprintf('psc_sess%d_%s.nii',ss+1,nameSeq{f}));
                        column_name{f} = fullfile(sprintf('psc_sessTr_%d-%d_%s.nii',ss,ss+1,name{f}));
                        G1  = surf_vol2surf(C1.vertices,C2.vertices,s1,'anatomicalStruct',hemname{h});
                        G2  = surf_vol2surf(C1.vertices,C2.vertices,s2,'anatomicalStruct',hemname{h});
                        data(:,f) = G2.cdata - G1.cdata;
                    end;      
                    outfile         = fullfile(subjDir,sprintf('%s.%s.psc_sessTr_%d-%d.func.gii',subj_name{s},hemI{h},ss,ss+1));
                    G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
                    % filter out any crazy values (on the border of mask)
                    for c=1:size(G.cdata,2)
                        G.cdata(G.cdata(:,c)>10,c)=0;
                        G.cdata(G.cdata(:,c)<-10,c)=0;
                    end
                    save(G, outfile);
                    fprintf('Sess transition: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session transition
    case 'SURF_wb:map_psc_zscore_individ'
        % map psc as zscores
        sessN = 1:4;
        sn=[5:9,11:31];
        vararginoptions(varargin,{'sn','sessTr'});
        name = {'psc_trained','psc_untrained'};  
        for s=sn
            for h=1:2
                for ss=sessN
                    subjDir = fullfile(wbDir,subj_name{s});
                    G1  = gifti(fullfile(subjDir,sprintf('%s.%s.psc.sess-%d.func.gii',subj_name{s},hemI{h},ss)));
                    for f = 1:length(name)
                        % fist subtract the mean
                        X = bsxfun(@minus,G1.cdata(:,f),nanmean(G1.cdata(:,f),1));
                        % And make standard deviation equal to 1, per column:
                        data(:,f) = bsxfun(@rdivide,X,nanstd(X,[],1));
                        column_name{f} = sprintf('Sess%d_Z_%s',ss,name{f});
                    end;
                    outfile         = fullfile(subjDir,sprintf('%s.%s.Z_psc.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
                    % filter out any crazy values (on the border of mask)
                    for c=1:size(G.cdata,2)
                        G.cdata(G.cdata(:,c)>15,c)=0;
                        G.cdata(G.cdata(:,c)<-15,c)=0;
                    end
                    save(G, outfile);
                    fprintf('Session: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % session
            end; % hemi
        end; % subject
        
    case 'SURF_wb:map_psc_group'
        sessN=1:4;
        sn=[5:9,11:31];
        atlas = 'fs_LR_164'; % 164 or 42
        metric = 'psc'; % psc or Z_psc
        replaceNaN = 1;
        INname = {'psc_all','psc_trained','psc_untrained'}; % psc_trained or Z_psc_trained
        vararginoptions(varargin,{'sessN','sn','metric','INname','replaceNaN'});
        
        for ss=sessN
            % Loop over hemispheres.
            for h = 1:2
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = fullfile(wbDir,atlas);
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'INname' and make a group metric file
                for j=1:length(INname)
                    % Loop over subjects...
                    for i = 1:length(sn);
                        % ...and define the names of their metric files
                        infilenames{i} = fullfile(wbDir,subj_name{sn(i)},sprintf('%s.%s.%s.sess-%d.func.gii',subj_name{sn(i)},hemI{h},metric,ss));
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    outfilenames    = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                    summaryname     = fullfile(surfaceGroupDir,sprintf('%s.group.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                    surf_groupGiftis(infilenames,'outfilenames',{outfilenames},'inputcol',j,'groupsummary',summaryname,'replaceNaNs',replaceNaN);
                end
                fprintf('Done: %s - sess%d\n',hemI{h},ss);
            end;
        end;
    case 'SURF_wb:map_dist_group'
        sessN=1:4;
        sn=[5:9,11:31];
        atlas = 'fs_LR_164'; % 164 or 42
        vararginoptions(varargin,{'sessN','sn'});        
        INname    = {'dist_all','dist_trained','dist_untrained','dist_cross'};
        replaceNaN = 1;
        for ss=sessN
            % Loop over hemispheres.
            for h = 1:2
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = fullfile(wbDir,atlas);
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'INname' and make a group metric file
                for j=1:length(INname)
                    % Loop over subjects...
                    for i = 1:length(sn);
                        % ...and define the names of their metric files
                        infilenames{i} = fullfile(wbDir,subj_name{sn(i)},sprintf('%s.%s.dist.sess-%d.func.gii',subj_name{sn(i)},hemI{h},ss));
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    outfilenames    = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                    summaryname     = fullfile(surfaceGroupDir,sprintf('%s.group.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                    surf_groupGiftis(infilenames,'outfilenames',{outfilenames},'inputcol',j,'groupsummary',summaryname,'replaceNaNs',replaceNaN);     
                end
                fprintf('Done: %s - sess%d\n',hemI{h},ss);
            end;
        end;
    case 'SURF_wb:map_seqDiff_group' 
        sessN=1:4;
        sn=[5:9,11:31];
        atlas = 'fs_LR_164'; % 164 or 42
        metric = 'dist'; % dist or psc
        vararginoptions(varargin,{'sessN','sn','metric'});
        INname    = sprintf('%sDifference',metric);
        replaceNaN = 1;
        for ss=sessN
            % Loop over hemispheres.
            for h = 1:2
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = fullfile(wbDir,atlas);
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'INname' and make a group metric file
                % Loop over subjects...
                for i = 1:length(sn);
                    % ...and define the names of their metric files
                    infilenames{i} = fullfile(wbDir,subj_name{sn(i)},sprintf('%s.%s.%sDiff.sess-%d.func.gii',subj_name{sn(i)},hemI{h},metric,ss));
                    % Name the output filename for this group metric file in average surface folder
                end;
                outfilenames    = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname,ss));
                summaryname     = fullfile(surfaceGroupDir,sprintf('%s.group.%s.sess-%d.func.gii',hemI{h},INname,ss));
                surf_groupGiftis(infilenames,'outfilenames',{outfilenames},'inputcol',1,'groupsummary',summaryname,'replaceNaNs',replaceNaN);
                fprintf('Done: %s - sess%d\n',hemI{h},ss);
            end;
        end;
    case 'SURF_wb:map_sessTrans_group' 
        sessTr=1:3;
        sn=[5:9,11:31];
        atlas       = 'fs_LR_164'; % 164 or 42
        INname      = {'psc_all','psc_trained','psc_untrained'}; % psc or dist: dist, dist_trained, dist_untrained
        metric      = 'psc'; % psc or dist
        vararginoptions(varargin,{'sessN','sn','metric','INname'});              
        replaceNaN = 1;
        for ss=sessTr
            % Loop over hemispheres.
            for h = 1:2
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = fullfile(wbDir,atlas);
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'INname' and make a group metric file
                for j=1:length(INname)
                    % Loop over subjects...
                    for i = 1:length(sn);
                        % ...and define the names of their metric files
                        infilenames{i} = fullfile(wbDir,subj_name{sn(i)},sprintf('%s.%s.%s_sessTr_%d-%d.func.gii',subj_name{sn(i)},hemI{h},metric,ss,ss+1));
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    outfilenames    = fullfile(surfaceGroupDir,sprintf('%s.%s.sessTr_%d-%d.func.gii',hemI{h},INname{j},ss,ss+1));
                    summaryname     = fullfile(surfaceGroupDir,sprintf('%s.group.sessTr_%d-%d_%s.func.gii',hemI{h},ss,ss+1,INname{j}));
                    surf_groupGiftis(infilenames,'outfilenames',{outfilenames},'inputcol',j,'groupsummary',summaryname,'replaceNaNs',replaceNaN);     
                end
                fprintf('Done: %s - sessTr:%d-%d\n',hemI{h},ss,ss+1);
            end;
        end;
    case 'SURF_wb:map_dist_sessTrans_group'
        sessTr=1:3;
        sn=[5:9,11:31];
        atlas = 'fs_LR_164'; % 164 or 42
        vararginoptions(varargin,{'sessN','sn'});        
        INname    = {'dist_trained','dist_untrained'};
        replaceNaN = 1;
        for ss=sessTr
            % Loop over hemispheres.
            for h = 1:2
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = fullfile(wbDir,atlas);
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'INname' and make a group metric file
                for j=1:length(INname)
                    % Loop over subjects...
                    for i = 1:length(sn);
                        % ...and define the names of their metric files
                        infilenames{i} = fullfile(wbDir,subj_name{sn(i)},sprintf('%s.%s.dist_sessTr_%d-%d.func.gii',subj_name{sn(i)},hemI{h},ss,ss+1));
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    outfilenames    = fullfile(surfaceGroupDir,sprintf('%s.%s.sessTr_%d-%d.func.gii',hemI{h},INname{j},ss,ss+1));
                    summaryname     = fullfile(surfaceGroupDir,sprintf('%s.group.sessTr_%d-%d_%s.func.gii',hemI{h},ss,ss+1,INname{j}));
                    surf_groupGiftis(infilenames,'outfilenames',{outfilenames},'inputcol',j+1,'groupsummary',summaryname,'replaceNaNs',replaceNaN);     
                end
                fprintf('Done: %s - sessTr:%d-%d\n',hemI{h},ss,ss+1);
            end;
        end;
    case 'SURF_wb:cSPM_group'
        % calculate mean statistics (t-test)
        metric = 'psc'; % psc or dist
        sessN = 1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        vararginoptions(varargin,{'metric','sessN','atlas'});
        switch metric
            case 'psc'
                INname = {'psc_all','psc_trained','psc_untrained'};
            case 'dist'
                INname = {'dist_all','dist_trained','dist_untrained','dist_cross'};
        end
        for h=1:2
            indx=1;
            fprintf('\n%s: session:',hemI{h});
            for ss=sessN
                surfaceGroupDir = fullfile(wbDir,atlas);
                for n=1:length(INname)
                    inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname{n},ss));
                    G1 = gifti(inFile);
                    cSPM = surf_getcSPM('onesample_t','data',G1.cdata,'maskthreshold',0.7);
                    outFile = sprintf('%s.%s.sess-%d',hemI{h},INname{n},ss);
                    save(fullfile(surfaceGroupDir,sprintf('cSPM_%s.mat',outFile)),'cSPM');
                    data(:,indx)=cSPM.con(1).con; % mean
                    data(:,indx+length(sessN)*length(INname))=cSPM.con(1).Z; % T
                    column_name{indx}=['mean_' outFile];
                    column_name{indx+length(sessN)*length(INname)}=['T_' outFile];
                    indx=indx+1;
                end
                fprintf('%d.',ss);
            end
            % here save the new functional files (per hemisphere)
            G2 = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            summaryName = fullfile(surfaceGroupDir,sprintf('%s.summary_%s.func.gii',hemI{h},metric));
            save(G2,summaryName);
        end
    case 'SURF_wb:cSPM_diff_group'
        sessN = 1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric = 'dist'; % dist or psc
        vararginoptions(varargin,{'metric','sessN','atlas'});
        INname = sprintf('%sDifference',metric);
        for h=1:2
            indx=1;
            fprintf('\n%s: session:',hemI{h});
            for ss=sessN
                surfaceGroupDir = fullfile(wbDir,atlas);
                inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname,ss));
                G1 = gifti(inFile);
                cSPM = surf_getcSPM('onesample_t','data',G1.cdata,'maskthreshold',0.7);
                outFile = sprintf('%s.%s.sess-%d',hemI{h},INname,ss);
                save(fullfile(surfaceGroupDir,sprintf('cSPM_%s.mat',outFile)),'cSPM');
                data(:,indx)=cSPM.con(1).con; % mean
                data(:,indx+length(sessN))=cSPM.con(1).Z; % T
                column_name{indx}=['mean_' outFile];
                column_name{indx+length(sessN)}=['T_' outFile];
                indx=indx+1;
                
                fprintf('%d.',ss);
            end
            % here save the new functional files (per hemisphere)
            G2 = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            summaryName = fullfile(surfaceGroupDir,sprintf('%s.summary_%s.func.gii',hemI{h},INname));
            save(G2,summaryName);
        end
    case 'SURF_wb:cSPM_sessTrans_group'
        % calculate mean statistics (t-test)
        sessTr = 1:3;
        atlas = 'FS_LR_164'; % 164 or 42
        INname = {'dist_all','dist_trained','dist_untrained'}; % dist: dist, dist_trained, dist_untrained
        metric = 'dist'; % dist or psc
        useSmooth = 0;
        vararginoptions(varargin,{'metric','sessN','atlas','INname','useSmooth'});

        for h=1:2
            indx=1;
            fprintf('\n%s: session transition:',hemI{h});
            for ss=sessTr
                surfaceGroupDir = fullfile(wbDir,atlas);
                for n=1:length(INname)
                    if useSmooth
                        inFile = fullfile(surfaceGroupDir,sprintf('s%s.%s.sessTr_%d-%d.func.gii',hemI{h},INname{n},ss,ss+1));
                        outFile = sprintf('s%s.%s.sessTr_%d-%d',hemI{h},INname{n},ss,ss+1);
                    else
                        inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sessTr_%d-%d.func.gii',hemI{h},INname{n},ss,ss+1));
                        outFile = sprintf('%s.%s.sessTr_%d-%d',hemI{h},INname{n},ss,ss+1);
                    end
                    G1 = gifti(inFile);
                    cSPM = surf_getcSPM('onesample_t','data',G1.cdata,'maskthreshold',0.7);
                    save(fullfile(surfaceGroupDir,sprintf('cSPM_%s.mat',outFile)),'cSPM');
                    data(:,indx)=cSPM.con(1).con; % mean
                    data(:,indx+length(sessTr)*length(INname))=cSPM.con(1).Z; % T
                    column_name{indx}=['mean_' outFile];
                    column_name{indx+length(sessTr)*length(INname)}=['T_' outFile];
                    indx=indx+1;
                end
                fprintf('%d.',ss);
            end
            % here save the new functional files (per hemisphere)
            G2 = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            if useSmooth
                summaryName = fullfile(surfaceGroupDir,sprintf('s%s.summary_%s_sessTrans.func.gii',hemI{h},metric));
            else
                summaryName = fullfile(surfaceGroupDir,sprintf('%s.summary_%s_sessTrans.func.gii',hemI{h},metric));
            end
            save(G2,summaryName);
        end
    case 'SURF_wb:zscore'
        % zscore psc / dist maps
        metric = 'psc'; % psc or dist
        sessN = 1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        vararginoptions(varargin,{'metric','sessN','atlas'});
        switch metric
            case 'psc'
                INname = {'psc_trained','psc_untrained'};
            case 'dist'
                INname = {'dist_all','dist_trained','dist_untrained','dist_cross'};
        end
        for h=1:2
            indx=1;
            fprintf('\n%s: session:',hemI{h});
            for ss=sessN
                surfaceGroupDir = fullfile(wbDir,atlas);
                for n=1:length(INname)
                    inFile = fullfile(surfaceGroupDir,sprintf('%s.group.%s.sess-%d.func.gii',hemI{h},INname{n},ss));
                    G1 = gifti(inFile);
                    % fist subtract the mean
                    X = bsxfun(@minus,G1.cdata,nanmean(G1.cdata,1));
                    % And make standard deviation equal to 1, per column:
                    data(:,indx) = bsxfun(@rdivide,X,nanstd(X,[],1));
                    outName = sprintf('%s.%s.sess-%d',hemI{h},INname{n},ss);
                    column_name{indx}=['Z_' outName];
                    indx=indx+1;
                end
                fprintf('%d.',ss);
            end
            % here save the new functional files (per hemisphere)
            G2 = surf_makeLabelGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            %G2 = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            summaryName = fullfile(surfaceGroupDir,sprintf('%s.Zscore_%s.func.gii',hemI{h},metric));
            save(G2,summaryName);
        end
    case 'SURF_wb:zscore_diff'
        metric = 'psc'; % psc or dist
        sessN = 1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        vararginoptions(varargin,{'metric','sessN','atlas'});
        switch metric
            case 'psc'
                col = [1,2]; mu = 2;
            case 'dist'
                col = [2,3]; mu = 4;
        end
        surfaceGroupDir = fullfile(wbDir,atlas);
        
        for h=1:2
            indx=1;
            fprintf('\n%s: session:',hemI{h});
            inFile = fullfile(surfaceGroupDir,sprintf('%s.Zscore_%s.func.gii',hemI{h},metric));
            G1 = gifti(inFile);
            for ss=sessN
                data(:,indx) = G1.cdata(:,(ss-1)*mu+col(1))-G1.cdata(:,(ss-1)*mu+col(2));
                column_name{indx}=[sprintf('Z_difference_sess%d_%s',ss,metric)];
                indx=indx+1;
                fprintf('%d.',ss);
            end
            % here save the new functional files (per hemisphere)
            G2 = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            summaryName = fullfile(surfaceGroupDir,sprintf('%s.ZscoreDifference_%s.func.gii',hemI{h},metric));
            save(G2,summaryName);
        end
    case 'SURF_wb:smooth_individ'
        %metric      = 'psc'; % psc, dist, distDifference, dist_sessTrans
        metric      = {'sessTr_1-2','sessTr_2-3','sessTr_3-4'};
        atlas       = 'FS_LR_164'; % 164 or 42
        INname      = {'psc_all','psc_trained','psc_untrained'}; % psc or dist
        kernel      = 2;
        vararginoptions(varargin,{'metric','sessN','kernel','INname'});
        
        % for ss=sessN
        surfaceGroupDir = fullfile(wbDir,atlas);
        for f=1:length(INname)
            for h=1:2
                for m=1:length(metric)
                    %inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},metric,ss));
                    inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.%s.func.gii',hemI{h},INname{f},metric{m}));
                    surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
                    surf_smooth(inFile,'surf',surfFile,'kernel',kernel);
                end
            end
        end
    case 'SURF_wb:smooth_group'
        metric = 'summary_psc'; % psc, dist, distDifference, dist_sessTrans
        atlas = 'FS_LR_164'; % 164 or 42
        kernel = 2;
        vararginoptions(varargin,{'metric','sessN','kernel'});
        
        % for ss=sessN
        for h=1:2
            surfaceGroupDir = fullfile(wbDir,atlas);
            %inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},metric,ss));
            inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.func.gii',hemI{h},metric));
            surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
            surf_smooth(inFile,'surf',surfFile,'kernel',kernel);
        end
        %  end
    case 'SURF_wb:map_group_diff'
        sessN=1:4;
        atlas = 'FS_LR_164'; % 164 or 42
        metric = 'dist';
        vararginoptions(varargin,{'metric'});
        switch metric
            case 'dist'
                m=4;
                a=2;
            case 'psc'
                m=2;
                a=1;
        end
        % Loop over hemispheres.
        for h = 1:2
            for ss=sessN
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = fullfile(wbDir,atlas);
                inFile     = fullfile(surfaceGroupDir,sprintf('%s.summary_%s.func.gii',hemI{h},metric));
                I = gifti(inFile);
                Data(:,ss) = I.cdata(:,(ss-1)*m+a)-I.cdata(:,(ss-1)*m+a+1); % trained - untrained
                column_name{ss}=sprintf('%s_sess-%d',metric,ss);
            end;
            if strcmp(metric,'dist') % just so it can be modulated in workbench (view); otherwise numbers too small
                Data = Data.*100;
            end
            G = surf_makeFuncGifti(Data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            outFile = fullfile(surfaceGroupDir,sprintf('%s.group.%sDiff.func.gii',hemI{h},metric));
            save(G,outFile);
            % smooth in addition
            surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
            surf_smooth(outFile,'surf',surfFile,'kernel',2);
        end;
    case 'SURF_wb:map_group_transDiff' % this in use for session difference / transition !!!
        atlas = 'FS_LR_164'; % 164 or 42
        metric = 'dist';
        seqType = {'all','trained','untrained'};
        vararginoptions(varargin,{'metric'});
        switch metric
            case 'dist'
                m=4;
                a=1;
            case 'psc'
                m=2;
                a=1;
        end
        % Loop over hemispheres.
        % Go to the directory where the group surface atlas resides
        surfaceGroupDir = fullfile(wbDir,atlas);
        for h = 1:2
            inFile     = fullfile(surfaceGroupDir,sprintf('%s.summary_%s.func.gii',hemI{h},metric));
            I = gifti(inFile);
            indx=1;
            for ss=1:3 % transition
                for st =1:3
                    % dist
                    Data(:,indx) = I.cdata(:,(ss)*m+a+(st-1))-I.cdata(:,(ss-1)*m+a+(st-1)); 
                    column_name{indx}=sprintf('%s_sess-%d-%d_seqType-%s',metric,ss,ss+1,seqType{st});
                    % T_values
                    if strcmp(metric,'dist')
                        Data(:,indx+1) = I.cdata(:,(ss)*m+a+(st-1)+16)-I.cdata(:,(ss-1)*m+a+(st-1)+16); 
                        column_name{indx+1}=sprintf('T-value_%s_sess-%d-%d_seqType-%s',metric,ss,ss+1,seqType{st});
                        indx = indx+2;
                    else
                        indx=indx+1;
                    end 
                end
            end;
          %  if strcmp(metric,'dist') % just so it can be modulated in workbench (view); otherwise numbers too small
          %      Data = Data.*100;
           % end
            G = surf_makeFuncGifti(Data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            outFile = fullfile(surfaceGroupDir,sprintf('%s.group.%sDiff_transitions_groupVersion.func.gii',hemI{h},metric));
            save(G,outFile);
            % smooth in addition
            surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
            surf_smooth(outFile,'surf',surfFile,'kernel',1);
        end;
    case 'SURF_wb:make_distMask' %old
        % make a mask based on group T-value of distances
        mask_thres = 1.71; % mask threshold (here so that t-map is p<.5)
        atlas = 'FS_LR_164'; % 164 or 42
        for h=1:2
            % first load the surface
            surfaceGroupDir = fullfile(wbDir,atlas);
            surfFile = fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
            inFile = fullfile(surfaceGroupDir,sprintf('%s.summary_dist.func.gii',hemI{h}));
            G = gifti(inFile);
            for ss=1:4 % across regions
                DD(:,ss)=G.cdata(:,(4*4)+(ss-1)*4+1); % get the correct t-map
            end
            maxDist(:,1)    = max(DD,[],2); % here the maximum distance 
            % first save the maxDist mask
            outFile         = fullfile(surfaceGroupDir,sprintf('%s.maskDist.func.gii',hemI{h}));
            G = surf_makeFuncGifti(maxDist,'anatomicalStruct',hemname{h},'columnNames',{'maxDist_acrossSess'});
            save(G,outFile);
            % now remove clusters + threshold
            outFile2 = fullfile(surfaceGroupDir,sprintf('%s.maskDist_clusters.func.gii',hemI{h}));
            comm=sprintf('wb_command -metric-find-clusters %s %s %f %f %s',...
                surfFile,outFile,mask_thres,7,outFile2);
            fprintf('%s\n',comm)
            [err,out]=system(comm);
            % last smooth
            surf_smooth(outFile2,'surf',surfFile,'kernel',1);
            mask = gifti(fullfile(surfaceGroupDir,sprintf('s%s.maskDist_clusters.func.gii',hemI{h})));
            save(fullfile(surfaceGroupDir,sprintf('%s.maskDist.mat',hemI{h})),'mask'); % save binary values as a mask (for roi definitions)
        end
    case 'SURF_wb:defineTessels' %old
        nNode = 642; % 42, 162, 362, 642, 1002, 1442
        thres = 0.5; % how much of the tessel needs to have significant distance
        atlas = 'FS_LR_164'; % 164 or 42
        vararginoptions(varargin,{'nNode','thres'});
        tessel = cell(2,1);
        for h=1:2
             surfaceGroupDir = fullfile(wbDir,atlas);
             M=load(fullfile(surfaceGroupDir,sprintf('%s.maskDist.mat',hemI{h})));
             T = gifti(fullfile(surfaceGroupDir,sprintf('fs_LR.164k.%s.icosahedron-%d.label.gii',hemI{h},nNode)));
           %  tessel{h} = unique(M.mask.cdata.*double(T.cdata));
             countAll=zeros(nNode,1); countMask=countAll;
             for i=1:numel(unique(T.cdata))
                 countAll(i)=sum(T.cdata==i);
                 countMask(i)=sum(double(T.cdata).*M.mask.cdata==i);
             end
             prop=countMask./countAll;
             tessel{h}=find(prop>thres);
        end
        varargout{1} = tessel; % per hemisphere region definitions
    case 'SURF_wb:speed_mask'
        % here mask the speed surface files 
        % just for sessions 3 and 4
        % used in the paper
        metric = 'summary_distDifference'; % summary_distDifference or group.distDiff_transitions_groupVersion
        atlas = 'FS_LR_164'; % 164 or 42
        Tcut  = 1.7; % what T-stat cut-off to use
        vararginoptions(varargin,{'metric','Tcut'});
        switch metric
            case 'summary_distDifference'
                INCol = [3,4];
                colName = {'sess3_distDifference','sess4_distDifference'};
            case 'group.distDiff_transitions_groupVersion'
                INCol = [13,15,17];
                colName = {'all_sess3-4','trained_sess3-4','untrained_sess3-4'};
        end
        maskCol = [27,26,27]; % third session
        for h=1:2
            surfaceGroupDir = fullfile(wbDir,atlas);
            inFile = fullfile(surfaceGroupDir,sprintf('s%s.%s.func.gii',hemI{h},metric));
            maskFile = fullfile(surfaceGroupDir,sprintf('s%s.summary_dist.func.gii',hemI{h}));
            M = gifti(maskFile);
            IN = gifti(inFile);
            
            maskIdx = sum(M.cdata(:,maskCol)>Tcut,2)>0;
            for i=1:length(INCol);
                data(:,i) = IN.cdata(:,INCol(i)).*maskIdx;
                column_name{i} = colName{i};
            end
            outfile         = fullfile(surfaceGroupDir,sprintf('s%s.%s_masked.func.gii',hemI{h},metric));
            G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            save(G, outfile);
        end
        
        
    case 'BETA_get'                                                     % STEP 5.6   :  Harvest betas from rois (raw, univ, multiv prewhit)    
        sessN = 1:4;
        sn  = [5:9,11:31];    
        %roi = 1:16; 
        roi = [1:3,7,8];
        roiDefine = 'subset'; % determine all regions from region file R, or 'tesselSelect' - based on distance mask
        parcelType = 'Brodmann'; % 162tessels, Brodmann, cortex_buckner, BG-striatum, thalamus, tesselsWB_162, tesselWB_362
        numError = 1;
        vararginoptions(varargin,{'sn','sessN','roi','parcelType','roiDefine','numError'});
       
        if strcmp(parcelType,'Yokoi_clusters')
            regType = [1 2 3 4 5 6 7 8 9 10 1 2 3 6 8 9 10];
            regSide = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2];
        end
        idxRoi=0;
        for ss=sessN
            % harvest
            for s=sn % for each subj
                T=[];
                fprintf('\nSubject: %d\n',s) % output to user
                % here first check if the beta file already exists - run
                % only if not
                %if exist(fullfile(betaDir,subj_name{s},sprintf('betas_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'file')
                %    fprintf('Beta file already exists....skipping.\n');
                %else
                    fprintf('Starting to extract betas...\n');
                    tElapsed=tic;
                    % load files
                    load(fullfile(glmErrorDir{numError}{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                    %load(fullfile(glmSessDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                    load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)]));  % load subject's region parcellation (R)
                    
                    if idxRoi==0 % determine for first subject
                        if strcmp(roiDefine,'all')
                            roi=1:size(R,2);
                            nReg_hemi = size(R,2)/2;
                        elseif strcmp(roiDefine,'tesselSelect')
                            nReg_hemi = size(R,2)/2;
                            roi = [];
                            for h=1:2
                                tessels{h}=sml_connect('TESSEL:select','hemi',h,'nTessel',642);
                                roi = [roi nReg_hemi*(h-1)+tessels{h}];
                            end
                        elseif strcmp(roiDefine,'subset')
                            roi = roi;
                            nReg_hemi = max(roi)+1;
                        end
                        idxRoi = 1;
                    end
                    % this used to ensure the tessel / region allocation is
                    % correct
%                     for i=roi
%                         if i<nReg_hemi
%                             hemT = tessels{1};
%                             regI = i;
%                         else
%                             hemT = tessels{2};
%                             regI = i-nReg_hemi;
%                         end
%                         t = R{i}.name;       
%                         idx=t(10:end)==num2str(regI);
%                         if sum(idx==0)>0
%                             keyboard;
%                         end
%                     end
                    %cd (fullfile(glmSessDir{ss},subj_name{s})); % maybe change back when remove glm           
                    cd (fullfile(glmErrorDir{numError}{ss},subj_name{s})); % maybe change back when remove glm           
                    %P=SPM.Vbeta(SPM.xX.iC);
                 %   P=SPM.Vbeta;
                    % Add a few extra images
                    %----task against rest
                    O{1}=sprintf('psc_sess%d_TrainSeq.nii',ss); %psc trained
                    O{2}=sprintf('psc_sess%d_UntrainSeq.nii',ss); %psc untrained
                   % if any(strcmp(parcelType,{'tesselsWB_1062','tesselsWB_642'}))
                   %     O{3}=sprintf('%s_sess%d_dist.nii',subj_name{s},ss); %dist trained
                   %     O{4}=sprintf('%s_sess%d_dist_trained.nii',subj_name{s},ss); %dist trained
                   %     O{5}=sprintf('%s_sess%d_dist_untrained.nii',subj_name{s},ss); %dist untrained
                   %     O{6}=sprintf('%s_sess%d_dist_cross.nii',subj_name{s},ss); %dist cross
                   % end
                    oP=spm_vol(char(O));
                    
                    V = SPM.xY.VY;
                    % new
                    TI = load('SPM_info.mat');
                    idx_noErr = [TI.isError==0;ones(numel(unique(TI.run)),1)];
                    for r = roi % for each region
                        % get raw data for voxels in region
                        % determine if any voxels for that parcel
                        idx=[];
                        if size(R{r}.data,1)==0 % no voxels
                            % make data into NaN
                            S.betaW             = {NaN};
                            S.betaUW            = {NaN};
                            S.betaRAW           = {NaN};
                            S.resMS             = {NaN};
                            S.psc_train         = NaN;
                            S.psc_untrain       = NaN;
                           % if any(strcmp(parcelType,{'tesselsWB_1062','tesselsWB_642'}))
                           %     S.dist_all      = NaN;
                           %     S.dist_train    = NaN;
                           %     S.dist_untrain  = NaN;
                           %     S.dist_cross    = NaN;
                           % end
                        else
                          %  if strcmp(parcelType,'BG-striatumbla')
                                for ss1=1:4
                                    %load(fullfile(glmSessDir{ss1}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                                    load(fullfile(glmErrorDir{numError}{ss1}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                                    V_test = SPM.xY.VY;
                                    % different idea, just load psc map,
                                    % not sure it works
                                    %V_test = spm_vol(fullfile(glmErrorDir{numError}{ss1},subj_name{s},sprintf('psc_sess%d_TrainSeq.nii',ss1)));
                                    dataTest{ss1} = region_getdata(V_test,R{r});
                                end
                                i1 = find(dataTest{1}(1,:)); i2 = find(dataTest{2}(1,:)); i3 = find(dataTest{3}(1,:)); i4 = find(dataTest{4}(1,:));  
                                iA = intersect(i1,i2); iA = intersect(iA,i3); idx = intersect(iA,i4);
                         %   end
                            
                            Y = region_getdata(V,R{r});  % Data Y is N x P
                            data = region_getdata(oP,R{r}); % from added images
                       %     betaRaw = region_getdata(P,R{r});
                            % exclude any missing data in voxels
                           % idx = find(Y(1,:));
                            % estimate region betas
                            [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y(:,idx),SPM,'normmode','runwise');
                            %[betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y(:,idx),SPM,'normmode','overall');
                            %S.betaW                   = {betaW};                             % multivariate pw
                            %S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                            %S.betaRAW                 = {beta};
                            S.betaW     = {betaW(idx_noErr==1,:)};
                            S.betaUW    = {bsxfun(@rdivide,beta(idx_noErr==1,:),sqrt(resMS))};
                            S.betaRAW   = {beta(idx_noErr==1,:)};
                         %   S.betaRAW2  = {betaRaw(idx_noErr==1,:)};
                        %    S.betaUW2   = {bsxfun(@rdivide,betaRaw(idx_noErr==1,:),sqrt(resMS))};
                            S.resMS     = {resMS};
                            % info from maps for surface
                            S.psc_train         = {data(1,idx)};
                            S.psc_untrain       = {data(2,idx)};
                       %     if any(strcmp(parcelType,{'tesselsWB_1062','tesselsWB_642'}))
                       %         S.dist_all      = {data(3,idx)};
                       %         S.dist_train    = {data(4,idx)};
                       %         S.dist_untrain  = {data(5,idx)};
                       %         S.dist_cross    = {data(6,idx)};
                       %     end
                            clear Y data betaW beta resMS;
                        end
                        % voxel position
                        S.volcoord = {R{r}.data(idx,:)'};
                        S.SN = s;
                        S.region = r;
                        %if r<(numel(roi)/2)+1
                        if r<nReg_hemi
                            S.regSide = 1;
                            S.regType = S.region;
                        else
                            S.regSide = 2;
                            %S.regType = S.region-(numel(roi)/2);
                            S.regType = S.region-nReg_hemi;
                        end
                        if any(strcmp(parcelType,{'BG-striatum','thalamus'}))
                            S.regName = {R{r}.name};
                        elseif strcmp(parcelType,{'Yokoi_clusters'})
                            S.regSide = regSide(r);
                            S.regType = regType(r); 
                        end
                        T = addstruct(T,S);
                        fprintf('%d.',r);
                        clear idx;
                        %fprintf('elapsed %d\n',telapsed);
                    end
                    dircheck(fullfile(betaDir,subj_name{s}));
                    %save(fullfile(betaDir,subj_name{s},sprintf('betas_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'-struct','T');
                    save(fullfile(betaDir,subj_name{s},sprintf('betasError%d_%s_%s_sess%d.mat',numError,parcelType,subj_name{s},ss)),'-struct','T');
                    fprintf('\nDone beta extraction for sess%d-%s\n',ss,subj_name{s}); toc(tElapsed);
                    clear idx_noErr;
             %   end
            end
        end
    case 'BETA_voxels_check'
        sn = [5:9,11:31];
        parcelType = 'Brodmann';
        sessN = 1:4;
        numError = 1;
        problSubj = [];
        for s=sn
            voxNum = zeros(4,5);
            for ss=sessN
                T = load(fullfile(betaDir,subj_name{s},sprintf('betasError%d_%s_%s_sess%d.mat',numError,parcelType,subj_name{s},ss)));
                for r=1:5
                    voxNum(ss,r) = size(T.betaW{r},2);
                end
            end
            if numel(unique(voxNum))>5
                problSubj = [problSubj s];
            end
        end
        keyboard;
    case 'BETA_combineGroup'
        % combine individual subject beta structures into the whole
        % structure
        sessN=1:4;
        sn=[5:9,11:31];
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        numError = 1;
        parcelType='Brodmann'; % 162tessels, Brodmann, cortex_buckner, BG-striatum, thalamus, tesselsWB_162, tesselsWB_362
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type','parcelType','numError'});
        
        for ss=sessN
            switch(type)
                case 'new'
                    T=[];
                case 'add'
                    T=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
                  %  T=getrow(T,T.SN<29);
            end
            fprintf('subjects added for sess-%d:\n',ss);
            for s=sn
                %S=load(fullfile(betaDir,subj_name{s},sprintf('betas_%s_%s_sess%d',parcelType,subj_name{s},ss)));
                S=load(fullfile(betaDir,subj_name{s},sprintf('betasError%d_%s_%s_sess%d',numError,parcelType,subj_name{s},ss)));
                T=addstruct(T,S);
                
                fprintf('%d.',s);
            end
            %save(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)),'-struct','T','-v7.3');
            save(fullfile(betaDir,'group',sprintf('betasError%d_%s_sess%d.mat',numError,parcelType,ss)),'-struct','T','-v7.3');
        end
    case 'BETA_check'
        %check if searchlights exist
        sn = [5:9,11:31];
        sessN = 4;
        parcelType = 'tesselsWB_362';
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        for ss=sessN
            idx=0;
            for s=sn
                % go to subject's glm directory
                 if exist(fullfile(betaDir,subj_name{s},sprintf('betas_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'file')
                    fprintf('Beta file: %s sess%d already exists!\n',subj_name{s},ss);
                else
                    fprintf('Beta file MISSING: %s sess%d!!!!\n',subj_name{s},ss);
                    idx=idx+1;
                end
            end
            fprintf('\nMissing %d beta files for sess-%d!!\n',idx,ss);
        end
    case 'BETA_splitHalf'
          % crossvalidated version - separately for even and odd runs
        sn=[5:9,11:31];
        sessN=1:4;
        parcelType='162tessels'; % Brodmann or 162tessels, tesselsWB_162, tesselsWB_362
        betaChoice='multi';
        vararginoptions(varargin,{'sn','sessN','parcelType','betaChoice'});
       
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
            regSide(roi>158)=2;
            regType = roi;
            regType(regSide==2)=regType(regSide==2)-158;
        elseif strcmp(parcelType,'tesselsWB_162') || strcmp(parcelType,'tesselsWB_362')
            roi = sml1_imana_dist('CLUSTER_choose','sessN',sessN,'parcelType',parcelType)';
            T1=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,1)));
            T1 = getrow(T1,T1.SN==5);
            % choose only clusters where group dist >0
            regSide = T1.regSide(roi);
            regType = T1.regType(roi);
            %regSide = ones(size(roi));
           % regSide(roi>158)=2;
            %regType = roi;
            %regType(regSide==2)=regType(regSide==2)-158;
        elseif strcmp(parcelType,'tesselsWB_642')
            BB=load(fullfile(betaDir,'s05','betas_tesselsWB_642_s05_sess1')); % used for 642
            roi = BB.region;
            regSide = BB.regSide;
            regType = BB.regType;
        end

        for ss=sessN
            TT = [];
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
                        T.G_partA         = ones(1,12*12)*NaN;
                        T.G_partB         = ones(1,12*12)*NaN;
                    else
                        distA = rsa_squareRDM(rsa.distanceLDC(beta(idxA,:),pVec(idxA),cVec));
                        distB = rsa_squareRDM(rsa.distanceLDC(beta(idxB,:),pVec(idxB),cVec));
                        T.partA_train     = rsa_vectorizeRDM(distA(1:6,1:6));
                        T.partA_untrain   = rsa_vectorizeRDM(distA(7:12,7:12));
                        T.partA_all       = rsa_vectorizeRDM(distA);
                        T.partB_train     = rsa_vectorizeRDM(distB(1:6,1:6));
                        T.partB_untrain   = rsa_vectorizeRDM(distB(7:12,7:12));
                        T.partB_all       = rsa_vectorizeRDM(distB);
                        T.G_partA         = rsa_vectorizeIPMfull(pcm_estGCrossval(beta(idxA,:),pVec(idxA),cVec));      
                        T.G_partB         = rsa_vectorizeIPMfull(pcm_estGCrossval(beta(idxB,:),pVec(idxB),cVec));
                    end   
                    T.sn = s;
                    T.sessN = ss;
                    T.roi = roi(r);
                    T.regType = regType(r); % refer to D
                    T.regSide = regSide(r);
                    fprintf('%d.',r);
                    TT=addstruct(TT,T);
                end
            end
            save(fullfile(betaDir,'group',sprintf('betas_partition_%s_%sPW_sess%d',parcelType,betaChoice,ss)),'-struct','TT');
        end
    case 'BETA_splithalf_combine_crossSess'
        % here combine beta structures across sessions
        sessN = 1:4;
        parcelType = 'tesselsWB_162';
        TT = [];
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('betas_partition_%s_sess%d.mat',parcelType,ss)));
            keyboard;
            TT = addstruct(TT,T);
        end
        save(fullfile(betaDir,'group',sprintf('betas_partition_%s.mat',parcelType)));
    case 'BETA_combine_splitHalf'
        % combine different split-half beta structures
        % Brodmann, BG-striatum, thalamus
        sessN=[1:4];
        betaType='multi';
        parcelType = {'Brodmann','BG-striatum','thalamus'};
        
        TT=[];
        for ss=sessN
            idx = [0 0];
            for i=1:length(parcelType)
                T = load(fullfile(betaDir,'group',sprintf('betas_partition_%sPW_%s.mat',betaType,parcelType{i})));
                T = getrow(T,T.sessN==ss);
                T.regStruct = ones(size(T.sn))*i;
                T.roi       = T.roi + idx(1);
                T.regType   = T.regType + idx(2);
                idx = [(max(T.roi)) (max(T.regType))];
                
                TT = addstruct(TT,T);
            end
        end
    save(fullfile(betaDir,'group',sprintf('betas_partition_%sPW_combined',betaType)),'-struct','TT');
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
                open s
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
        sessN = 1:4;
        sn=[5:9,11:31];
        roi = 1:16;
        roiDefine = 'all'; % determine from region file
        betaChoice = 'multi'; % uni, multi or raw
        checkG=0; % test across all regions / subjects 
        simulations=0;
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        parcelType = 'Brodmann'; % or Brodmann, tesselsWB_162
        numError = 2;
        excludeError = 1;
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','checkG','simulations','type','parcelType','roiDefine','numError','excludeError'});
        
        for ss=sessN            
            if excludeError
                T = load(fullfile(betaDir,'group',sprintf('betasError%d_%s_sess%d.mat',numError,parcelType,ss))); % loads region data (T)
            else
                T = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            end
            
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
                D = load(fullfile(glmErrorDir{numError}{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                D = getrow(D,D.isError==0);
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
                        So.IPMfull = ones(1,144)*NaN;
                        So.psc_train = NaN;
                        So.psc_untrain = NaN;
                        Do.dist_train=NaN;
                        Do.dist_untrain=NaN;
                        Do.dist_cross=NaN;
                        Do.dist_all=NaN;
                        Do.eigTrain=ones(1,6)*NaN;
                        Do.eigUntrain=ones(1,6)*NaN;
                        Do.RDM_train=ones(1,15)*NaN;
                        Do.RDM_untrain=ones(1,15)*NaN;
                    else % run distance / RDM analyses
                        % crossval second moment matrix
                        [G,Sig]     = pcm_estGCrossval(betaW(1:size(D.run,1),:),D.run,D.seqNumb);
                        So.IPM      = rsa_vectorizeIPM(G);
                        So.Sig      = rsa_vectorizeIPM(Sig);
                        So.IPMfull  = rsa_vectorizeIPMfull(G);
                        % squared distances
                       % So.RDM      = rsa.distanceLDC(betaW(1:size(D.run,1),:),D.run,D.seqNumb);
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
                    So.regSide  = S.regSide;
                    So.regType  = S.regType;
%                     if r<(numel(roi)/2)+1
%                         So.regSide = 1;
%                         So.regType = So.region;
%                     else
%                         So.regSide = 2;
%                         So.regType = So.region-(numel(roi)/2);
%                     end
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
            
            if excludeError
                save(fullfile(betaDir,'group',sprintf('statsError%d_%s_%sPW_sess%d.mat',numError,parcelType,betaChoice,ss)),'-struct','To');
            else
                save(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)),'-struct','To');
            end
            
            if simulations==1
                save(fullfile(regDir,sprintf('stats_%s_%sPW_SIMULATIONS_sess%d.mat',parcelType,betaChoice,ss)),'-struct','Po');
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
            parcelType = 'Brodmann';
            vararginoptions(varargin,{'sessN','roi','betaChoice'});
            
            % load data and simulations
            D = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,sessN)));
            S = load(fullfile(regDir,sprintf('stats_%s_%sPW_SIMULATIONS_sess%d.mat',parcelType,betaChoice,sessN)));
            
            reg_indx=unique(D.region);
            
            for r=1:length(reg_indx)
              figure
                %  figure(1);
              %  subplot(1,length(reg_indx),r);
                histplot(ssqrt(S.dist_cross),'subset',S.region==r); % simulations
                drawline(mean(ssqrt(D.dist_cross(D.region==r))),'dir','vert','color',[1 0 0]);
                title(sprintf('%s',regname_cortex{r}));
               % if r==1
                    ylabel('Simulation count'); xlabel('Trained vs. Untrained dist');
               % end
                
%                 figure(2);
%                 subplot(1,length(reg_indx),r);
%                 histplot(S.dist_train,'subset',S.region==r); % simulations
%                 drawline(mean(D.dist_train(D.region==r)),'dir','vert','color',[1 0 0]);
%                 title(sprintf('%s',regname{r}));
%                 if r==1
%                     ylabel('Simulation count'); xlabel('Trained dist');
%                 end
%                 
%                 figure(3);
%                 subplot(1,length(reg_indx),r);
%                 histplot(S.dist_untrain,'subset',S.region==r); % simulations
%                 drawline(mean(D.dist_untrain(D.region==r)),'dir','vert','color',[1 0 0]);
%                 title(sprintf('%s',regname{r}));
%                 if r==1
%                     ylabel('Simulation count'); xlabel('Untrained dist');
%                 end
                
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
            
    case 'SURF_flatmap_rgb'
        sessN=[1:4];
        smooth=0;
        vararginoptions(varargin,{'sessN','smooth'});
        
        for ss=sessN
            for h=1:2  
                surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
                cd(surfaceGroupDir)
                if smooth == 1
                    Cp=caret_load(['s' hem{h} sprintf('.summary_psc_sess%d.metric',ss)]);
                    Cd=caret_load([hem{h} sprintf('.summary_dist_sess%d.metric',ss)]);
                else
                    Cp=caret_load([hem{h} sprintf('.summary_psc_sess%d.metric',ss)]);
                    Cd=caret_load([hem{h} sprintf('.summary_dist_sess%d.metric',ss)]);
                end
                flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
                topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
                
                %----PSC
                low_1=0.9;
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
                
                name={sprintf('psc_sess%d',ss),sprintf('dist_sess%d',ss)};
                %name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN),sprintf('dist_cross_sess%d',sessN)};
                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1,sc2},'column_name',name);
                
                if smooth == 1
                    caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.sess%d.RGB_paint',ss)],C);
                else
                    caret_save([surfaceGroupDir filesep hem{h} sprintf('.sess%d.RGB_paint',ss)],C);
                end
            end;
        end
    case 'SURF_flatmap_trainedPsc_rgb'
        sessN=[1:3];
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
                flat = caret_load(fullfile(surfaceGroupDir,[hem{h},'.FLAT.coord']));
                topo = caret_load(fullfile(surfaceGroupDir,[hem{h},'.CLOSED.topo']));
                
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
                
                name={sprintf('psc_trained_sess%d',ss)};
                %name={sprintf('psc_sess%d',sessN),sprintf('dist_sess%d',sessN),sprintf('dist_cross_sess%d',sessN)};
                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1},'column_name',name);
                
                if smooth == 1
                    caret_save([surfaceGroupDir filesep 's' hem{h} sprintf('.PSC_trained_sess%d.RGB_paint',ss)],C);
                else
                    caret_save([surfaceGroupDir filesep hem{h} sprintf('.PSC_trained_sess%d.RGB_paint',ss)],C);
                end
            end;
        end
    case 'SURF_mapTessels'
        sn=[4:9,11:31];
        sessN=[1:4];
        parcelType='162tessels';
        betaChoice='multi';
        thres=[0 0.01];
        vararginoptions(varargin,{'sn','sessN','parcelType','betaChoice','thres'});
        for ss=sessN
            T=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
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
                    fprintf('Done sess-%d %s\n',ss,subj_name{s});
                end
            end
        end
    case 'SURF_mapTessels_group'
        %sn=[5,7,8,9,11,13,14,15,16];
        sn=[4:9,11:31]; % 1-3: different resolution, 6, 12- excessive movement, 10- no session 4
        sessN=[1:4];
        regType='162tessels';
        thres=[-0.01 0.01];
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
        sn=[4:9,11:31];
        sessTr=[1:3];
        regType='162tessels';
        thres=[-0.01 0.01];
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
        
    case 'SURF_wb:mapTessels_mahalanobis'
        sessN = 1:4;
        nNodes = 162;
        betaChoice = 'multi';
        vararginoptions(varargin,{'sessN','nNodes','betaChoice'});
        
        for h=1:2 % hemisphere
            for ss=sessN
                T = load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d.mat',nNodes,betaChoice,ss)));
                G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nNodes,hemI{h})));
                nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
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
            end
            save(G, outfile);
            fprintf('Hemisphere: %d\n',h);
        end
    case 'SURF_wb:mapTessels_mahalanobis_searchlight'
        sessN = 1:4;
        nNodes = 162;
        vararginoptions(varargin,{'sessN','nNodes','betaChoice'});
        
        for h=1:2 % hemisphere
            for ss=sessN
                T = load(fullfile(betaDir,'group',sprintf('betas_tesselsWB_%d_sess%d.mat',nNodes,ss)));
                G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nNodes,hemI{h})));
                nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
                for r=1:nReg
                    t = getrow(T,T.regType==r&T.regSide==h);
                    dist_T = zeros(size(t.dist_train,1),1);
                    dist_UT = dist_T;
                    for s=1:size(t.dist_train,1)
                        if isnan(t.dist_train{s})
                            dist_T(s) = NaN;
                            dist_UT(s) = NaN;
                        else
                            dist_T(s) = nanmean(t.dist_train{s});
                            dist_UT(s) = nanmean(t.dist_untrain{s});
                        end
                    end
                    idx = find(G.cdata(:,1)==r);
                    if sum(isnan(dist_T))/length(dist_T)>0.7% ensure more than 70% have data
                        data(idx,(ss-1)*2+1) = NaN;
                        data(idx,(ss-1)*2+2) = NaN;
                    else
                        data(idx,(ss-1)*2+1) = nanmean(dist_T);
                        data(idx,(ss-1)*2+2) = nanmean(dist_UT);
                    end
                    column_name{(ss-1)*2+1} = sprintf('dist_trained-sess%d',ss);
                    column_name{(ss-1)*2+2} = sprintf('dist_untrained-sess%d',ss);
                    fprintf('%d.',r);
                end
            end
            outfile         = fullfile(wbDir,'FS_LR_164',sprintf('%s.tesselsWB_%d-search_dist.func.gii',hemI{h},nNodes));
            G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            % filter out any crazy values (on the border of mask)
            for c=1:size(G.cdata,2)
                G.cdata(G.cdata(:,c)>5,c)=0;
                G.cdata(G.cdata(:,c)<-5,c)=0;
            end
            save(G, outfile);
            fprintf('Hemisphere: %d\n',h);
        end
    case 'SURF_wb:mapTessels_mahalanobis_crossSess'
        sessTr = 1:3;
        nNodes = 162;
        betaChoice = 'multi';
        vararginoptions(varargin,{'sessTr','nNodes','betaChoice'});
        
        for h=1:2 % hemisphere
            for ss=sessTr
                T1 = load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d.mat',nNodes,betaChoice,ss)));
                T2 = load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d.mat',nNodes,betaChoice,ss+1)));
                G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nNodes,hemI{h})));
                nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
                for r=1:nReg
                    t1 = getrow(T1,T1.regType==r&T1.regSide==h);
                    t2 = getrow(T2,T2.regType==r&T2.regSide==h);
                    idx = find(G.cdata(:,1)==r);
                    if sum(isnan(t1.dist_train))/length(t1.dist_train)>0.7% ensure more than 70% have data
                        data(idx,(ss-1)*2+1) = NaN;
                        data(idx,(ss-1)*2+2) = NaN;
                    else
                        data(idx,(ss-1)*2+1) = nanmean(t2.dist_train)-nanmean(t1.dist_train);
                        data(idx,(ss-1)*2+2) = nanmean(t2.dist_untrain)-nanmean(t1.dist_untrain);
                    end
                    column_name{(ss-1)*2+1} = sprintf('dist_trained-sessTr_%d-%d',ss,ss+1);
                    column_name{(ss-1)*2+2} = sprintf('dist_untrained-sessTr_%d-%d',ss,ss+1);
                end
            end
            outfile         = fullfile(wbDir,'FS_LR_164',sprintf('%s.tesselsWB_%d_sessTrans-dist.func.gii',hemI{h},nNodes));
            G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            % filter out any crazy values (on the border of mask)
            for c=1:size(G.cdata,2)
                G.cdata(G.cdata(:,c)>5,c)=0;
                G.cdata(G.cdata(:,c)<-5,c)=0;
            end
            save(G, outfile);
            fprintf('Hemisphere: %d\n',h);
        end
    case 'SURF_wb:mapTessels_mahalanobis_difference'
        % trained - untrained
        sessN = 1:4;
        nNodes = 162;
        betaChoice = 'multi';
        vararginoptions(varargin,{'sessN','nNodes','betaChoice'});
        
        for h=1:2 % hemisphere
            for ss=sessN
                T = load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d.mat',nNodes,betaChoice,ss)));
                G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nNodes,hemI{h})));
                nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
                for r=1:nReg
                    t = getrow(T,T.regType==r&T.regSide==h);
                    idx = G.cdata(:,1)==r;
                    if sum(isnan(t.dist_train))/length(t.dist_train)>0.7% ensure more than 70% have data
                        data(idx,ss) = NaN;
                    else
                        data(idx,ss) = nanmean(t.dist_train)-nanmean(t.dist_untrain);
                    end
                    column_name{ss} = sprintf('distDiff-sess%d',ss);
                end
            end
            outfile         = fullfile(wbDir,'FS_LR_164',sprintf('%s.tesselsWB_%d-distDifference.func.gii',hemI{h},nNodes));
            G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            % filter out any crazy values (on the border of mask)
            for c=1:size(G.cdata,2)
                G.cdata(G.cdata(:,c)>5,c)=0;
                G.cdata(G.cdata(:,c)<-5,c)=0;
            end
            save(G, outfile);
            fprintf('Hemisphere: %d\n',h);
        end
    case 'SURF_wb:mapTessels_cosine'
        sessN = 1:4;
        nNodes = 162;
        vararginoptions(varargin,{'sessN','nNodes'});
        
        for h=1:2 % hemisphere
            for ss=sessN
                T = load(fullfile(distPscDir,sprintf('corrDist_tesselsWB_%d_ROI.mat',nNodes)));
                G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nNodes,hemI{h})));
                nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
                for r=1:nReg
                    t = getrow(T,T.regType==r&T.regSide==h&T.sessN==ss);
                    idx = find(G.cdata(:,1)==r);
                    if sum(isnan(t.corrDist))/length(t.corrDist)>0.7% ensure more than 70% have data
                        data(idx,(ss-1)*2+1) = NaN;
                        data(idx,(ss-1)*2+2) = NaN;
                    else
                        data(idx,(ss-1)*2+1) = nanmean(t.corrDist(t.seqType==1))*100;
                        data(idx,(ss-1)*2+2) = nanmean(t.corrDist(t.seqType==2))*100;
                    end
                    column_name{(ss-1)*2+1} = sprintf('dist_corr_trained-sess%d',ss);
                    column_name{(ss-1)*2+2} = sprintf('dist_corr_untrained-sess%d',ss);
                end
            end
            outfile         = fullfile(wbDir,'FS_LR_164',sprintf('%s.tesselsWB_%d-corr-dist.func.gii',hemI{h},nNodes));
            G               = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
            % filter out any crazy values (on the border of mask)
            for c=1:size(G.cdata,2)
                G.cdata(G.cdata(:,c)>5,c)=0;
                G.cdata(G.cdata(:,c)<-5,c)=0;
            end
            save(G, outfile);
            fprintf('Hemisphere: %d\n',h);
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
    case 'PLOT_eig_dim2'
        roi = 1:8;
        sessN=1:4;
        vararginoptions(varargin,{'roi','sessN','hemi'});
        seqType={'trained','untrained'};
        
        T=load(fullfile(distPscDir,'eigenvalues.mat'));
        T.numDim=sum(T.eig(:,1:5)>0,2); % number of non-negative eigenvalues

        for r=roi
            style.use('SeqShade_small');
            figure
            for ss=sessN
                subplot(3,5,ss)
               % histogram(T.numDim(T.roi==r & T.sessN==ss&T.seqType==2));
               % hold on; histogram(T.numDim(T.roi==r & T.sessN==ss&T.seqType==1));
                plt.hist(T.numDim,'subset',T.roi==r & T.sessN==ss,'split',T.seqType,'leg',seqType);
                title(sprintf('%s - sess%d',regname_cortex{r},ss)); 
                if ss==1
                    ylabel('N'); xlabel('dimension');
                else
                    ylabel('');
                end
            end
            style.use('Seq');
            for d=1:5 % per dimension
                subplot(3,5,5+d)
                plt.line([T.sessN>3 T.sessN],T.eig(:,d),'subset',T.roi==r & T.numDim>=d,'split',T.seqType,'leg',seqType);
                title(sprintf('dimension-%d',d));
                if d==1
                    ylabel('variance explained'); xlabel('session');
                else
                    ylabel('');
                end
            end
            subplot(3,5,11)
            plt.line([T.sessN>3 T.sessN],sum(T.eig(:,1:T.numDim),2),'subset',T.roi==r,'split',T.seqType);
            title('sum of all dimensions');
            ylabel('variance explained'); xlabel('session');
            style.use('Sess');
            subplot(3,5,12)
            plt.hist(T.numDim,'subset',T.roi==r & T.seqType==1,'split',T.sessN,'leg',{'1','2','3','4'}); title('Trained'); xlabel('dimension');
            subplot(3,5,13)
            plt.hist(T.numDim,'subset',T.roi==r & T.seqType==2,'split',T.sessN,'leg',{'1','2','3','4'}); title('Control'); xlabel('dimension');
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
        roi=[1:16];
        sn=[4:9,11:31];
        sessN=[1:4];
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','sessN','parcelType'});
        
        Dim=[];
        for ss=sessN
            T=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            fprintf('Reg:\t');
            for r=roi
                for s=sn
                    Tt=getrow(T,T.region==r & T.SN==s);
                    % get SPM info
                    TI=load(fullfile(glmSessDir{ss},subj_name{s},'SPM_info'));
                    for st=1:2 % trained / untrained
                        % extract betas for seqType in question (st)
                        tmp=Tt.betaW{:};
                        betas=tmp(TI.seqType==st,:);
                        % extract partition vectors
                        partVec=TI.run(TI.seqType==st);
                        % calculate identity matrix for 6 sequences
                        ind=repmat(indicatorMatrix('identity',num_train),8,1);
                        
                        [corrC_ms,dC_ms]=crossval_estDimCorr(betas,ind,num_train,partVec); % subtracting the mean across conditions
                        [corrC_wm,dC_wm]=crossval_estDimCorr(betas,ind,num_train,partVec,'meanSubtract',0); % keep the mean pattern
                        % create structure
                        N.corr_ms           = corrC_ms(1:5)';
                        N.corr_ms_norm      = (corrC_ms(1:5)./corrC_ms(5))'; % normalised version
                        N.corrDiff_ms       = dC_ms(1:5)';
                        N.corrDiff_ms_norm  = (dC_ms(1:5)./dC_ms(1))';
                        N.corr_wm           = corrC_wm(1:5)';
                        N.corr_wm_norm      = (corrC_wm(1:5)./corrC_wm(5))'; % normalised version
                        N.corrDiff_wm       = dC_wm(1:5)';
                        N.corrDiff_wm_norm  = (dC_wm(1:5)./dC_wm(1))';
                        N.indDim            = (1:5)';
                        N.sn                = ones(size(N.corr_ms))*s;
                        N.roi               = ones(size(N.corr_ms))*r;
                        N.regType           = ones(size(N.corr_ms))*Tt.regType;
                        N.regSide           = ones(size(N.corr_ms))*Tt.regSide;
                        N.sessN             = ones(size(N.corr_ms))*ss;
                        N.seqType           = ones(size(N.corr_ms))*st;
                        
                        Dim=addstruct(Dim,N);
                    end
                end
                fprintf('%d.',r);
            end
            fprintf('\t...done sess-%d.\n',ss);
        end

        % save dimensionality
        save(fullfile(distPscDir,sprintf('Dimensionality_%s.mat',parcelType)),'-struct','Dim');
    case 'CALC_dimensions_overall'
         % calculate dimensions overall for the 12 sequences
        roi=1:16;
        sn=[4:9,11:31];
        sessN=1:4;
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','sessN','parcelType'});
        
        Dim=[];
        for ss=sessN
            T=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            fprintf('Reg:\t');
            for r=roi
                for s=sn
                    Tt=getrow(T,T.region==r & T.SN==s);
                    % get SPM info
                    TI=load(fullfile(glmSessDir{ss},subj_name{s},'SPM_info'));
                    % extract betas for seqType in question (st)
                    betas=Tt.betaW{:};
                    % extract partition vectors
                    partVec=TI.run;
                    % calculate identity matrix for all 12 sequences
                    ind=repmat(indicatorMatrix('identity',1:12),8,1);
                    
                    [corrC_ms,dC_ms]=crossval_estDimCorr(betas,ind,1:12,partVec); % subtracting the mean across conditions
                    [corrC_wm,dC_wm]=crossval_estDimCorr(betas,ind,1:12,partVec,'meanSubtract',0); % keep the mean pattern
                    % create structure
                    N.corr_ms           = corrC_ms(1:11)';
                    N.corr_ms_norm      = (corrC_ms(1:11)./corrC_ms(11))'; % normalised version
                    N.corrDiff_ms       = dC_ms(1:11)';
                    N.corrDiff_ms_norm  = (dC_ms(1:11)./dC_ms(1))';
                    N.corr_wm           = corrC_wm(1:11)';
                    N.corr_wm_norm      = (corrC_wm(1:11)./corrC_wm(11))'; % normalised version
                    N.corrDiff_wm       = dC_wm(1:11)';
                    N.corrDiff_wm_norm  = (dC_wm(1:11)./dC_wm(1))';
                    N.indDim            = (1:11)';
                    N.sn                = ones(size(N.corr_ms))*s;
                    N.roi               = ones(size(N.corr_ms))*r;
                    N.regType           = ones(size(N.corr_ms))*Tt.regType;
                    N.regSide           = ones(size(N.corr_ms))*Tt.regSide;
                    N.sessN             = ones(size(N.corr_ms))*ss;
                    
                    Dim=addstruct(Dim,N);
                end
                fprintf('%d.',r);
            end
            fprintf('\t...done sess-%d.\n',ss);
        end
        
        % save dimensionality
        save(fullfile(distPscDir,sprintf('Dimensionality_overallSeq_%s.mat',parcelType)),'-struct','Dim');
    case 'PLOT_dimensions'
        roi=3;
        parcelType='Brodmann';
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'roi','parcelType'});
        
        D=load(fullfile(distPscDir,sprintf('Dimensionality_%s.mat',parcelType)));
        D=getrow(D,D.roi==roi);
            
        figure
        for st=1:2
            subplot(1,2,st)
            plt.line(D.indDim,D.corr_ms,'split',D.sessN,'subset',D.seqType==st,'style',stySess);
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
        
        D1 = tapply(D,{'sessN','sn','seqType'},{'corr_ms','sum'});
        figure
        plt.line([D1.sessN>3 D1.sessN],D1.corr_ms,'split',D1.seqType);
        
        figure
        for st=1:2
            subplot(1,2,st)
            plt.line(D.indDim,D.corr_ms_norm,'split',D.sessN,'subset',D.seqType==st,'style',stySess);
            drawline(0,'dir','horz');
            plt.match('y');
            hold on;
            drawline(1,'dir','horz');
            title(sprintf('%s-%s',regname{roi},seqLabel{st}));
            if st==1
                ylabel('Correlations');
                xlabel('Possible dimensions');
            else
                ylabel('');
            end
        end
    case 'PLOT_dimensions_overall'
        roi=3;
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','parcelType'});
        
        D=load(fullfile(distPscDir,sprintf('Dimensionality_overallSeq_%s.mat',parcelType)));
        D=getrow(D,D.roi==roi);
        
        figure
        subplot(121)
        plt.line(D.indDim,D.corr_ms,'split',D.sessN,'style',stySess);
        drawline(0,'dir','horz');
        plt.match('y');
        title(sprintf('%s',regname{roi}));
        ylabel('Correlations');
        xlabel('Possible dimensions');
        
        
        subplot(122)
        plt.line(D.indDim,D.corr_ms_norm,'split',D.sessN,'style',stySess);
        drawline(0,'dir','horz');
        hold on;
        drawline(1,'dir','horz');
        title(sprintf('%s',regname{roi}));
        ylabel('Correlations');
        xlabel('Possible dimensions');
        
    case 'STATS_dimensions'
        parcelType='Brodmann';
        roi=1:8;
        vararginoptions(varargin,{'parcelType','roi'});
        D=load(fullfile(distPscDir,sprintf('Dimensionality_%s.mat',parcelType)));
     
        for r=roi
            T=getrow(D,D.roi==r); 
            fprintf('%s - trained\n',regname{r});
            anovaMixed(T.corr_ms,T.sn,'within',[T.sessN T.indDim],{'session','dimension'},'subset',T.seqType==1);
            fprintf('%s - untrained\n',regname{r});
            anovaMixed(T.corr_ms,T.sn,'within',[T.sessN T.indDim],{'session','dimension'},'subset',T.seqType==2);
            T1 = tapply(T,{'sessN','sn','seqType'},{'corr_ms','sum'});
            T2 = tapply(T,{'sessN','sn','seqType'},{'corr_ms','sum'},'subset',T.indDim~=1);
            anovaMixed(T1.corr_ms,T1.sn,'within',[T1.sessN,T1.seqType],{'session','seqType'},'subset',T1.sessN~=4);
            anovaMixed(T2.corr_ms,T2.sn,'within',[T2.sessN,T2.seqType],{'session','seqType'},'subset',T2.sessN~=4);
            ttestDirect(T2.corr_ms,[T2.seqType T2.sn],2,'paired','split',T2.sessN);
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
        sn = [5:9,11:31];
        roi = [1:3,7,8];
        sessN = 1:4;
        betaChoice='multiPW';
        parcelType='Brodmann';
        hemi = 1;
        numError = 1;
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice','parcelType','hemi','numError'});
        
        Stats = [];
        for ss=sessN
            %T = load(fullfile(betaDir,'group',sprintf('stats_%s_%s_sess%d.mat',parcelType,betaChoice,ss))); % loads region data (D)
            T = load(fullfile(betaDir,'group',sprintf('statsError%d_%s_%s_sess%d.mat',numError,parcelType,betaChoice,ss))); % loads region data (D)
            for s=sn
                for h=hemi
                    for r=roi
                        R = getrow(T,T.regType==r & T.SN==s & T.regSide==h);
                        S.dist_train    = R.dist_train;
                        S.dist_untrain  = R.dist_untrain;
                        S.dist_cross    = R.dist_cross;
                        S.sn            = s;
                        S.roi           = r;
                        S.regType       = R.regType;
                        S.regSide       = R.regSide;
                        S.sessN         = ss;
                        Stats = addstruct(Stats,S);
                    end
                end
            end
        end
        % save structure
        %save(fullfile(distPscDir,sprintf('dist_%s_ROI',parcelType)),'-struct','Stats');
        save(fullfile(distPscDir,sprintf('distError%d_%s_ROI',numError,parcelType)),'-struct','Stats');
    case 'PLOT_dist'
        roi=[1:8];
        sessN=[1:4];
        hemi=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        T=load(fullfile(distPscDir,sprintf('dist_%s_ROI.mat',parcelType)));
        switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        % create structure with trained / untrained together, seqType
        D.dist=[T.dist_train;T.dist_untrain];
        D.seqType=[ones(size(T.dist_train));ones(size(T.dist_train))*2];
        D.sn=[T.sn;T.sn];
        D.roi=[T.roi;T.roi];
        D.regType=[T.regType;T.regType];
        D.regSide=[T.regSide;T.regSide];
        D.sessN=[T.sessN;T.sessN];
        
        indx=1;
        figure
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.dist),'split',D.seqType,'subset',D.regType==r & D.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.dist),'split',D.seqType,'subset',D.regType==r & D.regSide==hemi & D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
             end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regLab{r}));
            if r==1
                ylabel('Mahalanobis distance');
                xlabel('Session');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'PLOT_dist_individDiff'
        roi=[1:8];
        hemi=1;
        parcelType='Brodmann';
        var='x1_mean';
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi','var'});
        
        T=load(fullfile(distPscDir,sprintf('corrDist_seqType_%s_ROI.mat',parcelType)));
        % load behaviour;
        B=load(fullfile(baseDir,'behavioral_data','analyze','exponentialFits'));
        B=getrow(B,B.sn~=4);
        for r=roi
            figure
            for ss=1:4
                subplot(1,4,ss)
              %  t=getrow(T,T.roi==r&T.regSide==hemi&T.sn~=4&T.sessN==ss);
             %   plt.scatter(B.(var),t.dist_train);
                t=getrow(T,T.roi==r & T.regSide==hemi & T.sn~=4 & T.sessN==ss & T.seqType==1);
                plt.scatter(B.(var),t.corrDist);
                title(sprintf('%s - sess%d',regname{r},ss));
            end
        end
    case 'PLOT_dist_sessions'
        sessN=[1:3];
        regExcl=6;
        seqType='trained';
        parcelType='Brodmann';
        vararginoptions(varargin,{'sessN','regExcl','seqType'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        switch (seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        T=load(fullfile(distPscDir,sprintf('dist_%s_ROI.mat',parcelType)));
        
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
        parcelType='Brodmann';
        hemi=1;
        
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        D=load(fullfile(distPscDir,sprintf('dist_%s_ROI.mat',parcelType)));
         switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.dist_cross),'subset',D.regType==r & D.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.dist_cross),'subset',D.regType==r & D.regSide==hemi & D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==2
                plt.line(D.sessN,ssqrt(D.dist_cross),'subset',D.regType==r & D.regSide==hemi&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regLab{r}));
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
        parcelType='Brodmann';
        vararginoptions(varargin,{'sessN','regExcl','seqType','parcelType'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(distPscDir,sprintf('dist_%s_ROI.mat',parcelType)));
        
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
        roi=1:8;
        sessN=1:4;
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
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.dist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r&T.regSide==hemi);
        end
        keyboard;
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.dist,[T.seqType T.sn],2,'paired','subset',T.regType==r&T.regSide==hemi&T.sessN==ss);
                fprintf('trained seq vs. 0\n');
                ttestDirect(T.dist,[T.sn],2,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==1);
                fprintf('untrained seq vs. 0\n');
                ttestDirect(T.dist,[T.sn],2,'onesample','subset',T.regType==r&T.regSide==hemi&T.sessN==ss&T.seqType==2);
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
        
    case 'CALC_corrDist'
        reg = 1:8;
        sn  = [5:9,11:31];
        sessN = 1:4;
        parcelType='Brodmann'; %'tesselsWB_162'
        regType = 'all';
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        hemi = 1;
        numError=2;
        vararginoptions(varargin,{'sn','reg','sessN','subtract_mean','parcelType','hemi','numError'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
            D   = load(fullfile(betaDir,'group',sprintf('betasError%d_%s_sess%d.mat',numError,parcelType,ss)));
            %D   = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for h=hemi
                if strcmp(regType,'all')
                    reg = unique(D.regType(D.regSide==h))';
                end
                for roi = reg;
                    sprintf('Done regions:');
                    for s = sn;
                        %SI = load(fullfile(glmSessDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                        SI = load(fullfile(glmErrorDir{numError}{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                        t   = getrow(D,D.regType==roi & D.regSide==h & D.SN==s);
                        if size(t.betaW{1},2)~=1 % only for NaNs
                            data = t.betaW{1}(1:size(SI.SN,1),:);
                            numRuns    = 1:numruns_task_sess;
                            numConds   = num_seq;
                            conds   = repmat(numConds,1,length(numRuns))';
                            runNums = kron(numRuns,ones(1,length(numConds)))';
                            if subtract_mean
                                for r=numRuns
                                    data(runNums==r,:) = bsxfun(@minus,data(runNums==r,:),mean(data(runNums==r,:)));
                                end
                            else
                                data=data;
                            end
                            % crossval second moment matrix
                            [G,~]     = pcm_estGCrossval(data,SI.run,SI.seqNumb);
                            C=corr_crossval(G,'reg','minvalue');
                            C=rsa_squareRDM(C);
                            
                            % average trained dist
                            S.corrDist(1,:) = sum(sum(triu(C(1:6,1:6))))/(6*5/2);
                            % average untrained dist
                            S.corrDist(2,:) = sum(sum(triu(C(7:12,7:12))))/(6*5/2);
                            S.seqType   = [1;2]; % trained, untrained
                            S.sessN     = [ss;ss];
                            S.roi       = [roi;roi];
                            S.regType   = repmat(t.regType,2,1);
                            S.regSide   = repmat(t.regSide,2,1);
                            S.sn        = [s;s];
                            SAll=addstruct(SAll,S);
                            
                            ST.corr_seqType = sum(sum(triu(C(7:12,1:6))))/(6*5);
                            ST.sessN        = ss;
                            ST.roi          = roi;
                            ST.regType      = t.regType;
                            ST.regSide      = t.regSide;
                            ST.sn           = s;
                            STAll=addstruct(STAll,ST);
                        end
                    end
                    fprintf('%d.',roi);
                end
            end
            fprintf('\nDone sess%d.\n',ss);
        end
        save(fullfile(distPscDir,sprintf('corrDistError%d_%s_ROI.mat',numError,parcelType)),'-struct','SAll');
        save(fullfile(distPscDir,sprintf('corrDistError%d_seqType_%s_ROI.mat',numError,parcelType)),'-struct','STAll');
        %save(fullfile(distPscDir,sprintf('corrDist_%s_ROI.mat',parcelType)),'-struct','SAll');
        %save(fullfile(distPscDir,sprintf('corrDist_seqType_%s_ROI.mat',parcelType)),'-struct','STAll');
    case 'CALC_corrDist_wholeStriatum'
            % combine caudate and pallidum
        reg = 1:8;
        sn  = [5:9,11:31];
        sessN = 1:4;
        parcelType='Brodmann'; %'tesselsWB_162'
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        vararginoptions(varargin,{'sn','reg','sessN','subtract_mean','parcelType'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
            D   = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for h=1:2
                    sprintf('Done regions:');
                    for s = sn;
                        SI = load(fullfile(glmSessDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                        t   = getrow(D,ismember(D.regType,1:2) & D.regSide==h & D.SN==s);
                        if size(t.betaW{1},2)~=1 || size(t.betaW{2},2)~=1 % only for NaNs
                            data = [t.betaW{1}(1:size(SI.SN,1),:) t.betaW{2}(1:size(SI.SN,1),:)];
                            numRuns    = 1:numruns_task_sess;
                            numConds   = num_seq;
                            conds   = repmat(numConds,1,length(numRuns))';
                            runNums = kron(numRuns,ones(1,length(numConds)))';
                            if subtract_mean
                                for r=numRuns
                                    data(runNums==r,:) = bsxfun(@minus,data(runNums==r,:),mean(data(runNums==r,:)));
                                end
                            else
                                data=data;
                            end
                            % crossval second moment matrix
                            [G,~]     = pcm_estGCrossval(data,SI.run,SI.seqNumb);
                            C=corr_crossval(G,'reg','minvalue');
                            C=rsa_squareRDM(C);
                            
                            % average trained dist
                            S.corrDist(1,:) = sum(sum(triu(C(1:6,1:6))))/(6*5/2);
                            % average untrained dist
                            S.corrDist(2,:) = sum(sum(triu(C(7:12,7:12))))/(6*5/2);
                            S.seqType   = [1;2]; % trained, untrained
                            S.sessN     = [ss;ss];
                            S.regType   = [1;1];
                            S.regSide   = [h;h];
                            S.sn        = [s;s];
                            SAll=addstruct(SAll,S);
                            
                            ST.corr_seqType = sum(sum(triu(C(7:12,1:6))))/(6*5);
                            ST.sessN        = ss;
                            ST.regType      = 1;
                            ST.regSide      = h;
                            ST.sn           = s;
                            STAll=addstruct(STAll,ST);
                    end
                end
            end
            fprintf('\nDone sess%d.\n',ss);
        end
        save(fullfile(distPscDir,sprintf('corrDist_wholeStriatum_ROI.mat')),'-struct','SAll');
        save(fullfile(distPscDir,sprintf('corrDist_seqType_wholeStriatum_ROI.mat')),'-struct','STAll');
    case 'PLOT_corrDist'
        roi=[1:8];
        sessN=[1:4];
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        T=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI.mat',parcelType)));
        switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        figure
        indx=1;
        for r=roi
            for h=hemi
                subplot(1,numel(roi),indx)
                if numel(sessN)==4
                    plt.line([T.sessN>3 T.sessN],T.corrDist,'split',T.seqType,'subset',T.regType==r & T.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
                elseif numel(sessN)==3
                    plt.line(T.sessN,T.corrDist,'split',T.seqType,'subset',T.regType==r & T.regSide==hemi&T.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
                end
                drawline(0,'dir','horz');
                plt.match('y');
                title(sprintf('%s',regLab{r}));
                if r==1
                    ylabel('Correlation distance');
                    xlabel('Session');
                else
                    ylabel('');
                end
                indx=indx+1;
            end
        end
    case 'PLOT_corrDist_sessions'
        sessN=[1:3];
        regExcl=6;
        seqType='trained';
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'sessN','regExcl','seqType','parcelType','hemi'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        switch (seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        T=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI.mat',parcelType)));
        
        T=getrow(T,~ismember(T.regType,regExcl)&ismember(T.sessN,sessN)&T.regSide==hemi); % exclude V12

        % list in terms of anterior -> posterior
        reg=[3 4 5 1 2 7 8];
        T.reg=T.regType;
        for i=1:length(T.regType) 
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
        parcelType='Brodmann';
        seqType={'trained','untrained'};
        hemi=1;
        vararginoptions(varargin,{'sessN','regExcl','seqType','parcelType','hemi'});
        
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI.mat',parcelType)));
        
        T=getrow(T,~ismember(T.regType,regExcl)&ismember(T.sessN,sessN)&T.regSide==hemi); % exclude V12

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
        parcelType='Brodmann';
        hemi=1;
        
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        D=load(fullfile(distPscDir,sprintf('corrDist_seqType_%s_ROI.mat',parcelType)));
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],D.corr_seqType,'subset',D.regType==r&D.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==3
                plt.line(D.sessN,D.corr_seqType,'subset',D.regType==r&D.regSide==hemi&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==2
                plt.line(D.sessN,D.corr_seqType,'subset',D.regType==r&D.regSide==hemi&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
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
        parcelType='Brodmann';
        hemi=1;
        vararginoptions(varargin,{'roi','sessN','parcelType','hemi'});
        
        T=load(fullfile(distPscDir,sprintf('corrDist_%s_ROI.mat',parcelType)));

        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.regType==r&T.regSide==hemi);
            %anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r&(T.sessN==1|T.sessN==4));
        end
        
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.corrDist,[T.seqType T.sn],2,'paired','subset',T.regType==r&T.sessN==ss&T.regSide==hemi);
            end
        end
    case 'STATS_corrDist_seqType'
        roi=[1:8];
        sessN=[1:4];
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'roi','sessN'});
        
        T=load(fullfile(distPscDir,sprintf('corrdist_seqType_%s_ROI.mat',parcelType)));
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
        sn = [4:9,11:31];
        roi = [1:8];
        sessN = [1:4];
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
                        
                        Stats=addstruct(Stats,S);
                    end
                end
            end
        end
           
        % save structure
        save(fullfile(distPscDir,sprintf('psc_%s_ROI',parcelType)),'-struct','Stats');
    case 'PLOT_psc'
        roi=[1:8];
        hemi=1;
        sessN=[1:4];
        parcelType='Brodmann';
        
        vararginoptions(varargin,{'roi','sessN','hemi','parcelType'});
        
        T=load(fullfile(distPscDir,sprintf('psc_%s_ROI.mat',parcelType)));
        
        switch parcelType
            case 'Brodmann'
                regLab = regname_cortex;
            case 'BG-striatum'
                regLab = regname_BG;
        end
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([T.sessN>3 T.sessN],T.psc,'split',T.seqType,'subset',T.regType==r&T.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(T.sessN,T.psc,'split',T.seqType,'subset',T.regType==r&T.sessN<4&T.regSide==hemi,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
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
    case 'PCM_constructModelFamily_simple'
        % here simpler version of model family:
        % includes first finger, all fingers, seqType, sequences (T+U)
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31]; 
        sessN = 1:4;
        reg = 1:8;
        hemi = 1;
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllRegSess=[];
        KKK=[]; GGG=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN,'parcelType',parcelType)';
            end
        elseif strcmp(parcelType,'tesselsWB_642') % for 642 parcels
            reg = sml1_imana_dist('CLUSTER_choose_all','sessN',sessN,'parcelType',parcelType)';
            BB = load(fullfile(betaDir,'s05','betas_tesselsWB_642_s05_sess1')); % used for 642
            BB = getrow(BB,ismember(BB.regSide,hemi) & ismember(BB.region,reg));
            reg = BB.region';
            regType = BB.regType';
            regSide = BB.regSide';
        elseif strcmp(parcelType,'Yokoi_clusters')
            BB = load(fullfile(betaDir,'s05','betas_Yokoi_clusters_s05_sess1')); % 
            reg = BB.region';
            regType = BB.regType';
            regSide = BB.regSide';
        end    
        for ss = sessN
            AllReg=[];
            KK=[]; GG=[];
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
%             if ~strcmp(parcelType,'162tessels')
%                 reg=unique(B.region)';
%                 regSide=(reg>size(reg,2)/2)+1;
%                 regType=reg-((regSide-1)*size(reg,2)/2);
%             else % only for 162tessels
%                 regType=reg;
%                 regType(regType>158)=regType(regType>158)-158;
%                 regSide=ones(size(regType));
%                 regSide(reg>158)=2;
%             end
            for r = 1:length(reg)
                for p=1:length(sn)          
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));               
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==reg(r))}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==reg(r))}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                    end  
                    partVec{p} = D.run;  % runs/partitions
                    if ~naturalStats
                        m = pcm_defineSequenceModels_fixed_noChunks(Seq,sn(p));
                    else
                        load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                        NatStat = NatStats.G_cent;
                        m = pcm_defineSequenceModels_fixed_noChunks_natStats(Seq,sn(p),NatStat); % make
                    end
                    [M,Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                    [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                    condVec{p} = repmat(Z,max(D.run),1);
                    indx = ones(size(D.run));
                    Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                end;
                % fit models
                T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                T.roi = ones(size(T.SN))*reg(r);
                T.regType = ones(size(T.SN))*regType(r);
                T.regSide = ones(size(T.SN))*regSide(r);
                T.sessN = ones(size(T.SN))*ss;
                T=rmfield(T,{'reg','thetaCr'});
                AllReg=addstruct(AllReg,T);
                
                % calculations - posterior probability, knockIN/OUT
                K = pcm_calc(T.cross_likelihood,Comb);
                K.roi = ones(length(K.indx),1)*reg(r);
                K.regType = ones(length(K.indx),1)*regType(r);
                K.regSide = ones(length(K.indx),1)*regSide(r);
                K.sessN = ones(length(K.indx),1)*ss;
                KK = addstruct(KK,K);
                idxST = Comb(:,3)==1; % always including SeqType
                G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,4,5]));
                G.roi = ones(length(G.indx),1)*reg(r);
                G.regType = ones(length(G.indx),1)*regType(r);
                G.regSide = ones(length(G.indx),1)*regSide(r);
                G.sessN = ones(length(G.indx),1)*ss;
                GG = addstruct(GG,G);
                fprintf('Done sess-%d reg-%d/%d\n',ss,r,length(reg));
            end
            save(fullfile(pcmDir,sprintf('ModelFamily_Fit_simple_%s_naturalStats-%d_sess-%d.mat',parcelType,naturalStats,ss)),'-struct','AllReg');
            save(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_%s_naturalStats-%d_sess-%d.mat',parcelType,naturalStats,ss)),'-struct','KK');
            save(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_noSeqType_%s_naturalStats-%d_sess-%d.mat',parcelType,naturalStats,ss)),'-struct','GG');
            KKK = addstruct(KKK,KK);
            GGG = addstruct(GGG,GG);
            AllRegSess = addstruct(AllRegSess,AllReg);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing','SeqType','Trained','Untrained'};
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_simple_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllRegSess');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KKK');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_noSeqType_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','GGG');
    case 'PCM_combineSess'
        parcelType = 'tesselsWB_642';
        sessN = 1:4;
        naturalStats = 1;
        vararginoptions(varargin,{'sessN','naturalStats','parcelType'});
        GG = []; KK = []; AllRegSess = [];
        for ss=sessN
            AllReg = load(fullfile(pcmDir,sprintf('ModelFamily_Fit_simple_%s_naturalStats-%d_sess-%d.mat',parcelType,naturalStats,ss)));
            K = load(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_%s_naturalStats-%d_sess-%d.mat',parcelType,naturalStats,ss)));
            G = load(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_noSeqType_%s_naturalStats-%d_sess-%d.mat',parcelType,naturalStats,ss)));
            GG = addstruct(GG,G);
            KK = addstruct(KK,K);
            AllRegSess = addstruct(AllRegSess,AllReg);         
        end
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllRegSess');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KK');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_noSeqType_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','GG');
    case 'PCM_PLOT_modelFamily_simple'
        parcelType='Brodmann';
        reg=[2,3];
        modelType = 'seqType'; % seqType of noSeqType
        metric='logBayes'; % logBayes, knockIN
        naturalStats = 0;
        vararginoptions(varargin,{'reg','parcelType','modelType','metric','naturalStats'});
        switch modelType
            case 'seqType'
                KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)));
                name = {'FF','AF','ST','TS','US'};
            case 'noSeqType'
                KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_noSeqType_%s_naturalStats-%d.mat',parcelType,naturalStats)));
                name = {'FF','AF','TS','US'};
        end
        for r = reg
            % posterior Log
            figure;
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.(metric));
                set(gca,'XTickLabel',name);
                drawline(0,'dir','horz');
                plt.match('y');
                if s==1
                    ylabel('Log Bayes');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
        end    
    case 'PCM_noiseCeilings'
        % construct noise ceilings
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31]; 
        sessN = 1:4;
        reg = 1:8;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType'})
        AllReg=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        % here models - null and noise ceiling
        M{1}.type       = 'feature';
        M{1}.numGparams = 1;
        M{1}.name       = 'null';
        M{1}.Ac(:,1:12,1)  = zeros(12);
        M{2}.type       = 'freedirect';
        M{2}.numGparams = 0;
        M{2}.theta0     = [];
        M{2}.name       = 'noice_ceiling';
        % two subject groups - split by even / odd subject assignment
        sn_group{1} = sn(rem(sn,2)==1);
        sn_group{2} = sn(rem(sn,2)==0);
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for g = 1:2 % two groups of subjects
                    for p=1:length(sn_group{g})
                        glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn_group{g}(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn_group{g}(p)&B.region==reg(r))}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn_group{g}(p)&B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn_group{g}(p)&B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                        end
                        partVec{p} = D.run;  % runs/partitions;
                        condVec{p} = D.seqNumb;
                        Data{p} = beta(:,1:length(D.run))';  % Data is N x P (cond x voxels) - no intercept
                    end;
                    % fit models
                    T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T=rmfield(T,{'reg','thetaCr'});
                    if g==2
                        T.SN = T.SN + length(sn_group{1});
                    end
                    AllReg=addstruct(AllReg,T);
                    clear partVec condVec Data;
                end
            end
            fprintf('Done sess-%d.\n\n\n',ss);
        end
        % save variables;
        save(fullfile(pcmDir,sprintf('NoiseCeilings_allSeq_%s.mat',parcelType)),'-struct','AllReg');
    case 'PCM_noiseCeilings_seqType'
        % construct noise ceilings
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31]; 
        sessN = 1:4;
        reg = 1:8;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType'})
        AllReg=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        % here models - null and noise ceiling
        M{1}.type       = 'feature';
        M{1}.numGparams = 1;
        M{1}.name       = 'null';
        M{1}.Ac(:,1:6,1)  = zeros(6);
        M{2}.type       = 'freedirect';
        M{2}.numGparams = 0;
        M{2}.theta0     = [];
        M{2}.name       = 'noice_ceiling';
        % two subject groups - split by even / odd subject assignment
        sn_group{1} = sn(rem(sn,2)==1);
        sn_group{2} = sn(rem(sn,2)==0);
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for st = 1:2 % sequence type
                    for g = 1:2 % two groups of subjects
                        for p=1:length(sn_group{g})
                            glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn_group{g}(p)});
                            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                            switch (beta_choice)
                                case 'uw'
                                    beta = B.betaUW{(B.sn==sn_group{g}(p)&B.region==reg(r))}';
                                case 'mw'
                                    beta = B.betaW{(B.SN==sn_group{g}(p)&B.region==reg(r))}';
                                case 'raw'
                                    beta = B.betaRAW{(B.sn==sn_group{g}(p)&B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                            end
                            partVec{p} = D.run(D.seqType==st);  % runs/partitions;
                            condVec{p} = D.seqNumb(D.seqType==st);
                            Data{p} = beta(:,(D.seqType==st))';  % Data is N x P (cond x voxels) - no intercept
                        end;
                        % fit models
                        T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                        T.roi = ones(size(T.SN))*reg(r);
                        T.regType = ones(size(T.SN))*regType(r);
                        T.regSide = ones(size(T.SN))*regSide(r);
                        T.sessN = ones(size(T.SN))*ss;
                        T.seqType = ones(size(T.SN))*st;
                        T=rmfield(T,{'reg','thetaCr'});
                        if g==2
                            T.SN = T.SN + length(sn_group{1});
                        end
                        AllReg=addstruct(AllReg,T);
                        clear partVec condVec Data;
                    end % group
                end % seqType
            end
            fprintf('Done sess-%d.\n\n\n',ss);
        end
        % save variables;
        save(fullfile(pcmDir,sprintf('NoiseCeilings_seqType_%s.mat',parcelType)),'-struct','AllReg');
    case 'PCM_PLOT_surface_modelFamily_simple'
        % here project results to the surface
        parcelType = 'Yokoi_clusters';
        naturalStats = 1;
        hemi = 1:2;
        sessN = 1:4;
        modelType = 'noSeqType';
        metric = 'logBayes';
        nTessel = 642;
        vararginoptions(varargin,{'parcelType','naturalStats','modelType','metric','hemi','sessN'});
        
        switch modelType
            case 'seqType'
                KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)));
                name = {'FirstFinger','AllFingers','SeqType','TrainedSeq','UntrainedSeq'};
            case 'noSeqType'
                KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_noSeqType_%s_naturalStats-%d.mat',parcelType,naturalStats)));
                name = {'FirstFinger','AllFingers','TrainedSeq','UntrainedSeq'};
        end
        
        for h=hemi
            switch parcelType
                case 'Yokoi_clusters'
                    G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('%s.Yokoi2019.10cluster.label.gii',hemI{h})));
                case 'tesselsWB_642'
                    G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{h})));
            end
            reg = unique(KK.regType(KK.regSide==h));
            %reg = reg(reg~=0);
            % initialize data structure, counter
            idx = 1;
            data = zeros(size(G.cdata,1),numel(sessN)*numel(unique(KK.indx)));
            for ss=sessN
                for i=unique(KK.indx)'
                    for r=reg'
                        t = getrow(KK,KK.sessN==ss&KK.regType==r&KK.regSide==h&KK.indx==i);
                        data(G.cdata==r,idx) = mean(t.(metric));
                    end
                    columnName{idx} = sprintf('sess-%d_%s_%s',ss,metric,name{i});
                    idx = idx+1;
                end
            end
            outfile = fullfile(wbDir,'FS_LR_164',sprintf('%s.PCM.%s.naturalStats-%d.%s_modelFamily_%s.func.gii',hemI{h},parcelType,naturalStats,metric,modelType));
            G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',columnName);
            save(G, outfile);
        end
        
    case 'PCM_noiseCeilings_plot'
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        sessN = 1:4;
        reg = 1:8;
        naturalStats = 0;
        vararginoptions(varargin,{'parcelType','reg','naturalStats'});
        % plot PCM noise ceilings
        T = load(fullfile(pcmDir,sprintf('NoiseCeilings_allSeq_%s.mat',parcelType)));
        K = load(fullfile(pcmDir,sprintf('ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('ModelFamilyComb_simple_%s.mat',parcelType)));
        for r=reg
            for ss=sessN
                figure(r)
                t=getrow(T,T.roi==r&T.sessN==ss);
                k=getrow(K,K.roi==r&T.sessN==ss);
                [~,bestM] = max(mean(k.cross_likelihood));
             %   bestM = mode(j); % determine the best model
                % form a new structure
                data = k.cross_likelihood(:,[1:6,bestM]);
                data = bsxfun(@minus,data,data(:,1));
                M.data = data(:);
                M.indx = kron([1:7]',ones(length(t.SN),1));
                M.sessN = repmat(k.sessN,7,1);
                M.SN = repmat(k.SN,7,1);
                M = normData(M,'data');
                subplot(1,4,ss)
                plt.bar(M.indx,M.normdata,'subset',M.indx~=1); hold on; 
                drawline(mean(t.cross_likelihood(:,2)-t.cross_likelihood(:,1)),'dir','horz'); 
                drawline(0,'dir','horz'); 
                set(gca,'XTickLabel',{'FF','AF','ST','TS','US','bestM'});
                title(sprintf('%s - sess%d ',regname_cortex{r},ss));
                if ss==1
                    ylabel('logBayes');
                else
                    ylabel('');
                end
            end
        end
    case 'PCM_constructModelFamily'  % old
        % includes seqType, chunks, transitions etc.
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        parcelType='162tessels'; % Brodmann or 162tessels
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
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
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
                for p=1:length(sn)          
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));               
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==reg(r))}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==reg(r))}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                    end
                    
                    partVec{p} = D.run;  % runs/partitions
                    switch(runEffect)
                        case 'fixed'
                            %m = pcm_defineSequenceModels_fixed(Seq,p);
                            m = pcm_defineSequenceModels_fixed_new(Seq,SeqChunks,p);
                            [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                            [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        case 'random'
                            m = pcm_defineSequenceModels_random(Seq,p);
                            [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                            [Mf,Comb] = pcm_constructModelFamily(M,'alwaysInclude',1,'fullModel',1);
                    end
                    condVec{p} = repmat(Z,max(D.run),1);
                    indx = ones(size(D.run));
                    Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                end;
                
                % fit models
                T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                T.roi = ones(size(T.SN))*reg(r);
                T.regType = ones(size(T.SN))*regType(r);
                T.regSide = ones(size(T.SN))*regSide(r);
                T.sessN = ones(size(T.SN))*ss;
                T=rmfield(T,{'reg','thetaCr'});
                AllReg=addstruct(AllReg,T);
                
                % calculations - posterior probability, knockIN/OUT
                K = pcm_calc(T.cross_likelihood,Comb);
                K.roi = ones(length(K.indx),1)*reg(r);
                K.regType = ones(length(K.indx),1)*regType(r);
                K.regSide = ones(length(K.indx),1)*regSide(r);
                K.sessN = ones(length(K.indx),1)*ss;
                KK=addstruct(KK,K); 
                fprintf('Done sess-%d reg-%d/%d\n',ss,r,length(reg));
            end
        end
        
        % save variables;
        dircheck(fullfile(pcmDir));
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_%s_new.mat',parcelType)),'Comb');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_%s_fing%s_new.mat',parcelType,fingType)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_fing%s_new.mat',parcelType,fingType)),'-struct','KK');
    
    
    case 'PCM_PLOT_modelFamily'
        %plots knock-IN / OUT / posterior probability
        parcelType='Brodmann';
        reg=1;
        fingType='Count';
        vararginoptions(varargin,{'reg','parcelType'});
        KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_fing%s.mat',parcelType,fingType)));
        name = {'FF','AF','FT','Chu','ST','TS','US'};
        fig_num=0;
        
        for r = reg
            % knock-in
            figure(fig_num+1);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s);
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
                a=getrow(KK,KK.roi==r&KK.sessN==s);
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
                    ylabel('Log Bayes');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
            fig_num=fig_num+5;  
        end    
    case 'PCM_PLOT_modelFamily_allSess'
        %plots knock-IN / OUT / posterior probability
        parcelType='Brodmann';
        reg=1;
        fingType='Count';
        var='logBayes';
        vararginoptions(varargin,{'reg','parcelType','var'});
        KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_fing%s.mat',parcelType,fingType)));
        for r=reg
            figure
            plt.bar(KK.indx,KK.(var),'split',KK.sessN,'style',stySess,'subset',KK.roi==r);
            ylabel(var);
            % set(gca,'XTickLabel',name);
            title(regname{r});
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
        % PCM modelling separately for trained and untrained sequences
        % in use Dec 2019
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        parcelType='Brodmann'; % Brodmann or 162tessels
        modelType = 'allComponents';
        regSelect='all'; % all or subset
        sn=[5:9,11:31]; % change
        sessN=1:4;
        
        vararginoptions(varargin,{'sn','reg','sessN','algorithm','modelType','parcelType','hemi','regSelect','fingType'})
        if strcmp(parcelType,'162tessels')
            reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
        end
        
        AllReg=[];
        K=[];
        KK=[];
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            if strcmp(regSelect,'all')
                reg=unique(B.region)';
                regSide=(reg>size(reg,2)/2)+1;
                regType=reg-((regSide-1)*size(reg,2)/2);
            else
                regType=reg;
                regType(regType>158)=regType(regType>158)-158;
                if strcmp(parcelType,'162tessels')
                    regSide=ones(size(regType));
                    regSide(reg>158)=2;
                end
            end
            
            for r = 1:length(reg)
                for st=1:2
                    for p=1:length(sn)  
                        glmDirSubj=fullfile(glmSessDir{ss},subj_name{sn(p)});     
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));                     
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p)&B.region==r)}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p)&B.region==r)}'; % no intercept - use T.betaRAWint otherwise
                        end                  
                        partVec{p} = D.run(D.seqType==st);  % runs/partitions
                        switch modelType
                            case 'allComponents'
                                m = pcm_defineSequenceModels_seqType(Seq,SeqChunks,p,st);
                                modelNames = {'FirstFing','AllFing','FingTransition','Chunk','Sequence'};
                            case 'fewComponents'
                                m = pcm_defineSequenceModels_seqType_simple(Seq,p,st);
                                modelNames = {'FirstFing','AllFing','Sequence'};
                        end
                        [M Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component','name','SeqTypeFeatures');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);               
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = zeros(size(D.seqType));
                        indx(D.seqType==st)=1; 
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept                        
                    end;
                    
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.seqType = ones(size(T.SN))*st;
                    T.sessN = ones(size(T.SN))*ss;
                    T=rmfield(T,{'reg','thetaCr'});
                    AllReg=addstruct(AllReg,T);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.seqType = ones(length(K.indx),1)*st;
                    K.sessN = ones(length(K.indx),1)*ss;
                    
                    KK=addstruct(KK,K);
                    fprintf('Done seqType %d/%d - reg %d/%d - sess %d/%d.\n',st,2,r,length(reg),2,ss,max(sessN));
                end
            end
        end
        fprintf('************************************* Done all %s. ***************************************\n\n\n\n\n',modelType);
        % save variables;
        dircheck(fullfile(pcmDir));
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_seqType_%s_%s.mat',modelType,parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_seqType_%s_%s.mat',modelType,parcelType)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_seqType_%s_%s.mat',modelType,parcelType)),'-struct','KK');
    
    case 'PCM_FingSeq_modelFamily'
        % only first finger, all fingers, transitions, seqType
        %(Feb 22nd 2019)
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        parcelType='162tessels'; % Brodmann or 162tessels
        reg = 1:8; % for tessels, 1:8 for Brodmann
        regSelect='all';
        sn=[5:9,11:31]; % change
        sessN=1:4;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi'})
        
        if strcmp(parcelType,'162tessels')
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            AllReg=[];
            KK=[];
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
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
                for p=1:length(sn)
                    glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p)&B.region==reg(r))};
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p)&B.region==reg(r))};
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p)&B.regoin==reg(r))}; % no intercept - use T.betaRAWint otherwise
                    end
                    partVec{p} = D.run;  % runs/partitions
                    %condVec{p} = D.seqNumb;
                    indx = ones(size(D.run));
                    Data{p} = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                    m = pcm_defineSequenceFingerModels_new(Seq,p);
                    [M,Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                    [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                    condVec{p} = repmat(Z,max(D.run),1);
                end;
                % fit models
                T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                T.roi = ones(size(T.SN))*reg(r);
                T.regType = ones(size(T.SN))*regType(r);
                T.regSide = ones(size(T.SN))*regSide(r);
                T.sessN = ones(size(T.SN))*ss;
                T=rmfield(T,{'reg','thetaCr'});
                AllReg=addstruct(AllReg,T);
                % calculations - posterior probability, knockIN/OUT -
                % removed
                K = pcm_calc(T.cross_likelihood,Comb);
                K.roi = ones(length(K.indx),1)*reg(r);
                K.regType = ones(length(K.indx),1)*regType(r);
                K.regSide = ones(length(K.indx),1)*regSide(r);
                K.sessN = ones(length(K.indx),1)*ss;
                KK=addstruct(KK,K);
                fprintf('Done sess:%d  region:%d/%d\n\n',ss,r,numel(reg));
            end
            % save variables;
            if ss==1
                save(fullfile(pcmDir,sprintf('ModelFamilyComb_SeqFing_%s.mat',parcelType)),'Comb');
            end
            save(fullfile(pcmDir,sprintf('ModelFamily_Fit_SeqFing_%s_sess-%d.mat',parcelType,ss)),'-struct','AllReg');
            save(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_%s_sess-%d.mat',parcelType,ss)),'-struct','KK');
            fprintf('Done sess-%d\n\n\n',ss);
        end
        

        %name = {'TS','US','FF','FA'};   
    case 'PCM_FingSeq_modelFamily_seqType'
        % only first finger, all fingers, transitions, seqType
        %(Feb 22nd 2019)
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        parcelType='162tessels'; % Brodmann or 162tessels
        reg = 1:8; % for tessels, 1:8 for Brodmann
        regSelect='all';
        sn=[5:9,11:31]; % change
        sessN=1:4;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi'})
        
        if strcmp(parcelType,'162tessels')
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            AllReg=[];
            KK=[];
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
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
                for st=1:2 % trained / untrained
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmSessDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p)&B.region==reg(r))};
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p)&B.region==reg(r))};
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p)&B.regoin==reg(r))}; % no intercept - use T.betaRAWint otherwise
                        end
                        partVec{p} = D.run(D.seqType==st);  % runs/partitions
                        %condVec{p} = D.seqNumb;
                        indx = ones(size(D.run));
                        Data{p} = beta(indx==1&D.seqType==st,:);  % Data is N x P (cond x voxels) - no intercept
                        m = pcm_defineSequenceFingerModels_seqType(Seq,p,st);
                        [M,Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        condVec{p} = repmat(Z,max(D.run),1);
                    end;
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.seqType = ones(size(T.SN))*st;
                    T.sessN = ones(size(T.SN))*ss;
                    T=rmfield(T,{'reg','thetaCr'});
                    AllReg=addstruct(AllReg,T);
                    % calculations - posterior probability, knockIN/OUT -
                    % removed
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.seqType = ones(length(K.indx),1)*st;
                    K.sessN = ones(length(K.indx),1)*ss;
                    KK=addstruct(KK,K);
                    fprintf('Done sess:%d sequence-%d region:%d/%d\n\n',ss,st,r,numel(reg));
                end
            end
            % save variables;
            if ss==1
                save(fullfile(pcmDir,sprintf('ModelFamilyComb_SeqFing_seqType_%s.mat',parcelType)),'Comb');
            end
            save(fullfile(pcmDir,sprintf('ModelFamily_Fit_SeqFing_seqType_%s_sess-%d.mat',parcelType,ss)),'-struct','AllReg');
            save(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_seqType_%s_sess-%d.mat',parcelType,ss)),'-struct','KK');
            fprintf('Done sess-%d\n\n\n',ss);
        end
        
    
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
    case 'PLOT_PCM_SeqFing_ROI_component'
        parcelType='Brodmann';
        roi=1:8;
        
        vararginoptions(varargin,{'parcelType','sessN','roi'});
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_%s.mat',parcelType)));
        
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
    case 'PLOT_PCM_SeqFing_seqType_ROI_component'
        parcelType='Brodmann';
        sessN=1:4;
        roi=1:8;
        vararginoptions(varargin,{'parcelType','sessN','roi'});
        
        A=[];
        for ss=sessN
            T = load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_seqType_%s_sess-%d.mat',parcelType,ss)));
            A=addstruct(A,T);
        end
        
        indx=1;
        for r = roi
            figure
            subplot(421);
            plt.bar([A.sessN>3 A.sessN],A.logBayes,'subset',A.roi==r&A.seqType==1,'split',[A.indx],'leg',{'Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s-trained',regname{r}));
            if indx==1
                ylabel('logBayes factor');
            else
                ylabel('');
            end
            subplot(422);
            plt.bar([A.sessN>3 A.sessN],A.logBayes,'subset',A.roi==r&A.seqType==2,'split',A.indx,'leg',{'Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s-untrained',regname{r}));
            if indx==1
                ylabel('logBayes factor');
            else
                ylabel('');
            end
            
            subplot(423);
            plt.bar([A.sessN>3 A.sessN],A.knockIN,'subset',A.roi==r&A.seqType==1,'split',A.indx,'leg',{'Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockIN factor - trained');
            else
                ylabel('');
            end
            
            subplot(424);
            plt.bar([A.sessN>3 A.sessN],A.knockIN,'subset',A.roi==r&A.seqType==2,'split',A.indx,'leg',{'Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockIN factor - untrained');
            else
                ylabel('');
            end
            
            subplot(425);
            plt.bar([A.sessN>3 A.sessN],A.knockOUT,'subset',A.roi==r&A.seqType==1,'split',A.indx,'leg',{'Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockOUT factor - trained');
            else
                ylabel('');
            end
            
            subplot(426);
            plt.bar([A.sessN>3 A.sessN],A.knockOUT,'subset',A.roi==r&A.seqType==2,'split',A.indx,'leg',{'Seq','First Fing','All Fing'},'leglocation','northeast');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('knockOUT factor - untrained');
            else
                ylabel('');
            end
            
            subplot(427);
            plt.bar([A.sessN>3 A.sessN],A.postProp,'subset',A.roi==r&A.seqType==1,'split',A.indx,'leg',{'Seq','First Fing','All Fing'},'leglocation','northeast');
            hold on;
            drawline(0.5,'dir','horz');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('Posterior probability - trained');
            else
                ylabel('');
            end
            
            subplot(428);
            plt.bar([A.sessN>3 A.sessN],A.postProp,'subset',A.roi==r&A.seqType==2,'split',A.indx,'leg',{'Trained Seq','Untrained Seq','First Fing','All Fing'},'leglocation','northeast');
            hold on;
            drawline(0.5,'dir','horz');
            title(sprintf('%s',regname{r}));
            if indx==1
                ylabel('Posterior probability - untrained');
            else
                ylabel('');
            end
            
           indx=indx+1;
        end
    case 'PLOT_PCM_SeqFing_ROI_logBayes'
        parcelType='Brodmann';
        roi=1:8;
        vararginoptions(varargin,{'parcelType','sessN','roi'});
        A=load(fullfile(pcmDir,sprintf('ModelFamily_Fit_SeqFing_%s.mat',parcelType)));
        
        for r=roi
            figure
            barplot(A.sessN,A.bayesEst(:,1:5),'subset',A.roi==r);
            ylabel('logBayes');
            title(regname{r});
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
    case 'PLOT_PCM_modelFamily_seqType' % this relevant (Feb 14 2019)
        %plots knock-IN / OUT / posterior probability
        reg=1;
        modelType = 'fewComponents'; % fewComponents or allComponents
        parcelType='Brodmann'; % Brodmann or 162tessels
        vararginoptions(varargin,{'reg','parcelType','modelType'});
        
        KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_seqType_%s_%s.mat',modelType,parcelType)));        
        nameFig = {'FF','AF','FT','Chu','Seq'};
        fig_num=0;
        for r = reg
            % knock-in
            figure(fig_num+1);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s);
                subplot(1,4,s)
                plt.bar(a.indx,a.knockIN,'split',a.seqType,'style',stySeq,'leg',{'trained','untrained'},'leglocation','north');
                
                if s==1
                    ylabel('Log Bayes - knock IN');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
                set(gca,'XTick',[1.5 4.5 8 11.5 15],'XTickLabel',nameFig);
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
                set(gca,'XTick',[1.5 4.5 8 11.5 15],'XTickLabel',nameFig);
                drawline(0,'dir','horz');
                plt.match('y');
            end
            % posterior probability
            figure(fig_num+3);
            for s=1:4
                a=getrow(KK,KK.roi==r&KK.sessN==s&KK.sn>3);
                subplot(1,4,s)
                plt.bar(a.indx,a.postProp,'split',a.seqType,'style',stySeq,'leg',{'trained','untrained'},'leglocation','north');             
                set(gca,'XTick',[1.5 4.5 8 11.5 15],'XTickLabel',nameFig);
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
                set(gca,'XTick',[1.5 4.5 8 11.5 15],'XTickLabel',nameFig);
                drawline(0,'dir','horz');
                plt.match('y');
                if s==1
                    ylabel('Log Bayes - component');
                    title(sprintf('Sess%d %s',s, regname{r}));
                else
                    ylabel('');
                    title(sprintf('Sess%d',s));
                end
            end
        end    
        % save combined structure of trained / untrained
       % save(fullfile(pcmDir,sprintf('ModelFamily_Stats_%s_bothSeqType.mat',parcelType)),'-struct','KK');
    case 'PLOT_SURFACE_PCM_tessels_seqType'
        parcelType='162tessels';
        sessN=1:4;
        thres=[0 1]; % [0 10] for knockIN [0 1 for knockOUT]
        featureSet=1:4;
        var='logBayes'; % logBayes, knockOUT, knockIN
        featureName={'FirstFing','AllFing','FingTrans','Chunk','Seq'};
        vararginoptions(varargin,{'parcelType','sessN','thres'});
        KK=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_seqType_%s_fingCount.mat',parcelType)));

        
        tesselNum=unique(KK.roi);

        for h=1:2;    
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            for ss=sessN
                for f=featureSet
                    RGBdata=zeros(size(C.data,1),3);
                    for r=tesselNum'
                        S=getrow(KK,KK.roi==r & KK.hemi==h & KK.sessN==ss & KK.indx==f);
                        % data from that hemi, tessel, session - per feature
                        
                        indx=C.index(C.data==r);    % indicate the right tessel
                        
                        RGBdata(indx,1)=mean(S.(var)(S.seqType==1)*(-1)); % trained
                        RGBdata(indx,3)=mean(S.(var)(S.seqType==2)*(-1)); % untrained
                        RGBdata(indx,2)=mean(S.(var)(S.seqType==1)*(-1))-mean(S.(var)(S.seqType==2)*(-1)); % difference between the two
                    end
                    RGBdata(RGBdata(:,1)<thres(1))=0;
                    RGBdata(RGBdata(:,3)<thres(1))=0;
                    
                    scale=[thres;thres;thres];
                    name={sprintf('%s_%s_sess%d',var,featureName{f},ss)};
                    R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                    caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.%s_%s_sess%d.RGB_paint',var,hem{h},featureName{f},ss)),R);
                    %caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.logBayes_%s_sess%d.RGB_paint',hem{h},featureName{f},ss)),R);
                end
            end
        end
    case 'PLOT_SURFACE_PCM_tessels'     
        parcelType='162tessels';
        sessN=1:4;
        thres=[0 3]; % [0 10] for knockIN [0 1 for knockOUT]
        var='logBayes'; %knockIN / knockOUT / logBayes
        featureSet=1:4;
         featureName={'Trained','Untrained','FirstFing','AllFing'};
        % colours: red - blue - green - green
        vararginoptions(varargin,{'parcelType','sessN','thres','var'});

        for h=1:2;    % only contralateral
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            for ss=sessN
                A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_%s_sess-%d.mat',parcelType,ss)));    
                tesselNum=unique(A.regType(A.regSide==h));
                for f=featureSet
                    RGBdata=zeros(size(C.data,1),3);
                    for r=tesselNum'
                        S=getrow(A,A.regType==r & A.regSide==h & A.sessN==ss & A.indx==f);
                        % data from that hemi, tessel, session - per feature
                        
                        indx=C.index(C.data==r);    % indicate the right tessel
                        TI=S.(var);
                        if f==1
                        %if f==2 % red
                            RGBdata(indx,1)=mean(TI); % trained
                        elseif f==2 % blue
                            RGBdata(indx,3)=mean(TI); % untrained
                        else % always green
                            RGBdata(indx,2)=mean(TI); % finger models
                        end
                    end
                    RGBdata(RGBdata(:,1)<thres(1))=0; % threshold
                    RGBdata(RGBdata(:,2)<thres(1))=0;
                    RGBdata(RGBdata(:,3)<thres(1))=0;
                    
                    scale=[thres;thres;thres];
                    % save
                    name={sprintf('PCM_%s_%s_sess%d',var,featureName{f},ss)};
                    R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                    caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.PCM_%s_%s_sess%d.RGB_paint',hem{h},var,featureName{f},ss)),R);
                end
            end
        end
    case 'PLOT_SURFACE_PCM_trained-untrained'
        parcelType='162tessels';
        sessN=1:4;
        thres=[0 3]; 
        var='logBayes'; %knockIN / knockOUT / logBayes
        vararginoptions(varargin,{'parcelType','sessN','thres','var'});
        
        for h=1:2;    % only contralateral
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            for ss=sessN
                A=load(fullfile(pcmDir,sprintf('ModelFamily_Stats_SeqFing_%s_sess-%d.mat',parcelType,ss)));
                tesselNum=unique(A.regType(A.regSide==h));
                data=zeros(size(C.data,1),3);
                for r=tesselNum'
                    S_tr=getrow(A,A.regType==r & A.regSide==h & A.sessN==ss & A.indx==1);
                    S_utr=getrow(A,A.regType==r & A.regSide==h & A.sessN==ss & A.indx==2);
                    % data from that hemi, tessel, session - per feature
                    indx=C.index(C.data==r);    % indicate the right tessel
                    data(indx,1)=mean(S_tr.(var)-S_utr.(var));
                end
                scale=[thres;thres;thres];
                % save
                column_name{1}=fullfile(sprintf('PCM_%s_trained-untrained_sess%d',var,ss));
              %  R=caret_struct('metric','data',data,'scales',{scale},'column_name',column_name);
                R=caret_struct('metric','data',data,'column_name',column_name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.PCM_%s_trained-untrained_sess%d.metric',hem{h},var,ss)),R);
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
        sn=[4:9,11:31];
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
    
    case 'CLUSTER_choose'               % CLUSTER start
        % choose based on significance
        sessN=1;
        betaChoice='multi';
        parcelType = '162tessels';
        vararginoptions(varargin,{'sessN','betaChoice','parcelType'});
        
        for ss=sessN       
            indx=[];
            T=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,ss)));
            for r=unique(T.region)'
                T1 = getrow(T,T.region==r);
                if sum(isnan(T1.dist_train))==0 % make sure all subjects have data there
                 %   [t,p]=ttestDirect(T.dist_train,[T.SN],1,'onesample','subset',T.region==r);
                 %   if (p<0.1 && t>0) % include into selection
                        indx=[indx;r];
                 %   end
                end
            end
            if ss==sessN(1)
                choice=indx;
            else
                choice=intersect(indx, choice);
            end
        end
        varargout={choice}; 
    case 'CLUSTER_choose_all'
        % choose all clusters where data present in all subj
        % use also for clustering
        sessN=1;
        betaChoice='multi';
        parcelType = '162tessels';
        vararginoptions(varargin,{'sessN','betaChoice','parcelType'});
        
        for ss=sessN       
            indx=[];
            T=load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,ss)));
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
    case 'CLUSTER_interSubj_RDM'        % RDM LEVEL - per tessel
        % estimate the consistency of RDMs in each tessel
        % both within and across subjects
        sn=[4:9,11:31];
        sessN=[1:4];
        parcelType='162tessels'; % Brodmann or 162tessels
        sessType = 'within';
        betaType='multi'; % multi,uni,RAW
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','betaType'});
        
        T = load(fullfile(betaDir,'group',sprintf('betas_partition_%sPW_%s',betaType,parcelType)));
        % split participants by groups
        S{1}=sn(mod(sn,2)==1);
        S{2}=sn(mod(sn,2)==0);
        for ss=sessN
            if strcmp(sessType,'within')
                T = getrow(T,T.sessN==ss);
            end
            RR=[];
            roi=unique(T.roi);
            
            if strcmp(parcelType,'Brodmann')
                regSide = [ones(1,8) ones(1,8)*2];
                regType = [1:8 1:8];
            elseif strcmp(parcelType,'162tessels')
                regSide = ones(size(roi));
                regSide(roi>158)=2;
                regType = roi;
                regType(regSide==2)=regType(regSide==2)-158;
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
            save(fullfile(distPscDir,sprintf('consist_interSubjRDM_%s_sess%d',parcelType,ss)),'-struct','RR');
        end
    case 'PLOT_interSubj_RDM'
        sessN=[1:4];
        parcelType='162tessels';
        vararginoptions(varargin,{'sessN','parcelType'});
        
        TT=[];
        for ss=sessN
            T = load(fullfile(distPscDir,sprintf('consist_interSubj_%s_sess%d',parcelType,ss)));
            T.sessN = ones(size(T.roi))*ss;
            TT=addstruct(TT,T);
        end
        figure
        plt.bar(TT.sessN,TT.corrDist,'subset',TT.sn1~=TT.sn2);
        xlabel('session'); ylabel('inter-subject RDM consistency');
        hold on;
        drawline(mean(TT.corrDist(TT.sn1==TT.sn2)),'dir','horz','linestyle','--');
    case 'PLOT_interSubj_RDM_surface'
        sessN=[1:4]; 
        parcelType='162tessels';
        type = {'within','between'};
        vararginoptions(varargin,{'sessN','parcelType'});
        
        for ss=sessN
            T = load(fullfile(distPscDir,sprintf('consist_interSubj_%s_sess%d',parcelType,ss)));
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
                        column_name{t} = fullfile(sprintf('consistRDM_%sSubj_sess%d.nii',type{t},ss));
                    end
                end
                R=caret_struct('metric','data',data,'column_name',column_name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.consistRDM_sess%d.metric',hem{h},ss)),R);
            end
        end
   
    case 'CLUSTER_calcAlpha'            % ALPHA LEVEL
        % calculate the alpha matrix
        % estimate the distance between RDMs across tessels
        sn=[4:9,11:31];
        sessN=[1:4];
        hemi=[1,2]; % 1- contra, 2 - ipsi, [1,2] - both
        parcelType='162tessels'; % Brodmann or 162tessels; or combined
        sessType = 'across'; % within or across
        seqType='all'; % all, trained, untrained
        crossType='all'; % all or splithalf
        distType='cosine'; % distance type
        betaType='multi'; % multi,uni,RAW
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','hemi','crossType','distType','seqType','betaType'});
        
        T = load(fullfile(betaDir,'group',sprintf('betas_partition_%sPW_%s',betaType,parcelType)));
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
            regNum = size(T1.roi,1);
            switch seqType
                case 'all'
                    t = [T1.partA_all;T1.partB_all];
                case 'trained' 
                    t = [T1.partA_train;T1.partB_train];
                case 'untrained'
                    t = [T1.partA_untrain;T1.partB_untrain];
            end
            a_all = squareform(pdist(t,distType)); % between all regions / parcels
           %A_all = rsa_calcdist(t,distType);
            switch crossType
                case 'all'
                    % take the mean of partA-A partB-B (not cross-part)
                    tmp(:,:,1) = a_all(1:regNum,1:regNum);
                    tmp(:,:,2) = a_all(regNum+1:end,regNum+1:end);
                    A = mean(tmp,3);
                case 'splithalf'
                    % consider distance of partA-B
                    tmp(:,:,1) = a_all(1:regNum,regNum+1:end);
                    tmp(:,:,2) = a_all(regNum+1:end,1:regNum);
                    A = mean(tmp,3);
            end
            % make confidence - 1-distance (within-region)
            % only really works for splithalf, otherwise 1
            Conf_all(:,:,s)=1-diag(A);    
            % set diagonal to 1 (for splithalf)
            % apply Gaussian similarity (for across region only)
            % create a matrix of similarity across regions 
            thres = quantile(rsa_vectorizeRDM(A),0.05);
            A     = exp(-A.^2 ./ (2*thres^2));
            A(1:regNum+1:end)=ones(size(A,1),1);
            A_all(:,:,s) = A;
        end
        A = nanmean(A_all,3);
        Conf=nanmean(Conf_all,3);
        Roi=T1.roi;
        RegType=T1.regType;
        RegSide=T1.regSide;
        sessN=T1.sessN;
       % save(fullfile(clusterDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)),'A','A_all','Conf','Conf_all','Roi','RegSide','RegType','sessN');
       switch sessType
           case 'within'
                save(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d_%s',seqType,distType,crossType,sessType,sessN(1),parcelType)),'A','A_all','Conf','Conf_all','Roi','RegSide','RegType','sessN');
           case 'across'
                save(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),parcelType)),'A','A_all','Conf','Conf_all','Roi','RegSide','RegType','sessN');
       end
    case 'CLUSTER_KLdivergence'
        % calculate the alpha matrix
        % using KL divergence
        sn=[4:9,11:31];
        sessN=4;
        hemi=[1,2]; % 1- contra, 2 - ipsi, [1,2] - both
        parcelType='162tessels'; % Brodmann or 162tessels; or combined
        sessType = 'within'; % within or across
        seqType='all'; % all, trained, untrained
        crossType='all'; % all or splithalf
        distType='cosine'; % distance type
        betaType='multi'; % multi,uni,RAW
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','hemi','crossType','distType','seqType','betaType'});
        
        T = load(fullfile(betaDir,'group',sprintf('betas_partition_%sPW_%s',betaType,parcelType)));
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
            regNum = size(T1.roi,1);
            switch seqType
                case 'all'
                    t = [T1.partA_all;T1.partB_all];
                case 'trained'
                    t = [T1.partA_train;T1.partB_train];
                case 'untrained'
                    t = [T1.partA_untrain;T1.partB_untrain];
            end
            a_all = squareform(pdist(t,distType)); % between all regions / parcels
            %A_all = rsa_calcdist(t,distType);
            switch crossType
                case 'all'
                    % take the mean of partA-A partB-B (not cross-part)
                    tmp(:,:,1) = a_all(1:regNum,1:regNum);
                    tmp(:,:,2) = a_all(regNum+1:end,regNum+1:end);
                    A = mean(tmp,3);
                case 'splithalf'
                    % consider distance of partA-B
                    tmp(:,:,1) = a_all(1:regNum,regNum+1:end);
                    tmp(:,:,2) = a_all(regNum+1:end,1:regNum);
                    A = mean(tmp,3);
            end
            % make confidence - 1-distance (within-region)
            % only really works for splithalf, otherwise 1
            Conf_all(:,:,s)=1-diag(A);
            % set diagonal to 1 (for splithalf)
            % apply Gaussian similarity (for across region only)
            % create a matrix of similarity across regions
            thres = quantile(rsa_vectorizeRDM(A),0.05);
            A     = exp(-A.^2 ./ (2*thres^2));
            A(1:regNum+1:end)=ones(size(A,1),1);
            A_all(:,:,s) = A;
        end
        A = nanmean(A_all,3);
        Conf=nanmean(Conf_all,3);
        Roi=T1.roi;
        RegType=T1.regType;
        RegSide=T1.regSide;
        sessN=T1.sessN;
        % save(fullfile(clusterDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)),'A','A_all','Conf','Conf_all','Roi','RegSide','RegType','sessN');
        switch sessType
            case 'within'
                save(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d_%s',seqType,distType,crossType,sessType,sessN(1),parcelType)),'A','A_all','Conf','Conf_all','Roi','RegSide','RegType','sessN');
            case 'across'
                save(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),parcelType)),'A','A_all','Conf','Conf_all','Roi','RegSide','RegType','sessN');
        end
    case 'CLUSTER_interSubj_old'
        parcelType='162tessels';
        sessN=[1:4];
        seqType='all'; % all, trained, untrained
        sn=[1:25];
        vararginoptions(varargin,{'var','clustN','parcelType','sessN','sessType','seqType'});
        
        CC=[];
        T=load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)));
        % try MDS plot
        [Y,l] =rsa_classicalMDS(T.A);
        figure
        scatterplot3(Y(:,1),Y(:,2),Y(:,3),'split',T.sessN,'markersize',5,'leg',{'sess1','sess2','sess3','sess4'});
        % loads as A
        M=zeros(numel(sn));
        for ss=sessN
            idx = find(T.sessN==ss);
            for s1=1:numel(sn)
                for s2=s1:numel(sn)
                    conn1 = rsa_vectorizeRDM(T.A_all(idx,idx,s1));
                    conn2 = rsa_vectorizeRDM(T.A_all(idx,idx,s2));
                    C.corr = corr(conn1',conn2');
                    C.sn1  = s1;
                    C.sn2  = s2;
                    C.sessN = ss;
                    CC=addstruct(CC,C);
                    M(s2,s1)=corr(conn1',conn2');
                end
            end
        end
        figure
        plt.bar(CC.sessN,CC.corr,'subset',CC.sn1~=CC.sn2);
        xlabel('session'); ylabel('inter-subject Alpha consistency');
    case 'CLUSTER_interSubj_Alpha'
        parcelType='162tessels';
        sessN=[1:4];
        seqType='all'; % all, trained, untrained
        distType='cosine';
        crossType='splithalf';% splithalf or all
        sessType='across';
        sn=[1:27];
        vararginoptions(varargin,{'var','clustN','parcelType','sessN','crossType','seqType'});
        
        CC=[];
        T = load(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),parcelType)));
        % loads as A
        M=zeros(numel(sn));
        for ss=sessN
            idx = find(T.sessN==ss);
            for s1=1:numel(sn)
                for s2=s1:numel(sn)
                    conn1 = rsa_vectorizeRDM(T.A_all(idx,idx,s1));
                    conn2 = rsa_vectorizeRDM(T.A_all(idx,idx,s2));
                    C.corr = corr(conn1',conn2');
                    C.sn1  = s1;
                    C.sn2  = s2;
                    C.sessN = ss;
                    CC=addstruct(CC,C);
                    M(s2,s1)=corr(conn1',conn2');
                end
            end
        end
        figure
        plt.bar(CC.sessN,CC.corr,'subset',CC.sn1~=CC.sn2);
        xlabel('session'); ylabel('inter-subject Alpha consistency');
    case 'CLUSTER_interSess_Alpha'
        parcelType='162tessels';
        sessN=[1:4];
        seqType='all'; % all, trained, untrained
        sn=[1:27];
        distType='cosine';
        crossType='splithalf';
        sessType='across';
        vararginoptions(varargin,{'crossType','distType','parcelType','sessN','sessType','seqType'});
        
        CC=[];
        T = load(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),parcelType)));        % loads as A
        for s=1:numel(sn)
            for s1=1:numel(sessN)
                idx1 = find(T.sessN==s1);
                for s2=s1:numel(sessN)
                    idx2 = find(T.sessN==s2);
                    conn1 = rsa_vectorizeRDM(T.A_all(idx1,idx1,s));
                    conn2 = rsa_vectorizeRDM(T.A_all(idx2,idx2,s));
                    C.corr = corr(conn1',conn2');
                    C.sn     = s;
                    C.sess1  = s1;
                    C.sess2  = s2;
                    CC=addstruct(CC,C);
                end
            end
        end
        figure
        subplot(211)
        plt.bar(CC.sess1,CC.corr,'subset',CC.sess2==CC.sess1+1);
        xlabel('session transitions'); title('inter-session Alpha consistency');
        ylabel('');
        subplot(212)
        plt.bar(CC.sess2,CC.corr,'subset',CC.sess1==1 & CC.sess2~=1);
        ylabel('');
    case 'PLOT_confidenceAlpha_surface' 
        % plot the confidence in splithalf dist per tessel
        parcelType='162tessels';
        sessN=[1:4];
        seqType='all'; % all, trained, untrained
        distType='cosine';
        crossType='splithalf';
        sessType='across';
        vararginoptions(varargin,{'crossType','distType','parcelType','sessN','sessType','seqType'});
        
        T = load(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),parcelType)));        % loads as A
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
                column_name{1} = fullfile(sprintf('confidence_alpha_sess%d_%sSeq.nii',ss,seqType));
                
                R=caret_struct('metric','data',data,'column_name',column_name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.confidenceAlpha_%sSeq_sess%d.metric',hem{h},seqType,ss)),R);
            end
        end
        
    case 'CLUSTER_make_group'           % CLUSTER LEVEL
        parcelType  = '162tessels'; % Brodmann or 162tessels; or combined
        seqType     = 'all'; % all, trained, untrained
        sessType    = 'within'; % within or across
        crossType   = 'all'; % all or splithalf
        distType    = 'cosine'; % distance type
        sessN       = 4;
        maxClust    = 10;
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','seqType','distType','crossType','maxClust','seqType'});
        switch sessType
            case 'within'
                T1 = load(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d_%s',seqType,distType,crossType,sessType,sessN(1),parcelType)));
            case 'across'
                T1 = load(fullfile(clusterDir,sprintf('alpha_%sSeq_%sDist_%s_%sSess-%d-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),parcelType)));
        end
        [C L U]=SpectralClustering(T1.A,maxClust,3);
       % Clust = kmeans(T1.A,maxClust);
        T1.community = community_louvain(T1.A);
        T1.cluster   = C;
        T1.laplace   = L;
        T1.eigenv    = U;
        switch sessType
            case 'within'
                save(fullfile(clusterDir,sprintf('cluster_%sSeq_%sDist_%s_%sSess-%d_clustN-%d_%s',seqType,distType,crossType,sessType,sessN(1),maxClust,parcelType)),'-struct','T1');
            case 'across'
                save(fullfile(clusterDir,sprintf('cluster_%sSeq_%sDist_%s_%sSess-%d-%d_clustN-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),maxClust,parcelType)),'-struct','T1');
        end
    case 'CLUSTER_make_subj' % TO DO
       % cluster separately for each subject
        parcelType='162tessels'; % Brodmann or 162tessels; or combined
        seqType='all'; % all, trained, untrained
        maxClust=7;
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','hemi','maxClust','seqType'});
        T1 = load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)));
        
        for s=1:size(T1.A_all,3);
            [C L U]=SpectralClustering(T1.A_all(:,:,s),maxClust,3);
            % Clust = kmeans(T1.A,maxClust);
            A{s}.A         = T1.A_all(:,:,s);
            A{s}.community = community_louvain(T1.A);
            A{s}.cluster   = C;
            A{s}.laplace   = L;
            A{s}.eigenv    = U;
            A{s}.Roi       = T1.Roi;
            A{s}.RegType   = T1.RegType;
            A{s}.RegSide   = T1.RegSide;
            A{s}.sessN     = T1.sessN;
        end
        save(fullfile(distPscDir,sprintf('consist_individ_crossval_%sSeq_%s_%d',seqType,parcelType,maxClust)),'A');
    
    
    case 'CLUSTER_extractROI_search'    % old
        % extract searchlight from ROI
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
    case 'CLUSTER_extractROI_dist'      % old
        % extract searchlight from ROI (non-crossval)
        sn=[4:9,11:28,30];
        sessN=[1:4];
        parcelType='Brodmann'; % Brodmann or 162tessels
        betaChoice='multi';
        sessType = 'within'; % within one session only / across sessions
        crossType = 'all'; % all the data or splithalf (even / odd)
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
                regSide(roi>158)=2;
                regType = roi;
                regType(regSide==2)=regType(regSide==2)-158;
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
            save(fullfile(clusterDir,sprintf('RDM_dist_%s_%sSess%d_%s',crossType,sessType,sessN,parcelType)),'-struct','TT');
        else
            save(fullfile(clusterDir,sprintf('RDM_dist_%s_%sSess%d-%d_%s',crossType,sessType,sessN(1),sessN(end),parcelType)),'-struct','TT');
        end

    case 'CLUSTER_acrossROI'        % REMOVE
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
    case 'CLUSTER_extractRDM'
        % extract RDMs and betas from each cluster (per subject /
        % hemisphere)
        sessType='within';
        sessN=4;
        distType='cosine';
        clustN = 10;
        parcelType='162tessels';
        seqType='trained';
        sn=[4:9,11:31];
        crossType='all'; % all or splithalf
        betaType='multi';
        vararginoptions(varargin,{'sessType','sessN','distType','clustN','crossType','betaType'});
        
         switch sessType
            case 'within'
                R=load(fullfile(clusterDir,sprintf('cluster_%sSeq_%sDist_%s_withinSess-%d_clustN-%d_%s',seqType,distType,crossType,sessN,clustN,parcelType)));
             case 'across'
            %    R=load(fullfile(distPscDir,sprintf('clusterResults_%s_%sDist_%dclusters_%sSeq_%sSess%d-%d',parcelType,distType,clustN,seqType,sessType,sessN(1),sessN(end))));
            %    T=load(fullfile(distPscDir,sprintf('cluster_RDM_dist_%sSess%d-%d_%s',sessType,sessN(1),sessN(end),parcelType)));
         end
         DD=[];
         for ss=sessN
             % get betas and rdms
             B = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss)));
             S = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaType,ss)));
             for s=sn
                 for c=1:clustN % extract data for each cluster
                     C = getrow(R,R.cluster==c);
                     % determine if regions span both hemispheres
                     hemiIdx = unique(B.regSide(ismember(B.region,C.Roi)));
                     for h=hemiIdx'
                         B1 = getrow(B,ismember(B.region,C.Roi) & B.SN==s & B.regSide==h);
                         S1 = getrow(S,ismember(S.region,C.Roi) & S.SN==s & S.regSide==h);
                         C_h = getrow(C,C.RegSide==h);
                         % concatenate betas into one region
                         beta=[];
                         switch seqType
                             case 'trained'
                                 D.rdm=nanmean(S1.RDM_train,1);
                             case 'untrained'
                                 D.rdm=nanmean(S1.RDM_untrain,1);
                             case 'all'
                                 D.rdm=nanmean(S1.RDM,1);
                         end
                         for i=1:length(unique(C_h.Roi))
                             if sum(sum(isnan(B1.betaW{i})))==0
                                 beta = [beta B1.betaW{i}];
                             end
                         end
                         D.beta = {beta};
                         D.cluster = c;
                         D.sn = s;
                         D.hemi = h;
                         DD=addstruct(DD,D);
                     end
                 end
                 fprintf('Done subject %d/%d\n',find(s==sn),length(sn));
             end
         end

         switch sessType
             case 'within'
                 save(fullfile(clusterDir,sprintf('clusterData_%sSeq_%sDist_%s_withinSess-%d_clustN-%d_%s',seqType,distType,crossType,sessN,clustN,parcelType)),'-struct','DD');
             case 'between'
                 save(fullfile(clusterDir,sprintf('clusterData_%sDist-%dcluster_%sSess_%d-%d',distType,clustN,sessType,sessN(1),sessN(end))),'-struct','DD');
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
        
    case 'CLUSTER_color' % assign colors, dendogram
        parcelType  = '162tessels'; % Brodmann or 162tessels; or combined
        seqType     = 'all'; % all, trained, untrained
        sessType    = 'within'; % within or across
        crossType   = 'all'; % all or splithalf
        distType    = 'cosine'; % distance type
        sessN       = 4;
        maxClust    = 10;
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','seqType','distType','crossType','maxClust','seqType'});
        switch sessType
            case 'within'
                T = load(fullfile(clusterDir,sprintf('cluster_%sSeq_%sDist_%s_%sSess-%d_clustN-%d_%s',seqType,distType,crossType,sessType,sessN(1),maxClust,parcelType)));
            case 'across'
                T = load(fullfile(clusterDir,sprintf('cluster_%sSeq_%sDist_%s_%sSess-%d-%d_clustN-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),maxClust,parcelType)));
        end
        % Evaluate similarity graph (W)
        [Csort,idx]=sort(T.cluster);
        A = full(T.A);
        Asort=A(idx,:);
        Asort=Asort(:,idx);
        Ncluster = numel(unique(T.cluster));
        
        % Laplacian eigenvector
        Lsort = T.laplace(idx,:); Lsort=Lsort(:,idx);
        Usort = T.eigenv(idx,:);
        
        % Merge similarities according to clusters
        Asort(logical(eye(size(Asort)))) = NaN;
        for i=1:Ncluster
            for j=1:Ncluster
                idxi=Csort==i;
                idxj=Csort==j;
                rWsrot(i,j) = nanmean(vec(Asort(idxi,idxj)));
                rLsort(i,j) = nanmean(vec(Lsort(idxi,idxj)));
            end
            rUsort(i,:) = nanmean(Usort(idxi,:),1);
        end
        
        % inspect sorted diagram
        figure
        imagesc(Asort);
        % get dendogram linkage
        Z = linkage(real(rUsort),'ward','euclidean');
        Col = colorDendrogram(Z,size(rUsort,1),'colorspace','rgb','order',[1,2,3],'fig',0,'weight',1);

        figure
        [h,t,per] = dendrogram(Z,'Orientation','right');
        hold on;
        xlim=get(gca,'xlim');
        for i=1:numel(per)
            plot(xlim(1),i,'o','markersize',8,...
                'markeredgecolor','k',...
                'markerfacecolor',Col(per(i),:));
            hold on;
        end
        for i=1:numel(h)
            set(h(i),'color',[0 0 0],'linewidth',1);
        end
        cd(clusterDir);
         switch sessType
            case 'within'
                dlmwrite(sprintf('my%dColors_%sSess-%d_%sSeq_%s.txt',Ncluster,sessType,sessN(1),seqType,crossType),Col);
            case 'across'
                dlmwrite(sprintf('my%dColors_%sSess-%d-%d_%sSeq_%s.txt',Ncluster,sessType,sessN(1),sessN(end),seqType,crossType),Col);
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
        figure
        for s=1:numel(sn)
            fprintf('Subject %d/%d\n',s,numel(sn));
            T1 = getrow(T,T.sn==sn(s));
            switch seqType
                case 'all'
                    WA = corrN(T1.partA_all',T1.partB_all');
                    WB = corrN(T1.partB_all',T1.partA_all');
                case 'trained' 
                    WA = corrN(T1.partA_train',T1.partB_train');
                    WB = corrN(T1.partB_train',T1.partA_train');
                case 'untrained'
                    WA = corrN(T1.partA_untrain',T1.partB_untrain');
                    WB = corrN(T1.partB_untrain',T1.partA_untrain');
            end
            W = (WA + WB)./2;
            Conf_all(:,:,s)=diag(W);
            W2 = 1-W;
            % alternatively here make diagonal 0, apply gaussian similarity
            W3=W2;
            W3(1:size(W,1)+1:end)=zeros(size(W,1),1);
            thres = quantile(rsa_vectorizeRDM(W3),0.05);
            A = exp(-W3.^2./(2*thres^2));
          %  W3 = W2/2;
          %  W4 = 1-W3;
          %  W5 = W4;
          %  W5(1:size(W,1)+1:end)=ones(size(W,1),1);
          %%%  W_sim = W./max(max(W));
          %%%  W_sim = W./2;
          %%%  W_sim = 1-W_sim;            
            %A_all(:,:,s) = W5;
            A_all(:,:,s) = A;
            subplot(5,5,s)
            imagesc(A);
        end
        A = nanmean(A_all,3);
        Conf=nanmean(Conf_all,3);
        Roi=T1.roi;
        RegType=T1.regType;
        RegSide=T1.regSide;
        sessN=T1.sessN;
        save(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s',seqType,parcelType)),'A','A_all','Conf','Conf_all','Roi','RegSide','RegType','sessN');
    case 'CLUSTER_svd'
        % estimate the consistency of RDMs in each tessel
        % both within and across subjects
        sn=[4:9,11:28,30];
        sessN=[1:4];
        hemi=[1,2]; % 1- contra, 2 - ipsi, [1,2] - both
        parcelType='162tessels'; % Brodmann or 162tessels; or combined
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'sn','sessN','parcelType','sessType','hemi','maxClust','seqType'});
        
        T = load(fullfile(betaDir,'group',sprintf('betas_partition_%s',parcelType))); 
        for ss=sessN
            for s=1:numel(sn)
                fprintf('Subject %d/%d\n',s,numel(sn));
                T1 = getrow(T,T.sn==sn(s)&T.sessN==ss);
                switch seqType
                    case 'all'
                        W = corrN(T1.partA_all',T1.partB_all');
                    case 'trained'
                        W = corrN(T1.partA_train',T1.partB_train');
                    case 'untrained'
                        W = corrN(T1.partA_untrain',T1.partB_untrain');
                end
                W2 = 1-W;
                % alternatively here make diagonal 0, apply gaussian similarity
                W3=W2;
                W3(1:size(W,1)+1:end)=zeros(size(W,1),1);
                thres = quantile(rsa_vectorizeRDM(W3),0.05);
                A = exp(-W3.^2./(2*thres^2));
                A_all{ss}(s,:) = rsa_vectorizeRDM(A);
            end
        end
        keyboard;
           
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
        for ss=sessN
            if strcmp(sessType,'within')
                T = getrow(T,T.sessN==ss);
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
            save(fullfile(distPscDir,sprintf('consist_interSubj_%s_sess%d',parcelType,ss)),'-struct','RR');
        end
    case 'PLOT_consist'
        sessN=4; 
        parcelType='162tessels';
        type = {'within','between'};
        vararginoptions(varargin,{'sessN','parcelType'});
        
        for ss=sessN
            T = load(fullfile(distPscDir,sprintf('consist_interSubj_%s_sess%d',parcelType,ss)));
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
                        column_name{t} = fullfile(sprintf('consistRDM_%sSubj_sess%d.nii',type{t},ss));
                    end
                end
                R=caret_struct('metric','data',data,'column_name',column_name);
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.consistRDM_sess%d.metric',hem{h},ss)),R);
            end
        end
        
    case 'CLUSTER_surface'
        distType='cosine';
        parcelType='162tessels';
        sessN=[1:2];
        maxClust=10;
        sessType='within';
        seqType='all'; % all, trained, untrained
        var='cluster';
        crossType='all'; % all or splithalf
        vararginoptions(varargin,{'distType','maxClust','parcelType','sessN','sessType','seqType','crossType'});
        
        switch sessType
            case 'within'
                T = load(fullfile(clusterDir,sprintf('cluster_%sSeq_%sDist_%s_%sSess-%d_clustN-%d_%s',seqType,distType,crossType,sessType,sessN(1),maxClust,parcelType)));
                clustN = length(unique(T.(var)));
                colMap = dlmread(sprintf('my%dColors_%sSess-%d_%sSeq_%s.txt',clustN,sessType,sessN(1),seqType,crossType));
            case 'across'
                T = load(fullfile(clusterDir,sprintf('cluster_%sSeq_%sDist_%s_%sSess-%d-%d_clustN-%d_%s',seqType,distType,crossType,sessType,sessN(1),sessN(end),maxClust,parcelType)));
                clustN = length(unique(T.(var)));
                colMap = dlmread(sprintf('my%dColors_%sSess-%d-%d_%sSeq_%s.txt',clustN,sessType,sessN(1),sessN(end),seqType,crossType));
        end

        for h=1:2
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            RGBdata=zeros(size(C.data,1),3);
            for ss=sessN
                for k=1:clustN
                    S=getrow(T,T.RegSide==h & T.(var)==k & T.sessN==ss);
                    % data from that hemi, tessel, per cluster
                    indx=C.index(ismember(C.data,S.RegType));    % indicate the right tessel
                    RGBdata(indx,1)=colMap(k,1);
                    RGBdata(indx,2)=colMap(k,2);
                    RGBdata(indx,3)=colMap(k,3);
                    scale=[1;1;1];
                end
                scale(1:clustN,1)=0;
                scale(1:clustN,2)=1;
                
                name={sprintf('clusters-%d_%sSeq_%sSess%d_%s',clustN,seqType,sessType,ss,crossType)};
                R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.clusters-%d_%sSeq_%sSess%d_%s.RGB_paint',hem{h},clustN,seqType,sessType,ss,crossType)),R);
            end
        end    
    case 'CLUSTER_surface_crossval'
        var='cluster'; % community or cluster
        parcelType='162tessels';
        sessN=[1:4];
        clustN=10;
        sessType='across';
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'var','clustN','parcelType','sessN','sessType','seqType'});
        
        T=load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));

        clustN = length(unique(T.(var)));
        colMap = [155 0 95; 10 41 92; 26 169 166; 152 95 153; 230 57 70;...
                    17 86 95; 193 100 94; 44 98 99; 32 164 243; 249 189 205];
        colMap = colMap./255;
        
        for h=1
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            RGBdata=zeros(size(C.data,1),3);
            for ss=sessN
                for k=1:clustN
                    S=getrow(T,T.RegSide==h & T.(var)==k & T.sessN==ss);
                    % data from that hemi, tessel, per cluster
                    indx=C.index(ismember(C.data,S.RegType));    % indicate the right tessel
                    RGBdata(indx,1)=colMap(k,1);
                    RGBdata(indx,2)=colMap(k,2);
                    RGBdata(indx,3)=colMap(k,3);
                    scale=[1;1;1];
                end
                scale(1:clustN,1)=0;
                scale(1:clustN,2)=1;
                
                name={sprintf('clusters%d_%s_%sSeq_%sSess%d',clustN,var,seqType,sessType,ss)};
                R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                
                caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.clusters-%d_%s_%sSeq_%sSess%d.RGB_paint',hem{h},clustN,var,seqType,sessType,ss)),R);
            end
        end
    case 'CLUSTER_surface_crossval_individ'
        % separately for each individual
        var='cluster'; % community or cluster
        parcelType='162tessels';
        sessN=[1:4];
        clustN=10;
        sessType='across';
        seqType='all'; % all, trained, untrained
        sn=[1:25];
        vararginoptions(varargin,{'var','clustN','parcelType','sessN','sessType','seqType'});
        
        load(fullfile(distPscDir,sprintf('consist_individ_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));
        % loads as A
        
        clustN = length(unique(A{1}.(var)));
        colMap = [155 0 95; 10 41 92; 26 169 166; 152 95 153; 230 57 70;...
            17 86 95; 193 100 94; 44 98 99; 32 164 243; 249 189 205];
        colMap = colMap./255;
        
        for h=1:2
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            RGBdata=zeros(size(C.data,1),3);
            for ss=sessN
                for s=sn
                    for k=1:clustN
                        T=A{s};
                        S=getrow(T,T.RegSide==h & T.(var)==k & T.sessN==ss);
                        % data from that hemi, tessel, per cluster
                        indx=C.index(ismember(C.data,S.RegType));    % indicate the right tessel
                        RGBdata(indx,1)=colMap(k,1);
                        RGBdata(indx,2)=colMap(k,2);
                        RGBdata(indx,3)=colMap(k,3);
                        scale=[1;1;1];
                    end
                    scale(1:clustN,1)=0;
                    scale(1:clustN,2)=1;
                    
                    name={sprintf('clusters%d_%s_%sSeq_%sSess%d_s%d',clustN,var,seqType,sessType,ss,s)};
                    R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                    
                    caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.clusters-%d_%s_%sSeq_%sSess%d_s%d.RGB_paint',hem{h},clustN,var,seqType,sessType,ss,s)),R);
                    fprintf('Done sess-%d subject-%d\n',ss,s);
                end
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
    
    case 'CLUSTER_sort'
        var='Clust'; % cluster or community
        parcelType='162tessels';
        clustN=9;
        hemi=[1,2];
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'var','clustN','parcelType','sessType','seqType','hemi'});
        
        T=load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));
        
        T1 = getrow(T,ismember(T.RegSide,hemi));
        % Evaluate similarity graph (W)
        [Asort,idx]=sort(T1.(var));
        W = full(T.A);
        Wsort=W(idx,:);
        Wsort=Wsort(:,idx);
        figure
        imagesc(Wsort);
    case 'CLUSTER_sort_subj'
        var='Clust'; % cluster or community
        parcelType='162tessels';
        clustN=9;
        hemi=[1,2];
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'var','clustN','parcelType','sessType','seqType','hemi'});
        
        load(fullfile(distPscDir,sprintf('consist_individ_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));
        sn=1:size(A,2);
        figure
        for s=sn
            T=A{s};
            T1 = getrow(T,ismember(T.RegSide,hemi));
            % Evaluate similarity graph (W)
            [Asort,idx]=sort(T1.(var));
            W = full(T.A);
            Wsort=W(idx,:);
            Wsort=Wsort(:,idx);
            subplot(5,5,s);
            imagesc(Wsort);
        end
    case 'CLUSTER_hierarchy'
        % calculates similarity within / across clusters for group
        var='Clust'; % Clust or Community
        parcelType='162tessels';
        sessN=[1:4];
        clustN=9;
        hemi=1;
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'var','clustN','parcelType','sesN','sessType','seqType','hemi'});
        
        T=load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));
        
        AA=[];
        for ss=sessN
            T1 = getrow(T,T.sessN==ss&ismember(T.RegSide,hemi));
            indx = find(T.sessN==ss&ismember(T.RegSide,hemi));
            T1.A = T.A(indx,indx);
            T1.A(1:size(T1.A,1)+1:end) = ones(size(T1.A,1),1)*NaN;
            clust = unique(T1.(var));
            M{ss} = zeros(length(clust));
            for c1=1:length(clust)
                idx1=find(T1.(var)==clust(c1)); % index for all tessels in that cluster
                if size(idx1,1)>1
                    for c2=c1:length(clust)
                        idx2=find(T1.(var)==clust(c2)); % all other tessels
                        dist11 = T1.A(idx1,idx1);
                        dist22 = T1.A(idx2,idx2);
                        dist12 = T1.A(idx1,idx2);
                        M{ss}(c1,c1) = nanmean(dist11(:));
                        M{ss}(c2,c2) = nanmean(dist22(:));
                        M{ss}(c1,c2) = nanmean(dist12(:));
                        M{ss}(c2,c1) = nanmean(dist12(:));
                    end
                end
            end
        end
        keyboard;
    case 'CLUSTER_hierarchy_acrossSess'
          % calculates similarity within / across clusters for group
        var='cluster'; % cluster or community
        parcelType='162tessels';
        clustN=9;
        hemi=[1,2];
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'var','clustN','parcelType','sessType','seqType','hemi'});
        
        T=load(fullfile(distPscDir,sprintf('consist_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));
        
        AA=[];
        T1 = getrow(T,ismember(T.RegSide,hemi));
        indx = find(ismember(T.RegSide,hemi));
        T1.A = T.A(indx,indx);
        T1.A(1:size(T1.A,1)+1:end) = ones(size(T1.A,1),1)*NaN;
        for c=unique(T1.(var))'
            idxc=find(T1.(var)==c); % index for all tessels in that cluster
            if size(idxc,1)>1
                idxnc=find(T1.(var)~=c); % all other tessels
                distClust = T1.A(idxc,idxc);
                distOut   = T1.A(idxc,idxnc);
                t1=rsa_vectorizeRDM(distClust);
                t2=distOut(:);
                A.similarity(1,:)    = nanmean(t1);
                A.similarity(2,:)    = nanmean(t2);
                A.withinOut          = [1;2];
                A.numTesselIn        = repmat(length(idxc),2,1);
                A.clustN             = repmat(c,2,1);
                AA=addstruct(AA,A);
            end
            clear idxc idxnc distClust distOut;
        end
        keyboard;
    case 'CLUSTER_hierarchy_perSubj'
           % calculates similarity within / across clusters for group
        var='cluster'; % cluster or community
        parcelType='162tessels';
        clustN=9;
        hemi=[1,2];
        seqType='all'; % all, trained, untrained
        vararginoptions(varargin,{'var','clustN','parcelType','sessType','seqType','hemi'});
        
        load(fullfile(distPscDir,sprintf('consist_individ_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));
        
        RR=[];
        sn=1:size(A,2);
        for s=sn
        T = A{s};
        T1 = getrow(T,ismember(T.RegSide,hemi));
        indx = find(ismember(T.RegSide,hemi));
        T1.A = T.A(indx,indx);
        T1.A(1:size(T1.A,1)+1:end) = ones(size(T1.A,1),1)*NaN;
        for c=unique(T1.(var))'
            idxc=find(T1.(var)==c); % index for all tessels in that cluster
            if size(idxc,1)>1
                idxnc=find(T1.(var)~=c); % all other tessels
                distClust = T1.A(idxc,idxc);
                distOut   = T1.A(idxc,idxnc);
                t1=rsa_vectorizeRDM(distClust);
                t2=distOut(:);
                R.similarity(1,:)    = nanmean(t1);
                R.similarity(2,:)    = nanmean(t2);
                R.withinOut          = [1;2];
                R.numTesselIn        = repmat(length(idxc),2,1);
                R.clustN             = repmat(c,2,1);
                R.sn                 = repmat(s,2,1);
                RR=addstruct(RR,R);
            end
            clear idxc idxnc distClust distOut;
        end
        end
        figure
        plt.bar(RR.clustN,RR.similarity,'split',RR.withinOut);

   
    case 'CLUSTER_surface_crossval_individ_old' % TO DO
        % separately for each individual
        var='cluster'; % community or cluster
        parcelType='162tessels';
        sessN=[1:4];
        clustN=10;
        sessType='across';
        seqType='all'; % all, trained, untrained
        sn=[1:25];
        vararginoptions(varargin,{'var','clustN','parcelType','sessN','sessType','seqType'});
        
        load(fullfile(distPscDir,sprintf('consist_individ_crossval_%sSeq_%s_%d',seqType,parcelType,clustN)));
        % loads as A
        
        clustN = length(unique(A{1}.(var)));
        colMap = [155 0 95; 10 41 92; 26 169 166; 152 95 153; 230 57 70;...
            17 86 95; 193 100 94; 44 98 99; 32 164 243; 249 189 205];
        colMap = colMap./255;
        
        for h=1:2
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
            % per hemisphere
            RGBdata=zeros(size(C.data,1),3);
            for ss=sessN
                for s=sn
                    for k=1:clustN
                        T=A{s};
                        S=getrow(T,T.RegSide==h & T.(var)==k & T.sessN==ss);
                        % data from that hemi, tessel, per cluster
                        indx=C.index(ismember(C.data,S.RegType));    % indicate the right tessel
                        RGBdata(indx,1)=colMap(k,1);
                        RGBdata(indx,2)=colMap(k,2);
                        RGBdata(indx,3)=colMap(k,3);
                        scale=[1;1;1];
                    end
                    scale(1:clustN,1)=0;
                    scale(1:clustN,2)=1;
                    
                    name={sprintf('clusters%d_%s_%sSeq_%sSess%d_s%d',clustN,var,seqType,sessType,ss,s)};
                    R=caret_struct('RGBpaint','data',RGBdata,'scales',{scale},'column_name',name);
                    
                    caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.clusters-%d_%s_%sSeq_%sSess%d_s%d.RGB_paint',hem{h},clustN,var,seqType,sessType,ss,s)),R);
                    fprintf('Done sess-%d subject-%d\n',ss,s);
                end
            end
        end
    case 'CLUSTER_consist_old' % old - make a community projection instead
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
    
    case 'TRANSFORM_rdms'
        % use this case to transform G from reg1->reg2
        % save in a subfolder of the clusterDir
        sn=[4:9,11:31];
        parcelType='Brodmann'; %do combination - Brodmann + BG-striatum + thalamus
        betaChoice='multi'; % multi, uni
        sessN=[1:4];
        seqType='all'; %all, trained, untrained
        vararginoptions(varargin,{'sn','parcelType','sessN','betaChoice','seqType'});
        
        TT=[];
        q=zeros(size(sessN));
        for ss=sessN
            if ~iscell(parcelType)
                S = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,ss)));
            else
                for i=1:size(parcelType,2)
                    Struct{i} = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType{i},betaChoice,ss)));
                end
                    S = sml1_imana_dist('combine_reg','Struct',Struct);
            end
            % pre-allocate alpha matrix
            A{ss}=zeros(length(unique(S.region)),length(unique(S.region)),length(unique(sn)));
            for s=1:length(sn)
                S1 = getrow(S,S.SN==sn(s));
                for r1=unique(S.region)'
                    for r2=unique(S.region)'
                        % extract Gs (square)
                        IPM1=S1.IPM(S1.region==r1,:);
                        IPM2=S1.IPM(S1.region==r2,:);
                      %  IPM1=S1.RDM(S1.region==r1,:);
                      %  IPM2=S1.RDM(S1.region==r2,:);
                        G1=rsa_squareIPM(IPM1);
                        G2=rsa_squareIPM(IPM2);
                       % G1=rsa_squareRDM(IPM1);
                       % G2=rsa_squareRDM(IPM2);
                       % new
                       H=eye(12)-ones(12,12)./12;  % centering matrix!
                       G1 = H*G1*H';  % double centered G matrix - rows and columns
                       G2 = H*G2*H'; 
                       IPM1 = rsa_vectorizeIPM(G1);
                       IPM2 = rsa_vectorizeIPM(G2);
                       %IPM1=diag(G1)';
                       %IPM2=diag(G2)';
                        % extract subsets of Gs for trained / untrained
                        if strcmp(seqType,'trained')
                            T.G1 = {G1(1:6,1:6)};
                            T.G2 = {G2(1:6,1:6)};
                        elseif strcmp(seqType,'untrained')
                            T.G1 = {G1(7:12,7:12)};
                            T.G2 = {G2(7:12,7:12)};
                        else
                            T.G1 = {G1};
                            T.G2 = {G2};
                        end
                        % calculate cosine and KL divergence between G1/2
                        T.cosDist = pdist([IPM1;IPM2],'cosine'); 
                    %    T.KL = KLdivergence(T.G1{1},T.G2{1});
                        % calculate transformation matrix G1 -> G2
                        [Trans,predG]=calcTransformG(T.G1{1},T.G2{1});
                        T.T={round(Trans,3)}; % round
                        T.predG={predG};
                        T.sn=sn(s);
                        T.sn_idx=s;
                        T.sessN=ss;
                        T.reg1=r1;
                        T.reg2=r2;
                        T.regType1=S1.regType(S1.region==r1);
                        T.regType2=S1.regType(S1.region==r2);
                        T.regSide1=S1.regSide(S1.region==r1);
                        T.regSide2=S1.regSide(S1.region==r2);
                        TT=addstruct(TT,T);
                        % construct alpha matrix with cosine distances
                        A_subj{ss}(r1,r2,s)=T.cosDist;
                    end
                end
                fprintf('Done sess-%d\tsubj-%d/%d\n',ss,s,length(unique(sn)));
            end
            % save 0.05 for constructing similarity matrix
            A_group{ss}=nanmean(A_subj{ss},3);
            q(ss)=quantile(rsa_vectorizeRDM(A_group{ss}),0.05);
        end
        % make a similarity alpha matrix - using the same threshold
        thres = max(q);
        for i=1:numel(sessN)
            A{i} = exp(-A_group{i}.^2./(2*thres^2));
            %A{i} = 1-A_group{i};
        end
        if ~iscell(parcelType)
            save(fullfile(transformDir,sprintf('alpha_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)),'A');
            save(fullfile(transformDir,sprintf('transform_cosine_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)),'-struct','TT');
        else
            numParcel = size(parcelType,2);
            for i=1:numParcel
                if i==1
                    nameParcel = sprintf('%s',parcelType{i});
                else
                    nameParcel = sprintf('%s-%s',nameParcel,parcelType{i});
                end
            end
            save(fullfile(transformDir,sprintf('alpha_%s_%sSeq_%sPW',nameParcel,seqType,betaChoice)),'A');
            save(fullfile(transformDir,sprintf('transform_cosine_%s_%sSeq_%sPW',nameParcel,seqType,betaChoice)),'-struct','TT');
        end
    case 'TRANSFORM_sess'
        sn=[4:9,11:31];
        parcelType='Brodmann'; %do combination - Brodmann + BG-striatum + thalamus
        betaChoice='multi'; % multi, uni
        sessN=[1:4];
        seqType='all'; %all, trained, untrained
        vararginoptions(varargin,{'sn','parcelType','sessN','betaChoice','seqType'});
        
        TT=[];
        for s1=1:3 % session
            S1 = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,s1)));
            S2 = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,s1+1)));
            for s=1:length(sn)
                K1 = getrow(S1,S1.SN==sn(s));
                K2 = getrow(S2,S2.SN==sn(s));
                for r=unique(K1.region)'
                    % extract RDMs
                    switch seqType
                        case 'trained'
                            dist1 = K1.RDM_train(r,:);
                            dist2 = K2.RDM_train(r,:);
                        case 'untrained'
                            dist1 = K1.RDM_untrain(r,:);
                            dist2 = K2.RDM_untrain(r,:);
                        case 'all'
                            dist1 = K1.RDM(r,:);
                            dist2 = K2.RDM(r,:);
                    end
                    % calculate cosine and KL divergence between dist1/2
                    T.cosDist = pdist([dist1;dist2],'cosine');
                    %    T.KL = KLdivergence(T.G1{1},T.G2{1});
                    % calculate transformation matrix G1 -> G2
                    %  [Trans,predG]=calcTransformG(T.G1{1},T.G2{1});
                    %  T.T={round(Trans,3)}; % round
                    %  T.predG={predG};
                    T.sn=sn(s);
                    T.sn_idx=s;
                    T.sess1=s1;
                    T.sess2=s1+1;
                    T.sessTrans = s1;
                    T.reg=r;
                    T.regType=K1.regType(K1.region==r);
                    T.regSide=K1.regSide(K1.region==r);
                    TT=addstruct(TT,T);    
                end
            end
        end
        
        save(fullfile(transformDir,sprintf('transform_sessTrans_cosine_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)),'-struct','TT');

    case 'TRANSFORM_uni'
        % use this case to transform G from reg1->reg2
        % save in a subfolder of the clusterDir
        sn=[4:9,11:31];
        parcelType='Brodmann'; %do combination - Brodmann + BG-striatum + thalamus
        betaChoice='multi'; % multi, uni
        sessN=[1:4];
        seqType='all'; %all, trained, untrained
        vararginoptions(varargin,{'sn','parcelType','sessN','betaChoice','seqType'});
        
        TT=[];
        condVec = repmat([1:12]',8,1);
        X=indicatorMatrix('identity_p',condVec);
        q=zeros(size(sessN));
        for ss=sessN
            
            S = load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d',parcelType,ss)));
            
            % pre-allocate alpha matrix
            A{ss}=zeros(length(unique(S.region)),length(unique(S.region)),length(unique(sn)));
            for s=1:length(sn)
                S1 = getrow(S,S.SN==sn(s));
                for r1=unique(S.region)'
                    for r2=unique(S.region)'
                        % extract Gs (square)
                        t1=S1.betaW{r1};
                        t2=S1.betaW{r2};
                        d1=nanmean(pinv(X)*t1(1:size(X,1),:),2);
                        d2=nanmean(pinv(X)*t2(1:size(X,1),:),2);
                        % calculate cosine and KL divergence between G1/2
                        T.cosDist = corr(d1,d2);
                        T.sn=sn(s);
                        T.sn_idx=s;
                        T.sessN=ss;
                        T.reg1=r1;
                        T.reg2=r2;
                        T.regType1=S1.regType(S1.region==r1);
                        T.regType2=S1.regType(S1.region==r2);
                        T.regSide1=S1.regSide(S1.region==r1);
                        T.regSide2=S1.regSide(S1.region==r2);
                        TT=addstruct(TT,T);
                        % construct alpha matrix with cosine distances
                        A_subj{ss}(r1,r2,s)=T.cosDist;
                    end
                end
                fprintf('Done sess-%d\tsubj-%d/%d\n',ss,s,length(unique(sn)));
            end
            % save 0.05 for constructing similarity matrix
            A_group{ss}=nanmean(A_subj{ss},3);
            q(ss)=quantile(rsa_vectorizeRDM(A_group{ss}),0.05);
        end
        % make a similarity alpha matrix - using the same threshold
        thres = max(q);
        for i=1:numel(sessN)
            A{i} = 1-exp(-A_group{i}.^2./(2*thres^2));
            %A{i} = 1-A_group{i};
        end
        save(fullfile(transformDir,sprintf('alpha_%s_%sSeq_%sPW_uni',parcelType,seqType,betaChoice)),'A');
        save(fullfile(transformDir,sprintf('transform_cosine_%s_%sSeq_%sPW_uni',parcelType,seqType,betaChoice)),'-struct','TT');
    case 'TRANSFORM_rdms_crossval'
        % use this case to transform G from reg1->reg2
        % save in a subfolder of the clusterDir
        sn=[4:9,11:31];
        parcelType='Brodmann'; %do combination - Brodmann + BG-striatum + thalamus
        betaChoice='multi'; % multi, uni
        sessN=[1:4];
        seqType='all'; %all, trained, untrained
        vararginoptions(varargin,{'sn','parcelType','sessN','betaChoice','seqType'});
        
        TT=[];
        q=zeros(size(sessN));
        
        if ~iscell(parcelType)
            SS = load(fullfile(betaDir,'group',sprintf('betas_partition_%sPW_%s',betaChoice,parcelType)));
            SS=rmfield(SS,{'G_partA','G_partB'});
        else
            for i=1:size(parcelType,2)
                Struct{i} = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType{i},betaChoice,ss)));
            end
            S = sml1_imana_dist('combine_reg','Struct',Struct);
        end
        for ss=sessN
            S = getrow(SS,SS.sessN==ss);
            % pre-allocate alpha matrix
            A{ss}=zeros(length(unique(S.roi)),length(unique(S.roi)),length(unique(sn)));
            for s=1:length(sn)
                S1 = getrow(S,S.sn==sn(s));
                for r1=unique(S.roi)'
                    for r2=unique(S.roi)'
                        % extract RDMs
                        switch seqType
                            case 'trained'
                                dist1 = [S1.partA_train(S1.roi==r1,:);S1.partB_train(S1.roi==r2,:)];
                                dist2 = [S1.partB_train(S1.roi==r1,:);S1.partA_train(S1.roi==r2,:)];
                            case 'untrained'
                                dist1 = [S1.partA_untrain(S1.roi==r1,:);S1.partB_untrain(S1.roi==r1,:)];
                                dist2 = [S1.partA_untrain(S1.roi==r2,:);S1.partB_untrain(S1.roi==r2,:)];
                            case 'all'
                                dist1 = [S1.partA_all(S1.roi==r1,:);S1.partB_all(S1.roi==r1,:)];
                                dist2 = [S1.partA_all(S1.roi==r2,:);S1.partB_all(S1.roi==r2,:)];
                        end               
                        % calculate cosine and KL divergence between G1/2
                        T.cosDist = mean([pdist(dist1,'cosine') pdist(dist2,'cosine')]); 
                    %    T.KL = KLdivergence(T.G1{1},T.G2{1});
                        % calculate transformation matrix G1 -> G2
                     %   [Trans,predG]=calcTransformG(T.G1{1},T.G2{1});
                     %   T.T={round(Trans,3)}; % round
                     %   T.predG={predG};
                        T.sn=sn(s);
                        T.sn_idx=s;
                        T.sessN=ss;
                        T.reg1=r1;
                        T.reg2=r2;
                        T.regType1=S1.regType(S1.roi==r1);
                        T.regType2=S1.regType(S1.roi==r2);
                        T.regSide1=S1.regSide(S1.roi==r1);
                        T.regSide2=S1.regSide(S1.roi==r2);
                        TT=addstruct(TT,T);
                        % construct alpha matrix with cosine distances
                        A_subj{ss}(r1,r2,s)=T.cosDist;
                    end
                end
                fprintf('Done sess-%d\tsubj-%d/%d\n',ss,s,length(unique(sn)));
            end
            % save 0.05 for constructing similarity matrix
            A_group{ss}=nanmean(A_subj{ss},3);
            q(ss)=quantile(rsa_vectorizeRDM(A_group{ss}),0.05);
        end
        % make a similarity alpha matrix - using the same threshold
        thres = max(q);
        for i=1:numel(sessN)
            A_group{i}(1:size(A_group{1},1)+1:end)=zeros(size(A_group{1},1),1);
            A{i} = exp(-A_group{i}.^2./(2*thres^2));
            %A{i} = 1-A_group{i};
        end
        if ~iscell(parcelType)
            save(fullfile(transformDir,sprintf('alpha_%s_%sSeq_%sPW_crossval',parcelType,seqType,betaChoice)),'A');
            save(fullfile(transformDir,sprintf('transform_cosine_%s_%sSeq_%sPW_crossval',parcelType,seqType,betaChoice)),'-struct','TT');
        else
            numParcel = size(parcelType,2);
            for i=1:numParcel
                if i==1
                    nameParcel = sprintf('%s',parcelType{i});
                else
                    nameParcel = sprintf('%s-%s',nameParcel,parcelType{i});
                end
            end
            save(fullfile(transformDir,sprintf('alpha_%s_%sSeq_%sPW',nameParcel,seqType,betaChoice)),'A');
            save(fullfile(transformDir,sprintf('transform_cosine_%s_%sSeq_%sPW',nameParcel,seqType,betaChoice)),'-struct','TT');
        end
    case 'combine_reg'
        vararginoptions(varargin,{'Struct'});
        Gs=[];
        % transform per hemisphere
        rIdx=1; rtF(1)=1; rtS(1)=1;
        for h=1:2
            for i=1:size(Struct,2);
                %  S = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType{i},betaChoice,sessN)));
                S = Struct{i};
                S = getrow(S,S.regSide==h);
                % what to add
                rAdd = rIdx-min(S.region);
                if h==1
                    rtAdd = rtF(i)-min(S.regType);
                else
                    rtAdd = rtS(i)-min(S.regType);
                end
                S.region = S.region+rAdd;
                S.regType = S.regType+rtAdd;
                %S.regType=S.regType+max(Gs.regType)
                Gs = addstruct(Gs,S);
                rIdx=max(S.region)+1; rtF(i+1)=max(S.regType)+1; rtS(i)=min(S.regType);
            end
        end
        varargout={Gs}; 
    case 'TRANSFORM_predict'
        % predict new G based on transformation matrix
        parcelType='Brodmann';
        seqType='all';
        betaChoice='multi';
        sn=[4:9,11:31];
        vararginoptions(varargin,{'parcelType','hemi','seqType','betaChoice'});
        
        S=load(fullfile(transformDir,sprintf('transform_cosine_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)));
        
        % initialize
        PP=[];
        reg = unique(S.reg1)';
        regComb = indicatorMatrix('allpairs',reg);
        D{1} = zeros(length(reg),length(reg),length(unique(S.sn)));
        D{2} = D{1}; D{3}=D{1}; D2=D;
        for s=1:length(unique(S.sn))
            for r=1:size(regComb,1) % go through all region combinations
                rIdx = find(regComb(r,:));
                % extract the right row and transform matrix T from sess1
                S1 = getrow(S,S.sn==sn(s) & S.reg1==rIdx(1) & S.reg2==rIdx(2));
                T = S1.T(S1.sessN==1);
                T = T{1};
                % predict for session 2-4
                for ss=2:4
                    [predG,P.corr,P.cosDist]=predictGfromTransform(S1.G1{ss},T,'G2',S1.G2{ss});
                    P.predG={predG};
                    P.sn=sn(s);
                    P.sessN=ss;
                    P.reg1=rIdx(1);
                    P.reg2=rIdx(2);
                    PP=addstruct(PP,P);
                    % fill deviation matrix D
                    D{ss-1}(rIdx(2),rIdx(1),s)=P.cosDist;
                    D2{ss-1}(rIdx(2),rIdx(1),s)=P.corr;
                end
            end
            fprintf('Done subject %d/%d\n',s,length(sn));
        end
        for i=1:3
            D_group{i}=nanmean(D{i},3);
            q(i)=quantile(rsa_vectorizeRDM(D_group{i}),0.05);
        end
        thres = max(q);
        for i=1:3
            D_group_cos{i} = 1-exp(-D_group{i}.^2./(2*thres^2));
            D_group_corr{i} = 1-nanmean(D2{i},3);
        end
        
        save(fullfile(transformDir,sprintf('predict_cosMatrix_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)),'D_group_cos');
        save(fullfile(transformDir,sprintf('predict_corrMatrix_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)),'D_group_corr');
        save(fullfile(transformDir,sprintf('predict_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)),'-struct','PP');       
    case 'TRANSFORM_feature'
        % assess features in the transformation matrix
         % use this case to transform G from reg1->reg2
        % save in a subfolder of the clusterDir
        sn=[4:9,11:31];
        parcelType='Brodmann'; %do combination - Brodmann + BG-striatum + thalamus
        betaChoice='multi'; % multi, uni
        sessN=[1:4];
        seqType='all'; %all, trained, untrained
        featureName = {'identity','scaling','firstFing','allFing','transitions','chunk','seq'};
        vararginoptions(varargin,{'sn','parcelType','sessN','betaChoice','seqType'});
        
        TT=[];
        
        for ss=sessN
            S = load(fullfile(betaDir,'group',sprintf('stats_%s_%sPW_sess%d',parcelType,betaChoice,ss)));
            for s=1:length(sn)
                S1 = getrow(S,S.SN==sn(s));
                for r1=unique(S.region)'
                    for r2=(r1+1):length(unique(S.region))
                        % extract Gs (square)
                        IPM1=S1.IPM(S1.region==r1,:);
                        IPM2=S1.IPM(S1.region==r2,:);
                        G1=rsa_squareIPM(IPM1);
                        G2=rsa_squareIPM(IPM2);
                        % extract subsets of Gs for trained / untrained
                        if strcmp(seqType,'trained')
                            G1 = G1(1:6,1:6);
                            G2 = G2(1:6,1:6);
                        elseif strcmp(seqType,'untrained')
                            G1 = G1(7:12,7:12);
                            G2 = G2(7:12,7:12);
                        end
                        % get the feature matrices
                        % identity, scaling, first fing, all fing, chunk, seq
                        F = defineFeatures(G1,G2,seqType,Seq,SeqChunks,sn(s));
                        nFeat = size(F,2);
                        for f=1:nFeat
                            % calculate transformation based on features G1 -> G2
                            [predG,T.cor(f,:),T.cosDist(f,:)]=predictGfromTransform(G1,F{f},'G2',G2);
                            T.predG(f,:)={predG};
                            T.featureLabel(f,:)=featureName(f);
                            T.feature(f,:)=f;
                        end
                        % save everything correctly
                        T.G1        = repmat({G1},nFeat,1);
                        T.G2        = repmat({G2},nFeat,1);
                        T.sn        = repmat(sn(s),nFeat,1);
                        T.sn_idx    = repmat(s,nFeat,1);
                        T.sessN     = repmat(ss,nFeat,1);
                        T.reg1      = repmat(r1,nFeat,1);
                        T.reg2      = repmat(r2,nFeat,1);
                        T.regType1  = repmat(S1.regType(S1.region==r1),nFeat,1);
                        T.regType2  = repmat(S1.regType(S1.region==r2),nFeat,1);
                        T.regSide1  = repmat(S1.regSide(S1.region==r1),nFeat,1);
                        T.regSide2  = repmat(S1.regSide(S1.region==r2),nFeat,1);
                        TT = addstruct(TT,T);
                    end
                end
                fprintf('Done sess-%d\tsubj: %d/%d\n',ss,s,length(unique(sn)));
            end
        end
        
        save(fullfile(transformDir,sprintf('transform_Features_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)),'-struct','TT');
    case 'TRANSFORM_plotAlpha'
        % plot alpha
        parcelType='Brodmann';
        hemi=1; % 1, 2 or [1,2]
        seqType='all';
        betaChoice='multi';
        sessNum=4;
        textLabel=regname;
        vararginoptions(varargin,{'parcelType','hemi','seqType','betaChoice'});
        
        if ~iscell(parcelType)
            load(fullfile(transformDir,sprintf('alpha_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)));
            C = load(fullfile(regDir,sprintf('region_%s_centroids_group',parcelType))); % load coordinates
        else
            numParcel = size(parcelType,2);
            for i=1:numParcel
                if i==1
                    nameParcel = sprintf('%s',parcelType{i});
                else
                    nameParcel = sprintf('%s-%s',nameParcel,parcelType{i});
                end
                C1{i} = load(fullfile(regDir,sprintf('region_%s_centroids_group',parcelType{i}))); % load coordinates
            end
            load(fullfile(transformDir,sprintf('alpha_%s_%sSeq_%sPW',nameParcel,seqType,betaChoice)));
            % create a new coordinates structure
            C=sml1_imana_dist('combine_reg','Struct',C1);
            C=rmfield(C,'flatCentr');
        end
        
        if size(hemi,2)==1
            idx=find(C.regSide==hemi)';
            C=getrow(C,C.regSide==hemi);
            for i=1:sessNum
                A{i}=A{i}(idx,idx,:);
            end
        else
            idx=unique(C.region)';
        end
        
        % all pairs of connections
        pIdx=indicatorMatrix('allpairs',[1:length(idx)]);
        figure
        for i=1:sessNum
            A_vec=rsa_vectorizeRDM(A{i});
            subplot(2,sessNum,i)
            imagesc(A{i});
            caxis([0 1]);
            title(sprintf('session-%d',i))
            for k=1:size(pIdx,1)
                p=find(pIdx(k,:));
                subplot(2,sessNum,i+sessNum)
                hold on;
                plot3([C.volCentr(p(1),1),C.volCentr(p(2),1)],...
                    [C.volCentr(p(1),2),C.volCentr(p(2),2)],...
                    [C.volCentr(p(1),3),C.volCentr(p(2),3)],'-o','LineWidth',A_vec(k)*4,'Color',[0 0 0],'MarkerSize',6);
              %  plot([C.flatCentr(p(1),1),C.flatCentr(p(2),1)],...
              %      [C.flatCentr(p(1),2),C.flatCentr(p(2),2)],'-o','LineWidth',A_vec(k)*4,'Color',[0 0 0],'MarkerSize',6);
            end
            for r=1:length(idx)
                if size(hemi,2)==1
                    text(C.volCentr(r,1),C.volCentr(r,2),C.volCentr(r,3),textLabel(r));
                else
                    if r<length(idx)/2+1
                        text(C.volCentr(r,1),C.volCentr(r,2),C.volCentr(r,3),textLabel(r));
                    else
                        text(C.volCentr(r,1),C.volCentr(r,2),C.volCentr(r,3),textLabel(r-length(idx)/2));
                    end
                end
            end
        end
    case 'TRANSFORM_plotPred'
        % plot deviation from prediction for ses2-4 based on T in ses1
        parcelType='Brodmann';
        hemi=1; % 1, 2 or [1,2]
        seqType='all';
        betaChoice='multi';
        metric='corr'; % corr or cos
        textLabel=regname;
        vararginoptions(varargin,{'parcelType','hemi','seqType','betaChoice','metric'});
        
        
        load(fullfile(transformDir,sprintf('predict_%sMatrix_%s_%sSeq_%sPW',metric,parcelType,seqType,betaChoice)));
        C = load(fullfile(regDir,sprintf('region_%s_centroids_group',parcelType))); % load coordinates
        D = eval(sprintf('D_group_%s',metric)); % change the input name to D
        sessNum=size(D,2);
        if size(hemi,2)==1
            idx=find(C.regSide==hemi)';
            C=getrow(C,C.regSide==hemi);
            for i=1:sessNum
                D{i}=D{i}(idx,idx,:);
            end
        end
        
        % all pairs of connections
        pIdx=indicatorMatrix('allpairs',[1:length(idx)]);
        figure
        for i=1:sessNum
            D_vec=rsa_vectorizeRDM(D{i});
            subplot(2,sessNum,i)
            imagesc(D{i});
            caxis([0 0.9]);
            title(sprintf('session-%d',i+1))
            for k=1:size(pIdx,1)
                p=find(pIdx(k,:));
                subplot(2,sessNum,i+sessNum)
                hold on;
                plot3([C.volCentr(p(1),1),C.volCentr(p(2),1)],...
                    [C.volCentr(p(1),2),C.volCentr(p(2),2)],...
                    [C.volCentr(p(1),3),C.volCentr(p(2),3)],'-o','LineWidth',D_vec(k)*3,'Color',[0 0 0],'MarkerSize',6);
            end
            for r=1:length(idx)
                text(C.volCentr(r,1),C.volCentr(r,2),C.volCentr(r,3),textLabel(r));
            end
        end
    case 'TRANSFORM_plotFeatures'
        % plot alpha
        parcelType='Brodmann';
        hemi=1; % 1, 2 or [1,2]
        seqType='all';
        betaChoice='multi';
        sessN=[1:4];
        textLabel=regname;
        vararginoptions(varargin,{'parcelType','hemi','seqType'});
        
        T=load(fullfile(transformDir,sprintf('transform_Features_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)));
        T=getrow(T,ismember(T.regSide1,hemi)&ismember(T.regSide2,hemi));
        % all pairs of connections
        regIdx=unique([unique(T.reg1);unique(T.reg2)]);
        pIdx=indicatorMatrix('allpairs',regIdx');
        
        for k=1:size(pIdx,1)
            idx=find(pIdx(k,:));
            S=getrow(T,T.reg1==idx(1)&T.reg2==idx(2));
            figure
            % subplot per feature
            for i=1:max(S.feature)
                subplot(1,max(S.feature),i);
                plt.bar(S.sessN,S.cosDist,'subset',S.feature==i);
                fLabel=S.featureLabel(S.feature==i);
                if i==1
                    xlabel('Session'); ylabel('Distance from true G')
                else
                    ylabel('');
                end
                title(sprintf('reg %d-%d %s pair',idx(1),idx(2),fLabel{1}));
            end
        end
    case 'TRANSFORM_plotTrans'
        % plot transitions
        parcelType='Brodmann';
        seqType='all';
        betaChoice='multi';
        vararginoptions(varargin,{'parcelType','seqType','betaChoice'});
        
        T = load(fullfile(transformDir,sprintf('transform_sessTrans_cosine_%s_%sSeq_%sPW',parcelType,seqType,betaChoice)));
        
        figure
        plt.bar(T.reg,1-T.cosDist,'split',T.sessTrans,'subset',T.regSide==1 & T.reg~=6);
        xlabel('session transition per region');
        ylabel('similarity');

    case 'job1'
      sml1_imana_dist('PCM_noiseCeilings_seqType');  
      sml1_imana_dist('PCM_noiseCeilings');
    case 'test_dist'
        sml1_imana_dist('SEARCH_dist_map','sn',[26:31],'sessN',4);
        sml1_imana_dist('SURF_wb:map_dist_individ','sessN',1:4); 
        sml1_imana_dist('SURF_wb:map_dist_group','sessN',1:4);
        sml1_imana_dist('SURF_wb:map_dist_individ','sessN',1:4); 
        sml1_imana_dist('SURF_wb:map_dist_sess_diff_individ');
        sml1_imana_dist('SURF_wb:map_sessTrans_group','metric','dist','INname',{'dist_all','dist_trained','dist_untrained','dist_cross'});
        sml1_imana_dist('SURF_wb:cSPM_sessTrans_group');
        sml1_imana_dist('SURF_wb:smooth_group','metric','summary_dist_sessTrans','kernel',1);
    case 'run_job1'
        sml1_imana_dist('HOUSEKEEPING:renameSPM','numError',3);
        sml1_imana_dist('PSC_create','errorNum',3);
        sml1_imana_dist('BETA_get','sn',[5:9,11:27,29:31],'sessN',1:4,'numError',3);
        
    case 'run_job2'
       % sml1_imana_dist('PCM_constructModelFamily_simple','parcelType','Brodmann','reg',1:8);
        sml1_imana_dist('BETA_get','parcelType','Yokoi_clusters','sn',[19:31],'sessN',4);
    
    case 'HOUSEKEEPING:renameSPM' 	% rename SPM directories - to do - new FUNC
        sn = [5:9,11:31];
        sessN = 1:4;
        numError = 1;
        vararginoptions(varargin,{'sn','sessN','numError'});
        fprintf('Renaming SPMs:\n');
        for ss=sessN
            for s=sn
                newGLMDir   = fullfile(glmErrorDir{numError}{ss},subj_name{s});
                cd(newGLMDir);
                load SPM;
                newRawDir   = fullfile(imagingDir,subj_name{s});
                SPM = spmj_move_rawdata(SPM,newRawDir);
                SPM.swd = fullfile(newGLMDir);
                save(fullfile(newGLMDir,'SPM.mat'),'SPM','-v7.3');
                fprintf('Done sess-%d sn-%d\n',ss,s);
            end
        end

    case 'get_ts'
        sessN=1;
        sn = 11:20;
        reg = 1:3;
        vararginoptions(varargin,{'sessN','sn','reg'});
        TT = []; OO = []; YY = [];
        for s=sn
            % load SPM and all regions
            load(fullfile(glmSessDir{sessN},subj_name{s},'SPM.mat'));
            SPM=spmj_move_rawdata(SPM,fullfile(imagingDir, subj_name{s}));
            load(fullfile(regDir,sprintf('%s_Brodmann_regions.mat',subj_name{s})));
            % extract data and onsets
           data = region_getdata(SPM.xY.VY,R(reg)); % here get the data
           O = spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
           O.sn = ones(size(O.block))*sn;
            % maybe consider filtering here to get y_adj
            for r=1:size(data,2)
           %for r=reg
                Y.y_raw     = {data{r}};
                Y.y_filt    = {spm_filter(SPM.xX.K,SPM.xX.W*Y.y_raw{:})};
                Y.B         = {SPM.xX.pKX*Y.y_filt{:}};
                Y.y_res     = {spm_sp('r',SPM.xX.xKXs,Y.y_filt{:})};
                Y.y_hat     = {SPM.xX.xKXs.X*Y.B{:}};
                Y.y_adj     = {Y.y_hat{:} + Y.y_res{:}};
                Y.sn        = s;
                Y.reg       = r;
                T.xyz       = {R{r}.data};
                T.flatcoord = {R{r}.flatcoord};
                T.depth     = {R{r}.depth};
                T.sn        = s;
                T.reg       = r;
                TT = addstruct(TT,T); %YY = addstruct(YY,Y);
            end
           % OO = addstruct(OO,O); 
            fprintf('Done %s\n',subj_name{s});
        end
        fprintf('Done all subjects.\n\n\n');
        % save structure
        keyboard;
        save(fullfile(distPscDir,'connect_timeseries'),'-struct','TT');
        
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
function M = pcm_defineSequenceModels_random(Seq,sn)
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
    
    M{5}(:,:)=zeros(5); % first finger
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
function M = pcm_defineSequenceModels_fixed(Seq,Chunks,sn)
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
    % Model4: Chunk model
    nChunks=max(max(Chunks));
    A=zeros(nSeq,nChunks);
    for i=1:nSeq
        idxChunk=Chunks(i,:);
        A(i,idxChunk)=1;
    end
    M{4} = A;
    % ---------------------------------------
    % Model5: Model with trained & untrained labels - 
    M{5}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{5}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model6: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{6}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model7: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{7}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
end
function M = pcm_defineSequenceModels_fixed_noChunks(Seq,sn)
% written November 29th 2019
%  No fixed / run component added
% specific sequence + overall component modelled together
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
    % ---------------------------------------
    % Model3: Model with trained & untrained labels - 
    M{3}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{3}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model4: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{4}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model5: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
end
function M = pcm_defineSequenceModels_fixed_noChunks_natStats(Seq,sn,NatStat)
% written November 29th 2019
%  No fixed / run component added
% incorporating natural statistics
% specific sequence + overall component modelled together
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
    M{1} = M{1}*NatStat*M{1}';
    % ---------------------------------------
    % Model2: All fingers
    M{2}(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{2}(u,j) = length(placenumb);
        end
    end
    M{2} = M{2}*NatStat*M{2}';
    % ---------------------------------------
    % Model3: Model with trained & untrained labels - 
    M{3}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{3}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model4: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{4}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model5: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
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
function M = pcm_defineSequenceModels_seqType(Seq,Chunks,sn,seqType)
%  Model for one seqType - trained or untrained sequences
% Run separately for trained and untrained
% seqType is index whether trained (1) or untrained (2)

    % --------------------------------------

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
    % Model1: First finger model
    M{1}(:,:)=zeros(6,5); % first finger
    for u = 1:size(SeqSet,1)    % 6 seq
        firstfing = SeqSet(:,1);
        M{1}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model2: All fingers
    M{2}(:,:)=zeros(6,5); % all fingers
    for u = 1:size(SeqSet,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(SeqSet(u,:)==j);
            M{2}(u,j) = length(placenumb);
        end
    end
    %---------------------------------------
    % Model3: Finger transitions
    Trans = [nchoosek([1:5],2); fliplr(nchoosek([1:5],2))];
    M{3}(:,:)=zeros(6,size(Trans,1)); % all transitions (no repetition)
    for k = 1:size(SeqSet,1)    % set of 12 seq
        for p = 1:size(Trans,1) % all possible transitions
            if any(ismember(SeqSet(k,1:8),Trans(p,1)))
                ind = find(ismember(SeqSet(k,1:8),Trans(p,1)));    % find matching digits to first finger in doublet
                if any(ismember(SeqSet(k,ind+1),Trans(p,2)))       % compare the second finger in doublet
                    M{3}(k,p) = M{3}(k,p) + sum(ismember(SeqSet(k,ind+1),Trans(p,2)));
                end
            end
        end  
    end
    tmp=M{3};
    % remove columns with 0s only
    tmp(:,~any(tmp,1))=[];
    M{3}=tmp;
    %----------------------------------------
    % Model 4: Chunks
    ChunkSet = Chunks([((seqType-1)*6)+1:((seqType-1)*6)+6],:);
    if seqType==2
        ChunkSet = ChunkSet-7; % 7 unique chunks
    end
    A=zeros(6,7);
    for i=1:6
        idxChunk=ChunkSet(i,:);
        A(i,idxChunk)=1;
    end
    M{4} = A;
    % --------------------------------------
    % Model5: Model for each specific sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}(:,1:6)    = A;       % Trained sequence patterns
end
function M = pcm_defineSequenceModels_seqType_simple(Seq,sn,seqType)
% simpler version with no chunks and transitions
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
    % Model1: First finger model
    M{1}(:,:)=zeros(6,5); % first finger
    for u = 1:size(SeqSet,1)    % 6 seq
        firstfing = SeqSet(:,1);
        M{1}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model2: All fingers
    M{2}(:,:)=zeros(6,5); % all fingers
    for u = 1:size(SeqSet,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(SeqSet(u,:)==j);
            M{2}(u,j) = length(placenumb);
        end
    end
    % --------------------------------------
    % Model3: Model for each specific sequence
    A=zeros(3);
    for i=1:6
        A(i,i)=1;
    end;
    M{3}(:,1:6)    = A;       % Trained sequence patterns
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
function M = pcm_defineSequenceFingerModels_seqType(Seq,sn,seqType)
%  No fixed / run component added
% specific sequence + overall component modelled together
    % --------------------------------------
     if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew(1:6,:)=Seq(7:12,:);
        SeqNew(7:12,:)=Seq(1:6,:);
        Seq=SeqNew;
    end
    if seqType==1
        SeqSet=Seq(1:6,:);
    elseif seqType==2
        SeqSet=Seq(7:12,:);
    end
    % Model1: Model for each specific SPECIFIC sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{1}(:,1:6)    = A;       % Trained sequence patterns
    % ---------------------------------------
    % Model2: First finger model
    M{2}(:,:)=zeros(6,5); % first finger
    for u = 1:size(SeqSet,1)    % 6 seq
        firstfing = SeqSet(:,1);
        M{2}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model3: All fingers
    M{3}(:,:)=zeros(6,5); % all fingers
    for u = 1:size(SeqSet,1)
        for j = 1:5                 % 5 fingers
            placenumb = find(SeqSet(u,:)==j);
            M{3}(u,j) = length(placenumb);
        end
    end
end
function F = defineFeatures(G1,G2,seqType,Seq,Chunks,sn)
%  defines feature matrices
% currently used in estimation of G1 -> G2 transformation
    numSeq=size(G1,1);
    % ---------------------------------------
    % F1 - identity
    F{1} = eye(numSeq);
    % ---------------------------------------
    % F2 - scaling
    % optimise weight
    w = fminsearch(@(w) abs(sum(sum(G2-(pinv(F{1}*w)*G1*pinv((F{1}*w)'))))),1);
    F{2} = F{1}*w;
    % ---------------------------------------
    % sequence-specific things
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end

    switch seqType
        case 'trained'
            SeqSet = Seq([1:6],:);
            ChunkSet = Chunks([1:6],:);
        case 'untrained'
            SeqSet = Seq([7:12],:);
            ChunkSet = Chunks([7:12],:);
            ChunkSet = ChunkSet-7; % 7 unique chunks
        case 'all'
            SeqSet = Seq;
            ChunkSet = Chunks;
    end
    % ---------------------------------------
    % F3: First finger model
    F{3}(:,:)=zeros(numSeq,5); % first finger
    for u = 1:numSeq   % 12 seq
        firstfing = SeqSet(:,1);
        F{3}(u,firstfing(u)) = 1;
    end
    F{3}=F{3}*F{3}';
    w = fminsearch(@(w) abs(sum(sum(G2-(pinv(F{3}*w)*G1*pinv((F{3}*w)'))))),1);
    F{3}=F{3}*w;
    % ---------------------------------------
    % F4: All fingers
    F{4}(:,:)=zeros(numSeq,5); % all fingers
    for u = 1:numSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(SeqSet(u,:)==j);
            F{4}(u,j) = length(placenumb);
        end
    end
    F{4}=F{4}*F{4}';
    w = fminsearch(@(w) abs(sum(sum(G2-(pinv(F{4}*w)*G1*pinv((F{4}*w)'))))),1);
    F{4}=F{4}*w;
    %---------------------------------------
    % F5: Finger transitions
    Trans = [nchoosek([1:5],2); fliplr(nchoosek([1:5],2))];
    F{5}=zeros(numSeq,size(Trans,1)); % all transitions (no repetition)
    for k = 1:size(SeqSet,1)    % set of 12 seq
        for p = 1:size(Trans,1) % all possible transitions
            if any(ismember(SeqSet(k,1:8),Trans(p,1)))
                ind = find(ismember(SeqSet(k,1:8),Trans(p,1)));    % find matching digits to first finger in doublet
                if any(ismember(SeqSet(k,ind+1),Trans(p,2)))       % compare the second finger in doublet
                    F{5}(k,p) = F{5}(k,p) + sum(ismember(SeqSet(k,ind+1),Trans(p,2)));
                end
            end
        end
    end
    tmp=F{5};
    % remove columns with 0s only
    tmp(:,~any(tmp,1))=[];
    F{5}=tmp*tmp';
    w = fminsearch(@(w) abs(sum(sum(G2-(pinv(F{5}*w)*G1*pinv((F{5}*w)'))))),1);
    F{5}=F{5}*w;
    %----------------------------------------
    % F6: Chunks
    A=zeros(numSeq,max(max(ChunkSet)));
    for i=1:numSeq
        idxChunk=ChunkSet(i,:);
        A(i,idxChunk)=1;
    end
    F{6} = A*A';
    w = fminsearch(@(w) abs(sum(sum(G2-(pinv(F{6}*w)*G1*pinv((F{6}*w)'))))),1);
    F{6}=F{6}*w;
    % --------------------------------------
    % F7: Sequence + seqType
    A=zeros(numSeq,numSeq);
    for i=1:numSeq
        A(i,i)=1;
    end;
    if strcmp(seqType,'all')
        A(:,i+1)=[ones(6,1);zeros(6,1)];
        A(:,i+2)=[zeros(6,1);ones(6,1)];
    else
        A(:,i+1)=[ones(6,1)];
    end
    F{7} = A*A';
    w = fminsearch(@(w) abs(sum(sum(G2-(pinv(F{7}*w)*G1*pinv((F{7}*w)'))))),1);
    F{7}=F{7}*w;
end

function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Crossvalidated model comparision:

    [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Group fit with %s algorithm done.\n',algorithm);
  %  [T1,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize');
  %  fprintf('Group fit with minimize algorithm done.\n');
  %  [T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR');
  %  fprintf('Group fit with NR algorithm done.\n');
    
  %  [T3,theta_hat3,G_pred3,theta03] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR','theta0',theta_hat);
  %  [T4,theta_hat4,G_pred4,theta04] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize','theta0',theta_hat2);
    
    [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);

    T.cross_likelihood = Tcross.likelihood;
    T.bayesEst = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
    T.thetaCr = thetaCr;
   % T.modelNum = [1:length(T.cross_likelihood)];
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