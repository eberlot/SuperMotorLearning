function varargout=sml1_imana_repsup(what,varargin)

% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
%baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
behavDir        =[baseDir '/behavioral_data/data'];  
anaDir          =[baseDir '/behavioral_data/analyze'];          
%imagingDir      =[baseDir '/imaging_data'];                      
imagingDir      ='/Volumes/Eva_data/SuperMotorLearning/imaging_data';
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia_new'];
pcmDir          =[baseDir '/pcm_stats'];
repSupDir       =[baseDir '/repsup_stats'];
wbDir           =[baseDir '/surfaceWB'];
codeDir         ='/Users/eberlot/Documents/MATLAB/projects/SuperMotorLearning';
% update glmDir when adding new glms
glmFoSExDir     ={[baseDir '/glmFoSEx/glmFoSEx1'],[baseDir '/glmFoSEx/glmFoSEx2'],[baseDir '/glmFoSEx/glmFoSEx3'],[baseDir '/glmFoSEx/glmFoSEx4']};    
glmFoSExErrorDir ={[baseDir '/glmFoSExError/glmFoSExError1'],[baseDir '/glmFoSExError/glmFoSExError2'],[baseDir '/glmFoSExError/glmFoSExError3'],[baseDir '/glmFoSExError/glmFoSExError4']};    
%glmFoSExDir     ={[baseDir '/glmFoSEx_old/glmFoSEx1'],[baseDir '/glmFoSEx_old/glmFoSEx2'],[baseDir '/glmFoSEx_old/glmFoSEx3'],[baseDir '/glmFoSEx_old/glmFoSEx4']};    
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
hem        = {'lh','rh'};      
hemI       = {'L','R'};      % left & right hemi folder names/prefixes
hemName    = {'LeftHem','RightHem'};
hemname    = {'CortexLeft','CortexRight'};
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
ms=6;
stySeq=style.custom({'red','blue'},'markersize',ms);
stySeqType_exe=style.custom({black,lightgray},'markersize',ms);
styTrained_exe=style.custom({red,lightred},'markersize',ms);
styUntrained_exe=style.custom({blue,lightblue},'markersize',ms);
stySess=style.custom({black,lightgray,silver,black},'markersize',ms);
styRS = style.custom({gray,silver,lightgray});
styRSbeh = style.custom({red,lightred,blue,lightblue},'markersize',ms);
sty_exe1 = style.custom({red,blue},'markersize',ms);
sty_exe2 = style.custom({lightred,lightblue},'markersize',ms);
stySeqShade = style.custom({'red','blue'},'errorbars','shade','markertype','.');
stySeqReg = style.custom({'blue','orange','magenta'},'errorbars','shade','markertype','.');
%% --------------------------------------- Analysis Cases -----------------------------------------------
switch(what)


    case 'BEH_scanner_MT' %----------------------- QUANTIFY BEHAVIOUR IN SCANNER -----------------------
        sn=[5:9,11:31];
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
        save(fullfile(repSupDir,'beh_scanner_new.mat'),'-struct','B_fun');   
    case 'PLOT_behScanner'
        
        B = load(fullfile(repSupDir,'beh_scanner_new.mat'));
        
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
        figure % here violin plot
        N1 = getrow(N,N.sn~=29);
        colblue=[49,130,189]/255;
        collightblue=[158,202,225]/255;
        colred=[222,45,38]/255;
        collightred=[251,177,168]/255;
 
        for ss=1:4
            subplot(1,4,ss)
            distributionPlot(N1.MT(N1.FoSEx==1 & N1.sessN==ss & N1.seqType==1),'histOri','left','color',colred,'widthDiv',[2 1],'xValues',1,'showMM',3,'globalNorm',0);
            distributionPlot(N1.MT(N1.FoSEx==2 & N1.sessN==ss & N1.seqType==1),'histOri','right','color',collightred,'widthDiv',[2 2],'xValues',1,'showMM',3,'globalNorm',0);
            distributionPlot(N1.MT(N1.FoSEx==1 & N1.sessN==ss & N1.seqType==2),'histOri','left','color',colblue,'widthDiv',[2 1],'xValues',2,'showMM',3,'globalNorm',0);
            distributionPlot(N1.MT(N1.FoSEx==2 & N1.sessN==ss & N1.seqType==2),'histOri','right','color',collightblue,'widthDiv',[2 2],'xValues',2,'showMM',3,'globalNorm',0);
            hold on;
            title(sprintf('session-%d',ss));
        end
        plt.match('y');
    case 'PSC_create_RepSup'
        % calculate psc for trained and untrained sequences (1st/2nd) - based on betas
        sn = [5:9,11:31];
        sessN = 1:4;
        vararginoptions(varargin,{'sn','sessN'});
        name={'AllSeq_1st','AllSeq_2nd','TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
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
                end
                for con=1:length(name)    % 4 contrasts
                    P{numB+1}=sprintf('con_%s.nii',name{con});
                    outname=sprintf('psc_sess%d_%s.nii',ss,name{con}); % ,subj_name{s}
                    formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);
                    spm_imcalc_ui(P,outname,formula,...
                        {0,[],spm_type(16),[]});        % Calculate percent signal change
                end
                fprintf('Subject %d sess %d: %3.3f\n',s,ss,h);
            end
        end
    case 'PSC_create_RepSup_noErr'
        % create psc volume niftis for glm with no errors
        sn = [5:9,11:31];
        sessN = 4;
        vararginoptions(varargin,{'sn','sessN'});
        name={'AllSeq_1st','AllSeq_2nd','TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
        for ss=sessN
            for s=sn
                cd(fullfile(glmFoSExErrorDir{ss}, subj_name{s}));
                load SPM;
                T=load('SPM_info.mat');
                X=(SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
                h=median(max(X));               % Height of response;
                P={};
                numB=length(SPM.xX.iB);         % Partitions - runs
                for p=SPM.xX.iB
                    P{end+1}=sprintf('beta_%4.4d.nii',p);       % get the intercepts and use them to calculate the baseline (mean images)
                end
                for con=1:length(name)    % 4 contrasts
                    P{numB+1}=sprintf('con_%s.nii',name{con});
                    outname=sprintf('psc_sess%d_%s.nii',ss,name{con}); % ,subj_name{s}
                    formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);
                    spm_imcalc_ui(P,outname,formula,...
                        {0,[],spm_type(16),[]});        % Calculate percent signal change
                end
                fprintf('Subject %d sess %d: %3.3f\n',s,ss,h);
            end
        end
    case 'PSC_surface_RepSup' %-------------------- CALCULATE PSC, PROJECT ON SURFACE - OLD CARET -------------------
        % calculate psc for trained and untrained sequences (1st/2nd) - based on betas    
        vararginoptions(varargin,{'sn','sessN'});
        name={'TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
        for seq=1:12
            name=[name sprintf('Seq%d-exe1',seq),sprintf('Seq%d-exe2',seq)];
        end
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
     for seq=1:12
         name=[name sprintf('Seq%d-exe1',seq),sprintf('Seq%d-exe2',seq)];
     end
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
    
    case 'SURF_wb:map_individ'  %------------------------ PROJECT ON SURFACE - WORKBENCH -----------------------
        % projects individual percent signal change volume files to WorkBench surface
        sessN=4;
        sn=[5:9,11:31];
      %  name={'AllSeq_1st','AllSeq_2nd','TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'}; % for psc
         name = {'exe1_dist','exe2_dist','exe1_dist_trained','exe2_dist_trained','exe1_dist_untrained',...
        'exe2_dist_untrained','exe1_dist_cross','exe2_dist_cross'}; % for dist
       % colName={'psc_all_exe1','psc_all_exe2','psc_trained_exe1','psc_trained_exe2','psc_untrained_exe1','psc_untrained_exe2'}; % for psc
        colName={'dist_all_exe1','dist_all_exe2','dist_trained_exe1','dist_trained_exe2','dist_untrained_exe1',...
            'dist_untrained_exe2','dist_cross_exe1','dist_cross_exe2'}; % for dist
        metric = 'dist'; % psc or dist
        vararginoptions(varargin,{'sn','sessN','name','colName','metric'});

        for ss=sessN
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    fprintf('\nSess: %d, Subj: %d, Hem: %d\t',ss,s,h);
                    white   = fullfile(subjDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},hemI{h}));
                    pial    = fullfile(subjDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},hemI{h}));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    for f=1:length(name)
                        if strcmp(metric,'psc')
                            images{f}=fullfile(glmFoSExDir{ss},subj_name{s},sprintf('psc_sess%d_%s.nii',ss,name{f}));
                        else
                            images{f}=fullfile(glmFoSExDir{ss},subj_name{s},sprintf('%s_sess%d_%s.nii',subj_name{s},ss,name{f}));
                        end
                        column_name{f} = fullfile(sprintf('Sess%d_%s.nii',ss,colName{f}));
                    end;
                    outfile         = fullfile(subjDir,sprintf('%s.%s.%s_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss));
                    G               = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names',column_name, ...
                        'anatomicalStruct',hemname{h});
                    if strcmp(metric,'dist')
                        for c=1:size(G.cdata,2)
                            G.cdata(G.cdata(:,c)>5,c)=0;
                            G.cdata(G.cdata(:,c)<-5,c)=0;
                        end
                    end
                    save(G, outfile);
                end; % hemi
            end; % sn
        end; % session
    case 'SURF_wb:map_RS_absolute'
        % projects repetition suppresion (1st-2nd exe)
        sessN=1:4;
        sn=[5:9,11:31];
        metric = 'psc'; % psc or dist
        vararginoptions(varargin,{'sn','sessN','name','metric'});
        
        colName={sprintf('%s_all_RS_abs',metric),sprintf('%s_trained_RS_abs',metric),sprintf('%s_untrained_RS_abs',metric)};

        for ss=sessN
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    G = gifti(fullfile(subjDir,sprintf('%s.%s.%s_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss)));
                    for f=1:length(colName)
                        %data(:,f) = G.cdata(:,(f-1)*2+1)-G.cdata(:,(f-1)*2+2); % columns 3-4 or 5-6 for trained  untrained
                        data(:,f) = G.cdata(:,(f-1)*2+2)-G.cdata(:,(f-1)*2+1); % columns 3-4 or 5-6 for trained  untrained
                        column_name{f} = sprintf('Sess%d_%s',ss,colName{f});
                    end;
                    outfile         = fullfile(subjDir,sprintf('%s.%s.%s_RepSup_absolute.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss));
                    G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
                    save(G, outfile);
                    fprintf('Sess: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session
    case 'SURF_wb:map_RS_relative'
        % projects repetition suppresion ((1st-2nd)/2nd exe)
        sessN=1:4;
        sn=[5:9,11:31];
        metric = 'psc'; % psc or dist
        smooth = 1; % here already using in calculation - otherwise some values can be very erratic
        vararginoptions(varargin,{'sn','sessN','name','colName','metric'});

        colName={sprintf('%s_all_RS_rel',metric),sprintf('%s_trained_RS_rel',metric),sprintf('%s_untrained_RS_rel',metric)};
        for ss=sessN
            for s=sn
                subjDir = fullfile(wbDir,subj_name{s});
                for h=1:2
                    if smooth
                        G = gifti(fullfile(subjDir,sprintf('s%s.%s.%s_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss)));
                    else
                        G = gifti(fullfile(subjDir,sprintf('%s.%s.%s_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss)));
                    end
                    for f=1:length(colName)
                        data(:,f) = (G.cdata(:,(f-1)*2+2)./G.cdata(:,(f-1)*2+1));
                        %data(:,f) = (G.cdata(:,(f-1)*2+3)-G.cdata(:,(f-1)*2+4))./-G.cdata(:,(f-1)*2+3);
                        column_name{f} = sprintf('Sess%d_%s',ss,colName{f});
                    end;
                    if smooth
                        outfile = fullfile(subjDir,sprintf('s%s.%s.%s_RepSup_relative.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss));
                    else
                        outfile = fullfile(subjDir,sprintf('%s.%s.%s_RepSup_relative.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss));
                    end
                    G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',column_name);
                    save(G, outfile);
                    fprintf('Sess: %d, Subj: %d, Hem: %d\n',ss,s,h);
                end; % hemi
            end; % sn
        end; % session
    case 'SURF_wb:smooth_individ'
        metric = 'psc_RepSup_absolute'; % dist_RepSup, psc_RepSup, dist_RepSup_absolute, psc_RepSup_absolute
        kernel = 2;
        sessN = 1:4;
        sn = [5:9,11:31];
        atlas = 'FoSEx_FS_LR_164'; % 164 or 42
        vararginoptions(varargin,{'metric','sessN','kernel','sn'});
        
        for ss=sessN
            for s=sn
                for h=1:2
                    subjDir = fullfile(wbDir,subj_name{s});
                    inFile = fullfile(subjDir,sprintf('%s.%s.%s.sess-%d.func.gii',subj_name{s},hemI{h},metric,ss));
                    surfFile = fullfile(wbDir,atlas,sprintf('fs_LR.164k.%s.inflated.surf.gii',hemI{h}));
                    surf_smooth(inFile,'surf',surfFile,'kernel',kernel);
                end
            end
            fprintf('Done all subjects sess-%d.\n',ss);
        end
    case 'SURF_wb:map_group'
        sessN=1:4;
        sn=[5:9,11:31];
        %metric = 'psc_RepSup'; % psc_RepSup / psc_RepSup_absolute /
        %psc_RepSup_relative / dist_RepSup / dist_RepSup_absolute /
        %dist_RepSup_relative / RepSup_fatigueModel
        metric = 'psc_RepSup_absolute';
        %metric = 'RepSup_fatigueModel_ico-642';
        replaceNaN = 1;
        smooth = 1;
        %INname = {'psc_all_exe1','psc_all_exe2','psc_trained_exe1','psc_trained_exe2','psc_untrained_exe1','psc_untrained_exe2'};  
        INname = {'psc_all_RS_abs','psc_trained_RS_abs','psc_untrained_RS_abs'};
        %INname = {'psc_all_RS_rel','psc_trained_RS_rel','psc_untrained_RS_rel'};
        %INname ={'dist_all_exe1','dist_all_exe2','dist_trained_exe1','dist_trained_exe2','dist_untrained_exe1','dist_untrained_exe2','dist_cross_exe1','dist_cross_exe2'};
        %INname = {'dist_all_RS_abs','dist_trained_RS_abs','dist_untrained_RS_abs'};
        %INname = {'dist_all_RS_rel','dist_trained_RS_rel','dist_untrained_RS_rel'};
        %INname = {'all_ratio_deviation_ico-642','trained_ratio_deviation_ico-642','untrained_ratio_deviation_ico-642'};
        vararginoptions(varargin,{'sessN','sn','metric','INname','replaceNaN','smooth'});
        
        for ss=sessN
            % Loop over hemispheres.
            for h = 1:2
                % Go to the directory where the group surface atlas resides
                surfaceGroupDir = fullfile(wbDir,'FoSEx_FS_LR_164');
                cd(surfaceGroupDir);
                % Loop over each input metric file in 'INname' and make a group metric file
                for j=1:length(INname)
                    for i = 1:length(sn);
                        if smooth
                            infilenames{i} = fullfile(wbDir,subj_name{sn(i)},sprintf('s%s.%s.%s.sess-%d.func.gii',subj_name{sn(i)},hemI{h},metric,ss));
                        else
                            infilenames{i} = fullfile(wbDir,subj_name{sn(i)},sprintf('%s.%s.%s.sess-%d.func.gii',subj_name{sn(i)},hemI{h},metric,ss));
                        end
                    end;

                    if smooth
                        outfilenames    = fullfile(surfaceGroupDir,sprintf('s%s.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                        summaryname     = fullfile(surfaceGroupDir,sprintf('s%s.group.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                    else
                        outfilenames    = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                        summaryname     = fullfile(surfaceGroupDir,sprintf('%s.group.%s.sess-%d.func.gii',hemI{h},INname{j},ss));
                    end
                    surf_groupGiftis(infilenames,'outfilenames',{outfilenames},'inputcol',j,'groupsummary',summaryname,'replaceNaNs',replaceNaN);
                end
                fprintf('Done: %s - sess%d\n',hemI{h},ss);
            end;
        end;
    case 'SURF_wb:cSPM_group'
        % calculate mean statistics (t-test)
        metric = 'psc_abs'; 
        smooth = 1;
        sessN = 1:4;
        vararginoptions(varargin,{'metric','sessN','atlas','smooth'});
        switch metric
            case 'psc'
                INname = {'psc_all_exe1','psc_all_exe2','psc_trained_exe1','psc_trained_exe2','psc_untrained_exe1','psc_untrained_exe2'}; 
            case 'psc_abs' 
                INname = {'psc_all_RS_abs','psc_trained_RS_abs', 'psc_untrained_RS_abs'};
            case 'psc_rel'
                INname = {'psc_all_RS_rel','psc_trained_RS_rel', 'psc_untrained_RS_rel'};
            case 'dist'
                INname = {'dist_all_exe1','dist_all_exe2','dist_trained_exe1','dist_trained_exe2','dist_untrained_exe1',...
                    'dist_untrained_exe2','dist_cross_exe1','dist_cross_exe2'}; 
            case 'dist_abs'
                INname = {'dist_all_RS_abs','dist_trained_RS_abs', 'dist_untrained_RS_abs'};
            case 'dist_rel'
                INname = {'dist_all_RS_rel','dist_trained_RS_rel', 'dist_untrained_RS_rel'};
            case 'fatigueModel'
                INname = {'trained_ratio_deviation','untrained_ratio_deviation'};
            case 'fatigueModel_ico-642'
                INname = {'all_ratio_deviation_ico-642','trained_ratio_deviation_ico-642','untrained_ratio_deviation_ico-642'};
        end
        for h=1:2
            indx=1;
            fprintf('\n%s: session:',hemI{h});
            for ss=sessN
                surfaceGroupDir = fullfile(wbDir,'FoSEx_FS_LR_164');
                for n=1:length(INname)
                    if smooth
                        inFile = fullfile(surfaceGroupDir,sprintf('s%s.%s.sess-%d.func.gii',hemI{h},INname{n},ss));
                    else
                        inFile = fullfile(surfaceGroupDir,sprintf('%s.%s.sess-%d.func.gii',hemI{h},INname{n},ss));
                    end
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
            if smooth
                summaryName = fullfile(surfaceGroupDir,sprintf('s%s.summary_%s_RepSup.func.gii',hemI{h},metric));
            else
                summaryName = fullfile(surfaceGroupDir,sprintf('%s.summary_%s_RepSup.func.gii',hemI{h},metric));
            end
            save(G2,summaryName);
        end
    case 'SURF_wb:group_smooth' % not to be used 
        metric = 'summary_psc_RepSup'; % summary_psc_RepSup, summary_psc_abs_RepSup
        atlas = 'FoSEx_FS_LR_164'; % 164 or 42
        kernel = 1;
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
    case 'SURF_wb:calculate_RepSup_model'
        % here calculate the repetition suppression model relationship
        % to what extent does activation vs. representation collapse
        % scaling, fatigue, sharpening
        sn = [5:9,11:31];
        sessN = 4;
        smooth = 1;
        vararginoptions(varargin,{'sn','sessN','smooth'});
        
        inCol = [3,5]; % trained / untrained
        colName = {'trained_ratio_deviation','untrained_ratio_deviation'};
        for ss=sessN
            for h=1:2
                for s=sn
                    subjDir = fullfile(wbDir,subj_name{s});
                    M = gifti(fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('s%s.summary_dist_RepSup.func.gii',hemI{h})));
                    mask = M.cdata(:,1)>.01;
                    if smooth
                        pscFile = fullfile(subjDir,sprintf('s%s.%s.psc_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                        distFile = fullfile(subjDir,sprintf('s%s.%s.dist_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    else
                        pscFile = fullfile(subjDir,sprintf('%s.%s.psc_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                        distFile = fullfile(subjDir,sprintf('%s.%s.dist_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    end
                    P = gifti(pscFile);
                    D = gifti(distFile);
                    data = zeros(size(P.cdata,1),size(inCol,2));
                    for i=1:size(inCol,2)
                        scale = (P.cdata(:,inCol(i))-P.cdata(:,inCol(i)+1))./P.cdata(:,inCol(i));
                        expDist = D.cdata(:,inCol(i)).*scale;
                        data(:,i) = D.cdata(:,inCol(i)+1)-expDist;
                        %p_ratio = P.cdata(:,inCol(i)+1)./P.cdata(:,inCol(i));
                        %d_ratio = D.cdata(:,inCol(i)+1)./D.cdata(:,inCol(i));
                        %data(:,i) = d_ratio-p_ratio;
                        %data(:,i) = data(:,i).*mask;
                    end
                    G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',colName);
                    if smooth
                        outfile = fullfile(subjDir,sprintf('s%s.%s.RepSup_fatigueModel.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    else
                        outfile = fullfile(subjDir,sprintf('%s.%s.RepSup_fatigueModel.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    end
                    save(G,outfile);
                    fprintf('sess-%d %s %s\n',ss,subj_name{s},hemName{h});
                end
            end
        end
    case 'SURF_wb:calculate_RepSup_model_ico'
        % calculate the repetition suppression model per icosahedron 
        % not per surface node as above
        sn = [5:9,11:31];
        sessN = 4;
        smooth = 1;
        nTessel = 642; % 164, 362, 642
        vararginoptions(varargin,{'sn','sessN','smooth','nTessel'});
        
        inCol = [1,3,5]; % all / trained / untrained
        thres = .3; % half of the icosahedron needs to have sufficiently large distances
        colName = {sprintf('all_ratio_deviation_ico-%d',nTessel),sprintf('trained_ratio_deviation_ico-%d',nTessel),sprintf('untrained_ratio_deviation_ico-%d',nTessel)};
        for ss=sessN
            for h=1:2
                Ico = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{h})));
                M = gifti(fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('s%s.summary_dist_RepSup.func.gii',hemI{h})));
                mask = M.cdata(:,1)>.01;
                for s=sn
                    subjDir = fullfile(wbDir,subj_name{s});

                    if smooth
                        pscFile = fullfile(subjDir,sprintf('s%s.%s.psc_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                        distFile = fullfile(subjDir,sprintf('s%s.%s.dist_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    else
                        pscFile = fullfile(subjDir,sprintf('%s.%s.psc_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                        distFile = fullfile(subjDir,sprintf('%s.%s.dist_RepSup.sess-%d.func.gii',subj_name{s},hemI{h},ss));
                    end
                    P = gifti(pscFile);
                    D = gifti(distFile);
                    data = zeros(size(P.cdata,1),size(inCol,2));
                    for i=1:max(Ico.cdata)
                        if sum(mask(Ico.cdata==i,:))/length(mask(Ico.cdata==i,:))>thres          
                            for j=1:size(inCol,2)
                                p_ratio = mean(P.cdata(Ico.cdata==i,inCol(j)+1))/mean(P.cdata(Ico.cdata==i,inCol(j)));
                                %d_ratio = mean(D.cdata(Ico.cdata==i,inCol(j)+1))/mean(D.cdata(Ico.cdata==i,inCol(j)));
                                %data(Ico.cdata==i,j) = d_ratio-p_ratio;
                                d_expect = p_ratio*mean(D.cdata(Ico.cdata==i,inCol(j)));
                                data(Ico.cdata==i,j) = mean(D.cdata(Ico.cdata==i,inCol(j)+1))-d_expect; % deviation from expected scaling
                            end
                        end
                    end
                    G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',colName);
                    if smooth
                        outfile = fullfile(subjDir,sprintf('s%s.%s.RepSup_fatigueModel_ico-%d.sess-%d.func.gii',subj_name{s},hemI{h},nTessel,ss));
                    else
                        outfile = fullfile(subjDir,sprintf('%s.%s.RepSup_fatigueModel_ico-%d.sess-%d.func.gii',subj_name{s},hemI{h},nTessel,ss));
                    end
                    save(G,outfile);
                    fprintf('sess-%d %s %s\n',ss,subj_name{s},hemName{h});
                end
            end
        end
    case 'SURF_wb:calculate_relative_ico'
        metric='dist'; % dist or psc
        % calculate per icosahedron
        sessN = 1:4;
        smooth = 1;
        nTessel = 642; % 164, 362, 642
        vararginoptions(varargin,{'sn','sessN','smooth','nTessel','metric'});
        
        thres = .5; % half of the icosahedron needs to have sufficiently large distances
        seqType = {'all','trained','untrained'};
        colName = {sprintf('trained_ratio_deviation_ico-%d',nTessel),sprintf('untrained_ratio_deviation_ico-%d',nTessel)};
        for h=1:2
            Ico = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{h})));
            M = gifti(fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('s%s.summary_dist_RepSup.func.gii',hemI{h})));
            idx=1;
            mask = M.cdata(:,1)>.005;
            data = zeros(size(M.cdata,1),size(seqType,2)*numel(sessN));
            for ss=sessN
                for st=1:size(seqType,2)
                    IN1 = gifti(fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('s%s.%s_%s_exe1.sess-%d.func.gii',hemI{h},metric,seqType{st},ss)));
                    IN2 = gifti(fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('s%s.%s_%s_exe2.sess-%d.func.gii',hemI{h},metric,seqType{st},ss)));
                    for i=1:max(Ico.cdata)
                        if sum(mask(Ico.cdata==i,:))/length(mask(Ico.cdata==i,:))>thres
                            data_ratio = nanmean(nanmean(IN2.cdata(Ico.cdata==i,:)))/nanmean(nanmean(IN1.cdata(Ico.cdata==i,:)));
                            data(Ico.cdata==i,idx) = data_ratio;
                        end
                    end
                    colName{idx} = sprintf('%s_%s_sess-%d',metric,seqType{st},ss);
                    idx=idx+1;
                end
                fprintf('Done %s, sess%d.\n',hemName{h},ss);    
            end
            G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',colName);
            if smooth
                outfile = fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('s%s.RepSup_%s_ratio_ico-%d.func.gii',hemI{h},metric,nTessel));
            else
                outfile = fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('%s.RepSup_%s_ratio_ico-%d.func.gii',hemI{h},metric,nTessel));
            end
            save(G,outfile);
        end
    case 'SURF_wb:create_mask'
        % here create a mask with distances >0.02, sessN = 4 (for further analyses)
        surfaceGroupDir = fullfile(wbDir,'FoSEx_FS_LR_164');
        G = gifti(fullfile(surfaceGroupDir,'sL.summary_dist_RepSup.func.gii'));
        mask = G.cdata(:,1)>0.02;
        outfile = fullfile(surfaceGroupDir,'sL.mask_distance.func.gii');
        G2 = surf_makeFuncGifti(mask,'anatomicalStruct','CortexLeft','columnNames',{'distance_mask'});
        save(G2,outfile);
 
    case 'BETA_get'                                                     % STEP 5.6   :  Harvest betas from rois (raw, univ, multiv prewhit)    
        sessN = 1:4;
        sn  = [5:9,11:31];    
        roi = 1:16; % [1:3,7,8]
        roiDefine = 'all'; % determine all regions from region file R - 'all', 'subset', 'tesselSelect'
        parcelType = 'Brodmann'; % 162tessels, Brodmann, cortex_buckner, BG-striatum, thalamus, tesselsWB_162
        nTessel = 162; % change that - used for tessels
        excludeError = 1; % whether to take error glm or not
        vararginoptions(varargin,{'sn','sessN','roi','parcelType','roiDefine','excludeError','nTessel'});
       
        if strcmp(parcelType,'Yokoi_clusters')
            regType = [1 2 3 4 5 6 7 8 9 10 1 2 3 6 8 9 10];
            regSide = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2];
        end
        % which glm to consider
        if excludeError
            glmDir = glmFoSExErrorDir;
        else
            glmDir = glmFoSExDir;
        end
        idxRoi = 0;
        for ss=sessN            
            % harvest
            for s=sn % for each subj
                tElapsed = tic;
                T=[];
                fprintf('\nSubject: %d\n',s) % output to user
                % load files
                load(fullfile(glmDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)]));  % load subject's region parcellation (R)
                
                if idxRoi==0 % determine for first subject
                    if strcmp(roiDefine,'all')
                        roi=1:size(R,2);
                        nReg_hemi = size(R,2)/2;
                    elseif strcmp(roiDefine,'tesselSelect') % used for WB tessels
                        nReg_hemi = size(R,2)/2;
                        roi = [];
                        for h=1%:2
                            tessels{h}=sml1_imana_repsup('TESSEL:select','hemi',h,'nTessel',nTessel);
                            roi = [roi nReg_hemi*(h-1)+tessels{h}];
                        end
                        nReg_hemi = max(roi)+1;
                    elseif strcmp(roiDefine,'mask')
                        nReg_hemi = 2;
                        roi = 1;
                    elseif strcmp(roiDefine,'subset')
                        roi = roi;
                        nReg_hemi = max(roi)+1;
                    end
                    idxRoi = 1;
                end
                cd(fullfile(glmDir{ss},subj_name{s})); % maybe change back when remove glm
                P=SPM.Vbeta(SPM.xX.iC);
                
                % Add a few extra images
                %----task against rest
                O{1}=sprintf('psc_sess%d_TrainSeq_1st.nii',ss); %psc trained - 1st execution
                O{2}=sprintf('psc_sess%d_UntrainSeq_1st.nii',ss); %psc untrained - 1st execution
                O{3}=sprintf('psc_sess%d_TrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                O{4}=sprintf('psc_sess%d_UntrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                oP=spm_vol(char(O));
                
                if ~exist(fileparts(SPM.xY.VY(1).fname))
                    SPM = sml1_imana_repsup('HOUSEKEEPING:renameSPM','sn',s,'sessN',ss);
                end
                V = SPM.xY.VY;
                if excludeError % new
                    TI = load('SPM_info.mat');
                    idx_noErr = [TI.isError==0;ones(numel(unique(TI.run)),1)];
                end
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
                        if strcmp(parcelType,'BG-striatum')
                            for ss1=1:4
                                load(fullfile(glmFoSExDir{ss1},subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                                V_test = SPM.xY.VY;
                                dataTest{ss1} = region_getdata(V_test,R{r});
                            end
                            i1 = find(dataTest{1}(1,:)); i2 = find(dataTest{2}(1,:)); i3 = find(dataTest{3}(1,:)); i4 = find(dataTest{4}(1,:));
                            iA = intersect(i1,i2); iA = intersect(iA,i3); idx = intersect(iA,i4);
                        else
                            idx = find(Y(1,:));
                        end
                        % estimate region betas
                        if excludeError
                            [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y(:,idx),SPM,'normmode','runwise');
                            S.betaW                   = {betaW(idx_noErr==1,:)}; % exclude errors, keep intercept multivariate pw
                            S.betaUW                  = {bsxfun(@rdivide,beta(idx_noErr==1,:),sqrt(resMS))}; % univariate pw
                            S.betaRAW                 = {beta(idx_noErr==1,:)};
                        else
                            [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y(:,idx),SPM,'normmode','overall');
                            S.betaW                   = {betaW}; % exclude errors, keep intercept multivariate pw
                            S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                            S.betaRAW                 = {beta};
                        end
                        S.resMS                   = {resMS};
                        % info from maps for surface
                        S.psc_train_1st   = {data(1,idx)};
                        S.psc_untrain_1st = {data(2,idx)};
                        S.psc_train_2nd   = {data(3,idx)};
                        S.psc_untrain_2nd = {data(4,idx)};
                    end
                    % voxel position
                    S.volcoord = {R{r}.data'};
                    S.SN                      = s;
                    S.region                  = r;
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
                    fprintf('%d.',r)
                    clear idx;
                end
                dircheck(fullfile(betaDir,subj_name{s}));
                if excludeError
                    save(fullfile(betaDir,subj_name{s},sprintf('betas_ErrFoSEx_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'-struct','T');
                else
                    save(fullfile(betaDir,subj_name{s},sprintf('betas_FoSEx_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'-struct','T');
                end
                clear idx_noErr;
                fprintf('\nDone beta extraction for sess%d-%s.\t',ss,subj_name{s}); toc(tElapsed);
            end            
        end
    case 'BETA_combineGroup'
        % combine individual subject beta structures into the whole
        % structure
        sessN=4;
        sn=[5:9,11:31];
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        excludeError = 1; % whether to take error glm or not
        parcelType='Brodmann'; % Brodmann, BG-striatum, tesselsWB_162
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type','parcelType','excludeError'});
        
        for ss=sessN
            switch(type)
                case 'new'
                    T=[];
                case 'add'
                    T=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            end
            fprintf('subjects added for sess-%d:\n',ss);
            for s=sn
                if excludeError
                    S=load(fullfile(betaDir,subj_name{s},sprintf('betas_ErrFoSEx_%s_%s_sess%d',parcelType,subj_name{s},ss)));
                else
                    S=load(fullfile(betaDir,subj_name{s},sprintf('betas_FoSEx_%s_%s_sess%d',parcelType,subj_name{s},ss)));
                end
                % for BG-striatum
%                 if s>11
%                     S.regType(S.region==2)=2;
%                     S.regSide(S.region==2)=1;
%                     % save per subject
%                     save(fullfile(betaDir,subj_name{s},sprintf('betas_FoSEx_%s_%s_sess%d',parcelType,subj_name{s},ss)),'-struct','S');
%                     % end of BG-striatum
%                 end
                T=addstruct(T,S);
                fprintf('%d.',s);
            end
            if excludeError
                save(fullfile(betaDir,'group',sprintf('betas_ErrFoSEx_%s_sess%d.mat',parcelType,ss)),'-struct','T','-v7.3');
            else
                save(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)),'-struct','T','-v7.3');
            end
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
        sessN = 4;
        sn  = [5:9,11:31];
        roi = [1:14];
        roiDefine = 'all'; % determine from region file
        betaChoice = 'multi'; % uni, multi or raw
        type='new';
        parcelType = 'Brodmann'; % tesselsWB_162, Brodmann, BG-striatum
        excludeError = 1;
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type','parcelType','roiDefine','excludeError'});
            
        for ss=sessN    
            if excludeError
                T = load(fullfile(betaDir,'group',sprintf('betas_ErrFoSEx_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
                glmDir = glmFoSExErrorDir;
            else 
                T = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
                glmDir = glmFoSExDir;
            end
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
                D = load(fullfile(glmDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                if excludeError
                    D = getrow(D,D.isError==0);
                end
                fprintf('\nSubject: %d session: %d\n',s, ss)
                num_run = numruns_task_sess;
                
                for r = roi % for each region
                    S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                    fprintf('%d.',r)
                    
                    for exe = 1:2   % FoSEx
                        % check first if betas are only nan
                        if size(S.betaW{1},1)==1 && isnan(S.betaW{1})
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
                            if excludeError
                                % for betas - exclude intercept(last 8 betas)
                                [G,Sig] = pcm_estGCrossval(betaW,D_exe.run,D_exe.seqNumb);
                                [G_train, Sig_train] = pcm_estGCrossval(betaW(D_exe.seqType==1,:),...
                                    D_exe.run(D_exe.seqType==1),D_exe.seqNumb(D_exe.seqType==1));
                                [G_untrain, Sig_untrain] = pcm_estGCrossval(betaW(D_exe.seqType==2,:),...
                                    D_exe.run(D_exe.seqType==2),D_exe.seqNumb(D_exe.seqType==2));
                                % squared distances
                                So.RDM_nocv = distance_euclidean(betaW',D_exe.seqNumb)';
                                ind12 = indicatorMatrix('allpairs',[1:12]); % All
                                So.RDM = (diag(ind12*G*ind12'))';
                            else
                                [G,Sig] = pcm_estGCrossval(betaW(1:(12*num_run),:),D_exe.run,D_exe.seqNumb);
                                [G_train, Sig_train] = pcm_estGCrossval(betaW(D_exe.seqType==1,:),D.run(D_exe.seqType==1),D.seqNumb(D_exe.seqType==1));
                                [G_untrain, Sig_untrain] = pcm_estGCrossval(betaW(D_exe.seqType==2,:),D.run(D_exe.seqType==2),D.seqNumb(D_exe.seqType==2));
                                % squared distances
                                So.RDM_nocv = distance_euclidean(betaW',D_exe.seqNumb)';
                                So.RDM      = rsa.distanceLDC(betaW,D_exe.run,D_exe.seqNumb);
                            end
                            So.IPM      = rsa_vectorizeIPM(G);
                            So.Sig      = rsa_vectorizeIPM(Sig);
                            So.IPMfull  = rsa_vectorizeIPMfull(G);

                            % trained seq
                            H=eye(6)-ones(6,6)./6;  % centering matrix! 
                            G_trainCent = H*G_train*H;  % double centered G matrix - rows and columns
                            So.eigTrain = sort(eig(G_trainCent)','descend');    % sorted eigenvalues
                            % untrained seq   
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
                            So.regType  = S.regType;
                            So.regSide  = S.regSide;
                            %if r<(numel(roi)/2)+1
                            %    So.regSide = 1;
                            %    So.regType = So.region;
                            %else
                           %     So.regSide = 2;
                          %      So.regType = So.region-(numel(roi)/2);
                           % end
                            So.FoSEx    = exe;
                            To          = addstruct(To,So);
                        end
                    end; % FoSEx
                end; % each region
            end; % each subject
            
            % % save
            if excludeError
                save(fullfile(betaDir,'group',sprintf('stats_ErrFoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)),'-struct','To');
            else
                save(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)),'-struct','To');
            end
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
        sn  = [5:9,11:31];
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
    case 'TESSEL:select_old'
        % here select based on the distance mask
        hemi=1;
        nTessel = 362;
        vararginoptions(varargin,{'hemi','nTessel','betaChoice','sessN'});
        % load in distances and icosahedron
        I = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{hemi})));
        G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('s%s.summary_dist.func.gii',hemI{hemi})));
        maskDist = G.cdata(:,17)>2.48; % all vertices where distances (all) in sess1 significant with p<.01
        tessels = unique(I.cdata)';
        tessels = tessels(tessels~=0); % exclude 0 - medial wall
        choice = [];
        for i=tessels
            numAll = sum(I.cdata==i);
            distPres = sum(maskDist(I.cdata==i)); % how many times significant distance (dist presence) 
            if distPres>(numAll*.25) % needs to be present in 25% of the tessel
                choice = [choice,i];
            end
        end
        varargout={double(choice)};
    case 'TESSEL:select'
        hemi=1;
        nTessel = 362;
        vararginoptions(varargin,{'hemi','nTessel','betaChoice','sessN'});
        % load in distances and icosahedron
        I = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{hemi})));
        G = gifti(fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('s%s.mask_distance.func.gii',hemI{hemi})));
        maskDist = G.cdata; % all vertices where distances (all) in sess1 significant with p<.01
        tessels = unique(I.cdata)';
        tessels = tessels(tessels~=0); % exclude 0 - medial wall
        choice = [];
        for i=tessels
            numAll = sum(I.cdata==i);
            distPres = sum(maskDist(I.cdata==i)); % how many times significant distance (dist presence) 
            if distPres>(numAll*.27) % needs to be present in 27% of the tessel
                choice = [choice,i];
            end
        end
        varargout={double(choice)};

    case 'CALC_act_dist' %-------------------------- ACTIVATION, MAHALANOBIS DIST  ---------------------- this final code 
        % summary structure for 1st / 2nd execution - act, mahalanobis dist
        sn = [5:9,11:31];
        roi = [1:3,7,8];
        sessN = 1:4;
        betaChoice = 'multiPW';
        parcelType = 'Brodmann';  % 162tessels, Brodmann, cortex_buckner
        subjfig = 0;
        excludeError = 1;
        vararginoptions(varargin,{'sn','roi','seq','sessN','betaChoice','fig','parcelType','excludeError'});

        Stats = [];
        
        for ss = sessN % do per session number
            if excludeError
                D = load(fullfile(betaDir,'group',sprintf('stats_ErrFoSEx_%s_%s_sess%d.mat',parcelType,betaChoice,ss))); % loads region data (D)
                T = load(fullfile(betaDir,'group',sprintf('betas_ErrFoSEx_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
                glmDir = glmFoSExErrorDir;
            else
                D = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%s_sess%d.mat',parcelType,betaChoice,ss))); % loads region data (D)
                T = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
                glmDir = glmFoSExDir;
            end

            runs=1:numruns_task_sess;

            indx_train = 1:6;
            indx_untrain = 7:12;
            
            for s=1:length(sn)
                SI = load(fullfile(glmDir{ss}, subj_name{sn(s)}, 'SPM_info.mat'));   % load subject's trial structure
                for r=roi
                    for exe=1:2;    % FoSEx
                        D_exe = getrow(D,D.FoSEx==exe & D.region==r & D.SN==sn(s));
                        % if empty - skip (Buckner network 5 for some subjects)
                        if size(D_exe.IPM,1)~=0
                            if excludeError
                                switch (betaChoice)
                                    case 'uniPW'
                                        beta = T.betaUW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe & SI.isError==0,:);
                                    case 'multiPW'
                                        beta = T.betaW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe & SI.isError==0,:);
                                    case 'raw'
                                        beta = T.betaRAW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe & SI.isError==0,:);
                                end
                                conditionVec = SI.seqNumb(SI.FoSEx==exe & SI.isError==0);
                            else
                                switch (betaChoice)
                                    case 'uniPW'
                                        beta = T.betaUW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe,:);
                                    case 'multiPW'
                                        beta = T.betaW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe,:);
                                    case 'raw'
                                        beta = T.betaRAW{T.SN==sn(s) & T.region==r}(SI.FoSEx==exe,:);
                                end
                                conditionVec = kron(ones(numel(runs),1),[1:12]');
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
                fprintf('Done sess%d %s\n',ss,subj_name{sn(s)});
            end; % sn
        end

        % % save 
        if excludeError
            save(fullfile(repSupDir,sprintf('RepSup_noError_stats_%s_%s.mat',parcelType,betaChoice)),'-struct','Stats');
        else
            save(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)),'-struct','Stats');
        end
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
       excludeError=1;
       sn=[5:9,11:31];
       sessN=1:4;
       vararginoptions(varargin,{'betaChoice','roi','parcelType','hemi','excludeError','sn','sessN'});
       if excludeError
           S = load(fullfile(repSupDir,sprintf('RepSup_noError_stats_%s_%s.mat',parcelType,betaChoice)));
       else
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
       end
       S = getrow(S,ismember(S.sn,sn)&ismember(S.sessN,sessN));
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
               %plt.line([S.speed S.sessN],var,'split',S.FoSEx,'subset',S.regType==f & S.regSide==hemi,'style',style,'leg',{'1st','2nd'});
               plt.line(S.sessN,var,'split',S.FoSEx,'subset',S.regType==f & S.regSide==hemi,'style',style,'leg',{'1st','2nd'});
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
       parcelType = 'Brodmann';
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
           
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
       
       %S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
       S = load(fullfile(repSupDir,sprintf('RepSup_noError_stats_%s_%s.mat',parcelType,betaChoice)));
       figure % trained
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.dist_train,'split',S.FoSEx,'subset',S.regType==roi(f) & S.regSide==hemi,'style',styTrained_exe,'leg',{'1st','2nd'});
           plt.match('y');
           if f==1
               ylabel('Distances trained'); xlabel('Session');
           else
               ylabel('');
           end
           drawline(0,'dir','horz');
           title(sprintf('%s',regname{f}))
       end
      
       figure % untrained
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line([S.speed S.sessN],S.dist_untrain,'split',S.FoSEx,'subset',S.regType==roi(f)&S.regSide==hemi,'style',styUntrained_exe,'leg',{'1st','2nd'});
           plt.match('y');
           drawline(0,'dir','horz');
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
           plt.line([S.speed S.sessN],S.dist_cross,'split',S.FoSEx,'subset',S.regType==roi(f)&S.regSide==hemi,'style',stySeqType_exe,'leg',{'1st','2nd'});
           plt.match('y');
           if f==1
               ylabel('Distance between seq sets'); xlabel('Session');
           else
               ylabel('');
           end
           drawline(0,'dir','horz');
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
       parcelType = 'Brodmann';
       roi=[1:8];
       exe=[1:2];
       vararginoptions(varargin,{'betaChoice','roi','exe'});
       
       exe_label={'1st','2nd'};
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
       
       for r=1:numel(roi)
           for e = exe
               if e==1
                   style=sty_exe1;
               else
                   style=sty_exe2;
               end
               figure(e)
               subplot(1,numel(roi),r);
               plt.line([[S.speed;S.speed] [S.sessN;S.sessN]],[S.psc_train;S.psc_untrain],'split',[ones(size(S.sn));ones(size(S.sn))*2],'subset',[S.FoSEx;S.FoSEx]==e & [S.roi;S.roi]==roi(r),'style',style,'leg',{'trained','untrained'});
               hold on; drawline(0,'dir','horz');
                plt.match('y');
               if r==1
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
           for e=exe
               if e==1
                   style=sty_exe1;
               else
                   style=sty_exe2;
               end
               figure(e)
               subplot(1,numel(roi),r);
               plt.line([[S.speed;S.speed] [S.sessN;S.sessN]],[S.dist_train;S.dist_untrain],'split',[ones(size(S.sn));ones(size(S.sn))*2],'subset',[S.FoSEx;S.FoSEx]==e & [S.roi;S.roi]==roi(r),'style',style,'leg',{'trained','untrained'});
               hold on; drawline(0,'dir','horz');
                plt.match('y');
               if r==1
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
        parcelType = 'Brodmann';
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        seqType={'trained','untrained'};
        
        D = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
        
        T.psc = [D.psc_train; D.psc_untrain];
        T.seqType = [ones(size(D.psc_train));ones(size(D.psc_train))*2];
        T.sn = [D.sn; D.sn];
        T.roi = [D.roi; D.roi];
        T.sessN = [D.sessN; D.sessN];
        T.FoSEx = [D.FoSEx; D.FoSEx];

        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n PSC in %s \n',regname{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.FoSEx T.seqType],{'session','repsup','seqType'},'subset',T.roi==r & ismember(T.sessN,sessN));
            fprintf('\n trained \n');
            for ss=sessN
                fprintf('\nSession %d 1st>2nd \n',ss)
                ttestDirect(T.psc,[T.FoSEx,T.sn],2,'paired','subset',T.roi==r&T.sessN==ss,'split',T.seqType);
                fprintf('\n Trained vs. untrained for each exe \n')
                ttestDirect(T.psc,[T.seqType,T.sn],2,'paired','subset',T.roi==r&T.sessN==ss,'split',T.FoSEx);
            end
        end
    case 'STATS_dist'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        D = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
        
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
                fprintf('\n%s - sess%d: post-hoc t-test on the effect of seqType on RS ratio \n',regname{r},ss);
                ttestDirect(T.dist,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
                fprintf('\n%s - sess%d: post-hoc t-test 1st vs 2nd for each sequence type\n',regname{r},ss);
                ttestDirect(T.dist,[T.FoSEx T.sn],2,'paired','subset',T.roi==r&T.sessN==ss,'split',T.seqType);      
                fprintf('\n%s - sess%d: post-hoc t-test for trained vs. untrained for each execution\n',regname{r},ss);
                ttestDirect(T.dist,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss,'split',T.FoSEx);    
            end
        end        
    case 'STATS_dist_seqType'
        betaChoice = 'multiPW';
        parcelType = 'Brodmann';
        roi=[1:8];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
        
           % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.dist_cross,T.sn,'within',[T.sessN T.FoSEx],{'session','repsup'},'subset',T.roi==r);
        end
        
    case 'PLOT_repsup_dist' % just trying different plots out
        % here plot the repetition suppression in activation
        % and distances (overall)
        parcelType = 'Brodmann';
        roi = [2,3,7];
        sessN = [1,4];
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'parcelType','roi','sessN'});
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
        idx=1;
        figure
        S = normData(S,'psc_train');
        S = normData(S,'psc_untrain');
        ylimA = zeros(size(roi,2)*size(sessN,2),1);
        for r=roi
            for ss=sessN
                subplot(1,numel(roi)*2,idx)
                t = getrow(S,S.sessN==ss&S.roi==r);
                style=styTrained_exe;
                plt.line(t.FoSEx,t.normpsc_train,'style',style);
                hold on;
                style=styUntrained_exe;
                plt.line(t.FoSEx,t.normpsc_untrain,'style',style);
                title(sprintf('%s - sess%d',regname_cortex{r},ss));
                A = gca;
                ylimA(idx) = A.YLim(2);
                ylim([0 max(ylimA)]);
                idx=idx+1;
            end
        end
        figure
        idx=1;
        colblue=[49,130,189]/255;
        collightblue=[158,202,225]/255;
        colred=[222,45,38]/255;
        collightred=[251,177,168]/255;
        for r=roi
            for ss=sessN
                subplot(1,numel(roi)*2,idx)
                t = getrow(S,S.sessN==ss&S.roi==r);
                distributionPlot(t.normpsc_train(t.FoSEx==1),'histOri','left','color',colred,'widthDiv',[2 1],'xValues',1,'showMM',2,'globalNorm',0);
                distributionPlot(t.normpsc_train(t.FoSEx==2),'histOri','right','color',collightred,'widthDiv',[2 2],'xValues',1,'showMM',2,'globalNorm',0);
                distributionPlot(t.normpsc_train(t.FoSEx==1),'histOri','left','color',colblue,'widthDiv',[2 1],'xValues',2,'showMM',2,'globalNorm',0);
                distributionPlot(t.normpsc_untrain(t.FoSEx==2),'histOri','right','color',collightblue,'widthDiv',[2 2],'xValues',2,'showMM',2,'globalNorm',0);
                
                idx=idx+1;
            end
        end
        figure
        idx=1;
        for r=roi
            for ss=sessN
                subplot(1,numel(roi)*2,idx)
                t = getrow(S,S.sessN==ss&S.roi==r);
                distributionPlot(t.normpsc_train(t.FoSEx==1),'histOri','left','color',colred,'widthDiv',[2 1],'xValues',1,'showMM',2,'globalNorm',0);
                distributionPlot(t.normpsc_train(t.FoSEx==2),'histOri','left','color',collightred,'widthDiv',[2 1],'xValues',2,'showMM',2,'globalNorm',0);
                distributionPlot(t.normpsc_train(t.FoSEx==1),'histOri','right','color',colblue,'widthDiv',[2 2],'xValues',1,'showMM',2,'globalNorm',0);
                distributionPlot(t.normpsc_untrain(t.FoSEx==2),'histOri','right','color',collightblue,'widthDiv',[2 2],'xValues',2,'showMM',2,'globalNorm',0);
                
                idx=idx+1;
            end
        end
        D = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));
        d = getrow(D,ismember(D.roi,roi)&ismember(D.sessN,sessN));
        idx=1;
        figure
        for r=roi
            for ss=sessN
                subplot(2,numel(roi)*2,idx)
                plt.box(d.seq,d.subtr_psc,'subset',d.roi==r&d.sessN&ss);
                idx=idx+1;
            end
        end
        
        
        D=load(fullfile(baseDir,'dist_psc_stats',sprintf('dist_%s_ROI',parcelType)));
        T.dist      = [D.dist_train; D.dist_untrain];
        T.sn        = [D.sn; D.sn];
        T.regType   = [D.regType; D.regType];
        T.regSide   = [D.regSide; D.regSide];
        T.sessN     = [D.sessN; D.sessN];
        T.seqType   = [ones(size(D.dist_train));ones(size(D.dist_train))*2];
        T = normData(T,'dist');
        figure
        idx=1;
        for r=roi
            subplot(1,numel(roi),idx)
            plt.bar(T.sessN,ssqrt(T.normdist),'split',T.seqType,'subset',T.regType==r & T.regSide==1 & ismember(T.sessN,sessN),'style',stySeq);
            ylabel('Distance'); xlabel('Session'); title(regname_cortex{r});
            idx=idx+1;
        end
        keyboard;
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
        sn=[5:9,11:31];
        sessN=[1:4];
        exe=[1,2];      
        vararginoptions(varargin,{'sessN'});
        
        missFile = 0;
        for ss=sessN
            for s=sn
                for e=exe
                    LDC_file = fullfile(glmFoSExDir{ss},subj_name{s},sprintf('%s_sess%d_exe%d_LDC.nii',subj_name{s},ss,e)); % searchlight nifti
                    if exist(LDC_file)
                        fprintf('Searchlight exists:\t%s-sess-%d-exe%d\n',subj_name{s},ss,e);
                    else
                        fprintf('Searchlight missing:\t%s-sess-%d-exe%d\t!!!\n',subj_name{s},ss,e);
                        missFile=missFile+1;
                    end
                end
            end
        end
        fprintf('altogether missing %d searchlights\n',missFile);
    case 'SEARCH_copy'
        % copy from the server locally
        sn=[5:9,11:31];
        sessN=1:3;
        exe=[1,2];      
        vararginoptions(varargin,{'sessN'});
        
        for ss=sessN
            fprintf('Working on session %d:\n',ss);
            for s=sn
                for e=exe
                    inFile = fullfile('/Volumes/MotorControl/data/SuperMotorLearning/glmFoSEx',sprintf('glmFoSEx%d',ss),...
                        subj_name{s},sprintf('%s_sess%d_exe%d_LDC.nii',subj_name{s},ss,e)); % searchlight nifti
                    outFile = fullfile(glmFoSExDir{ss},subj_name{s},sprintf('%s_sess%d_exe%d_LDC.nii',subj_name{s},ss,e)); % searchlight nifti
                    copyfile(inFile,outFile);
                    fprintf('Done %s, exe-%d.\n',subj_name{s},e);
                end
            end
        end
    case 'SEARCH_map'
        sn  = [5:9,11:31];
        sessN = 1:3;
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
    
    case 'SEARCH_dist_surf'  % old searchlights - to caret    
        % call dist function for mapping subject contrasts to surfaces
        sn  = [54:9,11:15,18:23,25];
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
        reg = 1:16;
        sn  = [5:9,11:31];
        sessN = 1:4;
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
        roi=[1:5,7,8];
        subtract_mean=0;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'roi','subtract_mean'});
        
        S=load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)));
        
        Exe={'1st exe','2nd exe'};

        for r=1:numel(roi)
            for exe=1:2
                if exe==1
                   style=sty_exe1;
               else
                   style=sty_exe2;
               end
               figure(exe)
               subplot(1,numel(roi),r);
                plt.line([S.sessN>3 S.sessN],S.corrDist,'style',style,'split',S.seqType,'subset',S.FoSEx==exe&S.roi==roi(r),...
                    'leg',{'trained','untrained'},'leglocation','southeast');
                plt.match('y');
                title(sprintf('%s-%s',regname{roi(r)},Exe{exe}));
                if r==1
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
            plt.line([S.sessN>3 S.sessN],S.corr_seqType,'style',stySeqType_exe,'split',S.FoSEx,'subset',S.regType==r&S.regSide==hemi,...
                'leg',{'1st','2nd'},'leglocation','southeast');
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
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'betaChoice','roi','sessN','subtract_mean'});
        
        T=load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)));
        
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.corrDist,T.sn,'within',[T.sessN T.seqType T.FoSEx],{'session','seqType','repsup'},'subset',T.roi==r);
        end
        
        for r=roi
            for ss=sessN
                fprintf('\n Distance in %s session%d - larger than 0?',regname{r},ss);
                fprintf('\n trained \n');
                ttestDirect(T.corrDist,[T.sn],2,'onesample','subset',T.roi==r&T.seqType==1&T.sessN==ss);
                fprintf('\n untrained \n');
                ttestDirect(T.corrDist,[T.sn],2,'onesample','subset',T.roi==r&T.seqType==2&T.sessN==ss);
                fprintf('\n 1st>2nd \n')
                ttestDirect(T.corrDist,[T.FoSEx,T.sn],2,'paired','subset',T.roi==r&T.sessN==ss,'split',T.seqType);
                fprintf('\n Trained vs. untrained for each exe \n')
                ttestDirect(T.corrDist,[T.seqType,T.sn],2,'paired','subset',T.roi==r&T.sessN==ss,'split',T.FoSEx);
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
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'parcelType'});
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));

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
        parcelType = 'Brodmann';
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));

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
        parcelType = 'Brodmann';
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));

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
           
    case 'CALC_subtract_actDist' %--------------------------- REPSUP AS SUBTRACTION: 1st - 2nd --------------------- this needed final code
       % Grafton version
       parcelType='Brodmann';
       betaChoice = 'multiPW';
       excludeError = 1;
       vararginoptions(varargin,{'betaChoice','roi','parcelType','excludeError'});
       
       if excludeError
           S = load(fullfile(repSupDir,sprintf('RepSup_noError_stats_%s_%s.mat',parcelType,betaChoice))); 
       else
           S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
       end
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
       if excludeError
           save(fullfile(repSupDir,sprintf('RepSup_subtract_noError_actDist_%s_%s.mat',parcelType,betaChoice)),'-struct','P');    
       else
           save(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)),'-struct','P');    
       end
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
    case 'PLOT_subtract_psc' % this needed final
        betaChoice = 'multiPW';
        roi=[1:3,7,8];
        parcelType = 'Brodmann';
        hemi=1;
        excludeError = 1;
        sessN=1:4;
        vararginoptions(varargin,{'betaChoice','roi','hemi','excludeError','sessN'});
        if excludeError
            S = load(fullfile(repSupDir,sprintf('RepSup_subtract_noError_actDist_%s_%s.mat',parcelType,betaChoice)));
        else
            S = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));
        end
        S=getrow(S,ismember(S.sessN,sessN));
       figure
       for f = 1:numel(roi)
           subplot(1,numel(roi),f);
           plt.line(S.sessN,S.subtr_psc,'split',S.seq,'subset',S.regType==roi(f) & S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
           %plt.line([S.speed S.sessN],S.subtr_psc,'split',S.seq,'subset',S.regType==roi(f) & S.regSide==hemi,'style',stySeq,'leg',{'Trained','Untrained'});
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
        sessN=1:3;
        seqType='trained';
        parcelType = 'Brodmann';
        roi=[1:5,7,8];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice))); 
        
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
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});

        T=load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice))); 
        
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
       parcelType = 'Brodmann';
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice))); 
 
       A = getrow(S,S.FoSEx==1);
       
       C.sn=A.sn;
       C.roi=A.roi;
       C.psc_distcross = [S.dist_cross(S.FoSEx==2)-S.dist_cross(S.FoSEx==1)];
       C.sessN=A.sessN;
       C.speed=A.speed;
       
       % % save
       save(fullfile(repSupDir,sprintf('RepSup_subtract_crossseq_%s_%s.mat',parcelType,betaChoice)),'-struct','C');
    case 'PLOT_subtract_distCross'
        % distance between sets of trained and untrained sequences
        betaChoice = 'multiPW';
        roi=[1:8];
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(repSupDir,sprintf('RepSup_subtract_crossseq_%s_%s.mat',parcelType,betaChoice)));
        
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
        parcelType = 'Brodmann';
        subtract_mean = 0;
        roi=[1:8];
        vararginoptions(varargin,{'roi'});
        
        S = load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)));
        
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
        parcelType = 'Brodmann';
        seqLabel={'trained','untrained'};
        excludeError = 1;
        vararginoptions(varargin,{'betaChoice','roi','sessN','excludeError'});
        
        if excludeError
            T = load(fullfile(repSupDir,sprintf('RepSup_subtract_noError_actDist_%s_%s.mat',parcelType,betaChoice)));
        else
            T = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));
        end
        % main ANOVA - session x sequence
        for r=roi
            fprintf('\n PSC in %s \n',regname{r});
            anovaMixed(T.subtr_psc,T.sn,'within',[T.sessN T.seq],{'session','seqType'},'subset',T.roi==r & ismember(T.sessN,sessN));
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
        parcelType = 'Brodmann';
        sessN=[1:4];
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));
        
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
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        
        T=load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice))); 
        
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
        parcelType = 'Brodmann';
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice)));
        
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
    case 'STATS_ratio_dist'
        betaChoice = 'multiPW';
        roi=[1:8];
        sessN=[1:4];
        parcelType = 'Brodmann';
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice)));
        
        % main ANOVA - session x sequence
        for r=roi
            fprintf('\n Dist in %s \n',regname{r});
            anovaMixed(T.psc_dist,T.sn,'within',[T.sessN T.seq],{'session','seqType'},'subset',T.roi==r);
            ttestDirect(T.psc_dist,T.sn,2,'onesample','subset',T.roi==r,'split',[T.seq,T.sessN]);
        end
        
        % t-tests comparing sequence type per session
        for r=roi
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS ratio in sess %d in %s \n',ss,regname{r});
                ttestDirect(T.psc_dist,[T.seq T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end
        
        % t-tests comparing across sessions - separately for T / UT seq
        for st=1:2
            fprintf('\n %s sequence\n',seqLabel{st});
            for sessTr=1:3
                for r=roi
                    fprintf('\n post-hoc t-test on comparison of sessions %d-%d in %s \n',sessTr,sessTr+1,regname{r});
                    ttestDirect(T.psc_dist,[T.sessN T.sn],2,'paired','subset',T.roi==r&ismember(T.sessN,[sessTr,sessTr+1])&T.seq==st);
                end               
            end
        end
    
    case 'SCATTER_plot' % xyplots - final figures (supplementary)!!!
        % plot first vs. second execution
        betaChoice = 'multiPW';
        parcelType = 'Brodmann';
        roi=1:8;
        metric = 'psc';
        vararginoptions(varargin,{'betaChoice','roi','exe','metric','plotSeq','parcelType'});
            
        switch metric
            case 'psc'
                margin = 0.08;
            case 'dist'
                margin = 0.009;
        end
        
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
        S1.psc      = [S.psc_train;S.psc_untrain];
        S1.seqType  = [ones(size(S.psc_train));ones(size(S.psc_train))*2];
        S1.FoSEx    = [S.FoSEx;S.FoSEx];
        S1.dist     = [S.dist_train;S.dist_untrain];
        S1.roi      = [S.roi; S.roi];
        S1.regSide  = [S.regSide; S.regSide];
        S1.regType  = [S.regType; S.regType];
        S1.sessN    = [S.sessN; S.sessN];
        S1.sn       = [S.sn; S.sn];
        
        figure
        for r=1:numel(roi)
            subplot(1,numel(roi),find(roi==roi(r)))
            T = getrow(S1,S1.roi==roi(r));
          %  t = tapply(T,{'FoSEx','sessN','seqType'},{metric,'mean'});
            
            T1 = getrow(T,T.FoSEx==1);
            T2 = getrow(T,T.FoSEx==2);
            T1 = normData(T1,(metric));
            T1.seqType(T1.seqType==1)=3; % change so that trained is red
            T2 = normData(T2,(metric)); %eval(sprintf('T1.norm%s',metric))
            xyplot(eval(sprintf('T1.norm%s',metric)),eval(sprintf('T2.norm%s',metric)),T1.sessN,'subset',T1.sessN<4,'errorbars','ellipse','style_thickline','split',T1.seqType)
            hold on;
            xyplot(eval(sprintf('T1.norm%s',metric)),eval(sprintf('T2.norm%s',metric)),T1.sessN,'subset',T1.sessN==4,'errorbars','ellipse','style_thickline','split',T1.seqType)
            g = gca;
            maxLim = max(g.XLim(2),g.YLim(2))+margin;
            minLim = min(g.XLim(1),g.YLim(1))-margin;
            xlim([minLim maxLim]); ylim([minLim maxLim]);
            x = linspace(minLim,maxLim);
            plot(x,x,'-k');
            title(regname_cortex{roi(r)});
            if r==1
                xlabel('Execution 1'); ylabel('Execution 2');
            end
        end
    case 'SCATTER_plot_twoSess'
        % plot ony session 1 and 4
        betaChoice = 'multiPW';
        parcelType = 'Brodmann';
        roi=[2,3,7];
        metric = 'psc';
        vararginoptions(varargin,{'betaChoice','roi','exe','metric','plotSeq','parcelType'});
            
        switch metric
            case 'psc'
                margin = 0.08;
            case 'dist'
                margin = 0.009;
        end
        
        S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
        S1.psc      = [S.psc_train;S.psc_untrain];
        S1.seqType  = [ones(size(S.psc_train));ones(size(S.psc_train))*2];
        S1.FoSEx    = [S.FoSEx;S.FoSEx];
        S1.dist     = [S.dist_train;S.dist_untrain];
        S1.roi      = [S.roi; S.roi];
        S1.regSide  = [S.regSide; S.regSide];
        S1.regType  = [S.regType; S.regType];
        S1.sessN    = [S.sessN; S.sessN];
        S1.sn       = [S.sn; S.sn];
        
        figure
        for r=1:numel(roi)
            subplot(1,numel(roi),find(roi==roi(r)))
            T = getrow(S1,S1.roi==roi(r));            
            T1 = getrow(T,T.FoSEx==1);
            T2 = getrow(T,T.FoSEx==2);
            T1 = normData(T1,(metric));
            T1.seqType(T1.seqType==1)=3; % change so that trained is red
            T2 = normData(T2,(metric)); %eval(sprintf('T1.norm%s',metric))
            xyplot(eval(sprintf('T1.norm%s',metric)),eval(sprintf('T2.norm%s',metric)),T1.sessN,'subset',T1.sessN==1,'errorbars','ellipse','style_thickline','split',T1.seqType)
            hold on;
            xyplot(eval(sprintf('T1.norm%s',metric)),eval(sprintf('T2.norm%s',metric)),T1.sessN,'subset',T1.sessN==4,'errorbars','ellipse','style_thickline','split',T1.seqType)
            g = gca;
            maxLim = max(g.XLim(2),g.YLim(2))+margin;
            minLim = min(g.XLim(1),g.YLim(1))-margin;
            xlim([minLim maxLim]); ylim([minLim maxLim]);
            x = linspace(minLim,maxLim);
            plot(x,x,'-k');
            title(regname_cortex{roi(r)});
            if r==1
                xlabel('Execution 1'); ylabel('Execution 2');
            end
        end
        keyboard;
    case 'SCATTER_corrDist_plot'
        % TO DO: do corrDist scatterplot
        parcelType = 'Brodmann';
        roi=1:8;
        vararginoptions(varargin,{'betaChoice','roi','exe','metric','plotSeq','parcelType'});
          
        margin = 0.004;        
        S1 = load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract0.mat',parcelType)));        
        figure
        for r=1:numel(roi)
            subplot(1,numel(roi),find(roi==roi(r)))
            T = getrow(S1,S1.roi==roi(r));            
            T1 = getrow(T,T.FoSEx==1);
            T2 = getrow(T,T.FoSEx==2);
            T1 = normData(T1,'corrDist');
            T1.seqType(T1.seqType==1)=3; % change so that trained is red
            T2 = normData(T2,'corrDist'); %eval(sprintf('T1.norm%s',metric))
            xyplot(T1.normcorrDist,T2.normcorrDist,T1.sessN,'subset',T1.sessN<4,'errorbars','ellipse','style_thickline','split',T1.seqType)
            hold on;
            xyplot(T1.normcorrDist,T2.normcorrDist,T1.sessN,'subset',T1.sessN==4,'errorbars','ellipse','style_thickline','split',T1.seqType)
            g = gca;
            maxLim = max(g.XLim(2),g.YLim(2))+margin;
            minLim = min(g.XLim(1),g.YLim(1))-margin;
            xlim([minLim maxLim]); ylim([minLim maxLim]);
            x = linspace(minLim,maxLim);
            plot(x,x,'-k');
            title(regname_cortex{roi(r)});
            if r==1
                xlabel('Execution 1'); ylabel('Execution 2');
            end
        end
        
    case 'DIST_scaling' % in final figures!!!!
        % here calculate how far off the distances are from 'scaling' model
        betaChoice = 'multiPW';
        parcelType = 'Brodmann';
        roi=[1:3,7,8];
        sessN=4;
        excludeError = 0;
        vararginoptions(varargin,{'betaChoice','roi','exe','metric','plotSeq','parcelType','excludeError'});
        if excludeError
            S = load(fullfile(repSupDir,sprintf('RepSup_noError_stats_%s_%s.mat',parcelType,betaChoice)));
        else
            S = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_%s.mat',parcelType,betaChoice)));
        end
        %S1.psc      = [S.psc_train;S.psc_untrain];
        S1.psc = [S.act_vecLength_train;S.act_vecLength_untrain];
        S1.seqType  = [ones(size(S.psc_train));ones(size(S.psc_train))*2];
        S1.FoSEx    = [S.FoSEx;S.FoSEx];
        S1.dist     = [ssqrt(S.dist_train);ssqrt(S.dist_untrain)];
        %S1.dist     = [S.dist_train;S.dist_untrain];
        S1.roi      = [S.roi; S.roi];
        S1.sessN    = [S.sessN; S.sessN];
        S1.sn       = [S.sn; S.sn];
        % here calculate the scaling
        R = getrow(S1,S1.FoSEx==1);
        R.scale = S1.psc(S1.FoSEx==2)./S1.psc(S1.FoSEx==1);
        R.expDist = S1.dist(S1.FoSEx==1).*R.scale;
        R = rmfield(R,{'dist','FoSEx','psc'});
        S1 = normData(S1,'dist');
        for r=roi
            figure
            subplot(121)
            plt.bar(S1.FoSEx,S1.normdist,'subset',S1.roi==r & S1.seqType==1 & S1.sessN==sessN,'style',styTrained_exe);
            hold on; drawline(mean(R.expDist(R.roi==r&R.sessN==sessN&R.seqType==1)),'dir','horz');
            title(sprintf('%s - trained',regname_cortex{r})); ylabel('Distance');
            subplot(122)
            plt.bar(S1.FoSEx,S1.normdist,'subset',S1.roi==r & S1.seqType==2 & S1.sessN==sessN,'style',styUntrained_exe);
            hold on; drawline(mean(R.expDist(R.roi==r&R.sessN==sessN&R.seqType==2)),'dir','horz');
            plt.match('y'); title(sprintf('%s - untrained',regname_cortex{r})); ylabel('Distance');
            figure
            plt.bar(S1.FoSEx,S1.normdist,'subset',S1.roi==r & S1.sessN==sessN,'split',S1.FoSEx,'style',styRS);
            hold on; drawline(mean(R.expDist(R.roi==r&R.sessN==sessN)),'dir','horz');
            title(sprintf('%s',regname_cortex{r})); ylabel('Distance');
        end
        % compare regions
        S2 = getrow(S1,S1.FoSEx==2);
        S2.diffDist = R.expDist - S2.dist;
        S2 = normData(S2,'diffDist');
        figure
        plt.bar(S2.roi,S2.normdiffDist.*(-1),'split',S2.seqType,'subset',ismember(S2.roi,roi),'style',stySeq,'leg',{'trained','untrained'});
        hold on; drawline(0,'dir','horz');
        ylabel('Difference from expected distance');
        % effect of region: difference M1-PMd, M1-SPLa
        anovaMixed(S2.diffDist,S2.sn,'within',[S2.roi S2.seqType],{'region','seqType'},'subset',ismember(S2.roi,roi)&S2.sessN==sessN);
        anovaMixed(S2.diffDist,S2.sn,'within',[S2.roi S2.seqType],{'region','seqType'},'subset',ismember(S2.roi,[2,3])&S2.sessN==sessN);
        anovaMixed(S2.diffDist,S2.sn,'within',[S2.roi S2.seqType],{'region','seqType'},'subset',ismember(S2.roi,[2,7])&S2.sessN==sessN);
        anovaMixed(S2.diffDist,S2.sn,'within',[S2.roi S2.seqType],{'region','seqType'},'subset',ismember(S2.roi,[3,7])&S2.sessN==sessN);
        ttestDirect(S2.diffDist,[S2.seqType S2.sn],2,'paired','subset',ismember(S2.roi,roi)&S2.sessN==sessN,'split',S2.roi);
        % effect of region: difference M1-PMd, M1-SPLa
        ttestDirect(S2.diffDist,[S2.roi,S2.sn],2,'paired','subset',ismember(S2.roi,[2,3])&S2.sessN==sessN);
        ttestDirect(S2.diffDist,[S2.roi,S2.sn],2,'paired','subset',ismember(S2.roi,[2,7])&S2.sessN==sessN);
        ttestDirect(S2.diffDist,[S2.roi,S2.sn],2,'paired','subset',ismember(S2.roi,[3,7])&S2.sessN==sessN);
        keyboard;
        ttestDirect(S2.diffDist,[S2.sn],2,'onesample','subset',S2.roi==2&S2.sessN==sessN);
        ttestDirect(S2.diffDist,[S2.sn],2,'onesample','subset',S2.roi==3&S2.sessN==sessN);
        ttestDirect(S2.diffDist,[S2.sn],2,'onesample','subset',S2.roi==7&S2.sessN==sessN);
    case 'DIST_scaling_ico'
        % here for the surface
        sessN = 4;
        nTessel = 642; % 164, 362, 642
        vararginoptions(varargin,{'sessN','nTessel'});
        
        parcelType = sprintf('tesselsWB_%d',nTessel);
        betaChoice = 'multi';
        colName = {sprintf('all_ratio_deviation_ico-%d',nTessel),sprintf('trained_ratio_deviation_ico-%d',nTessel),sprintf('untrained_ratio_deviation_ico-%d',nTessel)};
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
            for h=1
                Ico = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{h})));
                data = zeros(size(Ico.cdata,1),size(colName,2));
                for i=unique(T.region)'
                    S = getrow(T,T.region==i);
                    for j=1:size(S.RDM,1)
                        RDM = rsa_squareRDM(S.RDM(j,:));
                        RDM_train = RDM(1:6,1:6);
                        RDM_untrain = RDM(7:12,7:12);
                        S.dist_train(j,1) = mean(rsa_vectorizeRDM(RDM_train));
                        S.dist_untrain(j,1) = mean(rsa_vectorizeRDM(RDM_untrain));
                    end
                    S1.psc      = [S.act_vecLength_train;S.act_vecLength_untrain];
                    S1.seqType  = [ones(size(S.psc_train));ones(size(S.psc_train))*2];
                    S1.FoSEx    = [S.FoSEx;S.FoSEx];
                    S1.dist     = [S.dist_train;S.dist_untrain];
                    S1.region   = [S.region; S.region];
                    S1.SN       = [S.SN; S.SN];
                    scalePsc    = S1.psc(S1.FoSEx==2)./S1.psc(S1.FoSEx==1);
                    expDist     = S1.dist(S1.FoSEx==1).*scalePsc; % expected distance under scaling
                    diffDist    = -(expDist - S1.dist(S1.FoSEx==2));
                    data(Ico.cdata==i,1) = mean(diffDist);  % overall
                    data(Ico.cdata==i,2) = mean(diffDist(1:26));  % trained
                    data(Ico.cdata==i,3) = mean(diffDist(27:end));  % untrained
                end
                G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',colName);
                
                outfile = fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('%s.RepSup_scalingModel_ico-%d.sess-%d.func.gii',hemI{h},nTessel,ss));
                save(G,outfile);
                fprintf('sess-%d %s\n',ss,hemName{h});
            end
        end
        
    case 'RS_DIST' % this final figure !!!!
        % here plot repetition suppression and distances
        reg = [1:3,7,8];
        parcelType = 'Brodmann';
        betaChoice = 'multiPW';
        excludeError = 1;
        vararginoptions(varargin,{'reg','parcelType','betaChoice','excludeError'});
        if excludeError
            S = load(fullfile(repSupDir,sprintf('RepSup_subtract_noError_actDist_%s_%s.mat',parcelType,betaChoice)));
        else
            S = load(fullfile(repSupDir,sprintf('RepSup_subtract_actDist_%s_%s.mat',parcelType,betaChoice)));
        end
        D = load(fullfile(baseDir,'dist_psc_stats',sprintf('dist_%s_ROI',parcelType)));
        D2.dist     = [ssqrt(D.dist_train); ssqrt(D.dist_untrain)];
        D2.sn       = [D.sn; D.sn];
        D2.roi      = [D.roi; D.roi];
        D2.sessN    = [D.sessN; D.sessN];
        D2.seqType  = [ones(size(D.sn));ones(size(D.sn))*2];
        D2.regSide  = [D.regSide;D.regSide];
        idx=1;
        figure;
        for r=reg
            s = getrow(S,S.roi==r & S.regSide==1);
            d = getrow(D2,D2.roi==r & D2.regSide==1);
            s = normData(s,'subtr_psc');
            d = normData(d,'dist');
            subplot(2,3,idx);
            plt.line(s.sessN,-1*(s.normsubtr_psc),'subset',ismember(s.sessN,[1,4]),'split',s.seq,'leg',{'trained','untrained'},'style',stySeq);
            title(sprintf('%s repsup',regname_cortex{r})); xlabel('Session'); ylabel('repsup');
            subplot(2,3,idx+numel(reg))
            plt.line(d.sessN,d.normdist,'subset',ismember(d.sessN,[1,4]),'split',d.seqType,'leg',{'trained','untrained'},'style',stySeq);
            title(sprintf('%s distance',regname_cortex{r})); xlabel('Session'); ylabel('crossnobis');
            idx=idx+1;
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
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sessN','roi','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(repSupDir,sprintf('RepSup_percent_actDist_%s_%s.mat',parcelType,betaChoice))); 
        
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
     
    case 'firstFinger:distances'
        % here compare sequences with same first finger
        % distances amongst them
        % on exe 1 & 2
        parcelType = 'Brodmann';
        sn = [5:9,11:31];
        sessN = 4;
        reg = 1:8;
        vararginoptions(varargin,{'sn','parcelType','sessN','reg'});
        
        FF = [];
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,ss)));
            for r=reg
                for exe=1:2
                    for s=sn
                        t1 = getrow(T,T.region==r&T.SN==s&T.FoSEx==exe);
                        D = rsa_squareRDM(t1.RDM);
                        sameFing_train = [D(1,2),D(3,4),D(5,6)];
                        diffFing_train = [D(1,3),D(1,4),D(1,5),D(1,6),D(2,4),D(2,5),D(2,6),D(3,5),D(3,6),D(4,5),D(4,6)];
                        sameFing_untrain = [D(7,8),D(9,10),D(11,12)];
                        diffFing_untrain = [D(7,9),D(7,10),D(7,11),D(7,12),D(8,10),D(8,11),D(8,12),D(9,11),D(9,12),D(10,11),D(10,12)];
                        F.dist = [mean(sameFing_train);mean(diffFing_train);mean(sameFing_untrain);mean(diffFing_untrain)];
                        F.seqType = [1;1;2;2];
                        F.fingType = [1;2;1;2];
                        F.sn = [s;s;s;s];
                        F.exe = repmat(exe,[4,1]);
                        F.sessN = repmat(ss,[4,1]);
                        F.region = repmat(t1.region,[4,1]);
                        FF = addstruct(FF,F);
                    end
                end
                f1 = getrow(FF,FF.region==r);
                figure
                subplot(121)
                plt.bar(f1.sessN,f1.dist,'split',[f1.seqType,f1.fingType],'subset',f1.exe==1,'style',styRSbeh,'leg',{'trained-same','trained-diff','untrain-same','untrain-diff'});
                title(sprintf('%s - exe1',regname_cortex{r})); ylabel('distance');
                subplot(122)
                plt.bar(f1.sessN,f1.dist,'split',[f1.seqType,f1.fingType],'subset',f1.exe==2,'style',styRSbeh,'leg',{'trained-same','trained-diff','untrain-same','untrain-diff'});
                title(sprintf('%s - exe2',regname_cortex{r})); ylabel('distance');
                fprintf('Exe1:\n');
                anovaMixed(f1.dist,f1.sn,'within',[f1.seqType f1.fingType],{'seqType1','fingType'},'subset',f1.exe==1);
                ttestDirect(f1.dist,[f1.fingType f1.sn],2,'paired','subset',f1.exe==1,'split',f1.seqType);
                fprintf('Exe2:\n');
                anovaMixed(f1.dist,f1.sn,'within',[f1.seqType f1.fingType],{'seqType1','fingType'},'subset',f1.exe==2);
                ttestDirect(f1.dist,[f1.fingType f1.sn],2,'paired','subset',f1.exe==2,'split',f1.seqType);
            end
        end
        keyboard;
%         
%         plt.scatter(f1.dist(f1.seqType==1&f1.fingType==1&f1.exe==1),f1.dist(f1.seqType==1&f1.fingType==2&f1.exe==1)); hold on;
%         g = gca;
%         maxLim = max(g.XLim(2),g.YLim(2));
%         minLim = min(g.XLim(1),g.YLim(1));
%         xlim([minLim maxLim]); ylim([minLim maxLim]);
%         x = linspace(minLim,maxLim);
%         plot(x,x,'-k');
%         subplot(222)
%         plt.scatter(f1.dist(f1.seqType==1&f1.fingType==1&f1.exe==2),f1.dist(f1.seqType==1&f1.fingType==2&f1.exe==2)); hold on;
%         
    
    case 'trial_analysis'
        % here analyse pairs of trials (1-2, 1-2)
        % take into account whether the two pairs of trials start with the
        % same finger
        % specifically - analyse the third trial
        sessN = 4;
        sn = [5:9,11:31];
        reg = 1:8;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sessN','reg','parcelType'});
        for s=sn
            [t_info, t_mismatch] = sml1_imana_repsup('trials_retrieveTrials','sessN',sessN,'sn',s);
            [y_raw,y_filt,y_res,y_hat,y_adj] = sml1_imana_repsup('trials_retrieveBetas','sn',s,'sessN',sessN,'trialIdx',t_info.trial,'reg',reg,'parcelType',parcelType);
            save(fullfile(baseDir,'repsup_trial',sprintf('%s_%s_matching.mat',subj_name{s},parcelType)),'y_raw','y_filt','y_res','y_hat','y_adj','t_info');
            [y_raw,y_filt,y_res,y_hat,y_adj] = sml1_imana_repsup('trials_retrieveBetas','sn',s,'sessN',sessN,'trialIdx',t_mismatch.trial,'reg',reg,'parcelType',parcelType);
            t_info = t_mismatch;
            save(fullfile(baseDir,'repsup_trial',sprintf('%s_%s_mismatching.mat',subj_name{s},parcelType)),'y_raw','y_filt','y_res','y_hat','y_adj','t_info');
            fprintf('Done %s.\n',subj_name{s});
        end     
    case 'trials_retrieveTrials'
        % find subsequent paris of trials starting with same finger or not
        % e.g. 351...-351... 312... vs 351...-351... 132...
        sessN = 4;
        sn = [5:9,11:31];
        vararginoptions(varargin,{'sessN','sn'});
        
        for s=1:numel(sn)
            T = load(fullfile(anaDir,sprintf('sml1_%s.mat',subj_name{sn(s)})));
            t = getrow(T,T.ScanSess==sessN & T.blockType==9);
            t.inst = [1:size(t.day,1)]';
            t1 = getrow(t,t.FoSEx==1);
            t1.run = t1.BN;
            uniRun = unique(t1.run);
            % add run number
            for r=1:8 % 8 functional runs 
                t1.run(t1.run==uniRun(r))=r;
            end
            % determine which subsequent trials are matching
            f_match = find(diff(t1.press0)==0)+1;
            f_mismatch = find(diff(t1.press0))+1;
            TS.trial = t1.inst(f_match); % trial same first finger
            TD.trial = t1.inst(f_mismatch); % trial different first finger
            % add other info
            TS.seqNumb = t1.seqNumb(f_match);
            TS.seqType = t1.seqType(f_match);
            TS.seqType_before = t1.seqType(f_match-1);
            TS.isError = t1.isError(f_match);
            TS.MT = t1.MT(f_match);
            TS.run = t1.run(f_match);
            TD.seqNumb = t1.seqNumb(f_mismatch);
            TD.seqType = t1.seqType(f_mismatch);
            TD.seqType_before = t1.seqType(f_mismatch-1);
            TD.isError = t1.isError(f_mismatch);
            TD.MT = t1.MT(f_mismatch);
            TD.run = t1.run(f_mismatch);
        end
        varargout{1} = TS;
        varargout{2} = TD;
    case 'trials_retrieveBetas'
        % here get the activity
        sn=5;
        sessN=4;
        parcelType = 'Brodmann';
        trialIdx = [];
        pre = 0;
        post = 8;
        reg = 1:8;
        vararginoptions(varargin,{'sn','sessN','parcelType','trialIdx','reg'});
        
        load(fullfile(glmFoSExDir{sessN},subj_name{sn},'SPM.mat'));
        load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{sn},parcelType)));      % This is made in case 'ROI_define'
        % extract data and onsets
        Y = region_getdata(SPM.xY.VY,R(reg)); % here get the data
        O = spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
        % maybe consider filtering here to get y_adj
        for r=1:size(Y,2)
            y_raw = Y{r};
            y_filt = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
            B = SPM.xX.pKX*y_filt;                             
            y_res = spm_sp('r',SPM.xX.xKXs,y_filt);         
            y_hat = SPM.xX.xKXs.X*B;
            y_adj = y_hat + y_res;         
            for i=1:length(trialIdx)
                Y_raw{r}(i,:) = mean(cut(y_raw,pre,round(O.ons(trialIdx(i))),post,'padding','nan'));
                Y_filt{r}(i,:) = mean(cut(y_filt,pre,round(O.ons(trialIdx(i))),post,'padding','nan'));
                Y_res{r}(i,:) = mean(cut(y_res,pre,round(O.ons(trialIdx(i))),post,'padding','nan'));
                Y_hat{r}(i,:) = mean(cut(y_hat,pre,round(O.ons(trialIdx(i))),post,'padding','nan'));
                Y_adj{r}(i,:) = mean(cut(y_adj,pre,round(O.ons(trialIdx(i))),post,'padding','nan'));
            end
            clear B y_raw y_filt y_res y_hat y_adj;
        end
        varargout{1} = Y_raw;
        varargout{2} = Y_filt;
        varargout{3} = Y_res;
        varargout{4} = Y_hat;
        varargout{5} = Y_adj;
    case 'trials_meanAct'
        % here plot mean activation
        sn = [5:9,11:31];
        parcelType = 'Brodmann';
        metric = 'y_adj';
        roi = 1:8;
        vararginoptions(varargin,{'sn','parcelType','roi','metric'});
        
        DD = [];
        for s=sn
            T1 = load(fullfile(baseDir,'repsup_trial',sprintf('%s_%s_matching',subj_name{s},parcelType)));
            T2 = load(fullfile(baseDir,'repsup_trial',sprintf('%s_%s_mismatching',subj_name{s},parcelType)));
            for r=roi
%                 D.psc = [mean(mean(T1.(metric){r}(T1.t_info.isError==0 & T1.t_info.seqType_before==1,:)));mean(mean(T2.(metric){r}(T2.t_info.isError==0 & T2.t_info.seqType_before==1,:)))];
%                 D.psc = [D.psc;mean(mean(T1.(metric){r}(T1.t_info.isError==0 & T1.t_info.seqType_before==2,:)));mean(mean(T2.(metric){r}(T2.t_info.isError==0 & T2.t_info.seqType_before==2,:)))];
%                 D.repType = [1;1;2;2];
%                 D.seqType = [1;2;1;2];
%                 D.sn = [s;s;s;s];
%                 D.roi = [r;r;r;r];
               D.psc = [mean(mean(T1.(metric){r}(T1.t_info.isError==0 & T1.t_info.seqType==1 & T1.t_info.seqType_before==1,:)));...
                    mean(mean(T1.(metric){r}(T1.t_info.isError==0 & T1.t_info.seqType==1 & T1.t_info.seqType_before==2,:)));...
                    mean(mean(T1.(metric){r}(T1.t_info.isError==0 & T1.t_info.seqType==2 & T1.t_info.seqType_before==1,:)));...
                    mean(mean(T1.(metric){r}(T1.t_info.isError==0 & T1.t_info.seqType==2 & T1.t_info.seqType_before==2,:)));...
                    mean(mean(T2.(metric){r}(T2.t_info.isError==0 & T2.t_info.seqType==1 & T2.t_info.seqType_before==1,:)));...
                    mean(mean(T2.(metric){r}(T2.t_info.isError==0 & T2.t_info.seqType==1 & T2.t_info.seqType_before==2,:)));...
                    mean(mean(T2.(metric){r}(T2.t_info.isError==0 & T2.t_info.seqType==2 & T2.t_info.seqType_before==1,:)));...
                    mean(mean(T2.(metric){r}(T2.t_info.isError==0 & T2.t_info.seqType==2 & T2.t_info.seqType_before==2,:)))];
                D.repType = [1;1;1;1;2;2;2;2];
                D.seqType = [1;1;2;2;1;1;2;2];
                D.seqType_before = [1;2;1;2;1;2;1;2];
                D.roi = repmat(r,8,1);
                D.sn = repmat(s,8,1);
                DD=addstruct(DD,D);
            end
        end
        for r=roi
            figure
%             subplot(211)
%             plt.bar(DD.seqType_before,DD.psc,'subset',DD.roi==r&DD.seqType_before==DD.seqType,'split',DD.repType);
%             subplot(212)
%             plt.bar(DD.seqType_before,DD.psc,'subset',DD.roi==r&DD.seqType_before~=DD.seqType,'split',DD.repType);
%             title(regname_cortex{r});
%             anovaMixed(DD.psc,DD.sn,'within',[DD.repType DD.seqType DD.seqType_before],{'repsup_trial','seqType_now','seqType_bef'},'subset',DD.roi==r);
%             ttestDirect(DD.psc,[DD.repType DD.sn],2,'paired','split',DD.seqType,'subset',DD.seqType==DD.seqType_before & DD.roi==r);
%             ttestDirect(DD.psc,[DD.repType DD.sn],2,'paired','split',DD.seqType_before,'subset',DD.seqType==1 & DD.roi==r);
%             ttestDirect(DD.psc,[DD.repType DD.sn],2,'paired','split',DD.seqType_before,'subset',DD.seqType==2 & DD.roi==r);
           % plt.bar(DD.repType,DD.psc,'subset',DD.roi==r,'split',DD.seqType);
           % anovaMixed(DD.psc,DD.sn,'within',[DD.repType DD.seqType],{'repsup_trial','seqType_now'},'subset',DD.roi==r);
           % ttestDirect(DD.psc,[DD.repType DD.sn],2,'paired','subset',DD.roi==r,'split',DD.seqType);
           % title(regname_cortex{r}); ylabel('activation');
           plt.bar(DD.seqType_before,DD.psc,'subset',DD.roi==r,'split',DD.repType,'style',stySeqType_exe,'leg',{'same-firstFinger','diff-firstFinger'});
           title(regname_cortex{r}); ylabel('activation');
           set(gca,'XTick',[1.6,4.8],'XTickLabel',{'trained','untrained'});
           anovaMixed(DD.psc,DD.sn,'within',[DD.repType DD.seqType_before],{'repsup_trial','seqType_bef'},'subset',DD.roi==r);
           ttestDirect(DD.psc,[DD.repType DD.sn],2,'paired','subset',DD.roi==r,'split',DD.seqType_before);
        end
        keyboard;
    
    case 'trialwise_act_ROI' %------------------------------- TRIALWISE ANALYSIS (CLEAN-UP)--------------------
        % to check the model quality of the glm
        sn=[5:9,11:29];
        sessN=1:4;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        
        pre=0;         % How many TRs before the trial onset
        post=8;        % How many TRs after the trial onset
        for ss=sessN
            T=[];  
            for s=sn
                fprintf('Extracting the onsets and events for subject %s, and session %d:\n',subj_name{s},ss);
                load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM.mat'));
                SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s}));
                load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},parcelType)));      % This is made in case 'ROI_define'
                [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
                
                % Create a structure with trial onset and trial type (event)
                D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs        
                % D.event - conditions (seqNumb: 1-12 - first exe; 13-24 - second exe)
                % D.block - run
                % add other indicators - FoSEx, seqNumb, seqType
                D.FoSEx=D.event;
                D.FoSEx(D.FoSEx<13)=1;
                D.FoSEx(D.FoSEx>12)=2;
                D.seqNumb=D.event;
                D.seqNumb(D.seqNumb>12)=D.seqNumb(D.seqNumb>12)-12;
                D.seqType=D.event;
                D.seqType(D.seqNumb<7)=1;
                D.seqType(D.seqNumb>6)=2;
                if  sum(D.FoSEx==1)~=sum(D.FoSEx==2)
                    keyboard;
                end
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
                    fprintf('%d.',r);
                end;
                fprintf('...done.\n');
            end;
            cd(repSupDir);
            save(sprintf('trialwise_ACT_sess%d.mat',ss),'-struct','T');
        end
    case 'trialwise_BEH'
        sn=[5:9,11:29];
        sessN=1:4;
        vararginoptions(varargin,{'sn','sessN'});
        
        for ss=sessN
            RS2=[];
            for s=sn
                % load behavioural data
                B=dload(fullfile(behavDir,sprintf('sml1_%s.dat',subj_name{s})));
                % only that scan = functional runs
                BB = getrow(B,B.ScanSess==ss & (B.blockType==9 | B.blockType==3));
                for r = 1:numruns_task_sess
                    uniqrun=unique(BB.BN);
                    R = getrow(BB,BB.BN==uniqrun(r)); % 1-8 func runs of the session      
                    for c = 1:numel(num_seq) % each sequence
                        idx = find(R.seqNumb==c & R.FoSEx==1);   % find indx of all FIRST trials in run - 1:6 trained; 7-12 untrained
                        % difference in MT - 1st-2nd execution
                        RS.MT_dif=R.MT(idx)-R.MT(idx+1);
                        for i=1:numel(idx)
                            % determine if no response in one of two trials
                            if (R.MT(idx(i))==0 || R.MT(idx(i)+1)==0)
                                RS.noResponse(i,:)=1;
                            else
                                RS.noResponse(i,:)=0;
                            end
                            % determine error type
                            % 1 - no error
                            % 2 - 1st correct, 2nd error
                            % 3 - 1st error, 2nd correct
                            % 4 - both errors
                            if R.isError(idx(i))==0 && R.isError(idx(i)+1)==0
                                RS.errorType(i,:)=1;
                            elseif R.isError(idx(i))==0 && R.isError(idx(i)+1)==1
                                RS.errorType(i,:)=2;
                            elseif R.isError(idx(i))==1 && R.isError(idx(i)+1)==0
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
                        RS.sessN=ones(size(RS.seqType))*ss;
                        RS.run=ones(size(RS.seqType))*r;
                        RS2=addstruct(RS2,RS);
                    end
                end
            end  
            save(fullfile(repSupDir, sprintf('trialwise_BEH_sess%d.mat',ss)),'-struct','RS2');
            fprintf('Done sess-%d.\n',ss);
        end
    case 'trialwise_calc_RS'
        sessN=1:4;
        vararginoptions(varargin,{'sn','sessN'});
        for ss=sessN
            D=load(fullfile(repSupDir,sprintf('trialwise_ACT_sess%d.mat',ss)));
            B = [];
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
            save(fullfile(repSupDir, sprintf('trialwise_RS_sess%d.mat',ss)),'-struct','B');
        end
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
        sn=[5:9,11:29];
        sessN=1:4;
        roi=1:8;
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
    
    case 'PLOT:angle_length'
        % here plot the length as activation, angle as cosine
        roi=[2,3];
        subtract_mean=0;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'roi','subtract_mean'});
        
        S1=load(fullfile(repSupDir,sprintf('corrDist_RS_%s_meanSubtract%d.mat',parcelType,subtract_mean)));
        S2 = load(fullfile(repSupDir,sprintf('RepSup_stats_%s_multiPW.mat',parcelType)));
        for r=roi
            T1 = getrow(S1,S1.roi==r);
            T2 = getrow(S2,S2.roi==r);
            t1 = tapply(T1,{'FoSEx','sessN','seqType'},{'corrDist'});
            t1.corrPlot = t1.corrDist*100;
            t_tr = tapply(T2,{'FoSEx','sessN'},{'psc_train'});
            t_utr = tapply(T2,{'FoSEx','sessN'},{'psc_untrain'});
            t2.psc = [t_tr.psc_train;t_utr.psc_untrain];
            t2.FoSEx = [t_tr.FoSEx;t_utr.FoSEx];
            t2.sessN = [t_tr.sessN;t_utr.sessN];
            t2.seqType = [ones(size(t_tr.sessN));ones(size(t_tr.sessN))*2];
            figure
            for ss=1:4
                subplot(4,2,(ss-1)*2+1) % exe1
                % trained
                plot(-(t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==1)/2),t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==1),'ro');
                hold on;
                plot((t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==1)/2),t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==1),'ro');
                plot([0;-(t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==1)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==1)],'r-');
                plot([0;(t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==1)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==1)],'r-');
                % untrained
                plot(-(t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==2)/2),t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==2),'bo'); 
                plot((t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==2)/2),t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==2),'bo'); 
                plot([0;-(t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==2)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==2)],'b-');
                plot([0;(t1.corrDist(t1.sessN==ss&t1.FoSEx==1&t1.seqType==2)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==1&t2.seqType==2)],'b-');
                title(sprintf('%s - sess%d - exe1',regname_cortex{r},ss));
                
                subplot(4,2,(ss-1)*2+2) % exe2
                % trained
                plot(-(t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==1)/2),t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==1),'ro');
                hold on;
                plot((t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==1)/2),t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==1),'ro');
                plot([0;-(t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==1)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==1)],'r--');
                plot([0;(t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==1)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==1)],'r--');
                % untrained
                plot((t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==2)/2),t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==2),'bo'); 
                plot((-t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==2)/2),t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==2),'bo'); 
                plot([0;(t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==2)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==2)],'b--');
                plot([0;-(t1.corrDist(t1.sessN==ss&t1.FoSEx==2&t1.seqType==2)/2)],[0;t2.psc(t2.sessN==ss&t2.FoSEx==2&t2.seqType==2)],'b--');
                title(sprintf('%s - sess%d - exe2',regname_cortex{r},ss));
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
    
    case 'PCM:correlation_firstFing'
        % correlation with respect to first fingers
        runEffect   = 'random';
        beta_choice = 'mw';
        algorithm   = 'NR'; % minimize or NR
        reg         = [2,3,7];
        parcelType  = 'Brodmann';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessN       = 1:4;      % sessions
        corrLim     = [-1 1];    % bounds for lower / upper correlation of models - change this to [-1 1]!!!
        nModel      = 41;       % number of correlation models (determines how fine grained the corr estimates are)
        sn = [5:9,11:31];
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','withinCov','nModel'});
        
        corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
        withinCov = zeros(6,3);
        withinCov(1,1,1) = 1;
        withinCov(2,2,1) = 1;
        withinCov(3,3,2) = 1;
        withinCov(4,4,2) = 1;
        withinCov(5,5,3) = 1;
        withinCov(6,6,3) = 1;
        AllSess=[];
        for c=1:length(corrS)
            M{c} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',2,'numItems',6,'r',corrS(c));
        end
        for ss = sessN % session transition
            AllReg = [];
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = reg
                for st=1:2 % per seqType
                    tstart = tic;
                    for p=1:length(sn)
                        for rp = 1:2 % repetition
                            glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
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
                            indx = (D.seqType==st & D.FoSEx==rp);
                            if rp == 1
                                condVec{p}  = D.seqNumb(indx==1,:); % conditions
                                partVec{p}  = D.run(indx==1,:);
                                Data{p}     = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                            else
                                condVec{p}  = [condVec{p}; D.seqNumb(indx==1,:)+6]; % treat 2nd session as additional conditions
                                partVec{p}  = [partVec{p}; D.run(indx==1,:)+8];  % runs/partitions of 2nd session as additional runs
                                Data{p}     = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                            end;
                        end; % repetition
                    end; % subj
                    % group fit
                    T = sml1_imana_repsup('PCM:corrModel_run','Data',Data,'model',M,'partV',partVec,'condV',condVec,'algorithm',algorithm,'runEffect',runEffect);
                    % other variables of interest
                    T.roi                   = ones(size(T.SN))*r;
                    T.regType               = ones(size(T.SN))*t.regType;
                    T.regSide               = ones(size(T.SN))*t.regSide;
                    T.seqType               = ones(size(T.SN))*st;
                    T.sessN                 = ones(size(T.SN))*ss;
                    T = rmfield(T,'reg');
                    AllReg = addstruct(AllReg,T);
                    fprintf('Done seqType %d sessN %d- reg: %d/%d\n\n',st,ss,r,length(reg));
                    toc(tstart);
                end; % seqType
                fprintf('Done reg %d/%d:\t \tsess: %d\n\n',r,numel(reg),ss);
            end; % region
            save(fullfile(repSupDir,sprintf('PCM_corrModels_firstFing_%s_nonlinear_session-%d.mat',parcelType,ss)),'-struct','AllReg');
            fprintf('Done all:\tsess: %d-%d\n\n\n\n',ss);
            AllSess = addstruct(AllSess,AllReg);
        end; % session
        % save output
        save(fullfile(repSupDir,sprintf('PCM_corrModels_firstFing_%s_nonlinear.mat',parcelType)),'-struct','AllSess');
    case 'PCM:correlation_shape' % ------------------------- PCM model of all possible correlations ----------------------
        % estimate the likelihood that the data correspond to specific
        % correlations
        % here assessing models of different correlation values
        % across sessions
        runEffect   = 'random';
        beta_choice = 'mw';
        algorithm   = 'NR'; % minimize or NR
        reg         = 1:8;
        parcelType  = 'Brodmann';
        %modelType   = 'wSess';  % options: wSess, noSess, meanPattern
        sessN       = 1:4;      % sessions
        corrLim     = [-1 1];    % bounds for lower / upper correlation of models - change this to [-1 1]!!!
        nModel      = 41;       % number of correlation models (determines how fine grained the corr estimates are)
        withinCov   = 'individual'; % covariance type: individual or iid (or meanPattern)
        sn = [5:9,11:31];
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','withinCov','nModel'});
        AllSess=[];
        % construct models
        corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
        for c=1:length(corrS)
            if strcmp(withinCov,'meanPattern')
                M{c} = pcm_buildCorrModel('type','nonlinear','withinCov','iid','numCond',2,'numItems',6,'r',corrS(c),'condEffect',0);
            else
                M{c} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',2,'numItems',6,'r',corrS(c));
            end
        end
        for ss = sessN % session transition
            AllReg = [];
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = reg
                for st=1:2 % per seqType
                    tstart = tic;
                    condVec = cell(1,length(sn));
                    partVec = condVec; Data = condVec; % initialise
                    for p=1:length(sn)
                        for rp = 1:2 % repetition
                            glmDirSubj=fullfile(glmFoSExDir{sessN(ss)}, subj_name{sn(p)});
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
                            indx = (D.seqType==st & D.FoSEx==rp);
                            if rp == 1
                                condVec{p}  = D.seqNumb(indx==1,:); % conditions
                                partVec{p}  = D.run(indx==1,:);
                                Data{p}     = beta(indx==1,:);  % Data is N x P (cond x voxels) - no intercept
                            else
                                condVec{p}  = [condVec{p}; D.seqNumb(indx==1,:)+6]; % treat 2nd session as additional conditions
                                partVec{p}  = [partVec{p}; D.run(indx==1,:)+8];  % runs/partitions of 2nd session as additional runs
                                Data{p}     = [Data{p}; beta(indx==1,:)];  % Data is N x P (cond x voxels) - no intercept
                            end;
                        end; % repetition
                    end; % subj
                    % group fit
                    T = sml1_imana_repsup('PCM:corrModel_run','Data',Data,'model',M,'partV',partVec,'condV',condVec,'algorithm',algorithm,'runEffect',runEffect);
                    % other variables of interest
                    T.roi                   = ones(size(T.SN))*r;
                    T.regType               = ones(size(T.SN))*t.regType;
                    T.regSide               = ones(size(T.SN))*t.regSide;
                    T.seqType               = ones(size(T.SN))*st;
                    T.sessN                 = ones(size(T.SN))*ss;
                    T = rmfield(T,'reg');
                    AllReg = addstruct(AllReg,T);
                    fprintf('Done modelType: %s seqType %d sessN %d- reg: %d/%d\n\n',withinCov,st,ss,r,length(reg));
                    toc(tstart);
                end; % seqType
                fprintf('Done reg %d/%d:\tmodelType: %s \tsess: %d\n\n',r,numel(reg),withinCov,ss);
            end; % region
            save(fullfile(repSupDir,sprintf('PCM_corrModels_%s_nonlinear_%s_session-%d.mat',parcelType,withinCov,ss)),'-struct','AllReg');
            fprintf('Done all:\tmodelType: %s \tsess: %d-%d\n\n\n\n',withinCov,ss);
            AllSess = addstruct(AllSess,AllReg);
        end; % session
        % save output
        save(fullfile(repSupDir,sprintf('PCM_corrModels_%s_nonlinear_%s.mat',parcelType,withinCov)),'-struct','AllSess');
    case 'PCM:corrModel_run_old'
        algorithm='NR';
        vararginoptions(varargin,{'Data','model','partV','condV','runEffect','algorithm'});
        T = pcm_fitModelGroup(Data,model,partV,condV,'runEffect',runEffect,'fitAlgorithm',algorithm,'fitScale',1); % only group fit
        T.bayesEst = bsxfun(@minus,T.likelihood,T.likelihood(:,1)); % now relative to model with 0 correlation
        T.mean_likelihood = mean(T.likelihood,2);
        %T.norm_likelihood = bsxfun(@minus,T.likelihood,T.mean_likelihood); % normalised to mean
        T.posterior = exp(T.bayesEst);
        T.posterior = bsxfun(@rdivide,T.posterior,sum(T.posterior,2));
        varargout{1}=T;
    case 'PCM:corrModel_plot'
         modelType = 'nonlinear_individual';
        % for nonlinear models: nonlinear_indiviudal / nonlinear_individual, nonlinear_meanPattern
        parcelType = 'Brodmann';
        reg = [2,3,7];
        hemi = 1;
        sessN = 1:4;
        metric = 'bayesEst'; % bayesEst, bayesEst_cross or bayesEst_individ
        vararginoptions(varargin,{'modelType','parcelType','reg','metric','sessType','sessN'});
        
       % style.file(fullfile(codeDir,'sml_style.m'));
      %  style.use('default');
        T=load(fullfile(repSupDir,sprintf('PCM_corrModels_firstFing_%s_nonlinear.mat',parcelType)));

        for ss=sessN
            for r=reg
                t = getrow(T,T.regType==r & T.regSide==hemi & T.sessN==ss);
                t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
                t.metric3 = bsxfun(@minus,t.(metric),max(t.(metric),2));
                % reshape here
                nModel      = size(t.(metric),2);
                D.(metric)  = t.(metric)(:);
                D.metric2   = t.metric2(:);
                D.metric3   = t.metric3(:);
                D.SN        = repmat(t.SN,nModel,1);
                D.sessN    = repmat(t.sessN,nModel,1);
                D.seqType   = repmat(t.seqType,nModel,1);
                D.model     = kron((1:nModel)',ones(size(t.(metric),1),1));
                figure
                %    style.use('stySeq')
                plt.line(D.model,D.metric2,'split',D.seqType,'leg',{'trained','control'},'leglocation','southeast','style',stySeqShade);
                c = gca;
                corrTick = linspace(min(c.XTick),max(c.XTick),21);
                set(gca,'XTickLabel',(-1:.1:1),'XTick',corrTick); hold on; drawline(0,'dir','horz');
                xlabel('Correlation'); ylabel('log-BF'); title(sprintf('%s - sess%d',regname_cortex{r},ss));
                %   for ss=1:4
                %                 subplot(2,2,ss)
                %                 style.use('SeqShade_small')
                %                 plt.line(D.model,D.metric2,'split',D.seqType,'subset',D.sessN==ss,'leg',{'trained','control'},'leglocation','southeast');
                %                 hold on; drawline(0,'dir','horz'); title(sprintf('%s - session %d',regname_cortex{r},ss));
                %                 c = gca;
                %                 corrTick = linspace(min(c.XTick),max(c.XTick),21);
                %                 set(gca,'XTickLabel',(-1:.1:1),'XTick',corrTick);
                %          %   end
                %             figure
                %             subplot(2,2,1)
                %             style.use('Trained_shade')
                %             plt.line(D.model,D.metric2,'split',D.sessN,'subset',D.sessN<4&D.seqType==1,'leg',{'sess-1','sess-2','sess-3'},'leglocation','southeast');
                %             hold on; drawline(0,'dir','horz');
                %             title(sprintf('%s - trained 1-3',regname_cortex{r}));
                %             subplot(2,2,3)
                %             plt.line(D.model,D.metric2,'split',D.sessN,'subset',D.sessN>2&D.seqType==1,'leg',{'sess-3','sess-4'},'leglocation','southeast');
                %             hold on; drawline(0,'dir','horz');
                %             title(sprintf('%s - trained 3-4',regname_cortex{r}));
                %             style.use('Untrained_shade')
                %             subplot(2,2,2)
                %             plt.line(D.model,D.metric2,'split',D.sessN,'subset',D.sessN<4&D.seqType==2,'leg',{'sess-1','sess-2','sess-3'},'leglocation','southeast');
                %             hold on; drawline(0,'dir','horz');
                %             title(sprintf('%s - untrained 1-3',regname_cortex{r}));
                %             subplot(2,2,4)
                %             plt.line(D.model,D.metric2,'split',D.sessN,'subset',D.sessN>2&D.seqType==2,'leg',{'sess-3','sess-4'},'leglocation','southeast');
                %             hold on; drawline(0,'dir','horz');
                %             title(sprintf('%s - untrained 3-4',regname_cortex{r}));
            end
        end
    case 'PCM_simulate_repsupModels' %---------------------- PCM on CORRELATION 1st-2nd REPETITION - old ----------------
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
        reg = [2,3];
        sn=[5:9,11:31];
        sessN=1:4; % need to be two sessions at the time
        models={'generic','specific'};
        parcelType='Brodmann';
        M_type=2; % 1 - generic, 2 - specific
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','M_type','parcelType'})

        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            AllReg=[];
            for st=1:2
                for r = reg
                    for p=1:length(sn)
                        b=getrow(B,B.SN==sn(p)&B.region==r);
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        
                        switch (beta_choice)
                            case 'uw'
                                beta = b.betaUW{:}';
                            case 'mw'
                                beta = b.betaW{:}';
                            case 'raw'
                                beta = b.betaRAW{:}'; % no intercept - use T.betaRAWint otherwise
                        end
                        
                        indx = ones(size(D.run));
                        % conditions - 1-6 for first exe, 7-12 for second
                        cond = D.seqNumb;
                        if st==1
                            % make condVec for second exe into 7-12
                            cond(D.FoSEx==2)=cond(D.FoSEx==2)+6;
                        else
                            % make condVec for first exe into 1-6
                            cond(D.FoSEx==1)=cond(D.FoSEx==1)-6;
                        end
                        condVec{p} = cond(D.seqType==st); % conditions
                        partVec{p} = D.run(D.seqType==st);
                        Data{p} = beta(:,D.seqType==st)';  % Data is N x P (cond x voxels) - no intercept
                    end; % subj
                    
                    % construct models
                    switch M_type
                        case 1 % generic
                            M = pcm_repsupModel_generic;
                        case 2 % specific
                            M = pcm_repsupModel_specific;
                    end
                    T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                    C = pcm_correlation(Data,partVec,condVec,M{4},runEffect,M_type);
                    T.roi = ones(size(T.SN))*r;
                    T.regSide   = ones(size(T.SN))*b.regSide;
                    T.regType   = ones(size(T.SN))*b.regType;
                    T.seqType   = ones(size(T.SN))*st;
                    
                    AllReg=addstruct(AllReg,T);
                    AllReg=addstruct(AllReg,C);
                    AllReg=rmfield(AllReg,{'reg','theta','theta_hat','thetaCr'});
                    fprintf('Done: sess-%d reg-%d/%d\n',ss,r,length(reg));
                end % region
            end; % seqType
            % save output
            save(fullfile(pcmDir,sprintf('PCM_repsup_models_%s_sess%d_%s.mat',models{M_type},ss,parcelType)),'-struct','AllReg');
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
        reg = 1:16;
        parcelType='Brodmann'; % Brodmann
        sn=[5:9,11:31];
        sessN=1:4; % need to be two sessions at the time
        modelType='specific'; %generic or specific
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','seqType','modelType','parcelType'})
        
        % construct models
        switch modelType
            case 'generic'
                M = pcm_corrModel_sharedSeqType;
                mf=1;
            case 'specific'
                M = pcm_corrModel_indSeq_sharedSeqType;
                mf=2;
        end
        for ss=sessN
            CC=[];
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d',parcelType,ss)));
            for stIndx = 1:2 % seqType
                for r = reg
                    for p=1:length(sn)
                        b=getrow(B,B.SN==sn(p)&B.region==r);
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        
                        switch (beta_choice)
                            case 'uw'
                                beta = b.betaUW{:};
                            case 'mw'
                                beta = b.betaW{:}';
                            case 'raw'
                                beta = b.betaRAW{:}'; % no intercept - use T.betaRAWint otherwise
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
                        Data{p} = beta(D.seqType==stIndx,:);
                        fprintf('Extracted data for %s...\n',subj_name{sn(p)});
                    end; % subj
                    C           = pcm_correlation(Data,partVec,condVec,M{2},runEffect,mf);
                    C.SN        = [1:size(C.r_naive,1)]';
                    C.roi       = ones(size(C.r_naive))*r;
                    C.regSide   = ones(size(C.r_naive))*b.regSide;
                    C.regType   = ones(size(C.r_naive))*b.regType;
                    C.seqType   = ones(size(C.r_naive))*stIndx;
                    C           = rmfield(C,'theta');
                    CC          = addstruct(CC,C);
                    fprintf('Done calculating corr for region-%d\n\n\n',r);
                end
            end
            % save output
            save(fullfile(pcmDir,sprintf('PCM_repsup_corr_%s_sess%d_%s.mat',parcelType,ss,modelType)),'-struct','CC');
            % save(fullfile(pcmDir,sprintf('PCM_cov_%s_sess%d_%s_%s.mat',parcelType,ss,seqType{stIndx},modelType)),'-struct','CC');
            fprintf('Done all regions sess-%d\n\n',ss);
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
    
    case 'PCM_constructModelFamily'  % used in the paper
        % November 30th 2019
         % here simpler version of model family:
        % includes first finger, all fingers, seqType, sequences:T, sequences:U
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum, distanceMask
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 1:4;
        reg = 1:8;
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; GG=[]; PP=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe=1:2
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                        end
                       partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions
                       if ~naturalStats
                           m = pcm_defineSequenceModels_fixed_noChunks(Seq,sn(p));
                           %m = pcm_defineSequenceModels_test(Seq,sn(p));
                       else
                           load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                           NatStat = NatStats.G_cent;
                           m = pcm_defineSequenceModels_fixed_noChunks_natStats(Seq,sn(p),NatStat); % make
                           %m = pcm_defineSequenceModels_test_natStats(Seq,sn(p),NatStat); % make
                       end                        
                        [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = D.FoSEx==exe;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    end;
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T.exe = ones(size(T.SN))*exe;
                    
                    P.thetaCr = T.thetaCr;
                    P.G_pred = T.G_pred;
                    P.bayesEst = {T.bayesEst};
                    P.cross_likelihood = {T.cross_likelihood};
                    P.roi = reg(r);
                    P.regType = regType(r);
                    P.regSide = regSide(r);
                    P.sessN = ss;
                    P.exe = exe;
                    T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                    AllReg=addstruct(AllReg,T);
                    PP=addstruct(PP,P);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.sessN = ones(length(K.indx),1)*ss;
                    K.exe = ones(length(K.indx),1)*exe;
                    KK = addstruct(KK,K);
                    idxST = Comb(:,3)==1; % always including SeqType
                    G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,4,5]));
                    %idxST = Comb(:,4)==1; 
                    %G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,3,5,6]));
                    G.roi = ones(length(G.indx),1)*reg(r);
                    G.regType = ones(length(G.indx),1)*regType(r);
                    G.regSide = ones(length(G.indx),1)*regSide(r);
                    G.sessN = ones(length(G.indx),1)*ss;
                    G.exe = ones(length(G.indx),1)*exe;
                    GG = addstruct(GG,G);
                    fprintf('Done sess-%d reg-%d/%d - exe-%d\n',ss,r,length(reg),exe);
                end
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing','SeqType','Trained','Untrained'};
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_simple_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','PP');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_modelFamily_stats_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KK');
    case 'PCM_constructModelFamily_noError'
        % no error glm
        % includes first finger, all fingers, seqType, sequences:T, sequences:U
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum, distanceMask
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 4;
        reg = [1:3,7,8];
        naturalStats = 1;
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; GG=[]; PP=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_ErrFoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe=1:2
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExErrorDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        D=getrow(D,D.isError==0);
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                        end
                       partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions
                       if ~naturalStats
                           m = pcm_defineSequenceModels_fixed_noChunks(Seq,sn(p));
                       else
                           load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                           NatStat = NatStats.G_cent;
                           m = pcm_defineSequenceModels_fixed_noChunks_natStats(Seq,sn(p),NatStat); % make
                       end                        
                        [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        %condVec{p} = repmat(Z,max(D.run),1);
                        cVec = [];
                        for rr=1:max(D.run)
                            DR = getrow(D,D.run==rr & D.FoSEx==exe);
                            cVec = [cVec;Z(DR.seqNumb,:)];
                        end
                        condVec{p} = cVec;
                        indx = D.FoSEx==exe;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    end;
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T.exe = ones(size(T.SN))*exe;
                    
                    P.thetaCr = T.thetaCr;
                    P.G_pred = T.G_pred;
                    P.bayesEst = {T.bayesEst};
                    P.cross_likelihood = {T.cross_likelihood};
                    P.roi = reg(r);
                    P.regType = regType(r);
                    P.regSide = regSide(r);
                    P.sessN = ss;
                    P.exe = exe;
                    T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                    AllReg=addstruct(AllReg,T);
                    PP=addstruct(PP,P);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.sessN = ones(length(K.indx),1)*ss;
                    K.exe = ones(length(K.indx),1)*exe;
                    KK = addstruct(KK,K);
                    idxST = Comb(:,3)==1; % always including SeqType
                    G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,4,5]));
                    G.roi = ones(length(G.indx),1)*reg(r);
                    G.regType = ones(length(G.indx),1)*regType(r);
                    G.regSide = ones(length(G.indx),1)*regSide(r);
                    G.sessN = ones(length(G.indx),1)*ss;
                    G.exe = ones(length(G.indx),1)*exe;
                    GG = addstruct(GG,G);
                    fprintf('Done sess-%d reg-%d/%d - exe-%d\n',ss,r,length(reg),exe);
                end
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing','SeqType','Trained','Untrained'};
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_noError_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_noError_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_noError_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','PP');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_modelFamily_stats_noError_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KK');
    case 'PCM_constructModelFamily_new'
        % this is a simple sequence model
        % trained / untrained together
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum, distanceMask
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 1:4;
        reg = [1:3,7,8];
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; PP=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe=1:2
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                        end
                        partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions    
                        load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                        NatStat = NatStats.G_cent;
                        m = pcm_defineSequenceModels_new(Seq,sn(p),NatStat);
                        [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = D.FoSEx==exe;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    end;
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T.exe = ones(size(T.SN))*exe;
                    
                    P.thetaCr = T.thetaCr;
                    P.G_pred = T.G_pred;
                    P.bayesEst = {T.bayesEst};
                    P.cross_likelihood = {T.cross_likelihood};
                    P.roi = reg(r);
                    P.regType = regType(r);
                    P.regSide = regSide(r);
                    P.sessN = ss;
                    P.exe = exe;
                    T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                    AllReg=addstruct(AllReg,T);
                    PP=addstruct(PP,P);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.sessN = ones(length(K.indx),1)*ss;
                    K.exe = ones(length(K.indx),1)*exe;
                    KK = addstruct(KK,K);
                    fprintf('Done sess-%d reg-%d/%d - exe-%d\n',ss,r,length(reg),exe);
                end
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing','Sequence'};
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_new_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_new_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_new_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','PP');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_modelFamily_new_stats_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KK');
    case 'PCM_constructModelFamily_FirstFing'
        % first fing separately for trained and untrained 
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 1:4;
        reg = 1:8;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; GG=[]; PP=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe=1:2
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                        end
                        partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions
                        
                        load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                        NatStat = NatStats.G_cent;
                        m = pcm_defineSequenceModel_fixed_splitFirstFing(Seq,sn(p),NatStat); % make
                        [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = D.FoSEx==exe;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    end;
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T.exe = ones(size(T.SN))*exe;
                    
                    P.thetaCr = T.thetaCr;
                    P.G_pred = T.G_pred;
                    P.bayesEst = {T.bayesEst};
                    P.cross_likelihood = {T.cross_likelihood};
                    P.roi = reg(r);
                    P.regType = regType(r);
                    P.regSide = regSide(r);
                    P.sessN = ss;
                    P.exe = exe;
                    T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                    AllReg=addstruct(AllReg,T);
                    PP=addstruct(PP,P);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.sessN = ones(length(K.indx),1)*ss;
                    K.exe = ones(length(K.indx),1)*exe;
                    KK = addstruct(KK,K);
                    idxST = Comb(:,4)==1; % always including SeqType
                    G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,3,5,6]));
                    %idxST = Comb(:,4)==1; 
                    %G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,3,5,6]));
                    G.roi = ones(length(G.indx),1)*reg(r);
                    G.regType = ones(length(G.indx),1)*regType(r);
                    G.regSide = ones(length(G.indx),1)*regSide(r);
                    G.sessN = ones(length(G.indx),1)*ss;
                    G.exe = ones(length(G.indx),1)*exe;
                    GG = addstruct(GG,G);
                    fprintf('Done sess-%d reg-%d/%d - exe-%d\n',ss,r,length(reg),exe);
                end
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing-Trained','FirstFing-Untrained','AllFing','SeqType','Trained','Untrained'};
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_firstFing_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_firstFing_%s_naturalStats-1.mat',parcelType)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_firstFing_%s_naturalStats-1.mat',parcelType)),'-struct','KK');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_seqType_firstFing_%s_naturalStats-1.mat',parcelType)),'-struct','GG');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_firstFing_thetas_%s_naturalStats-1.mat',parcelType)),'-struct','PP');
    case 'PCM_constructModelFamily_Fing'
        % model family with finger effect split
        % here simpler version of model family:
        % includes first finger, all fingers, seqType, sequences (T+U)
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 4;
        reg = [2,3,7];
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; GG=[]; PP=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe=1:2
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                        end
                        partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions
                        
                        load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                        NatStat = NatStats.G_cent;
                        m = pcm_defineFingModels_natStats(Seq,sn(p),NatStat); % make
                        [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = D.FoSEx==exe;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    end;
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T.exe = ones(size(T.SN))*exe;
                    
                    P.thetaCr = T.thetaCr;
                    P.G_pred = T.G_pred;
                    P.bayesEst = {T.bayesEst};
                    P.cross_likelihood = {T.cross_likelihood};
                    P.roi = reg(r);
                    P.regType = regType(r);
                    P.regSide = regSide(r);
                    P.sessN = ss;
                    P.exe = exe;
                    T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                    AllReg=addstruct(AllReg,T);
                    PP=addstruct(PP,P);
                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.sessN = ones(length(K.indx),1)*ss;
                    K.exe = ones(length(K.indx),1)*exe;
                    KK = addstruct(KK,K);
                    G = pcm_calc(T.cross_likelihood,Comb);
                    G.roi = ones(length(G.indx),1)*reg(r);
                    G.regType = ones(length(G.indx),1)*regType(r);
                    G.regSide = ones(length(G.indx),1)*regSide(r);
                    G.sessN = ones(length(G.indx),1)*ss;
                    G.exe = ones(length(G.indx),1)*exe;
                    GG = addstruct(GG,G);
                    fprintf('Done sess-%d reg-%d/%d - exe-%d\n',ss,r,length(reg),exe);
                end
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing'};
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_Fing_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_Fing_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_Fing_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','PP');
    case 'PCM_constructModelFamily_bothRep'
        % here average across both repetitions
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 1:4;
        reg = [1:3,7,8];
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; GG=[]; PP=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for p=1:length(sn)
                    glmDirSubj=fullfile(baseDir,'glmSess',sprintf('glmSess%d',ss),subj_name{sn(p)});
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                    end
                    partVec{p} = D.run;  % runs/partitions
                    if ~naturalStats
                        m = pcm_defineSequenceModels_fixed_noChunks(Seq,sn(p));
                        %m = pcm_defineSequenceModels_test(Seq,sn(p));
                    else
                        load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                        NatStat = NatStats.G_cent;
                        m = pcm_defineSequenceModels_fixed_noChunks_natStats(Seq,sn(p),NatStat); % make
                        %m = pcm_defineSequenceModels_test_natStats(Seq,sn(p),NatStat); % make
                    end
                    [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                    [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                    condVec{p} = repmat(Z,max(D.run),1);
                    Data{p} = beta(:,1:size(D.SN,1))';  % Data is N x P (cond x voxels) - no intercept
                end;
                % fit models
                T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                T.roi = ones(size(T.SN))*reg(r);
                T.regType = ones(size(T.SN))*regType(r);
                T.regSide = ones(size(T.SN))*regSide(r);
                T.sessN = ones(size(T.SN))*ss;
                
                P.thetaCr = T.thetaCr;
                P.G_pred = T.G_pred;
                P.bayesEst = {T.bayesEst};
                P.cross_likelihood = {T.cross_likelihood};
                P.roi = reg(r);
                P.regType = regType(r);
                P.regSide = regSide(r);
                P.sessN = ss;
                T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                AllReg=addstruct(AllReg,T);
                PP=addstruct(PP,P);
                
                % calculations - posterior probability, knockIN/OUT
                K = pcm_calc(T.cross_likelihood,Comb);
                K.roi = ones(length(K.indx),1)*reg(r);
                K.regType = ones(length(K.indx),1)*regType(r);
                K.regSide = ones(length(K.indx),1)*regSide(r);
                K.sessN = ones(length(K.indx),1)*ss;
                KK = addstruct(KK,K);
                idxST = Comb(:,3)==1; % always including SeqType
                G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,4,5]));
                %idxST = Comb(:,4)==1;
                %G = pcm_calc(T.cross_likelihood(:,idxST),Comb(idxST,[1,2,3,5,6]));
                G.roi = ones(length(G.indx),1)*reg(r);
                G.regType = ones(length(G.indx),1)*regType(r);
                G.regSide = ones(length(G.indx),1)*regSide(r);
                G.sessN = ones(length(G.indx),1)*ss;
                GG = addstruct(GG,G);
                fprintf('Done sess-%d reg-%d/%d\n',ss,r,length(reg));
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing','SeqType','Trained','Untrained'};
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_simple_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','PP');
    case 'PCM_constructModelFamily_bothRep_new'
        % here both repetitions with a simpler model
        % first finger, all fingers, sequence
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 1:4;
        reg = [1:3,7,8];
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; GG=[]; PP=[];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for p=1:length(sn)
                    glmDirSubj=fullfile(baseDir,'glmSess',sprintf('glmSess%d',ss),subj_name{sn(p)});
                    D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                    switch (beta_choice)
                        case 'uw'
                            beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                        case 'mw'
                            beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                        case 'raw'
                            beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                    end
                    partVec{p} = D.run;  % runs/partitions
                    load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                    NatStat = NatStats.G_cent;
                    m = pcm_defineSequenceModels_new(Seq,sn(p),NatStat);
                    [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                    [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                    condVec{p} = repmat(Z,max(D.run),1);
                    Data{p} = beta(:,1:size(D.SN,1))';  % Data is N x P (cond x voxels) - no intercept
                end;
                % fit models
                T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                T.roi = ones(size(T.SN))*reg(r);
                T.regType = ones(size(T.SN))*regType(r);
                T.regSide = ones(size(T.SN))*regSide(r);
                T.sessN = ones(size(T.SN))*ss;
                
                P.thetaCr = T.thetaCr;
                P.G_pred = T.G_pred;
                P.bayesEst = {T.bayesEst};
                P.cross_likelihood = {T.cross_likelihood};
                P.roi = reg(r);
                P.regType = regType(r);
                P.regSide = regSide(r);
                P.sessN = ss;
                T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                AllReg=addstruct(AllReg,T);
                PP=addstruct(PP,P);
                
                % calculations - posterior probability, knockIN/OUT
                K = pcm_calc(T.cross_likelihood,Comb);
                K.roi = ones(length(K.indx),1)*reg(r);
                K.regType = ones(length(K.indx),1)*regType(r);
                K.regSide = ones(length(K.indx),1)*regSide(r);
                K.sessN = ones(length(K.indx),1)*ss;
                KK = addstruct(KK,K);
                fprintf('Done sess-%d reg-%d/%d\n',ss,r,length(reg));
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing','Sequence'};
        save(fullfile(pcmDir,sprintf('ModelFamilyComb_bothRep_new_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('ModelFamily_Fit_bothRep_new_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('ModelFamily_Stats_bothRep_new_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','PP');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_modelFamily_bothRep_new_stats_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KK');
    case 'PCM_constructModelFamily_seqType'
        % january 29th 2020
        % separately for trained / untrained sequences
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum, distanceMask
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31];
        sessN = 1:4;
        reg = 1:8;
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; PP = [];
        if strcmp(parcelType,'162tessels') % this to do for new tesselation
            switch regSelect
                case 'subset'
                    reg=sml1_imana_dist('CLUSTER_choose','sessN',sessN)';
                case 'all'
                    reg=sml1_imana_dist('CLUSTER_choose_all','sessN',sessN)';
            end
        end
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for st=1:2 % sequence type
                    for exe=1:2
                        for p=1:length(sn)
                            glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                            switch (beta_choice)
                                case 'uw'
                                    beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                                case 'mw'
                                    beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                                case 'raw'
                                    beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                            end
                            partVec{p} = D.run(D.FoSEx==exe & D.seqType==st);  % runs/partitions
                            if ~naturalStats
                                m = pcm_defineSequenceModels_seqType(Seq,sn(p),st);
                            else
                                load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                                NatStat = NatStats.G_cent;
                                m = pcm_defineSequenceModels_seqType_natStats(Seq,sn(p),NatStat,st); % make
                            end
                            [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                            [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                            condVec{p} = repmat(Z,max(D.run),1);
                            indx = D.FoSEx==exe & D.seqType==st;
                            Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                        end;
                        % fit models
                        T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                        T.roi = ones(size(T.SN))*reg(r);
                        T.regType = ones(size(T.SN))*regType(r);
                        T.regSide = ones(size(T.SN))*regSide(r);
                        T.sessN = ones(size(T.SN))*ss;
                        T.exe = ones(size(T.SN))*exe;
                        T.seqType = ones(size(T.SN))*st;
                        
                        P.thetaCr = T.thetaCr;
                        P.G_pred = T.G_pred;
                        P.bayesEst = {T.bayesEst};
                        P.cross_likelihood = {T.cross_likelihood};
                        P.roi = reg(r);
                        P.regType = regType(r);
                        P.regSide = regSide(r);
                        P.sessN = ss;
                        P.exe = exe;
                        P.seqType = st;
                        T=rmfield(T,{'reg','thetaCr','theta_hat'});
                        AllReg=addstruct(AllReg,T);
                        
                        % calculations - posterior probability, knockIN/OUT
                        K = pcm_calc(T.cross_likelihood,Comb);
                        K.roi = ones(length(K.indx),1)*reg(r);
                        K.regType = ones(length(K.indx),1)*regType(r);
                        K.regSide = ones(length(K.indx),1)*regSide(r);
                        K.sessN = ones(length(K.indx),1)*ss;
                        K.exe = ones(length(K.indx),1)*exe;
                        K.seqType = ones(length(K.indx),1)*st;
                        KK = addstruct(KK,K);
                        PP = addstruct(PP,P);
                        fprintf('Done sess-%d reg-%d seqType-%d/%d - exe-%d\n',ss,r,length(reg),st,exe);
                    end
                end
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
       % modelNames = {'FirstFing','AllFing','Sequence'};
        modelNames = {'FirstFing','AllFing'};
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_seqType_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_seqType_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_seqType_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KK');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_seqType_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','PP');
    case 'PCM_constructFingModel'
        % here only first finger and all finger
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum
        sn = [5:9,11:31];
        sessN = 1:4;
        reg = 1:3;
        naturalStats = 1;
        
        vararginoptions(varargin,{'beta_choice','sn','reg','sessN','algorithm','parcelType','hemi','regSelect','fingType','naturalStats'})
        AllReg=[];
        KK=[]; 
        for ss = sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe=1:2
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        switch (beta_choice)
                            case 'uw'
                                beta = B.betaUW{(B.sn==sn(p) & B.region==reg(r))}';
                            case 'mw'
                                beta = B.betaW{(B.SN==sn(p) & B.region==reg(r))}';
                            case 'raw'
                                beta = B.betaRAW{(B.sn==sn(p) & B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                        end
                       partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions
                       if ~naturalStats
                           m = pcm_defineFingModels(Seq,sn(p));
                       else
                           load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
                           NatStat = NatStats.G_cent;
                           m = pcm_defineFingModels_natStats(Seq,sn(p),NatStat); % make
                       end                        
                        [M, Z] = pcm_buildModelFromFeatures(m,'style','encoding_style','type','component');
                        [Mf,Comb] = pcm_constructModelFamily(M,'fullModel',1);
                        condVec{p} = repmat(Z,max(D.run),1);
                        indx = D.FoSEx==exe;
                        Data{p} = beta(:,indx==1)';  % Data is N x P (cond x voxels) - no intercept
                    end;
                    % fit models
                    T = pcm_fitModels(Data,Mf,partVec,condVec,runEffect,algorithm);
                    T.roi = ones(size(T.SN))*reg(r);
                    T.regType = ones(size(T.SN))*regType(r);
                    T.regSide = ones(size(T.SN))*regSide(r);
                    T.sessN = ones(size(T.SN))*ss;
                    T.exe = ones(size(T.SN))*exe;
                    T=rmfield(T,{'reg','thetaCr','theta_hat'});
                    AllReg=addstruct(AllReg,T);                    
                    % calculations - posterior probability, knockIN/OUT
                    K = pcm_calc(T.cross_likelihood,Comb);
                    K.roi = ones(length(K.indx),1)*reg(r);
                    K.regType = ones(length(K.indx),1)*regType(r);
                    K.regSide = ones(length(K.indx),1)*regSide(r);
                    K.sessN = ones(length(K.indx),1)*ss;
                    K.exe = ones(length(K.indx),1)*exe;
                    KK = addstruct(KK,K);
                    
                    fprintf('Done sess-%d reg-%d/%d - exe-%d\n',ss,r,length(reg),exe);
                end
            end
            fprintf('******************** Done sess-%d ********************\n\n',ss);
        end
        % save variables;
        modelNames = {'FirstFing','AllFing'};
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_finger_%s.mat',parcelType)),'Comb','modelNames');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_finger_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','AllReg');
        save(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_finger_%s_naturalStats-%d.mat',parcelType,naturalStats)),'-struct','KK');
    case 'PCM:create_scaling_model' % construct scaling model
        runEffect = 'random';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels, BG-striatum
        sn = [5:9,11:31];
        sessN = 4;
        reg = 3;
        scaleLim = [0 2];
        nModel = 20;
        
        A(:,:,1) = [[ones(6); zeros(6)] zeros(12,6)];
        A(:,:,2) = [zeros(12,6) [zeros(6); ones(6)]];
        A(:,:,3) = [zeros(12,6) [ones(6); zeros(6)]];
        A(:,:,4) = [[eye(6); zeros(6)] zeros(12,6)]; % common exe1
        A(:,:,5) = [zeros(12,6) [zeros(6); eye(6)]];

        
        scaleL = linspace(scaleLim(1),scaleLim(2),nModel); % correlation models to assess

        for r=1:length(scaleL);
            M{r} = pcm_buildCorrModel('type','nonlinear','withinCov','iid','numCond',2,'numItems',6,'r',scaleL(r),'condEffect',1);
        end
        for ss=sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r=reg
                for st=1:2 % trained, untrained
                    for p=1:length(sn)
                        glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn(p)});
                        D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                        beta = B.betaW{(B.SN==sn(p) & B.region==r)}';
                        partVec{p} = D.run(D.seqType==st);  % runs/partitions
                        condVec{p} = D.seqNumb(D.seqType==st);  % runs/partitions
                        Data{p} = beta(:,D.seqType==st)';
                    end
                    keyboard;
                end
            end
        end
    case 'PCM_compareFing'
        % compare PCM finger with .. without natStats
        reg = 1:3;
        sessN = 4;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sessN','reg'});
        TT = [];
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_finger_%s_naturalStats-%d.mat',parcelType,1)));
        T.natStat = ones(size(T.SN));
        T1 = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_finger_%s_naturalStats-%d.mat',parcelType,0)));
        T1.natStat = zeros(size(T.SN));
        TT = addstruct(TT,T);
        TT = addstruct(TT,T1);
        for ss=sessN
            for r=reg
                t = getrow(TT,TT.roi==r & TT.sessN==ss);
                figure
                plt.bar(t.natStat,t.bayesEst,'leg',{'null','firstF','allF','first+all'});
                xlabel('naturalStats'); ylabel('log-BF')
                title(regname{r});
            end
        end

    case 'PCM_variance_decomposition' % used in the paper
        % here decompose variance - from thetas estimated
        parcelType = 'Brodmann'; % Brodmann or distanceMask
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_simple_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,sessN)));
             
        for r=roi
            for e=1:2 % execution
                t = getrow(T,T.roi==r & T.sessN==sessN & T.exe==e);
                d = getrow(D,D.region==r & D.FoSEx==e);
                logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
                prob = exp(logLike);
                prob = bsxfun(@rdivide,prob,sum(prob,2));
                for s=1:size(d.SN,1)
                    G = rsa_squareIPMfull(d.IPMfull(s,:));
                    H=eye(12)-ones(12)./12;  % centering matrix!
                    G = H*G*H';  % double centered G matrix - rows and columns
                    Gcv_trace(s) = trace(G);
                end
                for c=1:size(Comb,2)
                    mIdx = Comb(:,c)==1;
                    %mIdx = [Comb(:,c)==1 & Comb(:,3)==1];
                    thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                    thIdx2 = thIdx(thIdx>0);
                    thetas = t.thetaCr(mIdx);
                    Gpred = t.G_pred(mIdx);
                    LL = logLike(:,mIdx);
                    PP = prob(:,mIdx);
                    ev = 0;
                    for m=1:size(LL,2)
                        ev = ev + mean(PP(:,m)) * mean(exp(thetas{m}(:,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace';
                    end
                    K.roi = r;
                    K.exe = e;
                    K.comb = c;
                    K.var = sum(ev);
                    KK = addstruct(KK,K);
                    M.var = ev;
                    M.exe = ones(size(M.var))*e;
                    M.comb = ones(size(M.var))*c;
                    M.roi = ones(size(M.var))*r;
                    M.Gcv_trace = Gcv_trace';
                    M.sn = [1:26]';
                    MM = addstruct(MM,M);
                end
            end
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            plt.bar(MM.comb,MM.var,'split',MM.exe,'subset',MM.roi==r,'leg',{'exe1','exe2'},'style',stySess);
            ylabel('variance'); title(regname{r});
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            G = pivottable(MM.comb,MM.exe,MM.var,'robustmean','subset',MM.roi==r);
            G(:,1) = G(:,1)./sum(G(:,1));
            G(:,2) = G(:,2)./sum(G(:,2));
            mosaic_plot(G); title(sprintf(regname_cortex{r}));
        end
        keyboard;
        % here quantify stats: % reduction with repetition for each feature across regions
        RR = [];
        for r=unique(roi)
            for i=unique(MM.comb)'
                t=getrow(MM,MM.roi==r&MM.comb==i);
                R.percent = 1-[t.var(t.exe==2)./t.var(t.exe==1)];
                R.sn = t.sn(t.exe==1);
                R.comb = t.comb(t.exe==1);
                R.roi = t.roi(t.exe==1);
                RR = addstruct(RR,R);
            end
        end
        % amount of reduction
        % first finger: M1,PMd,SPLa
        mean(RR.percent(RR.comb==1&RR.roi==2));
        mean(RR.percent(RR.comb==1&RR.roi==3));
        mean(RR.percent(RR.comb==1&RR.roi==7));
        % sequence type
        mean(RR.percent(RR.comb==3&RR.roi==3));
        mean(RR.percent(RR.comb==3&RR.roi==7));
        % trained sequence
        mean(RR.percent(RR.comb==4&RR.roi==3));
        mean(RR.percent(RR.comb==4&RR.roi==7));
        % compare first-all finger in M1
        ttestDirect(RR.percent,[RR.comb,RR.sn],2,'paired','subset',RR.roi==2&ismember(RR.comb,[1,2]));
        % compare first finger - trained sequence in PMd/SPLa
        ttestDirect(RR.percent,[RR.comb,RR.sn],2,'paired','subset',RR.roi==3&ismember(RR.comb,[1,4]));
        ttestDirect(RR.percent,[RR.comb,RR.sn],2,'paired','subset',RR.roi==7&ismember(RR.comb,[1,4]));
    case 'PCM_variance_decomposition_noError'
        % here PCM without any error
        parcelType = 'Brodmann'; % Brodmann or distanceMask
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_noError_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_noError_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_ErrFoSEx_%s_multiPW_sess%d.mat',parcelType,sessN)));
             
        for r=roi
            for e=1:2 % execution
                t = getrow(T,T.roi==r & T.sessN==sessN & T.exe==e);
                d = getrow(D,D.region==r & D.FoSEx==e);
                logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
                prob = exp(logLike);
                prob = bsxfun(@rdivide,prob,sum(prob,2));
                for s=1:size(d.SN,1)
                    G = rsa_squareIPMfull(d.IPMfull(s,:));
                    H=eye(12)-ones(12)./12;  % centering matrix!
                    G = H*G*H';  % double centered G matrix - rows and columns
                    Gcv_trace(s) = trace(G);
                end
                for c=1:size(Comb,2)
                    mIdx = Comb(:,c)==1;
                    %mIdx = [Comb(:,c)==1 & Comb(:,3)==1];
                    thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                    thIdx2 = thIdx(thIdx>0);
                    thetas = t.thetaCr(mIdx);
                    Gpred = t.G_pred(mIdx);
                    LL = logLike(:,mIdx);
                    PP = prob(:,mIdx);
                    ev = 0;
                    for m=1:size(LL,2)
                        ev = ev + mean(PP(:,m)) * mean(exp(thetas{m}(:,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace';
                    end
                    K.roi = r;
                    K.exe = e;
                    K.comb = c;
                    K.var = sum(ev);
                    KK = addstruct(KK,K);
                    M.var = ev;
                    M.exe = ones(size(M.var))*e;
                    M.comb = ones(size(M.var))*c;
                    M.roi = ones(size(M.var))*r;
                    M.Gcv_trace = Gcv_trace';
                    M.sn = [1:26]';
                    MM = addstruct(MM,M);
                end
            end
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            plt.bar(MM.comb,MM.var,'split',MM.exe,'subset',MM.roi==r,'leg',{'exe1','exe2'},'style',stySess);
            ylabel('variance'); title(regname{r});
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            G = pivottable(MM.comb,MM.exe,MM.var,'robustmean','subset',MM.roi==r);
            G(:,1) = G(:,1)./sum(G(:,1));
            G(:,2) = G(:,2)./sum(G(:,2));
            mosaic_plot(G); title(sprintf(regname_cortex{r}));
        end
        keyboard;
        % here quantify stats: first fing vs seqType compression: PMd, SPLa
        
        figure % here consider 'noise ceilings'
        T = tapply(MM,{'exe','roi','sn'},{'var','sum'},{'Gcv_trace','mean'});
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            plt.bar(T.exe,T.var,'subset',T.roi==r);
            hold on;
            drawline(mean(T.Gcv_trace(T.roi==r&T.exe==1)),'dir','horz');
            drawline(mean(T.Gcv_trace(T.roi==r&T.exe==2)),'dir','horz');
        end
    case 'PCM_variance_decomposition_new' % not in use
        parcelType = 'Brodmann'; % Brodmann or distanceMask
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_new_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_new_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,sessN)));
             
        for r=roi
            for e=1:2 % execution
                t = getrow(T,T.roi==r & T.sessN==sessN & T.exe==e);
                d = getrow(D,D.region==r & D.FoSEx==e);
                logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
                prob = exp(logLike);
                prob = bsxfun(@rdivide,prob,sum(prob,2));
                for s=1:size(d.SN,1)
                    G = rsa_squareIPMfull(d.IPMfull(s,:));
                    H=eye(12)-ones(12)./12;  % centering matrix!
                    G = H*G*H';  % double centered G matrix - rows and columns
                    Gcv_trace(s) = trace(G);
                end
                for c=1:size(Comb,2)
                    mIdx = Comb(:,c)==1;
                    %mIdx = [Comb(:,c)==1 & Comb(:,3)==1];
                    thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                    thIdx2 = thIdx(thIdx>0);
                    thetas = t.thetaCr(mIdx);
                    Gpred = t.G_pred(mIdx);
                    LL = logLike(:,mIdx);
                    PP = prob(:,mIdx);
                    ev = 0;
                    for m=1:size(LL,2)
                        ev = ev + mean(PP(:,m)) * mean(exp(thetas{m}(:,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace';
                    end
                    K.roi = r;
                    K.exe = e;
                    K.comb = c;
                    K.var = sum(ev);
                    KK = addstruct(KK,K);
                    M.var = ev;
                    M.exe = ones(size(M.var))*e;
                    M.comb = ones(size(M.var))*c;
                    M.roi = ones(size(M.var))*r;
                    MM = addstruct(MM,M);
                end
            end
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            plt.bar(MM.comb,MM.var,'split',MM.exe,'subset',MM.roi==r,'leg',{'exe1','exe2'},'style',stySess);
            ylabel('variance'); title(regname{r});
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            G = pivottable(MM.comb,MM.exe,MM.var,'robustmean','subset',MM.roi==r);
            G(:,1) = G(:,1)./sum(G(:,1));
            G(:,2) = G(:,2)./sum(G(:,2));
            mosaic_plot(G); title(sprintf(regname_cortex{r}));
        end
    case 'PCM_variance_decomposition_firstFing'
        % here split first finger - trained / untrained
         parcelType = 'Brodmann';
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_firstFing_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_firstFing_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,sessN)));
             
        for r=roi
            for e=1:2 % execution
                t = getrow(T,T.roi==r & T.sessN==sessN & T.exe==e);
                d = getrow(D,D.region==r & D.FoSEx==e);
                logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
                prob = exp(logLike);
                prob = bsxfun(@rdivide,prob,sum(prob,2));
                for s=1:size(d.SN,1)
                    G = rsa_squareIPMfull(d.IPMfull(s,:));
                    H=eye(12)-ones(12)./12;  % centering matrix!
                    G = H*G*H';  % double centered G matrix - rows and columns
                    Gcv_trace(s) = trace(G);
                end
                for c=1:size(Comb,2)
                    mIdx = Comb(:,c)==1;
                    %mIdx = [Comb(:,c)==1 & Comb(:,3)==1];
                    thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                    thIdx2 = thIdx(thIdx>0);
                    thetas = t.thetaCr(mIdx);
                    Gpred = t.G_pred(mIdx);
                    LL = logLike(:,mIdx);
                    PP = prob(:,mIdx);
                    ev = 0;
                    for m=1:size(LL,2)
                        ev = ev + mean(PP(:,m)) * mean(exp(thetas{m}(:,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace';
                    end
                    K.roi = r;
                    K.exe = e;
                    K.comb = c;
                    K.var = sum(ev);
                    KK = addstruct(KK,K);
                    M.var = ev;
                    M.exe = ones(size(M.var))*e;
                    M.comb = ones(size(M.var))*c;
                    M.roi = ones(size(M.var))*r;
                    MM = addstruct(MM,M);
                end
            end
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            plt.bar(MM.comb,MM.var,'split',MM.exe,'subset',MM.roi==r,'leg',{'exe1','exe2'},'style',stySess);
            ylabel('variance'); title(regname{r});
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            G = pivottable(MM.comb,MM.exe,MM.var,'robustmean','subset',MM.roi==r);
            G(:,1) = G(:,1)./sum(G(:,1));
            G(:,2) = G(:,2)./sum(G(:,2));
            mosaic_plot(G); title(sprintf(regname_cortex{r}));
        end
    case 'PCM_variance_decomposition_finger'
        parcelType = 'Brodmann';
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_Fing_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_Fing_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,sessN)));
             
        for r=roi
            for e=1:2 % execution
                t = getrow(T,T.roi==r & T.sessN==sessN & T.exe==e);
                d = getrow(D,D.region==r & D.FoSEx==e);
                logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
                prob = exp(logLike);
                prob = bsxfun(@rdivide,prob,sum(prob,2));
                for s=1:size(d.SN,1)
                    G = rsa_squareIPMfull(d.IPMfull(s,:));
                    H=eye(12)-ones(12)./12;  % centering matrix!
                    G = H*G*H';  % double centered G matrix - rows and columns
                    Gcv_trace(s) = trace(G);
                end
                for c=1:size(Comb,2)
                    mIdx = Comb(:,c)==1;
                    %mIdx = [Comb(:,c)==1 & Comb(:,3)==1];
                    thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                    thIdx2 = thIdx(thIdx>0);
                    thetas = t.thetaCr(mIdx);
                    Gpred = t.G_pred(mIdx);
                    LL = logLike(:,mIdx);
                    PP = prob(:,mIdx);
                    ev = 0;
                    for m=1:size(LL,2)
                        ev = ev + mean(PP(:,m)) * mean(exp(thetas{m}(:,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace';
                    end
                    K.roi = r;
                    K.exe = e;
                    K.comb = c;
                    K.var = sum(ev);
                    KK = addstruct(KK,K);
                    M.var = ev;
                    M.exe = ones(size(M.var))*e;
                    M.comb = ones(size(M.var))*c;
                    M.roi = ones(size(M.var))*r;
                    M.sn = d.SN;
                    MM = addstruct(MM,M);
                end
            end
        end
        figure
        for r=roi
            subplot(1,numel(roi),find(roi==r));
            plt.bar(MM.comb,MM.var,'split',MM.exe,'subset',MM.roi==r,'leg',{'exe1','exe2'},'style',stySess);
            ylabel('variance'); title(regname{r});
        end
        
        figure
        for r=roi
            subplot(1,numel(roi),find(roi==r))
            G = pivottable(MM.comb,MM.exe,MM.var,'robustmean','subset',MM.roi==r);
            G(:,1) = G(:,1)./sum(G(:,1));
            G(:,2) = G(:,2)./sum(G(:,2));
            mosaic_plot(G); title(sprintf(regname_cortex{r}));
        end
    case 'PCM_variance_decomposition_seqType'
        parcelType = 'Brodmann';
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_seqType_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_seqType_%s.mat',parcelType)));
        %T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_finger_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        %load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_finger_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,sessN)));
        
        for r=roi
            for st = 1:2 % seqType
                for e=1:2 % execution
                    t = getrow(T,T.roi==r & T.sessN==sessN & T.exe==e & T.seqType==st);
                    d = getrow(D,D.region==r & D.FoSEx==e);
                    logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
                    prob = exp(logLike);
                    prob = bsxfun(@rdivide,prob,sum(prob,2));
                    for s=1:size(d.SN,1)
                        G = rsa_squareIPMfull(d.IPMfull(s,:));
                        if st==1
                            G = G(1:6,1:6);
                        else
                            G = G(7:12,7:12);
                        end
                        H=eye(6)-ones(6)./6;  % centering matrix!
                        G = H*G*H';  % double centered G matrix - rows and columns
                        Gcv_trace(s) = trace(G);
                    end
                    for c=1:size(Comb,2)
                        mIdx = Comb(:,c)==1;
                        thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                        thIdx2 = thIdx(thIdx>0);
                        thetas = t.thetaCr(mIdx);
                        Gpred = t.G_pred(mIdx);
                        LL = logLike(:,mIdx);
                        PP = prob(:,mIdx);
                        ev = 0;
                        sub_idx = Gcv_trace>0; % index subjects - only those with some evidence
                        for m=1:size(LL,2)
                            ev = ev + mean(PP(sub_idx,m)) * mean(exp(thetas{m}(sub_idx,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace(sub_idx)';
                        end
                        K.roi = r;
                        K.exe = e;
                        K.comb = c;
                        K.var = sum(ev);
                        KK = addstruct(KK,K);
                        M.var = ev;
                        M.exe = ones(size(M.var))*e;
                        M.comb = ones(size(M.var))*c;
                        M.roi = ones(size(M.var))*r;
                        M.seqType = ones(size(M.var))*st;
                        M.Gcv = Gcv_trace(sub_idx)';
                        MM = addstruct(MM,M);
                    end
                end
            end
        end
        
        for r=roi
            figure
            subplot(121)
            plt.bar(MM.comb,MM.var,'split',MM.exe,'subset',MM.roi==r & MM.seqType==1,'leg',{'exe1','exe2'},'style',stySess);
            ylabel('variance'); title(sprintf('trained - %s',regname{r}));
            subplot(122)
            plt.bar(MM.comb,MM.var,'split',MM.exe,'subset',MM.roi==r & MM.seqType==2,'leg',{'exe1','exe2'},'style',stySess);
            ylabel('variance'); title(sprintf('untrained - %s',regname{r}));
            figure
            plt.bar(MM.comb,MM.var,'split',[MM.seqType,MM.exe],'subset',MM.roi==r,'leg',{'exe1','exe2'},'style',styRSbeh);
            ylabel('variance'); title(regname{r});
        end
        
        for r=roi
            figure
            for st=1:2
                subplot(1,2,st)
                G = pivottable(MM.comb,MM.exe,MM.var,'robustmean','subset',MM.roi==r & MM.seqType==st);
                for i=1:size(G,2)
                    G(:,i) = G(:,i)./sum(G(:,i));
                end
                mosaic_plot(G); title(sprintf(regname_cortex{r}));
            end
        end
    case 'PCM_variance_decomposition_bothRep'
        % both repetitions
        parcelType = 'Brodmann';
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('ModelFamily_Stats_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('ModelFamilyComb_simple_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_%s_multiPW_sess%d.mat',parcelType,sessN)));
             
        for r=roi
            t = getrow(T,T.roi==r & T.sessN==sessN);
            d = getrow(D,D.region==r);
            logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
            prob = exp(logLike);
            prob = bsxfun(@rdivide,prob,sum(prob,2));
            for s=1:size(d.SN,1)
                G = rsa_squareIPMfull(d.IPMfull(s,:));
                H=eye(12)-ones(12)./12;  % centering matrix!
                G = H*G*H';  % double centered G matrix - rows and columns
                Gcv_trace(s) = trace(G);
            end
            for c=1:size(Comb,2)
                mIdx = Comb(:,c)==1;
                %mIdx = [Comb(:,c)==1 & Comb(:,3)==1];
                thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                thIdx2 = thIdx(thIdx>0);
                thetas = t.thetaCr(mIdx);
                Gpred = t.G_pred(mIdx);
                LL = logLike(:,mIdx);
                PP = prob(:,mIdx);
                ev = 0;
                for m=1:size(LL,2)
                    ev = ev + mean(PP(:,m)) * mean(exp(thetas{m}(:,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace';
                end
                K.roi = r;
                K.comb = c;
                K.var = sum(ev);
                KK = addstruct(KK,K);
                M.var = ev;
                M.comb = ones(size(M.var))*c;
                M.roi = ones(size(M.var))*r;
                MM = addstruct(MM,M);
            end
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            plt.bar(MM.comb,MM.var,'subset',MM.roi==r,'style',stySess);
            ylabel('variance'); title(regname{r});
        end
        figure
        G = pivottable(MM.comb,MM.roi,MM.var,'robustmean');
        G(:,1) = G(:,1)./sum(G(:,1));
        G(:,2) = G(:,2)./sum(G(:,2));
        G(:,3) = G(:,3)./sum(G(:,3));
        mosaic_plot(G); title(sprintf(regname_cortex{r}));
    case 'PCM_variance_decomposition_bothRep_new'
        parcelType = 'Brodmann';
        roi=[2,3,7];
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        KK=[]; MM = [];
        T = load(fullfile(pcmDir,sprintf('ModelFamily_Stats_bothRep_new_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('ModelFamilyComb_bothRep_new_%s.mat',parcelType)));
        D = load(fullfile(betaDir,'group',sprintf('stats_%s_multiPW_sess%d.mat',parcelType,sessN)));
             
        for r=roi
            t = getrow(T,T.roi==r & T.sessN==sessN);
            d = getrow(D,D.region==r);
            logLike = bsxfun(@minus,t.cross_likelihood{1},max(t.cross_likelihood{1},[],2))+50;
            prob = exp(logLike);
            prob = bsxfun(@rdivide,prob,sum(prob,2));
            for s=1:size(d.SN,1)
                G = rsa_squareIPMfull(d.IPMfull(s,:));
                H=eye(12)-ones(12)./12;  % centering matrix!
                G = H*G*H';  % double centered G matrix - rows and columns
                Gcv_trace(s) = trace(G);
            end
            for c=1:size(Comb,2)
                mIdx = Comb(:,c)==1;
                %mIdx = [Comb(:,c)==1 & Comb(:,3)==1];
                thIdx = (sum(Comb(:,1:c-1),2)+1).*mIdx; % which theta to take into account per model
                thIdx2 = thIdx(thIdx>0);
                thetas = t.thetaCr(mIdx);
                Gpred = t.G_pred(mIdx);
                LL = logLike(:,mIdx);
                PP = prob(:,mIdx);
                ev = 0;
                for m=1:size(LL,2)
                    ev = ev + mean(PP(:,m)) * mean(exp(thetas{m}(:,thIdx2(m)))*trace(Gpred{m}))./Gcv_trace';
                end
                K.roi = r;
                K.comb = c;
                K.var = sum(ev);
                KK = addstruct(KK,K);
                M.var = ev;
                M.comb = ones(size(M.var))*c;
                M.roi = ones(size(M.var))*r;
                MM = addstruct(MM,M);
            end
        end
        figure
        for r=unique(roi)
            subplot(1,numel(roi),find(r==roi))
            plt.bar(MM.comb,MM.var,'subset',MM.roi==r,'style',stySess);
            ylabel('variance'); title(regname{r});
        end
        figure
        G = pivottable(MM.comb,MM.roi,MM.var,'robustmean');
        G(:,1) = G(:,1)./sum(G(:,1));
        G(:,2) = G(:,2)./sum(G(:,2));
        G(:,3) = G(:,3)./sum(G(:,3));
        mosaic_plot(G); title(sprintf(regname_cortex{r}));
        
    case 'PCM_plotG_pred_est'
        % here plot the G - estimated and predicted
        sessN=4;
        reg=[2,7];
        parcelType = 'Brodmann';
        betaChoice = 'multiPW';
        naturalStats = 1;
        numGroup = 2;
        vararginoptions(varargin,{'sn','sessN','reg','numGroup'});
        
        D = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        H = eye(12) - 1/12;
        
        load(fullfile(pcmDir,'naturalstatisticmodel.mat'));
        NatStat = NatStats.G_cent;

        idx=1;
        for g=numGroup % group
            m = pcm_defineSequenceModels_fixed_noChunks_natStats(Seq,g,NatStat); % make
            figure(1)
            for i=1:size(m,2)
                [M, Z] = pcm_buildModelFromFeatures(m(i),'style','encoding_style','type','component');
                subplot(1,size(m,2),i)
                imagesc(Z*Z');
                hold on; drawline(6.5,'dir','horz'); drawline(6.5,'dir','vert');
            end
            [M, Z] = pcm_buildModelFromFeatures(m(1),'style','encoding_style','type','component');
            M = Z*Z';
            for ss=sessN
                T = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%s_sess%d',parcelType,betaChoice,ss)));
                for e=1:2
                    for r=reg
                        if g==1
                            t = getrow(T,T.region==r & T.FoSEx==e & rem(T.SN,2)==g);
                        else
                            t = getrow(T,T.region==r & T.FoSEx==e & rem(T.SN,2)==0);
                        end
                        %t = getrow(T,T.region==r & T.FoSEx==e & T.SN==s);
                        %RDM = rsa_squareRDM(nanmean(t.RDM,1));
                        RDM = ssqrt(rsa_squareRDM(nanmean(t.RDM,1)));
                        G = -0.5*H*RDM*H';
                        if e==1&&r==2
                            G = G+0.18*M; %0.003 for non sqrt, 0.1 for sqrt
                            %G = G+0.08*M; %0.003 for non sqrt, 0.1 for sqrt
                        end
                        figure(3)
                        subplot(2,2,idx)
                        imagesc(G);
                        title(sprintf('%s, group %d, ssqrt, exe-%d',regname_cortex{r},g,e));
                        hold on; drawline(6.5,'dir','horz'); drawline(6.5,'dir','vert');
                        subplot(2,2,idx)
                        RDM = rsa_squareRDM(nanmean(t.RDM,1));
                        G = -0.5*H*RDM*H';
                        if e==1&&r==2
                            %G = G+0.01*M; %0.003 for non sqrt, 0.1 for sqrt
                            G = G+0.01*M; %0.003 for non sqrt, 0.1 for sqrt
                        end
                        %imagesc(G);
                        title(sprintf('%s, group %d, exe-%d',regname_cortex{r},g,e));
                        %hold on; drawline(6.5,'dir','horz'); drawline(6.5,'dir','vert');
                        idx=idx+1;
                    end
                end
            end; % session
        end
    case 'PCM_modelFamily_fing_seqType'
        parcelType = 'Brodmann';
        roi=[2,3,7];
        var='logBayes'; % knockOUT or logBayes;
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'roi','var','sessN'});
        
        T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_seqType_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        for r=roi
            for ss=sessN
                figure
                subplot(121)
                plt.bar(T.indx,T.(var),'split',T.exe,'subset',T.roi==r & T.sessN==ss & T.seqType==1,'leg',{'exe1','exe2'},'leglocation','northeast');
                ylabel(var);
                title(sprintf('%s - sess%d trained',regname{r},ss));
                subplot(122)
                plt.bar(T.indx,T.(var),'split',T.exe,'subset',T.roi==r & T.sessN==ss & T.seqType==2,'leg',{'exe1','exe2'},'leglocation','northeast');
                xlabel('modelType');
                ylabel(var);
                title(sprintf('%s - sess%d untrained',regname{r},ss));
            end
        end
    case 'PCM_modelFamily_firstFing_split'
        % splitting first finger into trained and untrained
        parcelType = 'Brodmann';
        roi=[2,3,7];
        var='logBayes'; % knockOUT or logBayes;
        sessN = 4;
        naturalStats = 1;
        seqType = 1; % whether only considering models with seqType or not
        vararginoptions(varargin,{'roi','var','sessN','seqType'});
        if ~seqType
            T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_firstFing_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        else
            T = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_seqType_firstFing_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        end
        for r=roi
            for ss=sessN
                figure
                plt.bar(T.indx,T.(var),'split',T.exe,'subset',T.roi==r & T.sessN==ss,'leg',{'exe1','exe2'},'leglocation','northeast');
                ylabel(var);
                title(sprintf('%s - sess%d',regname{r},ss));
                xlabel('modelType');
            end
        end
    case 'PCM_plotModelFamily'  % this needed for stats
        parcelType = 'Brodmann';
        roi=[2,3];
        var='logBayes'; % knockOUT or logBayes;
        type = 'seqType';
        sessN = 1:4;
        naturalStats = 0;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        KK = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        switch type
            case 'seqType'
                TT = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_modelFamily_stats_%s_naturalStats-%d.mat',parcelType,naturalStats)));        
                %TT = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_test_%s_naturalStats-%d.mat',parcelType,naturalStats)));
                % features: 1) first finger, 2) all fingers, 3) seqType, 4) trained Seq, 5) untrained seq
            case 'noSeqType'
                TT = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_simple_noSeqType_%s_naturalStats-%d.mat',parcelType,naturalStats)));
                %TT = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_test_noSeqType_%s_naturalStats-%d.mat',parcelType,naturalStats)));
                % features: 1) first finger, 2) all fingers, 3) trained Seq, 4) untrained seq
        end
        for r=roi
            for ss=sessN
                figure
                subplot(2,2,1)
                plt.bar(TT.indx,TT.(var),'split',TT.exe,'subset',TT.roi==r&TT.sessN==ss,'leg',{'exe1','exe2'},'leglocation','northeast');
                hold on;
                drawline(0,'dir','horz');
                xlabel('modelType');
                ylabel(var);
                %plt.match('y');
                title(sprintf('%s - sess%d',regname{r},ss));
                % here add stats - bayesfactors against 0
                fprintf('%s, exe-1:\n',regname_cortex{r});
                ttestDirect(TT.(var),[TT.sn],2,'onesample','subset',TT.roi==r&TT.sessN==ss&TT.exe==1,'split',TT.indx);
                fprintf('%s, exe-2:\n',regname_cortex{r});
                ttestDirect(TT.(var),[TT.sn],2,'onesample','subset',TT.roi==r&TT.sessN==ss&TT.exe==2,'split',TT.indx);
                subplot(2,2,3)
                k = getrow(KK,KK.roi==r&KK.sessN==ss);
                [~,j] = max(k.cross_likelihood);
                plt.bar(k.exe,k.bayesEst(:,mode(j)));
                ylabel('log-BF'); title('Best fitting model'); xlabel('execution');

            end
        end
    case 'PCM_plotModelFamily_seqType'
        parcelType = 'Brodmann';
        roi=[2,3];
        var='logBayes'; % knockOUT or logBayes;
        sessN = 4;
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','roi','var','type','naturalStats','sessN'});
        
        TT = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_seqType_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        % features: 1) first finger, 2) all fingers, 3) sequence
        for r=roi
            for ss=sessN
                figure
                subplot(121)
                plt.bar(TT.indx,TT.(var),'split',TT.seqType,'subset',TT.roi==r&TT.sessN==ss&TT.exe==1,'leg',{'trained','untrained'},'leglocation','northeast');
                hold on;
                drawline(0,'dir','horz');
                xlabel('modelType');
                ylabel(var);
                plt.match('y'); title(sprintf('%s - sess%d, exe1',regname{r},ss));
                subplot(122)
                plt.bar(TT.indx,TT.(var),'split',TT.seqType,'subset',TT.roi==r&TT.sessN==ss&TT.exe==2,'leg',{'trained','untrained'},'leglocation','northeast');
                hold on;
                drawline(0,'dir','horz');
                xlabel('modelType');
                ylabel(var); title(sprintf('%s - sess%d, exe2',regname{r},ss));
            end
        end
    case 'PCM_noiseCeilings_repsup' % this needed in final figures
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
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe = 1:2 % sequence type
                    for g = 1:2 % two groups of subjects
                        for p=1:length(sn_group{g})
                            glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn_group{g}(p)});
                            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                            switch (beta_choice)
                                case 'uw'
                                    beta = B.betaUW{(B.sn==sn_group{g}(p)&B.region==reg(r))}';
                                case 'mw'
                                    beta = B.betaW{(B.SN==sn_group{g}(p)&B.region==reg(r))}';
                                case 'raw'
                                    beta = B.betaRAW{(B.sn==sn_group{g}(p)&B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                            end
                            partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions;
                            condVec{p} = D.seqNumb(D.FoSEx==exe);
                            Data{p} = beta(:,(D.FoSEx==exe))';  % Data is N x P (cond x voxels) - no intercept
                        end;
                        % fit models
                        T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                        T.roi = ones(size(T.SN))*reg(r);
                        T.regType = ones(size(T.SN))*regType(r);
                        T.regSide = ones(size(T.SN))*regSide(r);
                        T.sessN = ones(size(T.SN))*ss;
                        T.exe = ones(size(T.SN))*exe;
                        T=rmfield(T,{'reg','thetaCr','theta_hat'});
                        if g==2
                            T.SN = T.SN + length(sn_group{1});
                        end
                        AllReg=addstruct(AllReg,T);
                        clear partVec condVec Data;
                    end % group
                end % repetition
            end
            fprintf('Done sess-%d.\n\n\n',ss);
        end
        % save variables;
        save(fullfile(pcmDir,sprintf('NoiseCeilings_repsup_%s.mat',parcelType)),'-struct','AllReg');   
    case 'PCM_noiseCeilings_repsup_noError'
        % construct noise ceilings - from glm without errors
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm = 'NR'; % minimize or NR
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        regSelect = 'all'; % all or subset
        sn = [5:9,11:31]; 
        sessN = 4;
        reg = [1:3,7,8];
        
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
            B=load(fullfile(betaDir,'group',sprintf('betas_ErrFoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe = 1:2 % sequence type
                    for g = 1:2 % two groups of subjects
                        for p=1:length(sn_group{g})
                            glmDirSubj=fullfile(glmFoSExErrorDir{ss}, subj_name{sn_group{g}(p)});
                            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                            D=getrow(D,D.isError==0);
                            switch (beta_choice)
                                case 'uw'
                                    beta = B.betaUW{(B.sn==sn_group{g}(p)&B.region==reg(r))}';
                                case 'mw'
                                    beta = B.betaW{(B.SN==sn_group{g}(p)&B.region==reg(r))}';
                                case 'raw'
                                    beta = B.betaRAW{(B.sn==sn_group{g}(p)&B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                            end
                            partVec{p} = D.run(D.FoSEx==exe);  % runs/partitions;
                            condVec{p} = D.seqNumb(D.FoSEx==exe);
                            Data{p} = beta(:,(D.FoSEx==exe))';  % Data is N x P (cond x voxels) - no intercept
                        end;
                        % fit models
                        T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                        T.roi = ones(size(T.SN))*reg(r);
                        T.regType = ones(size(T.SN))*regType(r);
                        T.regSide = ones(size(T.SN))*regSide(r);
                        T.sessN = ones(size(T.SN))*ss;
                        T.exe = ones(size(T.SN))*exe;
                        T=rmfield(T,{'reg','thetaCr','theta_hat'});
                        if g==2
                            T.SN = T.SN + length(sn_group{1});
                        end
                        T = rmfield(T,'G_pred');
                        AllReg=addstruct(AllReg,T);
                        clear partVec condVec Data;
                    end % group
                end % repetition
            end
            fprintf('Done sess-%d.\n\n\n',ss);
        end
        % save variables;
        save(fullfile(pcmDir,sprintf('NoiseCeilings_repsup_noError_%s.mat',parcelType)),'-struct','AllReg');   
    case 'PCM_noiseCeilings_plot'
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        sessN = 1:4;
        reg = [2:3,7];
        naturalStats = 0;
        vararginoptions(varargin,{'parcelType','reg','naturalStats','sessN'});
        % plot PCM noise ceilings
        T = load(fullfile(pcmDir,sprintf('NoiseCeilings_repsup_%s.mat',parcelType)));
        %K = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_modelFamily_stats_%s_naturalStats-%d',parcelType,naturalStats)));
        Th = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        K = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_simple_%s.mat',parcelType)));
        idx=1;
        figure
        for r=reg
            for ss=sessN
                t=getrow(T,T.roi==r&T.sessN==ss);
                k=getrow(K,K.roi==r&K.sessN==ss);
                th=getrow(Th,Th.roi==r&Th.sessN==ss);
                [~,bestM1] = max(mean(k.bayesEst(k.exe==1,:)));
                [~,bestM2] = max(mean(k.bayesEst(k.exe==2,:)));
             %   bestM = mode(j); % determine the best model
                % form a new structure
                data = k.cross_likelihood(:,[1:6,bestM1,bestM2]);
                data = bsxfun(@minus,data,data(:,1));
                M.data = data(:);
                M.indx = kron([1:8]',ones(length(t.SN),1));
                M.sessN = repmat(k.sessN,8,1);
                M.exe = repmat(k.exe,8,1);
                M.SN = repmat(k.SN,8,1);
                M = normData(M,'data');
                subplot(1,numel(reg),idx)
                plt.bar(M.indx,M.normdata,'subset',M.indx~=1,'split',M.exe,'leg',{'exe1','exe2'},'leglocation','northwest','style',stySess); hold on; 
                drawline(mean(t.cross_likelihood(t.exe==1,2)-t.cross_likelihood(t.exe==1,1)),'dir','horz');  % exe1
                drawline(mean(t.cross_likelihood(t.exe==2,2)-t.cross_likelihood(t.exe==2,1)),'dir','horz','linestyle','--');  % exe2
                drawline(0,'dir','horz'); 
                set(gca,'XTickLabel',{'FF','','AF','','ST','','TS','','US','','bestM1','','bestM2',''});
                title(sprintf('%s - sess%d ',regname_cortex{r},ss));
                if ss==1
                    ylabel('logBayes');
                else
                    ylabel('');
                end
                idx=idx+1;
            end
        end
    case 'PCM_noiseCeiling_plot2' % used in the paper
        % here plot the marginal log bayes for the different factors
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        sessN = 4;
        reg = [2:3,7];
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','reg','naturalStats','sessN'});
        % plot PCM noise ceilings
        T = load(fullfile(pcmDir,sprintf('NoiseCeilings_repsup_%s.mat',parcelType)));
        %Th = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Stats_simple_thetas_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        K = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_simple_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_simple_%s.mat',parcelType)));
        idx=1;
        figure
        for r=reg
            for ss=sessN
                t=getrow(T,T.roi==r&T.sessN==ss);
                k=getrow(K,K.roi==r&K.sessN==ss);
                % transform bayes factors into pseudo-r2
                t.bayesEst_cross = bsxfun(@minus,t.cross_likelihood,t.cross_likelihood(:,1));
                k.bayesEst_cross = bsxfun(@minus,k.cross_likelihood,k.cross_likelihood(:,1));
                % form a new structure
                nModel = size(k.cross_likelihood,2);
                M.data_cross = k.bayesEst_cross(:);
                M.indx = kron([1:nModel]',ones(length(t.SN),1));
                M.sessN = repmat(k.sessN,nModel,1);
                M.exe = repmat(k.exe,nModel,1);
                M.SN = repmat(k.SN,nModel,1);
                M = normData(M,'data_cross');
                
                % go through all models - determine which significantly
                % different from noise ceiling
                % exe 1
                mIdx1 = []; pIdx1 = [];
                for i=1:size(k.time,2)
                    [tstat,p]=ttest(k.bayesEst(k.exe==1,i),t.bayesEst(t.exe==1,2),2,'paired');
                    if p>.05
                        mIdx1 = [mIdx1 i];
                        pIdx1 = [pIdx1 p(1)];
                    end
                end
                % exe 2
                 mIdx2 = []; pIdx2 = [];
                for i=1:size(k.time,2)
                    [tstat,p]=ttest(k.bayesEst(k.exe==2,i),t.bayesEst(t.exe==2,2),2,'paired');
                    if p>.05
                        mIdx2 = [mIdx2 i];
                        pIdx2 = [pIdx2 p(1)];
                    end
                end
                % out of the models which fit noise ceiling
                % here determine which model is best fitting
                % which models are significantly worse than that model
                [~,i1]=max(k.cross_likelihood(k.exe==1,mIdx1)');
                [~,i2]=max(k.cross_likelihood(k.exe==1,mIdx2)');
                BM1 = mIdx1(mode(i1));
                BM2 = mIdx2(mode(i2));
                mIdxW1 = []; pIdxW1 = [];
                for i=1:length(mIdx1) % exe1
                    [tstat,p]=ttest(k.bayesEst(k.exe==1,mIdx1(i)),k.bayesEst(k.exe==1,BM1),2,'paired');
                    if p>.01
                        mIdxW1 = [mIdxW1 mIdx1(i)];
                        pIdxW1 = [pIdxW1 p(1)];
                    end
                end
                mIdxW2 = []; pIdxW2 = [];
                for i=1:length(mIdx2) % exe2
                    [tstat,p]=ttest(k.bayesEst(k.exe==2,mIdx2(i)),k.bayesEst(k.exe==2,BM2),2,'paired');
                     if p>.01
                        mIdxW2 = [mIdxW2 mIdx2(i)];
                        pIdxW2 = [pIdxW2 p(1)];
                    end
                end
                keyboard;
                for exe=1:2
                    subplot(numel(reg),2,idx)
                    plt.bar(M.indx,M.data_cross,'subset',M.exe==exe); hold on;
                    drawline(mean(t.bayesEst_cross(t.exe==exe,2)),'dir','horz');
                    drawline(mean(t.bayesEst_cross(t.exe==exe,2)+stderr(t.bayesEst_cross(t.exe==exe,2))),'dir','horz');
                    drawline(mean(t.bayesEst_cross(t.exe==exe,2)-stderr(t.bayesEst_cross(t.exe==exe,2))),'dir','horz');
                    title(sprintf('%s - sess%d exe-%d',regname_cortex{r},ss,exe));
                    ylabel('Log-BF');
                    idx=idx+1;
                end                
            end
        end
    case 'PCM_noiseCeiling_plot2_noError'
        % glm without errors
        parcelType = 'Brodmann'; % Brodmann, 162tessels
        sessN = 4;
        reg = [2:3,7];
        naturalStats = 1;
        vararginoptions(varargin,{'parcelType','reg','naturalStats','sessN'});
        % plot PCM noise ceilings
        T = load(fullfile(pcmDir,sprintf('NoiseCeilings_repsup_noError_%s.mat',parcelType)));
        K = load(fullfile(pcmDir,sprintf('FoSEx_ModelFamily_Fit_noError_%s_naturalStats-%d.mat',parcelType,naturalStats)));
        load(fullfile(pcmDir,sprintf('FoSEx_ModelFamilyComb_simple_%s.mat',parcelType)));
        idx=1;
        figure
        for r=reg
            for ss=sessN
                t=getrow(T,T.roi==r&T.sessN==ss);
                k=getrow(K,K.roi==r&K.sessN==ss);
                % transform bayes factors into pseudo-r2
                t.bayesEst_cross = bsxfun(@minus,t.cross_likelihood,t.cross_likelihood(:,1));
                k.bayesEst_cross = bsxfun(@minus,k.cross_likelihood,k.cross_likelihood(:,1));
                % form a new structure
                nModel = size(k.cross_likelihood,2);
                M.data_cross = k.bayesEst_cross(:);
                M.indx = kron([1:nModel]',ones(length(t.SN),1));
                M.sessN = repmat(k.sessN,nModel,1);
                M.exe = repmat(k.exe,nModel,1);
                M.SN = repmat(k.SN,nModel,1);
                M = normData(M,'data_cross');
                
                % go through all models - determine which significantly
                % different from noise ceiling
                % exe 1
                mIdx1 = []; pIdx1 = [];
                for i=1:size(k.time,2)
                    [tstat,p]=ttest(k.bayesEst(k.exe==1,i),t.bayesEst(t.exe==1,2),2,'paired');
                    if p>.05
                        mIdx1 = [mIdx1 i];
                        pIdx1 = [pIdx1 p(1)];
                    end
                end
                % exe 2
                 mIdx2 = []; pIdx2 = [];
                for i=1:size(k.time,2)
                    [tstat,p]=ttest(k.bayesEst(k.exe==2,i),t.bayesEst(t.exe==2,2),2,'paired');
                    if p>.05
                        mIdx2 = [mIdx2 i];
                        pIdx2 = [pIdx2 p(1)];
                    end
                end
                % out of the models which fit noise ceiling
                % here determine which model is best fitting
                % which models are significantly worse than that model
                [~,i1]=max(k.cross_likelihood(k.exe==1,mIdx1)');
                [~,i2]=max(k.cross_likelihood(k.exe==1,mIdx2)');
                BM1 = mIdx1(mode(i1));
                BM2 = mIdx2(mode(i2));
                mIdxW1 = []; pIdxW1 = [];
                for i=1:length(mIdx1) % exe1
                    [tstat,p]=ttest(k.bayesEst(k.exe==1,mIdx1(i)),k.bayesEst(k.exe==1,BM1),2,'paired');
                    if p>.01
                        mIdxW1 = [mIdxW1 mIdx1(i)];
                        pIdxW1 = [pIdxW1 p(1)];
                    end
                end
                mIdxW2 = []; pIdxW2 = [];
                for i=1:length(mIdx2) % exe2
                    [tstat,p]=ttest(k.bayesEst(k.exe==2,mIdx2(i)),k.bayesEst(k.exe==2,BM2),2,'paired');
                     if p>.01
                        mIdxW2 = [mIdxW2 mIdx2(i)];
                        pIdxW2 = [pIdxW2 p(1)];
                    end
                end
                keyboard;
                for exe=1:2
                    subplot(numel(reg),2,idx)
                    plt.bar(M.indx,M.data_cross,'subset',M.exe==exe); hold on;
                    drawline(mean(t.bayesEst_cross(t.exe==exe,2)),'dir','horz');
                    drawline(mean(t.bayesEst_cross(t.exe==exe,2)+stderr(t.bayesEst_cross(t.exe==exe,2))),'dir','horz');
                    drawline(mean(t.bayesEst_cross(t.exe==exe,2)-stderr(t.bayesEst_cross(t.exe==exe,2))),'dir','horz');
                    title(sprintf('%s - sess%d exe-%d',regname_cortex{r},ss,exe));
                    ylabel('Log-BF');
                    idx=idx+1;
                end                
            end
        end
    case 'PCM_noiseCeilings_seqType'
        % construct noise ceilings - separately for trained and untrained
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
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            for r = 1:length(reg)
                for exe = 1:2 % sequence type
                    for st=1:2 % seqType
                        for g = 1:2 % two groups of subjects
                            for p=1:length(sn_group{g})
                                glmDirSubj=fullfile(glmFoSExDir{ss}, subj_name{sn_group{g}(p)});
                                D=load(fullfile(glmDirSubj,'SPM_info.mat'));
                                switch (beta_choice)
                                    case 'uw'
                                        beta = B.betaUW{(B.sn==sn_group{g}(p)&B.region==reg(r))}';
                                    case 'mw'
                                        beta = B.betaW{(B.SN==sn_group{g}(p)&B.region==reg(r))}';
                                    case 'raw'
                                        beta = B.betaRAW{(B.sn==sn_group{g}(p)&B.region==reg(r))}'; % no intercept - use T.betaRAWint otherwise
                                end
                                partVec{p} = D.run(D.FoSEx==exe & D.seqType==st);  % runs/partitions;
                                condVec{p} = D.seqNumb(D.FoSEx==exe & D.seqType==st);
                                Data{p} = beta(:,(D.FoSEx==exe & D.seqType==st))';  % Data is N x P (cond x voxels) - no intercept
                            end;
                            % fit models
                            T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                            T.roi = ones(size(T.SN))*reg(r);
                            T.regType = ones(size(T.SN))*regType(r);
                            T.regSide = ones(size(T.SN))*regSide(r);
                            T.sessN = ones(size(T.SN))*ss;
                            T.seqType = ones(size(T.SN))*st;
                            T.exe = ones(size(T.SN))*exe;
                            T=rmfield(T,{'reg','thetaCr','theta_hat','G_pred'});
                            if g==2
                                T.SN = T.SN + length(sn_group{1});
                            end
                            AllReg=addstruct(AllReg,T);
                            clear partVec condVec Data;
                        end % group
                    end % seqType
                end % repetition
            end
            fprintf('Done sess-%d.\n\n\n',ss);
        end
        % save variables;
        save(fullfile(pcmDir,sprintf('NoiseCeilings_seqType_repsup_%s.mat',parcelType)),'-struct','AllReg');   
    case 'PCM_noiseCeilings_seqType_plot'
        sessN = 1:4;
        reg = 1:8;
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sessN','reg','parcelType'});
        T = load(fullfile(pcmDir,sprintf('NoiseCeilings_seqType_repsup_%s.mat',parcelType)));
        for r=reg
            figure
            for ss=sessN
                t = getrow(T,T.roi==r & T.sessN==ss);
                subplot(2,4,ss)
                plt.bar(t.seqType,t.likelihood(:,2)-t.likelihood(:,1),'split',t.exe,'leg',{'exe1','exe2'},'leglocation','northwest');
                title(sprintf('High noise ceiling %s sess-%d',regname_cortex{r},ss)); ylabel('log-BF');
                subplot(2,4,ss+4)
                plt.bar(t.seqType,t.cross_likelihood(:,2)-t.cross_likelihood(:,1),'split',t.exe,'leg',{'exe1','exe2'},'leglocation','northwest');
                title(sprintf('Low noise ceiling %s sess-%d',regname_cortex{r},ss)); ylabel('log-BF');
            end
        end
        
    case 'PCM:scaling_model'
        % here create a scaling model of execution 2 with respect to exe1
        parcelType = 'Brodmann';
        runEffect = 'fixed';
        sessN = 1:4;
        reg = [1:3,7,8];
        sn = [5:9,11:31];
        scaleLim = [0 1];
        regType = 'clusters'; % all or subset
        nModel = 20;
        vararginoptions(varargin,{'parcelType','sessN','reg','scaleLim','nModel'});
        scaleL = linspace(scaleLim(1),scaleLim(2),nModel); % correlation models to assess
        H=eye(12)-ones(12)./12;  % centering matrix!
        DD = [];
        if strcmp(regType,'clusters')
            regID = 0;
        end
        for ss=sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            if strcmp(regType,'clusters') || regID == 0
                regAll = unique(B.region)';
                idx=1;
                for i=regAll
                    B1 = getrow(B,B.region==i);
                    for s=1:size(B1.betaW,1)
                        hasData(s) = 1-sum(any(isnan(B1.betaW{s}))); % whether participant has data (not nan)
                    end
                    if sum(hasData) == size(B1.betaW,1)
                        reg(idx) = i;
                        idx=idx+1;
                    end
                end
                regID = 1;
            end
            TI = load(fullfile(glmFoSExDir{ss},subj_name{sn(1)},'SPM_info.mat'));
            for r=reg
                for p=1:length(sn)
                    beta = B.betaW{(B.SN==sn(p) & B.region==r)}';
                    [G,~] = pcm_estGCrossval(beta(:,TI.FoSEx==1)',TI.run(TI.FoSEx==1),TI.seqNumb(TI.FoSEx==1));
                    partVec{1} = TI.run(TI.FoSEx==2);  % runs/partitions
                    condVec{1} = TI.seqNumb(TI.FoSEx==2);  % runs/partitions
                    Data{1} = beta(:,TI.FoSEx==2)';
                    G_centr = H*G*H';
                    for ii=1:nModel
                        M{ii}.type           = 'fixed';
                        M{ii}.numGparams     = 0;
                        M{ii}.name           = sprintf('scale-%1.2f',scaleL(ii));
                        M{ii}.Gc = G_centr.*~eye(size(G_centr))*scaleL(ii) + G_centr.*eye(size(G_centr));
                    end
                    [D,~,~] = pcm_fitModelIndivid(Data(1),M,partVec(1),condVec(1),'runEffect',runEffect,'verbose',0);
                    D.SN = p;
                    D.roi=r;
                    D.sessN=ss;
                    DD=addstruct(DD,D);
                    clear Data beta;
                end
                fprintf('Done reg-%d.\n\n\n',r);
            end
            fprintf('Done all regions, session %d.\n\n\n\n',ss);
        end
        save(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s',parcelType)),'-struct','DD');
    case 'PCM:scaling_model_seqType'
        % here split into trained / untrained separately
        parcelType = 'Brodmann';
        runEffect = 'fixed';
        sessN = 1:4;
        reg = [1:3,7];
        sn = [5:9,11:31];
        scaleLim = [0 1];
        nModel = 20;
        vararginoptions(varargin,{'parcelType','sessN','reg','scaleLim','nModel'});
        scaleL = linspace(scaleLim(1),scaleLim(2),nModel); % correlation models to assess
        H=eye(6)-ones(6)./6;  % centering matrix!
        DD = [];
        for ss=sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            TI = load(fullfile(glmFoSExDir{ss},subj_name{sn(1)},'SPM_info.mat'));
            for r=reg
                for st=1:2
                    for p=1:length(sn)
                        beta = B.betaW{(B.SN==sn(p) & B.region==r)}';
                        [G,~] = pcm_estGCrossval(beta(:,TI.FoSEx==1 & TI.seqType==st)',TI.run(TI.FoSEx==1 & TI.seqType==st),TI.seqNumb(TI.FoSEx==1 & TI.seqType==st));
                        partVec{1} = TI.run(TI.FoSEx==2 & TI.seqType==st);  % runs/partitions
                        condVec{1} = TI.seqNumb(TI.FoSEx==2 & TI.seqType==st);  % runs/partitions
                        Data{1} = beta(:,TI.FoSEx==2 & TI.seqType==st)';
                        G_centr = H*G*H';
                        for ii=1:nModel
                            M{ii}.type           = 'fixed';
                            M{ii}.numGparams     = 0;
                            M{ii}.name           = sprintf('scale-%1.2f',scaleL(ii));
                            M{ii}.Gc = G_centr.*~eye(size(G_centr))*scaleL(ii) + G_centr.*eye(size(G_centr));
                        end
                        [D,~,~] = pcm_fitModelIndivid(Data(1),M,partVec(1),condVec(1),'runEffect',runEffect);
                        D.SN = p;
                        D.roi=r;
                        D.seqType=st;
                        D.sessN=ss;
                        DD=addstruct(DD,D);
                        clear Data beta;
                    end
                    fprintf('Done reg-%d, seqType-%d.\n\n\n',r,st);
                end
                fprintf('Done reg-%d.\n\n\n',r);
            end
        end
        save(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s_seqType',parcelType)),'-struct','DD');
    case 'PCM:scaling_model_seq'
        % here separate scalar for trained / untrained
        % here create a scaling model of execution 2 with respect to exe1
        parcelType = 'Brodmann'; % Brodmann, tesselsWB_162
        runEffect = 'fixed';
        sessN = 1:4;
        reg = [1:3,7,8];
        regType = 'clusters'; % all or subset
        sn = [5:9,11:31];
        scaleLim = [0 1];
        nModel = 20;
        vararginoptions(varargin,{'parcelType','sessN','reg','scaleLim','nModel'});
        scaleL = linspace(scaleLim(1),scaleLim(2),nModel); % correlation models to assess
        H=eye(12)-ones(12)./12;  % centering matrix!
        DD = [];
        for ss=sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            if strcmp(regType,'clusters')
                regAll = unique(B.region)';
                idx=1;
                for i=regAll
                    B1 = getrow(B,B.region==i);
                    for s=1:size(B1.betaW,1)
                        hasData(s) = 1-sum(any(isnan(B1.betaW{s}))); % whether participant has data (not nan)
                    end
                    if sum(hasData) == size(B1.betaW,1)
                        reg(idx) = i;
                        idx=idx+1;
                    end
                end
            end
            TI = load(fullfile(glmFoSExDir{ss},subj_name{sn(1)},'SPM_info.mat'));
            for r=reg
                for p=1:length(sn)
                    beta = B.betaW{(B.SN==sn(p) & B.region==r)}';
                    [G,~] = pcm_estGCrossval(beta(:,TI.FoSEx==1)',TI.run(TI.FoSEx==1),TI.seqNumb(TI.FoSEx==1));
                    partVec{1} = TI.run(TI.FoSEx==2);  % runs/partitions
                    condVec{1} = TI.seqNumb(TI.FoSEx==2);  % runs/partitions
                    Data{1} = beta(:,TI.FoSEx==2)';
                    G_centr = H*G*H';
                    idx=1;
                    for i=1:nModel
                        for j=1:nModel
                            M{idx}.type           = 'fixed';
                            M{idx}.numGparams     = 0;
                            M{idx}.name           = sprintf('scale-trained%1.2f-untrained%1.2f',scaleL(i),scaleL(j));
                            % separately construct G from trained /
                            % untrained / correlation
                            Gc = G_centr;
                            Gc(1:6,1:6) = G_centr(1:6,1:6).*~eye(6)*scaleL(i);
                            Gc(7:12,7:12) = G_centr(7:12,7:12).*~eye(6)*scaleL(j);
                            corTU = sqrt(scaleL(i)*scaleL(j));
                            Gc(1:6,7:12) = G_centr(1:6,7:12).*corTU;
                            Gc(7:12,1:6) = G_centr(7:12,1:6).*corTU;
                            M{idx}.Gc = Gc;
                            trComb(idx) = i;
                            utrComb(idx) = j;
                            idx=idx+1;
                        end
                    end
                    [D,~,~] = pcm_fitModelIndivid(Data(1),M,partVec(1),condVec(1),'runEffect',runEffect,'verbose',0);
                    D.trainedScale = trComb;
                    D.untrainedScale = utrComb;
                    D.SN = p;
                    D.roi=r;
                    D.sessN=ss;
                    DD=addstruct(DD,D);
                    clear Data beta;
                end
                fprintf('Done reg-%d.\n\n\n',r);
            end
            fprintf('Done all regions, session %d.\n\n\n\n',ss);
        end
        save(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s_seq',parcelType)),'-struct','DD');
    case 'PCM:scaling_check'
        parcelType = 'Brodmann';
        reg = [1:3,7];
        sn = [5:9,11:31];
        sessN = 4;
        vararginoptions(varargin,{'parcelType','sessN','reg','scaleLim','nModel'});
        TT = [];
        for ss=sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            TI = load(fullfile(glmFoSExDir{ss},subj_name{sn(1)},'SPM_info.mat'));
            H=eye(12)-ones(12)./12;  % centering matrix!
            for r=reg
                for p=1:length(sn)
                    beta = B.betaW{(B.SN==sn(p) & B.region==r)}';
                    [G1,~] = pcm_estGCrossval(beta(:,TI.FoSEx==1)',TI.run(TI.FoSEx==1),TI.seqNumb(TI.FoSEx==1));
                    [G2,~] = pcm_estGCrossval(beta(:,TI.FoSEx==2)',TI.run(TI.FoSEx==2),TI.seqNumb(TI.FoSEx==2));
                    B1 = mean(mean(beta(:,TI.FoSEx==1)));
                    B2 = mean(mean(beta(:,TI.FoSEx==2)));
                    D1 = rsa.distanceLDC(beta(:,TI.FoSEx==1)',TI.run(TI.FoSEx==1),TI.seqNumb(TI.FoSEx==1));
                    D2 = rsa.distanceLDC(beta(:,TI.FoSEx==2)',TI.run(TI.FoSEx==2),TI.seqNumb(TI.FoSEx==2));
                    GG1 = H*G1*H';
                    GG2 = H*G2*H';
                    scaleG = GG1./GG2;
                    scaleG2 = G1./G2;
                    T.sn = sn(p);
                    T.reg = r;
                    T.sessN = ss;
                    T.offDiag = mean(scaleG(triu(scaleG,1)>0));
                    T.diag = mean(diag(scaleG));
                    T.offDiag2 = mean(scaleG2(triu(scaleG2,1)>0));
                    T.diag2 = mean(diag(scaleG2));
                    T.distScale = mean(D1./D2);
                    T.actScale = B1/B2;
                    TT = addstruct(TT,T);
                end
            end
            keyboard;
        end
        figure
        scatterplot(TT.diag,TT.offDiag,'subset',TT.reg==3);hold on;
        g = gca;
        maxLim = max(g.XLim(2),g.YLim(2))+0.01;
        minLim = min(g.XLim(1),g.YLim(1))-0.01;
        xlim([minLim maxLim]); ylim([minLim maxLim]);
        x = linspace(minLim,maxLim);
        plot(x,x,'-k');
    case 'PCM:scaling_plot'
        % here plot the scaling results
        parcelType = 'Brodmann';
        reg = [1:3,7];
        sessN = 1:4;
        T = load(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s',parcelType)));
        T.bayesEst = bsxfun(@minus,T.likelihood,mean(T.likelihood,2));
        DD = [];
        for r=reg
            for ss=sessN
                T1 = getrow(T,T.roi==r & T.sessN==ss);
                % reshape
                nModel      = size(T1.bayesEst,2);
                D.bayesEst  = T1.bayesEst(:);
                D.SN        = repmat(T1.SN,nModel,1);
                D.model     = kron((1:nModel)',ones(size(T1.bayesEst,1),1));
                D.reg       = ones(size(D.model))*r;
                D.sessN     = ones(size(D.model))*ss;
                DD = addstruct(DD,D);
                figure(r)
                subplot(1,numel(sessN),find(ss==sessN))
                plt.line(D.model,D.bayesEst,'style',stySeqShade);
                hold on;
                drawline(0,'dir','horz');
                xlabel('Scaling model'); ylabel('Evidence'); title(sprintf('%s-session %d',regname_cortex{r},ss));
                t1 = gca;
                corrTick = linspace(min(t1.XTick),max(t1.XTick),11);
                set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
            end
        end
        D1 = getrow(DD,ismember(DD.reg,[2,3,7]));
        figure
        subplot(122)
        plt.line(D1.model,D1.bayesEst,'split',D1.reg,'subset',D1.sessN==4,'style',stySeqReg,'leg',{'M1','PMd','SPLa'});
        hold on; drawline(0,'dir','horz'); xlabel('Scaling model'); ylabel('Evidence'); title('session 4');
        t1 = gca;
        corrTick = linspace(min(t1.XTick),max(t1.XTick),11);
        set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
        subplot(121)
        plt.line(D1.model,D1.bayesEst,'split',D1.reg,'subset',D1.sessN==1,'style',stySeqReg,'leg',{'M1','PMd','SPLa'});
        hold on; drawline(0,'dir','horz'); xlabel('Scaling model'); ylabel('Evidence'); title('session 1');
        t1 = gca;
        corrTick = linspace(min(t1.XTick),max(t1.XTick),11);
        set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
    case 'PCM:scaling_plot_seqType'
        % here plot the scaling results, split across sequences
        parcelType = 'Brodmann';
        reg = [1:3,7];
        sessN = 1:4;
        T = load(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s_seqType',parcelType)));
        for r=reg
            for ss=sessN
                T1 = getrow(T,T.roi==r & T.sessN==ss);
                T1.bayesEst = bsxfun(@minus,T1.likelihood,mean(T1.likelihood,2));
                % reshape
                nModel      = size(T1.bayesEst,2);
                D.bayesEst  = T1.bayesEst(:);
                D.SN        = repmat(T1.SN,nModel,1);
                D.model     = kron((1:nModel)',ones(size(T1.bayesEst,1),1));
                D.seqType   = repmat(T1.seqType,nModel,1);
                figure(r)
                subplot(1,numel(sessN),find(ss==sessN))
                plt.line(D.model,D.bayesEst,'split',D.seqType,'style',stySeqShade);
                hold on;
                drawline(0,'dir','horz');
                xlabel('Scaling model'); ylabel('Evidence'); title(sprintf('%s-session %d',regname_cortex{r},ss));
                t1 = gca;
                corrTick = linspace(min(t1.XTick),max(t1.XTick),11);
                set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
            end
        end
    case 'PCM:scaling_stats'
        % here compare the best fitting model across regions
        parcelType = 'Brodmann';
        reg = [2:3,7];
        sessN = 4;
        vararginoptions(varargin,{'reg','sessN'});
        T = load(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s',parcelType)));
        T.bayesEst = bsxfun(@minus,T.likelihood,mean(T.likelihood,2));
        DD = [];
        for r=reg
            for ss=sessN
                T1 = getrow(T,T.roi==r & T.sessN==ss);
                [~,maxModel] = max(T1.bayesEst'); % identify max evidence model
                D.maxModel = maxModel';
                D.SN = [1:size(maxModel,2)]';
                D.reg = ones(size(D.SN))*r;
                D.sessN = ones(size(D.SN))*ss;
                DD = addstruct(DD,D);
            end
        end
        ttestDirect(DD.maxModel,[DD.reg DD.SN],2,'paired','subset',ismember(DD.reg,[2,3]));
        ttestDirect(DD.maxModel,[DD.reg DD.SN],2,'paired','subset',ismember(DD.reg,[2,7]));
        ttestDirect(DD.maxModel,[DD.reg DD.SN],2,'paired','subset',ismember(DD.reg,[3,7]));
    case 'PCM:scaling_seq'
        % scaling models
        reg = [1:3,7,8];
        parcelType = 'Brodmann';
        sessN = 1:4;
        vararginoptions(varargin,{'reg','parcelType','sessN'});
        
        T = load(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s_seq',parcelType)));
        T.bayesEst = bsxfun(@minus,T.likelihood,mean(T.likelihood,2));
        KK=[];
        for r=reg
            figure
            for ss=sessN
                subplot(2,2,ss)
                t = getrow(T,T.sessN==ss & T.roi==r);
                meanScale = mean(t.bayesEst);
                for i=1:size(meanScale,2)
                    R(t.trainedScale(1,i),t.untrainedScale(1,i))=meanScale(i);
                end
                imagesc(R);
                [y,~]=max(max(R));
                [k,l]=find(R==y);
                hold on; plot(l,k,'x');
                title(sprintf('%s-sess%d',regname_cortex{r},ss));
                % per subject
                [~,i1] = max(t.bayesEst');
                K.max = [t.trainedScale(1,i1)';t.untrainedScale(1,i1)']/20;
                K.seqType = [ones(size(t.SN,1),1);ones(size(t.SN,1),1)*2];
                K.sn = [t.SN;t.SN];
                K.reg = ones(size(K.sn))*r;
                K.sessN = ones(size(K.sn))*ss;
                KK=addstruct(KK,K);
            end
        end
        keyboard;
    case 'PCM:scaling_seq_surface'
        % project to the surface
        sessN = 4;
        parcelType = 'tesselsWB_162';
        vararginoptions(varargin,{'sessN'});
        
        T = load(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s_seq',parcelType)));
        T.bayesEst = bsxfun(@minus,T.likelihood,mean(T.likelihood,2));
        colIdx = 1;
        for h=1
            Ico = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-162.164k.%s.label.gii',hemI{h})));
            data = zeros(size(Ico.cdata));
            reg = unique(T.roi(T.sessN==1))';
            for ss=sessN
                T1 = getrow(T,T.sessN==ss);
                for i=reg
                    t = getrow(T1, T1.roi==i);
                    meanScale = mean(t.bayesEst);
                    [~,idx]=max(meanScale);
                    data(Ico.cdata==i,colIdx) = t.trainedScale(1,idx)/max(max(t.trainedScale));
                    data(Ico.cdata==i,colIdx+1) = t.untrainedScale(1,idx)/max(max(t.trainedScale));
                    data(Ico.cdata==i,colIdx+2) = t.trainedScale(1,idx)/max(max(t.trainedScale)) - t.untrainedScale(1,idx)/max(max(t.trainedScale));
                end
                colName{colIdx} = sprintf('trained_scale-sess%d',ss);
                colName{colIdx+1} = sprintf('untrained_scale-sess%d',ss);
                colName{colIdx+2} = sprintf('difference_scale-sess%d',ss);
                colIdx=colIdx+3;
            end
            G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',colName);
            outfile = fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('%s.PCM_scalingModel_seq.func.gii',hemI{h}));
            save(G,outfile);
        end
    case 'PCM:scaling_surface'
        sessN = 4;
        parcelType = 'tesselsWB_162';
        vararginoptions(varargin,{'sessN'});
        
        T = load(fullfile(pcmDir,sprintf('FoSEx_scalingModel_%s',parcelType)));
        T.bayesEst = bsxfun(@minus,T.likelihood,mean(T.likelihood,2));
        colIdx = 1;
        for h=1
            Ico = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-162.164k.%s.label.gii',hemI{h})));
            data = zeros(size(Ico.cdata));
            reg = unique(T.roi(T.sessN==1))';
            for ss=sessN
                T1 = getrow(T,T.sessN==ss);
                for i=reg
                    t = getrow(T1, T1.roi==i);
                    meanScale = mean(t.bayesEst);
                    [~,idx]=max(meanScale);
                    data(Ico.cdata==i,colIdx) = idx/size(t.likelihood,2);
                end
                colName{colIdx} = sprintf('scale-sess%d',ss);
                colIdx=colIdx+1;
            end
            G = surf_makeFuncGifti(data,'anatomicalStruct',hemname{h},'columnNames',colName);
            outfile = fullfile(wbDir,'FoSEx_FS_LR_164',sprintf('%s.PCM_scalingModel.func.gii',hemI{h}));
            save(G,outfile);
        end
        
    case 'PCM:correlation'
        % correlate executions 1 & 2
        parcelType = 'Brodmann';
        runEffect = 'random';
        sessN = 1:4;
        reg = [1:3,7,8];
        sn = [5:9,11:31];
        corrLim = [0 1];
        nModel = 30;
        withinCov = 'iid'; % iid or individual
        algorithm = 'NR';
        type = 'overall'; % overall - treat trained/untrained together; seqType - split
        vararginoptions(varargin,{'parcelType','sessN','reg','corrLim','nModel','type','withinCov','runEffect'});
        corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
        for c=1:length(corrS)
            switch type
                case 'overall'
                    M{c} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',2,'numItems',12,'r',corrS(c));
                case 'seqType'
                    M{c} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',4,'numItems',6,'r',corrS(c));
                    % seqType has different thetas for trained / untrained,
                    % but correlation is still the same
                case 'overall_flexible' % flexibly determine best r
                    M{1} = pcm_buildCorrModel('type','nonlinear','withinCov',withinCov,'numCond',2,'numItems',12,'r','flexible');
            end
        end
        DD = [];
        for ss=sessN
            B=load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            TI = load(fullfile(glmFoSExDir{ss},subj_name{sn(1)},'SPM_info.mat'));
            for r=reg
                for p=1:length(sn)
                    t = getrow(B,B.SN==sn(p) & B.region==r);
                    beta = B.betaW{(B.SN==sn(p) & B.region==r)};
                    partVec{p} = TI.run;  % runs/partitions
                    cond = TI.seqNumb;  % runs/partitions
                    cond(TI.FoSEx==2) = cond(TI.FoSEx==2)+max(cond);
                    condVec{p} = cond;
                    Data{p} = beta(1:size(TI.seqNumb),:);
                end
                D = sml1_imana_repsup('PCM:corrModel_run','Data',Data,'model',M,'partV',partVec,'condV',condVec,'algorithm',algorithm,'runEffect',runEffect);
                D.roi       = ones(size(D.SN))*r;
                D.regType   = ones(size(D.SN))*t.regType;
                D.regSide   = ones(size(D.SN))*t.regSide;
                D.sessN     = ones(size(D.SN))*ss;
                D = rmfield(D,'reg');
                DD=addstruct(DD,D);
                fprintf('Done reg-%d, sess-%d.\n\n\n',r,ss);
            end
            fprintf('Done all regions, session %d.\n\n\n\n',ss);
        end
        save(fullfile(pcmDir,sprintf('FoSEx_correlationModel_%s_covariance-%s_%s',parcelType,withinCov,type)),'-struct','DD');
    case 'PCM:corrModel_run'
        algorithm='NR';
        S = [];
        vararginoptions(varargin,{'Data','model','partV','condV','runEffect','algorithm','S'});
        if ~isempty(S)
            [T,~] = pcm_fitModelGroup(Data,model,partV,condV,'runEffect',runEffect,'fitAlgorithm',algorithm,'fitScale',1,'S',S); % only group fit
        else
            [T,~] = pcm_fitModelGroup(Data,model,partV,condV,'runEffect',runEffect,'fitAlgorithm',algorithm,'fitScale',1); % only group fit
        end
        T.bayesEst = bsxfun(@minus,T.likelihood,T.likelihood(:,1)); % now relative to model with 0 correlation
        T.mean_likelihood = mean(T.likelihood,2);
        T.posterior = exp(T.bayesEst);
        T.posterior = bsxfun(@rdivide,T.posterior,sum(T.posterior,2));
        % individual fit
      %  if ~isempty(S)
       %     [I,~] = pcm_fitModelIndivid(Data,model,partV{1},condV{1},'runEffect',runEffect,'S',S,'verbose',0);
      %  else
      %      [I,~] = pcm_fitModelIndivid(Data,model,partV{1},condV{1},'runEffect',runEffect,'S',S,'verbose',0);
      %  end
      %   T.individ_likelihood    = I.likelihood;
      %   T.bayesEst_individ      = bsxfun(@minus,T.individ_likelihood,T.individ_likelihood(:,1));
         %varargout={T,thetaG,thetaI};
         varargout{1}=T;
    case 'PCM:correlation_plot'
        reg = [1:3,7,8];
        sessN = 1:4;
        parcelType = 'Brodmann';
        withinCov = 'iid';
        type = 'overall';
        metric = 'bayesEst'; % bayesEst or bayesEst_individ
        vararginoptions(varargin,{'reg','sessN','withinCov','type','metric'});
        T = load(fullfile(pcmDir,sprintf('FoSEx_correlationModel_%s_covariance-%s_%s',parcelType,withinCov,type)));
        
        for r=reg
            for ss=sessN
                t = getrow(T,T.roi==r & T.sessN==ss);
                t.metric2 = bsxfun(@minus,t.(metric),mean(t.(metric),2));
                % reshape
                nModel      = size(t.(metric),2);
                D.(metric)  = t.(metric)(:);
                D.metric2   = t.metric2(:);
                D.SN        = repmat(t.SN,nModel,1);
                D.model     = kron((1:nModel)',ones(size(t.(metric),1),1));
                figure(r)
                subplot(1,numel(sessN),find(ss==sessN))
                plt.line(D.model,D.metric2,'style',stySeqShade);
                drawline(0,'dir','horz');
                xlabel('Model correlation'); ylabel('Evidence'); title(sprintf('%s, sess-%d',regname_cortex{r},ss));
                t1 = gca;
                corrTick = linspace(min(t1.XTick),max(t1.XTick),11);
                set(gca,'XTickLabel',(0:.1:1),'XTick',corrTick);
            end
        end
    case 'CORR:acrSess_G'
        % here calculate correlation across 1st and 2exe
        % and across sessions
        % taking into account G, not PCM
        reg = [1:3,7,8];
        regType = 'all'; % all or subset
        sessN = 1:4;
        parcelType = 'Brodmann';
        sn = [5:9,11:31];
        vararginoptions(varargin,{'reg','sessN','parcelType'});
        
        DD = [];
        if strcmp(regType,'all');
            regIdx = 0;
        end
        for ss=sessN
            T{ss} = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,ss)));
            if regIdx==0
                reg = unique(T{ss}.region)';
                regIdx=1;
            end
        end
        H=eye(12)-ones(12)./12;
        for r=reg
             for ss1=sessN
                 for ss2=ss1:max(sessN)
                     for s=sn
                         t1 = getrow(T{ss1},T{ss1}.region==r & T{ss1}.SN==s);
                         t2 = getrow(T{ss2},T{ss2}.region==r & T{ss2}.SN==s);
                         % calculate Gs (+ double center)
                         G1_exe1 = H*rsa_squareIPMfull(t1.IPMfull(1,:))*H;
                         G1_exe2 = H*rsa_squareIPMfull(t1.IPMfull(2,:))*H;
                         G2_exe1 = H*rsa_squareIPMfull(t2.IPMfull(1,:))*H;
                         G2_exe2 = H*rsa_squareIPMfull(t2.IPMfull(2,:))*H;
                         % cross-correlate all
                         D.corr(1,:) = corr(rsa_vectorizeIPM(G1_exe1)',rsa_vectorizeIPM(G1_exe2)');
                         D.corr(2,:) = corr(rsa_vectorizeIPM(G2_exe1)',rsa_vectorizeIPM(G2_exe2)');
                         D.corr(3,:) = corr(rsa_vectorizeIPM(G1_exe1)',rsa_vectorizeIPM(G2_exe1)');
                         D.corr(4,:) = corr(rsa_vectorizeIPM(G1_exe2)',rsa_vectorizeIPM(G2_exe2)');
                         D.corr(5,:) = corr(rsa_vectorizeIPM(G1_exe1)',rsa_vectorizeIPM(G2_exe2)');
                         D.corr(6,:) = corr(rsa_vectorizeIPM(G1_exe2)',rsa_vectorizeIPM(G2_exe1)');
                         % other info
                         D.sess1 = [ss1 ss2 ss1 ss1 ss1 ss1]';
                         D.sess2 = [ss1 ss2 ss2 ss2 ss2 ss2]';
                         D.exe1 = [1 1 1 2 1 2]';
                         D.exe2 = [2 2 1 2 2 1]';
                         D.region = ones(size(D.exe1))*r;
                         D.sn = ones(size(D.exe1))*s;
                         DD = addstruct(DD,D);
                     end
                 end
             end
             fprintf('Done reg %d/%d.\n',find(r==reg),numel(reg));
        end
        save(fullfile(pcmDir,sprintf('FoSEx_correlationModel_G_%s',parcelType)),'-struct','DD');
    case 'CORR:acrSess_G_seqType'
        % split into trained / untrained
        % here calculate correlation across 1st and 2exe
        % and across sessions
        % taking into account G, not PCM
        reg = [1:3,7,8];
        sessN = 1:4;
        regType = 'all'; % all or subset
        parcelType = 'Brodmann';
        sn = [5:9,11:31];
        vararginoptions(varargin,{'reg','sessN','regType','parcelType'});
        
        DD = [];
        if strcmp(regType,'all');
            regIdx = 0;
        end
        for ss=sessN
            T{ss} = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_multiPW_sess%d.mat',parcelType,ss)));
            if regIdx==0
                reg = unique(T{ss}.region)';
                regIdx=1;
            end
        end
        H=eye(6)-ones(6)./6;
        for r=reg
            for st=1:2
                for ss1=sessN
                    for ss2=ss1:max(sessN)
                        for s=sn
                            t1 = getrow(T{ss1},T{ss1}.region==r & T{ss1}.SN==s);
                            t2 = getrow(T{ss2},T{ss2}.region==r & T{ss2}.SN==s);
                            % calculate Gs (+ double center)
                            G1_exe1 = rsa_squareIPMfull(t1.IPMfull(1,:));
                            G1_exe2 = rsa_squareIPMfull(t1.IPMfull(2,:));
                            G2_exe1 = rsa_squareIPMfull(t2.IPMfull(1,:));
                            G2_exe2 = rsa_squareIPMfull(t2.IPMfull(2,:));
                            if st==1
                                G1_exe1 = G1_exe1(1:6,1:6);
                                G1_exe2 = G1_exe2(1:6,1:6);
                                G2_exe1 = G2_exe1(1:6,1:6);
                                G2_exe2 = G2_exe2(1:6,1:6);
                            else
                                G1_exe1 = G1_exe1(7:12,7:12);
                                G1_exe2 = G1_exe2(7:12,7:12);
                                G2_exe1 = G2_exe1(7:12,7:12);
                                G2_exe2 = G2_exe2(7:12,7:12);
                            end
                            G1_exe1 = H*G1_exe1*H;
                            G1_exe2 = H*G1_exe2*H;
                            G2_exe1 = H*G2_exe1*H;
                            G2_exe2 = H*G2_exe2*H;
                            % cross-correlate all
                            D.corr(1,:) = corr(rsa_vectorizeIPM(G1_exe1)',rsa_vectorizeIPM(G1_exe2)');
                            D.corr(2,:) = corr(rsa_vectorizeIPM(G2_exe1)',rsa_vectorizeIPM(G2_exe2)');
                            D.corr(3,:) = corr(rsa_vectorizeIPM(G1_exe1)',rsa_vectorizeIPM(G2_exe1)');
                            D.corr(4,:) = corr(rsa_vectorizeIPM(G1_exe2)',rsa_vectorizeIPM(G2_exe2)');
                            D.corr(5,:) = corr(rsa_vectorizeIPM(G1_exe1)',rsa_vectorizeIPM(G2_exe2)');
                            D.corr(6,:) = corr(rsa_vectorizeIPM(G1_exe2)',rsa_vectorizeIPM(G2_exe1)');
                            % other info
                            D.sess1 = [ss1 ss2 ss1 ss1 ss1 ss1]';
                            D.sess2 = [ss1 ss2 ss2 ss2 ss2 ss2]';
                            D.exe1 = [1 1 1 2 1 2]';
                            D.exe2 = [2 2 1 2 2 1]';
                            D.seqType = ones(size(D.exe1))*st;
                            D.region = ones(size(D.exe1))*r;
                            D.sn = ones(size(D.exe1))*s;
                            DD = addstruct(DD,D);
                        end
                    end
                end
                fprintf('Done seqType %d/%d, reg %d/%d.\n',st,2,find(r==reg),numel(reg));
            end
        end
        save(fullfile(pcmDir,sprintf('FoSEx_correlationModel_G_seqType_%s',parcelType)),'-struct','DD');
    case 'CORR:acrSess_G_plot'
        % here plot the results
        parcelType = 'Brodmann';
        reg = [1:3,7,8];
        D = load(fullfile(pcmDir,sprintf('FoSEx_correlationModel_G_%s',parcelType)));
        
        for r=reg
            d1 = getrow(D,D.region==r);
            figure(1)
            subplot(2,numel(reg),find(reg==r))
            plt.line([d1.sess1>3 d1.sess1],d1.corr,'subset',d1.sess1==d1.sess2,'style',stySess); ylabel('corr');
            title(sprintf('%s - within session exe1-2 corr',regname_cortex{r}));
            subplot(2,numel(reg),find(reg==r)+numel(reg))
            plt.line([d1.sess1>3 d1.sess1],d1.corr,'split',d1.exe1,'subset',d1.sess2==4 & d1.sess1~=4 & d1.exe1==d1.exe2,'style',stySess);
            title('relative to sess 4'); ylabel('corr');
            figure
            subplot(321)
            plt.line(d1.sess1,d1.corr,'subset',d1.sess1==d1.sess2); ylabel('corr');
            title(sprintf('%s - within session exe1-2 corr',regname_cortex{r}));
            subplot(323)
            plt.line(d1.sess2,d1.corr,'split',d1.exe1,'subset',d1.sess1==1 & d1.sess2~=1 & d1.exe1==d1.exe2);
            title('relative to sess 1'); ylabel('corr');
            subplot(324)
            plt.line(d1.sess1,d1.corr,'split',d1.exe1,'subset',d1.sess2==4 & d1.sess1~=4 & d1.exe1==d1.exe2);
            title('relative to sess 4'); ylabel('corr');
            subplot(325)
            plt.line(d1.sess1,d1.corr,'subset',d1.sess2==d1.sess1+1 & d1.exe1==1 & d1.exe2==2);
            title('subsequent sess, exe1-2'); ylabel('corr');
            subplot(326)
            plt.line(d1.sess1,d1.corr,'subset',d1.sess2==d1.sess1+1 & d1.exe1==2 & d1.exe2==1);
            title('subsequent sess, exe2-1'); ylabel('corr');
            figure
            plt.line(d1.sess1,d1.corr,'split',[d1.exe1 d1.exe2],'subset',d1.sess2==4 & d1.sess1~=4); title(sprintf(regname_cortex{r}));
            figure
            subplot(221)
            [F,~,~]=pivottable(d1.sess1, d1.sess2,d1.corr,'robustmean','subset',d1.exe1==1 & d1.exe2==1 & d1.sess1~=d1.sess2);
            imagesc(F); colorbar; title(sprintf('%s, exe1-1',regname_cortex{r}));
            subplot(222)
            [F,~,~]=pivottable(d1.sess1, d1.sess2,d1.corr,'robustmean','subset',d1.exe1==2 & d1.exe2==2 & d1.sess1~=d1.sess2);
            imagesc(F); colorbar; title(sprintf('%s, exe2-2',regname_cortex{r}));
            subplot(223)
            [F,~,~]=pivottable(d1.sess1, d1.sess2,d1.corr,'robustmean','subset',d1.exe1==1 & d1.exe2==2 & d1.sess1~=d1.sess2);
            imagesc(F); colorbar; title(sprintf('%s, exe1-2',regname_cortex{r}));
            subplot(224)
            [F,~,~]=pivottable(d1.sess1, d1.sess2,d1.corr,'robustmean','subset',d1.exe1==2 & d1.exe2==1 & d1.sess1~=d1.sess2);
            imagesc(F); colorbar; title(sprintf('%s, exe2-1',regname_cortex{r}));
        end
    case 'CORR:acrSess_G_seqType_plot'
        % here plot the results, split for sequence type
        parcelType = 'Brodmann';
        reg = [1:3,7,8];
        D = load(fullfile(pcmDir,sprintf('FoSEx_correlationModel_G_seqType_%s',parcelType)));
        
        for r=reg
            for st=1:2
                d1 = getrow(D,D.region==r & D.seqType==st);
                figure
                subplot(221)
                plt.line(d1.sess1,d1.corr,'subset',d1.sess1==d1.sess2); ylabel('corr');
                title(sprintf('%s - within session exe1-2 corr',regname_cortex{r}));
                subplot(222)
                plt.line(d1.sess1,d1.corr,'split',d1.exe1,'subset',d1.sess2==4 & d1.sess1~=4 & d1.exe1==d1.exe2);
                title(sprintf('relative to sess 4 - %s sequence',SeqType{st})); ylabel('corr');
                subplot(223)
                plt.line(d1.sess1,d1.corr,'subset',d1.sess2==d1.sess1+1 & d1.exe1==1 & d1.exe2==2);
                title('subsequent sess, exe1-2'); ylabel('corr');
                subplot(224)
                plt.line(d1.sess1,d1.corr,'subset',d1.sess2==d1.sess1+1 & d1.exe1==2 & d1.exe2==1);
                title('subsequent sess, exe2-1'); ylabel('corr');
            end
        end
    case 'CORR:acrSess_G_stats'
        parcelType = 'Brodmann';
        reg = [1:3,7,8];
        D = load(fullfile(pcmDir,sprintf('FoSEx_correlationModel_G_%s',parcelType)));
        for r=reg
            d1 = getrow(D,D.region==r);
            fprintf('%s\n',regname_cortex{r});
            ttestDirect(d1.corr,[d1.exe1,d1.sn],2,'paired','subset',d1.sess2==d1.sess1+1 & d1.exe1~=d1.exe2,'split',d1.sess1);
        end
    case 'CORR:acrSess_G_seqType_stats'
        parcelType = 'Brodmann';
        reg = [1:3,7,8];
        D = load(fullfile(pcmDir,sprintf('FoSEx_correlationModel_G_seqType_%s',parcelType)));
        for st=1:2
            for r=reg
                d1 = getrow(D,D.region==r);
                fprintf('%s - %s\n',regname_cortex{r},SeqType{st});
                ttestDirect(d1.corr,[d1.exe1,d1.sn],2,'paired','subset',d1.sess2==d1.sess1+1 & d1.exe1~=d1.exe2 & d1.seqType==st,'split',d1.sess1);
            end
        end
    case 'CORR:acrSess_pattern'
        % here correlate the patterns, not the G matrices
        reg = [1:3,7,8];
        regType = 'subset'; % all or subset
        sessN = 1:4;
        parcelType = 'Brodmann';
        sn = [5:9,11:27,29:31];
        vararginoptions(varargin,{'reg','sessN','parcelType'});
        
        DD = []; SS = [];
        if strcmp(regType,'all');
            regIdx = 0;
        else
            regIdx = 1;
        end
        for ss=sessN
            T{ss} = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d.mat',parcelType,ss)));
            if regIdx==0
                reg = unique(T{ss}.region)';
                regIdx=1;
            end
        end
        condVec = repmat([1:24]',8,1);
        for r=reg
             for ss1=sessN
                 for ss2=ss1:max(sessN)
                     for s=sn
                         t1 = getrow(T{ss1},T{ss1}.region==r & T{ss1}.SN==s);
                         t2 = getrow(T{ss2},T{ss2}.region==r & T{ss2}.SN==s);
                         % calculate mean patterns
                         beta1 = pinv(indicatorMatrix('identity',condVec))*t1.betaW{1}(1:size(condVec,1),:);
                         beta2 = pinv(indicatorMatrix('identity',condVec))*t2.betaW{1}(1:size(condVec,1),:);
                         beta1(1:12,:) = bsxfun(@minus,beta1(1:12,:),mean(beta1(1:12,:),1));
                         beta1(13:24,:) = bsxfun(@minus,beta1(13:24,:),mean(beta1(13:24,:),1));
                         beta2(1:12,:) = bsxfun(@minus,beta2(1:12,:),mean(beta2(1:12,:),1));
                         beta2(13:24,:) = bsxfun(@minus,beta2(13:24,:),mean(beta2(13:24,:),1));
                         seqCorr_w1 = zeros(12,1); seqCorr_w2 = zeros(12,1); seqCorr_acr12 = zeros(12,1); seqCorr_acr21 = zeros(12,1);
                         seqCorr_acr11 = zeros(12,1); seqCorr_acr22 = zeros(12,1);
                         for i=1:12 % each sequence
                             seqCorr_w1(i) = corr(beta1(i,:)',beta1(i+12,:)');
                             seqCorr_w2(i) = corr(beta2(i,:)',beta2(i+12,:)');
                             seqCorr_acr11(i) = corr(beta1(i,:)',beta2(i,:)');
                             seqCorr_acr22(i) = corr(beta1(i+12,:)',beta2(i+12,:)');
                             seqCorr_acr12(i) = corr(beta1(i,:)',beta2(i+12,:)');
                             seqCorr_acr21(i) = corr(beta1(i+12,:)',beta2(i,:)');
                         end
                         % averaged across all sequences
                         D.corr(1,:) = mean(seqCorr_w1);
                         D.corr(2,:) = mean(seqCorr_w2);
                         D.corr(3,:) = mean(seqCorr_acr11);
                         D.corr(4,:) = mean(seqCorr_acr22);
                         D.corr(5,:) = mean(seqCorr_acr12);
                         D.corr(6,:) = mean(seqCorr_acr21);
                         DD = addstruct(DD,D);
                         % trained / untrained
                         S.corr(1,:) = mean(seqCorr_w1(1:6));
                         S.corr(2,:) = mean(seqCorr_w2(1:6));
                         S.corr(3,:) = mean(seqCorr_acr11(1:6));
                         S.corr(4,:) = mean(seqCorr_acr22(1:6));
                         S.corr(5,:) = mean(seqCorr_acr12(1:6));
                         S.corr(6,:) = mean(seqCorr_acr21(1:6));
                         S.corr(7,:) = mean(seqCorr_w1(7:12));
                         S.corr(8,:) = mean(seqCorr_w2(7:12));
                         S.corr(9,:) = mean(seqCorr_acr11(7:12));
                         S.corr(10,:) = mean(seqCorr_acr11(7:12));
                         S.corr(11,:) = mean(seqCorr_acr12(7:12));
                         S.corr(12,:) = mean(seqCorr_acr21(7:12));
                         SS = addstruct(SS,S);
                         % other info
                         I.sess1 = [ss1 ss2 ss1 ss1 ss1 ss1]';
                         I.sess2 = [ss1 ss2 ss2 ss2 ss2 ss2]';
                         I.exe1 = [1 1 1 2 1 2]';
                         I.exe2 = [2 2 1 2 2 1]';
                         I.region = ones(size(I.exe1))*r;
                         I.sn = ones(size(I.exe1))*s;
                         DD = addstruct(DD,I);
                         II = [];
                         II = addstruct(II,I);
                         II = addstruct(II,I);
                         II.seqType = [ones(6,1);ones(6,1)*2];
                         SS = addstruct(SS,II);
                     end
                 end
             end
             fprintf('Done reg %d/%d.\n',find(r==reg),numel(reg));
        end
        save(fullfile(pcmDir,sprintf('FoSEx_correlationModel_pattern_%s',parcelType)),'-struct','DD');
        save(fullfile(pcmDir,sprintf('FoSEx_correlationModel_pattern_seqType_%s',parcelType)),'-struct','SS');
    case 'CORR:acrSess_pattern_plot'
        parcelType = 'Brodmann';
        reg = [1:3,7,8];
        vararginoptions(varargin,{'parcelType'});
        
        D = load(fullfile(pcmDir,sprintf('FoSEx_correlationModel_pattern_%s',parcelType)));
        S = load(fullfile(pcmDir,sprintf('FoSEx_correlationModel_pattern_seqType_%s',parcelType)));
        
        for r=reg
            d = getrow(D,D.region==r);
            s = getrow(S,S.region==r);
            figure(1)
            subplot(2,numel(reg),find(reg==r));
            plt.line([d.sess1>3 d.sess1],d.corr,'subset',d.sess1==d.sess2 & d.exe1==1 & d.exe2==2);
            xlabel('Session'); ylabel('Pattern correlation');
            title(regname_cortex{r})
            subplot(2,numel(reg),find(reg==r)+numel(reg));
            plt.line([s.sess1>3 s.sess1],s.corr,'subset',s.sess1==s.sess2 & s.exe1==1 & s.exe2==2,'split',s.seqType,'style',stySeq,'leg',{'trained','untrained'});
            xlabel('Session'); ylabel('Pattern correlation');
            figure
            subplot(221)
            plt.line(d.sess1,d.corr,'subset',d.sess1==d.sess2 & d.exe1==1 & d.exe2==2);
            title(regname_cortex{r})
            subplot(222)
            plt.line(d.sess1,d.corr,'subset',d.sess2==4 & d.sess1~=4 & d.exe1==d.exe2,'split',d.exe1);
            subplot(223)
            plt.line(s.sess1,s.corr,'subset',s.sess1==s.sess2 & s.exe1==1 & s.exe2==2,'split',s.seqType);
            subplot(224)
            plt.line(s.sess1,s.corr,'subset',s.sess2==4 & s.sess1~=4 & s.exe1==s.exe2,'split',[s.exe1,s.seqType]);
        end
        
    case 'PCM_models_Stats'
        parcelType = 'Brodmann';
        roi=[1:8];
        var='knockOUT'; % knockOUT or logBayes;
        vararginoptions(varargin,{'parcelType','roi','var'});
        sessN=[1:4];
        TT=[];
        % concatenate across all sessions
        for ss=sessN
            T = load(fullfile(pcmDir,sprintf('ModelFamily_FoSEx_Stats_%s_sess%d.mat',parcelType,ss)));     
            TT=addstruct(TT,T);
        end
        % features: 1) first finger, 2) all fingers, 3) transitions, 4)
        % chunk trained, 5) chunk untrained, 6) seqType, 7) trained Seq, 8)
        % untrained seq
        for r=roi
            for ss=1:4
                for i=1:8
                    fprintf('reg-%d session-%d model-%d\n',r,ss,i);
                    fprintf('one-sample 1st exe\n');
                    ttestDirect(TT.(var),[TT.sn],2,'onesample','subset',TT.roi==r & TT.indx==i & TT.sessN==ss & TT.exe==1);
                    fprintf('one-sample 2nd exe\n');
                    ttestDirect(TT.(var),[TT.sn],2,'onesample','subset',TT.roi==r & TT.indx==i & TT.sessN==ss & TT.exe==2);
                    fprintf('paired 1st vs. 2nd\n');
                    ttestDirect(TT.(var),[TT.exe TT.sn],2,'paired','subset',TT.roi==r & TT.indx==i & TT.sessN==ss);
                end 
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
    case 'PLOT_pcm_corr_allSess'        % use this
        reg = [1:8];
        sessN=[1:4];
        modelType='specific';
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','parcelType'});
        
        switch parcelType
            case 'Brodmann'
                lab = regname_cortex;
            case 'striatum'
                lab = {'Caudate','Putamen'};
            case 'striatum_networks'
                lab = {'Network-1','Network-2','Network-3','Network-4',...
                    'Network-5','Network-6','Network-7'};
        end
        
        T=[];
        for t=sessN
            R=load(fullfile(pcmDir,sprintf('PCM_repsup_corr_%s_sess%d_%s.mat',parcelType,t,modelType)));
            %R=load(fullfile(pcmDir,sprintf('PCM_repsup_models_%s_sess%d_%s.mat',modelType,t,parcelType)));
            R.sessN=ones(size(R.SN))*t;
            T=addstruct(T,R);
        end
        
        corrType={'naive','crossval','pcm','mean_naive','mean_crossval'};
        corrVar=[T.r_naive, T.r_crossval, T.r_model2, T.r_naive_meanOnly, T.r_crossval_mean];
        for ct=1:length(corrVar) % three types of correlation 
            figure(ct)
            sub=1;
            for r=reg
                subplot(1,numel(reg),sub)
                plt.line([T.sessN>3 T.sessN],corrVar(:,ct),'subset',T.roi==r,'split',T.seqType,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
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
    case 'PLOT_pcm_models'          %repsup models (Feb 16)
        reg = [1:8];
        sessN=[1:4];
        modelType='specific';
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','parcelType'});
        T=[];
        for t=sessN
            R=load(fullfile(pcmDir,sprintf('PCM_repsup_models_%s_sess%d_%s.mat',modelType,t,parcelType)));
            R.sessN=ones(size(R.SN))*t;
            T=addstruct(T,R);
        end
        model_noRep = bsxfun(@minus,T.bayesEst,T.bayesEst(:,2));
        for r=reg
            figure
            subplot(221)
            barplot(T.sessN,T.bayesEst,'subset',T.roi==r & T.seqType==1);
            xlabel('session');
            ylabel('logBayes');
            title(sprintf('trained %s',regname{r}));
            subplot(223)
            barplot(T.sessN,model_noRep(:,[2:end]),'subset',T.roi==r & T.seqType==1);
            xlabel('session');
            ylabel('logBayes - relative to repModel');
            subplot(222)
            barplot(T.sessN,T.bayesEst,'subset',T.roi==r & T.seqType==2);
            xlabel('session');
            ylabel('logBayes');
            title(sprintf('untrained %s',regname{r}));
            subplot(224)
            barplot(T.sessN,model_noRep(:,[2:end]),'subset',T.roi==r & T.seqType==2);
            xlabel('session');
            ylabel('logBayes - relative to repModel');
        end
    case 'STATS_pcm_models'
        reg = [1:8];
        sessN=[1:4];
        modelType='specific';
        parcelType = 'Brodmann';
        runEffect  = 'fixed';
        vararginoptions(varargin,{'reg','sessN','seqType','modelType','parcelType'});
        T=[];
        for t=sessN
            R=load(fullfile(pcmDir,sprintf('PCM_repsup_models_%s_sess%d_%s.mat',modelType,t,parcelType)));
            R.sessN=ones(size(R.SN))*t;
            T=addstruct(T,R);
        end
        model_noRep = bsxfun(@minus,T.bayesEst,T.bayesEst(:,2));
        % restructure back into a structure
        S.bayesEst  = model_noRep(:);
        S.modelInd  = [ones(size(T.SN));ones(size(T.SN))*2;ones(size(T.SN))*3;ones(size(T.SN))*4;ones(size(T.SN))*5];
        S.SN        = repmat(T.SN,5,1);
        S.sessN     = repmat(T.sessN,5,1);
        S.roi       = repmat(T.roi,5,1);
        S.seqType   = repmat(T.seqType,5,1);
        for r=reg
            fprintf('reg %s - model 3 vs. 2 - trained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==1 & ismember(S.modelInd,[2,3]),'split',S.sessN);
            fprintf('reg %s - model 3 vs. 2 - untrained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==2 & ismember(S.modelInd,[2,3]),'split',S.sessN);
            fprintf('reg %s - model 4 vs. 2 - trained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==1 & ismember(S.modelInd,[2,4]),'split',S.sessN);
            fprintf('reg %s - model 4 vs. 2 - untrained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==2 & ismember(S.modelInd,[2,4]),'split',S.sessN);
            fprintf('reg %s - model 4 vs. 3 - trained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==1 & ismember(S.modelInd,[3,4]),'split',S.sessN);
            fprintf('reg %s - model 4 vs. 3 - untrained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==2 & ismember(S.modelInd,[3,4]),'split',S.sessN);
            fprintf('reg %s - model 5 vs. 2 - trained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==1 & ismember(S.modelInd,[2,5]),'split',S.sessN);
            fprintf('reg %s - model 5 vs. 2 - untrained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==2 & ismember(S.modelInd,[2,5]),'split',S.sessN);
            fprintf('reg %s - model 5 vs. 4 - trained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==1 & ismember(S.modelInd,[4,5]),'split',S.sessN);
            fprintf('reg %s - model 5 vs. 4 - untrained\n',regname{r});
            ttestDirect(S.bayesEst,[S.modelInd S.SN],2,'paired','subset',S.roi==r & S.seqType==2 & ismember(S.modelInd,[4,5]),'split',S.sessN);
            
            fprintf('trained vs. untrained - model 3\n');
            ttestDirect(S.bayesEst,[S.seqType S.SN],2,'paired','subset',S.roi==r & S.modelInd==3,'split',S.sessN);
            fprintf('trained vs. untrained - model 4\n');
            ttestDirect(S.bayesEst,[S.seqType S.SN],2,'paired','subset',S.roi==r & S.modelInd==4,'split',S.sessN);
            fprintf('trained vs. untrained - model 5\n');
            ttestDirect(S.bayesEst,[S.seqType S.SN],2,'paired','subset',S.roi==r & S.modelInd==5,'split',S.sessN);
        end
        
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
    case 'STATS_modelPCM' % OLD
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
        sn=[4:9,11:31];
        sessN=1:4;
        betaChoice='multiPW';
        parcelType = 'Brodmann';
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice'});
        
        AllCorr=[];
        % Make the components for estimation of loadings in G
     %   Gc{1} = [ones(6,6)];   % seqType
     %   Gc{2} = [eye(6)];      % sequence
%         GX=[];
%         for i=1:2
%             GX=[GX Gc{i}(:)];
%         end;
        for ss=sessN
            %T = load(fullfile(regDir,sprintf('betas_FoSEx_sess%d',ss)));
            T = load(fullfile(betaDir,'group',sprintf('betas_FoSEx_%s_sess%d',parcelType,ss)));
            for s=sn
                D = load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM_info'));
                for r=1:max(T.region)
                    for i=1:2 %repetition 1
                        for j=1:2 % repetition 2
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
                                    K = splitHalfCorr(betas1,D1.run,D1.seqNumb,'withinSes');
                                else
                                    % corr across exe
                                    % define condVec, partVec
                                    partVec = [{D1.run} {D2.run}];
                                    condVec = [{D1.seqNumb} {D2.seqNumb}];
                                    data    = [{betas1} {betas2}];
                                    K = splitHalfCorr(data,partVec,condVec,'acrossSes');
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
                              %  K.G          = {G};
                              %  param        = (pinv(GX)*G(:)); % get the two parameters from regression on G
                              %  K.hcov       = param(1); % overall mean
                              %  K.dcov       = param(2); % digit-specific variance
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
        parcelType='Brodmann';
        metric='corr_vox'; % corr_vox; corr_RDM
        vararginoptions(varargin,{'roi','metric','betaChoice'});
       % T = load(fullfile(repSupDir,sprintf('corr_splitHalf_NEW_exe_%s',betaChoice)));
        T = load(fullfile(repSupDir,sprintf('corr_splitHalf_%s_exe_%s',parcelType,betaChoice)));
        
        if strcmp(metric,'corr_RDMAll')
            for r=roi
                R = getrow(T,T.regNum==r);
                figure(1)
                subplot(1,numel(roi),r)
                plt.line([R.sessN>3 R.sessN],R.(metric),'split',R.seqType,'subset',R.exe1~=R.exe2,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
                drawline(0,'dir','horz');
                title(regname{r});
                plt.match('y');
                if r==1
                    ylabel('RDM correlation 1st-2nd');
                else
                    ylabel('');
                end
                figure(2)
                subplot(1,numel(roi),r)
                plt.bar(R.sessN,R.(metric),'subset',R.seqType==1&ismember(R.sessN,[1,4]));
                title(regname{r});
                plt.match('y');
                ylabel('');
            end
        else
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
        end
        % split into within and across
    case 'CORR_splithalf_corrected' % voxel or RDM correlation
        exe = [1 2];
        betaChoice='multiPW';
        sessN=[1:4];
        sn=[4:9,11:31];
        parcelType='Brodmann';
        vararginoptions(varargin,{'betaChoice','parcelType'});
        CC=[];
        T = load(fullfile(repSupDir,sprintf('corr_splitHalf_%s_exe_%s',parcelType,betaChoice)));
        sn=unique(T.sn);
        for ss=sessN
            for s=1:length(sn)
                for r=1:max(T.regNum)
                    D = getrow(T,T.sn==sn(s) & T.regNum==r & T.sessN==ss);
                    for st=1:2 % seqType
                        for i=1:length(exe)
                            for j=i:length(exe)
                                t1=getrow(D,D.seqType==st & D.exe1==i & D.exe2==i);
                                t2=getrow(D,D.seqType==st & D.exe1==j & D.exe2==j);
                                tcross=getrow(D,D.seqType==st & D.exe1==i & D.exe2==j);
                                if size(tcross.corr_vox,1)~=0 % if there is data for both sessions
                                    C.corr_vox = tcross.corr_vox/sqrt(t1.corr_vox*t2.corr_vox);
                                    C.corr_RDM = tcross.corr_RDM/sqrt(t1.corr_RDM*t2.corr_RDM);
                                    if ~isreal(C.corr_vox)
                                        C.corr_vox=0; % or NaN
                                    end
                                    if ~isreal(C.corr_RDM)
                                        C.corr_RDM=0; % or NaN
                                    end
                                    % other info
                                    C.sn = sn(s);
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
        save(fullfile(repSupDir,sprintf('corr_splitHalf_%s_exe_%s_corrected',parcelType,betaChoice)),'-struct','CC');
    case 'PLOT_splithalf_corr_corrected' % voxel or RDM correlation
        roi=[1:8];
        betaChoice='multiPW';
        parcelType='Brodmann';
        metric='corr_vox'; % corr_vox or corr_RDM
        vararginoptions(varargin,{'roi','metric','betaChoice'});
        
        T = load(fullfile(repSupDir,sprintf('corr_splitHalf_%s_exe_%s_corrected',parcelType,betaChoice)));
        
        figure
        for r=roi
            subplot(1,length(roi),r);
            R = getrow(T,T.regNum==r);
            plt.line([R.sessN>3 R.sessN],R.(metric),'split',R.seqType,'subset',R.exe1~=R.exe2,'style',stySeq,'leg',{'trained','untrained'},'leglocation','northeast');
            title(sprintf('%s - across exe',regname{r}));  
        end  
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
                        diff_psc=t1.psc-t2.psc;
                        
                        T2.psc=[t1.psc;t2.psc;diff_psc];
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
                        diff_psc=t2.psc/t1.psc;
                        
                        T2.psc=[t1.psc;t2.psc;diff_psc];
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
     
    case 'RDM_plot'
        % plot the average RDM per region - exe1/2
        sessN=4;
        reg=[2,3,7];
        parcelType = 'Brodmann';
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'sn','sessN','reg'});
        
        H = eye(12) - 1/12; 
        for ss=sessN
            T = load(fullfile(betaDir,'group',sprintf('stats_FoSEx_%s_%s_sess%d',parcelType,betaChoice,ss)));
            idx=1;
            for e=1:2
                for r=reg
                    t = getrow(T,T.region==r&T.FoSEx==e);
                    RDM = rsa_squareRDM(nanmean(t.RDM,1));
                    %G = rsa_squareIPM(nanmean(t.IPM,1));
                    G = -0.5*H*RDM*H';
                    figure(1)
                    subplot(2,numel(reg),idx)
                    imagesc(RDM);
                    figure(2)
                    subplot(2,numel(reg),idx)
                    imagesc(G);
                    idx=idx+1;
                end
            end
        end; % session
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
                     diff_m = K1.(metric)-K2.(metric);
                     subplot(2,numel(sessN),indx);
                     imagesc(rsa_squareRDM(mean(diff_m,1)));
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
    
    case 'job2'
        sml1_imana_repsup('PCM_constructModelFamily','reg',[1:8],'naturalStats',1);
        sml1_imana_repsup('PCM_constructModelFamily','reg',[1:8]);
    case 'job3'
        sml1_imana_dist('SEARCH_dist_runLDC','sn',[13:25],'sessN',3);
    case 'run_job'
        sml1_imana_repsup('SURF_wb:map_RS_absolute');
        sml1_imana_repsup('SURF_wb:smooth_individ');
        sml1_imana_repsup('SURF_wb:map_group');
        sml1_imana_repsup('SURF_wb:cSPM_group');
        
    case 'HOUSEKEEPING:renameSPM' 	% rename SPM directories - to do - new FUNC
        sn = [5:9,11:31];
        sessN = 4;
        vararginoptions(varargin,{'sn','sessN'});
        fprintf('Renaming SPMs:\n');
        for ss=sessN
            for s=sn
                newGLMDir   = fullfile(glmFoSExDir{ss},subj_name{s});
                cd(newGLMDir);
                load SPM;
                newRawDir   = fullfile(imagingDir,subj_name{s});
                SPM = spmj_move_rawdata(SPM,newRawDir);
                SPM.swd = fullfile(newGLMDir);
                %save(fullfile(newGLMDir,'SPM.mat'),'SPM','-v7.3');
                fprintf('Done sess-%d sn-%d\n',ss,s);
            end
        end
        varargout{1} = SPM;
        
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
function M = pcm_corrModel_indSeq_sharedSeqType % OLD
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
    M{4}.Ac(:,13,4)   = [ones(6,1);zeros(6,1)]; % Exe1 mean
    M{4}.Ac(:,14,5)   = [zeros(6,1);ones(6,1)]; % Exe2 mean
    
    % Model 5: Execution + sequence specific + PERFECT correlation in sequences
    M{5}.type         = 'feature';
    M{5}.numGparams   = 4;
    M{5}.name         = 'RepSup+Seq+PerfectCorr';
    M{5}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{5}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    M{5}.Ac(:,3:8,3)  = [A;zeros(6)];       % Unique sess1 sequence patterns
    M{5}.Ac(:,3:8,4)  = [zeros(6);A];        % Identical sess2 sequence pattterns

    
end
function M = pcm_repsupModel_specific % TO USE
    
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
    
   % Model 5: Execution + sequence specific + PERFECT correlation in session
    M{5}.type         = 'feature';
    M{5}.numGparams   = 14;
    M{5}.name         = 'RepSup+Seq+PerfectCorr';
    M{5}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{5}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{5}.Ac(:,3:8,2+i)   = [A;zeros(6)];       % Unique exe1 sequence patterns
        M{5}.Ac(:,3:8,8+i)   = [zeros(6);A];       % Same exe2 sequence pattterns
    end;

end
function M = pcm_defineSequenceModels_G(Seq,sn,NatStat)
% define sequence models completely
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
    M{1}.Gc(:,:)=zeros(nSeq,5); % first finger
    for u = 1:nSeq    % 12 seq
        firstfing = Seq(:,1);
        M{1}.Gc(u,firstfing(u)) = 1;
    end
    M{1}.Gc = M{1}.Gc*NatStat*M{1}.Gc';
    M{1}.name = 'FirstFing';
    M{1}.numGparams = 1;
    M{1}.type = 'component';
    % ---------------------------------------
    % Model2: All fingers
    M{2}.Gc(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{2}.Gc(u,j) = length(placenumb);
        end
    end
    M{2}.Gc = M{2}.Gc*NatStat*M{2}.Gc'; 
    M{2}.name = 'AllFing';
    M{2}.numGparams = 1;
    M{2}.type = 'component';
    % ---------------------------------------
%     % Model3: Model with trained & untrained labels - 
%     M{3}.Gc(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
%     M{3}.Gc(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
%     M{3}.Gc         = M{3}.Gc*M{3}.Gc';
%     M{3}.name = 'SeqType';
%     M{3}.numGparams = 1;
%     M{3}.type = 'component';
    M{3}.Gc = ones(12);
    M{3}.name = 'Mean_activity';
    M{3}.numGparams = 1;
    M{3}.type = 'component';
    % --------------------------------------
    % Model4: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{4}.Gc(:,1:6)  = [A;zeros(6)];       % Trained sequence patterns
    M{4}.Gc         = M{4}.Gc*M{4}.Gc';
    M{4}.name = 'TrainedSeq';
    M{4}.numGparams = 1;
    M{4}.type = 'component';
    % --------------------------------------
    % Model5: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}.Gc(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
    M{5}.Gc         = M{5}.Gc*M{5}.Gc';
    M{5}.name = 'UntrainedSeq';
    M{5}.numGparams = 1;
    M{5}.type = 'component';
 
    % Model 6: mean activity
%     M{6}.Gc = ones(12);
%     M{6}.name = 'MeanActivity';
%     M{6}.numGparams = 1;
%     M{6}.type = 'component';
end
function M = pcm_defineSequenceModels_G_simple(Seq,sn,NatStat)
% define sequence models completely
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
    M{1}.Gc(:,:)=zeros(nSeq,5); % first finger
    for u = 1:nSeq    % 12 seq
        firstfing = Seq(:,1);
        M{1}.Gc(u,firstfing(u)) = 1;
    end
    M{1}.Gc = M{1}.Gc*NatStat*M{1}.Gc';
    M{1}.name = 'FirstFing';
    M{1}.numGparams = 1;
    M{1}.type = 'component';
    % ---------------------------------------
    % Model2: All fingers
    M{2}.Gc(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{2}.Gc(u,j) = length(placenumb);
        end
    end
    M{2}.Gc = M{2}.Gc*NatStat*M{2}.Gc'; 
    M{2}.name = 'AllFing';
    M{2}.numGparams = 1;
    M{2}.type = 'component';
    % ---------------------------------------
%     % Model3: Sequence model 
    M{3}.Gc    = eye(12);
    M{3}.name = 'Sequence';
    M{3}.numGparams = 1;
    M{3}.type = 'component';

    % ------
    % Model4: Intercept
    M{4}.Gc = ones(12);
    M{4}.name = 'Mean_activity';
    M{4}.numGparams = 1;
    M{4}.type = 'component';
end
function M = pcm_defineSequenceModels_G_seqType(Seq,sn,NatStat,st)
% define sequence models completely
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
    M{1}.Gc(:,:)=zeros(nSeq/2,5); % first finger
    for u = 1:6    % 12 seq
        if st==1
            firstfing = Seq(1:6,1);
        else
            firstfing = Seq(7:12,1);
        end
        M{1}.Gc(u,firstfing(u)) = 1;
    end
    M{1}.Gc = M{1}.Gc*NatStat*M{1}.Gc';
    M{1}.name = 'FirstFing';
    M{1}.numGparams = 1;
    M{1}.type = 'component';
    % ---------------------------------------
    % Model2: All fingers
    M{2}.Gc(:,:)=zeros(nSeq/2,5); % all fingers
    if st==1
        seqL = 1:6;
    else
        seqL = 7:12;
    end
    for u = 1:6
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(seqL(u),:)==j);
            M{2}.Gc(u,j) = length(placenumb);
        end
    end
    M{2}.Gc = M{2}.Gc*NatStat*M{2}.Gc'; 
    M{2}.name = 'AllFing';
    M{2}.numGparams = 1;
    M{2}.type = 'component';
    % ---------------------------------------
    % Model3: Model for each specific  sequence
    M{3}.Gc = eye(6);   
    M{3}.Gc = M{3}.Gc*M{3}.Gc';
    M{3}.name = 'Seq';
    M{3}.numGparams = 1;
    M{3}.type = 'component';
    % --------------------------------------
    % Model4: Mean activity
    M{4}.Gc = ones(6);
    M{4}.name = 'Mean_activity';
    M{4}.numGparams = 1;
    M{4}.type = 'component';

end


function M = pcm_defineSequenceModels_fixed_noChunks_natStats(Seq,sn,NatStat)
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
function M = pcm_defineSequenceModels_new(Seq,sn,NatStat)
% written April 17th 2020
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
    % Model3: Trained / Untrained sequences
    %M{3} = [ones(6) zeros(6); zeros(6) ones(6)];
    M{3}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{3}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{3}(:,3:8)    = [A;zeros(6)];       % Trained sequence patterns
     % each specific UNTRAINED sequence
    M{3}(:,9:14)   = [zeros(6);A];
end

function M = pcm_defineSequenceModel_fixed_splitFirstFing(Seq,sn,NatStat)
% written November 29th 2019
%  No fixed / run component added
% specific sequence + overall component modelled together
    % --------------------------------------
    % Model1: First finger model for trained
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    nSeq=size(Seq,1);
    M{1}(:,:)=zeros(nSeq,5); % first finger
    trainedSq = [6,5,4,3,2,1];
    for u = 1:6    % trained seq
    %for u = trainedSq    % trained seq
        firstfing = Seq(:,1);
        M{1}(u,firstfing(u)) = 1;
        %M{1}(u,firstfing(trainedSq(u))) = 1;
    end
    M{1} = M{1}*NatStat*M{1}';
    % --------------------------------------
    % Model2: First finger model for untrained
    M{2}(:,:)=zeros(nSeq,5); % first finger
    for u = 7:12 % untrained seq
        firstfing = Seq(:,1);
        M{2}(u,firstfing(u)) = 1;
    end
    M{2} = M{2}*NatStat*M{2}';
    % ---------------------------------------
    % Model3: All fingers
    M{3}(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{3}(u,j) = length(placenumb);
        end
    end
    M{3} = M{3}*NatStat*M{3}'; 
    % ---------------------------------------
    % Model4: Model with trained & untrained labels - 
    M{4}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{4}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model5: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model6: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{6}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
end 
function M = pcm_defineSequenceModels_Fing(Seq,sn,NatStat)
% written February 12 2020
% splitting finger components
    % --------------------------------------
    % Model1: First finger trained model
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    nSeq=size(Seq,1);
    M{1}(:,:)=zeros(nSeq,5); % first finger trained
    for u = 1:6    % trained
        firstfing = Seq(1:6,1);
        M{1}(u,firstfing(u)) = 1;
    end
    M{1} = M{1}*NatStat*M{1}';
    % ---------------------------------------
    % Model First finger untrained
    M{2}(:,:)=zeros(nSeq,5); % first finger trained
    for u = 1:6    %untrained
        firstfing = Seq(7:12,1);
        M{2}(u+6,firstfing(u)) = 1;
    end
    M{2} = M{2}*NatStat*M{2}';
    % ---------------------------------------
    % Model first finger trained - untrained correspondence
%      corrF = M{2}(7:12,7:12);
%      M{3}=zeros(12);
%      M{3}(1:6,7:12) = fliplr(corrF);
    % ---------------------------------------
    % Model: All fingers
    M{3}(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{3}(u,j) = length(placenumb);
        end
    end
    M{3} = M{3}*NatStat*M{3}'; 
    % ---------------------------------------
    % Model5: Model with trained & untrained labels - 
    M{4}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{4}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model6: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model7: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{6}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
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
function M = pcm_defineFingModels_natStats(Seq,sn,NatStat)
% just all fingers and first finger
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
end
function M = pcm_defineFingModels(Seq,sn)
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
end
function M = pcm_defineSequenceModels_test(Seq,sn)
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    nSeq=size(Seq,1);
    % ---------------------------------------
    firstfing = Seq(:,1);
    % Model1: First finger for trained
    M{1}(:,:)=zeros(nSeq,5); % first finger
    for u = 1:6
        M{1}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model2: First finger for untrained
    M{2}(:,:)=zeros(nSeq,5); % first finger
    for u = 7:12
        M{2}(u,firstfing(u)) = 1;
    end
    % ---------------------------------------
    % Model3: All fingers
    M{3}(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{3}(u,j) = length(placenumb);
        end
    end
    % ---------------------------------------
    % Model4: Model with trained & untrained labels - 
    M{4}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{4}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model5: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model6: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{6}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
end
function M = pcm_defineSequenceModels_test_natStats(Seq,sn,NatStat)
    if rem(sn,2)==0 % for group 2 - switch around sequences 1-6, 7-12
        SeqNew = zeros(12,9);
        SeqNew([1:6],:)=Seq([7:12],:);
        SeqNew([7:12],:)=Seq([1:6],:);
        Seq=SeqNew;
    end
    nSeq=size(Seq,1);
    % ---------------------------------------
    firstfing = Seq(:,1);
    % Model1: First finger for trained
    M{1}(:,:)=zeros(nSeq,5); % first finger
    for u = 1:6
        M{1}(u,firstfing(u)) = 1;
    end
    %M{1} = M{1}*NatStat*M{1}';
    % ---------------------------------------
    % Model2: First finger for untrained
    M{2}(:,:)=zeros(nSeq,5); % first finger
    for u = 7:12
        M{2}(u,firstfing(u)) = 1;
    end
    %M{2} = M{2}*NatStat*M{2}';
    % ---------------------------------------
    % Model3: All fingers
    M{3}(:,:)=zeros(nSeq,5); % all fingers
    for u = 1:nSeq
        for j = 1:5                 % 5 fingers
            placenumb = find(Seq(u,:)==j);
            M{3}(u,j) = length(placenumb);
        end
    end
    M{3} = M{3}*NatStat*M{3}';
    % ---------------------------------------
    % Model4: Model with trained & untrained labels - 
    M{4}(:,1)    = [ones(6,1);zeros(6,1)]; % Common component to trained
    M{4}(:,2)    = [zeros(6,1);ones(6,1)]; % Common component to untrained
    % --------------------------------------
    % Model5: Model for each specific TRAINED sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{5}(:,1:6)    = [A;zeros(6)];       % Trained sequence patterns
    % --------------------------------------
    % Model6: Model for each specific UNTRAINED sequence    
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{6}(:,1:6)  = [zeros(6);A];       % Untrained sequence patterns
 
end
function M = pcm_defineSequenceModels_seqType_natStats(Seq,sn,NatStat,seqType)
% written January 29th 2020
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
    if seqType == 1 % trained
        Seq = Seq(1:6,:);
    else % untrained
        Seq = Seq(7:12,:);
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
    % Model3: Model for each specific sequence
%      A=zeros(6);
%      for i=1:6
%         A(i,i)=1;
%      end;
%      M{3} = A;    
end
function M = pcm_defineSequenceModels_seqType(Seq,sn,seqType)
% written January 29th 2020
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
    if seqType == 1 % trained
        Seq = Seq(1:6,:);
    else % untrained
        Seq = Seq(7:12,:);
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
    % --------------------------------------
    % Model3: Model for each specific sequence
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    M{3} = A;    
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
function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm,S)
    % --------------------------------------
    % Crossvalidated model comparision:
    %[T,theta_hat,G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm,'S',S);
    [T,theta_hat,T.G_pred,theta0] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    [Tcross,thetaCr,Tcross.G_predCv] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
 %   [Tcross2,thetaCr2] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm',algorithm);
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);

    T.cross_likelihood = Tcross.likelihood;
    T.G_pred_cross = Tcross.G_predCv;
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
    Gcv(:,:,p) = crossval_estG(data,Z,condVec{p});
    b = pinv(Z)*data;
   % G2 = cov(b');
    C.r_crossval(p,1) = calcCorr(pcm_makePD(Gcv(:,:,p))); % this is the correct one
    % calculate correlation of mean pattern
    % make even and odd mean patterns for 1st / 2nd
    tmp11 = mean(meanP(1,:,[1:2:8]),3);
    tmp12 = mean(meanP(1,:,[2:2:8]),3);
    tmp21 = mean(meanP(2,:,[1:2:8]),3);
    tmp22 = mean(meanP(2,:,[2:2:8]),3);
    cor1 = corr(tmp11',tmp22');
    cor2 = corr(tmp12',tmp21');
    C.r_crossval_mean(p,1) = mean([cor1 cor2]); 
    % calculate correlation of Gs for 1st / 2nd
    RDM = rsa_squareRDM(rsa.distanceLDC(data,partVec{p},condVec{p}));
    Exe{1}=RDM(1:6,1:6);
    Exe{2}=RDM(7:12,7:12);
  %  C.RDM_corrRSA(p,1)  = rsa_calcCorrRDMs(Exe);
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
nElem = size(G,1)/2;
v1 = d0(1:nElem)';    % Variances first exe
v2 = d0(nElem+1:end)';   % Variances 2nd exe
cv=diag(G,nElem);     % Covariance
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
        G = crossval_estG(Y,X,ses);
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
