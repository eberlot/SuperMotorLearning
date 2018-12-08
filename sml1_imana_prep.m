function varargout=sml1_imana_prep(what,varargin)

% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
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

% total - per subject (the total in the end will always be 40)
numruns           = [40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40];

sess = [repmat(1,1,10),repmat(2,1,10),repmat(3,1,10),repmat(4,1,10)];   % all sessions

sess_sn = [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4];    % per subject

run_task   = [1:3 5:7 9:10;
              11:13 15:17 19:20;
              21:23 25:27 29:30;
              31:33 35:37 39:40];    % task
run_loc    = [4 8;
              14 18;
              24 28;
              34 38];             % localizer

run_num{1} = [1:10];
run_num{2} = [11:20];
run_num{3} = [21:30];
run_num{4} = [31:40];

runs{1}    = {'_01','_02','_03','_04','_05','_06','_07','_08','_09','_10'};
runs{2}    = {'_11','_12','_13','_14','_15','_16','_17','_18','_19','_20'};
runs{3}    = {'_21','_22','_23','_24','_25','_26','_27','_28','_29','_30'};
runs{4}    = {'_31','_32','_33','_34','_35','_36','_37','_38','_39','_40'};  
          
TRlength   = 1000;      % in ms
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
% The variables in this section must be updated for every new subject.
%       DicomName  :  first portion of the raw dicom filename
%       NiiRawName :  first protion of the nitfi filename (get after 'PREP_4d_nifti')
%       fscanNum   :  series # for corresponding functional runs. Enter in run order
%       anatNum    :  series # for anatomical scans (~208 or so imgs/series)
%       loc_AC     :  location of the anterior commissure. 
%
% The values of loc_AC should be acquired manually prior to the preprocessing
%   Step 1: get .nii file of anatomical data by running "spmj_tar2nii(TarFileName,NiiFileName)"
%   Step 2: open .nii file with MRIcron and manually find AC and read the xyz coordinate values
%           (note: there values are not [0 0 0] in the MNI coordinate)
%   Step 3: set those values into loc_AC (subtract from zero)

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
              's21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'};  


% different sessions denoted with {}
DicomName{1}  = {'2017_07_18_S01.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_07_25_S02.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_07_25_S03.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_08_S04.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_09_05_S05.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_09_19_S06.MR.Diedrichsen_LongSeqLearn',...
                 '2017_09_12_S07.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_10_03_S08.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_17_S09.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_19_S10.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_01_S11.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_01_S12.MR.Diedrichsen_LongSeqLearn',...
                 '2018_01_16_S13.MR.Diedrichsen_LongSeqLearn',...
                 '2018_01_19_S14.MR.Diedrichsen_LongSeqLearn',...
                 '2018_01_25_S15.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_06_S16.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_28_S17.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_01_S18.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_13_S19.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_23_S20.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_06_S21.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_05_S22.MR.Diedrichsen_LongSeqLearn',...
                 '2018_05_15_S23.MR.Diedrichsen_LongSeqLearn',...
                 '2018_05_24_S24.MR.Diedrichsen_LongSeqLearn',...
                 '2018_05_25_S025.MR.Diedrichsen_LongSeqLearn',...
                 '2018_07_06_S26.MR.Diedrichsen_LongSeqLearn',...
                 '2018_09_11_S27.MR.Diedrichsen_LongSeqLearn',...
                 '2018_07_16_S28.MR.Diedrichsen_LongSeqLearn',...
                 '2018_09_20_S29.MR.Diedrichsen_LongSeqLearn',...
                 '2018_09_20_S30.MR.Diedrichsen_LongSeqLearn',...
                 '2018_09_19_S31.MR.Diedrichsen_LongSeqLearn'};             
DicomName{2}  = {'2017_07_25_S01.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_01_S02.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_01_S03.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_15_S04.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_09_12_S05.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_09_29_S06.MR.Diedrichsen_LongSeqLearn',...
                 '2017_09_19_S07.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_10_S08.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_24_S09.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_03_S10.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_10_S11.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_07_S12.MR.Diedrichsen_LongSeqLearn',...
                 '2018_01_23_S13.MR.Diedrichsen_LongSeqLearn',...
                 '2018_01_29_S14.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_05_S15.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_13_S16.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_05_S17.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_08_S18.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_21_S19.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_29_S20.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_12_S21.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_10_S22.MR.Diedrichsen_LongSeqLearn',...
                 '2018_05_22_S23.MR.Diedrichsen_LongSeqLearn',...
                 '2018_06_01_S24.MR.Diedrichsen_LongSeqLearn',...
                 '2018_05_31_S25.MR.Diedrichsen_LongSeqLearn',...
                 '2018_07_13_S26.MR.Diedrichsen_LongSeqLearn',...
                 '2018_09_18_S27.MR.Diedrichsen_LongSeqLearn',...
                 '2018_07_23_S28.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_04_S29.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_01_S30.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_01_S31.MR.Diedrichsen_LongSeqLearn'};
DicomName{3}  = {'2017_08_15_S01.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_22_S02.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_23_S03.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_09_06_S04.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_10_03_S05.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_17_S06.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_10_S07.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_31_S08.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_13_S09.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_23_S10.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_28_S11.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_28_S12.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_05_S13.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_12_S14.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_13_S15.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_06_S16.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_27_S17.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_28_S18.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_09_S19.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_19_S20.MR.Diedrichsen_LongSeqLearn',...
                 '2018_05_03_S21.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_26_S22.MR.Diedrichsen_LongSeqLearn',...
                 '2018_06_14_S23.MR.Diedrichsen_LongSeqLearn',...
                 '2018_06_21_S24.MR.Diedrichsen_LongSeqLearn',...
                 '2018_06_19_S25.MR.Diedrichsen_LongSeqLearn',...
                 '2018_07_26_S26.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_09_S27.MR.Diedrichsen_LongSeqLearn',...
                 '2018_08_16_S28.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_16_S29.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_29_S30.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_15_S31.MR.Diedrichsen_LongSeqLearn'};             
DicomName{4}  = {'2017_08_16_S01.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_23_S02.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_08_25_S03.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_09_12_S04.MR.DIEDRICHSEN_LONGSEQLEARN',...
                 '2017_10_04_S05.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_18_S06.MR.Diedrichsen_LongSeqLearn',...
                 '2017_10_11_S07.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_01_S08.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_15_S09.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_24_S10.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_29_S11.MR.Diedrichsen_LongSeqLearn',...
                 '2017_11_29_S12.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_06_S13.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_15_S14.MR.Diedrichsen_LongSeqLearn',...
                 '2018_02_15_S15.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_07_S16.MR.Diedrichsen_LongSeqLearn',...
                 '2018_03_28_S17.MR.Diedrichsen_LongSeqLearn'....
                 '2018_03_29_S18.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_10_S19.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_20_S20.MR.Diedrichsen_LongSeqLearn',...
                 '2018_05_04_S21.MR.Diedrichsen_LongSeqLearn',...
                 '2018_04_27_S22.MR.Diedrichsen_LongSeqLearn',...
                 '2018_06_15_S23.MR.Diedrichsen_LongSeqLearn',...
                 '2018_06_22_S24.MR.Diedrichsen_LongSeqLearn',...
                 '2018_06_22_S25.MR.Diedrichsen_LongSeqLearn',...
                 '2018_08_01_S26.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_10_S27.MR.Diedrichsen_LongSeqLearn',...
                 '2018_08_17_S28.MR.Diedrichsen_LongSeqLearn',...
                 '2018_11_08_S29.MR.Diedrichsen_LongSeqLearn',...
                 '2018_10_31_S30.MR.Diedrichsen_LongSeqLearn',...
                 '2018_11_09_S31.MR.Diedrichsen_LongSeqLearn'};


NiiRawName{1} = {'170718114530DST131221107523418932',...
                 '170725105405DST131221107523418932',...
                 '170725090958DST131221107523418932',...
                 '170808113354DST131221107523418932',...
                 '170905134655DST131221107523418932',...
                 '170919130243DST131221107523418932',...
                 '170912103955DST131221107523418932',...
                 '171003152437DST131221107523418932',...
                 '171017144222DST131221107523418932',...
                 '171019144204DST131221107523418932',...
                 '171101104440DST131221107523418932',...
                 '171101123744DST131221107523418932',...
                 '2018_01_16_S13',...
                 '2018_01_19_S14',...
                 '2018_01_25_S15',...
                 '2018_02_06_S16',...
                 '2018_02_28_S17',...
                 '2018_03_01_S18',...
                 '2018_03_13_S19',...
                 '2018_03_23_S20',...
                 '2018_04_06_S21',...
                 '2018_04_05_S22',...
                 '2018_05_15_S23',...
                 '2018_05_24_S24',...
                 '180525120459DST131221107523418932',...
                 '2018_07_06_S26',...
                 '2018_09_11_S27',...
                 '2018_07_16_S28',...
                 '2018_09_20_S29',...
                 '2018_09_20_S30',...
                 '2018_09_19_S31'};
NiiRawName{2} = {'170725124629DST131221107523418932',...
                 '170801131823DST131221107523418932',...
                 '170801113315DST131221107523418932',...
                 '170815122423DST131221107523418932',...
                 '170912143950DST131221107523418932',...
                 '170929141715DST131221107523418932',...
                 '170919112438DST131221107523418932',...
                 '171010084845DST131221107523418932',...
                 '171024144829DST131221107523418932',...
                 'S10',...
                 '171110104631STD131221107523418932',...
                 '171107103831STD131221107523418932',...
                 '2018_01_23_S13',...
                 '2018_01_29_S14',...
                 '2018_02_05_S15',...
                 '2018_02_13_S16',...
                 '2018_03_05_S17',...
                 '2018_03_08_S18',...
                 '2018_03_21_S19',...
                 '2018_03_29_S20',...
                 '2018_04_12_S21',...
                 '2018_04_10_S22',...
                 '2018_05_22_S23',...
                 '2018_06_01_S24',...
                 '2018_05_31_S25',...
                 '2018_07_13_S26',...
                 '2018_09_18_S27',...
                 '2018_07_23_S28',...
                 '2018_10_04_S29',...
                 '2018_10_01_S30',...
                 '2018_10_01_S31'};
NiiRawName{3}  = {'170815144204DST131221107523418932',...
                  '170822110530DST131221107523418932',...
                  '170823110924DST131221107523418932',...
                  '170906090741DST131221107523418932',...
                  '171003133651DST131221107523418932',...
                  '171017120505DST131221107523418932',...
                  '171010145029DST131221107523418932',...
                  '171031090438DST131221107523418932',...
                  '171113145300STD131221107523418932',...
                  's10',...
                  '171128144155STD131221107523418932',...
                  '171128101440STD131221107523418932',...
                  '2018_02_05_S13',...
                  '2018_02_12_S14',...
                  '2018_02_13_S15',...
                  '2018_03_06_S16',...
                  '2018_03_27_S17',...
                  '2018_03_28_S18',...
                  '2018_04_09_S19',...
                  '2018_04_19_S20',...
                  '2018_05_03_S21',...
                  '2018_04_26_S22',...
                  '2018_06_14_S23',...
                  '2018_06_21_S24',...
                  '2018_06_19_S25',...
                  '2018_07_26_S26',...
                  '2018_10_09_S27',...
                  '2018_08_16_S28',...
                  '2018_10_16_S29',...
                  '2018_10_29_S30',...
                  '2018_10_15_S31'};  
NiiRawName{4}  = {'170816090254DST131221107523418932',...
                  '170823091559DST131221107523418932',...
                  '170825085040DST131221107523418932',...
                  's04',...
                  '171004134839DST131221107523418932',...
                  '171018114933DST131221107523418932',...
                  '171011103434DST131221107523418932',...
                  '171101085911DST131221107523418932',...
                  '171115143146STD131221107523418932',...
                  '171124120324STD131221107523418932',...
                  '171129100600STD131221107523418932',...
                  '171129133110STD131221107523418932',...
                  '2018_02_06_S13',...
                  '2018_02_15_S14',...
                  '2018_02_15_S15',...
                  '2018_03_07_S16',...
                  '2018_03_28_S17',...
                  '2018_03_29_S18',...
                  '2018_04_10_S19',...
                  '2018_04_20_S20',...
                  '2018_05_04_S21',...
                  '2018_04_27_S22',...
                  '2018_06_15_S23',...
                  '2018_06_22_S24',...
                  '2018_06_22_S25',...
                  '2018_08_01_S26',...
                  '2018_10_10_S27',...
                  '2018_08_17_S28',...
                  '2018_11_08_S29',...
                  '2018_10_31_S30',...
                  '2018_11_09_S31'};

fscanNum{1}   = {[16 18 20 22 24 26 28 30 32 34],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [17 19 21 23 25 27 29 31 33 35],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [19 21 23 25 27 29 31 33 35 37],...
                 [10 12 14 16 18 20 22 24 26 28],...
                 [19 21 23 25 27 29 31 33 35 37],...
                 [17 19 23 25 27 29 31 33 35 37],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [17 19 21 23 25 27 29 31 33 35],...
                 [17 19 21 23 25 27 29 31 33 35],...
                 [18 20 22 24 26 28 30 32 34 36],...
                 [18 20 22 24 26 30 32 34 36 38],...
                 [17 19 21 23 25 27 29 31 33 35],...
                 [18 20 22 24 26 28 30 32 34 36],...
                 [23 25 27 29 31 33 35 37 39 41],...
                 [17 19 21 23 25 27 29 31 33 35],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [17 19 21 23 25 29 31 33 35 37],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [17 19 21 23 25 27 29 31 33 35],...
                 [16 18 20 22 24 26 28 30 34 36],...
                 [16 18 20 22 28 30 32 34 36 38],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [18 20 22 24 26 28 32 34 36 38],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [16 18 20 22 24 26 28 30 32 34],...
                 [18 20 22 24 26 28 30 32 34 38],...
                 [20 22 24 26 28 30 32 34 36 38]};             
fscanNum{2}   = {[11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 23 27 31 25 33 35],...
                 [14 16 18 28 22 24 26 34 30 32],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [12 14 16 18 20 22 24 28 30 32],...
                 [16 18 20 24 26 28 30 32 34 36],...
                 [16 18 20 22 24 26 28 30 32 36],...
                 [17 19 24 26 28 30 32 34 36 38],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [11 15 17 19 21 23 25 27 29 31],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [13 15 17 19 21 25 27 29 31 33],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [11 13 15 17 19 21 25 27 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [11 13 15 17 19 21 25 27 29 31],...
                 [15 17 19 21 23 25 27 29 31 33],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 21 23 25 27 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 23 25 27 29 31],...
                 [11 15 17 19 21 23 25 27 29 31],...
                 [16 18 20 22 24 26 28 30 32 34]};    
fscanNum{3}   = {[11 13 15 17 19 21 23 25 27 29],...
                 [11 13 19 17 21 23 27 25 29 31],... 
                 [11 13 15 17 19 21 23 25 27 29],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [15 17 19 21 23 25 27 29 31 33],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [13 15 17 19 21 23 27 29 31 35],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [13 15 17 19 21 25 27 31 33 35],...
                 [11 13 15 17 23 25 27 29 31 33],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 31 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [12 14 18 20 22 24 26 28 30 32],...
                 [11 13 15 17 19 21 25 27 29 31],...
                 [11 15 17 19 21 23 25 27 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [13 15 17 19 21 23 25 27 29 33],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [11 13 15 17 19 21 23 25 27 29]};  
fscanNum{4}   = {[11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 21 23 27 25 29 31],...  
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 29 31],...
                 [14 16 18 20 22 24 26 28 30 32],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 19 21 23 25 27 29 31 33 35],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 21 23 25 27 29 31 33],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [12 14 16 18 20 22 24 26 28 30],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 27 25 29 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [14 16 18 20 22 24 26 28 30 32],...
                 [11 13 15 17 19 21 23 25 27 31],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [13 15 17 19 21 23 25 27 29 31],...
                 [14 16 18 20 22 24 26 28 30 32],...
                 [14 16 18 20 22 24 26 28 30 32],...
                 [13 15 17 19 23 25 27 33 37 39],...
                 [11 13 15 17 19 21 23 25 27 29],...
                 [15 17 23 25 27 29 31 33 39 41]};    

fieldNum{1}   = {[35,36],...
                 [35,36],...
                 [36,37],...
                 [35,36],...
                 [38,39],...
                 [29,30],...
                 [38,39],...
                 [38,39],...
                 [35,36],...
                 [35,36],...
                 [35,36],...
                 [36,37],...
                 [36,37],...
                 [37,38],...
                 [39,40],...
                 [36,37],...
                 [37,38],...
                 [42,43],...
                 [36,37],...
                 [35,36],...
                 [38,39],...
                 [35,36],...
                 [36,37],...
                 [37,38],...
                 [39,40],...
                 [35,36],...
                 [39,40],...
                 [35,36],...
                 [35,36],...
                 [39,40],...
                 [39,40]};
fieldNum{2}   = {[30,31],...
                 [30,31],...
                 [30,31],...
                 [36,37],...
                 [35,36],...
                 [31,32],...
                 [31 32],...
                 [33,34],...
                 [37,38],...
                 [37,38],...
                 [39,40],...
                 [33,34],...
                 [32,33],...
                 [32,33],...
                 [34,35],...
                 [31,32],...
                 [32,33],...
                 [30,31],...
                 [32,33],...
                 [32,33],...
                 [34,35],...
                 [30,31],...
                 [31,32],...
                 [30,31],...
                 [30,31],...
                 [30,31],...
                 [32,33],...
                 [30,31],...
                 [32,33],...
                 [32,33],...
                 [35,36]};
fieldNum{3}   = {[30,31],...
                 [32,33],...
                 [30,31],...
                 [31,32],...
                 [34 35],...
                 [32,33],...
                 [36,37],...
                 [31,32],...
                 [36,37],...
                 [34,35],...
                 [30,31],...
                 [30,31],...
                 [30,31],...
                 [32,33],...
                 [30,31],...
                 [31,32],...
                 [30,31],...
                 [30,31],...
                 [30,31],...
                 [32,33],...
                 [30,31],...
                 [30,31],...
                 [33,34],...
                 [32,33],...
                 [32,33],...
                 [30,31],...
                 [32,33],...
                 [30,31],...
                 [34,35],...
                 [31,32],...
                 [30,31]};
fieldNum{4}   = {[30,31],...
                 [32,33],...
                 [30,31],...
                 [32,33],...
                 [33,34],...
                 [31,32],...
                 [32 33],...
                 [30,31],...
                 [36,37],...
                 [32,33],...
                 [30,31],...
                 [30,31],...
                 [34,35],...
                 [30,31],...
                 [30,31],...
                 [31,32],...
                 [32,33],...
                 [30,31],...
                 [32,33],...
                 [30,31],...
                 [30,31],...
                 [30,31],...
                 [33,34],...
                 [32,33],...
                 [30,31],...
                 [32,33],...
                 [33,34],...
                 [33,34],...
                 [40,41],...
                 [30,31],...
                 [42,43]};

anatNum    = {[11:15],...
              [10:14],...
              [10:14],...
              [11:15],...
              [10:14],...
              [11:15],...
              [4:8],...
              [10:14],...
              [11:15],...
              [10:14],... % 2nd session
              [10:14],... % 2nd session
              [11:15],... % 2nd session
              [11:15],...
              [11:15],...
              [10:14],...
              [12:16],...
              [10:14],...
              [17:21],...
              [11:15],...
              [10:14],...
              [11:15],...
              [10:14],...
              [11:15],...
              [10:14],...
              [10:14],...
              [10:14],...
              [10:14],...
              [10:14],...
              [10:14],...
              [10:14],...
              [12:16]};
          
loc_AC     = {[-112 -165 -176],...
              [-106 -173 -163],...
              [-107 -178 -163],...
              [-103 -162 -167],...
              [-105 -163 -163],...
              [-107 -169 -156],...
              [-105 -163 -162],...
              [-103 -161 -175],...
              [-103 -159 -159],...
              [-107 -158 -191],...
              [-105 -163 -169],...
              [-106 -161 -178],...
              [-105 -160 -169],...
              [-108 -158 -181],...
              [-107 -167 -171],...
              [-104 -166 -175],...
              [-104 -161 -168],...
              [-105 -165 -176],...
              [-101 -163 -166],...
              [-102 -165 -178],...
              [-105 -154 -183],...
              [-107 -169 -174],...
              [-109 -163 -170],...
              [-105 -167 -168],...
              [-104 -169 -166],...
              [-108 -167 -174],...
              [-108 -175 -180],...
              [-108 -161 -171],...
              [-105 -170 -158],...
              [-105 -167 -182],...
              [-107 -169 -153]};

% Other random notes %

% Dicom system changed with s05 (sessN 3,4), s06 (sessN 1-4) and s07 (sessN 2-4)


% ------------------------------ Analysis Cases --------------------------------
   switch(what)
    
    case '0_MISC' % ------------ MISC: some aux. things ------------------------
    case 'MISC_check_time'                                                  % Check alignment of scanner and recorded time (sanity check): enter sn
        vararginoptions(varargin,{'sn'});
        %sn = Opt.sn;
        cd(behavDir);
        
        D=dload(sprintf('sml1_%s.dat',subj_name{sn}));
        figure('Name',sprintf('Timing of Task Onsets vs. TR onsets for Subj %d',sn),'NumberTitle','off')
        
        % plot alignment of TR time and trial onset time
        subplot(2,1,1); plot(D.startTimeReal/1000,(D.startTR-1)*TRlength/1000+D.startTRTime/1000)
        title('Alignment of TR time and Trial onset time')
        
        % plot difference of TR time and trial onset time
        subplot(2,1,2); plot((D.startTimeReal-(D.startTR-1)*TRlength+D.startTRTime)/1000)
        title('Difference of Trial onset time and TR time')
        xlabel('trial')
        ylabel('ms')
    case 'MISC_check_movement'                                              % Check movement of subject. Requires GLM 3 for specified subject.
        vararginoptions(varargin,{'sn','sessN'});
        
        glmSubjDir = [glmSessDir{sessN} filesep subj_name{sn}];
        cd(glmSubjDir);
        load SPM;
        spm_rwls_resstats(SPM)    
    case 'QC_motion'                                                        % Quality Control: motion
        vararginoptions(varargin,{'sn','sessN'});
        Q=[];
        for ss=sessN
            for s=sn
                subjDir = [glmSessDir{ss} filesep subj_name{s}];
                cd(subjDir);
                load SPM;
                M=motion_QC(SPM,'subjDir',fullfile(imagingDir,subj_name{s}));
                M.sn=s;
                M.sess=ss;
                Q=addstruct(Q,M);
            end
        end

        dircheck(QXDir);
        % save  
        save(fullfile(QCDir,'QC_motion'),'-struct','Q');
    case 'QC_betas'                                                         % Quality Control: contrast to noise ratio, pattern consistency
        betaChoice = 'raw';
        rmean=1;
        roi=[1:16];
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','rmean'});
        C=[];CC=[];
        
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_Brodmann_sess%d.mat',ss))); % loads region data (T)
            for s=sn
                for r=roi          
                    switch (betaChoice)
                        case 'uniPW'
                            beta = T.betaUW{T.SN==s & T.region==r};
                        case 'multiPW'
                            beta = T.betaW{T.SN==s & T.region==r};
                        case 'raw'
                            beta = T.betaRAW{T.SN==s & T.region==r};
                    end
                    
                    res=T.resMS{T.SN==s & T.region==r};
                    numCond=12;
                    numRun=8;
                    % make vectors for pattern consistency func
                    partVec=kron([1:numRun]',ones(numCond,1));
                    condVec=kron(ones(numRun,1),[1:numCond]');
                    
                    % calculate contrast to noise ratio
                    C.cnr = cnr_QC(beta,res,numCond,numRun);
                    % calculate pattern consistency
                    C.consist   = rsa_patternConsistency(beta,partVec,condVec,'removeMean',rmean);
                    [C.consist_crossval C.consist_voxel] = rsa_patternConsistency_crossval(beta,partVec,condVec,'removeMean',rmean,'voxSpecific',1);
                    C.sn=s;
                    C.roi=r;
                    C.sess=ss;     
                    CC=addstruct(CC,C);
                end
            end
        end
        dircheck(QCDir);
        % save  
        save(fullfile(QCDir,sprintf('QC_betas_%s',betaChoice)),'-struct','CC');
    case 'QC_plot'
       betaChoice = 'raw';
       roi=3;
       sessN=[1:4];
       vararginoptions(varargin,{'sessN','roi','betaChoice'});
       
       M=load(fullfile(QCDir,'QC_motion'));
       B=load(fullfile(QCDir,sprintf('QC_betas_%s.mat',betaChoice)));
       
       sub=1;
       for ss=sessN
           M2=getrow(M,M.sess==ss);
           B2=getrow(B,B.sess==ss & B.roi==roi);
           
           
           figure(1)
           subplot(1,numel(sessN),sub);
           plt.scatter(M2.fwd,B2.cnr);
           xlabel('framewise displacement');
           ylabel('pattern consistency');
           title(sprintf('sess%d',ss));
           
           figure(2)
           subplot(1,numel(sessN),sub);
           plt.scatter(M2.rms,B2.consist);
           xlabel('root mean square displacement');
           ylabel('pattern consistency');
           title(sprintf('sess%d',ss));
           
           sub=sub+1;
       end
       
    case 'PHYSIO_read'                                              % Process pulse and respiration!! DEPRECIATED   
        % Read raw physio and make regressors for both sessions
        vararginoptions(varargin,{'sn'});

        % indicator for number of loc/func runs (across sessions)
        indx_lc = 1;
        indx_fn = 1;
        
        for sessN = 1:sess_sn(sn)   % do for all sessions for that subject
            physioRawDir=[dicomDir,'/',sprintf('%s_%d',subj_name{sn},sessN),'/Physio'];
            
            snDir = fullfile(dicomDir,sprintf('%s_%d',subj_name{sn},sessN));
            ndicom = DicomName{sessN}{sn};
            nscans = fscanNum{sessN}{sn};
            
            D = physio_processSiemens(snDir,ndicom,nscans,'dummyTR',numDummys,'physioDir',physioRawDir);
            
            [~,W] = physio_getCardiacPhase(D,'fig',1);
            
            % split into localizer and functional runs - separate glms used!
            W_LOC = {W{4} W{8}};  % localizer
            W_FUNC = {W{1} W{2} W{3} W{5} W{6} W{7} W{9} W{10}};   % functional
            
            % ensure all the runs have equal length - cut the longer runs
            for i = 1: numruns_task_sess
                if length(W_FUNC{i}) ~= numTRs(1)-numDummys
                    temp=W_FUNC{i}(1:(numTRs(1)-numDummys)); % store the initial values into a temp var
                    W_FUNC{i}=zeros(1,(numTRs(1)-numDummys)); % clear the long run
                    W_FUNC{i}=temp;
                end
            end
            
            [CregsLOC,CnamesLOC] = physio_makeCardiacRegressors(W_LOC);     % localizer
            [CregsFUNC,CnamesFUNC] = physio_makeCardiacRegressors(W_FUNC);   % functional runs
            
            % re-arrange in a structure with cos/sin regressors per run
            
            % localizer
            
            for r = 1:numruns_loc_sess
                L=[];
                L.regSin=CregsLOC((r*2)-1,:); % use sth else, not r - increments after every loop
                L.regCos=CregsLOC((r*2),:);
                P_loc{indx_lc}=L;
                indx_lc=indx_lc+1;
            end
            
            for r = 1:numruns_task_sess
                T=[];
                T.regSin=CregsFUNC((r*2)-1,:);
                T.regCos=CregsFUNC((r*2),:);
                P_func{indx_fn}=T;
                indx_fn=indx_fn+1;
            end
        end
        
        % save regressors for localizer and functional runs
        dircheck(fullfile(physioDir,subj_name{sn}));
        save(fullfile(physioDir,subj_name{sn},'physioLoc.mat'),'P_loc');
        save(fullfile(physioDir,subj_name{sn},'physioFunc.mat'),'P_func');

    case 'BEH_scanner_er_mt' 
        % Analyse and plot ER and MT per subject
        sessN=4;
        vararginoptions(varargin,{'sn','sessN'});
        
        B_loc = [];
        B_fun = [];
        
        for ss = 1:sessN

            for s = sn

                D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));
                
                R = getrow(D,D.ScanSess==ss);
                F = getrow(R,R.blockType==3 | R.blockType==9);  % functional runs  
                L = getrow(R,R.blockType==4);   % localiser
                % calculate / save variables
                % Loc - MT, ER
                BL.MT = L.MT;
                BL.ER = L.isError;
                BL.ScanSess = L.ScanSess;
                BL.sn = s*ones(size(L.ScanSess));
                % Func - MT, ER (for trained and untrained separately)
                BF.MT = F.MT;
                BF.ER = F.isError;
                BF.ScanSess = F.ScanSess;
                BF.seqType = F.seqType;
                BF.sn = s*ones(size(F.ScanSess));
                % Add all var into structures B_loc and B_fun
                
                B_loc = addstruct(B_loc,BL);
                B_fun = addstruct(B_fun,BF);
            end
        end
        
        lab = {'Sess1','Sess2'};
        
        figure(1)
        subplot(2,2,1)
        barplot(B_loc.ScanSess,B_loc.MT); ylabel('Mov Time'); title('Localiser'); xlabel('Session');
        subplot(2,2,2)
        barplot(B_fun.ScanSess,B_fun.MT,'split',B_fun.seqType,'leg',{'train','untrain'},'leglocation','northeast'); ylabel('Mov Time'); xlabel('Session'); title('Functional runs');
        subplot(2,2,3)
        barplot(B_loc.ScanSess,B_loc.ER); ylabel('Error rate'); xlabel('Session');
        subplot(2,2,4)
        barplot(B_fun.ScanSess,B_fun.ER,'split',B_fun.seqType,'leg',{'train','untrain'},'leglocation','northeast'); ylabel('Error rate'); xlabel('Session');
        
        figure;
        for i=1:5
            subplot(1,5,i)
            barplot(B_fun.ScanSess,B_fun.MT,'split',B_fun.seqType,'subset',B_fun.sn==i,'leg',{'train','untrain'})
        end
        
        figure;
        sty = style.custom({'red','blue'},'markersize',12);
        plt.box(B_fun.ScanSess,B_fun.MT,'split',B_fun.seqType,'style',sty,'plotall',0,'leg',{'trained','untrained'});
        
        keyboard;
    case 'BEH_training'
        vararginoptions(varargin,{'sn'});
        
        AllSubj = [];
        figure        

        for s = sn
            
            D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));
            
            R = getrow(D,D.blockType==6);
            
            runs = 1:length(unique(R.BN));
            blockNum = kron(runs',ones(24,1)); % 24 trials
            subplot(1,numel(sn),s)
            lineplot(blockNum,R.MT,'style_thickline');
            
            R.sn = s*ones(size(R.timeTrial));
            AllSubj = addstruct(AllSubj,R);
            
        end
        
        keyboard;
 
    case '1_PREP' % ------------ PREP: preprocessing. Expand for more info. ----
        % The PREP cases are preprocessing cases.
        % You should run these in the following order:
        %       'PREP_dicom_import'*    :  call with 'series_type','functional',
        %                                  'series_type','anatomical'
        %                                   and 'series_type','fieldmap'.
        %       'PREP_process1_func'    :   Runs steps 1.3 - 1.8 (see below).
        %       'PREP_precess1_anat'    :   Runs steps 1.9 - 1.11 (see below).
        %       'PREP_coreg'*           :   Registers meanepi to anatomical img. (step 1.12)
        %       'PREP_process2'*        :   Runs steps 1.9 - 1.11 (see below).
        %
        %   * requires user input/checks after running BEFORE next steps.
        %       See corresponding cases for more info about required
        %       user input.
        %
        % When calling any case, you can submit an array of Subj#s as so:
        %       ('some_case','sn',[Subj#s])
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'PREP_dicom_import'                                                % STEP 1.1/2/3   :  Import functional/anatomical/field map dicom series: enter sn, series type
        % converts dicom to nifti files w/ spm_dicom_convert
        
        series_type = 'functional';
        sessN = 1;
        vararginoptions(varargin,{'sn','sessN','series_type'});
        cwd = pwd;
       
        
        for ss=sessN
            switch series_type
                case 'functional'
                    seriesNum = fscanNum{ss};
                case 'anatomical'
                    seriesNum = anatNum;
                case 'fieldmap'
                    seriesNum = fieldNum{ss};
            end
            % Loop through subjects
            for s = sn;
                dircheck(fullfile(dicomDir,[subj_name{s},sprintf('_%d',ss)]));
                cd(fullfile(dicomDir,[subj_name{s},sprintf('_%d',ss)]));
                
                % For each series number of this subject (in 'Subject Things')
                for i=1:length(seriesNum{s})
                    r     = seriesNum{s}(i);
                    % Get DICOM FILE NAMES
                    if (s<5 || s==5 && ss<3 || s==7 && ss<2)
                        DIR   = dir(sprintf('%s.%4.4d.*.IMA',DicomName{ss}{s},r));   % Get DICOM FILE NAMES
                    else
                        folder = fullfile(dicomDir,[subj_name{s},sprintf('_%d',ss)],[sprintf('%4.4d',r)]);
                        cd(folder)
                        DIR   = dir(sprintf('%s.%4.4d.*.dcm',DicomName{ss}{s},r));   % Get DICOM FILE NAMES
                    end
                    Names = vertcat(DIR.name);
                    % Convert the dicom files with these names.
                    if (~isempty(Names))
                        % Load dicom headers
                        HDR=spm_dicom_headers(Names,1);
                        
                        % Make a directory for series{r} for this subject.
                        % The nifti files will be saved here.
                        
                        dirname = fullfile(dicomDir,[subj_name{s},sprintf('_%d',ss)],sprintf('series%2.2d',r));
                        dircheck(dirname);
                        % Go to the dicom directory of this subject
                        cd(dirname);
                        % Convert the data to nifti
                        spm_dicom_convert(HDR,'all','flat','nii');
                        cd ..
                    end
                    display(sprintf('Series %d done \n',seriesNum{s}(i)))
                end
                % Display verbose messages to user.
                switch series_type
                    case 'functional'
                        fprintf('Subject %02d session %d functional runs imported. Copy the unique .nii name for subj files and place into ''Subject Things''.\n',s,ss)
                    case 'anatomical'
                        fprintf('Anatomical runs have been imported for subject %d session %d.\n',s,ss);
                        fprintf('Please locate the T1 weighted anatomical img. Copy it to the anatomical folder.\n')
                        fprintf('Rename this file to ''%s_anatomical_raw.nii'' in the anatomical folder.\n',subj_name{s});
                    case 'fieldmap'
                        fprintf('Subject %s session %d fieldmaps imported.\n',subj_name{s},ss);
                        fieldmapfold = fullfile(fieldmapDir,subj_name{s},sprintf('sess%d',ss));
                        dircheck(fieldmapfold);
                        fprintf('The subfolder ''sess%d'' in the subject fieldmap folder ''%s'' was created for you.\n',ss,subj_name{s});
                        fprintf('Please locate the magnitude and phase files (different series).\n');
                        fprintf('Rename the files into ''%s_magnitude.nii'' and ''%s_phase.nii''.\n',subj_name{s},subj_name{s});
                end
            end; % sn
        end; % sessN
        cd(cwd);

    case 'PREP_process1_func'                                               
        % need to have dicom_import done prior to this step.
        vararginoptions(varargin,{'sn','sessN'});
        
        for s = sn
            sml1_imana_prep('PREP_make_4dNifti','sn',s,'sessN',sessN);
            sml1_imana_prep('PREP_makefieldmap','sn',s,'sessN',sessN);
            sml1_imana_prep('PREP_make_realign_unwarp','sn',s,'sessN',sessN);
            sml1_imana_prep('PREP_move_data','sn',s,'sessN',sessN);
            sml1_imana_prep('PREP_meanimage_bias_correction','sn',s);
        end
    case 'PREP_make_4dNifti'                                                % STEP 1.4       :  Converts dicoms to 4D niftis out of your raw data files
        vararginoptions(varargin,{'sn','sessN'});
        for s = sn
            for ss = 1:sessN
                % For each functional run
                for i = 1:length(fscanNum{ss}{s})
                    outfilename = fullfile(imagingDirRaw,subj_name{s},sprintf('sess%d',ss),sprintf('%s_run_%2.2d.nii',subj_name{s},run_num{ss}(i)));
                    % Create a 4d nifti of all functional imgs in this run.
                    % Don't include the first few dummy scans in this 4d nifti.
                    P={};
                    for j = 1:(numTRs((ss-1)*10+i)-numDummys)
                        P{j}=fullfile(dicomDir,[subj_name{s},sprintf('_%d',ss)],sprintf('series%2.2d',fscanNum{ss}{s}(i)),...
                            sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{ss}{s},fscanNum{ss}{s}(i),j+numDummys,j+numDummys));
                    end;
                    dircheck(fullfile(imagingDirRaw,subj_name{s},sprintf('sess%d',ss)))
                    spm_file_merge(char(P),outfilename);
                    fprintf('Run %d in session %d done -> overall run %d\n',i,ss,run_num{ss}(i));
                end; % run
            end; % session - ss
        end; % sn
    case 'PREP_makefieldmap'                                                % STEP 1.5       :  Create field map
        prefix = '';
        vararginoptions(varargin,{'sn','sessN'});
        for ss=1:sessN
            subfolderRawdata    = sprintf('sess%d',ss);
            subfolderFieldmap   = sprintf('sess%d',ss);
            
            spmj_makefieldmap(baseDir, subj_name{sn}, runs{ss},'prefix',prefix,'subfolderRawdata',subfolderRawdata,'subfolderFieldmap',subfolderFieldmap);
        end
    case 'PREP_make_realign_unwarp'                                         % STEP 1.6       :  Realign + unwarp functional runs
        prefix  ='';
        vararginoptions(varargin,{'sn','sessN'});
                
        subfolderRawdata    = {'sess1','sess2','sess3','sess4'};
        subfolderFieldmap   = {'sess1','sess2','sess3','sess4'};
        subj_runs           = runs;
        
       % subfolderRawdata    = subfolderRawdata(sessN);
        subfolderRawdata    = subfolderRawdata(1:sessN);
       % subfolderFieldmap   = subfolderFieldmap(sessN);
        subfolderFieldmap   = subfolderFieldmap(1:sessN);
       % subj_runs           = subj_runs(sessN);
        subj_runs           = subj_runs(1:sessN);

        spmj_realign_unwarp_sess(baseDir, subj_name{sn}, subj_runs, numTRs, 'prefix',prefix, 'subfolderRawdata',subfolderRawdata,'subfolderFieldmap',subfolderFieldmap);      
    case 'PREP_plot_movementparameters'                                     % OPTIONAL       :  Investigate movement parameters
        vararginoptions(varargin,{'sn','sessN'});
        X=[];
        %r=[1:numruns(sn)];
        for r=1:10 % 10 functional runs
            x = dlmread (fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rp_' subj_name{sn},'_run' runs{sessN}{r},'.txt']));
            X = [X; x];
        end
       
        figure;
        clr = hsv(3);
        subplot(2,1,1); 
        for i = 1:3
            plot(X(:,i),'Color',clr(i,:))
            hold on;
        end
        legend('x', 'y', 'z', 'location' , 'EastOutside')

        subplot(2,1,2); 
        for j = 1:3
            plot(X(:,j+3)*180/pi,'Color',clr(j,:));
            hold on;
        end
        legend('pitch', 'roll', 'yaw', 'location' , 'EastOutside') 
    case 'PREP_move_data'                                                   % STEP 1.7       :  Moves subject data from raw directories to working directories
        % Moves image data from imaging_dicom_raw into a "working dir":
        % imaging_dicom.
        vararginoptions(varargin,{'sn','sessN'});

        prefix='';
        dircheck(fullfile(baseDir, 'imaging_data',subj_name{sn}));
        for ss=1:sessN;
            disp(['Sess' num2str(ss)]);
            for r=1:numruns_sess;
                
                source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},sprintf('sess%d',ss), ['u' prefix subj_name{sn},'_run',runs{ss}{r},'.nii']);
                dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['u' prefix subj_name{sn},'_run',runs{ss}{r},'.nii']);
                
                copyfile(source,dest);
                source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},sprintf('sess%d',ss), ['rp_' subj_name{sn},'_run',runs{ss}{r},'.txt']);
                dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rp_' subj_name{sn},'_run',runs{ss}{r},'.txt']);
                
                copyfile(source,dest);
            end;
            
            if ss == 1
                source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},sprintf('sess%d',ss), ['meanu' prefix subj_name{sn},'_run',runs{1}{1},'.nii']); %first run of first session used for mean
                dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);
                
                copyfile(source,dest);
            end
        end   

    %__________________________________________________________________
    case 'PREP_meanimage_bias_correction'                                   % STEP 1.8       :  Bias correct mean image prior to coregistration
        vararginoptions(varargin,{'sn'});
        
        % make copy of original mean epi, and work on that
        source  = fullfile(baseDir, 'imaging_data',subj_name{sn},['meanepi_' subj_name{sn} '.nii']);
        dest    = fullfile(baseDir, 'imaging_data',subj_name{sn},['bmeanepi_' subj_name{sn} '.nii']);
        copyfile(source,dest);
        
        % bias correct mean image for grey/white signal intensities 
        P{1}    = dest;
        spmj_bias_correct(P);
  
    case 'PREP_process1_anat'                                               
        % need to have dicom_import done prior to this step.
        % only run once per subject (after first session)
        vararginoptions(varargin,{'sn'});
        
        for s = sn
            sml1_imana_prep('PREP_reslice_LPI','sn',s);
            sml1_imana_prep('PREP_centre_AC','sn',s);
            sml1_imana_prep('PREP_segmentation','sn',s);
        end         
    case 'PREP_reslice_LPI'                                                 % STEP 1.9       :  Reslice anatomical image within LPI coordinate systems
        vararginoptions(varargin,{'sn'});
        
        % (1) Reslice anatomical image to set it within LPI co-ordinate frames
        source  = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical_raw','.nii']);
        dest    = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']);
        spmj_reslice_LPI(source,'name', dest);
        
        % (2) In the resliced image, set translation to zero
        V               = spm_vol(dest);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = [0 0 0];
        spm_write_vol(V,dat);
        display 'Manually retrieve the location of the anterior commissure (x,y,z) before continuing'
        
        
        %___________
    case 'PREP_centre_AC'                                                   % STEP 1.10      :  Re-centre AC in anatomical image
        % Set origin of anatomical to anterior commissure (must provide
        % coordinates in section (4)).
        vararginoptions(varargin,{'sn'});
        
        img    = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']);
        V               = spm_vol(img);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = loc_AC{sn};
        spm_write_vol(V,dat);
        display 'Done'
        
        
        %_____
    case 'PREP_segmentation'                                                % STEP 1.11      :  Segmentation & normalization
        vararginoptions(varargin,{'sn'});

        SPMhome=fileparts(which('spm.m'));
        J=[];
        for s=sn
            if s==27
                home=fullfile(anatomicalDir,subj_name{sn});
            else
                home=SPMhome;
            end
            J.channel.vols = {fullfile(anatomicalDir,subj_name{sn},[subj_name{sn},'_anatomical.nii,1'])};
            J.channel.biasreg = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write = [0 0];
            J.tissue(1).tpm = {fullfile(home,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm = {fullfile(home,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm = {fullfile(home,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm = {fullfile(home,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm = {fullfile(home,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm = {fullfile(home,'tpm/TPM.nii,6')};
            J.tissue(6).ngaus = 2;
            J.tissue(6).native = [0 0];
            J.tissue(6).warped = [0 0];
            J.warp.mrf = 1;
            J.warp.cleanup = 1;
            J.warp.reg = [0 0.001 0.5 0.05 0.2];
            J.warp.affreg = 'mni';
            J.warp.fwhm = 0;
            J.warp.samp = 3;
            J.warp.write = [1 1];   
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);
            fprintf('Check segmentation results for %s\n', subj_name{s})
        end;
    %__________________________________________________________________
    
    case 'PREP_coreg'                                                       % STEP 1.12      :  Coregister meanepi to anatomical image - only needs to be done for 1st session
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi image
        % Note: this only needs to be done for the 1st session, because..
        % .. mean epi is taken from the 1st run of the 1st session
        vararginoptions(varargin,{'sn'});
        
        cd(fullfile(anatomicalDir,subj_name{sn}));
        coregtool;
        keyboard();
        
        % (2) Automatically co-register functional and anatomical images
        %sn=varargin{1};
        
        J.ref = {fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii'])};
        J.source = {fullfile(imagingDir,subj_name{sn},['rbbmeanepi_' subj_name{sn} '.nii'])}; 
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % (3) Manually check again
        coregtool;
        keyboard();
        
        % NOTE:
        % Overwrites meanepi, unless you update in stsml1ep one, which saves it
        % as rmeanepi.
        % Each time you click "update" in coregtool, it saves current
        % alignment by appending the prefix 'r' to the current file
        % So if you continually update rmeanepi, you'll end up with a file
        % called r...rrrmeanepi.
      
        %__________________________________________________________________

    case 'PREP_process2'                                                                                                     
        vararginoptions(varargin,{'sn'});
        
        for s=sn
            sml1_imana_prep('PREP_make_samealign','sn',s);
            sml1_imana_prep('PREP_check_samealign','sn',s);
            sml1_imana_prep('PREP_make_maskImage','sn',s);
        end
    case 'PREP_make_samealign'                                              % STEP 1.13     :  Align to first image (rbmeanepi_* of first session)
        prefix  = 'u';
        vararginoptions(varargin,{'sn'});

        cd(fullfile(imagingDir,subj_name{sn}));

        % Select image for reference
        P{1} = fullfile(imagingDir,subj_name{sn},sprintf('rbbmeanepi_%s.nii',subj_name{sn}));

        % Select images to be realigned
        Q={};
        for r=1:numruns(sn)
            for i=1:numTRs(r)-numDummys;
                Q{end+1}    = fullfile(imagingDir,subj_name{sn},...
                    sprintf('%s%s_run_%2.2d.nii,%d',prefix, subj_name{sn},r,i));
            end;
        end;
        spmj_makesamealign_nifti(char(P),char(Q));
        % Run spmj_makesamealign_nifti to bring all functional runs into
        % same space as realigned mean epis  
    case 'PREP_check_samealign'                                             % OPTIONAL      :  Check if all functional scans are align to the anatomical
        prefix='u';
        vararginoptions(varargin,{'sn'});
        Q={};
        cd(fullfile(imagingDir,subj_name{sn}));
        for i=1:numruns(sn)
            r{i}=int2str(i);
        end;
        for r= 1:numel(r)
            for i=1:numTRs(r)-numDummys;
                %Q{end+1} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix, subj_name{sn},'_run',r(r),'.nii,',num2str(i)])];
                Q{end+1}    = fullfile(imagingDir,subj_name{sn},...
                    sprintf('%s%s_run_%2.2d.nii,%d',prefix, subj_name{sn},r,i));
            end
        end
        P{1}= fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rbbmeanepi_' subj_name{sn} '.nii']);
        spmj_checksamealign(char(P),char(Q))       
    case 'PREP_make_maskImage'                                              % STEP 1.14     :  Make mask images (noskull and gray_only) - only for the 1st session (using mean epi)
    vararginoptions(varargin,{'sn'});

    for sn=sn
        cd(fullfile(imagingDir,subj_name{sn}));

        nam{1}  = fullfile(imagingDir,subj_name{sn}, ['rbbmeanepi_' subj_name{sn} '.nii']);
        nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
        nam{3}  = fullfile(anatomicalDir, subj_name{sn}, ['c2' subj_name{sn}, '_anatomical.nii']);
        nam{4}  = fullfile(anatomicalDir, subj_name{sn}, ['c3' subj_name{sn}, '_anatomical.nii']);
        spm_imcalc_ui(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')

        nam={};
        nam{1}  = fullfile(imagingDir,subj_name{sn}, ['rbbmeanepi_' subj_name{sn} '.nii']);
        nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
        spm_imcalc_ui(nam, 'rmask_gray.nii', 'i1>1 & i2>0.4')

    end
    
    case '2_SURF' % ------------ SURF: Freesurfer funcs. Expand for more info. -
        % The SURF cases are the surface reconstruction functions. Surface
        % reconstruction is achieved via freesurfer.
        % All functions can be called with ('SURF_processAll','sn',[Subj#s]).
        % You can view reconstructed surfaces with Caret software.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'SURF_processAll'                                                 
        vararginoptions(varargin,{'sn'});
        % You can call this case to do all the freesurfer processing.
        % 'sn' can be an array of subjects because each processing case
        % contained within loops through the subject array submitted to the
        % case.
        sml1_imana_prep('SURF_freesurfer','sn',sn);
        sml1_imana_prep('SURF_xhemireg','sn',sn);
        sml1_imana_prep('SURF_map_ico','sn',sn);
        sml1_imana_prep('SURF_make_caret','sn',sn);
    case 'SURF_freesurfer'                                                  % STEP 2.1   :  Call recon-all in                                                                                                                  Freesurfer                                                   
        vararginoptions(varargin,{'sn'});
        for i=sn
            freesurfer_reconall(freesurferDir,subj_name{i},fullfile(anatomicalDir,subj_name{i},[subj_name{i} '_anatomical.nii']));
        end
    case 'SURF_xhemireg'                                                    % STEP 2.2   :  Cross-Register surfaces left / right hem
        vararginoptions(varargin,{'sn'});
        for i=sn
            freesurfer_registerXhem({subj_name{i}},freesurferDir,'hemisphere',[1 2]); % For debug... [1 2] orig
        end;
    case 'SURF_map_ico'                                                     % STEP 2.3   :  Align to the new atlas surface (map icosahedron)
        vararginoptions(varargin,{'sn'});
        for i=sn
            freesurfer_mapicosahedron_xhem(subj_name{i},freesurferDir,'smoothing',1,'hemisphere',[1:2]);
        end;
    case 'SURF_make_caret'                                                  % STEP 2.4   :  Translate into caret format
        vararginoptions(varargin,{'sn'});
        for i=sn
            caret_importfreesurfer(['x' subj_name{i}],freesurferDir,caretDir);
        end;
     
    case '3a_GLM_SESS' % ------------- GLM per session! ------------------
        % makes the glm per subject and per session
    case 'GLM_sess_all'
        glm=2;
        vararginoptions(varargin,{'sn','glm','sessN'});
        sml1_imana_prep('GLM_make_sess','sn',sn,'glm',glm,'sessN',sessN);
        sml1_imana_prep('GLM_estimate_sess','sn',sn,'sessN',sessN,'glm',glm);
        sml1_imana_prep('GLM_contrast_sess','sn',sn,'sessN',sessN,'glm',glm);
    case 'GLM_make_sess'
        % makes the GLM file for each subject, and a corresponding 
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        
        glm = 3;    %1/2/3
        vararginoptions(varargin,{'sn','glm','sessN'});
        % Set some constants.
        prefix		 = 'u';
        T			 = [];
        dur			 = 2.5;                                                 % secs (length of task dur, not trial dur)
        % adjusting hrf per subject & session based on extracted timeseries!  
        delay     = [0.5 1 1 0.5 1 0 0 0.5 1 1 0.5 1 1 1 1 0 1 1.5 1 1 1 1 1 1 0.5 0 1 1 1 1 1];  
        %delay     = [0.5 1 1 0.5 1 0 0 0.5 1 1 0.5 1 1 1 1 0 0.5 1.5];  
        
        announceTime = 0;                                                 % length of task announce time - currently not used
        % Gather appropriate GLM presets.
        switch glm
            case 1  % wls
                hrf_params = [5.5 12.5]; 
                hrf_cutoff = 128;
                cvi_type   = 'wls';
            case 2  % fast + hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = 128;
                cvi_type   = 'fast';
            case 3  % fast, no hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = inf;
                cvi_type   = 'fast';
        end
        for ss=sessN % session
            % Loop through subjects / all sessions and make SPM files.
            for s = sn
                D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));
                T=[];
                
                % Do some subject structure fields.
              %  dircheck(fullfile(glmSessDir{sessN}, subj_name{s}));
                dircheck(fullfile(glmSessDir{ss},subj_name{s})); % add glm type
              %  J.dir 			 = {fullfile(glmSessDir{ss}, subj_name{s})};
                J.dir 			 = {fullfile(glmSessDir{ss},subj_name{s})};
                J.timing.units   = 'secs';                                      % timing unit that all timing in model will be
                J.timing.RT 	 = 1.0;                                         % TR (in seconds, as per 'J.timing.units')
                J.timing.fmri_t  = 16;
                J.timing.fmri_t0 = 1;
                
                L = getrow(D,D.ScanSess==ss);    % only blocks of that scan session
                
                % Loop through runs.
                for r = 1:numruns_task_sess
                    if ss == 4
                        Rr = getrow(L,L.blockType==9); %blockType==9 - func imaging run without metronome
                    else
                        Rr = getrow(L,L.blockType==3); %blockType==3 - funct imaging run with metornome
                    end
                    
                    uniqrun=unique(Rr.BN);
                    R = getrow(Rr,Rr.BN==uniqrun(r)); % 1-8 func runs of the session
                    
                    for i = 1:(numTRs(run_task(r))-numDummys)                   % get nifti filenames, correcting for dummy scancs
                        
                        N{i} = [fullfile(baseDir, 'imaging_data',subj_name{s}, ...
                            [prefix subj_name{s},'_run',runs{ss}{run_task(1,r)},'.nii,',num2str(i)])];
                        
                    end;
                    J.sess(r).scans = N;                                        % number of scans in run
                    % Loop through conditions.
                    
                    for c = 1:numel(num_seq)
                        idx						   = find(R.seqNumb==c);             % find indx of all trials in run - 1:6 trained; 7-12 untrained
                        condName = sprintf('SeqNumb-%d',R.seqNumb(idx(1)));
                        J.sess(r).cond(c).name 	   = condName;
                        % Correct start time for numDummys removed & convert to seconds
                        J.sess(r).cond(c).onset    = [R.startTimeReal(idx)/1000 - J.timing.RT*numDummys + announceTime + delay(s)];
                        J.sess(r).cond(c).duration = dur;                       % durations of task we are modeling (not length of entire trial)
                        
                        J.sess(r).cond(c).tmod     = 0;
                        J.sess(r).cond(c).orth     = 0;
                        J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                        
                        % Do some subject info for fields in SPM_info.mat.
                        S.SN    		= s;
                        S.run   		= r;    % 1-8: functional runs
                        S.runAll        = (ss-1)*8 + r;  % 1-32
                        S.seqNumb 		= R.seqNumb(idx(1));
                        S.seqType    	= R.seqType(idx(1));
                        S.isMetronome   = R.isMetronome(idx(1));
                        S.ScanSess      = R.ScanSess(idx(1));
                        T				= addstruct(T,S);
                    end;
                    
                    % Add any additional regressors here.
                    J.sess(r).multi 	= {''};
                    J.sess(r).regress 	= struct('name', {}, 'val', {});
                    J.sess(r).multi_reg = {''};
                    % Define high pass filter cutoff (in seconds): see glm cases.
                    J.sess(r).hpf 		= hrf_cutoff;
                end;
                J.fact 			   = struct('name', {}, 'levels', {});
                J.bases.hrf.derivs = [0 0];
                J.bases.hrf.params = hrf_params;    % make it subject specific
                J.volt 			   = 1;
                J.global 		   = 'None';
                J.mask 	           = {fullfile(baseDir, 'imaging_data',subj_name{s}, 'rmask_noskull.nii,1')};
                J.mthresh 		   = 0.05;
                J.cvi_mask 		   = {fullfile(baseDir, 'imaging_data',subj_name{s},'rmask_gray.nii')};
                J.cvi 			   = cvi_type;
                % Save the GLM file for this subject.
                spm_rwls_run_fmri_spec(J);
                % Save the aux. information file (SPM_info.mat).
                % This file contains user-friendly information about the glm
                % model, regressor types, condition names, etc.
                save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            end;    % sn
        end
    case 'GLM_estimate_sess'
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        vararginoptions(varargin,{'sn','sessN','glm'});
        for ss = sessN
            for s = sn
                % Load files
              %  load(fullfile(glmSessDir{ss},subj_name{s},'SPM.mat'));
                load(fullfile(glmSessDir{ss},subj_name{s},'SPM.mat')); % add glm type
                SPM.swd = fullfile(glmSessDir{ss},subj_name{s});
                % Run the GLM.
                spm_rwls_spm(SPM);
            end; % sn
        end; %sessN
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)
    case 'GLM_contrast_sess'
        % enter sn, sessN #
        % 1:   Trained seq vs. rest
        % 2:   Novel seq vs. rest
        
        vararginoptions(varargin,{'sn','sessN','glm'});
        cwd = pwd;
        for ss=sessN
            % Loop through subjects.
            for s = sn
               % glmSubjDir = [glmSessDir{ss} filesep subj_name{s}];
                glmSubjDir = fullfile(glmSessDir{ss},subj_name{s}); % add glm type
                cd(glmSubjDir);
                
                load SPM;
                SPM = rmfield(SPM,'xCon');
                T   = load('SPM_info.mat');
                
                nrun    = numel(SPM.nscan);
                tt      = repmat([1:numel(num_seq)],1,nrun);
                
                % (1) t contrasts each sequence against rest
                % (1:12)
                for c = 1:numel(num_seq)
                    con = zeros(1,size(SPM.xX.X,2));
                    con(tt==c) = 1;
                    con = con/sum(con);
                    SPM.xCon(c) = spm_FcUtil('Set',sprintf('Seq%d',c), 'T', 'c',con',SPM.xX.xKXs);
                end
                
                % (2) t contrast for trained seq vs. rest
                con                = zeros(1,size(SPM.xX.X,2));
                con(:,T.seqNumb<7) = 1;
                con                = con/sum(con);
                SPM.xCon(end+1)    = spm_FcUtil('Set',sprintf('TrainSeq'), 'T', 'c',con',SPM.xX.xKXs);
                
                % (3) t contrast for novel seq vs. rest
                con                = zeros(1,size(SPM.xX.X,2));
                con(:,T.seqNumb>6 & T.seqNumb<13)  = 1;
                con                = con/sum(con);
                SPM.xCon(end+1)    = spm_FcUtil('Set',sprintf('UntrainSeq'), 'T', 'c',con',SPM.xX.xKXs);
                
                %____do the constrasts
                SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
                save('SPM.mat','SPM');
                
                % rename contrast images and spmT images
                conName = {'con','spmT'};
                for i=1:length(SPM.xCon),
                    for n=1:numel(conName),
                        oldName{i} = fullfile(glmSubjDir,sprintf('%s_%2.4d.nii',conName{n},i));
                        newName{i} = fullfile(glmSubjDir,sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                        movefile(oldName{i},newName{i});
                    end
                end
            end; % sn
        end; %sessN
        cd(cwd);
        
    case '3b_GLM_RepSup' % ------- GLM with separate regressors for first / second execution --
        % fits regressors separately for FoSEx
    case 'GLM_RepSup_all'
        glm=3;
        vararginoptions(varargin,{'sn','glm','sessN'});
        sml1_imana_prep('GLM_make_RepSup','sn',sn,'glm',glm,'sessN',sessN);
        sml1_imana_prep('GLM_estimate_RepSup','sn',sn,'sessN',sessN);
        sml1_imana_prep('GLM_contrast_RepSup','sn',sn,'sessN',sessN);
    case 'GLM_make_RepSup'
        % functional runs - separate regressors for first / second execution
        % makes the GLM file for each subject, and a corresponding
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        
        glm = 3;
        vararginoptions(varargin,{'sn','glm','sessN'});
        % Set some constants.
        prefix		 = 'u';
        dur			 = 2.5;                                                 % secs (length of task dur, not trial dur)
        %delay     = [0.5 1 1 0.5 1 0 0 0.5 1 1 0.5 1 1 1 1];              % adjusting hrf per subject based on extracted timeseries!
        delay     = [0.5 1 1 0.5 1 0 0 0.5 1 1 0.5 1 1 0.5 0.5 0 0.5 1 0.5 1 1 1 1 1 1 1 1 1 1 1 1]; 
        announceTime = 0;                                                   % length of task announce time - currently not used
        % Gather appropriate GLM presets.
        switch glm
            case 1  % wls
                hrf_params = [5.5 12.5];
                hrf_cutoff = 128;
                cvi_type   = 'wls';
            case 2  % fast + hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = 128;
                cvi_type   = 'fast';
            case 3  % fast, no hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = inf;
                cvi_type   = 'fast';
        end
        
        for ss = sessN
            % Loop through subjects and make SPM files.
            for s = sn
                D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));
                T = [];
                % Do some subject structure fields.
                dircheck(fullfile(glmFoSExDir{ss}, subj_name{s}));
                J.dir 			 = {fullfile(glmFoSExDir{ss}, subj_name{s})};
                J.timing.units   = 'secs';                                      % timing unit that all timing in model will be
                J.timing.RT 	 = 1.0;                                         % TR (in seconds, as per 'J.timing.units')
                J.timing.fmri_t  = 16;
                J.timing.fmri_t0 = 1;
                
                L = getrow(D,D.ScanSess==ss);    % only blocks of that scan session

                % Loop through sessions
                % Loop through runs.
                for r = 1:numruns_task_sess       % 8 functional runs
                    if sessN == 4
                        Rr = getrow(L,L.blockType==9); %blockType==9 - func imaging run without metronome
                    else
                        Rr = getrow(L,L.blockType==3); %blockType==3 - funct imaging run with metornome
                    end
                    
                    uniqrun=unique(Rr.BN);
                    R = getrow(Rr,Rr.BN==uniqrun(r)); % 1-8 func runs of the session

                    for i = 1:(numTRs(run_task(r))-numDummys)                   % get nifti filenames, correcting for dummy scancs
                        
                        N{i} = [fullfile(baseDir, 'imaging_data',subj_name{s}, ...
                            [prefix subj_name{s},'_run',runs{ss}{run_task(1,r)},'.nii,',num2str(i)])];
                        
                    end;
                    J.sess(r).scans = N;                                        % number of scans in run
                    % Loop through conditions.
                    for FoSEx = 1:2     % first / second execution
                        for c = 1:numel(num_seq) % each sequence
                            idx = find(R.seqNumb==c & R.FoSEx==FoSEx);   % find indx of all trials in run - 1:6 trained; 7-12 untrained
                            condName = sprintf('SeqNumb-%d-Ex%d',R.seqNumb(idx(1)),R.FoSEx(idx(1)));
                            J.sess(r).cond((FoSEx-1)*numel(num_seq)+c).name 	   = condName;
                            % Correct start time for numDummys removed & convert to seconds
                            J.sess(r).cond((FoSEx-1)*numel(num_seq)+c).onset    = [R.startTimeReal(idx)/1000 - J.timing.RT*numDummys + announceTime + delay(s)];
                            J.sess(r).cond((FoSEx-1)*numel(num_seq)+c).duration = dur;                       % durations of task we are modeling (not length of entire trial)
                            
                            J.sess(r).cond((FoSEx-1)*numel(num_seq)+c).tmod     = 0;
                            J.sess(r).cond((FoSEx-1)*numel(num_seq)+c).orth     = 0;
                            J.sess(r).cond((FoSEx-1)*numel(num_seq)+c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            
                            % Do some subject info for fields in SPM_info.mat.
                            S.SN    		= s;
                            S.run   		= r;    % 1-8: functional runs
                            S.runAll        = (ss-1)*8 + r;  % 1-32
                            S.seqNumb 		= R.seqNumb(idx(1));
                            S.seqType    	= R.seqType(idx(1));
                            S.FoSEx         = R.FoSEx(idx(1));
                            S.isMetronome   = R.isMetronome(idx(1));
                            S.ScanSess      = R.ScanSess(idx(1));
                            T				= addstruct(T,S);
                        end;
                    end
                    % Add any additional regressors here.
                    J.sess(r).multi 	= {''};
                    J.sess(r).regress 	= struct('name', {}, 'val', {});
                    J.sess(r).multi_reg = {''};
                    % Define high pass filter cutoff (in seconds): see glm cases.
                    J.sess(r).hpf 		= hrf_cutoff;
                end;    % runs
                
                J.fact 			   = struct('name', {}, 'levels', {});
                J.bases.hrf.derivs = [0 0];
                J.bases.hrf.params = hrf_params;    % make it subject specific
                J.volt 			   = 1;
                J.global 		   = 'None';
                J.mask 	           = {fullfile(baseDir, 'imaging_data',subj_name{s}, 'rmask_noskull.nii,1')};
                J.mthresh 		   = 0.05;
                J.cvi_mask 		   = {fullfile(baseDir, 'imaging_data',subj_name{s},'rmask_gray.nii')};
                J.cvi 			   = cvi_type;
                % Save the GLM file for this subject.
                spm_rwls_run_fmri_spec(J);
                % Save the aux. information file (SPM_info.mat).
                % This file contains user-friendly information about the glm
                % model, regressor types, condition names, etc.
                save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
                
            end; % sn
        end; % sessN
    case 'GLM_estimate_RepSup'
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        vararginoptions(varargin,{'sn','sessN'});
        for ss = sessN
        for s = sn
            % Load files
            load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM.mat'));
            SPM.swd = fullfile(glmFoSExDir{ss},subj_name{s});
            % Run the GLM.
            spm_rwls_spm(SPM);
        end; % sn
        end; % sessN
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)
    case 'GLM_contrast_RepSup'  
     
        % enter sn
        % 1:   Trained seq (1st) vs. rest
        % 2:   Trained seq (2nd) vs. rest
        % 3:   Untrained seq (1st) vs. rest
        % 4:   Untrained seq (2nd) vs. rest
 
        vararginoptions(varargin,{'sn','sessN'});
        cwd = pwd;
        % Loop through subjects.
        for ss = sessN
        for s = sn
            glmSubjDir = [glmFoSExDir{ss} filesep subj_name{s}];
            cd(glmSubjDir);

            load SPM;
            SPM = rmfield(SPM,'xCon');
            T   = load('SPM_info.mat');

            %_____t contrast for trained seq 1st vs. rest
            con                                = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb<7 & T.FoSEx == 1) = 1;
            con                                = con/sum(con);
            SPM.xCon(1)                        = spm_FcUtil('Set',sprintf('TrainSeq_1st'), 'T', 'c',con',SPM.xX.xKXs);
            
            %_____t contrast for trained seq 2nd vs. rest
            con                                 = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb<7 & T.FoSEx == 2)   = 1;
            con                                 = con/sum(con);
            SPM.xCon(2)                         = spm_FcUtil('Set',sprintf('TrainSeq_2nd'), 'T', 'c',con',SPM.xX.xKXs);

            %_____t contrast for untrained seq 1st vs. rest
            con                                = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb>6 & T.FoSEx == 1) = 1;
            con                                = con/sum(con);
            SPM.xCon(3)                        = spm_FcUtil('Set',sprintf('UntrainSeq_1st'), 'T', 'c',con',SPM.xX.xKXs);
            
            %_____t contrast for untrained seq 2nd vs. rest
            con                                 = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb>6 & T.FoSEx == 2)   = 1;
            con                                 = con/sum(con);
            SPM.xCon(4)                         = spm_FcUtil('Set',sprintf('UntrainSeq_2nd'), 'T', 'c',con',SPM.xX.xKXs);
            %____do the constrasts
            SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat','SPM');
            
                % rename contrast images and spmT images
                conName = {'con','spmT'};
                for i=1:length(SPM.xCon),
                    for n=1:numel(conName),
                        oldName{i} = fullfile(glmSubjDir,sprintf('%s_%2.4d.nii',conName{n},i));
                        newName{i} = fullfile(glmSubjDir,sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                        movefile(oldName{i},newName{i});
                    end
                end
        end; % sn
        end; % sessN
        cd(cwd);
        
    case '3c_GLM_LOC'  %  ------- localizer GLM across the 4 sessions ------- % 
        % finger mapping 1-5
    case 'GLM_make_LOC'                                               % STEP 3.4a  :  Make the SPM.mat and SPM_info.mat files (prep the GLM - localizer)
        % localizer runs
        % makes the GLM file for each subject, and a corresponding 
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        
        glm = 2;    %1/2/3
        vararginoptions(varargin,{'sn','glm'});
        % Set some constants.
        prefix		 = 'u';
        T			 = [];
        dur			 = 2.5;                                                 % secs (length of task dur, not trial dur)
        delay        = [0.5 1 0 0];     
        announceTime = 1.0;                                                 % length of task announce time - currently not used
        % Gather appropriate GLM presets.
        switch glm
            case 1  % wls
                hrf_params = [5.5 12.5]; % change to 5.5 12.5
                hrf_cutoff = 128;
                cvi_type   = 'wls';
            case 2  % fast + hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = 128;
                cvi_type   = 'fast';
            case 3  % fast, no hpf
                hrf_params = [4.5 11]; % change to 5.5 12.5
                hrf_cutoff = inf;
                cvi_type   = 'fast';
        end

        % Loop through subjects and make SPM files.
        for s = sn
            D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));     
             
            % Do some subject structure fields.
            dircheck(fullfile(glmLocDir{glm}, subj_name{s}));
            J.dir 			 = {fullfile(glmLocDir{glm}, subj_name{s})};
            J.timing.units   = 'secs';                                      % timing unit that all timing in model will be
            J.timing.RT 	 = 1.0;                                         % TR (in seconds, as per 'J.timing.units')
            J.timing.fmri_t  = 16;
            J.timing.fmri_t0 = 1;
            
            % Loop through sessions
            for sessN = 1:4     % all 4 scanning sessions
                L = getrow(D,D.ScanSess==sessN);    % only blocks of that scan session
                uniqrun = unique(L.BN);
                % Loop through runs.
                for r = 1:numruns_loc_sess       % 2 localizer runs per session
                    
                    for i = 1:(numTRs(run_loc(r))-numDummys)                   % get nifti filenames, correcting for dummy scancs
                        
                        N{i} = [fullfile(baseDir, 'imaging_data',subj_name{s}, ...
                            [prefix subj_name{s},'_run',runs{sessN}{run_loc(1,r)},'.nii,',num2str(i)])];
                        
                    end;
                    J.sess((sessN-1)*2+r).scans = N;                                        % number of scans in run
                    % Loop through conditions.
                    
                    for c = 1:numel(num_fing)   % 5
                        R = getrow(L,L.BN==uniqrun(run_loc(1,r)));
                        idx						   = find(R.seqNumb==num_fing(c));    % find indx for fingers 1-5 (overall 13-17)
                        condName = sprintf('SeqNumb-%d',R.seqNumb(idx(1)));
                        J.sess((sessN-1)*2+r).cond(c).name 	   = condName;
                        % Correct start time for numDummys removed & convert to seconds
                        J.sess((sessN-1)*2+r).cond(c).onset    = [R.startTimeReal(idx)/1000 - J.timing.RT*numDummys + announceTime + delay(sn)];
                        J.sess((sessN-1)*2+r).cond(c).duration = dur;                       % durations of task we are modeling (not length of entire trial)
                        
                        J.sess((sessN-1)*2+r).cond(c).tmod     = 0;
                        J.sess((sessN-1)*2+r).cond(c).orth     = 0;
                        J.sess((sessN-1)*2+r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                        
                        % Do some subject info for fields in SPM_info.mat.
                        S.SN    		= s;
                        S.run           = (sessN-1)*2 + r;  % 1-8 overall func runs
                        S.seqNumb 		= R.seqNumb(idx(1));
                        S.seqType    	= R.seqType(idx(1));
                        S.isMetronome   = R.isMetronome(idx(1));
                        S.ScanSess      = R.ScanSess(idx(1));
                        T				= addstruct(T,S);
                    end;
                    
                    % Add any additional regressors here.
                    J.sess((sessN-1)*2+r).multi 	= {''};
                    J.sess((sessN-1)*2+r).regress 	= struct('name', {}, 'val', {});
                    J.sess((sessN-1)*2+r).multi_reg = {''};
                    % Define high pass filter cutoff (in seconds): see glm cases.
                    J.sess((sessN-1)*2+r).hpf 		= hrf_cutoff;
                end;    % runs
            end;        % scanning session
            
            J.fact 			   = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = hrf_params;    % make it subject specific
            J.volt 			   = 1;
            J.global 		   = 'None';
            J.mask 	           = {fullfile(baseDir, 'imaging_data',subj_name{s}, 'rmask_noskull.nii,1')};
            J.mthresh 		   = 0.05;
            J.cvi_mask 		   = {fullfile(baseDir, 'imaging_data',subj_name{s},'rmask_gray.nii')};
            J.cvi 			   = cvi_type;
            % Save the GLM file for this subject.
            spm_rwls_run_fmri_spec(J);        
            % Save the aux. information file (SPM_info.mat).
            % This file contains user-friendly information about the glm
            % model, regressor types, condition names, etc.
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            
        end;
    case 'GLM_estimate_LOC'                                           % STEP 3.5a  :  Run the GLM according to model defined by SPM.mat
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        vararginoptions(varargin,{'sn','glm'});
        glm=2;
        for s = sn
            % Load files
            load(fullfile(glmLocDir{glm},subj_name{s},'SPM.mat'));
            SPM.swd = fullfile(glmLocDir{glm},subj_name{s});
            % Run the GLM.
            spm_rwls_spm(SPM);
        end; % sn
        
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)    
    case 'GLM_contrast_LOC'                                           % STEP 3.6a  :  Make t-contsmlrasts - any / every finger vs. rest.
        % enter sn, glm #
        % 1:   Finger average vs. rest
        % 2-6: Single finger mapping (1-5)
 
        glm=2;
        vararginoptions(varargin,{'sn','glm'});
        cwd = pwd;
        % Loop through subjects.
        for s = sn
            glmSubjDir = [glmLocDir{glm} filesep subj_name{s}];
            cd(glmSubjDir);

            load SPM;
            SPM = rmfield(SPM,'xCon');
            T   = load('SPM_info.mat');

            %_____t contrast for single finger mapping (average) vs. rest
            con                = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb>12)= 1;
            con                = con/sum(con);
            SPM.xCon(1)        = spm_FcUtil('Set',sprintf('DigitAny'), 'T', 'c',con',SPM.xX.xKXs);
         
            %_____t contrast for single finger mapping vs. rest
            for d = 1:5
                con                = zeros(1,size(SPM.xX.X,2));
                con(:,T.seqNumb==d+12)  = 1;
                con                = con/sum(con);
                SPM.xCon(d+1)      = spm_FcUtil('Set',sprintf('Digit%d',d), 'T', 'c',con',SPM.xX.xKXs);
            end;
            
            %____do the constrasts
            SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat','SPM');
            
                % rename contrast images and spmT images
                conName = {'con','spmT'};
                for i=1:length(SPM.xCon),
                    for n=1:numel(conName),
                        oldName{i} = fullfile(glmSubjDir,sprintf('%s_%2.4d.nii',conName{n},i));
                        newName{i} = fullfile(glmSubjDir,sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                        movefile(oldName{i},newName{i});
                    end
                end
        end;
        cd(cwd);
    
    case '3d_GLM_LOC_sess'  % ------- construct GLM for localizer runs per session! ------- %
        % important for assessing consistencies of patterns etc.
    case 'GLM_LOC_sess_all'
        glm=2;
        vararginoptions(varargin,{'sn','glm','sessN'});
        sml1_imana_prep('GLM_make_LOC_sess','sn',sn,'glm',glm,'sessN',sessN);
        sml1_imana_prep('GLM_estimate_LOC_sess','sn',sn,'sessN',sessN);
        sml1_imana_prep('GLM_contrast_LOC_sess','sn',sn,'sessN',sessN);
    case 'GLM_make_LOC_sess'
        % localizer runs - per session
        % makes the GLM file for each subject, and a corresponding 
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        
        glm = 2;    %1/2/3
        sessN = 1;
        vararginoptions(varargin,{'sn','glm','sessN'});
        % Set some constants.
        prefix		 = 'u';
        T			 = [];
        dur			 = 2.5;                                                 % secs (length of task dur, not trial dur)
        delay        = [0.5 1 1 0.5 1 0 0];     
        announceTime = 0;                                                   % length of task announce time - currently not used
        % Gather appropriate GLM presets.
        switch glm
            case 1  % wls
                hrf_params = [5.5 12.5]; % change to 5.5 12.5
                hrf_cutoff = 128;
                cvi_type   = 'wls';
            case 2  % fast + hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = 128;
                cvi_type   = 'fast';
            case 3  % fast, no hpf
                hrf_params = [4.5 11]; % change to 5.5 12.5
                hrf_cutoff = inf;
                cvi_type   = 'fast';
        end

        % Loop through subjects and make SPM files.
        for s = sn
            D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));
          for ss = 1:sessN   
            T = [];
            % Do some subject structure fields.
            dircheck(fullfile(glmLocSessDir{ss}, subj_name{s}));
            J.dir 			 = {fullfile(glmLocSessDir{ss}, subj_name{s})};
            J.timing.units   = 'secs';                                      % timing unit that all timing in model will be
            J.timing.RT 	 = 1.0;                                         % TR (in seconds, as per 'J.timing.units')
            J.timing.fmri_t  = 16;
            J.timing.fmri_t0 = 1;
            
            L = getrow(D,D.ScanSess==ss);    % only blocks of that scan session
            %if (sn==2 & ss==3)
             %   uniqrun = [181,182,191,184,185,186,187,188,189,190];
            %elseif (sn==2 & ss==4)
             %   uniqrun = [192,193,194,195,202,197,198,199,200,201];
            %else
             %   uniqrun = unique(L.BN);
            %end
            % Loop through runs.
            for r = 1:numruns_loc_sess       % 2 localizer runs per session
                Rr = getrow(L,L.blockType==4); %blockType==4 - funct imaging localizer
                
                uniqrun=unique(Rr.BN);
                R = getrow(Rr,Rr.BN==uniqrun(r)); % 1-2 localizer runs
                
                for i = 1:(numTRs(run_loc(r))-numDummys)                   % get nifti filenames, correcting for dummy scancs
                    
                    N{i} = [fullfile(baseDir, 'imaging_data',subj_name{s}, ...
                        [prefix subj_name{s},'_run',runs{ss}{run_loc(1,r)},'.nii,',num2str(i)])];
                    
                end;
                J.sess(r).scans = N;                                        % number of scans in run
                % Loop through conditions.
                
                for c = 1:numel(num_fing)   % 5
                    idx						   = find(R.seqNumb==num_fing(c));    % find indx for fingers 1-5 (overall 13-17)
                    condName = sprintf('SeqNumb-%d',R.seqNumb(idx(1)));
                    J.sess(r).cond(c).name 	   = condName;
                    % Correct start time for numDummys removed & convert to seconds
                    J.sess(r).cond(c).onset    = [R.startTimeReal(idx)/1000 - J.timing.RT*numDummys + announceTime + delay(sn)];
                    J.sess(r).cond(c).duration = dur;                       % durations of task we are modeling (not length of entire trial)
                    
                    J.sess(r).cond(c).tmod     = 0;
                    J.sess(r).cond(c).orth     = 0;
                    J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    
                    % Do some subject info for fields in SPM_info.mat.
                    S.SN    		= s;
                    S.run           = r;
                    S.runAll        = (ss-1)*2 + r;  % 1-8 overall loc runs
                    S.seqNumb 		= R.seqNumb(idx(1))-12;     % 1-5
                    S.seqType    	= R.seqType(idx(1));
                    S.isMetronome   = R.isMetronome(idx(1));
                    S.ScanSess      = R.ScanSess(idx(1));
                    T				= addstruct(T,S);
                end;
                
                % Add any additional regressors here.
                J.sess(r).multi 	= {''};
                J.sess(r).regress 	= struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                % Define high pass filter cutoff (in seconds): see glm cases.
                J.sess(r).hpf 		= hrf_cutoff;
            end;    % runs
            
            J.fact 			   = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = hrf_params;    % make it subject specific
            J.volt 			   = 1;
            J.global 		   = 'None';
            J.mask 	           = {fullfile(baseDir, 'imaging_data',subj_name{s}, 'rmask_noskull.nii,1')};
            J.mthresh 		   = 0.05;
            J.cvi_mask 		   = {fullfile(baseDir, 'imaging_data',subj_name{s},'rmask_gray.nii')};
            J.cvi 			   = cvi_type;
            % Save the GLM file for this subject.
            spm_rwls_run_fmri_spec(J);        
            % Save the aux. information file (SPM_info.mat).
            % This file contains user-friendly information about the glm
            % model, regressor types, condition names, etc.
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
        end; % sessN
        end; % subject    
    case 'GLM_estimate_LOC_sess'
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        vararginoptions(varargin,{'sn','sessN'});
        for s = sn
            % Load files
            for ss = 1:sessN
                load(fullfile(glmLocSessDir{ss},subj_name{s},'SPM.mat'));
                SPM.swd = fullfile(glmLocSessDir{ss},subj_name{s});
                % Run the GLM.
                spm_rwls_spm(SPM);
            end; % session
        end; % subject
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)   
    case 'GLM_contrast_LOC_sess'
         % enter sn, glm #
        % 1:   Finger average vs. rest
        % 2-6: Single finger mapping (1-5)
 
        vararginoptions(varargin,{'sn','sessN'});
        cwd = pwd;
        % Loop through subjects.
        for s = sn
            for ss = 1:sessN
                glmSubjDir = [glmLocSessDir{ss} filesep subj_name{s}];
                cd(glmSubjDir);
                
                load SPM;
                SPM = rmfield(SPM,'xCon');
                T   = load('SPM_info.mat');
                
                %_____t contrast for single finger mapping (average) vs. rest
                con                = zeros(1,size(SPM.xX.X,2));
                con(:,T.seqNumb>0)= 1;
                con                = con/sum(con);
                SPM.xCon(1)        = spm_FcUtil('Set',sprintf('DigitAny'), 'T', 'c',con',SPM.xX.xKXs);
                
                %_____t contrast for single finger mapping vs. rest
                for d = 1:5
                    con                = zeros(1,size(SPM.xX.X,2));
                    con(:,T.seqNumb==d)  = 1;
                    con                = con/sum(con);
                    SPM.xCon(d+1)      = spm_FcUtil('Set',sprintf('Digit%d',d), 'T', 'c',con',SPM.xX.xKXs);
                end;
                
                %____do the constrasts
                SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
                save('SPM.mat','SPM');
                
                % rename contrast images and spmT images
                conName = {'con','spmT'};
                for i=1:length(SPM.xCon),
                    for n=1:numel(conName),
                        oldName{i} = fullfile(glmSubjDir,sprintf('%s_%2.4d.nii',conName{n},i));
                        newName{i} = fullfile(glmSubjDir,sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                        movefile(oldName{i},newName{i});
                    end
                end
                clear SPM;
            end; 
        end; % subject
        cd(cwd);
    
    case '3e_GLM_trialwise' % ---------- construct GLM estimating beta for each trial separately --------
     % fits regressor for each trial   
    case 'GLM_trial_all'
        glm=2;
        vararginoptions(varargin,{'sn','sessN'});
        sml1_imana_prep('GLM_make_trial','sn',sn,'sessN',sessN);
        sml1_imana_prep('GLM_estimate_trial','sn',sn,'sessN',sessN);
       % sml1_imana('GLM_contrast_trial','sn',sn,'sessN',sessN);
    case 'GLM_make_trial'
        % makes the GLM file for each subject, and a corresponding 
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        
        glm = 2;    %1/2/3
        vararginoptions(varargin,{'sn','glm','sessN'});
        % Set some constants.
        prefix		 = 'u';
        T			 = [];
        dur			 = 2.5;                                                 % secs (length of task dur, not trial dur)
        % adjusting hrf per subject & session based on extracted timeseries!  
        delay     = [0.5 1 1 0.5 1 0 0 0.5 1 1 0.5 1 1];  

        announceTime = 0;                                                 % length of task announce time - currently not used
        % Gather appropriate GLM presets.
        switch glm
            case 1  % wls
                hrf_params = [5.5 12.5]; 
                hrf_cutoff = 128;
                cvi_type   = 'wls';
            case 2  % fast + hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = 128;
                cvi_type   = 'fast';
            case 3  % fast, no hpf
                hrf_params = [5.5 12.5];
                hrf_cutoff = inf;
                cvi_type   = 'fast';
        end

        % Loop through subjects / all sessions and make SPM files.
            for s = sn
                D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));
                
                % Do some subject structure fields.
                dircheck(fullfile(glmTrialDir{sessN}, subj_name{s}));
                J.dir 			 = {fullfile(glmTrialDir{sessN}, subj_name{s})};
                J.timing.units   = 'secs';                                      % timing unit that all timing in model will be
                J.timing.RT 	 = 1.0;                                         % TR (in seconds, as per 'J.timing.units')
                J.timing.fmri_t  = 16;
                J.timing.fmri_t0 = 1;
                
                L = getrow(D,D.ScanSess==sessN);    % only blocks of that scan session

                % Loop through runs.
                for r = 1:numruns_task_sess
                    if sessN == 4
                        Rr = getrow(L,L.blockType==9); %blockType==9 - func imaging run without metronome
                    else
                        Rr = getrow(L,L.blockType==3); %blockType==3 - funct imaging run with metornome
                    end
                    
                    uniqrun=unique(Rr.BN);
                    R = getrow(Rr,Rr.BN==uniqrun(r)); % 1-8 func runs of the session

                    for i = 1:(numTRs(run_task(r))-numDummys)                   % get nifti filenames, correcting for dummy scancs
                        
                        N{i} = [fullfile(baseDir, 'imaging_data',subj_name{s}, ...
                            [prefix subj_name{s},'_run',runs{sessN}{run_task(1,r)},'.nii,',num2str(i)])];
                        
                    end;
                    J.sess(r).scans = N;                                        % number of scans in run
                    % Loop through conditions.
                    
                    for c = 1:size(R.TN,1) % all trials
                        idx						   = c;             % find indx of all trials in run - 1:6 trained; 7-12 untrained
                        condName = sprintf('Trial number-%d',c);
                        J.sess(r).cond(c).name 	   = condName;
                        % Correct start time for numDummys removed & convert to seconds
                        J.sess(r).cond(c).onset    = [R.startTimeReal(idx)/1000 - J.timing.RT*numDummys + announceTime + delay(sn)];
                        J.sess(r).cond(c).duration = dur;                       % durations of task we are modeling (not length of entire trial)
                        
                        J.sess(r).cond(c).tmod     = 0;
                        J.sess(r).cond(c).orth     = 0;
                        J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                        
                        % Do some subject info for fields in SPM_info.mat.
                        S.SN    		= s;
                        S.run   		= r;    % 1-8: functional runs
                        S.runAll        = (sessN-1)*8 + r;  % 1-32
                        S.seqNumb 		= R.seqNumb(idx);
                        S.seqType    	= R.seqType(idx);
                        S.MT            = R.MT(idx);
                        S.isError       = R.isError(idx);
                        S.isMetronome   = R.isMetronome(idx);
                        S.ScanSess      = R.ScanSess(idx);
                        T				= addstruct(T,S);
                    end;
                    
                    % Add any additional regressors here.
                    J.sess(r).multi 	= {''};
                    J.sess(r).regress 	= struct('name', {}, 'val', {});
                    J.sess(r).multi_reg = {''};
                    % Define high pass filter cutoff (in seconds): see glm cases.
                    J.sess(r).hpf 		= hrf_cutoff;
                end;
                J.fact 			   = struct('name', {}, 'levels', {});
                J.bases.hrf.derivs = [0 0];
                J.bases.hrf.params = hrf_params;    % make it subject specific
                J.volt 			   = 1;
                J.global 		   = 'None';
                J.mask 	           = {fullfile(baseDir, 'imaging_data',subj_name{s}, 'rmask_noskull.nii,1')};
                J.mthresh 		   = 0.05;
                J.cvi_mask 		   = {fullfile(baseDir, 'imaging_data',subj_name{s},'rmask_gray.nii')};
                J.cvi 			   = cvi_type;
                % Save the GLM file for this subject.
                spm_rwls_run_fmri_spec(J);
                % Save the aux. information file (SPM_info.mat).
                % This file contains user-friendly information about the glm
                % model, regressor types, condition names, etc.
                save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
                
            end;    % sn
    case 'GLM_estimate_trial'
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        vararginoptions(varargin,{'sn','sessN'});
        
            for s = sn
                % Load files
                load(fullfile(glmTrialDir{sessN},subj_name{s},'SPM.mat'));
                SPM.swd = fullfile(glmTrialDir{sessN},subj_name{s});
                % Run the GLM.
                spm_rwls_spm(SPM);
            end; % sn
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)
    case 'GLM_contrast_trial'   % DEPRECIATED - TO DO
        % enter sn, sessN #
        % 
        
        vararginoptions(varargin,{'sn','sessN'});
        cwd = pwd;

            % Loop through subjects.
            for s = sn
                glmSubjDir = [glmSessDir{sessN} filesep subj_name{s}];
                cd(glmSubjDir);
                
                load SPM;
                SPM = rmfield(SPM,'xCon');
                T   = load('SPM_info.mat');
                
                nrun    = numel(SPM.nscan);
                tt      = repmat([1:numel(num_seq)],1,nrun);
                
                % (1) t contrasts each sequence against rest
                % (1:12)
                for c = 1:numel(num_seq)
                    con = zeros(1,size(SPM.xX.X,2));
                    con(tt==c) = 1;
                    con = con/sum(con);
                    SPM.xCon(c) = spm_FcUtil('Set',sprintf('Seq%d',c), 'T', 'c',con',SPM.xX.xKXs);
                end
                
                % (2) t contrast for trained seq vs. rest
                con                = zeros(1,size(SPM.xX.X,2));
                con(:,T.seqNumb<7) = 1;
                con                = con/sum(con);
                SPM.xCon(end+1)    = spm_FcUtil('Set',sprintf('TrainSeq'), 'T', 'c',con',SPM.xX.xKXs);
                
                % (3) t contrast for novel seq vs. rest
                con                = zeros(1,size(SPM.xX.X,2));
                con(:,T.seqNumb>6 & T.seqNumb<13)  = 1;
                con                = con/sum(con);
                SPM.xCon(end+1)    = spm_FcUtil('Set',sprintf('UntrainSeq'), 'T', 'c',con',SPM.xX.xKXs);
                
                %____do the constrasts
                SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
                save('SPM.mat','SPM');
                
                % rename contrast images and spmT images
                conName = {'con','spmT'};
                for i=1:length(SPM.xCon),
                    for n=1:numel(conName),
                        oldName{i} = fullfile(glmSubjDir,sprintf('%s_%2.4d.nii',conName{n},i));
                        newName{i} = fullfile(glmSubjDir,sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                        movefile(oldName{i},newName{i});
                    end
                end
            end; % sn

        cd(cwd); 
        
    case '4_ROI'% ------------------- ROI analysis. Expand for more info. ------
                % The ROI cases are used to:
        %       - map ROIs to each subject
        %       - harvest timeseries from each roi for each condition
        %
        % There is no 'processAll' case here. However, the following cases
        % must be called to utilize other roi cases:
        %       'ROI_define'      :  Maps rois to surface of each subject-
        %                             requires paint files from above case.
        %       'ROI_timeseries'  :  Only if you wish to plot timeseries
        %       'ROI_getBetas'    :  Harvest patterns from roi
        %       'ROI_stats'       :  Estimate distances, etc. 
        %
        % You can view roi maps by loading paint files and subject surfaces
        % in Caret (software).
        % 
        % Most functionality is achieved with rsa toolbox by JDiedrichsen
        % and NEjaz (among others).
        %
        % See blurbs in each SEARCH case to understand what they do.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
        
    case 'SURF_162tessels'

        for h=1:length(hemName) % per hemisphere  
            % load in metric file
            F=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.metric',hem{h})));
            B=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.cerebral_cortex.paint',hem{h})));
            
            % get medial wall
            %                 find(F.data(B.data==6))=0;
            F.data(B.data==6)=163;
            
            I=setdiff([1:158],unique(F.data));
            N=[159,160,162,163];
            
            % relabel data (medial wall is now index 149)
            for i=1:length(N),
                F.data(F.data==N(i))=I(i);
            end
            
            % save out new metric file (medial wall removed)
            caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162_new.metric',hem{h})),F);
            
            % make and save new paintfile
            F=rmfield(F,'column_color_mapping');
            
            % loop over tessels
            for ii=1:length(unique(F.data)),
                F.paintnames{ii}=sprintf('tessel%2.2d',ii);
            end
            %                 F.paintnames{I(end)}='medial_wall';
            F.num_paintnames=length(unique(F.data));
            
            caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h})),F);
            
            % make new area colour file
            cmap=[1:1:255];
            for c=1:length(unique(F.data)),
                M.data(c,:)=randsample(cmap,3);
            end
            M.data(I(end),:)=[0 0 0];
            M.encoding={'BINARY'};
            M.column_name={'Column_01'};
            M.num_rows=F.num_cols;
            M.num_cols=F.num_cols;
            M.column_color_mapping=repmat([-5 5],length(unique(F.data)),1);
            M.paintnames=F.paintnames;
            caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('tessel162.areacolor')),M);
        end
        
    case 'MASK_combine'
        % combine glm masks of all runs, create one mask overall
        
        vararginoptions(varargin,{'sn'});
        for s = sn
            for sess = 1:sess_sn(s)
                file = fullfile(glmSessDir{sess},subj_name{s},'mask.nii');
                V = spm_vol(file);
                mask(:,:,:,sess)=spm_read_vols(V);
            end
            maskNew=mask(:,:,:,1).*mask(:,:,:,2).*mask(:,:,:,3).*mask(:,:,:,4);
            dest = fullfile(regDir, ['mask_' subj_name{s} '.nii']);
            V.fname=dest;
            V.dim=size(maskNew);
            V=rmfield(V,'descrip');
            V.descrip = 'combined_mask';
            spm_write_vol(V,maskNew);
            fprintf('Combined mask definition done: %s\n',subj_name{s});
        end
    case 'ROI_define_sess'  % only used if not complete dataset - mask per glm rather than overall mask
        
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','sessN','parcelType'});
        for s=sn
            R=[];
            for h=1:2
                C = caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},['ROI.paint']));
                caretSubjDir = fullfile(caretDir,['x' subj_name{s}]);  
               % suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{s});

               %file = fullfile(regDir,['mask_' subj_name{s} '.nii']); % mask constructed in mask_combine
               file = fullfile(glmSessDir{sessN},subj_name{s},'mask.nii'); % mask constructed in glm    
                
                 for i=1:numregions_surf
                    R{i+(h-1)*numregions_surf}.type='surf_nodes';
                    R{i+(h-1)*numregions_surf}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                    R{i+(h-1)*numregions_surf}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                    R{i+(h-1)*numregions_surf}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                    R{i+(h-1)*numregions_surf}.flat=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                    
                    R{i+(h-1)*numregions_surf}.linedef=[10,0,1];
                    R{i+(h-1)*numregions_surf}.image=file;
                    R{i+(h-1)*numregions_surf}.name=[subj_name{s} '_' regname_cortex{i} '_' hem{h}];
                    R{i+(h-1)*numregions_surf}.location=find(C.data(:,1)==i);
                 end            
                          
            end;
                 
            R=region_calcregions(R);
            dircheck(regDir);
            cd(regDir);
            save(fullfile(regDir,[subj_name{s},sprintf('_%s_regions.mat',parcelType)]),'R');
            
            fprintf('\nROIs have been defined for %s \n',subj_name{sn});
        end    
    case 'ROI_define'                                                       % STEP 5.1   :  Define ROIs
        sn=[1:4];
        parcelType='162tessels'; % Brodmann or 162tessels
        vararginoptions(varargin,{'sn','parcelType'});
        for s=sn
            R=[];
            idx=0;
            for h=1:2
                caretSubjDir = fullfile(caretDir,['x' subj_name{s}]);
                file = fullfile(regDir,['mask_' subj_name{s} '.nii']); % mask constructed in mask_combine
                switch (parcelType)
                    case 'Brodmann'
                        C = caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},['ROI.paint']));  
                        for i=1:numregions_surf
                            R{i+(h-1)*numregions_surf}.type='surf_nodes';
                            R{i+(h-1)*numregions_surf}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            R{i+(h-1)*numregions_surf}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            R{i+(h-1)*numregions_surf}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            R{i+(h-1)*numregions_surf}.flat=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                            
                            R{i+(h-1)*numregions_surf}.linedef=[10,0,1];
                            R{i+(h-1)*numregions_surf}.image=file;
                            R{i+(h-1)*numregions_surf}.name=[subj_name{s} '_' regname_cortex{i} '_' hem{h}];
                            R{i+(h-1)*numregions_surf}.location=find(C.data(:,1)==i);
                        end
                    case '162tessels'
                        C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
                        tesselNames=unique(C.data);
                        if h==2,
                            idx=numel(C.paintnames); % append to the first hemisphere indices
                        end
                        for i=1:numel(C.paintnames),
                            R{i+idx}.type='surf_nodes';
                            R{i+idx}.location=find(C.data(:,1)==tesselNames(i));
                            R{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            R{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            R{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                          %  R{i+idx}.flat=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                            R{i+idx}.linedef=[5,0,1];
                            R{i+idx}.image=file;
                            R{i+idx}.file=file;
                            R{i+idx}.name=[C.paintnames{i}];
                        end;
                        R{149}.name='medial_wall'; % tessel 149 is the medial wall
                        R{149*2}.name='medial_wall'; % tessel 149*2 is the other medial wall
                    case 'cortex_buckner'
                        C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},sprintf('%s.Yeo7.paint',hem{h}))); % freesurfer
                        tesselNames=unique(C.data);
                        if h==2,
                            idx=C.num_paintnames-1; % append to the first hemisphere indices
                        end
                        for i=1:(C.num_paintnames-1),% no medial wall: label=1
                            R{i+idx}.type='surf_nodes';
                            R{i+idx}.location=find(C.data(:,1)==tesselNames(i+1));
                            R{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            R{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            R{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                          %  R{i+idx}.flat=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                            R{i+idx}.linedef=[5,0,1];
                            R{i+idx}.image=file;
                            R{i+idx}.file=file;
                            R{i+idx}.name=sprintf('%s-%s',C.paintnames{i+1},hem{h});
                        end;
                end
            end;
            
            R=region_calcregions(R);
            
            dircheck(regDir);
            save(fullfile(regDir,[subj_name{s},sprintf('_%s_regions.mat',parcelType)]),'R');
            
            fprintf('\nROIs %s have been defined for %s \n',parcelType,subj_name{s});
        end
    case 'ROI_make_nii'                                                     % OPTIONAL   :  Convert ROI def (.mat) into multiple .nii files (to check!)
           parcelType='Brodmann';
           vararginoptions(varargin,{'sn','sessN','parcelType'});
           
        for s=sn
            % load ROI definition
            load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},parcelType)));
            % loop over rois
            for roi = 1:size(R,2)
                % mask volume
                mask = fullfile(regDir,sprintf('mask_%s.nii',subj_name{s}));           
                % Save region file as nifti
                cd(regDir);
                region_saveasimg(R{roi},mask);      
            end
            
        end
    case 'ROI_extract_coords'
     % extracts coordinate positions (centroids) of regions of interest
     % useful for plotting positions of different ROIs (connectivity
     % analysis)
     sn = [4:9,11:31];
     parcelType='Brodmann';
     vararginoptions(varargin,{'sn','parcelType'});
     
     CC=[];
     for s=sn
         load(fullfile(regDir,sprintf('%s_%s_regions',subj_name{s},parcelType)));
         for i=1:size(R,2)
             if strcmp(parcelType,'Brodmann');
                C.flatCentr = mean(R{i}.flatcoord,1);
             end
             C.volCentr  = mean(R{i}.data,1);
             C.roi       = i;
             C.hemi      = (i>size(R,2)/2)+1;
             C.regType   = i-((C.hemi-1)*size(R,2)/2);
             C.sn        = s;
             CC = addstruct(CC,C);
         end
         fprintf('Done extracting centroids - %s\n',subj_name{s});
     end
     if strcmp(parcelType,'Brodmann')
         CC.flatCentr(:,3)=[];
         G.flatCentr(:,1) = nanmean(pivottable(CC.sn,CC.roi,CC.flatCentr(:,1),'nanmean'),1)';
         G.flatCentr(:,2) = nanmean(pivottable(CC.sn,CC.roi,CC.flatCentr(:,2),'nanmean'),1)';
     end
     G.volCentr(:,1)  = nanmean(pivottable(CC.sn,CC.roi,CC.volCentr(:,1),'nanmean'),1)';
     G.volCentr(:,2)  = nanmean(pivottable(CC.sn,CC.roi,CC.volCentr(:,2),'nanmean'),1)';
     G.volCentr(:,2)  = nanmean(pivottable(CC.sn,CC.roi,CC.volCentr(:,3),'nanmean'),1)';
     G.roi = [1:size(R,2)]';
     G.hemi = (G.roi>size(R,2)/2)+1;
     G.regType = G.roi-((G.hemi-1)*size(R,2)/2);

     save(fullfile(regDir,sprintf('region_%s_centroids_individSubj',parcelType)),'-struct','CC');
     save(fullfile(regDir,sprintf('region_%s_centroids_group',parcelType)),'-struct','G');
     
    case 'ROI_timeseries'                                                   % STEP 5.2   :  Extract onsets and events/trials - hrf
        % to check the model quality of the glm
        glm=2;
        parcelType='Brodmann';
        sessN=1;
        vararginoptions(varargin,{'sn','sessN','glm','regType'});

        pre=4;          % How many TRs before the trial onset
        post=16;        % How many TRs after the trial onset
        T=[];
        for ss=sessN
            for s=sn
                fprintf('Extracting the onsets and events for subject %s, glm %d and session %d\n',subj_name{s},glm,ss);
                load(fullfile(glmSessDir{ss},subj_name{s},'SPM.mat'));
                
                SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s}));
                load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},parcelType)));      % This is made in case 'ROI_define'
                [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
                
                % Create a structure with trial onset and trial type (event)
                D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
                % D.event - conditions (seqNumb: 1-12)
                % D.block - run
                for r=1:size(y_raw,2)   % regions
                    S.block=D.block;
                    for i=1:(size(S.block,1));
                        S.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    end;
                    S.event=D.event;
                    S.sn=ones(length(S.event),1)*s;
                    S.region=ones(length(S.event),1)*r; % region
                    S.regType=regType(S.region)';
                    S.regSide=regSide(S.region)';
                    % indicate sequence type
                    S.seqType=zeros(length(S.event),1);
                    S.seqType(S.event<7)=1;
                    S.seqType(S.event>6)=2;
                    T=addstruct(T,S);
                end;
                cd(regDir);
                save(sprintf('hrf_%s_sess%d_glm%d.mat',subj_name{s},ss,glm),'-struct','T');
            end;
        end
    case 'ROI_plot_timeseries'                                              % STEP 5.3   :  Plot estimated hrf by sequence type (trained/untrained)                     
        sn=1; 
        glm=2; 
        sessN=1;
        reg=5;
        regS=1; % 1 - LH, 2 - RH

        vararginoptions(varargin,{'sn','sessN','glm','reg','regS','BG'});
        T=load(fullfile(regDir,sprintf('hrf_%s_sess%d_glm%d.mat',subj_name{sn},sessN,glm)));
        
        traceplot([-4:16],T.y_adj,'errorfcn','stderr','subset',T.regSide==regS & T.regType==reg,'split',T.seqType,'leg','auto');
        hold on;
        traceplot([-4:16],T.y_hat,'subset',T.regSide==regS & T.regType==reg,'linestyle',':','split',T.seqType);
        hold off;
        xlabel('TR');
        ylabel('activation');
        drawline(0); 
    case 'ROI_timeseries_FoSEx'                                                  
        % to check the model quality of the glm
        glm=2;
        parcelType='Brodmann';
        sessN=1;
        vararginoptions(varargin,{'sn','sessN','glm','parcelType','sessN'});
        
        pre=4;          % How many TRs before the trial onset
        post=16;        % How many TRs after the trial onset
        T=[];
        for ss=sessN
            for s=sn
                fprintf('Extracting the onsets and events for subject %s, glm %d and session %d\n',subj_name{s},glm,sessN);
                load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM.mat'));
                
                SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s}));
                load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},parcelType)));      % This is made in case 'ROI_define'
                [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
                
                % Create a structure with trial onset and trial type (event)
                D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
                % D.event - conditions (seqNumb: 1-12)
                % D.block - run
                for r=1:size(y_raw,2)   % regions
                    S.block=D.block;
                    for i=1:(size(S.block,1));
                        S.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                        S.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    end;
                    S.event=D.event;
                    S.sn=ones(length(S.event),1)*s;
                    S.region=ones(length(S.event),1)*r; % region
                    S.regType=regType(S.region)';
                    S.regSide=regSide(S.region)';
                    % indicate first / second execution
                    S.FoSEx=zeros(length(S.event),1);
                    S.FoSEx(S.event<12)=1;
                    S.FoSEx(S.event>6)=2;
                    T=addstruct(T,S);
                end;
                cd(regDir);
                save(sprintf('hrf_%s_sess%d_FoSEx.mat',subj_name{s},ss),'-struct','T');
            end;
        end
    case 'ROI_plot_timeseries_FoSEx'                
        sn=1; 
        sessN=1;
        reg=5;
        regS=1; % 1 - LH, 2 - RH

        vararginoptions(varargin,{'sn','sessN','glm','reg','regS','BG'});
        T=load(fullfile(regDir,sprintf('hrf_%s_sess%d_FoSEx.mat',subj_name{sn},sessN)));
        
        traceplot([-4:16],T.y_adj,'errorfcn','stderr','subset',T.regSide==regS & T.regType==reg,'split',T.FoSEx,'leg','auto');
        hold on;
        traceplot([-4:16],T.y_hat,'subset',T.regSide==regS & T.regType==reg,'linestyle',':','split',T.FoSEx);
        hold off;
        xlabel('TR');
        ylabel('activation');
        drawline(0); 
    case 'ROI_timeseries_localizer'                                         % STEP 5.4   :  Extract onsets and events/trials for localizer - hrf
        % to check the model quality of the glm
        sessN=1;
        parcelType='Brodmann';
        vararginoptions(varargin,{'sn','sessN','parcelType'});

        pre=4;          % How many TRs before the trial onset
        post=16;        % How many TRs after the trial onset
        T=[];
        for s=sn
            fprintf('Extracting the onsets and events for subject %s and session %d\n',subj_name{s},sessN);
            load(fullfile(glmLocSessDir{sessN},subj_name{s},'SPM.mat'));
            
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s})); % This accounts for shifting
            load(fullfile(regDir,fullfile('%s_%s_regions.mat',subj_name{s},parcelType)));   % This is made in case 'ROI_define'
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
            
            % Create a structure with trial onset and trial type (event)
            D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
            % D.event - conditions (seqNumb: 1-10)
            % D.block - run
            for r=1:size(y_raw,2)   % regions
                for i=1:size(D.block,1);  % runs
                    D.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                end;                
                D.sn=ones(size(D.event,1),1)*s;
                D.region=ones(size(D.event,1),1)*r; % region
                D.regType=regType(D.region)';
                D.regSide=regSide(D.region)';
                % indicate sequence type
                D.seqType=ones(size(D.event,1),1)*3; % localiser mapping - 3
                T=addstruct(T,D);
            end;
        end;
        cd(regDir);
        save(sprintf('hrf_%s_LOC_sess%d.mat',subj_name{sn},sessN),'-struct','T');
    case 'ROI_plot_timeseries_localizer'                                    % STEP 5.5   :  Plot estimated hrf for finger presses                    
        sn=7; 
        sessN=1; 
        reg=2;
        regS=1; % 1 - LH, 2 - RH
        vararginoptions(varargin,{'sn','sessN','reg','regS'});
        T=load(fullfile(regDir,sprintf('hrf_%s_LOC_sess%d.mat',subj_name{sn},sessN)));
        
        traceplot([-4:16],T.y_adj,'errorfcn','stderr','subset',T.regSide==regS & T.regType==reg,'leg','auto');
        hold on;
        traceplot([-4:16],T.y_hat,'subset',T.regSide==regS & T.regType==reg,'linestyle',':');
        hold off;
        xlabel('TR');
        ylabel('activation');
        drawline(0);
    
    case 'ROI_define_BG'
        vararginoptions(varargin,{'sn'});
        for s=sn
            R=[];
            for h=1:2
                C = caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},['ROI.paint']));
                
                file = fullfile(regDir,['mask_' subj_name{s} '.nii']); % mask constructed in mask_combine
                %file = fullfile(glmSessDir{sessN},subj_name{s},'mask.nii');
                
                for j=1:numregions_BG
                    % Get basal ganglia
                    fileBG = fullfile(BGDir,'FSL',subj_name{s},sprintf('%s_%s_%s.nii', subj_name{s},regname_BG{j}, hem{h}));
                    R{j+(h-1)*numregions_BG}.type = 'roi_image';
                    R{j+(h-1)*numregions_BG}.file= fileBG;
                    R{j+(h-1)*numregions_BG}.name = [subj_name{s} '_' regname_BG{j} '_' hem{h}];
                    R{j+(h-1)*numregions_BG}.value = 1;
                end
            end;
            
            R=region_calcregions(R);
            dircheck(regDir);
            cd(regDir);
            save([subj_name{s} '_BG_regions.mat'],'R');
            
            fprintf('\nBG ROIs have been defined for %s \n',subj_name{s});
        end
    case 'ROI_BG_timeseries'
        glm=2;
        vararginoptions(varargin,{'sn','sessN','glm'});

        pre=4;          % How many TRs before the trial onset
        post=16;        % How many TRs after the trial onset
        T=[];
        for s=sn
            fprintf('Extracting the onsets and events for subject %s, glm %d and session %d\n',subj_name{s},glm,sessN);
            load(fullfile(glmSessDir{sessN},subj_name{s},'SPM.mat')); 
            
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s}));
            load(fullfile(regDir,[subj_name{s} '_BG_regions.mat']));      % This is made in case 'ROI_define'
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
            
            % Create a structure with trial onset and trial type (event)
            D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
            % D.event - conditions (seqNumb: 1-12)
            % D.block - run
            for r=1:size(y_raw,2)   % regions
                S.block=D.block;
                for i=1:(size(S.block,1));
                    S.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    S.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    S.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    S.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                end;
                S.event=D.event;
                S.sn=ones(length(S.event),1)*s;
                S.region=ones(length(S.event),1)*r; % region
                % indicate sequence type
                S.seqType=zeros(length(S.event),1);
                S.seqType(S.event<7)=1;
                S.seqType(S.event>6)=2;
                T=addstruct(T,S);
            end;
            T.regType=T.region;
            T.regType(T.region>4)=T.regType(T.region>4)-4;
            T.regSide=ones(size(T.regType));
            T.regSide(T.region<5)=1;
            T.regSide(T.region>4)=2;
            cd(regDir);
            save(sprintf('hrf_BG_%s_sess%d.mat',subj_name{s},sessN),'-struct','T');
        end;
    case 'ROI_BG_plot_timeseries'
        sn=4; 
        sessN=1;
        reg=1;
        regS=1; % 1 - LH, 2 - RH

        vararginoptions(varargin,{'sn','sessN','glm','reg','regS','BG'});
        T=load(fullfile(regDir,sprintf('hrf_BG_%s_sess%d.mat',subj_name{sn},sessN)));
        
        traceplot([-4:16],T.y_adj,'errorfcn','stderr','subset',T.regSide==regS & T.regType==reg,'split',T.seqType,'leg','auto');
        hold on;
        traceplot([-4:16],T.y_hat,'subset',T.regSide==regS & T.regType==reg,'linestyle',':','split',T.seqType);
        hold off;
        xlabel('TR');
        ylabel('activation');
        drawline(0); 
 
        
    case '5_SUIT' % ------------- BG and SUIT preparation. Expand for more info. CURRENTLY NOT IN USE ------
        % The ROI cases are used to:
        %       - map ROIs to each subject
        %       - harvest timeseries from each roi for each condition
        %
        % There is no 'processAll' case here. However, the following cases
        % must be called to utilize other roi cases:
        %       'ROI_makePaint'   :  Creates roi paint files (see case)
        %       'ROI_define'      :  Maps rois to surface of each subject-
        %                             requires paint files from above case.
        %       'ROI_timeseries'  :  Only if you wish to plot timeseries
        %       'ROI_getBetas'    :  Harvest patterns from roi
        %       'ROI_stats'       :  Estimate distances, etc. 
        %                               This is the big kahuna as it is
        %                               often loaded by future cases.
        %
        % You can view roi maps by loading paint files and subject surfaces
        % in Caret (software).
        % 
        % Most functionality is achieved with rsa toolbox by JDiedrichsen
        % and NEjaz (among others).
        %
        % See blurbs in each SEARCH case to understand what they do.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
                                                      % STEP 4.1-2  :  Do segmentation for BG using FSL (step 1&2)
        % uses FSL to do segmentation for BG and map the regions to
        % individual subject anatomy
        % need to run MATLAB through terminal 
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        
        % 1) Run step 1: sml1_imana_prep('BG_FSLmask',1,1)
        % 2) Open the zip folder p01b_BG_all_fast_firstseg.nii.gz
        % 3) Run step 2
        
        sn = varargin{1};
        step = varargin{2};
        
        switch (step)
            case 1 % run FSL routine                
                for s= sn%1:numel(subj_name)
                    IN= fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']);
                    outDir = fullfile(baseDir, 'basal_ganglia', 'FSL');
                    dircheck(outDir);
                    OUT= fullfile(outDir, subj_name{s}, [subj_name{s}, '_BG.nii']);
                    %calc with FSL
                    comm=sprintf('run_first_all -i %s -o %s', IN, OUT);
                    fprintf('%s\n',comm);
                    system(comm);
                end
            case 2 % make the ROI images in subject space
                %         10 Left-Thalamus-Proper 40
                %         11 Left-Caudate 30
                %         12 Left-Putamen 40
                %         13 Left-Pallidum 40
                %         49 Right-Thalamus-Proper 40
                %         50 Right-Caudate 30
                %         51 Right-Putamen 40
                %         52 Right-Pallidum 40
                BGnumber= [11 13 12 10; 50 52 51 49];
                %'CaudateN' 'Pallidum', 'Putamen' 'Thalamus'
               
                for s= sn%1:numel(subj_name)
                    %----deform info for basal ganglia ROI to individual space
                    nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_to_std_sub.mat']);
                    mniDir = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', subj_name{s});
                    if ~exist(mniDir,'dir')
                        mkdir(mniDir);
                    end
                    
                    for h=1:2
                        for i=1:numregions_BG
                            fprintf('Working on subj: %i region: %s \n', s, [regname{i+numregions_surf},'_',hem{h}])
                            
                            %----get names!
                            IN= fullfile(baseDir, 'basal_ganglia', 'FSL', subj_name{s},...
                                [subj_name{s},'_BG_all_fast_firstseg.nii']);
                            
                            OUT{i}= fullfile(baseDir, 'basal_ganglia', 'FSL', subj_name{s},...
                                [subj_name{s},'_',regname{i+numregions_surf},'_',hem{h},'.nii']);
                            
                            OUT_MNI{i}= fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', subj_name{s},...
                                [subj_name{s},'_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                            
                            %----make subj specific ROI image
                            spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                        end
                        %----do deformation
                        % spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end
            case 3 %make the avrg mask image
                for h=1:2
                    for i=1:numregions_BG
                        for s = 1:numel(sn)%1:numel(subj_name)
                            IN{s} = fullfile(baseDir, 'basal_ganglia', 'FSL',subj_name{sn(s)},...
                                [subj_name{sn(s)},'_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                        end
                        outDir = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', 'avrg');
                        if ~exist(outDir, 'dir');
                            mkdir(outDir);
                        end
                        OUT = fullfile(outDir,...
                            ['avrg_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                        spmj_imcalc_mtx(IN,OUT,'mean(X)');
                    end
                end

        end
    case 'SUIT_process_all'  
        % example: sml1_imana('SUIT_process_all',1,2)'
        sn = varargin{1}; % subjNum
        glm = varargin{2}; % glmNum
        %         spm fmri - call first
        for s=sn,
            sml1h_imana('SUIT_isolate_segment','sn',s);
            sml1h_imana('SUIT_normalize','sn',s);
            sml1h_imana('SUIT_reslice',s,glm,'betas');
            sml1h_imana('SUIT_reslice',s,glm,'contrast');
            %sml1h_imana('SUIT_reslice_localizer',s,glm,'betas');
            %sml1h_imana('SUIT_reslice_localizer',s,glm,'contrast');
            sml1h_imana('SUIT_make_mask','sn',s,'glm',glm);
            sml1h_imana('SUIT_roi','sn',s);
            fprintf('Suit data processed for %s',subj_name{s});
        end
    case 'SUIT_isolate_segment'                                             % STEP 4.3   :  Isolate cerebellum using SUIT
        vararginoptions(varargin,{'sn'});     
        for s=sn,
            suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{sn});
            dircheck(suitSubjDir);

            source = fullfile(anatomicalDir,subj_name{s},[subj_name{sn}, '_anatomical','.nii']);
            dest = fullfile(suitSubjDir,'anatomical.nii');
            copyfile(source,dest);
            cd(fullfile(suitSubjDir));
            suit_isolate_seg({fullfile(suitSubjDir,'anatomical.nii')},'keeptempfiles',1);
        end     
    case 'SUIT_normalize'                                                   % STEP 4.4   :  Put cerebellum into SUIT space
        vararginoptions(varargin,{'sn'});     
        for s=sn
            cd(fullfile(suitDir,'anatomicals',subj_name{sn}));
            %----run suit normalize
            suit_normalize(['c_anatomical.nii'], 'mask',['c_anatomical_pcereb.nii'])
        end;
    case 'SUIT_reslice'                                                     % STEP 4.5   :  Reslice functional images (glm as input)                      
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        type=varargin{3}; % 'betas','contrast','mask'     
        
        for s=sn
            suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{s});
            glmSubjDir = [glmDir{glm} filesep subj_name{s}];
            suitGLMDir = fullfile(suitDir,sprintf('glm%d',glm),subj_name{s});
            dircheck(suitGLMDir);
            mask   = fullfile(suitSubjDir,['c_anatomical_pcereb.nii']);
            params = fullfile(suitSubjDir,['mc_anatomical_snc.mat']);
            
            switch type
               case 'mask'
                    source=dir(fullfile(glmSubjDir,'mask.nii')); % images to be resliced
                case 'betas'
                    images='beta_0';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                case 'contrast'
                    images='con';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
            end
            for i=1:size(source,1)
                P=fullfile(glmSubjDir,source(i).name);
                O=fullfile(suitGLMDir,[source(i).name]);
                cd(fullfile(suitGLMDir));
                suit_reslice(P,params,'outfilename',O,'mask',mask,'interp',1,'vox',[1 1 1]);
            end;
        end;
    case 'SUIT_reslice_localizer'
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        type=varargin{3}; % 'betas','contrast','mask'     
        
        for s=sn
            suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{s});
            glmSubjDir = [glmLocDir{glm} filesep subj_name{s}];
            suitGLMDir = fullfile(suitDir,sprintf('glm%d_loc',glm),subj_name{s});
            dircheck(suitGLMDir);
            mask   = fullfile(suitSubjDir,['c_anatomical_pcereb.nii']);
            params = fullfile(suitSubjDir,['mc_anatomical_snc.mat']);
            
            switch type
                case 'betas'
                    images='beta_0';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                case 'contrast'
                    images='con';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
            end
            for i=1:size(source,1)
                P=fullfile(glmSubjDir,source(i).name);
                O=fullfile(suitGLMDir,[source(i).name]);
                cd(fullfile(suitGLMDir));
                suit_reslice(P,params,'outfilename',O,'mask',mask,'interp',1,'vox',[1 1 1]);
            end;
        end;
    case 'SUIT_make_mask'                                                   % STEP 4.6   :  Make cerebellar mask (glm as input)
        vararginoptions(varargin,{'sn','glm'});
        for s=sn
            suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{s});
            mask_all = fullfile(baseDir,sprintf('glm%d',glm),subj_name{s}, 'mask.nii');
            mask_cerebellum = fullfile(suitSubjDir,['c_anatomical_pcereb.nii']);
            mask_suit = fullfile(suitSubjDir, 'mask_suit.nii'); % create
            spm_imcalc_ui({mask_all,mask_cerebellum}, mask_suit, 'i1>0 & i2>0');
        end;        
    case 'SUIT_roi'                                                         % STEP 4.7   :  Generate ROI from cerebellar atlas
        vararginoptions(varargin,{'sn'});
        for s=sn
            suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{s});
            cd(suitSubjDir);
            
            SUIT_num = [1 3 5 2 4 7]; % Lobules IV, V, VI - left/right
            V = spm_vol(fullfile(fileparts(which('suit_reslice')), 'atlas', 'Cerebellum-SUIT.nii'));
            R_atlas= spm_read_vols(V);
            R=zeros(size(R_atlas));
                     
            for i= 1:6
                R(R_atlas==SUIT_num(i))=i;
            end
            % lh: 1,2,3; rh: 4,5,6 (LIV-VI)
            V.fname= fullfile(['ROI_cerebellum_orig.nii']);
            spm_write_vol(V,R);
                       
            %----create the ROI image in subj space saved in the anatomical folder
            refImage= fullfile(baseDir, 'imaging_data',subj_name{s}, ['rmask_noskull.nii']);
            defMat = fullfile(suitSubjDir,'mc_anatomical_snc.mat');
            source = fullfile(suitSubjDir, 'ROI_cerebellum_orig.nii');
            suit_reslice_inv(source, defMat,'reference',refImage,'prefix', 'subjspace_'); %ROI_cerebellum_orig is now in subject space!!!
        end;    
    
    case 'job1'
      %  sml1_imana_prep('ROI_define','sn',[29,31],'parcelType','Brodmann');
        sml1_imana_BG_new('ROI_define_BG','sn',[29,31],'regType','BG-striatum');
        sml1_imana_BG_new('ROI_define_thalamus','sn',[29,31]);
    case 'job2'
        sml1_imana_prep('GLM_estimate_RepSup','sessN',3,'sn',[4:9,11:25]);
    case 'job3'
        sml1_imana_prep('GLM_contrast_RepSup','sn',[4:9,11:31],'sessN',[1,2,4]);

    case 'subjStruct'
        sn=[4:9,11:28,30];
        sessN = 1;
        vararginoptions(varargin,{'sn','sessN'});
        T=[];
        for ss = sessN
            for s = sn
                D = dload(fullfile(behavDir,['sml1_',subj_name{s},'.dat']));
                
                L = getrow(D,D.ScanSess==ss);    % only blocks of that scan session
                if ss == 4
                    Rr = getrow(L,L.blockType==9); %blockType==9 - func imaging run without metronome
                else
                    Rr = getrow(L,L.blockType==3); %blockType==3 - funct imaging run with metornome
                end   
                for c = 1:numel(num_seq)
                    t = getrow(Rr,Rr.seqNumb==c);
                    % Do some subject info for fields in SPM_info.mat.
                    S.SN    		= s;
                    S.sessN  		= ss;    
                    S.seqNumb       = c;
                    S.seqType       = unique(t.seqType);
                    S.presses       = [t.press0(1) t.press1(1) t.press2(1) t.press3(1) t.press4(1) t.press5(1) t.press6(1) t.press7(1) t.press8(1)];
                    S.pressID       = unique(t.cueP);
                    T				= addstruct(T,S);
                end;
                
            end
        end
        keyboard;
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

