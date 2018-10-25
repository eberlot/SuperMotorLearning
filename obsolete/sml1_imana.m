function varargout=sml1_imana(what,varargin)

% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
behavDir        =[baseDir '/behavioral_data'];            
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

% update glmDir when adding new glms
glmDir          ={[baseDir '/glm1'],[baseDir '/glm2'],[baseDir '/glm3']};

% ------------------------- Experiment Info -------------------------------
numDummys  = 4;        % per run
numTRs     = 480;      % per run (includes dummies)
numruns    = [7 6 4 8 8 8];  % per subject
run        = {'1','2','3','4','5','6','7','8'};
TRlength   = 1000;      % in ms
% seqNumb - all sequences: 1-19
% seqType - types of sequences 
    % 1-trained
    % 2-novel
    % 3-trained other hand
    % 4-random sequences (behavioural training only)
    % 5-single finger mapping 

% ------------------------- ROI things ------------------------------------
hem        = {'lh','rh'};                                                   % left & right hemi folder names/prefixes
hemName    = {'LeftHem','RightHem'};
regname         = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp','CaudateN' 'Pallidum', 'Putamen' 'Thalamus','CIV','CV','CVI'};
regname_cortex  = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp'};
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen', 'Thalamus'};
regname_cerebellum = {'LobIV','LobV','LobVI'};
numregions_surf = 7;
numregions_BG   = 4;
numregions_cerebellum = 3;
numregions = numregions_surf+numregions_BG+numregions_cerebellum;        
regSide=[ones(1,14) ones(1,14)*2]; % 1-left, 2-right
regType=[1:14  1:14]; % cortical areas: 1-7, BG: 8-11, cereb: 12-14


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

subj_name  = {'p01a','p01b','p01c','p02a','p02b','p03'};  % same pilot subject
loc_AC     = {[-106 -163 -172],...
    [-106 -163 -172],...
    [-106 -163 -172],...
    [-108 -166 -166],...
    [-108 -166 -166],...
    [-107 -168 -174]};

DicomName  = {'2017_03_15_TQ-001.MR.DIEDRICHSEN_LONGSEQLEARN',...
    '2017_03_16_TQ-02.MR.DIEDRICHSEN_LONGSEQLEARN',...
    '2017_03_22_TQ-03.MR.DIEDRICHSEN_LONGSEQLEARN',...
    '2017_04_27_EVA.MR.DIEDRICHSEN_LONGSEQLEARN',...
    '2017_05_11_EVA.MR.DIEDRICHSEN_LONGSEQLEARN',...
    '2017_05_31_DN.MR.DIEDRICHSEN_LONGSEQLEARN'};

NiiRawName = {'170315110129DST131221107523418932',...
    '170316091541DST131221107523418932',...
    '170322145028DST131221107523418932',...
    '170427143409DST131221107523418932',...
    '170511130632DST131221107523418932',...
    '170531143217DST131221107523418932'};

fscanNum   = {[16 18 20 22 24 28 30],...                                 
    [12 18 20 22 28 30],...
    [11 15 17 27],...
    [16 18 20 22 24 26 28 30],...
    [16 18 20 22 24 26 28 30],...
    [20 22 24 26 28 30 32 34]};

anatNum    = {[10:14],...
    [10:14],...
    [10:14],...
    [10:14],...
    [10:14],...
    [11:15]};   % 2 anatomicals - 2 subjects; 13 - T1 

fieldNum   = {[31,32],...
    [31,32],...
    [28,29],...
    [31,32],...
    [31,32],...
    [35,36]};

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
        vararginoptions(varargin,{'sn','glm'});
        
        glmSubjDir = [glmDir{glm} filesep subj_name{sn}];
        cd(glmSubjDir);
        load SPM;
        spm_rwls_resstats(SPM)    
        
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
        vararginoptions(varargin,{'sn','series_type'});
        cwd = pwd;
        
        switch series_type
            case 'functional'
                seriesNum = fscanNum;
            case 'anatomical'
                seriesNum = anatNum;
            case 'fieldmap'
                seriesNum = fieldNum;
        end
        
        % Loop through subjects
        for s = sn;
            cd(fullfile(dicomDir,subj_name{s}));
            
            % For each series number of this subject (in 'Subject Things')
            for i=1:length(seriesNum{s})
                r     = seriesNum{s}(i);
                % Get DICOM FILE NAMES
                DIR   = dir(sprintf('%s.%4.4d.*.IMA',DicomName{s},r));
                Names = vertcat(DIR.name);
                % Convert the dicom files with these names.
                if (~isempty(Names))
                    % Load dicom headers
                    HDR=spm_dicom_headers(Names,1);
                    % Make a directory for series{r} for this subject.
                    % The nifti files will be saved here.
                    dirname = fullfile(dicomDir,subj_name{s},sprintf('series%2.2d',r));
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
                    fprintf('Subject %02d functional runs imported. Copy the unique .nii name for subj files and place into ''Subject Things''.\n',s)
                case 'anatomical'
                    fprintf('Anatomical runs have been imported for subject %d.\n',s);
                    fprintf('Please locate the T1 weighted anatomical img. Copy it to the anatomical folder.\n')
                    fprintf('Rename this file to ''%s_anatomical_raw.nii'' in the anatomical folder.\n',subj_name{s});
                case 'fieldmap'
                    fprintf('Subject %s fieldmaps imported.\n',subj_name{s});
                    fprintf('Please locate the magnitude and phase files (different series) into the fieldmap folder.\n')
                    fprintf('Rename the files into ''%s_magnitude.nii'' and ''%s_phase.nii''.\n',subj_name{s},subj_name{s});
            end
        end
        cd(cwd);

    case 'PREP_process1_func'                                               
        % need to have dicom_import done prior to this step.
        vararginoptions(varargin,{'sn'});
        
        for s = sn
            sml1_imana('PREP_make_4dNifti','sn',s);
            sml1_imana('PREP_makefieldmap','sn',s);
            sml1_imana('PREP_make_realign_unwarp','sn',s);
            sml1_imana('PREP_move_data','sn',s);
            sml1_imana('PREP_meanimage_bias_correction','sn',s);
        end
    case 'PREP_make_4dNifti'                                                % STEP 1.4       :  Converts dicoms to 4D niftis out of your raw data files
        vararginoptions(varargin,{'sn'});
        for s = sn
            % For each functional run
            for i = 1:length(fscanNum{s})                                      
                outfilename = fullfile(imagingDirRaw,subj_name{s},sprintf('%s_run_%2.2d.nii',subj_name{s},i));
                % Create a 4d nifti of all functional imgs in this run.
                % Don't include the first few dummy scans in this 4d nifti.
                for j = 1:(numTRs-numDummys)                                        
                    P{j}=fullfile(dicomDir,subj_name{s},sprintf('series%2.2d',fscanNum{s}(i)),...
                        sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{s},fscanNum{s}(i),j+numDummys,j+numDummys));
                end;
                dircheck(fullfile(imagingDirRaw,subj_name{s}))
                spm_file_merge(char(P),outfilename);
                fprintf('Run %d done\n',i);
            end
        end
    case 'PREP_makefieldmap'                                                % STEP 1.5       :  Create field map
        prefix = '';
        vararginoptions(varargin,{'sn'});

        if sn == 1
            runs    = {'_01','_02','_03','_04','_05','_06','_07'};
        elseif sn == 2
            runs    = {'_01','_02','_03','_04','_05','_06'};
        elseif sn == 3
            runs    = {'_01','_02','_03','_04'};
        else
            runs    = {'_01','_02','_03','_04','_05','_06','_07','_08'};
        end
        spmj_makefieldmap(baseDir, subj_name{sn}, runs,'prefix',prefix);
    case 'PREP_make_realign_unwarp'                                         % STEP 1.6       :  Realign + unwarp functional runs
        prefix  ='';
        vararginoptions(varargin,{'sn'});
         if sn == 1
            runs    = {'_01','_02','_03','_04','_05','_06','_07'};
        elseif sn == 2
            runs    = {'_01','_02','_03','_04','_05','_06'};
         elseif sn == 3
             runs    = {'_01','_02','_03','_04'};
        else
            runs    = {'_01','_02','_03','_04','_05','_06','_07','_08'};
         end
        spmj_realign_unwarp_sess(baseDir, subj_name{sn}, {runs}, numTRs,'prefix',prefix);
    case 'PREP_realign'                                                     % STEP 1.5/6     :  Realign functional runs (if no unwarping used)
        % SPM realigns first volume in each run to first volume of first
        % run, and then registers each image in that run to the first
        % volume of that run. Hence also why it's often better to run
        % anatomical before functional scans.

        % SPM does this with 4x4 affine transformation matrix in nifti
        % header (see function 'coords'). These matrices convert from voxel
        % space to world space(mm). If the first image has an affine
        % transformation matrix M1, and image two has one (M2), the mapping
        % from 1 to 2 is: M2/M1 (map image 1 to world space-mm - and then
        % mm to voxel space of image 2).

        % Registration determines the 6 parameters that determine the rigid
        % body transformation for each image (described above). Reslice
        % conducts these transformations; resampling each image according
        % to the transformation parameters. This is for functional only!
        vararginoptions(varargin,{'sn'});

        cd(fullfile(imagingDirRaw,subj_name{sn}));
        for s=sn;
            data={};
            for i=1:length(fscanNum{sn});
                for j=1:numTRs-numDummys;
                    data{i}{j,1}=sprintf('%s_run_%2.2d.nii,%d',subj_name{sn},i,j);
                end;
            end;
            spmj_realign(data);
            fprintf('Subj %d realigned\n',s);
        end;


    %__________________________________________________________________
    case 'PREP_plot_movementparameters'                                     % OPTIONAL       :  Investigate movement parameters
        vararginoptions(varargin,{'sn'});
        X=[];
        %r=[1:numruns(sn)];
        for r=1:numruns(sn)
            x = dlmread (fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rp_' subj_name{sn},'_run_0',num2str(r),'.txt']));
            X = [X; x];
        end
       
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
        vararginoptions(varargin,{'sn'});

        prefix='';
        dircheck(fullfile(baseDir, 'imaging_data',subj_name{sn}))
        for r=1:numruns(sn);

            source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, ['u' prefix subj_name{sn},'_run_0',run{r},'.nii']);
            dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['u' prefix subj_name{sn},'_run_0',run{r},'.nii']);

            copyfile(source,dest);
            source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, ['rp_' subj_name{sn},'_run_0',run{r},'.txt']);
            dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rp_' subj_name{sn},'_run_0',run{r},'.txt']);

            copyfile(source,dest);
        end;
        source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, ['meanu' prefix subj_name{sn},'_run_0',run{1},'.nii']);
        dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);

        copyfile(source,dest);


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
        vararginoptions(varargin,{'sn'});
        
        for s = sn
            sml1_imana('PREP_reslice_LPI','sn',s);
            sml1_imana('PREP_centre_AC','sn',s);
            sml1_imana('PREP_segmentation','sn',s);
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
            J.channel.vols = {fullfile(anatomicalDir,subj_name{sn},[subj_name{sn},'_anatomical.nii,1'])};
            J.channel.biasreg = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write = [0 0];
            J.tissue(1).tpm = {fullfile(SPMhome,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm = {fullfile(SPMhome,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm = {fullfile(SPMhome,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm = {fullfile(SPMhome,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm = {fullfile(SPMhome,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm = {fullfile(SPMhome,'tpm/TPM.nii,6')};
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
        
    case 'PREP_coreg'                                                       % STEP 1.12      :  Coregister meanepi to anatomical image
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi
        %   image
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
        % Overwrites meanepi, unless you update in step one, which saves it
        % as rmeanepi.
        % Each time you click "update" in coregtool, it saves current
        % alignment by appending the prefix 'r' to the current file
        % So if you continually update rmeanepi, you'll end up with a file
        % called r...rrrmeanepi.
      
        %__________________________________________________________________

    case 'PREP_process2'                                                                                                     
        vararginoptions(varargin,{'sn'});
        
        for s=sn
            sml1_imana('PREP_make_samealign','sn',s);
            sml1_imana('PREP_check_samealign','sn',s);
            sml1_imana('PREP_make_maskImage','sn',s);
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
            for i=1:numTRs-numDummys;
                Q{end+1}    = fullfile(imagingDir,subj_name{sn},...
                    sprintf('%s%s_run_%2.2d.nii,%d',prefix, subj_name{sn},r,i));
            end;
        end;

        % Run spmj_makesamealign_nifti to bring all functional runs into
        % same space as realigned mean epis
        spmj_makesamealign_nifti(char(P),char(Q));
    case 'PREP_check_samealign'                                             % OPTIONAL      :  Check if all functional scans are align to the anatomical
        prefix='u';
        vararginoptions(varargin,{'sn'});
        Q={};
        cd(fullfile(imagingDir,subj_name{sn}));
        for i=1:numruns(sn)
            r{i}=int2str(i);
        end;
        for r= 1:numel(r)
            for i=1:numTRs-numDummys;
                %Q{end+1} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix, subj_name{sn},'_run',r(r),'.nii,',num2str(i)])];
                Q{end+1}    = fullfile(imagingDir,subj_name{sn},...
                    sprintf('%s%s_run_%2.2d.nii,%d',prefix, subj_name{sn},r,i));
            end
        end
        P{1}= fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rbbmeanepi_' subj_name{sn} '.nii']);
        spmj_checksamealign(char(P),char(Q))       
    case 'PREP_make_maskImage'                                              % STEP 1.14     :  Make mask images (noskull and gray_only)
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
        sml1_imana('SURF_freesurfer','sn',sn);
        sml1_imana('SURF_xhemireg','sn',sn);
        sml1_imana('SURF_map_ico','sn',sn);
        sml1_imana('SURF_make_caret','sn',sn);
    case 'SURF_freesurfer'                                                  % STEP 2.1   :  Call recon-all in Freesurfer                                                   
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
    
    case '3_GLM' % ------------- GLM: SPM GLM fitting. Expand for more info. ---
        % The GLM cases fit general linear models to subject data with 
        % SPM functionality.
        %
        % All functions can be called with ('GLM_processAll','sn',[Subj#s]).      
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    case 'GLM_make'                                                         % STEP 3.1   :  Make the SPM.mat and SPM_info.mat files (prep the GLM)
        % makes the GLM file for each subject, and a corresponding 
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        
        glm = 2;    %1/2/3
        vararginoptions(varargin,{'sn','glm'});
        % Set some constants.
        prefix		 = 'u';
        T			 = [];
        dur			 = 4.5;                                                   % secs (length of task dur, not trial dur)
        delay        = [2.0 2.0 2.0 2.0 2.0 2.0];                                         % adjusting hrf per subject based on extracted timeseries!
        announceTime = 2.0;                                                   % length of task announce time - currently not used
        % Gather appropriate GLM presets.
        switch glm
            case 1
                hrf_params = [5.5 12.5]; % change to 5.5 12.5
                hrf_cutoff = 128;
                cvi_type   = 'wls';
            case 2
                hrf_params = [5.5 12.5];
                hrf_cutoff = 128;
                cvi_type   = 'fast';
            case 3
                hrf_params = [4.5 11]; % change to 5.5 12.5
                hrf_cutoff = inf;
                cvi_type   = 'fast';
        end
        % Loop through subjects and make SPM files.
        for s = sn
            D = dload(fullfile(baseDir, 'behavioral_data',['sml1_',subj_name{s},'.dat']));     
            
            if s==1   % define how many runs per subject
                run={'1','2','3','4','5','6','7'}; 
            elseif s==2
                run={'1','2','3','4','5','6'}; 
            elseif s==3
                run={'1','2','3','4'}; 
            else
                run=run;
            end                
                 
            % Do some subject structure fields.
            J.dir 			 = {fullfile(glmDir{glm}, subj_name{s})};
            J.timing.units   = 'secs';                                      % timing unit that all timing in model will be
            J.timing.RT 	 = 1.0;                                         % TR (in seconds, as per 'J.timing.units')
            J.timing.fmri_t  = 16;
            J.timing.fmri_t0 = 1;
            
            % Loop through runs. 
            for r = 1:numel(run)                                            
                R = getrow(D,D.BN==r);
                for i = 1:(numTRs-numDummys)                                % get nifti filenames, correcting for dummy scancs
                    N{i} = [fullfile(baseDir, 'imaging_data',subj_name{s}, ...
                        [prefix subj_name{s},'_run_0',run{r},'.nii,',num2str(i)])];
                end;
                J.sess(r).scans = N;                                        % number of scans in run
                % Loop through conditions.
                           
                for c = 1:numel(unique(D.seqNumb))
                    idx						   = find(R.seqNumb==c);             % find indx of all trials in run of that sequence type (seqNumb 1-14: 3x, seqNumb 15-19: 1x -> finger mapping)
                    condName = sprintf('SeqNumb-%d',R.seqNumb(idx(1)));
                    J.sess(r).cond(c).name 	   = condName;
                    % Correct start time for numDummys removed & convert to seconds 
                    J.sess(r).cond(c).onset    = [R.startTimeReal(idx)/1000 - J.timing.RT*numDummys + delay(sn)];    
                    J.sess(r).cond(c).duration = dur;                       % durations of task we are modeling (not length of entire trial)
                    J.sess(r).cond(c).tmod     = 0;
                    J.sess(r).cond(c).orth     = 0;
                    J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
					
                    % Do some subject info for fields in SPM_info.mat.
                    S.SN    		= s;
                    S.run   		= r;
                    S.seqNumb 		= R.seqNumb(idx(1));
                    S.seqType    	= R.seqType(idx(1));
                    S.hand          = R.hand(idx(1));
                    S.isMetronome   = R.isMetronome(idx(1));
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
            spm_rwls_run_fmri_spec(J);          % problem - subj 1
            % Save the aux. information file (SPM_info.mat).
            % This file contains user-friendly information about the glm
            % model, regressor types, condition names, etc.
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            
        end;
    case 'GLM_estimate'                                                     % STEP 3.2   :  Run the GLM according to model defined by SPM.mat
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        vararginoptions(varargin,{'sn','glm'});
        glm=2;
        for s = sn
            % Load files
            load(fullfile(glmDir{glm},subj_name{s},'SPM.mat'));
            SPM.swd = fullfile(glmDir{glm},subj_name{s});
            % Run the GLM.
            spm_rwls_spm(SPM);
        end;
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)
    case 'GLM_contrast'                                                     % STEP 3.3   :  Make t-contrasts for specified GLM estimates.
        % enter sn, glm #
        % 1:   Trained seq vs. rest
        % 2:   Novel seq vs. rest
        % 3:   Trained seq vs. rest
        % 4:   Finger average vs. rest
        % 5-9: Single finger mapping (1-5)
 
        vararginoptions(varargin,{'sn','glm'});
        cwd = pwd;
        % Loop through subjects.
        for s = sn
            glmSubjDir = [glmDir{glm} filesep subj_name{s}];
            cd(glmSubjDir);

            load SPM;
            SPM = rmfield(SPM,'xCon');
            T   = load('SPM_info.mat');

            %_____t contrast for trained seq vs. rest
            con                = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb<7) = 1;
            con                = con/sum(con);
            SPM.xCon(1)        = spm_FcUtil('Set',sprintf('TrainedSeq_Rest'), 'T', 'c',con',SPM.xX.xKXs);
            
            %_____t contrast for novel seq vs. rest (digit mapping)
            con                = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb>6 & T.seqNumb<11)  = 1;
            con                = con/sum(con);
            SPM.xCon(2)        = spm_FcUtil('Set',sprintf('NovelSeq_Rest'), 'T', 'c',con',SPM.xX.xKXs);

            %_____t contrast for trained seq other hand vs. rest
            con                    = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb>10 & T.seqNumb<15)  = 1;
            con                    = con/sum(con);
            SPM.xCon(3)           = spm_FcUtil('Set',sprintf('OtherHand_Rest'), 'T', 'c',con',SPM.xX.xKXs);
            
            %_____t contrast for single finger mapping vs. rest
            con                = zeros(1,size(SPM.xX.X,2));
            con(:,T.seqNumb>14)= 1;
            con                = con/sum(con);
            SPM.xCon(4)        = spm_FcUtil('Set',sprintf('DigitAny_Rest'), 'T', 'c',con',SPM.xX.xKXs);
         
            %_____t contrast for single finger mapping vs. rest
            for d = 1:5
                con                = zeros(1,size(SPM.xX.X,2));
                con(:,T.seqNumb==d+14)  = 1;
                con                = con/sum(con);
                SPM.xCon(d+4)      = spm_FcUtil('Set',sprintf('Digit%d',d), 'T', 'c',con',SPM.xX.xKXs);
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
    
    case '4_BG_SUIT' % ------------- BG and SUIT preparation. Expand for more info. ------
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
    case 'BG_FSLmask'                                                       % STEP 4.1-2  :  Do segmentation for BG using FSL (step 1&2)
        % uses FSL to do segmentation for BG and map the regions to
        % individual subject anatomy
        % need to run MATLAB through terminal 
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        
        % 1) Run step 1: sml1_imana('BG_FSLmask',1,1)
        % 2) Open the zip folder p01b_BG_all_fast_firstseg.nii.gz
        % 3) Run step 2
        
        sn = varargin{1};
        step = varargin{2};
        
        switch (step)
            case 1 % run FSL routine                
                for s= sn%1:numel(subj_name)
                    IN= fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']);
                    outDir = fullfile(baseDir, 'basal_ganglia', 'FSL');
                    if ~exist(outDir,'dir')
                        mkdir(outDir);
                    end
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
                    nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_seg_sn.mat']);
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

                    end
                end

        end
    case 'SUIT_process_all'  
        % example: sml1_imana('SUIT_process_all',1,2)'
        sn = varargin{1}; % subjNum
        glm = varargin{2}; % glmNum
        %         spm fmri - call first
        for s=sn,
            sml1_imana('SUIT_isolate_segment','sn',s);
            sml1_imana('SUIT_normalize','sn',s);
            sml1_imana('SUIT_make_mask','sn',s,'glm',glm);
            sml1_imana('SUIT_reslice',s,glm,'betas');
            sml1_imana('SUIT_reslice',s,glm,'contrast');
            sml1_imana('SUIT_make_mask','sn',s,'glm',glm);
            sml1_imana('SUIT_roi','sn',s);
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
                O=fullfile(suitDir,sprintf('glm%d',glm),subj_name{s},[source(i).name]);
                cd(fullfile(suitDir,sprintf('glm%d',glm),subj_name{s}));
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
    
    case '5_ROI'% ------------------- ROI analysis. Expand for more info. ------
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
    case 'ROI_define'                                                       % STEP 5.1   :  Define ROIs
        glm=2;
        vararginoptions(varargin,{'sn','glm'});
        for s=sn
            R=[];
            for h=1:2
                C = caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},['ROI.paint']));
                caretSubjDir = fullfile(caretDir,['x' subj_name{s}]);   % same subject - anatomy
                suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{s});

                file = fullfile(glmDir{glm},subj_name{s},'mask.nii');
                
                indx = [1 2 3 4 5 7 8]; % index of what regions to take from ROI.paint
                
                 for i=1:numregions_surf
                    R{i+(h-1)*numregions}.type='surf_nodes';
                    R{i+(h-1)*numregions}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                    R{i+(h-1)*numregions}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                    R{i+(h-1)*numregions}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                    R{i+(h-1)*numregions}.flat=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                    
                    R{i+(h-1)*numregions}.linedef=[10,0,1];
                    R{i+(h-1)*numregions}.image=file;
                    R{i+(h-1)*numregions}.name=[subj_name{s} '_' regname{i} '_' hem{h}];
                    R{i+(h-1)*numregions}.location=find(C.data(:,1)==indx(i));
                 end;

                 for j=1:numregions_BG
                    % Get basal ganglia
                    fileBG = fullfile(BGDir,'FSL',subj_name{sn},sprintf('%s_%s_%s.nii', subj_name{sn},regname_BG{j}, hem{h}));
                    R{j+i+(h-1)*numregions}.type = 'roi_image';
                    R{j+i+(h-1)*numregions}.file= fileBG;
                    R{j+i+(h-1)*numregions}.name = [subj_name{s} '_' regname_BG{j} '_' hem{h}];
                    R{j+i+(h-1)*numregions}.value = 1;
                 end    
                

                C=gifti(fullfile(caretDir,'suit_flat','Cerebellum-lobules.label.gii'));
                file_suit = fullfile(suitSubjDir,'mask_suit.nii'); 
                
                % Cerebellum
                
                for c=1:numregions_cerebellum
                    R{c+j+i+(h-1)*numregions}=region('roi_image',fullfile(suitSubjDir, ['subjspace_ROI_cerebellum_orig.nii']), c+(numregions_cerebellum*(h-1)));
                    R{c+j+i+(h-1)*numregions}.name=[subj_name{s} '_' regname_cerebellum{c} '_' hem{h}];
                end;
                
                          
            end;
                 
            R=region_calcregions(R);
            cd(regDir);
            save([subj_name{s} '_glm' num2str(glm) '_regions.mat'],'R');
            
            
            fprintf('\nROIs have been defined for %s \n',subj_name{sn});
        end         
    case 'ROI_make_nii'                                                     % OPTIONAL   :  Convert ROI def (.mat) into multiple .nii files (to check!)
              
        glm = 2;
        vararginoptions(varargin,{'sn','glm'});        
        
        for s=sn
            glmSubjDir = [glmDir{glm} filesep subj_name{s}];
            suitSubjDir = fullfile(suitDir,'anatomicals',subj_name{s});

            cd(glmSubjDir);
            % load ROI definition
            load(fullfile(regDir,sprintf('%s_glm%d_regions.mat',subj_name{s},glm)));

            % loop over rois
            for roi = 1:size(R,2)
                % mask volume
                mask = fullfile(glmSubjDir,'mask.nii');           
                % Save region file as nifti
                cd(regDir);
                if sum(roi==[12:14 27:28])==1 % for cerebellar ROIs 
                    mask = fullfile(suitSubjDir,'mask_suit.nii');
                end
                region_saveasimg(R{roi},mask);      
            end
            
        end    
    case 'ROI_timeseries'                                                   % STEP 5.2   :  Extract onsets and events/trials - hrf
        % to check the model quality of the glm
        glm=2;
        vararginoptions(varargin,{'sn','glm'});

        pre=4;          % How many TRs before the trial onset
        post=16;        % How many TRs after the trial onset
        T=[];
        for s=sn
            fprintf('Extracting the onsets and events for subject %s and glm %d\n',subj_name{s},glm);
            load(fullfile(glmDir{glm},subj_name{s},'SPM.mat'));
            
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s})); % This accounts for shifting
            load(fullfile(regDir,[subj_name{s} '_glm' num2str(glm) '_regions.mat']));   % This is made in case 'ROI_define'
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);      % Gets the time series data for the data
            
            % Create a structure with trial onset and trial type (event)
            D=spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
            % D.event - conditions (seqNumb: 1-19)
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
                D.seqType=zeros(size(D.event,1),1);
                D.seqType(D.event<7)=1;
                D.seqType(D.event>6 & D.event<11)=2;
                D.seqType(D.event>10 & D.event<15)=3;
                D.seqType(D.event>14)=5;
                T=addstruct(T,D);
            end;
        end;
        cd(regDir);
        save(sprintf('hrf_%s_glm%d.mat',subj_name{sn},glm),'-struct','T');
    case 'ROI_plot_timeseries'                                              % STEP 5.3   :  Plot estimated hrf by sequence type                     
        sn=6; 
        glm=2; 
        reg=5;
        regS=1; % 1 - LH, 2 - RH
        vararginoptions(varargin,{'sn','glm','reg','regS'});
        T=load(fullfile(regDir,sprintf('hrf_%s_glm%d.mat',subj_name{sn},glm)));
        
        traceplot([-4:16],T.y_adj,'errorfcn','stderr','subset',T.regSide==regS & T.regType==reg,'split',T.seqType,'leg','auto');
        hold on;
        traceplot([-4:16],T.y_hat,'subset',T.regSide==regS & T.regType==reg,'linestyle',':','split',T.seqType);
        hold off;
        xlabel('TR');
        ylabel('activation');
        drawline(0);      
    case 'ROI_getBetas'                                                     % STEP 5.4   :  Harvest betas from rois (raw, univ, multiv prewhit)
        glm = 2;
        sn  = 1:6;
        roi = [1:28];
        vararginoptions(varargin,{'sn','glm','roi'});
        
        T=[];
            
        % harvest
        for s=sn % for each subj
            fprintf('\nSubject: %d\n',s) % output to user
            
            % load files
            load(fullfile(glmDir{glm}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            load(fullfile(regDir,[subj_name{s} '_glm' num2str(glm) '_regions.mat']));          % load subject's region parcellation & depth structure (R)

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
        save(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)),'-struct','T'); 
        fprintf('\n')
    case 'ROI_patternconsistency'                                           % OPTIONAL   :  Calculates pattern consistencies for each subject in roi across glms.
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
        glm = 2;
        sn  = 1;
        roi = 2; % default LH primary motor cortex
        beta = 1;  % RAW or UW betas produce the best effect
        removeMean = 'yes'; % are we removing pattern means for patternconsistency?
        vararginoptions(varargin,{'sn','glm','roi','beta','removeMean'});
        
        if strcmp(removeMean,'yes')
             keepmean = 1; % we are removing the mean
        else keepmean = 0; % we are keeping the mean (yeilds higher consistencies but these are biased)
        end
  
        Rreturn=[];
        %========%
        for g=glm
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',g))); % loads in struct 'T'
            for r=roi
                Rall=[]; %prep output variable
                for s=sn
                    S = getrow(T,(T.SN==s & T.region==r));
                    runs = 1:numruns(s);
                    switch(beta)
                        case 0
                            betaW  = S.betaW{1}; 
                        case 1
                            betaW  = S.betaUW{1}; 
                        case 2
                            betaW  = S.betaRAW{1}; 
                    end
                    
                    % make vectors for pattern consistency func
                    conditionVec = kron(ones(numel(runs),1),[1:19]');
                    partition    = kron(runs',ones(19,1));
                    % calculate the pattern consistency
                    R2   = rsa_patternConsistency(betaW,partition,conditionVec,'removeMean',keepmean);
                    Rall = [Rall,R2];
                end
                Rreturn = [Rreturn;Rall];
            end
        end
        varargout = {Rreturn};
        % output arranged such that each row is an roi, each col is subj
        
        %_______________    
    case 'ROI_stats'                                                        % STEP 5.5   :  Calculate stats/distances on activity patterns
        glm = 2;
        sn  = 1:6;
        roi = [1:28];
        vararginoptions(varargin,{'sn','glm','roi'});
        
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        
        % output structures
        Ts = [];
        To = [];
        
        % do stats
        for s = sn % for each subject
            D = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
            fprintf('\nSubject: %d\n',s)
            num_run = numruns(s);
            
            for r = roi % for each region
                S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                fprintf('%d.',r)
                
                betaW  = S.betaUW{1};
                
                % % TseqNumb structure stats (all seqNumb - 19 conditions)
                % crossval second moment matrix
                %[G,Sig]     = pcm_estGCrossval(betaW(1:(19*num_run),:),D.run,D.seqNumb);
                [G,Sig]     = pcm_estGCrossval(betaW(1:(19*num_run),:),D.run(1:(19*num_run),:),D.seqNumb(1:(19*num_run),:));
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
                
                % % TseqType structure stats- structure for each seq
                % type (4 conditions)
                for t=1:4 % for each sequenceType condition
                    type_indx = [1 2 3 5];  % index for condition types
                    % distances
                    indx_subt  = min(D.seqNumb(D.seqType==type_indx(t)))-1;
                    Ss.RDM     = {rsa.distanceLDC(betaW(D.seqType==type_indx(t)),D.run(D.seqType==type_indx(t)),D.seqNumb(D.seqType==type_indx(t))-indx_subt)};   % because different number of seq
                    Ss.act     = mean(mean(betaW(D.seqType==type_indx(t),:)));
                    act = mean(betaW(D.seqType==type_indx(t),:),2);
                    num_example = length(act)/num_run;
                    a=zeros(1,num_example);
                    for i = 1:num_example
                        a(i)=mean(act(i:num_example:length(act)));
                    end
                    Ss.actInd  = {a};
                    % indexing fields
                    Ss.SN      = s;
                    Ss.region  = r;
                    Ss.seqType = type_indx(t);
                    Ss.regSide = regSide(r);
                    Ss.regType = regType(r);
                    Ts         = addstruct(Ts,Ss);
                end
            end; % each region
        end; % each subject

        % imagesc(rsa_squareRDM(To.RDM(To.region==1 & To.SN==1,:)));
        % % save
        save(fullfile(regDir,sprintf('glm%d_reg_TseqType.mat',glm)),'-struct','Ts');
        save(fullfile(regDir,sprintf('glm%d_reg_TallSeq.mat',glm)),'-struct','To');
        fprintf('\nDone.\n')    
    
    case 'PLOT_FingDist'
        glm = 2;
        vararginoptions(varargin,{'glm'});
        
        Ts = load(fullfile(regDir,sprintf('glm%d_reg_TallSeq.mat',glm))); % loads region data (Ts)
        for j = 1:2
            for i = 1:6
                AllDist = rsa_squareRDM(Ts.RDM(Ts.region==j & Ts.SN==i,:));
                Digit = triu(AllDist(15:19,15:19));
                DigitStr = Digit(Digit~=0);
                
                figure(1)
                hold on
                subplot(6,1,i)
                plot(ssqrt(DigitStr),'linewidth',2)
                
                DigitStrNorm = DigitStr./mean(DigitStr);
                
                figure(2)
                hold on
                subplot(6,1,i)
                
                plot(ssqrt(DigitStrNorm),'linewidth',2)
                
                if j==2
                    figure(1)
                    legend('S1','M1');
                    figure(2)
                    legend('S1','M1');                  
                end
            end
        end
    case 'PLOT_SeqDist'
        
        roi = [3,6]; % PMd, PMv, SMA, SPLa, SPLp
        sn = 1:6;
        seq = 'new';
        glm=2;
        vararginoptions(varargin,{'sn','roi','seq','glm'});
        
        Ts = load(fullfile(regDir,sprintf('glm%d_reg_TallSeq.mat',glm))); % loads region data (Ts)

        switch(seq)
            case 'learnt'
                seq_indx = [1:6];
            case 'new'
                seq_indx = [7:10];
            case 'otherhand'
                seq_indx = [11:14];
        end
        
        for r=roi
            for s=sn
                AllDist = rsa_squareRDM(Ts.RDM(Ts.region==r & Ts.SN==s,:));
                Seq = AllDist(seq_indx,seq_indx);
                
                figure(r);
                subplot(max(sn),1,s);
                imagesc(Seq);     
            end
        end
    case 'CALC_SeqDist'
        
        roi = [1:7]; % S1, M1, PMd, PMv, SMA, SPLa, SPLp
        sn = 1:6;
        seq = 'trained';
        plotRDM = 0;
        glm = 2;
        vararginoptions(varargin,{'sn','roi','seq','plotRDM','glm'});
        
        Ts = load(fullfile(regDir,sprintf('glm%d_reg_TallSeq.mat',glm))); % loads region data (T)

        switch(seq)
            case 'trained'
                seq_indx = [1:6,1:6];
            case 'untrained'
                seq_indx = [7:10,7:10];
            case 'otherhand'
                seq_indx = [11:14,11:14];
            case 'finger'
                seq_indx = [15:19,15:19];
        end
        
        for r=roi
            for s=sn
                AllDist = rsa_squareRDM(Ts.RDM(Ts.region==r & Ts.SN==s,:));
                
                Seq = triu(AllDist(seq_indx,seq_indx));
                SeqAll = Seq(Seq~=0);
                
                Dist(r,s)=mean(ssqrt(SeqAll));  
                
                switch(plotRDM)
                    case 1
                        figure(r)
                        title(regname{r})
                        subplot(max(sn),1,s)
                        imagesc(AllDist(seq_indx,seq_indx));
                    case 0
                end
            end
        end 
        myboxplot([],Dist');
        keyboard;
        
    case 'CALC_SeqPairDist'
        
        roi = [1:7]; % S1, M1, PMd, PMv, SMA, SPLa, SPLp
        sn = 1:5;
        seqpair = 'trained_untrained';
        plotRDM = 0;
        glm = 2;
        
        vararginoptions(varargin,{'sn','roi','seq','plotRDM','glm'});
        
        Ts = load(fullfile(regDir,sprintf('glm%d_reg_TallSeq.mat',glm))); % loads region data (T)

        switch(seqpair)
            case 'trained_untrained'
                seq_indx1 = [1:6];
                seq_indx2 = [7:10];
            case 'trained_otherhand'
                seq_indx1 = [1:6];
                seq_indx2 = [11:14];
        end
        
        for r=roi
            for s=sn
                AllDist = rsa_squareRDM(Ts.RDM(Ts.region==r & Ts.SN==s,:));
                
                SeqPair = AllDist(seq_indx1,seq_indx2);
                
                Dist(r,s)=mean(ssqrt(SeqPair(:)));  
                
                switch(plotRDM)
                    case 1
                        figure(r)
                        title(regname{r})
                        subplot(max(sn),1,s)
                        imagesc(AllDist(seq_indx1,seq_indx2));
                    case 0
                end
            end
        end 
        
        keyboard;
    case 'Class_Seq'
        
        glm = 2;
        sn = 1:5;
        sequence = 'finger';
        roi = 2;
        vararginoptions(varargin,{'glm','sn','sequence','roi'});
        
        switch(sequence)
            case 'trained'
                seq = 1;
            case 'untrained'
                seq = 2;
            case 'otherhand'
                seq = 3;
            case 'finger'
                seq = 5;
        end
        
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        
        for reg = roi
            for s = sn
                D = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                
                betas = T.betaUW{T.SN==s & T.region==reg};
                
                % do classification
                for r=1:numruns(s)   % all runs/shuffle combinations of train+test
                    
                    % calculate mean pattern out of 7 training sets - cent
                    test = betas(D.seqType==seq & D.run==r,:);
                    train = betas(D.seqType==seq & D.run~=r,:);   % 5 centroids
                    
                    cond_num = numel(unique(D.seqNumb(D.seqType==seq))); % how many different conditions
                    
                    for t = 1:cond_num         % test betas
                        for c = 1:cond_num     % defined centroids
                            train_mean = mean(train(c:cond_num:size(train,1),:));
                            dist(t,c,r) = (test(t,:)-train_mean)*(test(t,:)-train_mean)';
                        end
                        
                        chosFing(r,t) = find(dist(t,:,r) == min(dist(t,:,r)));  % choose finger - classify
                    end
                end
                
                acc=0;
                for class = 1:cond_num
                    acc = acc + sum(chosFing(:,class)==class);
                end 
                acc_group(reg,s) = acc/length(chosFing(:));
            end
        end
        keyboard;
    case 'RDM_consistency'
        roi = 1; 
        sn = 1;
        seq = 'finger';
        plotRDM = 0;
        glm = 2;
        vararginoptions(varargin,{'sn','roi','seq','plotRDM','glm'});
        
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
         for reg = roi
            for s = sn
                D = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                betas = T.betaUW{T.SN==s & T.region==reg};
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