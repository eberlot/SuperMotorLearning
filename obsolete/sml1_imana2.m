function varargout=sml1_imana2(what,varargin)

% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         = '/Volumes/MotorControl/data/SuperMotorLearning/'; 
behavDir        =[baseDir '/behavioral_data'];            
imagingDir      =[baseDir '/imaging_data'];              
imagingDirRaw   =[baseDir '/imaging_data_raw'];           
dicomDir        =[baseDir '/imaging_data_dicom'];         
anatomicalDir   =[baseDir '/anatomicals'];       
fieldmapDir     =[baseDir '/fieldmaps/'];
freesurferDir   =[baseDir '/surfaceFreesurfer'];          
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           = fullfile(baseDir,'basal_ganglia');
% update glmDir when adding new glms
glmDir          ={[baseDir '/glm1'],[baseDir '/glm2'],[baseDir '/glm3']};

% ------------------------- Experiment Info -------------------------------
numDummys  = 4;        % per run
numTRs     = 480;      % per run (includes dummies)
numruns    = [7 6 5];  % per subject
run        = {'1','2','3','4','5','6','7'};
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

regname    = {'BA1','BA2','BA3','BA4','BA6','CaudateN' 'Pallidum', 'Putamen' 'Thalamus'}; % Cortical ROIs
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen' 'Thalamus'};
numregions_surf = 5;
numregions_BG   = 4;
%numregions_cerebellum - to be defined!!
numregions = numregions_surf+numregions_BG;         % sum - volume + surface
regSide=[1 1 1 1 1  2 2 2 2 2]; % 1-left, 2-right
regType=[1 2 3 4 5  1 2 3 4 5]; % Brodmann areas: 1,2,3,4,6


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

subj_name  = {'p01a','p01b','p01c'};  % same pilot subject
loc_AC     = {[-106 -163 -172],...
    [-106 -163 -172],...
    [-106 -163 -172]};

DicomName  = {'2017_03_15_TQ-001.MR.DIEDRICHSEN_LONGSEQLEARN',...
    '2017_03_16_TQ-02.MR.DIEDRICHSEN_LONGSEQLEARN',...
    '2017_03_22_TQ-03.MR.DIEDRICHSEN_LONGSEQLEARN'};

NiiRawName = {'170315110129DST131221107523418932',...
    '170316091541DST131221107523418932',...
    '170322145028DST131221107523418932'};

fscanNum   = {[16 18 20 22 24 28 30],...                                 
    [12 18 20 22 28 30],...
    [11 15 17 25 27]};

anatNum    = {[10:14]};   % same subject, 1 anatomical 

fieldNum   = {[31,32],...
    [31,32],...
    [28,29]};

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
        glm = 1;
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
                    fprintf('Rename this file to ''s%02d_anatomical_raw.nii'' in the anatomical folder.\n',s);
                case 'fieldmap'
                    fprintf('Subject %02d fieldmaps imported.\n',s)
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
            sml1_imana('PREP_realign','sn',s);
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
%         for i=1:numruns(sn)
%             runs{i}=int2str(i);
%         end;
        if sn == 1
            runs    = {'_01','_02','_03','_04','_05','_06','_07'};
        elseif sn == 2
            runs    = {'_01','_02','_03','_04','_05','_06'};
        elseif sn == 3
            runs    = {'_01','_02','_03','_04','_05'};
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
             runs    = {'_01','_02','_03','_04','_05'};
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
            J.warp.write = [0 0];
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
        delay        = 2.0;                                                  % adjusting hrf per subject based on extracted timeseries!
        announceTime = 2.0;                                                   % length of task announce time - currently not used
        % Gather appropriate GLM presets.
        switch glm
            case 1
                hrf_params = [4.5 11];
                hrf_cutoff = 128;
                cvi_type   = 'wls';
            case 2
                hrf_params = [4.5 11];
                hrf_cutoff = 128;
                cvi_type   = 'fast';
            case 3
                hrf_params = [4.5 11];
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
                run={'1','2','3','4','5'}; 
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
                    J.sess(r).cond(c).duration = dur(sn);                       % durations of task we are modeling (not length of entire trial)
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
            J.bases.hrf.params = hrf_params;
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
    
    case '4_ROI' % ------------- ROI: roi analyses. Expand for more info. ------
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
    case 'ROI_makepaint'                                                    % STEP 4.1   :  Make paint file for ROIs
        % Creates ROI boundaries on the group template/atlas (fsaverage).
        % ROIs are defined using:
        %       - probabilistic atlas
        %       - flatmap coordinates for X Y coordinate restriction
        
        for h=1:2
            groupDir=[caretDir filesep 'fsaverage_sym' filesep hemName{h} ];
            cd(groupDir);
            M=caret_load([hem{h} '.propatlas.metric']);
            C=caret_load([hem{h} '.FLAT.coord']);
            
            PROI(:,1)=M.data(:,1);
            PROI(:,2)=M.data(:,2);
            PROI(:,3)=sum(M.data(:,[3 4]),2);
            PROI(:,4)=sum(M.data(:,[7 8]),2);
            PROI(:,5)=M.data(:,9);
            [Prop,ROI]=max(PROI,[],2);
            ROI(Prop<0.2)=0;
            ROI(C.data(:,2)<-15 | C.data(:,2)>25)=0; % Hand area M1
            ROI(C.data(:,1)<-22)=0; % Take caudal PMd only
            % Save Paint file
            Paint=caret_struct('paint','data',ROI,'paintnames',regname,'column_name',{'ROI'});
            caret_save(['ROIsm.paint'],Paint);
        end;
    case 'ROI_define'                                                       % STEP 4.2   :  Define ROIs
        glm=3;
        vararginoptions(varargin,{'sn','glm'});
        for s=sn
            R=[];
            for h=1:2
                C = caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},['ROIsm.paint']));
                caretSubjDir = fullfile(caretDir,['x' subj_name{1}]);   % same subject - anatomy
                file = fullfile(glmDir{glm},subj_name{s},'mask.nii');
                for i=1:numregions_surf
                    R{i+(h-1)*numregions}.type='surf_nodes';
                    R{i+(h-1)*numregions}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                    R{i+(h-1)*numregions}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                    R{i+(h-1)*numregions}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                    R{i+(h-1)*numregions}.flat=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                    
                    R{i+(h-1)*numregions}.linedef=[10,0,1];
                    R{i+(h-1)*numregions}.image=file;
                    R{i+(h-1)*numregions}.name=[subj_name{s} '_' regname{i} '_' hem{h}];
                    R{i+(h-1)*numregions}.location=find(C.data(:,1)==i);
                end;

            end;
            R=region_calcregions(R);
            cd(regDir);
            save([subj_name{s} '_regions.mat'],'R');
            fprintf('ROIs have been defined for %s \n',subj_name{sn});
        end
         
    case 'ROI_MNI_normalization_write'%___________________________________

        vararginoptions(varargin,{'sn'});
        
        for s=sn
            defor= fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);

            % Do the anatomical
            sn_images={}; out_images={};
            sn_images=fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']);
            out_images=fullfile(anatomicalDir, subj_name{s},[subj_name{s},'_anatomical_MNI.nii']);
            spmj_normalization_write(defor, sn_images,'outimages',out_images);

        end;    
    case 'ROI_BG_MNImask'                     % Do segmentation for BG using FSL in MNI space
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
            case 2 % make the ROI images in subject space and do mni transform
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
                            
                            %----make subj specific ROI image; still in MNI space!
                            spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                        end
                        %----do deformation
                        spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end
            case 3 % make the avrg mask image
                for h=1:2
                    for i=1:numregions_BG
                        for s = 1:numel(sn)%1:numel(subj_name)
                            IN{s} = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', subj_name{sn(s)},...
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
            case 4 % make avrg mask image for all defined BG roi
                outDir = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', 'avrg');
                IN = {};
                for h=1:2
                    for i=1:numregions_BG                        
                        IN{end+1} = fullfile(outDir,...
                                    ['avrg_',regname{i+numregions_surf},'_',hem{h}, '.nii']);                        
                    end
                end
                OUT = fullfile(outDir,...
                            ['avrg_BG_all.nii']);
                spmj_imcalc_mtx(IN,OUT,'sum(X)');
        end
    case 'ROI_define_BG'                      % Direct specification of the BG ROIS
        %sl1_imana('ROI_define_BG', 1:16)
        sn      = varargin{1};                
        regtype = 'regAll_movement_1.mat';
        
        vararginoptions(varargin(2:end), {'regtype'});
        
        % -------------------------------
        % Read in BG-regions in MNI space: define possibly different
        for h=1:2
            for i=1:numregions_BG
                %----make subj specific ROI image; still in MNI space!
                %IN=fullfile(baseDir, 'basal_ganglia', [regname{i+numregions_surf},'_',hem{h}, '.nii']);
                %BGreg{h,i}=region('roi_image',IN,255);
                % defined through FSL BG
                IN=fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', 'avrg',...
                    ['avrg_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                BGreg{h,i}=region('image',IN,0.001);
                BGreg{h,i}=region_calcregions(BGreg{h,i});
            end;
        end;
        
        % -------------------------------
        % Transform each of the regions into the individual space
        cd(regDir);
        for s=sn            
            %nam_def = fullfile(anatomicalDir,subj_name{s},
            %[subj_name{s},'_anatomical_seg_sn.mat']); % this should be
            %inverse map
            nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_seg_inv_sn.mat']);
            for c=1:2
                load(fullfile('./',subj_name{s},[subj_name{s},'_',regtype]));
                
                for h=1:2
                    for i = 1:numregions_BG
                        num = numregions_surf*2 + numregions_BG*(h-1) + i; % Append BG ROI info
                        R{num} = region_deformation(BGreg{h,i},nam_def);
                        R{num}.xyz_mni = BGreg{h,i}.data;             % Use the original MNI space coordinates
                        R{num}.name = [subj_name{s},'_',regname_BG{i},'_',hem{h}];
                    end;
                end; %hem
                save(fullfile('./',subj_name{s},[subj_name{s},'_',regtype]),'R');
            end %con
            fprintf('Working on subj: %d\n',s);
        end; %subj
        varargout={R};    
        
    case 'ROI_timeseries'                                                   % STEP 4.3   :  Extract onsets and events/trials - hrf
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
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));   % This is made in case 'ROI_define'
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
    case 'ROI_plot_timeseries_reg'                                          % STEP 4.4   :  Plot estimated hrf by region
        vararginoptions(varargin,{'sn','glm'});                             
        T=load(fullfile(regDir,sprintf('hrf_%s_glm%d.mat',subj_name{sn},glm)));
        
        %T=getrow(T,isincluded(sn,T.sn) & isincluded(reg,T.regType));
        traceplot([-4:16],T.y_adj,'errorfcn','stderr','subset',T.regSide==1,'split',T.region,'leg','auto');
        hold on;
        traceplot([-4:16],T.y_hat,'subset',T.regSide==1 ,'linestyle',':','split',T.region);
        hold off;
        xlabel('TR');
        ylabel('activation');
        drawline(0);      
        keyboard;
    case 'ROI_plot_timeseries_cond'                                         % STEP 4.5   :  Plot estimated hrf by sequence type                     
        sn=1; 
        glm=2; 
        reg=3;
        regS=2; 
        vararginoptions(varargin,{'sn','glm','reg'});
        T=load(fullfile(regDir,sprintf('hrf_%s_glm%d.mat',subj_name{sn},glm)));
        
        traceplot([-4:16],T.y_adj,'errorfcn','stderr','subset',T.regSide==regS & T.regType==reg,'split',T.seqType,'leg','auto');
        hold on;
        traceplot([-4:16],T.y_hat,'subset',T.regSide==regS & T.regType==reg,'linestyle',':','split',T.seqType);
        hold off;
        xlabel('TR');
        ylabel('activation');
        drawline(0);      

    otherwise
        disp('there is no such case.')
end;    % switch(what)
end


%  % Local functions

function dircheck(dir)
% Checks existance of specified directory. Makes it if it does not exist.
% SArbuckle 01/2016
if ~exist(dir,'dir');
    %warning('%s didn''t exist, so this directory was created.\n',dir);
    mkdir(dir);
end
end