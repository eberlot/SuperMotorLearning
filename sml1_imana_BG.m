function varargout=sml1_imana_BG(what,varargin)

% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
behavDir        =[baseDir '/behavioral_data/data'];            
imagingDir      =[baseDir '/imaging_data'];                     
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia'];
suitDir         =[baseDir '/suit'];
physioDir       =[baseDir '/physio'];
pcmDir          =[baseDir '/pcm_stats'];

% update glmDir when adding new glms
glmLocDir       ={[baseDir '/glmLoc/glmL1'],[baseDir '/glmLoc/glmL2'],[baseDir '/glmLoc/glmL3']};   % localiser glm
glmLocSessDir   ={[baseDir '/glmLocSeiess/glmLocSess1'],[baseDir '/glmLocSess/glmLocSess2'],[baseDir '/glmLocSess/glmLocSess3'],[baseDir '/glmLocSess/glmLocSess4']}; % one glm for loc run per session
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session
glmFoSExDir     ={[baseDir '/glmFoSEx/glmFoSEx1'],[baseDir '/glmFoSEx/glmFoSEx2'],[baseDir '/glmFoSEx/glmFoSEx3'],[baseDir '/glmFoSEx/glmFoSEx4']};   

% ------------------------- Experiment Info -------------------------------

% Stimuli - numbers given in SeqNumb
num_train = 1:6;
num_untrain = 7:12;   
num_seq = 1:12;
num_fing = 13:17;
num_seqtype = 2;
SeqType={'trained','untrained'};

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
numruns           = [40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 20 20 0 10];
numruns_task      = 32;
numruns_loc       = 8;

sess = [repmat(1,1,10),repmat(2,1,10),repmat(3,1,10),repmat(4,1,10)];   % all sessions

sess_sn = [4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,2,2,0,1];    % per subject

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
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen', 'Thalamus'};
numregions_BG   = 4;
regSide=[ones(1,4) ones(1,4)*2]; % 1-left, 2-right
regType=[1:4  1:4]; 


% ------------------------- Freesurfer things -----------------------------         
atlasA    = 'x';                                                            % freesurfer filename prefix
atlasname = 'fsaverage_sym';                                                % freesurfer average atlas                                      % freesurfer hemisphere folder names    

% ------------------------- Subject things --------------------------------
% The variables in this section must be updated for every new subject.

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22'};  

% subjects where BG included
sn_BG=[4,5,6,7,8,9,11,12,13,14,15,16];
subj_name_BG = {'s04','s05','s06','s07','s08','s09','s11','s12','s13','s14','s15','s16'};

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
% ------------------------------ Analysis Cases --------------------------------

switch(what)
    
    % Define BG - normalisation via FSL - FLIRT / FNIRT, individual subject
    case 'BG_FSLmask'                                                       % STEP 4.1-2  :  Do segmentation for BG using FSL (step 1&2)
        % uses FSL to do segmentation for BG and map the regions to
        % individual subject anatomy
        % need to run MATLAB through terminal
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        
        % 1) Run step 1: sml1_imana_BG('BG_FSLmask',1,1) / sml1_imana_prep('BG_FSLmask',1,1)
        % 2) Open the zip folder s01_BG_all_fast_firstseg.nii.gz
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
                            fprintf('Working on subj: %i region: %s \n', s, [regname_BG{i},'_',hem{h}])
                            
                            %----get names!
                            IN= fullfile(baseDir, 'basal_ganglia', 'FSL', subj_name{s},...
                                [subj_name{s},'_BG_all_fast_firstseg.nii']);
                            
                            OUT{i}= fullfile(baseDir, 'basal_ganglia', 'FSL', subj_name{s},...
                                [subj_name{s},'_',regname_BG{i},'_',hem{h},'.nii']);
                            
                            OUT_MNI{i}= fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', subj_name{s},...
                                [subj_name{s},'_',regname_BG{i},'_',hem{h}, '.nii']);
                            
                            %----make subj specific ROI image
                            spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                        end
                        %----do deformation
                        % spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end
            case 3 % make the avrg mask image - CAN BE REMOVED
                for h=1:2
                    for i=1:numregions_BG
                        for s = 1:numel(sn)%1:numel(subj_name)
                            IN{s} = fullfile(baseDir, 'basal_ganglia', 'FSL',subj_name{sn(s)},...
                                [subj_name{sn(s)},'_',regname_BG{i},'_',hem{h}, '.nii']);
                        end
                        outDir = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', 'avrg');
                        if ~exist(outDir, 'dir');
                            mkdir(outDir);
                        end
                        OUT = fullfile(outDir,...
                            ['avrg_',regname_BG{i},'_',hem{h}, '.nii']);
                        spmj_imcalc_mtx(IN,OUT,'mean(X)');
                    end
                end
            case 4 % STEP1 MNI: register subject anatomical to MNI anatomical
                atlasDir='standard/MNI152_T1_1mm.nii.gz';
                for s=sn
                    subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                    IN=fullfile(anatomicalDir,[subj_name{s},filesep,sprintf('%s_anatomical.nii',subj_name{s})]);
                    OUT=fullfile(subjMNIDir,sprintf('%s_MNI_test.nii.gz',subj_name{s}));
                    REF=fullfile(BGDir,atlasDir);
                    mat=fullfile(subjMNIDir,'subjNATIVEtoMNI_test.mat');
                    %calc with FSL
                    comm=sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6',IN,REF,OUT,mat);
                    %comm=sprintf('flirt -in %s -ref %s -out %s -omat %s',IN,REF,OUT,mat);
                    fprintf('%s\n',comm);
                    system(comm);
                end
            case 5 % STEP1 MNI: register individual anatomical BG structures to MNI with transformation in 4
                atlasDir='standard/MNI152_T1_1mm.nii.gz';
                for s=sn
                    for h=1:2
                        for i=1:numregions_BG
                            subjDir=fullfile(BGDir,'FSL',subj_name{s});
                            subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                            IN=fullfile(subjDir,sprintf('%s_%s_%s.nii',subj_name{s},regname_BG{i},hem{h}));
                            OUT=fullfile(subjMNIDir,sprintf('%s_MNI_%s_%s.nii.gz',subj_name{s},regname_BG{i},hem{h}));
                            REF=fullfile(BGDir,atlasDir);
                            mat=fullfile(subjMNIDir,'subj2MNI.mat');
                            %calc with FSL
                            comm=sprintf('flirt -in %s -applyxfm -init %s -out %s -paddingsize 0.05 -interp trilinear -ref %s',IN,mat,OUT,REF);
                            fprintf('%s\n',comm);
                            system(comm);
                        end
                    end
                end
            case 6 % STEP2 MNI: register MNI anatomical stuctures to sn4's MNI anatomical BG structures
                for s=sn
                    for h=1:2
                        for i=1:numregions_BG
                            subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                            IN=fullfile(subjMNIDir,sprintf('%s_MNI_%s_%s.nii',subj_name{s},regname_BG{i},hem{h}));
                            OUT=fullfile(subjMNIDir,sprintf('%s_MNI_step2_%s_%s.nii.gz',subj_name{s},regname_BG{i},hem{h}));
                            REF=fullfile(BGDir,'FSL','MNI',subj_name{4},sprintf('s04_MNI_%s_%s.nii',regname_BG{i},hem{h}));
                            mat=fullfile(subjMNIDir,'subj2MNI_step2.mat');
                            %calc with FSL
                            comm=sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6',IN,REF,OUT,mat);
                            fprintf('%s\n',comm);
                            system(comm);
                        end
                    end
                end
            case 7 % STEP2 MNI: register MNI overall anatomy to sn4's MNI anatomy with transformation in 6
                for s=sn
                    subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                    IN=fullfile(subjMNIDir,sprintf('%s_MNI.nii',subj_name{s}));
                    OUT=fullfile(subjMNIDir,sprintf('%s_MNI_step2.nii.gz',subj_name{s}));
                    REF=fullfile(BGDir,'FSL','MNI','s04','s04_MNI.nii.gz');
                    mat=fullfile(subjMNIDir,'subj2MNI_step2.mat');
                    %calc with FSL
                    comm=sprintf('flirt -in %s -applyxfm -init %s -out %s -paddingsize 0.05 -interp trilinear -ref %s',IN,mat,OUT,REF);
                    fprintf('%s\n',comm);
                    system(comm);
                end
            case 8 % CHECK
                for s = sn%1:numel(subj_name)
                    IN{s} = fullfile(baseDir, 'basal_ganglia', 'FSL','MNI',subj_name{sn(s)},...
                        [subj_name{s},'_MNI_step2.nii']);
                end
                outDir = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', 'avrg');
                dircheck(outDir);
                OUT = fullfile(outDir,...
                    ['avrg_MNI_step2.nii']);
                spmj_imcalc_mtx(IN,OUT,'mean(X)');

            case 8 % functional - DEPRECIATED
                funcNames={'TrainSeq','UntrainSeq'};
                for ss=1%:4 %for all sessions                   
                    for s=sn
                        subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                        %REF=fullfile(BGDir,atlasDir);
                        REF=fullfile(anatomicalDir,[subj_name{s},filesep,sprintf('%s_anatomical.nii',subj_name{s})]);
                        for f=1:length(funcNames)
                            IN=fullfile(glmSessDir{ss},[subj_name{s},filesep,sprintf('psc_sess%d_%s.nii',ss,funcNames{f})]); % functional file
                            MASK=fullfile(regDir,sprintf('mask_%s.nii',subj_name{s}));
                            OUT=fullfile(subjMNIDir,sprintf('partial_%s_2anat.nii',funcNames{f})); % new naming of the reslices func file

                            mat=fullfile(subjMNIDir,sprintf('partial_%s_2anat.mat',funcNames{f}));
                            %calc with FSL
                             comm=sprintf('flirt -in %s -ref %s -dof 6 -nosearch -inweight %s -out %s -omat %s',IN,REF,MASK,OUT,mat);
                            %comm=sprintf('flirt -in %s -ref %s -dof 6 -nosearch -init /tmp/nudge_kFKUB8.xfm -inweight %s -out %s -omat %s',IN,REF,MASK,OUT,mat);
                            fprintf('%s\n',comm);
                            system(comm);
                        end
                    end
                end
            case 9 % DEPRECIATED - register functional images to MNI-aligned anatomical (whole-brain)
                funcNames={'TrainSeq','UntrainSeq'};
                atlasDir='standard/MNI152_T1_2mm.nii.gz';
                for ss=1%:4 %for all sessions
                    
                    for s=sn
                        subjBGDir=fullfile(BGDir,'FSL',subj_name{s});
                        subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                        REF=fullfile(BGDir,atlasDir);
                        %REF=fullfile(anatomicalDir,[subj_name{s},filesep,sprintf('%s_anatomical.nii',subj_name{s})]);
                        % first remove the skull - make anatomical BETted
                        %REFbet=fullfile(anatomicalDir,[subj_name{s},filesep,sprintf('%s_anatomical_BET.nii',subj_name{s})]);
                      % comm=sprintf('bet %s %s',REF,REFbet);
                       % fprintf('%s\n',comm);
                        %system(comm);
                        for f=1:length(funcNames)
                            IN=fullfile(glmSessDir{ss},[subj_name{s},filesep,sprintf('psc_sess%d_%s.nii',ss,funcNames{f})]); % functional file
                            %OUT=fullfile(subjMNIDir,sprintf('psc_sess%d_%s_MNI.nii',ss,funcNames{f})); % new naming of the reslices func file
                            OUT=fullfile(subjMNIDir,sprintf('test2_%s_MNI.nii',funcNames{f})); % new naming of the reslices func file
                            %REF=fullfile(subjMNIDir,sprintf('%s_MNI.nii.gz',subj_name{s}));
                            %mat=fullfile(subjMNIDir,'subj2MNI.mat');
                            %mat=fullfile(subjMNIDir,sprintf('func_%s_sess%d_toNATIVE.mat',funcNames{f},ss));
                            mat=fullfile(subjMNIDir,sprintf('test_%s_MNI.mat',funcNames{f}));
                            %calc with FSL
                            %comm=sprintf('epi_reg --epi=%s --t1=%s --t1brain=%s --out=%s',IN,REF,REFbet,OUT);
                            comm=sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6 -cost bbr',IN,REF,OUT,mat);
                         %   comm=sprintf('flirt -in %s -applyxfm -init %s -out %s -paddingsize 0.05 -interp trilinear -ref %s',IN,mat,OUT,REF);
                            fprintf('%s\n',comm);
                            system(comm);
                        end
                    end
                end
            case 10 % FLIRT + FNIRT
                atlasDir='standard/MNI152_T1_2mm.nii.gz';
                for s=sn
                    subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                    dircheck(subjMNIDir);
                    IN=fullfile(anatomicalDir,[subj_name{s},filesep,sprintf('%s_anatomical.nii',subj_name{s})]);
                    OUT_flirt=fullfile(subjMNIDir,sprintf('%s_MNI_FLIRT.nii.gz',subj_name{s}));
                    OUT_fnirt=fullfile(subjMNIDir,sprintf('%s_MNI_FNIRT.nii.gz',subj_name{s}));
                    REF=fullfile(BGDir,atlasDir);
                    matFLIRT=fullfile(subjMNIDir,'subjtoMNI_FLIRT.mat');
                    transFNIRT=fullfile(subjMNIDir,'subjtoMNI_FNIRT_trans.nii.gz');
                   % matFNIRT=fullfile(subjMNIDir,'subjtoMNI_FNIRT.mat.nii.gz');
                    %calc with FSL
                    comm1=sprintf('flirt -in %s -ref %s -out %s -omat %s',IN,REF,OUT_flirt,matFLIRT);
                    comm2=sprintf('fnirt --ref=%s --in=%s --aff=%s --cout=%s',REF,IN,matFLIRT,transFNIRT);
                    %comm2=sprintf('fnirt --in=%s --aff=%s --cout%s --config=%s',IN,matFLIRT,matFNIRT,REF);
                    %comm3=sprintf('applywarp --ref=%s --in=%s --out=%s --coef=%s --premat=%s',REF,IN,OUT_final,OUT_fnirt,matFLIRT);
                    comm3=sprintf('applywarp --ref=%s --in=%s --warp=%s --out=%s',REF,IN,transFNIRT,OUT_fnirt);
                    fprintf('%s\n',comm1);
                    system(comm1);
                    fprintf('%s\n',comm2);
                    system(comm2);
                    fprintf('%s\n',comm3);
                    system(comm3);
                end
        end      
    case 'COREG_anat_MNI' % coregister grey matter to subject aligned MNI anatomical (step2)
         % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi image
        
        % Note: this is done for aligning functional to subject MNI space!!
        
        vararginoptions(varargin,{'sn'});
        
        % anatomical file here is subject aligned to MNI - step2
        subjMNIDir=fullfile(BGDir,'FSL','MNI',subj_name{sn});
        cd(fullfile(subjMNIDir));
        coregtool;
        keyboard();
        
        % (2) Automatically co-register functional and anatomical images
        %sn=varargin{1};
        
        %J.ref = {fullfile(subjMNIDir,[subj_name{sn}, '_MNI','.nii'])};
        J.ref = {fullfile(subjMNIDir,[subj_name{sn}, '_MNI_step2','.nii'])};
        J.source = {fullfile(anatomicalDir,subj_name{sn},['c1' subj_name{sn} '_anatomical.nii'])}; 
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

%___________________________________________________________   
    case 'COREG_func_MNI' % DEPRECIATED coregister functional mean epi to grey matter aligned to MNI-subject anatomical
         % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi image
        
        % Note: this is done for aligning functional to subject MNI space!!
        
        vararginoptions(varargin,{'sn'});
        
        % first copy the mean epi image
        mean_epi=fullfile(imagingDir,subj_name{sn},['bbmeanepi_' subj_name{sn} '.nii']);
        mean_epi_vol=spm_vol(mean_epi);
        V=spm_read_vols(mean_epi_vol);
        mean_epi_BG=mean_epi_vol;
        mean_epi_BG.fname=fullfile(imagingDir,subj_name{sn},['bbmeanepi_MNI_' subj_name{sn} '.nii']);
        spm_write_vol(mean_epi_BG,V);
        
        % anatomical file here is subject aligned to MNI - step2
        subjAnatDIR=fullfile(anatomicalDir,subj_name{sn});
        cd(fullfile(subjAnatDIR));
        coregtool;
        keyboard();
        
        % (2) Automatically co-register functional and anatomical images
        %sn=varargin{1};
        
        %J.ref = {fullfile(subjMNIDir,[subj_name{sn}, '_MNI','.nii'])};
        J.ref = {fullfile(anatomicalDir,subj_name{sn},['rc1' subj_name{sn} '_anatomical.nii'])}; 
        J.source = {fullfile(imagingDir,subj_name{sn},['rbbmeanepi_MNI_' subj_name{sn} '.nii'])}; 
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
%___________________________________________________________
    case 'COREG_make_samealign' % DEPRECIATED reslice other functional images to aligned mean epi (rmeanepi_MNI)

        sn=[2:16];
        sessN=[1:4];
        seqType={'Train','Untrain'};
        vararginoptions(varargin,{'sn','sessN'});  

        Q={};
        for s=sn
            cd(fullfile(imagingDir,subj_name{s}));
            indx=1; % of images to align
            % Select image for reference
            P{1} = fullfile(imagingDir,subj_name{s},sprintf('rbbmeanepi_MNI_%s.nii',subj_name{s}));
            for ss=sessN
                for st=1:length(seqType)
                    % Select images to be realigned - original, save them
                    % separately
                    image    = fullfile(glmSessDir{ss},subj_name{s},...
                        sprintf('psc_sess%d_%sSeq.nii',ss,seqType{st}));
                    imageV=spm_vol(image);
                    V=spm_read_vols(imageV);
                    image_BG=imageV;
                    image_BG.fname=fullfile(glmSessDir{ss},subj_name{s},...
                        sprintf('psc_sess%d_%sSeq_MNI_BG.nii',ss,seqType{st}));
                    spm_write_vol(image_BG,V);
                    % Now select real images - the newly created one
                    Q{indx}     = image_BG.fname;
                    indx=indx+1;
                end
            end; % session
            spmj_makesamealign(char(P),char(Q));
            % Run spmj_makesamealign_nifti to bring all functional runs into
            % same space as realigned mean epis to MNI BG
        end; % subject
    
    case 'FNIRT_SEGMENT_BG'
        % uses FSL to do segmentation for BG and map the regions to
        % individual subject anatomy
        % need to run MATLAB through terminal
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        
        % 1) Run step 1: sml1_imana_BG('BG_FSLmask',1,1) / sml1_imana_prep('BG_FSLmask',1,1)
        % 2) Open the zip folder s01_BG_all_fast_firstseg.nii.gz
        % 3) Run step 2
        
        sn = varargin{1};
        step = varargin{2};
        
        switch (step)
            case 1 % run FSL routine
                for s= sn%1:numel(subj_name)
                    BGsubjDir=fullfile(BGDir,'FSL','MNI',subj_name{s});
                    IN= fullfile(BGsubjDir, [subj_name{s}, '_MNI_FNIRT.nii']);
                    outDir = fullfile(BGDir, 'FNIRT');
                    dircheck(outDir);
                    OUT= fullfile(outDir, subj_name{s}, [subj_name{s}, '_MNI_FNIRT_BG.nii']);
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
                  %  nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_to_std_sub.mat']);
                    mniSPMDir = fullfile(baseDir, 'basal_ganglia', 'FNIRT', subj_name{s});
                    
                    for h=1:2
                        for i=1:numregions_BG
                            fprintf('Working on subj: %i region: %s \n', s, [regname_BG{i},'_',hem{h}])
                            
                            %----get names!
                            IN= fullfile(mniSPMDir, [subj_name{s},'_MNI_FNIRT_BG_all_fast_firstseg.nii']);
                            
                            OUT{i}= fullfile(mniSPMDir,...
                                [subj_name{s},'_',regname_BG{i},'_',hem{h},'_MNI_FNIRT.nii']);                   
                            
                            %----make subj specific ROI image
                            spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                        end
                        %----do deformation
                        % spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end
        end
        
    case 'MNI_anat'     % SPM version
        % run SPM jobman to align anatomical to MNI space 
        vararginoptions(varargin,{'sn'});
        
        spm('defaults','fmri');
        spm_jobman('initcfg');
        for s=sn
            subjAnat=fullfile(anatomicalDir,subj_name{s},sprintf('%s_anatomical.nii,1',subj_name{s}));
            cd(fullfile(anatomicalDir,subj_name{s}));
            J.subj.vol = {subjAnat};
            J.subj.resample = {subjAnat};
            J.eoptions.biasreg = 0.0001;
            J.eoptions.biasfwhm = 60;
            J.eoptions.tpm = {'/Users/eberlot/Documents/MATLAB/spm12/tpm/TPM.nii'};
            J.eoptions.affreg = 'mni';
            J.eoptions.reg = [0 0.001 0.5 0.05 0.2];
            J.eoptions.fwhm = 0;
            J.eoptions.samp = 3;
            J.woptions.bb = [-78 -112 -70; 78 76 85];
            J.woptions.vox = [1 1 1]; % 1mm isotropic atlas
            J.woptions.interp = 4;
            
            matlabbatch{1}.spm.spatial.normalise.estwrite=J;
            spm_jobman('run',matlabbatch);
            fprintf('\n MNI alignment for %s done\n',subj_name{s});    
        end
        % creates MNI-aligned anatomical ws01_anatomical.nii and
        % transformation matrix y_s01_anatomical.nii -> used for coreg
    case 'MNI_coreg_func' % coregister different functional images (mask, contrasts, raw imaging runs) to MNI subject template
        
        sessN=[1:4];
        imageType='Tmaps';
        seqType={'Train','Untrain'};
        vararginoptions(varargin,{'sn','imageType','sessN'});
        % align functional images to MNI-aligned anatomical
        
        for s=sn
            subjMNItrans=fullfile(anatomicalDir,subj_name{s},sprintf('y_%s_anatomical.nii',subj_name{s}));
            switch (imageType)
                case 'meanEPI'
                    % first create a clean copy of the mean epi image
                    mean_epi=fullfile(imagingDir,subj_name{s},['rbbmeanepi_' subj_name{s} '.nii']);
                    mean_epi_vol=spm_vol(mean_epi);
                    V=spm_read_vols(mean_epi_vol);
                    mean_epi_BG=mean_epi_vol;
                    mean_epi_BG.fname=fullfile(imagingDir,subj_name{s},['rbbmeanepi_MNI_SPM_' subj_name{s} '.nii']);
                    spm_write_vol(mean_epi_BG,V);
                    
                    Q = sprintf('%s,1',mean_epi_BG.fname);
                case 'mask'
                    % first create a clean copy of the mask
                    maskFile=fullfile(regDir,['mask_' subj_name{s} '.nii']);
                    maskFile_vol=spm_vol(maskFile);
                    V=spm_read_vols(maskFile_vol);
                    maskFile_BG=maskFile_vol;
                    maskFile_BG.fname=fullfile(regDir,['mask_' subj_name{s} '_MNI_SPM.nii']);
                    spm_write_vol(maskFile_BG,V);
                    
                    Q = sprintf('%s,1',maskFile_BG.fname);
                case 'psc'
                    indx=1;
                    for ss=sessN
                        for st=1:length(seqType)
                            % Select images to be realigned - save them
                            % with a new name                          
                            image    = fullfile(glmSessDir{ss},subj_name{s},...
                                sprintf('psc_sess%d_%sSeq.nii',ss,seqType{st}));
                            imageV=spm_vol(image);
                            V=spm_read_vols(imageV);
                            image_BG=imageV;
                            image_BG.fname=fullfile(glmSessDir{ss},subj_name{s},...
                                sprintf('psc_sess%d_%sSeq_MNI_BG.nii',ss,seqType{st}));
                            spm_write_vol(image_BG,V);
                            
                            % Now select real the newly created one for realignment...
                            Q{indx}   = sprintf('%s,1',image_BG.fname);
                            indx=indx+1;
                        end; % seqType
                    end; % session
                case 'Tmaps'
                    indx=1;
                    for ss=sessN
                        for st=1:length(seqType)
                            % Select images to be realigned - save them
                            % with a new name
                            image    = fullfile(glmSessDir{ss},subj_name{s},...
                                sprintf('spmT_%sSeq.nii',seqType{st}));
                            imageV=spm_vol(image);
                            V=spm_read_vols(imageV);
                            image_BG=imageV;
                            image_BG.fname=fullfile(glmSessDir{ss},subj_name{s},...
                                sprintf('spmT_%sSeq_MNI_BG.nii',seqType{st}));
                            spm_write_vol(image_BG,V);
                            
                            % saved with w before!!!
                            % Now select real the newly created one for realignment...
                            Q{indx}   = sprintf('%s,1',image_BG.fname);
                            indx=indx+1;
                        end; % seqType
                    end; % session
                case 'imagingRuns'
                    indx=1;
                    indxRun=run_task';
                    images=sprintf('u%s_run_',subj_name{s});
                    source=dir(fullfile(imagingDir,subj_name{s},sprintf('*%s*',images))); % images to be resliced
                    nameFiles={source.name};
                    for imageNum=1:size(nameFiles,2)
                        % Select images to be realigned - save them
                        % with a new name
                        cd(fullfile(imagingDir,subj_name{s}));
                        image    = fullfile(imagingDir,subj_name{s},...
                            nameFiles{imageNum});
                        imageV=spm_vol(image);
                        V=spm_read_vols(imageV);
                        image_BG=imageV;
                        dirBG=fullfile(BGDir,'SPM','MNI_imaging_data',subj_name{s});
                        dircheck(dirBG);
                        
                        name=fullfile(dirBG,...
                            sprintf('MNI_%s',nameFiles{imageNum}));
                        for i=1:436
                            image_BG(i).fname=name;
                        end
                        spm_write_vol(image_BG,V);
                        
                        % Now select real the newly created one for realignment...
                        Q{indx}   = sprintf('%s,1',image_BG.fname);
                        indx=indx+1;
                    end; % number of images
            end
            
            % definitions for realignment
            
            J.subj.def = {subjMNItrans};
            J.woptions.bb = [-78 -112 -70; 78 76 85];
            J.woptions.vox = [2 2 2];
            J.woptions.interp = 4;
            % run realignment
            if strcmp(imageType,'psc') || strcmp(imageType,'Tmaps') || strcmp(imageType,'imagingRuns')
                for i=1:size(Q,2)
                    J.subj.resample = {Q{i}}; % image to realign
                    matlabbatch{1}.spm.spatial.normalise.write=J;
                    spm_jobman('run',matlabbatch);
                end
            else
            J.subj.resample = {Q}; 
            matlabbatch{1}.spm.spatial.normalise.write=J;
            spm_jobman('run',matlabbatch);
            end
            fprintf('\n Coregistration %s image for %s done\n',imageType,subj_name{s});    
        end
    case 'SEGMENT_BG_MNI_subj' % segment BG structures in SPM MNI space per subject
        % uses FSL to do segmentation for BG and map the regions to
        % individual subject anatomy
        % need to run MATLAB through terminal
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        
        % 1) Run step 1: sml1_imana_BG('BG_FSLmask',1,1) / sml1_imana_prep('BG_FSLmask',1,1)
        % 2) Open the zip folder s01_BG_all_fast_firstseg.nii.gz
        % 3) Run step 2
        
        sn = varargin{1};
        step = varargin{2};
        
        switch (step)
            case 1 % run FSL routine
                for s= sn%1:numel(subj_name)
                    IN= fullfile(anatomicalDir, subj_name{s}, ['w', subj_name{s}, '_anatomical.nii']);
                    outDir = fullfile(BGDir,'SPM',subj_name{s});
                    dircheck(outDir);
                    OUT= fullfile(outDir,[subj_name{s},'_MNI_BG.nii']);
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
                  %  nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_to_std_sub.mat']);
                    mniSPMDir = fullfile(baseDir, 'basal_ganglia', 'SPM', subj_name{s});
                    
                    for h=1:2
                        for i=1:numregions_BG
                            fprintf('Working on subj: %i region: %s \n', s, [regname_BG{i},'_',hem{h}])
                            
                            %----get names!
                            IN= fullfile(mniSPMDir, [subj_name{s},'_MNI_BG_all_fast_firstseg.nii']);
                            
                            OUT{i}= fullfile(mniSPMDir,...
                                [subj_name{s},'_',regname_BG{i},'_',hem{h},'_MNI_SPM.nii']);                   
                            
                            %----make subj specific ROI image
                            spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                        end
                        %----do deformation
                        % spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end
        end
    case 'SEGMENT_BG_MNI_avrg'
        step = varargin{1};
        
        switch (step)
            case 1 % run FSL routine
                BGDirGroup=fullfile(BGDir,'spm','avrg');
                IN= fullfile(BGDirGroup,'MNI152_t1_2mm.nii');
                OUT= fullfile(BGDirGroup,'MNI_avrg.nii');
                %calc with FSL
                comm=sprintf('run_first_all -i %s -o %s', IN, OUT);
                fprintf('%s\n',comm);
                system(comm);
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
                
                %----deform info for basal ganglia ROI to individual space
                %  nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_to_std_sub.mat']);
               BGDirGroup=fullfile(BGDir,'spm','avrg');
                
                for h=1:2
                    for i=1:numregions_BG
                        fprintf('Working on average MNI region: %s \n',[regname_BG{i},'_',hem{h}])
                        
                        %----get names!
                        IN= fullfile(BGDirGroup, 'MNI_avrg_all_fast_firstseg.nii');
                        
                        OUT{i}= fullfile(BGDirGroup,...
                           ['SPM_MNI_',regname_BG{i},'_',hem{h},'.nii']);
                        
                        %----make subj specific ROI image
                        spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                    end
                end
        end
    case 'CALC_MNI_subjAvrg'
        sn=[3:20];
        template='SPM'; % SPM or FNIRT
        vararginoptions(varargin,{'sn','type','template'});
        
        outDir = fullfile(baseDir, 'basal_ganglia', template, 'avrg');
        dircheck(outDir);
        for s = 1:numel(sn)
            IN{s} = fullfile(anatomicalDir,subj_name{sn(s)},...
                ['w',subj_name{sn(s)},'_anatomical.nii']);
        end
        OUT = fullfile(outDir,'subjAvrg_MNI.nii');
        spmj_imcalc_mtx(IN,OUT,'mean(X)');
    case 'CALC_MNI_BG_reg_subjAvrg'
        sn=[3:20];
        template='SPM'; % SPM or FNIRT
        vararginoptions(varargin,{'sn','type','template'});
        
        outDir = fullfile(baseDir, 'basal_ganglia', template, 'avrg');
        dircheck(outDir);
        for h=1:2
            for i=1:numregions_BG
                for s = 1:numel(sn) 
                    IN{s} = fullfile(BGDir,'regionMasks',sprintf('%s_%s_%s_SPM.nii',subj_name{sn(s)},regname_BG{i},hem{h}));
                   % IN{s} = fullfile(BGDir, template,subj_name{sn(s)},...
                    %    [subj_name{sn(s)},'_',regname_BG{i},'_',hem{h}, '_MNI_',template,'.nii']);   
                end
                OUT = fullfile(outDir,...
                    ['subjAvrg_anat_MNI_',regname_BG{i},'_',hem{h},'.nii']);
                spmj_imcalc_mtx(IN,OUT,'mean(X)');
            end
        end
    case 'CALC_MNI_BG_reg_groupMask'
        sn=[3:20];
        template='SPM'; % SPM or FNIRT
        vararginoptions(varargin,{'sn','template'});

        for h=1:2
            for i=1:numregions_BG
                for s = 1:numel(sn)
                    % accummulate regions for all subjects
                    subjFile=fullfile(BGDir,'regionMasks',[subj_name{sn(s)},'_',regname_BG{i},'_',hem{h},'_',template,'.nii']);
                  %  subjFile = fullfile(BGDir,template,subj_name{sn(s)},...
                   %     [subj_name{sn(s)},'_',regname_BG{i},'_',hem{h}, '_MNI_',template,'.nii']); 
                    Vol1=spm_vol(subjFile);
                    Vol2=spm_read_vols(Vol1);
                    numVox(s)=length(find(Vol2));
                end
                % calculate on average how many voxels per subject
                avrgVoxNum=round(mean(numVox));
                % load in group average 
                groupAvrg = fullfile(baseDir,'basal_ganglia',template,'avrg',sprintf('subjAvrg_anat_MNI_%s_%s.nii',regname_BG{i},hem{h}));
                VolGroup=spm_vol(groupAvrg);
                VolGroup2=spm_read_vols(VolGroup);
                voxIndx=find(VolGroup2);
                % sort voxels by numbers (how many subjects in common)
                [sortVox indx]=sort(VolGroup2(voxIndx),'descend');
                % choose the same number of voxels with most overlap
                newVox=voxIndx(indx(1:avrgVoxNum));
                % calculate mask - preserving the volume
                maskVol=zeros(size(VolGroup2));
                maskVol(newVox)=1;
                maskNii=VolGroup;
                dircheck(fullfile(BGDir,'regionMasks','avrg'));
                dest= fullfile(BGDir,'regionMasks','avrg',sprintf('subjAvrg_anat_MNI_%s_%s_MASK.nii',regname_BG{i},hem{h}));
                maskNii.fname=dest;
                maskNii=rmfield(maskNii,'descrip');
                maskNii.descrip='avrg_mask';
                spm_write_vol(maskNii,maskVol);
            end
        end
    case 'CALC_MNI_BG_whole_groupMask'
        % makes a mask for all of the BG regions (both hemisphere)
        for h=1:2
            for i=1:3
                mask=fullfile(BGDir,'regionMasks','avrg',sprintf('subjAvrg_anat_MNI_%s_%s_MASK.nii',regname_BG{i},hem{h}));
                Reg_mask=spm_vol(mask);
                Reg_mask2=spm_read_vols(Reg_mask);
                if h==1 && i==1
                    BG_mask=Reg_mask2;
                else
                    BG_mask=BG_mask+Reg_mask2;
                end
            end
        end
        % write the mask
        BG_writeMask=Reg_mask;
        BG_writeMask.fname=fullfile(BGDir,'regionMasks','avrg','subjAvrg_anat_MNI_BG_MASK.nii');
        BG_writeMask.descip='avrg_wholeBG_mask';
        spm_write_vol(BG_writeMask,BG_mask);
        
    case 'CALC_MNI_BG_subjAvrg_func'
        sn=[3:18];
        template='SPM'; % SPM
        sessN=[1:4];
        seqType={'Train','Untrain'};
        imageType='spmT';
        vararginoptions(varargin,{'sn','template','imageType','sessN'});
        
        outDir = fullfile(baseDir, 'basal_ganglia', template, 'avrg');
        dircheck(outDir);        
        for ss=sessN     
            for st=1:length(seqType)
                for  s=1:numel(sn)
                    % Select images to be realigned - save them
                    % with a new name
                    switch imageType
                        case 'spmT'
                            IN{s}    = fullfile(glmSessDir{ss},subj_name{sn(s)},...
                                sprintf('w%s_%sSeq_MNI_BG.nii',imageType,seqType{st}));
                        case 'psc'
                            IN{s}    = fullfile(glmSessDir{ss},subj_name{sn(s)},...
                                sprintf('w%s_sess%d_%sSeq_MNI_BG.nii',imageType,ss,seqType{st}));
                    end
                end; % subject
                OUT = fullfile(outDir,sprintf('subjAvrg_%s_%sSeq_sess%d.nii',imageType,seqType{st},ss));
                spmj_imcalc_mtx(IN,OUT,'nanmean(X)');
            end; % seqType
        end % session
    case 'MNI_BG_mask_func'
        sessN=[1:4];
        seqType={'Train','Untrain'};
        imageType='psc'; % psc or spmT
        template='SPM';
        vararginoptions(varargin,{'sn','template','imageType','sessN'});
        
        for ss=sessN
            for st=1:length(seqType)
                % Select images to be realigned - save them
                % with a new name
                IN  = fullfile(BGDir,template,'avrg',...
                    sprintf('subjAvrg_%s_%sSeq_sess%d.nii',imageType,seqType{st},ss));
                image=spm_vol(IN); image2=spm_read_vols(image);
                
                mask = fullfile(BGDir,'regionMasks','avrg','subjAvrg_anat_MNI_BG_MASK.nii');
                maskRead=spm_vol(mask); maskRead2=spm_read_vols(maskRead);
               
                newImage=image2.*maskRead2;
                BG_write=image;
                BG_write.fname=fullfile(BGDir,template,'avrg', sprintf('subjAvrg_%s_%sSeq_sess%d_MASKED.nii',imageType,seqType{st},ss));
                BG_write.descip=('masked_funcImage');
                spm_write_vol(BG_write,newImage);
                
            end; % seqType
        end % sess
        
        
    case 'ROI_define_BG_subj' 
        sn=1;
        regType='MNI_step2'; % native, or SPM (MNI, MNI_step2, FNIRT)
        vararginoptions(varargin,{'sn','regType'});
        for s=sn
            R=[];
            for h=1:2
                for j=1:numregions_BG
                    % Get basal ganglia
                    if strcmp(regType,'native') 
                        fileBG = fullfile(BGDir,'FSL',subj_name{s},sprintf('%s_%s_%s.nii', subj_name{s}, regname_BG{j}, hem{h}));
                    elseif strcmp(regType,'SPM')
                        fileBG = fullfile(BGDir,'SPM',subj_name{s},sprintf('%s_%s_%s_MNI_SPM.nii', subj_name{s}, regname_BG{j}, hem{h}));
                    elseif strcmp(regType,'FNIRT')
                        fileBG = fullfile(BGDir,'FNIRT',subj_name{s},sprintf('%s_%s_%s_MNI_FNIRT.nii', subj_name{s}, regname_BG{j}, hem{h}));
                    else
                        fileBG = fullfile(BGDir,'FSL','MNI',subj_name{s},sprintf('%s_%s_%s_%s.nii', subj_name{s}, regType,regname_BG{j}, hem{h}));
                    end
                    R{j+(h-1)*numregions_BG}.type = 'roi_image';
                    R{j+(h-1)*numregions_BG}.file= fileBG;
                    R{j+(h-1)*numregions_BG}.name = [subj_name{s} '_' regname_BG{j} '_' hem{h}];
                    R{j+(h-1)*numregions_BG}.value = 1;
                end
            end
            R=region_calcregions(R);
            dircheck(regDir);
            cd(regDir);
            save([subj_name{s} sprintf('_BG_%s_regions.mat',regType)],'R');

            fprintf('\nROIs in basal ganglia have been defined for %s \n',subj_name{s});
        end
    case 'ROI_define_BG_group'
        regType='subjAvrg'; % SPM or subjAvrg
        template='SPM';
        vararginoptions(varargin,{'sn','regType','template'});
        R=[];
        for h=1:2
            for j=1:numregions_BG
                % Get basal ganglia
                fileBG = fullfile(BGDir,template,'avrg',sprintf('%s_anat_MNI_%s_%s_MASK.nii', regType, regname_BG{j}, hem{h}));
                
                R{j+(h-1)*numregions_BG}.type = 'roi_image';
                R{j+(h-1)*numregions_BG}.file= fileBG;
                R{j+(h-1)*numregions_BG}.name = [regType '_' regname_BG{j} '_' hem{h}];
                R{j+(h-1)*numregions_BG}.value = 1;
            end
        end
        R=region_calcregions(R);
        dircheck(regDir);
        cd(regDir);
        save(sprintf('%s_MNI_BG_regions.mat',regType),'R');
        fprintf('\n Average ROIs - template %s in basal ganglia have been defined\n',regType);
    
    case 'ROI_make_nii'                                                     % MAKES MASK!!! :  Convert ROI def (.mat) into multiple .nii files (to check!)
        sn=3;
        regType='native'; % FNIRT, SPM, native
        imageType='func'; % func or anat
        vararginoptions(varargin,{'sn','regType'});
        
        for s=sn
            % load ROI definition
            load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
            % loop over rois
            for roi = 1:size(R,2)
                % mask volume - anatomical aligned to MNI
                switch(regType)
                    case 'SPM'
                        if strcmp(imageType,'anat')
                            maskAnat = fullfile(anatomicalDir,subj_name{s},sprintf('w%s_anatomical.nii',subj_name{s})); %anatomical aligned to MNI
                        elseif strcmp(imageType,'func')
                            maskAnat = fullfile(regDir,sprintf('wmask_%s_MNI_SPM.nii',subj_name{s}));
                        end
                    case 'FNIRT'
                        maskAnat = fullfile(BGDir,'FSL','MNI',subj_name{s},sprintf('%s_MNI_FNIRT.nii',subj_name{s})); %anatomical aligned to MNI
                    case 'native'
                        if strcmp(imageType,'anat')
                            maskAnat = fullfile(anatomicalDir,subj_name{s},sprintf('%s_anatomical.nii',subj_name{s})); %anatomical aligned to MNI
                        elseif strcmp(imageType,'func')
                            maskAnat = fullfile(regDir,sprintf('mask_%s.nii',subj_name{s})); %functional mask - used for searchlight
                        end
                end
                % Save region file as nifti
                dircheck(fullfile(BGDir,'regionMasks'));
                cd(fullfile(BGDir,'regionMasks'));
  
                region_saveasimg(R{roi},maskAnat,'name',sprintf('%s_%s.nii',R{roi}.name,regType));
            end
        end
        
    case 'DICE_intersubject_overlap'
        template={'SPM','FNIRT'}; 
        templateType=[1 2]; % both
        sn=[3:20];
        roi=[1:3];
        hemi=[1:2];
        vararginoptions(varargin,{'sn','templateType','roi','hemi'});
        
        DD=[];
        for t=templateType
        for h=hemi
            for r=roi
                % all combinations
                indMat=indicatorMatrix('allpairs',sn);
                BGsubjDir=fullfile(BGDir,sprintf('%s',template{t}));
                
                for c=1:size(indMat,1)
                    [i j]=find(indMat(c,:));
                    s1=i(2); s2=j(2);
                    subj1=sprintf(fullfile(BGsubjDir,subj_name{sn(s1)},sprintf('%s_%s_%s_MNI_%s.nii',subj_name{sn(s1)},regname_BG{r},hem{h},template{t})));
                    subj2=sprintf(fullfile(BGsubjDir,subj_name{sn(s2)},sprintf('%s_%s_%s_MNI_%s.nii',subj_name{sn(s2)},regname_BG{r},hem{h},template{t})));
                    subjV1=spm_vol(subj1);
                    subjRV1=spm_read_vols(subjV1);
                    subjV2=spm_vol(subj2);
                    subjRV2=spm_read_vols(subjV2);
                    % calculate overlap
                    ind1=find(subjRV1);
                    ind2=find(subjRV2);
                    overlap=intersect(ind1,ind2);
                    dice=2*length(overlap)/(length(ind1)+length(ind2));
                    D.dice=dice;
                    D.roi=r;
                    D.hemi=h;
                    D.subj1=s1;
                    D.subj2=s2;
                    D.templateType=t;
                    DD=addstruct(DD,D);
                end
            end
        end
        end
        % plot the overlap
        
        figure
        for f=1:2
            subplot(1,2,f)
            plt.line(DD.roi,DD.dice,'split',DD.templateType,'subset',DD.hemi==f,'leg',template,'leglocation','northeast');
            ylabel('Intersubject volume overlap - dice coefficient');
            set(gca,'XTickLabel',regname_BG(1:3));
            title(sprintf('%s DICE',hem{f}));
        end
        plt.match('y');  
    case 'DICE_group_overlap'
        template={'SPM','FNIRT'}; 
        templateType=[1,2]; % 1, 2 or both
        sn=[3:20];
        roi=[1:3];
        hemi=[1:2];
        vararginoptions(varargin,{'sn','templateType','roi','hemi'});
        
        DD=[];
        for t=templateType
        for h=hemi
            for r=roi
                % group file
                BGgroupDir=fullfile(BGDir,template{t},'avrg');
                group=fullfile(BGgroupDir,sprintf('subjAvrg_anat_MNI_%s_%s_MASK.nii',regname_BG{r},hem{h}));
                groupV=spm_vol(group);
                groupRV=spm_read_vols(groupV);
                % all subjects
                indMat=indicatorMatrix('identity',sn);
                BGsubjDir=fullfile(BGDir,sprintf('%s',template{t}));     
                % subject file
                for c=1:size(indMat,1)
                    i=find(indMat(c,:));
                    s1=i; 
                    subj=sprintf(fullfile(BGsubjDir,subj_name{sn(s1)},sprintf('%s_%s_%s_MNI_%s.nii',subj_name{sn(s1)},regname_BG{r},hem{h},template{t})));
                    subjV=spm_vol(subj);
                    subjRV=spm_read_vols(subjV);
                    % calculate overlap - group to subj
                    ind1=find(subjRV);
                    ind2=find(groupRV);
                    overlap=intersect(ind1,ind2);
                    dice=2*length(overlap)/(length(ind1)+length(ind2));
                    D.dice=dice;
                    D.roi=r;
                    D.hemi=h;
                    D.subj=s1;
                    D.templateType=t;
                    DD=addstruct(DD,D);
                end
            end
        end
        end
        % plot the overlap
        figure
        for f=1:2
            subplot(1,2,f)
            plt.line(DD.roi,DD.dice,'split',DD.templateType,'subset',DD.hemi==f,'leg',template,'leglocation','northeast');
            ylabel('Intersubject volume overlap - dice coefficient');
            set(gca,'XTickLabel',regname_BG(1:3));
            title(sprintf('%s DICE',hem{f}));
        end
        plt.match('y'); 
 
    case 'ROI_BG_timeseries' % ------------- timeseries from GLM
        glm=2;
        vararginoptions(varargin,{'sn','sessN','glm'});
        
        pre=4;          % How many TRs before the trial onset
        post=16;        % How many TRs after the trial onset
        T=[];
        for s=sn
            fprintf('Extracting the onsets and events for subject %s, glm %d and session %d\n',subj_name{s},glm,sessN);
            load(fullfile(glmSessDir{sessN},subj_name{s},'SPM.mat'));
            
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data', subj_name{s}));
            load(fullfile(regDir,[subj_name{s} '_BG_native_regions.mat']));      % This is made in case 'ROI_define'
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
    case 'PSC_create' % ------ percent signal change - create in case not defined in sml1_imana_dist/stability and sml1_imana_repsup
        % calculate psc for trained and untrained sequences - based on betas
        vararginoptions(varargin,{'sn','sessN'});
        name={'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','TrainSeq','UntrainSeq'};
        for s=sn
            cd(fullfile(glmSessDir{sessN}, subj_name{s}));
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
                outname=sprintf('psc_sess%d_%s.nii',sessN,name{con}); % ,subj_name{s}
                
                formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);
                
                spm_imcalc_ui(P,outname,formula,...
                    {0,[],spm_type(16),[]});        % Calculate percent signal change
            end;
            fprintf('Subject %d sess %d: %3.3f\n',s,sessN,h);
        end;
    case 'PSC_create_RepSup'
        % calculate psc for trained and untrained sequences (1st/2nd) - based on betas
        vararginoptions(varargin,{'sn','sessN'});
        name={'TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
        for s=sn
            cd(fullfile(glmFoSExDir{sessN}, subj_name{s}));
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
                outname=sprintf('psc_sess%d_%s.nii',sessN,name{con}); % ,subj_name{s}
                
                formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);
                
                spm_imcalc_ui(P,outname,formula,...
                    {0,[],spm_type(16),[]});        % Calculate percent signal change
            end;
            fprintf('Subject %d sess %d: %3.3f\n',s,sessN,h);
        end;
        
    case 'BG_betas' % -------- save betas per glm (per session) - for psc / dist
        sessN = 1;
        %sn  = [4:16];
        sn  = [4:9,11:16];
        %roi = [1:8];
        roi=[1,2,3,5,6,7]; % no thalamus
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        regType='subjAvrg'; % 1) native or group: 2)subjAvrg, 3)SPM (MNI, MNI_step2, SPM or FNIRT) - DO NATIVE OR GROUP_AVERAGE ONLY
        vararginoptions(varargin,{'sn','sessN','roi','type','regType'});
        
        for ss=sessN
            switch(type)
                case 'new'
                    T=[];
                case 'add'
                    T=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            end
            
            % harvest
            for s=sn % for each subj
                fprintf('\nSubject: %d\n',s) % output to user
                cd (fullfile(glmSessDir{ss},subj_name{s})); 
                
                if strcmp(regType,'native')
                    % load files
                    load(fullfile(glmSessDir{sessN}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                    load(fullfile(regDir,[subj_name{s} sprintf('_BG_%s_regions.mat',regType)]));        % load subject's region parcellation (R)
                    
                    % Add a few extra images
                    %----task against rest
                    O{1}=sprintf('psc_sess%d_TrainSeq.nii',ss); %psc trained
                    O{2}=sprintf('psc_sess%d_UntrainSeq.nii',ss); %psc untrained
                    O{3}='spmT_TrainSeq.nii';
                    O{4}='spmT_UntrainSeq.nii';
                else % subjAvrg (or SPM)
                    load(fullfile(regDir,'subjAvrg_MNI_BG_regions.mat'));        % load subject's region parcellation (R)
                    O{1}=sprintf('wpsc_sess%d_TrainSeq_MNI_BG.nii',ss); %psc trained
                    O{2}=sprintf('wpsc_sess%d_UntrainSeq_MNI_BG.nii',ss); %psc untrained
                    O{3}='wspmT_TrainSeq_MNI_BG.nii';
                    O{4}='wspmT_UntrainSeq_MNI_BG.nii';
                end
                
                %P=SPM.Vbeta(SPM.xX.iC);
                V = SPM.xY.VY; 
                oP=spm_vol(char(O));
                
                for r = roi % for each region
                    % get raw data for voxels in region
                    Y = region_getdata(V,R{r});  % Data Y is N x P
                    data = region_getdata(oP,R{r}); % from added images
                    % voxel position
                    S.volcoord = {R{r}.data'};   
                    %        % estimate region betas
                    if strcmp(regType,'native')
                            [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                            S.betaW                   = {betaW};                             % multivariate pw
                            S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                            S.betaRAW                 = {beta};
                            S.resMS                   = {resMS};
                    end
                    % info from maps for psc
                    S.psc_train    = {data(1,:)};
                    S.psc_untrain  = {data(2,:)};
                    S.Tmap_train   = {data(3,:)};
                    S.Tmap_untrain = {data(4,:)};
                    
                    S.SN                      = s;
                    S.region                  = r;
                    T = addstruct(T,S);
                    fprintf('%d.',r)
                end
            end
            % save T
            save(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)),'-struct','T');
            fprintf('\n');
        end
    case 'BG_stats'
        sessN = [1:3];
        sn  = [4:18];
        roi = [1:3,5:7];
        atlasType = 'MNI_step2';
        betaChoice = 'uni'; % uni, multi or raw
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type','atlasType'});
        
        for ss=sessN;
            T = load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',atlasType,ss))); % loads region data (T)
            
            % output structures / or load from before
            switch(type)
                case 'new'
                    To=[];
                case 'add'
                    To=load(fullfile(BGDir,sprintf('stats_BG_%sPW_sess%d.mat',betaChoice,ss)));
            end
            
            % do stats
            for s = sn % for each subject
                D = load(fullfile(glmSessDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                fprintf('\nSubject: %d session: %d\n',s,ss)
                num_run = numruns_task_sess;
                
                for r = roi % for each region
                    beta_seqType=[];
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
                    
                    % % To structure stats (all seqNumb - 12 conditions)
                    % crossval second moment matrix
                    
                    if strcmp(atlasType,'native')
                    [G,Sig]     = pcm_estGCrossval(betaW(1:(12*num_run),:),D.run,D.seqNumb);
                    
                    So.IPM      = rsa_vectorizeIPM(G);
                    So.Sig      = rsa_vectorizeIPM(Sig);
                    % squared distances
                    So.RDM      = rsa.distanceLDC(betaW,D.run,D.seqNumb);
                    
                    % calculate G and dist
                    Do = calcDist(D,betaW,G);
                    end
                    % stats from additional images
                    So.psc_train = nanmean(S.psc_train{:});
                    So.psc_untrain = nanmean(S.psc_untrain{:});
                    
                    % indexing fields
                    So.SN       = s;
                    So.region   = r;
                    So.regSide  = regSide(r);
                    So.regType  = regType(r);
                    
                    % data structure
                    To          = addstruct(To,So); % indexing fields, other images
                    To          = addstruct(To,Do); % distances
                    
                end; % each region
            end; % each subject
            
            % % save - stats data and simulations
            save(fullfile(BGDir,sprintf('stats_BG_%s_%sPW_sess%d.mat',atlasType,betaChoice,ss)),'-struct','To');
            
            fprintf('\nDone.\n')
        end
    case 'pattern_consist'
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 

        % enter sn, region,  beta: 0=betaW, 1=betaU, 2=raw betas
        % (1) Set parameters
        sn  = [4:18];
        roi = [1:3,5:7];
        sessN=1;
        removeMean = 'yes'; % are we removing pattern means for patternconsistency?
        betaChoice='multiPW'; % multiPW, uniPW, raw
        vararginoptions(varargin,{'sn','roi','removeMean','sessN','betaChoice'});
        
        numRun=numruns_task_sess;
        numSeq=numel(num_seq);
        
        if strcmp(removeMean,'yes')
             rmean = 1; % we are removing the mean
        else rmean = 0; % we are keeping the mean (yields higher consistencies but these are biased)
        end

        RR=[];
        %========%
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_BG_native_sess%d',ss))); % loads in struct 'T'
            for r=roi
                for s=sn
                    S = getrow(T,(T.SN==s & T.region==r));
                    
                    switch(betaChoice)
                        case 'raw'
                            beta  = S.betaRAW{1}; 
                        case 'uniPW' 
                            beta  = S.betaUW{1}; 
                        case 'multiPW'
                            beta  = S.betaW{1}; 
                    end
                    % make vectors for pattern consistency func
                    partVec=kron([1:numRun]',ones(numSeq,1));
                    condVec=kron(ones(numRun,1),[1:numSeq]');

                    % calculate the pattern consistency
                    R.consist   = rsa_patternConsistency(beta,partVec,condVec,'removeMean',rmean);
                    R.consist_crossval = rsa_patternConsistency_crossval(beta,partVec,condVec,'removeMean',rmean);
                    R.sn =s;
                    R.sessN=ss;
                    R.roi=r;
                    if r<4
                        R.regType=r;
                        R.hemi=1;
                    else
                        R.regType=r-4;
                        R.hemi=2;
                    end
                    RR=addstruct(RR,R);
                end
            end
        end
                
        figure
        subplot(2,3,[1:3]);
        plt.bar(RR.regType,RR.consist,'split',RR.hemi,'leg',{'Contra','Ipsi'},'leglocation','northeast');
        ylabel(sprintf('Pattern consistency %s',betaChoice));
        title('Pattern consistency');
        xlabel('ROIs');
        
        for rr=1:3
            subplot(2,3,rr+3)
            plt.hist(RR.consist,'subset',RR.regType==rr,'split',RR.hemi,'leg',{'contra','ipsi'},'leglocation','north');
            title(sprintf('%s',regname_BG{rr}));
        end
        
        figure
        subplot(2,3,[1:3]);
        plt.bar(RR.regType,RR.consist_crossval,'split',RR.hemi,'leg',{'Contra','Ipsi'},'leglocation','northeast');
        title('Pattern consistency - crossvalidated');
        ylabel(sprintf('Pattern consistency %s',betaChoice));
        xlabel('ROIs');
        
        for rr=1:3
            subplot(2,3,rr+3)
            plt.hist(RR.consist_crossval,'subset',RR.regType==rr,'split',RR.hemi,'leg',{'contra','ipsi'},'leglocation','north');
            title(sprintf('%s',regname_BG{rr}));
        end
        
        
         for p=1:2 % normal and crossvalidated
            if p==1
                var=RR.consist;
                type='normal';
            else
                var=RR.consist_crossval;
                type='crossval';
            end
            
            % split by hemisphere
            figure
            subplot(2,3,[1:3]);
            plt.bar(RR.regType,var,'split',RR.hemi,'leg',{'Contra','Ipsi'},'leglocation','northeast');
            ylabel(sprintf('Pattern consistency %s',betaChoice));
            title(sprintf('Pattern consistency %s',type));
            xlabel('ROIs');
            
            for rr=roi
                subplot(2,3,rr+3)
                plt.hist(var,'subset',RR.regType==rr,'split',RR.hemi,'leg',{'contra','ipsi'},'leglocation','north');
                title(sprintf('%s',regname_BG{rr}));
            end
            
            % split by session
            figure
            subplot(2,3,[1:3]);
            plt.bar(RR.regType,var,'subset',RR.hemi==1,'split',RR.sessN,'leg',{'sess1','sess2','sess3','sess4'},'leglocation','northeast');
            ylabel(sprintf('Pattern consistency %s',betaChoice));
            title(sprintf('Pattern consistency %s',type));
            xlabel('ROIs');
            
            for rr=roi
                subplot(2,3,rr+3)
                plt.hist(var,'subset',RR.hemi==1,'split',RR.sessN,'leg',{'sess1','sess2','sess3','sess4'},'leglocation','north');
                title(sprintf('%s',regname_BG{rr}));
            end
            
            
        end

    case 'pattern_consist_crossval'
            reg = [1:8]; 
        sn  = [1:9,11,12];
        sessN = [1:4];
        remove_mean = 1; % subtract
        betaChoice = 'mw'; % uw / mw / raw

        vararginoptions(varargin,{'roi','sn','remove_mean','sessN','betaChoice'});
                
        partitions = [1:2:numruns_task_sess; 2:2:numruns_task_sess];
        numRuns    = [1:numruns_task_sess];
        numConds   = num_seq;
        conds   = repmat([numConds],1,length(numRuns))';
        runNums = kron([numRuns],ones(1,length(numConds)))';
        
        C=[];
        for ss = sessN
            D   = load(fullfile(regDir,sprintf('betas_sess%d.mat',ss)));
            splitcorrs = [];
            for roi = reg;
                T   = getrow(D,D.region==reg(roi));
                
                for s = 1:numel(sn) % for each subject
                    sbeta=[]; prepBetas=[];
                    t = getrow(T,T.SN==sn(s));
                    
                    prepBetas = [];
                    switch (betaChoice)
                        case 'uw'
                            beta = t.betaUW{1};
                        case 'mw'
                            beta = t.betaW{1};
                        case 'raw'
                            beta = t.betaRAW{1};
                    end
                    
                    if remove_mean
                        for r=numRuns
                            sbeta(runNums==r,:) = bsxfun(@minus,beta(runNums==r,:),mean(beta(runNums==r,:)));
                        end
                    else
                        sbeta=beta;
                    end
                    % prep betas (harvest and subtract partition mean)
                    for i = 1:size(partitions,1)
                        partitionIdx = logical(ismember(runNums,partitions(i,:)));
                        condIdx{i}   = conds(partitionIdx);
                        prepBetas{i} = sbeta(partitionIdx,:);
                    end
                    
                    % correlate patterns across partitions, both within and across
                    % conditions
                    for c1 = numConds % for each condition
                        % condition mean activity pattern for the odd run partition
                        oddCon   = condIdx{1}==c1;
                        oddBetas = mean(prepBetas{1}(oddCon,:));
                        % condition mean activity pattern for the even run partition
                        evenCon   = condIdx{2}==c1;
                        evenBetas = mean(prepBetas{2}(evenCon,:));
                        % correlate condition patterns across partitions
                        tmp = corrcoef(evenBetas,oddBetas);
                        splitcorrs(c1) = tmp(1,2);
                    end
                    
                    Corr.sn = sn(s);
                    Corr.reg = roi;
                    Corr.trainCorr = mean(splitcorrs(:,[1:6]));
                    Corr.untrainCorr = mean(splitcorrs(:,[7:12]));
                    Corr.overallCorr = mean(splitcorrs);
                    Corr.sessN = ss;
                    
                    C = addstruct(C,Corr);         
                end
            end
        end
        
    case 'CALC_psc' % ------------------- PSC, dist (Mahalanobis, correlation) ------------
        sn = [4:9,11:18];
        roi = [1:3,5:7];
        sessN = [1:4];
        regType = 'native'; % subjAvrg

        vararginoptions(varargin,{'sn','roi','seq','sessN','regType'});

        Stats = [];
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss))); % loads region data (T)
            for s=1:numel(sn)
                for r=roi
                   
                    S.psc(1,:)=nanmean(T.psc_train{(T.region==r & T.SN==sn(s)),:});
                    S.psc(2,:)=nanmean(T.psc_untrain{(T.region==r & T.SN==sn(s)),:});
                    S.seqType=[1;2];
                    S.sn=ones(size(S.seqType))*sn(s);
                    S.roi=ones(size(S.seqType))*r;
                    S.sessN=ones(size(S.seqType))*ss;
                    
                    Stats=addstruct(Stats,S);
                end
            end
        end
           
        % save structure
        save(fullfile(BGDir,sprintf('psc_%s_BG.mat',regType)),'-struct','Stats');
    case 'PLOT_psc'
        roi=[1:3,5:7];
        sessN=[1:4];
        regType='subjAvrg'; % subjAvrg or native
        
        vararginoptions(varargin,{'roi','sessN','regType'});
        
        T=load(fullfile(BGDir,sprintf('psc_%s_BG.mat',regType)));
        
        sub=0;
        figure
        for r=roi
            sub=sub+1;
            subplot(1,numel(roi),sub)
            if numel(sessN)==4
                plt.line([T.sessN>3 T.sessN],T.psc,'split',T.seqType,'subset',T.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(T.sessN,T.psc,'split',T.seqType,'subset',T.roi==r&T.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            end
            plt.match('y');
            title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Percent signal change');
                xlabel('Session');
            else
                ylabel('');
            end
        end
    case 'STATS_psc'
        roi=[1:3];
        sessN=[1:4];
        regType='native'; % native or subjAvrg
        vararginoptions(varargin,{'roi','regType','sessN'});
        
        T=load(fullfile(BGDir,sprintf('psc_%s_BG.mat',regType)));
        
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n PSC in %s \n',regname_BG{r});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.seqType],{'session','seqType'},'subset',T.roi==r);
        end
        for r=roi
            for ss=sessN
                fprintf('\n post-hoc t-test on the effect of seqType in sess %d in %s \n',ss,regname_BG{r});
                ttestDirect(T.psc,[T.seqType T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);
            end
        end

    case 'CALC_dist'
        sn = [4:9,11:18];
        roi = [1:3,5:7];
        sessN = [1:4];
        betaChoice='multiPW';
        regType='native';
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice','fig','regType'});
        
        Stats = [];
        for ss=sessN
            D = load(fullfile(BGDir,sprintf('stats_BG_%s_%s_sess%d.mat',regType,betaChoice,ss))); % loads region data (D)
            for s=1:numel(sn)
                for r=roi
                    S.dist_train=D.dist_train(D.SN==sn(s) & D.region==r);
                    S.dist_untrain=D.dist_untrain(D.SN==sn(s) & D.region==r);
                    S.dist_cross=D.dist_cross(D.SN==sn(s) & D.region==r);
                    S.sn=sn(s);
                    S.roi=r;
                    S.sessN=ss;
                    Stats=addstruct(Stats,S); 
                end
            end
        end
        % save structure
        save(fullfile(BGDir,sprintf('dist_%s_BG.mat',regType)),'-struct','Stats');
    case 'PLOT_dist'
        roi=[1:3,5:7];
        sessN=[1:4];
        regType='native'; % native or subjAvrg
        
        vararginoptions(varargin,{'roi','sessN','regType'});
        
        T=load(fullfile(BGDir,sprintf('dist_%s_BG.mat',regType)));
        
        % create structure with trained / untrained together, seqType
        D.dist=[T.dist_train;T.dist_untrain];
        D.seqType=[ones(size(T.dist_train));ones(size(T.dist_train))*2];
        D.sn=[T.sn;T.sn];
        D.roi=[T.roi;T.roi];
        D.sessN=[T.sessN;T.sessN];
        
        sub=0;
        figure
        for r=roi
            sub=sub+1;
            subplot(1,numel(roi),sub)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.dist),'split',D.seqType,'subset',D.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.dist),'split',D.seqType,'subset',D.roi==r&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==2
                plt.line(D.sessN,ssqrt(D.dist),'split',D.seqType,'subset',D.roi==r&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Mahalanobis distance');
                xlabel('Session');
            else
                ylabel('');
            end
        end
    case 'PLOT_dist_seqType'
        roi=[1:3,5:7];
        sessN=[1:4];
        regType='native'; % native or subjAvrg
        
        vararginoptions(varargin,{'roi','sessN','regType'});
        
        D=load(fullfile(BGDir,sprintf('dist_%s_BG.mat',regType)));
        figure
        sub=0;
        for r=roi
            sub=sub+1;
            subplot(1,numel(roi),sub)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.dist_cross),'subset',D.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.dist_cross),'subset',D.roi==r&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==2
                plt.line(D.sessN,ssqrt(D.dist_cross),'subset',D.roi==r&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Mahalanobis distance');
                xlabel('Session');
            else
                ylabel('');
            end
        end
    
    case 'CALC_corrdist'
        reg = [1:3,5:7];
        sn  = [4:9,11:16];
        sessN = [1:4];
        regType='native'; % 1) native or group: 2)subjAvrg, 3)SPM (MNI, MNI_step2, SPM or FNIRT) - DO NATIVE OR GROUP_AVERAGE ONLY
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        vararginoptions(varargin,{'sn','reg','sessN','subtract_mean','regType'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
            D   = load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
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
                    
                    ST.corr_seqType = sum(sum(triu(C(7:12,1:6)))/(6*5));
                    ST.sessN=ss;
                    ST.roi=roi;
                    ST.sn=s;
                    STAll=addstruct(STAll,ST);
                end
            end
        end
        save(fullfile(BGDir,sprintf('corrDist_%s_ROI.mat',regType)),'-struct','SAll');
        save(fullfile(BGDir,sprintf('corrDist_seqType_%s_ROI.mat',regType)),'-struct','STAll');
    case 'PLOT_corrDist'
        roi=[1:3,5:7];
        sessN=[1:4];
        regType='native';
        
        vararginoptions(varargin,{'roi','sessN','regType'});
        
        T=load(fullfile(BGDir,sprintf('corrDist_%s_ROI.mat',regType)));
        
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([T.sessN>3 T.sessN],ssqrt(T.corrDist),'split',T.seqType,'subset',T.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            elseif numel(sessN)==3
                plt.line(T.sessN,ssqrt(T.corrDist),'split',T.seqType,'subset',T.roi==r&T.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Correlation distance');
                xlabel('Session');
            else
                ylabel('');
            end
        indx=indx+1;    
        end
    case 'PLOT_corrDist_seqType'
        roi=[1:3,5:7];
        sessN=[1:4];
        regType='native';
        
        vararginoptions(varargin,{'roi','sessN','regType'});
        
        D=load(fullfile(BGDir,sprintf('corrDist_seqType_%s_ROI.mat',regType)));
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.corr_seqType),'subset',D.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.corr_seqType),'subset',D.roi==r&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==2
                plt.line(D.sessN,ssqrt(D.corr_seqType),'subset',D.roi==r&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            end
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Correlation distance betweeen seqTypes');
                xlabel('Session');
            else
                ylabel('');
            end
        indx=indx+1;    
        end
    
    case 'BG_betas_RepSup' % --------- save betas for each FoSEx glm - repetition suppression
        sessN = 1;
        sn  = [4:16];    
        roi = [1:3,5:7];
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        atlasType='native';
        vararginoptions(varargin,{'sn','sessN','roi','type','atlasType'});
        
         switch(type)
            case 'new'
                T=[];
            case 'add'
                T=load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',sessN)));
        end
            
        % harvest
        for s=sn % for each subj
            fprintf('\nSubject: %d\n',s) % output to user
            
            % load files
            load(fullfile(glmFoSExDir{sessN}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            load(fullfile(regDir,[subj_name{s} '_BG_' atlasType '_regions.mat']));        % load subject's region parcellation (R)
            V = SPM.xY.VY; 
            glmSubjDir = fullfile(glmFoSExDir{sessN},subj_name{s});
            cd(glmSubjDir);
            
            % Add a few extra images
            %----task against rest
            O{1}=sprintf('psc_sess%d_TrainSeq_1st.nii',sessN); %psc trained - 1st execution
            O{2}=sprintf('psc_sess%d_UntrainSeq_1st.nii',sessN); %psc untrained - 1st execution
            O{3}=sprintf('psc_sess%d_TrainSeq_2nd.nii',sessN); %psc trained - 2nd execution
            O{4}=sprintf('psc_sess%d_UntrainSeq_2nd.nii',sessN); %psc trained - 2nd execution
            oP=spm_vol(char(O));

            
            for r = roi % for each region
                % get raw data for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P 
                data = region_getdata(oP,R{r}); % from added images
                
                % estimate region betas
                [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                S.betaW                   = {betaW};                             % multivariate pw
                S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw 
                S.betaRAW                 = {beta};
                S.resMS                   = {resMS};
                
                 % info from maps for surface
                S.psc_train_1st   = {data(1,:)}; 
                S.psc_untrain_1st = {data(2,:)};
                S.psc_train_2nd   = {data(3,:)}; 
                S.psc_untrain_2nd = {data(4,:)};
                
                S.SN                      = s;
                S.region                  = r;
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
        end
        % save T
        save(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',sessN)),'-struct','T'); 
        fprintf('\n');      
    case 'BG_stats_RepSup'
        sessN = 1;
        sn  = [4:16];
        roi = [1:3,5:7];
        betaChoice = 'multi'; % uni, multi or raw
        type='new';
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type'});
        
        T = load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',sessN))); % loads region data (T)
        
        % output structures
        switch(type)
            case 'new'
                To=[];
            case 'add'
                To=load(fullfile(BGDir,sprintf('stats_FoSEx_%sPW_sess%d.mat',betaChoice,sessN)));
        end
        
        % do stats
        for s = sn % for each subject
            D = load(fullfile(glmFoSExDir{sessN}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
            fprintf('\nSubject: %d session: %d\n',s, sessN)
            num_run = numruns_task_sess;
            
            for r = roi % for each region
                S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                fprintf('%d.',r)
                
                for exe = 1:2   % FoSEx
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
                    So.psc_train_1st = nanmean(S.psc_train_1st{:});
                    So.psc_train_2nd = nanmean(S.psc_train_2nd{:});
                    So.psc_untrain_1st = nanmean(S.psc_untrain_1st{:});
                    So.psc_untrain_2nd = nanmean(S.psc_untrain_2nd{:});
                    % indexing fields
                    So.SN       = s;
                    So.region   = r;
                    So.regSide  = regSide(r);
                    So.regType  = regType(r);
                    So.FoSEx    = exe;
                    To          = addstruct(To,So);
                    
                end; % FoSEx
            end; % each region
        end; % each subject

        % % save
        save(fullfile(BGDir,sprintf('stats_FoSEx_%sPW_sess%d.mat',betaChoice,sessN)),'-struct','To');
        fprintf('\nDone.\n') 
    
    case 'CALC_repsup' % ----------- save psc / dist for FoSEx - for cortex in sml1_imana_repsup
        
        sn = [4:9,11:22];
        roi = [1:3,5:7];
        sessN = 1:4;
        betaChoice = 'multiPW'; % uniPW or multiPW

        vararginoptions(varargin,{'sn','roi','seq','sessN','betaChoice'});

        Stats = [];
        
        for ss = sessN % do per session number
            D = load(fullfile(BGDir,sprintf('stats_FoSEx_%s_sess%d.mat',betaChoice,ss))); % loads region data (D)
            T = load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',ss))); % loads region data (T)
            
            runs=1:numruns_task_sess;
            conditionVec = kron(ones(numel(runs),1),[1:12]');       
            indx_train = 1:6;
            indx_untrain = 7:12;
            
            for s=1:length(sn)
                SI = load(fullfile(glmFoSExDir{ss}, subj_name{sn(s)}, 'SPM_info.mat'));   % load subject's trial structure
                for r=roi
                    for exe=1:2;    % FoSEx
                       D_exe = getrow(D,D.FoSEx==exe); 
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
                        
                        AllDist = ssqrt(rsa_squareRDM(D_exe.RDM(D_exe.region==r & D_exe.SN==sn(s),:)));
                        SeqTrain = triu(AllDist(indx_train,indx_train));
                        SeqUntrain = triu(AllDist(indx_untrain,indx_untrain));
                        SeqCross = triu(AllDist(indx_train,indx_untrain));
                        SeqTrainAll = SeqTrain(SeqTrain~=0);
                        SeqUntrainAll = SeqUntrain(SeqUntrain~=0);
                        SeqCrossAll = SeqCross(SeqCross~=0);
                        
                        S.sn=sn(s);
                        S.roi=r;
                        S.dist_train=mean(SeqTrainAll);
                        S.dist_untrain=mean(SeqUntrainAll);
                        S.dist_cross=mean(SeqCrossAll);
                        S.beta_train=mean(mean(C.beta_seq_train));
                        S.beta_untrain=mean(mean(C.beta_seq_untrain));
                        if exe==1
                            S.psc_train=D.psc_train_1st(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            S.psc_untrain=D.psc_untrain_1st(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                        else
                            S.psc_train=D.psc_train_2nd(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                            S.psc_untrain=D.psc_untrain_2nd(D.SN==sn(s)&D.region==r&D.FoSEx==exe);
                        end
                        S.sessN=ss;
                        S.FoSEx=exe;
                        Stats=addstruct(Stats,S);
                    end; % FoSEx
                end; % roi
            end; % sn
        end

       % % save
        save(fullfile(BGDir,sprintf('RepSup_stats_%s.mat',betaChoice)),'-struct','Stats');    
    case 'PLOT_repsup_psc'
        
        betaChoice = 'multiPW';
        roi=[1:3,5:7];
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(BGDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        for st=1:2
            if st==1
                style=styTrained_exe;
                var=S.psc_train;
            else
                style=styUntrained_exe;
                var=S.psc_untrain;
            end
            figure
            sub=0;
            for f = roi
                sub=sub+1;
                subplot(1,numel(roi),sub);
                plt.line([S.sessN>3 S.sessN],var,'split',S.FoSEx,'subset',S.roi==f,'style',style,'leg',{'1st','2nd'});
                plt.match('y');
                %  plt.match('y');
                if f==1
                    ylabel(sprintf('Psc %s',SeqType{st}));
                    xlabel('Session');
                else
                    ylabel('');
                end
                title(sprintf('%s',regname_BG{rem(f,4)}));
            end
        end
    case 'PLOT_repsup_dist'
       
       betaChoice = 'multiPW';
       roi=[1:3,5:7];
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(BGDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
 
       figure
       sub=0;
       for f = roi
           sub=sub+1;
           subplot(1,numel(roi),sub);
           plt.line([S.sessN>3 S.sessN],S.dist_train,'split',S.FoSEx,'subset',S.roi==f,'style',stySeqType_exe,'leg',{'1st','2nd'});
           drawline(0,'dir','horz');
      %     plt.match('y');
           if f==1
               ylabel('Distances trained'); xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname_BG{rem(f,4)}))
       end
      
       figure
       sub=0;
       for f = roi
           sub=sub+1;
           subplot(1,numel(roi),sub);
           plt.line([S.sessN>3 S.sessN],S.dist_untrain,'split',S.FoSEx,'subset',S.roi==f,'style',stySeqType_exe,'leg',{'1st','2nd'});
           drawline(0,'dir','horz');
        %   plt.match('y');
           if f==1
               ylabel('Distances untrained'); xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname_BG{rem(f,4)}))
       end
       
       figure
       sub=0;
       for f = roi
           sub=sub+1;
           subplot(1,numel(roi),sub);
           plt.line([S.sessN>3 S.sessN],S.dist_cross,'split',S.FoSEx,'subset',S.roi==f,'style',stySeqType_exe,'leg',{'1st','2nd'});
           drawline(0,'dir','horz');
        %   plt.match('y');
           if f==1
               ylabel('Distance between seq sets'); xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname_BG{rem(f,4)}))
       end
    case 'STATS_repsup_psc'
        betaChoice = 'multiPW';
        roi=[1:3,5:7];
        sessN=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        D = load(fullfile(BGDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        
        T.psc = [D.psc_train; D.psc_untrain];
        T.seqType = [ones(size(D.psc_train));ones(size(D.psc_train))*2];
        T.sn = [D.sn; D.sn];
        T.roi = [D.roi; D.roi];
        T.sessN = [D.sessN; D.sessN];
        T.FoSEx = [D.FoSEx; D.FoSEx];
        
        
        % ANOVA - session x seqType - per region
        for r=roi
            fprintf('\n PSC in %s \n',regname_BG{rem(r,4)});
            anovaMixed(T.psc,T.sn,'within',[T.sessN T.FoSEx T.seqType],{'session','repsup','seqType'},'subset',T.roi==r);
        end
    case 'PLOT_psc_sess1'
       betaChoice = 'multiPW';
       roi=[1,3];
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(BGDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
           
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
           set(gca,'XTick',[2 5 8 11 14.5 18 21],'XTickLabel',regname_BG(roi));
           xlabel('ROI');      
           
           for r=roi
               fprintf('%s \t',regname_BG{r});
               ttestDirect(SS.psc,[SS.FoSEx SS.sn],2,'paired','subset',SS.roi==r&SS.sessN==1);
           end
        
    case 'CALC_percent_repsup' % ------- calculate amount of RS 
        
        betaChoice = 'multiPW'; % or uniPW
        roi=[1:3,5:7];
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(BGDir,sprintf('RepSup_stats_%s.mat',betaChoice)));
        
        A = getrow(S,S.FoSEx==1);
        
        P.psc_beta = [S.beta_train(S.FoSEx==2)./S.beta_train(S.FoSEx==1);S.beta_untrain(S.FoSEx==2)./S.beta_untrain(S.FoSEx==1)];
        P.psc_psc  = [S.psc_train(S.FoSEx==2)./S.psc_train(S.FoSEx==1);S.psc_untrain(S.FoSEx==2)./S.psc_untrain(S.FoSEx==1)];
        P.psc_dist = [S.dist_train(S.FoSEx==2)./S.dist_train(S.FoSEx==1);S.dist_untrain(S.FoSEx==2)./S.dist_untrain(S.FoSEx==1)];
        
        P.sn = [A.sn; A.sn];
        P.roi = [A.roi; A.roi];
        P.seq = [ones(size(A.sn)); ones(size(A.sn))*2];  % 1 for trained, 2 for untrained
        P.sessN = [A.sessN; A.sessN];
        
        % % save
        save(fullfile(BGDir,sprintf('RepSup_percent_%s.mat',betaChoice)),'-struct','P');
    case 'PLOT_percent_repsup_psc'
        
        betaChoice = 'multiPW';
        roi=[1:3,5:7];
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(BGDir,sprintf('RepSup_percent_%s.mat',betaChoice)));
        
        figure
        sub=0;
        for f = roi
            sub=sub+1;
            subplot(1,numel(roi),sub);
            plt.line(S.sessN,S.psc_psc,'split',S.seq,'subset',S.roi==f,'style',stySeq,'leg',{'trained','untrained'});
            drawline(1,'dir','horz');
          %  plt.match('y');
            if f==1
                ylabel('Repetition suppression magnitude (activation)');
                xlabel('Session');
            else
                ylabel('');
            end
            title(sprintf('%s',regname_BG{rem(f,4)}));
        end      
    case 'PLOT_percent_repsup_dist'
       
       betaChoice = 'multiPW';
       roi=[1:3,5:7];
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(BGDir,sprintf('RepSup_percent_%s.mat',betaChoice))); 
 
       figure
       sub=0;
       for f = roi
           sub=sub+1;
           subplot(1,numel(roi),sub);
           plt.line(S.sessN,S.psc_dist,'split',S.seq,'subset',S.roi==f,'style',stySeq,'leg',{'trained','untrained'});
           drawline(1,'dir','horz');
      %     plt.match('y');
           if f==1
               ylabel('Repetition suppression magnitude (distances)'); xlabel('Session');
           else
               ylabel('');
           end
           title(sprintf('%s',regname_BG{rem(f,4)}))
       end
     
    case 'CALC_subtract_repsup'
        % Graffton version
       betaChoice = 'multiPW';
       vararginoptions(varargin,{'betaChoice','roi'});
       
       S = load(fullfile(BGDir,sprintf('RepSup_stats_%s.mat',betaChoice))); 
 
       A = getrow(S,S.FoSEx==1);
       
       P.subtr_beta = [S.beta_train(S.FoSEx==1)-S.beta_train(S.FoSEx==2);S.beta_untrain(S.FoSEx==1)-S.beta_untrain(S.FoSEx==2)];
       P.subtr_psc = [S.psc_train(S.FoSEx==1)-S.psc_train(S.FoSEx==2);S.psc_untrain(S.FoSEx==1)-S.psc_untrain(S.FoSEx==2)];
       P.subtr_dist = [S.dist_train(S.FoSEx==1)-S.dist_train(S.FoSEx==2);S.dist_untrain(S.FoSEx==1)-S.dist_untrain(S.FoSEx==2)];
       P.sn = [A.sn; A.sn];
       P.roi = [A.roi; A.roi];
       P.seq = [ones(size(A.sn)); ones(size(A.sn))*2];  % 1 for trained, 2 for untrained
       P.sessN = [A.sessN; A.sessN];
       
       % % save
        save(fullfile(BGDir,sprintf('RepSup_subtract_%s.mat',betaChoice)),'-struct','P');    
    case 'PLOT_subtract_repsup_psc'
        betaChoice = 'multiPW';
        roi=[1:3,5:7];
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(BGDir,sprintf('RepSup_subtract_%s.mat',betaChoice)));

       figure
       sub=0;
       for f = roi
           sub=sub+1;
           subplot(1,numel(roi),sub);
           plt.line(S.sessN,S.subtr_psc,'split',S.seq,'subset',S.roi==f,'style',stySeq,'leg',{'Trained','Untrained'});
           drawline(0,'dir','horz');
           if f==1
               ylabel('Subtract 1st - 2nd activation')
           else
               ylabel('')
           end
           title(sprintf('%s',regname_BG{rem(f,4)}))
       end
    case 'PLOT_subtract_repsup_dist'
        betaChoice = 'multiPW';
        roi=[1:3,5:7];
        vararginoptions(varargin,{'betaChoice','roi'});
        
        S = load(fullfile(BGDir,sprintf('RepSup_subtract_%s.mat',betaChoice)));

       figure
       sub=0;
       for f = roi
           sub=sub+1;
           subplot(1,numel(roi),sub);
           plt.line(S.sessN,S.subtr_dist,'split',S.seq,'subset',S.roi==f,'style',stySeq,'leg',{'Trained','Untrained'});
           drawline(0,'dir','horz');
           if f==1
               ylabel('Subtract 1st - 2nd dist')
           else
               ylabel('')
           end
           title(sprintf('%s',regname_BG{rem(f,4)}))
       end
    case 'PLOT_subtract_psc_sessions'
        sessN=[1:3];
        seqType='trained';
        betaChoice='multiPW';
        regExcl=[2,6];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        S = load(fullfile(BGDir,sprintf('RepSup_subtract_%s.mat',betaChoice)));
        T=getrow(S,~ismember(S.roi,regExcl)&ismember(S.sessN,sessN));
        
        T.hemi=T.roi;
        T.hemi(T.roi<4,:)=1;
        T.hemi(T.roi>3,:)=2;
        
        switch(seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        figure
        for h=1:2
            subplot(1,2,h)
            plt.bar(T.roi,T.subtr_psc,'split',T.sessN,'subset',T.seq==st&T.hemi==h,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
            drawline(0,'dir','horz');
            ylabel('PSC repetition suppression - 1st-2nd');
            set(gca,'XTickLabel',regname_BG(1:3));
            title(sprintf('%s %s',seqType,hemName{h}));
            xlabel('ROI');
        end
    case 'PLOT_subtract_psc_seqType'
        sessN=[1:3];
        betaChoice='multiPW';
        regExcl=[2,6];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        S = load(fullfile(BGDir,sprintf('RepSup_subtract_%s.mat',betaChoice)));
        T=getrow(S,~ismember(S.roi,regExcl)&ismember(S.sessN,sessN));
        
        T.hemi=T.roi;
        T.hemi(T.roi<4,:)=1;
        T.hemi(T.roi>3,:)=2;
        
        for ss=sessN
            figure
            for h=1
                plt.bar(T.roi,T.subtr_psc,'split',T.seq,'subset',T.sessN==ss&T.hemi==h,'leg',{'trained','untrained'},'leglocation','northeast','style',stySeq);
                drawline(0,'dir','horz');
                ylabel('PSC repetition suppression - 1st-2nd');
                set(gca,'XTickLabel',regname_BG(1:3));
                title(sprintf('sess-%d %s',ss,hemName{h}));
                xlabel('ROI');
            end
        end
    case 'STATS_subtract_psc'
        betaChoice = 'multiPW';
        roi=[1:3,5:7];
        sessN=[1:4];
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(BGDir,sprintf('RepSup_subtract_%s.mat',betaChoice)));
        
        % main ANOVA - session x sequence
        for r=roi
            if r<4
                hem='LH';
            else
                hem='RH';
            end
            fprintf('\n PSC in %s %s \n',regname_BG{rem(r,4)},hem);
            anovaMixed(T.subtr_psc,T.sn,'within',[T.sessN T.seq],{'session','seqType'},'subset',T.roi==r);
        end
        
        % t-tests comparing sequence type per session
        for r=roi
            if r<4
                hem='LH';
            else
                hem='RH';
            end
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS subtraction in sess %d in %s %s \n',ss,regname_BG{rem(r,4)},hem);
                ttestDirect(T.subtr_psc,[T.seq T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end
        
        % t-tests comparing across sessions - separately for T / UT seq
        for st=1:2
            fprintf('\n %s sequence\n',seqLabel{st});
            for sessTr=1:3
                for r=roi
                    if r<4
                        hem='LH';
                    else
                        hem='RH';
                    end
                    fprintf('\n post-hoc t-test on comparison of sessions %d-%d in %s %s\n',sessTr,sessTr+1,regname_BG{rem(r,4)},hem);
                    ttestDirect(T.subtr_psc,[T.sessN T.sn],2,'paired','subset',T.roi==r&ismember(T.sessN,[sessTr,sessTr+1])&T.seq==st);
                end
            end
        end
       
    case 'CALC_corrDist_repsup'
        reg = [1:3,5:7];
        sn  = [4:9,11:22];
        sessN = [1:4];
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        vararginoptions(varargin,{'sn','reg','sessN','subtract_mean'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
            D   = load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',ss)));
            for roi = reg;
                for s = sn;
                    SI = load(fullfile(glmFoSExDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                    
                    for exe = 1:2
                        t   = getrow(D,D.region==roi & D.SN==s);
                        SI_exe = getrow(SI,SI.FoSEx==exe);
                        
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
                        [G,Sig]     = pcm_estGCrossval(data,SI_exe.run,SI_exe.seqNumb);
                        C=corr_crossval(G,'reg','minvalue');
                        C=rsa_squareRDM(C);
                        
                       % average trained dist
                        S.corrDist(1,:) = sum(sum(C(1:6,1:6)))/(6*5);
                        % average untrained dist
                        S.corrDist(2,:) = sum(sum(C(7:12,7:12)))/(6*5);
                        S.seqType=[1;2]; % trained, untrained
                        S.sessN=[ss;ss];
                        S.roi=[roi;roi];
                        S.sn=[s;s];
                        S.FoSEx=[exe;exe];
                        SAll=addstruct(SAll,S);
                             
                        ST.corr_seqType = sum(sum(C(7:12,1:6)))/(6*5);
                        ST.sessN=ss;
                        ST.roi=roi;
                        ST.sn=s;
                        ST.FoSEx=exe;
                        STAll=addstruct(STAll,ST);
                    end
                end
            end
        end
        save(fullfile(BGDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)),'-struct','SAll');
        save(fullfile(BGDir,sprintf('corrDist_seqType_RS_meanSubtract%d.mat',subtract_mean)),'-struct','STAll');
    case 'PLOT_corrDist_RS'
        roi=[1:3,5:7];
        subtract_mean=0;
        distType='crossval';
        vararginoptions(varargin,{'roi','subtract_mean','distType'});
        
        S=load(fullfile(BGDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        
        SeqType={'trained','untrained'};
        
        switch(distType)
            case 'crossval'
                var=S.corrDist;
            case 'non-crossval'
                var=S.corrDist_nonCrossval;
        end
        for seq=1:2
            indx=1;
            figure
            if seq==1
                style=styTrained_exe;
            else
                style=styUntrained_exe;
            end
            for r=roi
                subplot(1,numel(roi),indx)
                plt.line([S.sessN>3 S.sessN],var,'style',style,'split',S.FoSEx,'subset',S.seqType==seq&S.roi==r,'leg',{'1st','2nd'});
                plt.match('y');
                title(sprintf('%s',regname_BG{rem(r,4)}));
                if r==1
                    xlabel('Session'); ylabel(sprintf('Corr distance %s',SeqType{seq}));
                else
                    ylabel('');
                end
                indx=indx+1;
            end
        end   
    case 'PLOT_corrDist_RS_seqType'
        roi=[1:3,5:7];
        subtract_mean=0;
        vararginoptions(varargin,{'roi','subtract_mean'});
        
        S=load(fullfile(BGDir,sprintf('corrDist_seqType_RS_meanSubtract%d.mat',subtract_mean)));
        
        Exe={'1st exe','2nd exe'};
        
        figure;
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.line([S.sessN>3 S.sessN],S.corr_seqType,'style',stySeqType_exe,'split',S.FoSEx,'subset',S.roi==r,'leg',{'1st','2nd'});
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                xlabel('Session'); ylabel('Corr distance between seqType');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'CALC_corrDist_subtract'
       roi=[1:3,5:7];
       subtract_mean=0;
       vararginoptions(varargin,{'roi','subtract_mean'});
       
       S=load(fullfile(BGDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
 
       A = getrow(S,S.FoSEx==1);
       
       C.sn=A.sn;
       C.roi=A.roi;
       C.corrDist = [S.corrDist(S.FoSEx==1)-S.corrDist(S.FoSEx==2)];
       C.sessN=A.sessN;
       C.seqType=A.seqType;
       
       % % save
       save(fullfile(BGDir,'RepSup_subtract_corrDist.mat'),'-struct','C');
    case 'PLOT_corrDist_subtract_sess'
        sessN=[1:3];
        seqType='trained';
        regExcl=[2,6];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(BGDir,'RepSup_subtract_corrDist.mat')); 
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN));
        
        T.hemi=T.roi;
        T.hemi(T.roi<4,:)=1;
        T.hemi(T.roi>3,:)=2;
        
        switch(seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        figure
        for h=1:2
            subplot(1,2,h)
            plt.bar(T.roi,T.corrDist,'split',T.sessN,'subset',T.seqType==st&T.hemi==h,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
            drawline(0,'dir','horz');
            ylabel('Distance repetition suppression');
            set(gca,'XTickLabel',regname_BG(1:3));
            title(sprintf('%s %s',seqType,hemName{h}));
            xlabel('ROI');
        end
    case 'PLOT_corrDist_subtract_seqType'
        sessN=[1:3];
        regExcl=[2,6];
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(BGDir,'RepSup_subtract_corrDist.mat')); 
        T=getrow(T,~ismember(T.roi,regExcl)&ismember(T.sessN,sessN));
        
        T.hemi=T.roi;
        T.hemi(T.roi<4,:)=1;
        T.hemi(T.roi>3,:)=2;
        
        for ss=sessN
            figure
            for h=1
                plt.bar(T.roi,T.corrDist,'split',T.seqType,'subset',T.sessN==ss&T.hemi==h,'leg',{'trained','untrained'},'leglocation','northeast','style',stySeq);
                drawline(0,'dir','horz');
                ylabel('Distance repetition suppression');
                set(gca,'XTickLabel',regname_BG(1:3));
                title(sprintf('sess-%d %s',ss,hemName{h}));
                xlabel('ROI');
            end
        end
    case 'STATS_corrDist_subtract'
        betaChoice = 'multiPW';
        roi=[1:3,5:7];
        sessN=[1:4];
        seqLabel={'trained','untrained'};
        vararginoptions(varargin,{'betaChoice','roi','sessN'});
        
        T = load(fullfile(BGDir,sprintf('RepSup_subtract_%s.mat',betaChoice)));
        
        % main ANOVA - session x sequence
        for r=roi
            if r<4
                hem='LH';
            else
                hem='RH';
            end
            fprintf('\n PSC in %s %s \n',regname_BG{rem(r,4)},hem);
            anovaMixed(T.subtr_psc,T.sn,'within',[T.sessN T.seq],{'session','seqType'},'subset',T.roi==r);
        end
        
        % t-tests comparing sequence type per session
        for r=roi
            if r<4
                hem='LH';
            else
                hem='RH';
            end
            for ss=sessN    
                fprintf('\n post-hoc t-test on the effect of seqType on RS subtraction in sess %d in %s %s \n',ss,regname_BG{rem(r,4)},hem);
                ttestDirect(T.subtr_psc,[T.seq T.sn],2,'paired','subset',T.roi==r&T.sessN==ss);               
            end
        end
        
        % t-tests comparing across sessions - separately for T / UT seq
        for st=1:2
            fprintf('\n %s sequence\n',seqLabel{st});
            for sessTr=1:3
                for r=roi
                    if r<4
                        hem='LH';
                    else
                        hem='RH';
                    end
                    fprintf('\n post-hoc t-test on comparison of sessions %d-%d in %s %s\n',sessTr,sessTr+1,regname_BG{rem(r,4)},hem);
                    ttestDirect(T.subtr_psc,[T.sessN T.sn],2,'paired','subset',T.roi==r&ismember(T.sessN,[sessTr,sessTr+1])&T.seq==st);
                end
            end
        end
    case 'PLOT_corrDist_sess1'    
        roi=[1,3];
        subtract_mean=0;
        vararginoptions(varargin,{'roi','subtract_mean','distType'});
        
        SS=load(fullfile(BGDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        
        figure
        plt.bar(SS.roi,SS.corrDist,'split',SS.FoSEx,'subset',SS.sessN==1&ismember(SS.roi,roi),'style',stySeqType_exe,'leg',{'1st','2nd'});
        ylabel('Cosine distance');
        set(gca,'XTick',[2 5 8 11 14.5 18 21],'XTickLabel',regname_BG(roi));
        xlabel('ROI');
        
        for r=roi
        fprintf('%s \t',regname_BG{r});
            ttestDirect(SS.corrDist,[SS.FoSEx SS.sn],2,'paired','subset',SS.roi==r & SS.sessN==1);
        end
    case 'CALC_expectedDist_2nd'
        roi=[1,3];
        sn=[4:9,11:22];
        subtract_mean=0;
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'roi','subtract_mean','distType'});
        
        D = load(fullfile(BGDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));        
        P = load(fullfile(BGDir,sprintf('RepSup_percent_%s.mat',betaChoice)));
        
        keyboard;
        
        R=[];
        for s=sn
            for st=1:2
                for r=roi
                    Pp=getrow(P,P.sessN==1&P.roi==r&P.seq==st&P.sn==s);
                    Dd=getrow(D,D.sessN==1&D.roi==r&D.seqType==st&D.sn==s);
                    
                    dist1=Dd.corrDist(Dd.FoSEx==1);
                    dist2=Dd.corrDist(Dd.FoSEx==2);
                    expDist=Pp.psc_psc.*dist1;
                    
                    D2.dist=[dist1;expDist;dist2];
                    D2.indx=[1;2;3]; % first, expected 2nd, 2nd
                    D2.roi=[r;r;r];
                    D2.seqType=ones(size(D2.roi)).*st;
                    
                    
                    R=addstruct(R,D2);
                end
            end
        end
        
        figure
        plt.bar(R.roi,R.dist,'split',R.indx,'leg',{'1st distance','Expected 2nd distance','2nd distance'},...
            'leglocation','northeast','style',stySess);
        ylabel('Cosine distance');
        set(gca,'XTick',[2 6.5],'XTickLabel',regname_BG(roi));    
    case 'CALC_expectedDist_2nd_learning'
        roi=[1,3];
        sessN=[1:3];
        subtract_mean=0;
        betaChoice = 'multiPW';
        vararginoptions(varargin,{'roi','subtract_mean','distType','sessN'});
        
        D = load(fullfile(BGDir,sprintf('corrDist_RS_meanSubtract%d.mat',subtract_mean)));
        P = load(fullfile(BGDir,sprintf('RepSup_percent_%s.mat',betaChoice)));
        Dnew=[];
        for ss=sessN
            for r=roi
                for st=1:2
                    Pp=getrow(P,P.sessN==ss & P.roi==r & P.seq==st);
                    Dd=getrow(D,D.sessN==ss & D.roi==r & D.seqType==st);
                    
                    dist1=Dd.corrDist(Dd.FoSEx==1);
                    dist2=Dd.corrDist(Dd.FoSEx==2);
                    expDist=Pp.psc_psc.*dist1;
                    
                    D2.dist=[dist1;expDist;dist2];
                    D2.indx=[ones(size(expDist));ones(size(expDist))*2;ones(size(expDist))*3]; % first, expected 2nd, 2nd
                    D2.roi=[Dd.roi(Dd.FoSEx==1);Dd.roi(Dd.FoSEx==1);Dd.roi(Dd.FoSEx==1)];
                    D2.sn=[Dd.sn(Dd.FoSEx==1);Dd.sn(Dd.FoSEx==1);Dd.sn(Dd.FoSEx==1)];
                    D2.seqType=ones(size(D2.roi))*st;
                    D2.sessN=ones(size(D2.roi))*ss;
                    
                    Dnew=addstruct(Dnew,D2);
                end
            end
            
        end

        save(fullfile(BGDir,'expectedDist_RS.mat'),'-struct','Dnew');
    case 'PLOT_expectedDist'
        roi=[1,3];
        sessN=[1:3];
        vararginoptions(varargin,{'roi','sessN'});
        
        D=load(fullfile(BGDir,'expectedDist_RS.mat'));
      
        for ss=sessN
            figure(ss)
            sub=1;
            for r=roi
                subplot(1,numel(unique(roi)),sub)
                plt.bar(D.seqType,D.dist,'split',D.indx,'subset',D.roi==r &D.sessN==ss,'leg',{'1st distance','Expected 2nd distance','2nd distance'},...
                    'leglocation','northeast','style',stySess);
                xlabel('SeqType'); title(sprintf('%s',regname_BG{r}));
                if sub==1
                    ylabel('Distance');
                else
                    ylabel('');
                end
                sub=sub+1;
            end
        end
    
    case 'PCM_data_repsupModels' % new - repsup models
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1,3];
        sn=[4:9,11:22];
        sessN=1; % need to be two sessions at the time
        seqType='trained';
        models={'generic','specific','spec+genCorr'};
        M_type=1; % 1 - generic, 2 - specific, 3 - specific with gen corr
        
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
                    B=load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',ss)));
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
                    case 3
                        M=pcm_repsupModel_specific_genCorr;
                end
                T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm);
                %  C = pcm_correlation(Data,partVec,condVec,M{2},runEffect,M_type);
                if M_type==3
                    C = pcm_correlation(Data,partVec,condVec,M{5},runEffect,M_type);
                else
                    C = pcm_correlation(Data,partVec,condVec,M{4},runEffect,M_type);
                end
                T.roi = ones(size(T.SN))*r;
                AllReg=addstruct(AllReg,C);
                AllReg=addstruct(AllReg,T);
            end % region
            
            %remove some fields
            a1='reg'; a2='theta'; a3='theta_hat'; a4='thetaCr';
            AllReg=rmfield(AllReg,a1); AllReg=rmfield(AllReg,a2); AllReg=rmfield(AllReg,a3); AllReg=rmfield(AllReg,a4);
            % save output
            save(fullfile(BGDir,sprintf('PCM_repsup_reliability_%s_%s_sess%d_%s.mat',models{M_type},seqType,ss,runEffect)),'-struct','AllReg');
        end % session    
    case 'PLOT_pcm_corr_allSess'
        reg = [1,3];
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
                R=load(fullfile(BGDir,sprintf('PCM_repsup_reliability_%s_%s_sess%d_%s.mat',modelType,seqType{st},t,runEffect)));
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
                title(sprintf('%s',regname_BG{rem(r,4)}));
                if sub==1
                    ylabel(sprintf('%s corr',corrType{ct}));
                    xlabel('Session');
                else
                    ylabel('');
                end
                sub=sub+1;
            end     
        end
    case 'PLOT_pcm_RS_logBayes'
         reg = [1:3];
         sessN=[1:4];
         seqType='trained';
         runEffect='fixed';
         modelType='specific'; % or specific
         vararginoptions(varargin,{'reg','sessN','seqType','runEffect','modelType'});
         
         T=[];
         for t=sessN
            R=load(fullfile(BGDir,sprintf('PCM_repsup_reliability_%s_%s_sess%d_%s.mat',modelType,seqType,t,runEffect)));
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
    case 'CROSSVAL_RS_corr'
        % looking both into simple model (generic / specific)
        % and model with added pattern across sequences (shared)
        beta_choice = 'mw';
        reg = [1:3,5:7];
        sn=[4:9,11:22];
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
                    B=load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',ss)));
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
                    d1a=Data1(rem(partVec1,2)==1,:);
                    part1a=partVec1(rem(partVec1,2)==1,:);
                    d2a=Data2(rem(partVec1,2)==0,:);
                    part2a=partVec2(rem(partVec2,2)==0,:);
                    dA=[d1a;d2a];
                    partA=[part1a;part2a];
                    d1b=Data1(rem(partVec1,2)==0,:);
                    part1b=partVec1(rem(partVec1,2)==0,:);
                    d2b=Data2(rem(partVec1,2)==1,:);
                    part2b=partVec2(rem(partVec2,2)==1,:);
                    dB=[d1b;d2b];
                    partB=[part1b;part2b];
                    
                    
                    %pattern consistency
                    R1=rsa_patternConsistency(Data1,partVec1,condVec1);
                    R2=rsa_patternConsistency(Data2,partVec2,condVec2);
                    R1c=rsa_patternConsistency_crossval(Data1,partVec1,condVec1);
                    R2c=rsa_patternConsistency_crossval(Data2,partVec2,condVec2);
                    % across
                    RAacr=rsa_patternConsistency_crossval(dA,partA,condVec1);
                    RBacr=rsa_patternConsistency_crossval(dB,partB,condVec1);
                    Racr=mean([RAacr,RBacr]);
                    
                    C.sn=[p;p];
                    C.pConsist=[R1;R2];
                    C.pConsistCross=[R1c;R2c];
                    C.FoSEx=[1;2];
                    C.roi=[r;r];
                    if r<4
                        C.hemi=[1;1];
                    else
                        C.hemi=[2;2];
                    end
                    C.sessN=[ss;ss];
                    
                    R.consistFoSEx=Racr;
                    R.corr=Racr/sqrt(R1c*R2c);
                    R.sn=p;
                    R.roi=r;
                    R.sessN=ss;
                    
                    CC=addstruct(CC,C);
                    RR=addstruct(RR,R);
                end; % session
            end; % subj
        end

        keyboard;
        % save structure
        save(fullfile(BGDir,sprintf('RS_patternConsist_%s',seqType)),'-struct','CC');
        save(fullfile(BGDir,sprintf('RS_patternConsistCorr_%s',seqType)),'-struct','RR');
    case 'PLOT_crossval_RS_corr'
        reg=[1:3,5:7];
        sessN=[1:4];
        vararginoptions(varargin,{'sessN','reg'});
        
        seqType={'trained','untrained'};
        
        TT=[];
        for st=1:2
            T=load(fullfile(BGDir,sprintf('RS_patternConsistCorr_%s',seqType{st})));
            T.seqType=ones(size(T.sn))*st;
            TT=addstruct(TT,T);
        end
        
        indx=1;
        figure
        for r=reg
            subplot(1,numel(reg),indx)
            plt.line([TT.sessN>3 TT.sessN],TT.consistFoSEx,'subset',TT.roi==r,'split',TT.seqType,'leg',seqType,'style',stySeq,'leglocation','north');
            plt.match('y');
            title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Pattern consistency');
                xlabel('Session');
            else
                ylabel('');
                xlabel('');
            end
            indx=indx+1;
        end

            
    case 'BG_plot_activation' %---------------- PLOTTING OF ACTIVATION / 2D 3D etc. - work on this!!!
        sn=15; 
        roi=[1:3,5:7];
        sessN=[1:4];
        regType='native';
        seqType='trained';
        vararginoptions(varargin,{'sn','roi','sessN','regType','seqType'});
        
        for ss=sessN
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            for s=sn
                % load subject region def
                load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                Ss=getrow(S,S.SN==s);
                switch seqType
                    case 'trained'
                        seqData=Ss.psc_train;
                    case 'untrained'
                        seqData=Ss.psc_untrain;
                end
                figure
                sub=0;
                for r=roi
                    sub=sub+1;
                    if r<4
                        side='lh';
                        reg=sub;
                      %  data=[R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3)];
                        data=[R{r}.data(:,2),R{r}.data(:,3),R{r}.data(:,1)];
                    else
                        side='rh';
                        reg=sub-3;
                       % data=[R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3)];
                        data=[R{r}.data(:,2),R{r}.data(:,3),R{r}.data(:,1)];
                    end
                    subplot(2,3,sub)
                  %  scatter(R{r}.data(:,2),R{r}.data(:,3),10,Ss.psc_train{sub});
                  
                    plot3k(data,'ColorData',seqData{sub}','ColorRange',[-1 1],'Marker',{'o' 4}); axis equal;
                  
                    if  r>3
                        set(gca,'Zdir','reverse');
                    end
                    title(sprintf('PSC %s %s sess-%d %s',side,regname_BG{reg},ss,seqType));
                end
            end
        end
    case 'STRIATUM_plot_activation_3D'
        sn=[5:9,11:16]; 
        hemi=1:2;
        sessN=[1:4];
        regType='native';
        typeFig='func'; % anat, func, topVox - only anatomy or activation overall, only top voxels
        seqType={'trained','untrained'}; % trained or untrained
        seqIndx=[1,2]; % take both - or just trained / untrained
        numVox=100;
        vararginoptions(varargin,{'sn','roi','sessN','typeFig','numVox','seqIndx','regType'});
   
        for ss=sessN
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            
            A=load(fullfile(BGDir,sprintf('BG_activation_locus_top%dvoxels.mat',numVox)));
            for s=sn
                % load subject region def
                load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                Ss=getrow(S,S.SN==s);
                for st=seqIndx
                    
                    switch seqType{st}
                        case 'trained'
                            psc_seq=Ss.psc_train;
                        case 'untrained'
                            psc_seq=Ss.psc_untrain;
                    end
                    Aa=getrow(A,A.sn==s&A.sessN==ss);
                    figure(2*s-st)
                    for h=hemi
                        if h==1
                            data = [[R{1}.data(:,3); R{3}.data(:,3)] [R{1}.data(:,2); R{3}.data(:,2)] [R{1}.data(:,1); R{3}.data(:,1)]];
                            switch(typeFig)
                                case 'topVox'
                                    psc1=psc_seq{1}';
                                    psc2=psc_seq{3}';
                                    psc_temp1=zeros(size(psc1));
                                    psc_temp2=zeros(size(psc2));
                                    psc_temp1(Aa.voxIndx(Aa.roi==1,:))=psc1(Aa.voxIndx(Aa.roi==1,:));
                                    psc_temp2(Aa.voxIndx(Aa.roi==3,:))=psc2(Aa.voxIndx(Aa.roi==3,:));
                                    psc=[psc_temp1;psc_temp2];
                                case 'func'
                                    psc = [psc_seq{1}'; psc_seq{3}'];
                            end
                            sub=2*ss-1;
                        else
                            data = [[R{5}.data(:,3); R{7}.data(:,3)] [R{5}.data(:,2); R{7}.data(:,2)] [R{5}.data(:,1); R{7}.data(:,1)]];
                            switch(typeFig)
                                case 'topVox'
                                    psc1=psc_seq{4}';
                                    psc2=psc_seq{6}';
                                    psc_temp1=zeros(size(psc1));
                                    psc_temp2=zeros(size(psc2));
                                    psc_temp1(Aa.voxIndx(Aa.roi==5,:))=psc1(Aa.voxIndx(Aa.roi==5,:));
                                    psc_temp2(Aa.voxIndx(Aa.roi==7,:))=psc2(Aa.voxIndx(Aa.roi==7,:));
                                    psc=[psc_temp1;psc_temp2];
                                case 'func'
                                    psc = [psc_seq{4}'; psc_seq{6}'];
                            end
                            sub=2*ss;
                        end
                        subplot(4,2,sub)
                        switch typeFig
                            case 'func'
                                plot3k(data,'ColorData',psc,'ColorRange',[-1 1],'Marker',{'o' 4}); axis equal;
                            case 'topVox'
                                plot3k(data,'ColorData',psc,'ColorRange',[-1 1],'Marker',{'o' 4}); axis equal;
                            case 'anat'
                                plot3k(data,'Marker',{'o' 4});
                        end
                        if h==2
                            set(gca,'Zdir','reverse');
                        end
                        title(sprintf('%s caudate subj %d %s',hemName{h},s,seqType{st}));
                    end
                end
            end
        end
    case 'BG_plot_activation_2Dimagesc'
        sn=11;
        hemi=[1:2];
        sessN=[1:4];
        regType='native'; % or native
        pivot_operation='nanmean';
        imageType='psc'; % psc or Tmap;
        seqType='trained'; % trained / untrained
        thres=[0 2]; % threshold for plotting - [0 1] for psc
        vararginoptions(varargin,{'sn','roi','sessN','regType','hemi','pivot_operation','imageType','seqType','thres'});
        
        for ss=sessN
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            for s=sn
                figure(s)
                % load subject region def
                load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                Ss=getrow(S,S.SN==s);
                switch imageType
                    case 'psc'
                        if strcmp(seqType,'trained')
                            data_seq=Ss.psc_train;
                        elseif strcmp(seqType,'untrained')
                            data_seq=Ss.psc_untrain;
                        end
                    case 'Tmap'
                        if strcmp(seqType,'trained')
                            data_seq=Ss.Tmap_train;
                        elseif strcmp(seqType,'untrained')
                            data_seq=Ss.Tmap_untrain;
                        end
                end
                for h=hemi
                    
                    data_caud=[]; data_put=[]; data_pall=[];
                    if h==1
                        data_caud = [data_caud;data_seq{1}];
                        data_put = [data_put;data_seq{3}];
                        data_pall = [data_pall;data_seq{2}];
                        C.x=R{1}.data(:,1);
                        C.y=R{1}.data(:,2);
                        C.z=R{1}.data(:,3);
                        C.data=nanmean(data_caud,1)';
                        Pu.x=R{3}.data(:,1);
                        Pu.y=R{3}.data(:,2);
                        Pu.z=R{3}.data(:,3);
                        Pu.data=nanmean(data_put,1)';
                        Pa.x=R{2}.data(:,1);
                        Pa.y=R{2}.data(:,2);
                        Pa.z=R{2}.data(:,3);
                        Pa.data=nanmean(data_pall,1)';
                    else
                        data_caud = [data_caud;data_seq{4}];
                        data_put = [data_put;data_seq{6}];
                        data_pall = [data_pall;data_seq{5}];
                        C.x=R{5}.data(:,1);
                        C.y=R{5}.data(:,2);
                        C.z=R{5}.data(:,3);
                        C.data=nanmean(data_caud,1)';
                        Pu.x=R{7}.data(:,1);
                        Pu.y=R{7}.data(:,2);
                        Pu.z=R{7}.data(:,3);
                        Pu.data=nanmean(data_put,1)';
                        Pa.x=R{6}.data(:,1);
                        Pa.y=R{6}.data(:,2);
                        Pa.z=R{6}.data(:,3);
                        Pa.data=nanmean(data_pall,1)';
                    end
                    
                    % plot
                    orient_names={'sagittal','coronal','axial'};
                    for o=1:size(orient_names,2)
                        if o==1
                            c=pivottable(C.z,C.y,C.data,pivot_operation);
                            put=pivottable(Pu.z,Pu.y,Pu.data,pivot_operation);
                            pall=pivottable(Pa.z,Pa.y,Pa.data,pivot_operation);
                        elseif o==2
                            c=pivottable(C.x,C.y,C.data,pivot_operation);
                            put=pivottable(Pu.x,Pu.y,Pu.data,pivot_operation);
                            pall=pivottable(Pa.x,Pa.y,Pa.data,pivot_operation);
                        else
                            c=pivottable(C.x,C.z,C.data,pivot_operation);
                            put=pivottable(Pu.x,Pu.y,Pu.data,pivot_operation);
                            pall=pivottable(Pa.x,Pa.y,Pa.data,pivot_operation);
                        end
                        
                        figure((ss-1)*3+o)
                        subplot(3,2,h)
                        imagesc(flipud(c));
                        caxis(thres);
                        title(sprintf('%s caudate sess-%d %s %s %s',hemName{h},ss,orient_names{o},seqType,imageType));
                        subplot(3,2,2+h)
                        imagesc(flipud(put));
                        caxis(thres);
                        title(sprintf('%s putamen sess-%d %s %s %s %s',hemName{h},ss,regType,orient_names{o},seqType,imageType));
                        subplot(3,2,4+h)
                        imagesc(flipud(pall));
                        caxis(thres);
                        title(sprintf('%s palldum sess-%d %s %s %s %s',hemName{h},ss,regType,orient_names{o},seqType,imageType));
                end
                end
            end
        end
    
    case 'STRIATUM_plotAnat_allSubj_3D'
        sn = [2:9,11:20];
        hemi=1:2;
        regType='MNI_step2';  % MNI, MNI_step2, native, SPM
        vararginoptions(varargin,{'sn','roi','regType'});
        
        figure;
        for h=hemi
            data=[];
            for s=sn
                % load subject region def
                 load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                if h==1
                    data = [data; [R{1}.data(:,3); R{3}.data(:,3)] [R{1}.data(:,2); R{3}.data(:,2)] [R{1}.data(:,1); R{3}.data(:,1)]];
                else
                    data = [data; [R{5}.data(:,3); R{7}.data(:,3)] [R{5}.data(:,2); R{7}.data(:,2)] [R{5}.data(:,1); R{7}.data(:,1)]];
                end             
            end
            
            subplot(1,2,h);
            plot3k(data,'Marker',{'o' 4});
            if h==2
                set(gca,'Zdir','reverse');
            end
            title(sprintf('%s caudate',hemName{h}));
        end
    case 'BG_plotAnat_allSubj_2D'
        sn = [2:9,11:20];
        hemi=1:2;
        regType='SPM'; % MNI, MNI_step2, native, SPM
        plotAvrg=1; % overlay average shape
        avrgType='subjAvrg'; % which method - SPM or subjAvrg
        vararginoptions(varargin,{'sn','roi','figType','regType','plotAvrg','avrgType'});
        colMarker=hsv(length(sn));
        
        for h=hemi
            colIndx=1;
            for s=sn
                % load subject region def
                load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                %load(fullfile(regDir,sprintf('%s_BG_regions.mat',subj_name{s})));
                
                if h==1
                    caud = [R{1}.data(:,1) R{1}.data(:,2) R{1}.data(:,3)];
                    put = [R{3}.data(:,1)  R{3}.data(:,2) R{3}.data(:,3)];
                    pall = [R{2}.data(:,1)  R{2}.data(:,2) R{2}.data(:,3)];
                else
                    caud = [R{5}.data(:,1) R{5}.data(:,2) R{5}.data(:,3)];
                    put = [R{7}.data(:,1)  R{7}.data(:,2) R{7}.data(:,3)];
                    pall = [R{6}.data(:,1)  R{6}.data(:,2) R{6}.data(:,3)];
                end
                
                orient=[2,3;1,3;1,2];
                orient_names={'sagittal','coronal','axial'};
                for o=1:size(orient,1)       
                    figure(o)
                    subplot(3,2,h)
                    k=boundary(caud(:,orient(o,1)),caud(:,orient(o,2))); % determine coordinates
                    plot(caud(k,orient(o,1)),caud(k,orient(o,2)),'color',colMarker(colIndx,:));
                    hold on;
                    title(sprintf('%s caudate %s %s',hemName{h},regType,orient_names{o}));
                    subplot(3,2,2+h)
                    k=boundary(put(:,orient(o,1)),put(:,orient(o,2))); % Y and Z
                    plot(put(k,orient(o,1)),put(k,orient(o,2)),'color',colMarker(colIndx,:));
                    hold on;
                    title(sprintf('%s putamen %s %s',hemName{h},regType,orient_names{o}));
                    subplot(3,2,4+h)
                    k=boundary(pall(:,orient(o,1)),pall(:,orient(o,2))); % Y and Z
                    plot(pall(k,orient(o,1)),pall(k,orient(o,2)),'color',colMarker(colIndx,:));
                    hold on;
                    title(sprintf('%s pallidum %s %s',hemName{h},regType,orient_names{o}));                  
                end 
                colIndx=colIndx+1;
            end
            
            if plotAvrg % plot an average trace on top in black
                load(fullfile(regDir,sprintf('%s_MNI_BG_regions.mat',avrgType)));
                if h==1
                    caud = [R{1}.data(:,1) R{1}.data(:,2) R{1}.data(:,3)];
                    put = [R{3}.data(:,1)  R{3}.data(:,2) R{3}.data(:,3)];
                    pall = [R{2}.data(:,1)  R{2}.data(:,2) R{2}.data(:,3)];
                else
                    caud = [R{5}.data(:,1) R{5}.data(:,2) R{5}.data(:,3)];
                    put = [R{7}.data(:,1)  R{7}.data(:,2) R{7}.data(:,3)];
                    pall = [R{6}.data(:,1)  R{6}.data(:,2) R{6}.data(:,3)];
                end
                for o=1:size(orient,1)
                    figure(o)
                    subplot(3,2,h)
                    k=boundary(caud(:,orient(o,1)),caud(:,orient(o,2))); % determine coordinates
                    plot(caud(k,orient(o,1)),caud(k,orient(o,2)),'color',[0 0 0],'linewidth',2);
                    hold on;
                    subplot(3,2,2+h)
                    k=boundary(put(:,orient(o,1)),put(:,orient(o,2))); % Y and Z
                    plot(put(k,orient(o,1)),put(k,orient(o,2)),'color',[0 0 0],'linewidth',2);
                    hold on;  
                    subplot(3,2,4+h)
                    k=boundary(pall(:,orient(o,1)),pall(:,orient(o,2))); % Y and Z
                    plot(pall(k,orient(o,1)),pall(k,orient(o,2)),'color',[0 0 0],'linewidth',2);
                    hold on;  
                end
            end
            
        end
    case 'STRIATUM_plotFunc'
        sn = [5:9,11:16];
        hemi=1;
        figType='2D';
        sessN=1;
        seqType='trained';
        regType='SPM'; % MNI or MNI_step2 or native
        vararginoptions(varargin,{'sn','roi','figType','regType','sessN','seqType'});
        colMarker=hsv(length(sn));
        
        for ss=sessN
            %S=load(fullfile(regDir,sprintf('betas_BG_sess%d.mat',ss)));
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            for h=hemi
                colIndx=1;
                for s=sn
                    % load subject region def
                    load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                    %load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                    Ss=getrow(S,S.SN==s);
                    switch seqType
                        case 'trained'
                            psc_seq=Ss.psc_train;
                        case 'untrained'
                            psc_seq=Ss.psc_untrain;
                    end
         
                    if h==1
                        data = [[R{1}.data(:,3); R{3}.data(:,3)] [R{1}.data(:,2); R{3}.data(:,2)] [R{1}.data(:,1); R{3}.data(:,1)]];
                        caud = [R{1}.data(:,1) R{1}.data(:,2) R{1}.data(:,3)];
                        pall = [R{3}.data(:,1)  R{3}.data(:,2) R{3}.data(:,3)];
                    else
                        data = [[R{5}.data(:,3); R{7}.data(:,3)] [R{5}.data(:,2); R{7}.data(:,2)] [R{5}.data(:,1); R{7}.data(:,1)]];
                        caud = [R{5}.data(:,1) R{5}.data(:,2) R{5}.data(:,3)];
                        pall = [R{7}.data(:,1)  R{7}.data(:,2) R{7}.data(:,3)];
                    end
                    
                    switch(figType)
                        case '3D'
                            figure(1)
                            subplot(1,2,h);
                            plot3k(data,'Marker',{'o' 4},'MarkerColor',colMarker(colIndx,:));
                            hold on;
                            title(sprintf('%s striatum',hemName{h}));
                            if h==2
                                set(gca,'Zdir','reverse');
                            end
                        case '2D'
                            figure(1)
                            subplot(4,2,2*(ss-1)+1)
                            k=boundary(caud(:,2),caud(:,3)); % Y and Z
                            plot(caud(k,2),caud(k,3),'color',colMarker(colIndx,:));
                            hold on;
                            scatter(caud(:,2),caud(:,3),14,psc_seq{1}','filled');
                            colorbar; caxis([-1 1]);
                            title(sprintf('%s caudate-sess%d',hemName{h},ss));
                            subplot(4,2,2*ss)
                            k=boundary(pall(:,2),pall(:,3)); % Y and Z
                            plot(pall(k,2),pall(k,3),'color',colMarker(colIndx,:));
                            hold on;
                            scatter(pall(:,2),pall(:,3),14,psc_seq{3}','filled');
                            colorbar; caxis([-1 1]);
                            title(sprintf('%s pallidum-sess%d',hemName{h},ss));
                            
                            figure(2)
                            subplot(4,2,2*(ss-1)+1)
                            k=boundary(caud(:,1),caud(:,3)); % Y and Z
                            plot(caud(k,1),caud(k,3),'color',colMarker(colIndx,:));
                            hold on;
                            scatter(caud(:,1),caud(:,3),14,psc_seq{1}','filled');
                            colorbar; caxis([-1 1]);
                            title(sprintf('%s caudate-sess%d',hemName{h},ss));
                            subplot(4,2,2*ss)
                            k=boundary(pall(:,1),pall(:,3)); % Y and Z
                            plot(pall(k,1),pall(k,3),'color',colMarker(colIndx,:));
                            hold on;
                            scatter(pall(:,1),pall(:,3),14,psc_seq{3}','filled');
                            colorbar; caxis([-1 1]);
                            title(sprintf('%s pallidum-sess%d',hemName{h},ss));
                            
                            
                            figure(3)
                            subplot(4,2,2*(ss-1)+1)
                            k=boundary(caud(:,1),caud(:,2)); % Y and Z
                            plot(caud(k,1),caud(k,2),'color',colMarker(colIndx,:));
                            hold on;
                            scatter(caud(:,1),caud(:,2),14,psc_seq{1}','filled');
                            colorbar; caxis([-1 1]);
                            title(sprintf('%s caudate-sess%d',hemName{h},ss));
                            subplot(4,2,2*ss)
                            k=boundary(pall(:,1),pall(:,2)); % Y and Z
                            plot(pall(k,1),pall(k,2),'color',colMarker(colIndx,:));
                            hold on;
                            scatter(pall(:,1),pall(:,2),14,psc_seq{3}','filled');
                            colorbar; caxis([-1 1]);
                            title(sprintf('%s pallidum-sess%d',hemName{h},ss));
                    end
                    colIndx=colIndx+1;
                end
            end
        end
        
    case 'STRIATUM_plotFunc_allSubj'
        sn=[4:9,11:16]; 
        hemi=1:2;
        sessN=[1:4];
        typeFig='func'; % anat, func, topVox - only anatomy or activation overall, only top voxels
        regType='MNI_step2';
        seqType='trained'; % trained or untrained
        numVox=100;
        vararginoptions(varargin,{'sn','roi','sessN','typeFig','numVox','seqType'});
        
        for ss=sessN
            figure
            %S=load(fullfile(regDir,sprintf('betas_BG_sess%d.mat',ss)));
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            
           % A=load(fullfile(BGDir,sprintf('BG_activation_locus_top%dvoxels.mat',numVox)));
            A=load(fullfile(BGDir,sprintf('BG_activation_locus_%s_top%dvoxels.mat',regType,numVox)));
            for h=hemi
                data=[];
                psc=[];
                for s=sn
                    % load subject region def
                    %load(fullfile(regDir,sprintf('%s_BG_regions.mat',subj_name{s})));
                    load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                    Ss=getrow(S,S.SN==s);
                    switch seqType
                        case 'trained'
                            psc_seq=Ss.psc_train;
                        case 'untrained'
                            psc_seq=Ss.psc_untrain;
                    end
                    Aa=getrow(A,A.sn==s&A.sessN==ss);
                    if h==1
                        data = [data; [R{1}.data(:,3); R{3}.data(:,3)] [R{1}.data(:,2); R{3}.data(:,2)] [R{1}.data(:,1); R{3}.data(:,1)]];
                        switch(typeFig)
                            case 'topVox'
                                psc1=psc_seq{1}';
                                psc2=psc_seq{3}';
                                psc_temp1=zeros(size(psc1));
                                psc_temp2=zeros(size(psc2));
                                psc_temp1(Aa.voxIndx(Aa.roi==1,:))=psc1(Aa.voxIndx(Aa.roi==1,:));
                                psc_temp2(Aa.voxIndx(Aa.roi==3,:))=psc2(Aa.voxIndx(Aa.roi==3,:));
                                psc=[psc; psc_temp1;psc_temp2];
                            case 'func'
                                psc = [psc;psc_seq{1}'; psc_seq{3}'];
                        end

                    else
                        data = [data; [R{5}.data(:,3); R{7}.data(:,3)] [R{5}.data(:,2); R{7}.data(:,2)] [R{5}.data(:,1); R{7}.data(:,1)]];
                        switch(typeFig)
                            case 'topVox'
                                psc1=psc_seq{4}';
                                psc2=psc_seq{6}';
                                psc_temp1=zeros(size(psc1));
                                psc_temp2=zeros(size(psc2));
                                psc_temp1(Aa.voxIndx(Aa.roi==5,:))=psc1(Aa.voxIndx(Aa.roi==5,:));
                                psc_temp2(Aa.voxIndx(Aa.roi==7,:))=psc2(Aa.voxIndx(Aa.roi==7,:));
                                psc=[psc;psc_temp1;psc_temp2];
                            case 'func'
                                psc = [psc;psc_seq{4}'; psc_seq{6}'];
                        end
                    end
                    subplot(1,2,h)
                    switch typeFig
                        case 'func'
                           plot3k(data,'ColorData',psc,'ColorRange',[-1 1],'Marker',{'o' 4});
                            hold on
                        case 'topVox'
                            plot3k(data,'ColorData',psc,'ColorRange',[-1 1],'Marker',{'o' 4});
                            hold on
                        case 'anat'
                            plot3k(data,'Marker',{'o' 4});
                    end
                    if h==2
                         set(gca,'Zdir','reverse');
                    end
                   title(sprintf('%s caudate - sess%d',hemName{h},ss));
                end
            end
        end  
    case 'STRIATUM_plotFunc_group_3D'
        sn=[4:9,11:16];      
        hemi=1:2;
        sessN=[1:4];
        regType='subjAvrg';
        seqType='trained'; % trained or untrained
        vararginoptions(varargin,{'sn','roi','sessN','seqType'});
        
        for ss=sessN
            figure
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            load(fullfile(regDir,sprintf('%s_MNI_BG_regions.mat',regType)));
            for h=hemi
                s_indx=1;
                for s=sn
                    psc=[];
                    Ss=getrow(S,S.SN==s);
                    switch seqType
                        case 'trained'
                            psc_seq=Ss.psc_train;
                        case 'untrained'
                            psc_seq=Ss.psc_untrain;
                    end
                    if h==1
                        psc(s_indx,:) = [psc;psc_seq{1} psc_seq{3}];
                    else
                        psc(s_indx,:) = [psc;psc_seq{4} psc_seq{6}];
                    end
                    s_indx=s_indx+1;
                end
                subplot(1,2,h)
                if h==1
                    data = [[R{1}.data(:,3); R{3}.data(:,3)] [R{1}.data(:,2); R{3}.data(:,2)] [R{1}.data(:,1); R{3}.data(:,1)]];
                else
                    data = [[R{5}.data(:,3); R{7}.data(:,3)] [R{5}.data(:,2); R{7}.data(:,2)] [R{5}.data(:,1); R{7}.data(:,1)]];
                end       
                plot3k(data,'ColorData',nanmean(psc,1),'ColorRange',[-0.2 0.2],'Marker',{'o' 4});
                %plot3k(data,'ColorData',max(psc),'ColorRange',[-1 1],'Marker',{'o' 4});
                if h==2
                    set(gca,'Zdir','reverse');
                end
                title(sprintf('%s caudate - sess%d',hemName{h},ss));
            end
        end
    case 'STRIATUM_plotFunc_group_2Dscatter'
         sn=[4:9,11:16];      
        hemi=1:2;
        sessN=[1:4];
        regType='subjAvrg';
        seqType='trained'; % trained or untrained
        vararginoptions(varargin,{'sn','roi','sessN','seqType'});
        
        for ss=sessN
            figure
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            load(fullfile(regDir,sprintf('%s_MNI_BG_regions.mat',regType)));
            for h=hemi
                s_indx=1;
                psc_caud=[]; psc_pall=[];
                for s=sn
                    Ss=getrow(S,S.SN==s);
                    switch seqType
                        case 'trained'
                            psc_seq=Ss.psc_train;
                        case 'untrained'
                            psc_seq=Ss.psc_untrain;
                    end
                    if h==1
                        psc_caud = [psc_caud;psc_seq{1}];
                        psc_pall = [psc_pall;psc_seq{3}];
                    else
                        psc_caud = [psc_caud;psc_seq{4}];
                        psc_pall = [psc_pall;psc_seq{6}];
                    end
                    s_indx=s_indx+1;
                end
                subplot(1,2,h)
                if h==1
                    caud = [R{1}.data(:,1) R{1}.data(:,2) R{1}.data(:,3)];
                    pall = [R{3}.data(:,1) R{3}.data(:,2) R{3}.data(:,3)];
                else
                    caud = [R{5}.data(:,1) R{5}.data(:,2) R{5}.data(:,3)];
                    pall = [R{7}.data(:,1) R{7}.data(:,2) R{7}.data(:,3)];
                end
                % plot
                orient=[2,3;1,3;1,2];
                orient_names={'sagittal','coronal','axial'};
                for o=1:size(orient,1)                 
                    figure(o)
                    subplot(2,2,h)
                    k=boundary(caud(:,orient(o,1)),caud(:,orient(o,2))); % determine coordinates
                    plot(caud(k,orient(o,1)),caud(k,orient(o,2)),'color',[0 0 0]);
                    hold on;

                    scatter(caud(:,orient(o,1)),caud(:,orient(o,2)),50,nanmean(psc_caud,1)','filled');
                    title(sprintf('%s caudate sess-%d %s',hemName{h},ss,orient_names{o}));
                    subplot(2,2,2+h)
                    k=boundary(pall(:,orient(o,1)),pall(:,orient(o,2))); % Y and Z
                    plot(pall(k,orient(o,1)),pall(k,orient(o,2)),'color',[0 0 0]);
                    hold on;
                    scatter(pall(:,orient(o,1)),pall(:,orient(o,2)),70,nanmean(psc_pall,1),'filled');
                    title(sprintf('%s pallidum %s %s',hemName{h},regType,orient_names{o}));            
                end
                
            end
        end   
    case 'STRIATUM_plotFunc_group_2Dimagesc'
        sn=[4:9,11:18];      
        hemi=1:2;
        sessN=[1:4];
        regType='subjAvrg';
        seqType='trained'; % trained or untrained
        imageType='psc'; % psc or Tmap
        pivot_operation='nanmean'; % nanmean or nanmax
        thres=[0 2]; % threshold for plotting - [0 1] for psc
        vararginoptions(varargin,{'sn','roi','sessN','seqType','imageType','thres','pivot_operation'});
        
        for ss=sessN
            figure
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            load(fullfile(regDir,sprintf('%s_MNI_BG_regions.mat',regType)));
            for h=hemi
                s_indx=1;
                data_caud=[]; data_put=[]; data_pall=[];
                for s=sn
                    Ss=getrow(S,S.SN==s);
                    switch imageType
                        case 'psc'
                            if strcmp(seqType,'trained')
                                data_seq=Ss.psc_train;
                            elseif strcmp(seqType,'untrained')
                                data_seq=Ss.psc_untrain;
                            end
                        case 'Tmap'
                            if strcmp(seqType,'trained')
                                data_seq=Ss.Tmap_train;
                            elseif strcmp(seqType,'untrained')
                                data_seq=Ss.Tmap_untrain;
                            end
                    end
                    if h==1
                        data_caud = [data_caud;data_seq{1}];
                        data_put = [data_put;data_seq{3}];
                        data_pall = [data_pall;data_seq{2}];
                    else
                        data_caud = [data_caud;data_seq{4}];
                        data_put = [data_put;data_seq{6}];
                        data_pall = [data_pall;data_seq{5}];
                    end
                    s_indx=s_indx+1;
                end
                subplot(1,2,h)
                if h==1
                    C.x=R{1}.data(:,1);
                    C.y=R{1}.data(:,2);
                    C.z=R{1}.data(:,3);
                    C.data=nanmean(data_caud,1)';
                    Pu.x=R{3}.data(:,1);
                    Pu.y=R{3}.data(:,2);
                    Pu.z=R{3}.data(:,3);
                    Pu.data=nanmean(data_put,1)';
                    Pa.x=R{2}.data(:,1);
                    Pa.y=R{2}.data(:,2);
                    Pa.z=R{2}.data(:,3);
                    Pa.data=nanmean(data_pall,1)';
                else
                    C.x=R{5}.data(:,1);
                    C.y=R{5}.data(:,2);
                    C.z=R{5}.data(:,3);
                    C.data=nanmean(data_caud,1)';
                    Pu.x=R{7}.data(:,1);
                    Pu.y=R{7}.data(:,2);
                    Pu.z=R{7}.data(:,3);
                    Pu.data=nanmean(data_put,1)';
                    Pa.x=R{6}.data(:,1);
                    Pa.y=R{6}.data(:,2);
                    Pa.z=R{6}.data(:,3);
                    Pa.data=nanmean(data_pall,1)';
                end
                % plot
                orient_names={'sagittal','coronal','axial'};
                for o=1:size(orient_names,2)                 
                    if o==1
                        c=pivottable(C.z,C.y,C.data,pivot_operation);
                        put=pivottable(Pu.z,Pu.y,Pu.data,pivot_operation);
                        pall=pivottable(Pa.z,Pa.y,Pa.data,pivot_operation);
                    elseif o==2
                        c=pivottable(C.x,C.y,C.data,pivot_operation);
                        put=pivottable(Pu.x,Pu.y,Pu.data,pivot_operation);
                        pall=pivottable(Pa.x,Pa.y,Pa.data,pivot_operation);
                    else
                        c=pivottable(C.x,C.z,C.data,pivot_operation);
                        put=pivottable(Pu.x,Pu.y,Pu.data,pivot_operation);
                        pall=pivottable(Pa.x,Pa.y,Pa.data,pivot_operation);
                    end
                                        
                    figure((ss-1)*3+o)
                    subplot(3,2,h)
                    imagesc(flipud(c));
                    caxis(thres);
                    title(sprintf('%s caudate sess-%d %s %s %s',hemName{h},ss,orient_names{o},seqType,imageType));
                    subplot(3,2,2+h)
                    imagesc(flipud(put));
                    caxis([0 1]);
                    title(sprintf('%s putamen sess-%d %s %s %s %s',hemName{h},ss,regType,orient_names{o},seqType,imageType));      
                    subplot(3,2,4+h)
                    imagesc(flipud(pall));
                    caxis([0 1]);
                    title(sprintf('%s palldum sess-%d %s %s %s %s',hemName{h},ss,regType,orient_names{o},seqType,imageType)); 
                end
                
            end
        end
        
    case 'BG_activation_topVoxels'
        sn=[5:9,11:16]; 
        hemi=1:2;
        sessN=[1:4];
        roi=[1:3,5:7];
        numVox=100;
        thres=2.5;
        regType='MNI_step2';
        vararginoptions(varargin,{'sn','roi','sessN','hemi','numVox','regType'});
        
        AA=[];
        for ss=sessN
            %S=load(fullfile(regDir,sprintf('betas_BG_sess%d.mat',ss)));
            S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            for s=sn
                % load subject region def
                %load(fullfile(regDir,sprintf('%s_BG_regions.mat',subj_name{s})));
                load(fullfile(regDir,sprintf('%s_BG_%s_regions.mat',subj_name{s},regType)));
                for r=roi
                    for st=1:2 % sequence type
                        Ss=getrow(S,S.SN==s & S.region==r);
                        if st==1         
                            psc=Ss.psc_train{:};
                        else
                            psc=Ss.psc_untrain{:};
                        end
                        psc(:,isnan(psc))=0;
                        psc(:,psc>thres)=0;
                        [sortedX, sortedInds1]=sort(psc,'descend');
                        topVox = sortedInds1(1:numVox);
                        
                        A.psc=sortedX(1:numVox)';
                        A.voxIndx=topVox';
                        A.x_dim=R{r}.data(topVox,1);
                        A.y_dim=R{r}.data(topVox,2);
                        A.z_dim=R{r}.data(topVox,3);
                        A.roi=ones(numVox,1)*r;
                        A.sn=ones(numVox,1)*s;
                        A.sessN=ones(numVox,1)*ss;
                        A.seqType=ones(numVox,1)*st;
                        AA=addstruct(AA,A);
                    end
                end
            end
        end
        
        save(fullfile(BGDir,sprintf('BG_activation_locus_%s_top%dvoxels.mat',regType,numVox)),'-struct','AA');
    case 'BG_plot_activation_shift'
        numVox=100;
        roi=[1:3];
        seqType='trained'; % trained or untrained
        vararginoptions(varargin,{'numVox','roi','seqType'});
        T=load(fullfile(BGDir,sprintf('BG_activation_locus_top%dvoxels.mat',numVox)));
        
        switch (seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        for r=roi
            figure
            subplot(1,3,1)
            scatterplot(T.x_dim,T.psc,'split',T.sessN,'subset',T.roi==r & T.sessN<4 & T.seqType==st,'leg',{'sess1','sess2','sess3'});
            title(sprintf('%s - X',regname_BG{r}));
            subplot(1,3,2)
            scatterplot(T.y_dim,T.psc,'split',T.sessN,'subset',T.roi==r & T.sessN<4 & T.seqType==st,'leg',{'Sess1','sess2','sess3'});
            title(sprintf('%s - Y',regname_BG{r}));
            subplot(1,3,3)
            scatterplot(T.z_dim,T.psc,'split',T.sessN,'subset',T.roi==r & T.sessN<4 & T.seqType==st,'leg',{'Sess1','sess2','sess3'});
            title(sprintf('%s - Z',regname_BG{r}));
            figure
            scatterplot(T.y_dim,T.z_dim,'split',T.sessN,'subset',T.roi==r & T.sessN<4 & T.seqType==st,'leg',{'Sess1','sess2','sess3'});
            title(sprintf('%s - YZ coordinates',regname_BG{r}));
            figure
            scatterplot(T.x_dim,T.y_dim,'split',T.sessN,'subset',T.roi==r & T.sessN<4 & T.seqType==st,'leg',{'Sess1','sess2','sess3'});
            title(sprintf('%s - XY coordinates',regname_BG{r}));
            figure
            scatterplot(T.x_dim,T.z_dim,'split',T.sessN,'subset',T.roi==r & T.sessN<4 & T.seqType==st,'leg',{'Sess1','sess2','sess3'});
            title(sprintf('%s - XZ coordinates',regname_BG{r}));
        end
        
    case 'BG_plot_newSurface'
        
        regType='subjAvrg';
        ss=1; % session
        roi=1;
        vararginoptions(varargin,'roi');
        S=load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
        load(fullfile(regDir,sprintf('%s_MNI_BG_regions.mat',regType)));
        
        for r=roi
            x=R{r}.data(:,1);
            y=R{r}.data(:,2);
            z=R{r}.data(:,3);
            
        end
        
        keyboard;

        % Option 1
        tri=delaunay(x,y,z);
        p.Faces=tri;
        p.Vertices=[x, y, z];
        p.FaceVertexCData=z;
        patch(p)
        view(3)
        shading interp
        lighting p
        camlight
        
        % Option 2
        dt = delaunayTriangulation(x,y) ;
        tri = dt.ConnectivityList ;
        xi = dt.Points(:,1) ;
        yi = dt.Points(:,2) ;
        F = scatteredInterpolant(x,y,z);
        zi = F(xi,yi) ;
        trisurf(tri,xi,yi,zi)
        view(2)
        shading interp
        
        % option 3
        x_indx=unique(x); 
        y_indx=unique(y); 
        z_indx=unique(z); 
        
        %optional - interpolate 
        x_indx=interp(x_indx,3);
        y_indx=interp(y_indx,3);
        z_indx=interp(z_indx,3);
        
        
        
        X=size(x_indx,1);
        Y=size(y_indx,1);
        Z=size(z_indx,1);
        
        
        Vol=zeros(X,Y,Z);
        
        for i=1:Z
            x_val=x(z==z_indx(i));
            y_val=y(z==z_indx(i));
            
            for f=1:length(x_val)
                x_entry=find(ismember(x_indx,x_val(f)));
                y_entry=find(ismember(y_indx,y_val(f)));
                Vol(x_entry,y_entry,i)=1;
            end
        end
        
        Vq=interpn(Vol,'linear');
        xq=interpn(x_indx,'linear');
        yq=interpn(y_indx,'linear');
        zq=interpn(z_indx,'linear');
        
        
        PATCH_3Darray(Vol,x_indx',y_indx',z_indx');
        PATCH_3Darray(Vt,xt,yt,zt);
        %PATCH_3Darray(Vol,x_indx',y_indx',z_indx','cmap',[1 0.5 0.5]);
        view(3);
        
        
        % option 4
         V1=interp3(Vol,'linear');
         V2=smooth3(V1,'gaussian');
         figure
         vol3d('cdata',V2);
        view(3);
        
    case 'SEARCH_BG_define' % ------------------------ SEARCHLIGHT ---------------------
        sn=11;
        fig=1;
        BG_roi=[1,2,3];
        sphere_size=20; % size of searchlight sphere
        voxNum=100; % number of voxels to use
        vararginoptions(varargin,{'sn','fig','BG_roi','sphere_size'});
        
        for s=sn
            fprintf('Defining searchlight on %s\n',subj_name{s});
            SL = struct('LI',[],'li',[],'voxel',[],'voxmin',[],'voxmax',[]);
            % get functional mask
            funcMask    = rsa.readMask({fullfile(regDir,sprintf('mask_%s.nii',subj_name{s}))});
            funcMask.data = funcMask.mask;
            funcMask.fname = 'FuncMask';
            
            for h = 1:numel(hem)
                for bg = BG_roi
                    fprintf(' - Extracting: %d voxels from %s %s...\n',sphere_size,regname_BG{bg},hemName{h});
                    
                    % read BG mask file (.nii) - defined based on the func
                    % mask from overlap of glm (mask in regDir)
                    searchBGDir=fullfile(BGDir,'searchlight',subj_name{s});
                    dircheck(searchBGDir);
                    cd(searchBGDir);
                    BGMask = rsa.readMask({fullfile(BGDir,'regionMasks',...
                        sprintf('%s_%s_%s_native.nii',subj_name{s},regname_BG{bg},hem{h}))});

                    BGMask.data = BGMask.mask;
                    BGMask.fname = 'AnatMask';
                    % calc searchlight definition
                    L = rsa_defineSearchlight_volume(BGMask,funcMask,'sphere',[sphere_size voxNum]);
                    
                    % add searchlight
                    SL.LI       = [SL.LI;L.LI];
                    SL.voxel    = [SL.voxel;L.voxel];
                    SL.voxmin   = [SL.voxmin;L.voxmin];
                    SL.voxmax   = [SL.voxmax;L.voxmax];
                    SL.li{h,bg} = L.voxel;
                           
%                     % save searchlight
%                     SLname = fullfile(searchBGDir,sprintf('%s_%s_%s_vol_searchlight_%dvox.mat',subj_name{s},regname_BG{bg},hem{h},voxNum)); 
%                     save(SLname,'-struct','L');
                    
                    if fig==1
                        % visualize searchlight to check
                        figure('name',sprintf('Searchlight-%s_%s',regname_BG{bg},hemName{h}));
                        centers = surfing_inds2subs(funcMask.dim,L.voxel);
                        [coords(:,1),coords(:,2),coords(:,3)] = spmj_affine_transform(centers(:,1),centers(:,2),centers(:,3),funcMask.mat);
                        plot3(coords(:,1),coords(:,2),coords(:,3),'sk'); axis equal;
                        hold on;clear coords
                        centers = surfing_inds2subs(funcMask.dim,L.LI{1});
                        centers = double(centers);
                        [coords(:,1),coords(:,2),coords(:,3)] = spmj_affine_transform(centers(:,1),centers(:,2),centers(:,3),funcMask.mat);
                        plot3(coords(:,1),coords(:,2),coords(:,3),'.r'); axis equal;
                        clear coords
                        centers = surfing_inds2subs(funcMask.dim,L.LI{end});
                        centers = double(centers);
                        [coords(:,1),coords(:,2),coords(:,3)] = spmj_affine_transform(centers(:,1),centers(:,2),centers(:,3),funcMask.mat);
                        plot3(coords(:,1),coords(:,2),coords(:,3),'.b'); axis equal;
                        clear coords
                    end
                    
                    fprintf('...done.\n');
                end
            end

            % save searchlight
            SLname = fullfile(searchBGDir,sprintf('%s_BG_vol_searchlight_%dvox.mat',subj_name{s},voxNum));
            save(SLname,'-struct','SL');
            
        end
    case 'SEARCH_BG_runLDC'
        sessN=1;
        sn=11;
        BG_roi=[1,2,3];
        voxNum=100;
        spaceType='native'; % MNI or native
        vararginoptions(varargin,{'sn','sessN','spaceType'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        for s=sn
            for ss=sessN
                %   for h=1:numel(hem)
                %   for bg=BG_roi
                runs = 1:numruns_task_sess;
                % make index vectors
                conditionVec  = kron(ones(numel(runs),1),[1:12]');      % 12 sequences
                partition     = kron(runs',ones(12,1));
                % go to subject's glm directory
                cd(fullfile(glmSessDir{ss},subj_name{s}));
                % load their searchlight definitions and SPM file
                searchBGDir=fullfile(BGDir,'searchlight',subj_name{s});
                
                name = sprintf('%s_sess%d_BG',subj_name{s},ss);
                
                switch(spaceType)
                    case 'native'
                        L=load(fullfile(searchBGDir,sprintf('%s_BG_vol_searchlight_%dvox.mat',subj_name{s},voxNum)));
                        load SPM;
                        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                        % run the searchlight
                        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name);
                    case 'MNI'
                        L=load(fullfile(searchBGDir,sprintf('%s_BG_vol_searchlight_%dvox_MNI.mat',subj_name{s},voxNum)));
                        dirImageData;
                        % run the searchlight
                        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name,'imageDataDir',dirImageData);
                end
                %L = load(fullfile(searchBGDir,sprintf('%s_%s_%s_vol_searchlight_100vox.mat',subj_name{s},regname_BG{bg},hem{h})));

                
                %name = sprintf('%s_sess%d_%s_%s',subj_name{s},sessN,regname_BG{bg},hem{h});

                % end
                %   end % hem
            end % session
        end % subject
        cd(cwd);
    case 'SEARCH_BG_dist_map'
        sn  = 11;
        sessN = 1;
       % BG_roi=[1:3];
        vararginoptions(varargin,{'sn','sessN'});

        cWD = cd;
        for s = sn
            for ss=sessN
                %   for h=1:numel(hem)
                %  for bg=BG_roi
                % Load subject surface searchlight results (1 vol per paired conds)
                % LDC_file            = fullfile(glmSessDir{sessN},subj_name{s},sprintf('%s_sess%d_%s_%s_LDC.nii',subj_name{s},sessN,regname_BG{bg},hem{h})); % searchlight nifti
                LDC_file            = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_BG_LDC.nii',subj_name{s},ss)); % searchlight nifti
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
                Y.fname   = sprintf('%s_sess%d_BG_dist.nii',subj_name{s},ss);
                %Y.fname   = sprintf('%s_sess%d_%s_%s_dist.nii',subj_name{s},sessN,regname_BG{bg},hem{h});
                
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
                %T.fname = sprintf('%s_sess%d_%s_%s_dist_trained.nii',subj_name{s},sessN,regname_BG{bg},hem{h});
                %U.fname = sprintf('%s_sess%d_%s_%s_dist_untrained.nii',subj_name{s},sessN,regname_BG{bg},hem{h});
                %Z.fname = sprintf('%s_sess%d_%s_%s_dist_cross.nii',subj_name{s},sessN,regname_BG{bg},hem{h});
                T.fname = sprintf('%s_sess%d_BG_dist_trained.nii',subj_name{s},ss);
                U.fname = sprintf('%s_sess%d_BG_dist_untrained.nii',subj_name{s},ss);
                Z.fname = sprintf('%s_sess%d_BG_dist_cross.nii',subj_name{s},ss);
                
                % save outputs
                spm_write_vol(Y,Y.LDC);
                spm_write_vol(T,T.LDC);
                spm_write_vol(U,U.LDC);
                spm_write_vol(Z,Z.LDC);
                fprintf('Done %s sess%d BG\n',subj_name{s},ss);
                %fprintf('Done %s_sess%d %s-%s\n',subj_name{s},sessN,regname_BG{bg},hem{h});
                
                clear vol vdat LDC Y T U Z
                %  end; % region
                %  end; % hemi
            end; % session
        end; % subject
        cd(cWD);  % return to working directory
    case 'SEARCH_copy_native'
        sn  = [4:18];
        sessN = [1:4];
        fileName={'trained','untrained','cross'};
       % BG_roi=[1:3];
        vararginoptions(varargin,{'sn','sessN'});

        for s=sn
            for ss=sessN
                for f=1:size(fileName,2)
                    file=sprintf('%s_sess%d_BG_dist_%s.nii',subj_name{s},ss,fileName{f});
                    source=fullfile(glmSessDir{ss},subj_name{s},file);
                    dest=fullfile(BGDir,'searchlight',subj_name{s},file);
                    copyfile(source,dest);
                end
                fprintf('Done %s session %d \n',subj_name{s},ss);
            end
        end     
    case 'SEARCH_make_MNI'
        sessN=[1:4];
        sn=[4:9,11:18];
        fileName={'trained','untrained','cross'};
        vararginoptions(varargin,{'sn','sessN'});
        % align functional images - searchlight results to MNI-aligned anatomical
        
        for s=sn
            subjMNItrans=fullfile(anatomicalDir,subj_name{s},sprintf('y_%s_anatomical.nii',subj_name{s}));
            subjBGMNIDir=fullfile(BGDir,'searchlight','MNI',subj_name{s});
            dircheck(subjBGMNIDir);
            indx=1;
            for ss=sessN
                for f=1:length(fileName)
                    % Select images to be realigned - save them with a new name
                    image    = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_BG_dist_%s.nii',subj_name{s},ss,fileName{f}));
                    imageV=spm_vol(image);
                    V=spm_read_vols(imageV);
                    image_MNI=imageV;
                    image_MNI.fname=fullfile(subjBGMNIDir,sprintf('MNI_%s_sess%d_BG_dist_%s.nii',subj_name{s},ss,fileName{f}));
                    spm_write_vol(image_MNI,V);
                    
                    % Now select real the newly created one for realignment...
                    Q{indx}   = sprintf('%s,1',image_MNI.fname);
                    Q_del{indx} = sprintf('%s',image_MNI.fname); % whole file to delete once copied over
                    indx=indx+1;
                end; % file
            end; % session
            
            % definitions for realignment
            J.subj.def = {subjMNItrans};
            J.woptions.bb = [-78 -112 -70; 78 76 85];
            J.woptions.vox = [2 2 2];
            J.woptions.interp = 4;
            % run realignment
            for i=1:size(Q,2)
                J.subj.resample = {Q{i}}; % image to realign
                matlabbatch{1}.spm.spatial.normalise.write=J;
                spm_jobman('run',matlabbatch);
                % now delete the copied files (only wMNI_ is the transformed one)
                delete(Q_del{i});
            end
            
            fprintf('\n MNI transformation for %s searchlight done\n',subj_name{s});
        end % subject
    case 'SEARCH_AVRG'
        sessN=[1:4];
        sn=[4:9,11:18];
        fileName={'trained','untrained','cross'};
        vararginoptions(varargin,{'sn','sessN'})
        
        for ss=sessN
            for f=1:size(fileName,2)
                for s = 1:numel(sn)%1:numel(subj_name)
                    file=sprintf('wMNI_%s_sess%d_BG_dist_%s.nii',subj_name{sn(s)},ss,fileName{f});
                    IN{s} = fullfile(BGDir, 'searchlight','MNI',subj_name{sn(s)},file);
                end
                
                outDir = fullfile(BGDir,'searchlight','MNI','avrg');
                dircheck(outDir);
                OUT = fullfile(outDir,sprintf('avrg_sess%d_BG_dist_%s.nii',ss,fileName{f}));
                spmj_imcalc_mtx(IN,OUT,'nanmean(X)');
            end
        end        
    
    case 'SEARCH_BG_hist'
        sn=11;
        sessN=1;
        BG_roi=[1:3];
        vararginoptions(varargin,{'sn','sessN'});
        
        for s=sn
            for ss=sessN
                indx=1;
                for h=1:numel(hem)
                    for bg=BG_roi
                        Dist = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_%s_%s_dist.nii',subj_name{s},ss,regname_BG{bg},hem{h}));
                        DistVol     = spm_vol(Dist);
                        Dist_VDat    = spm_read_vols(DistVol); % is searchlight data
                        figure((ss-1)*3+1)
                        subplot(2,3,indx)
                        histogram(Dist_VDat); drawline(0,'dir','vert');
                        title(sprintf('%s %s Distance overall - sess%d',regname_BG{bg},hem{h},ss));
                        
                        DistTrain = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_%s_%s_dist_trained.nii',subj_name{s},ss,regname_BG{bg},hem{h}));
                        DistTrainVol = spm_vol(DistTrain);
                        DistTrain_VDat = spm_read_vols(DistTrainVol);
                        figure((ss-1)*3+2)
                        subplot(2,3,indx)
                        histogram(DistTrain_VDat); drawline(0,'dir','vert');
                        title(sprintf('%s %s Distance trained - sess%d',regname_BG{bg},hem{h},ss));
                        
                        DistUntrain = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_%s_%s_dist_untrained.nii',subj_name{s},ss,regname_BG{bg},hem{h}));
                        DistUntrainVol = spm_vol(DistUntrain);
                        DistUntrain_VDat = spm_read_vols(DistUntrainVol);
                        figure((ss-1)*3+3)
                        subplot(2,3,indx)
                        histogram(DistUntrain_VDat); drawline(0,'dir','vert');
                        title(sprintf('%s %s Distance untrained - sess%d',regname_BG{bg},hem{h},ss));
                        
                        indx=indx+1;
                    end
                end
            end
        end
    case 'SEARCH_BG_3D'
        sn=11;
        sessN=1;
        BG_roi=[1:3];
        seqType='trained';
        spaceType='native';
        vararginoptions(varargin,{'sn','sessN','seqType','spaceType'});
        
        for s=sn
            for h=1:numel(hem)
                for bg=BG_roi
                    searchBGDir=fullfile(BGDir,'searchlight',subj_name{s});
                    L = load(fullfile(searchBGDir,sprintf('%s_%s_%s_vol_searchlight_100vox.mat',subj_name{s},regname_BG{bg},hem{h})));
                    anatMask = rsa.readMask({fullfile(BGDir,'regionMasks',...
                        sprintf('%s_%s_%s_native.nii',subj_name{s},regname_BG{bg},hem{h}))});
                    funcMask    = rsa.readMask({fullfile(regDir,sprintf('mask_%s.nii',subj_name{s}))});
                    mask = anatMask.mask.*funcMask.mask;
                    
                    for ss=sessN
                        Dist = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_%s_%s_dist_%s.nii',subj_name{s},ss,regname_BG{bg},hem{h},seqType));
                        DistVol = spm_vol(Dist);
                        Dist_VDat = spm_read_vols(DistVol);
                        % make sure you are including only voxels in region
                        voxRegion=Dist_VDat.*mask;
                        % only non-nans - in the region
                        voxDist=voxRegion(~isnan(voxRegion));
                        
                        % calculate anatomical coordinates of search values
                        centers = surfing_inds2subs(funcMask.dim,L.voxel);
                        [S.x,S.y,S.z] = spmj_affine_transform(centers(:,1),centers(:,2),centers(:,3),funcMask.mat);
                        S.distValue = voxDist;
                        
                        % plot
                        figure(ss)
                        subplot(2,3,(h-1)*3+bg)
                        plot3k([S.y S.z S.x],'ColorData',S.distValue,'Marker',{'o' 4}); axis equal;
                        if h==2
                            set(gca,'Zdir','reverse');
                        end
                        title(sprintf('SEARCHLIGHT %s %s sess-%d %s',hem{h},regname_BG{bg},ss,seqType));
                    end
                    
                end
            end
        end

        
    otherwise
        disp('there is no such case.')
end % switch(what)

end

function dircheck(dir)
% Checks existance of specified directory. Makes it if it does not exist.

if ~exist(dir,'dir');
    %warning('%s didn''t exist, so this directory was created.\n',dir);
    mkdir(dir);
end
end

function C = calcDist(D,betaW,G)
% calculates G and distances for all sequences, trained / untrained /
% between the two
% INPUT: D - structure with run / cond etc.
%        betaW - all betas
%        G - overall G structure (12x12 - all seq)
% OUTPUT: Do - containes distances, eigenvalues

    % calculate trained / untrained G
    [G_train, Sig_train] = pcm_estGCrossval(betaW(D.seqType==1,:),D.run(D.seqType==1),D.seqNumb(D.seqType==1));
    [G_untrain, Sig_untrain] = pcm_estGCrossval(betaW(D.seqType==2,:),D.run(D.seqType==2),D.seqNumb(D.seqType==2));
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

    % calculate distances
    dist_train = rsa.rdm.squareRDM(diag(ind6*G_train*ind6'));
    dist_untrain = rsa.rdm.squareRDM(diag(ind6*G_untrain*ind6'));
    dist_cross = diag(ind2*G_seqType*ind2');
    dist_all = rsa.rdm.squareRDM(diag(ind12*G*ind12'));
    % check - new!!!!! Mar8th 2018
    dist_train = triu(dist_train); dist_train(dist_train==0)=NaN; 
    dist_untrain = triu(dist_untrain); dist_untrain(dist_untrain==0)=NaN; 
    dist_cross = triu(dist_cross); dist_cross(dist_cross==0)=NaN; 
    dist_all = triu(dist_all); dist_all(dist_all==0)=NaN; 
    
    
    C.dist_train = nanmean(dist_train(:));
    C.dist_untrain = nanmean(dist_untrain(:));
    C.dist_cross = nanmean(dist_cross);
    C.dist_all = nanmean(dist_all(:));
 
    % calculate eigenvalues
    H=eye(6)-ones(6,6)./6;  % centering matrix!
    G_trainCent = H*G_train*H';  % double centered G matrix - rows and columns
    C.eigTrain = sort(eig(G_trainCent)','descend');    % sorted eigenvalues - trainedSeq
    G_untrainCent = H*G_untrain*H';  % double centered G matrix - rows and columns
    C.eigUntrain = sort(eig(G_untrainCent)','descend'); % sorted eigenvalues - untrainedSeq
end
function M = pcm_repsupModel_generic

    % for sequence-specific modelling - one parameter for all seq           
    A=zeros(6);
    for i=1:6
        A(i,i)=1;
    end;
    
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'RepSup';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 4;
    M{3}.name       = 'RepSup+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    M{3}.Ac(:,3:8,3)  = [A;zeros(6)];      % Unique exe1 sequence patterns
    M{3}.Ac(:,9:14,4)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    
    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 5;
    M{4}.name         = 'RepSup+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    M{4}.Ac(:,3:8,3)  = [A;zeros(6)];      % Unique exe1 sequence patterns
    M{4}.Ac(:,9:14,4) = [zeros(6);A];     % Unique exe2 sequence pattterns
    M{4}.Ac(:,3:8,5)  = [zeros(6);A];     % Correlation exe1-exe2
    
end
function M = pcm_repsupModel_specific
    
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'RepSup';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 14;
    M{3}.name       = 'RepSup+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,3:8,2+i)  = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{3}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    end;

    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 20;
    M{4}.name         = 'RepSup+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,3:8,2+i)   = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{4}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
        M{4}.Ac(:,3:8,14+i)  = [zeros(6);A];     % Correlation exe1-exe2
    end;

end
function M = pcm_repsupModel_specific_genCorr
    
    % Model 1: No sequence pattern
    M{1}.type       = 'feature';
    M{1}.numGparams = 1;
    M{1}.name       = 'null';
    M{1}.Ac(:,1:12 ,1)  = zeros(12);
    
    % Model 2: First vs. second execution
    M{2}.type       = 'feature';
    M{2}.numGparams = 2;
    M{2}.name       = 'RepSup';
    M{2}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{2}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    
    % Model 3: Execution + sequence specific
    M{3}.type       = 'feature';
    M{3}.numGparams = 14;
    M{3}.name       = 'RepSup+Seq';
    M{3}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
    M{3}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{3}.Ac(:,3:8,2+i)  = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{3}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
    end;

    % Model 4: Execution + sequence specific + correlation in repetition
    M{4}.type         = 'feature';
    M{4}.numGparams   = 20;
    M{4}.name         = 'RepSup+Seq+Corr';
    M{4}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
    M{4}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
    % for sequence-specific modelling- one parameter per sequence
    for i=1:6
        A=zeros(6);
        A(i,i)=1;
        M{4}.Ac(:,3:8,2+i)   = [A;zeros(6)];      % Unique exe1 sequence patterns
        M{4}.Ac(:,9:14,8+i)  = [zeros(6);A];     % Unique exe2 sequence pattterns
        M{4}.Ac(:,3:8,14+i)  = [zeros(6);A];     % Correlation exe1-exe2
    end;
  
    % Model 5: Repetition suppression + generic correlation
    M{5}.type       = 'feature';
    M{5}.numGparams = 3;
    M{5}.name       = 'RepSup+GenCorr';
    M{5}.Ac(:,1,1) = [ones(6,1);zeros(6,1)];
    M{5}.Ac(:,2,2) = [zeros(6,1);ones(6,1)];
    M{5}.Ac(:,1,3) = [zeros(6,1);ones(6,1)];
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
% 1a. Empirical correlation
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec{p});
    b = pinv(Z)*Data{p};           % Estimate mean activities
    b(1:6,:)  = bsxfun(@minus,b(1:6,:) ,mean(b(1:6,:))); % Subtract mean per condition - first exe
    b(7:12,:) = bsxfun(@minus,b(7:12,:),mean(b(7:12,:))); % second exe
    G=cov(b');
    C.r_naive(p,1) = calcCorr(G);
end;

% 1b. Empirical correlation - subtract the whole mean pattern together
% for p=sn
%     Z=pcm_indicatorMatrix('identity',condVec{p});
%     b = pinv(Z)*Data{p};           % Estimate mean activities
%     b(1:12,:)  = bsxfun(@minus,b(1:12,:),mean(b(1:12,:))); % Subtract mean for all sequences
%     G=cov(b');
%     C.r_naive_meanSubtract(p,1) = calcCorr(G);
% end;

% --------------------------------------
% 2a. Crossvalidated correlation - make PD
for p=sn
    Z=pcm_indicatorMatrix('identity',condVec{p});
    % Subtract mean for each condition and run
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>6)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    C.r_crossval(p,1)=calcCorr(pcm_makePD(Gcv(:,:,p)));
   % C.r_crossval_noMakePD(p,1)=calcCorr(Gcv(:,:,p));
end;

% % 2b. Crossvalidated correlation - use nearestSPD
% for p=sn
%     Z=pcm_indicatorMatrix('identity',condVec{p});
%     % Subtract mean for each condition and run
%     X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>6)-1);
%     R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
%     Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
%     C.r_crossval_spd(p,1)=calcCorr(nearestSPD(Gcv(:,:,p)));
%     %C.r_crossval(p,1)=calcCorr(Gcv(:,:,p));
% end;
% --------------------------------------
% 3. Fit model 2  and infer correlations from the parameters
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);
C.theta=theta;
% Get the correlations
switch M_type
    case 1
        var1       = (theta{1}(3,:).^2)';
        var2       = (theta{1}(4,:).^2+theta{1}(5,:).^2)';
        cov12      = (theta{1}(3,:).*theta{1}(5,:))';
       % var1       = (theta{1}(1,:).^2)';
       % var2       = (theta{1}(2,:).^2+theta{1}(3,:).^2)';
       % cov12      = (theta{1}(1,:).*theta{1}(3,:))';
        C.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
    case 2
        var1       = (theta{1}(3:8,:).^2)';
        var2       = (theta{1}(9:14,:).^2+theta{1}(15:20,:).^2)';
        cov12      = (theta{1}(3:8,:).*theta{1}(15:20,:))';
        %var1       = (theta{1}(1:6,:).^2)';
        %var2       = (theta{1}(7:12,:).^2+theta{1}(13:18,:).^2)';
        %cov12      = (theta{1}(1:6,:).*theta{1}(13:18,:))';
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