function varargout=sml1_imana_BG_new(what,varargin)

% ------------------------- Directories -----------------------------------
%baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
behavDir        =[baseDir '/behavioral_data/data'];   
betaDir         =[baseDir '/betas'];
imagingDir      =[baseDir '/imaging_data'];                     
anatomicalDir   =[baseDir '/anatomicals'];       
caretDir        =[baseDir '/surfaceCaret'];              
regDir          =[baseDir '/RegionOfInterest/']; 
BGDir           =[baseDir '/basal_ganglia_new'];
BGanatDir       =[BGDir   '/anat'];
BGfuncDir       =[BGDir   '/func'];
BGstatsDir      =[BGDir   '/stats'];
pcmDir          =[baseDir '/pcm_stats'];
QCDir           =[baseDir '/quality_control'];

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
numruns           = [40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 20 40 10 10 10];
numruns_task      = 32;
numruns_loc       = 8;

sess = [repmat(1,1,10),repmat(2,1,10),repmat(3,1,10),repmat(4,1,10)];   % all sessions

sess_sn = [4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,4,1,1,1];    % per subject

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
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen'};
numregions_BG   = 3;
regSide=[ones(1,3) ones(1,3)*2]; % 1-left, 2-right
regType=[1:3  1:3]; 
regname_striatum = {'CaudateN','Putamen'};
regname_striatum_networks = {'Network1','Network2','Network3','Network4','Network5','Network6','Network7'};

% ------------------------- Freesurfer things -----------------------------         
atlasA    = 'x';                                                            % freesurfer filename prefix
atlasname = 'fsaverage_sym';                                                % freesurfer average atlas                                      % freesurfer hemisphere folder names    

% ------------------------- Subject things --------------------------------
% The variables in this section must be updated for every new subject.

subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
                's11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
                's21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'};  

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
stySeqType=style.custom({'gray'},'markersize',ms);
stySeqType_exe=style.custom({gray,lightgray},'markersize',ms);
styTrained_exe=style.custom({red,lightred},'markersize',ms);
styUntrained_exe=style.custom({blue,lightblue},'markersize',ms);
stySess=style.custom({gray,lightgray,silver,black},'markersize',ms);
% ------------------------------ Analysis Cases --------------------------------

switch(what)
    
    case '0_ANAT'                           % ---------------- ANATOMY (segment, MNI) -------------
    case 'ANAT_mni_FNIRT'                   % define MNI template via FLIRT - FNIRT
        sn=varargin{1};
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        % sml1_imana_BG_new('BG_MNI_FNIRT',4);
        atlasDir='standard/MNI152_T1_2mm.nii.gz';
        for s=sn
            subjMNIDir=fullfile(BGanatDir,'FNIRT',subj_name{s});
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
    case 'ANAT_mni_SPM'                     % define MNI template in SPM
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
    case 'ANAT_bg_segment'                  % do segmentation for BG - in native / FNIRT / SPM space
        % uses FSL to do segmentation for BG and ?map the regions 
        % option3 - space: native / FNIRT or SPM space
        % need to run MATLAB through terminal
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        
        % 1) Run step 1: sml1_imana_BG('ANAT_bg_segment',1,1,'native') / sml1_imana_prep('BG_FSLmask',1,1)
        % 2) Open the zip folder s01_BG_all_fast_firstseg.nii.gz
        % 3) Run step 2
        
        sn = varargin{1};
        step = varargin{2};
        space = varargin{3}; % native, FNIRT or SPM
        
        switch (step)
            case 1 % run FSL routine
                for s=sn%1:numel(subj_name)
                    switch(space)
                        case 'native'
                            IN= fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']);
                        case 'FNIRT'
                            IN= fullfile(BGanatDir,space,subj_name{s}, [subj_name{s}, '_MNI_FNIRT.nii']);
                        case 'SPM'
                            IN= fullfile(anatomicalDir, subj_name{s}, ['w', subj_name{s}, '_anatomical.nii']);
                    end
                    outDir = fullfile(BGanatDir, space);
                    dircheck(outDir);
                    OUT= fullfile(outDir, subj_name{s}, [subj_name{s},'_',space,'_BG.nii']);
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
                
                for s=sn%1:numel(subj_name)
                    %----deform info for basal ganglia ROI to individual space
                    
                    for h=1:2
                        for i=1:numregions_BG
                            fprintf('Working on subj: %i region: %s \n', s, [regname_BG{i},'_',hem{h}])
                            
                            subjBGDir=fullfile(BGanatDir,space,subj_name{s});
                            %----get names!
                            IN= fullfile(subjBGDir,[subj_name{s},'_',space,'_BG_all_fast_firstseg.nii']);
                            
                            if ~strcmp(space,'native')
                                OUT{i}= fullfile(subjBGDir,[subj_name{s},'_',regname_BG{i},'_',hem{h},'_',space,'.nii']);
                            else
                                OUT{i}= fullfile(subjBGDir,[subj_name{s},'_',regname_BG{i},'_',hem{h},'.nii']);
                            end
                            
                            %----make subj specific ROI image
                            spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                        end
                        %----do deformation
                        % spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end           
        end    
    case 'ANAT_mni_brain_subjAvrg'          % make average MNI template from subjects - whole brain
        sn=[5:9,11:31];
        template='SPM'; 
        vararginoptions(varargin,{'sn','type','template'});
        
        outDir = fullfile(BGanatDir, template, 'avrg');
        dircheck(outDir);
        for s = 1:numel(sn)
            IN{s} = fullfile(anatomicalDir,subj_name{sn(s)},...
                ['w',subj_name{sn(s)},'_anatomical.nii']);
        end
        OUT = fullfile(outDir,'subjAvrg_anat_MNI_brain.nii');
        spmj_imcalc_mtx(IN,OUT,'mean(X)');
    case 'ANAT_mni_BG_subjAvrg'             % make average MNI template from subjects - per BG roi
        sn=[5:9,11:17,19:25];
        template='SPM'; % SPM or FNIRT, striatum_buckner
        vararginoptions(varargin,{'sn','type','template'});
        
        outDir = fullfile(BGanatDir, template, 'avrg');
        dircheck(outDir);
        for h=1:2
            if strcmp(template,'SPM') | strcmp(template,'FNIRT') % per BG region 
                for i=1:numregions_BG
                    for s = 1:numel(sn)
                        IN{s} = fullfile(BGanatDir,template,subj_name{sn(s)},...
                            [subj_name{sn(s)},'_',regname_BG{i},'_',hem{h},'_',template,'.nii']);
                    end
                    OUT = fullfile(outDir,...
                        ['subjAvrg_anat_MNI_',regname_BG{i},'_',hem{h},'.nii']);
                    spmj_imcalc_mtx(IN,OUT,'mean(X)');
                end
            else % whole mask - funcStriatum, striatum_buckner
                for s = 1:numel(sn)
                    IN{s} = fullfile(BGanatDir,template,subj_name{sn(s)},['w' subj_name{sn(s)},'_',template,'_1mm.nii']);
                end
                OUT = fullfile(outDir,['subjAvrg_anat_',template,'.nii']);
                    spmj_imcalc_mtx(IN,OUT,'mean(X)');
            end
        end     
    case 'ANAT_MASK_make_subjAvrg_roi'      % make mask images from MNI (BG roi) templates (0 or 1 value)
        % make mask images per ROI, threshold initial templates based on
        % average size of ROI across subjects
        sn=[4:9,11:17,19:25];
        template='SPM'; % SPM or FNIRT 
        vararginoptions(varargin,{'sn','template'});

        for h=1:2
            for i=1:numregions_BG
                for s = 1:numel(sn)
                    % accummulate regions for all subjects  
                    subjFile = fullfile(BGanatDir,template,subj_name{sn(s)},...
                        [subj_name{sn(s)},'_',regname_BG{i},'_',hem{h}, '_',template,'.nii']);                   
                    Vol1=spm_vol(subjFile);
                    Vol2=spm_read_vols(Vol1);
                    numVox(s)=length(find(Vol2));
                end
                % calculate on average how many voxels per subject
                avrgVoxNum=round(mean(numVox));
                % load in group average 
                groupAvrg = fullfile(BGanatDir,template,'avrg',sprintf('subjAvrg_anat_MNI_%s_%s.nii',regname_BG{i},hem{h}));
                % destination anatomical mask
                dest= fullfile(BGanatDir,template,'avrg',sprintf('subjAvrg_anat_MNI_%s_%s_MASK.nii',regname_BG{i},hem{h}));
                makeGroupMask(avrgVoxNum,groupAvrg,dest);
                
                fprintf('Done: group anatMASK %s %s \n',regname_BG{i},hemName{h}); 
            end
        end
    case 'ANAT_MASK_combine_subjAvrg'       % create a mask for all of BG or striatum - combine previous masks
        template='SPM';
        maskType='BGall'; % whole BG or striatum
        maskNum = 2; % 2 - per hemisphere, 1 - overall
        vararginoptions(varargin,{'maskType','maskNum'});
        
        switch(maskType)
            case 'BGall'
                reg=[1:3]; % striatum, pallidum, putamen
            case 'striatum'
                reg=[1,3]; % striatum and putamen only
        end
        % makes a mask for all of the BG regions (both hemisphere)
        indx=1;
        for h=1:2
            for i=reg
                mask{indx}=fullfile(BGanatDir,template,'avrg',sprintf('subjAvrg_anat_MNI_%s_%s_MASK.nii',regname_BG{i},hem{h}));
                indx=indx+1;
            end
            if maskNum==2
                M{h}=combineMask(mask);
                indx=1; % start again
            end
        end
        
        M{maskNum}=combineMask(mask);
        % write the mask(s)
        for m = 1:maskNum
            BG_writeMask=spm_vol(mask{1});
            if maskNum==1
                BG_writeMask.fname=fullfile(BGanatDir,template,'avrg',sprintf('subjAvrg_anat_MNI_%s_MASK.nii',maskType));
            else
                BG_writeMask.fname=fullfile(BGanatDir,template,'avrg',sprintf('subjAvrg_anat_MNI_%s_%s_MASK.nii',maskType,hem{m}));
            end
            BG_writeMask.descip='avrg_wholeBG_mask';
            spm_write_vol(BG_writeMask,M{m});
        end
    case 'DICE_intersubject_overlap'        % calculate DICE coefficient for inter-subject overlap (SPM/FNIRT)
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
                BGsubjDir=fullfile(BGanatDir,sprintf('%s',template{t}));
                
                for c=1:size(indMat,1)
                    [i j]=find(indMat(c,:));
                    s1=i(2); s2=j(2);
                    subj1=sprintf(fullfile(BGsubjDir,subj_name{sn(s1)},sprintf('%s_%s_%s_%s.nii',subj_name{sn(s1)},regname_BG{r},hem{h},template{t})));
                    subj2=sprintf(fullfile(BGsubjDir,subj_name{sn(s2)},sprintf('%s_%s_%s_%s.nii',subj_name{sn(s2)},regname_BG{r},hem{h},template{t})));
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
    case 'DICE_group_overlap'               % calculate DICE coefficient for avrg-subject overlap (SPM/FNIRT)
        template={'SPM','FNIRT'}; 
        templateType=[1,2]; % 1, 2 or both
        sn=[4:9,11:17,19:22];
        roi=[1:3];
        hemi=[1:2];
        vararginoptions(varargin,{'sn','templateType','roi','hemi'});
        
        DD=[];
        for t=templateType
        for h=hemi
            for r=roi
                % group file
                BGgroupDir=fullfile(BGanatDir,template{t},'avrg');
                group=fullfile(BGgroupDir,sprintf('subjAvrg_anat_MNI_%s_%s_MASK.nii',regname_BG{r},hem{h}));
                groupV=spm_vol(group);
                groupRV=spm_read_vols(groupV);
                % all subjects
                indMat=indicatorMatrix('identity',sn);
                BGsubjDir=fullfile(BGanatDir,sprintf('%s',template{t}));     
                % subject file
                for c=1:size(indMat,1)
                    i=find(indMat(c,:));
                    s1=i; 
                    subj=sprintf(fullfile(BGsubjDir,subj_name{sn(s1)},sprintf('%s_%s_%s_%s.nii',subj_name{sn(s1)},regname_BG{r},hem{h},template{t})));
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
     
    case 'ANAT_thalamus_segment' % same for thalamus
        % uses FSL to do segmentation for BG and ?map the regions 
        % option3 - space: native / FNIRT or SPM space
        % need to run MATLAB through terminal
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        
        % 1) Run step 1: sml1_imana_BG('ANAT_bg_segment',1,1,'native') / sml1_imana_prep('BG_FSLmask',1,1)
        % 2) Open the zip folder s01_BG_all_fast_firstseg.nii.gz
        % 3) Run step 2
        
        sn = varargin{1};
        step = varargin{2};
        space = varargin{3}; % native only for now
        
        switch (step)
            case 1 % run FSL routine - this is already run for BG
                for s=sn%1:numel(subj_name)
                    
                    IN= fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']);
                    
                    outDir = fullfile(BGanatDir, space);
                    dircheck(outDir);
                    OUT= fullfile(outDir, subj_name{s}, [subj_name{s},'_',space,'_BG.nii']);
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
                BGnumber= [10; 49];
                %'CaudateN' 'Pallidum', 'Putamen' 'Thalamus'
                
                for s=sn%1:numel(subj_name)
                    %----deform info for basal ganglia ROI to individual space
                    
                    for h=1:2
                        fprintf('Working on subj: %i region: Thalamus-%s \n', s, hem{h});
                        
                        subjBGDir=fullfile(BGanatDir,space,subj_name{s});
                        %----get names!
                     %   IN= fullfile(subjBGDir,[subj_name{s},'_',space,'_BG_all_fast_firstseg.nii']);
                        IN= fullfile(subjBGDir,[subj_name{s},'_BG_all_fast_firstseg.nii']);
                        %----make subj specific ROI image
                        OUT{h}= fullfile(subjBGDir,[subj_name{s},'_Thalamus_',hem{h},'.nii']);
                        spm_imcalc_ui(IN,OUT{h},sprintf('i1==%d',BGnumber(h)));
                        %----do deformation
                        % spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end
        end
        
    case '0_SEARCH'                         % ------------------- SEARCHLIGHT  -----------------
    case 'SEARCH_striatum_define'           % define searchlights in striatum (caudate, putamen)
        sn=11;
        fig=1;
        BG_roi=[1,3];
        sphere_size=20; % size of searchlight sphere
        voxNum=100; % number of voxels to use
        vararginoptions(varargin,{'sn','fig','BG_roi','sphere_size'});
        
        for s=sn
            fprintf('Defining searchlight on %s\n',subj_name{s});
            SL = struct('LI',[],'li',[],'voxel',[],'voxmin',[],'voxmax',[]);
            % get functional mask - from glms of all sessions  
            funcMask    = rsa.readMask({fullfile(BGfuncDir,'native',subj_name{s},sprintf('maskGLM_%s.nii',subj_name{s}))});
            funcMask.data = funcMask.mask;
            funcMask.fname = 'FuncMask';
            
            for h = 1:numel(hem)
                for bg = BG_roi
                    fprintf(' - Extracting: %d voxels from %s %s...\n',sphere_size,regname_BG{bg},hemName{h});
                    
                    % read BG mask file (.nii) - defined based on the func
                    % mask from overlap of glm (mask in regDir)
                    searchBGDir=fullfile(BGfuncDir,'native',subj_name{s});
                    dircheck(searchBGDir);
                    cd(searchBGDir);
                  %  anatMask = rsa.readMask({fullfile(BGanatDir,'native',subj_name{s},...
                   %     sprintf('%s_%s_%s.nii',subj_name{s},regname_BG{bg},hem{h}))});
                    anatMask = rsa.readMask({fullfile(BGfuncDir,'native',subj_name{s},...
                        sprintf('funcMask_%s_%s_%s.nii',subj_name{s},regname_BG{bg},hem{h}))});
                   
                    anatMask.data = anatMask.mask;
                    anatMask.fname = 'AnatMask';
                    % calc searchlight definition
                    L = rsa_defineSearchlight_volume(anatMask,funcMask,'sphere',[sphere_size voxNum]);
                    
                    % add searchlight
                    SL.LI       = [SL.LI;L.LI];
                    SL.voxel    = [SL.voxel;L.voxel];
                    SL.voxmin   = [SL.voxmin;L.voxmin];
                    SL.voxmax   = [SL.voxmax;L.voxmax];
                    SL.li{h,bg} = L.voxel;
                                      
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
            SLname = fullfile(searchBGDir,sprintf('%s_striatum_vol_searchlight_%dvox.mat',subj_name{s},voxNum));
            save(SLname,'-struct','SL');
            
        end
    case 'SEARCH_striatum_runLDC'           % run searchlights - overall distances
        sessN=1;
        sn=11;
        voxNum=100;
        vararginoptions(varargin,{'sn','sessN'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        for s=sn
            for ss=sessN
                runs = 1:numruns_task_sess;
                % make index vectors
                conditionVec  = kron(ones(numel(runs),1),[1:12]');      % 12 sequences
                partition     = kron(runs',ones(12,1));
                % go to subject's glm directory
                cd(fullfile(glmSessDir{ss},subj_name{s}));
                % load their searchlight definitions and SPM file
                searchBGDir=fullfile(BGfuncDir,'native',subj_name{s});
                
                name = sprintf('%s_sess%d_striatum',subj_name{s},ss);
                
                L=load(fullfile(searchBGDir,sprintf('%s_striatum_vol_searchlight_%dvox.mat',subj_name{s},voxNum)));
                load SPM;
                SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                % run the searchlight
                rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name);
                
            end % session
        end % subject
        cd(cwd);
    case 'SEARCH_striatum_runFoSEx'         % run searchlights - 1st / 2nd exe - Mahalanobis dist
        sessN=1;
        sn=11;
        voxNum=100;
        repInd=[1,2]; % do searchlight for first or second execution
        vararginoptions(varargin,{'sn','sessN','spaceType','repInd'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        for s=sn
            for ss=sessN
                % go to subject's glm directory
                cd(fullfile(glmFoSExDir{ss},subj_name{s}));
                for exe=repInd
                    % check first if the searchlight map already exists
                    if exist(sprintf('%s_sess%d_exe%d_striatum_LDC.nii',subj_name{s},ss,exe));
                        fprintf('Searchlight map: %s sess%d exe%d already exists - skipping\n',subj_name{s},ss,exe);
                    else
                        D=load('SPM_info.mat');
                        
                        % initialise vectors
                        runs = 1:numruns_task_sess;
                        conditionVec  = zeros(size(D.run));    % 12 sequences
                        partition     = zeros(size(D.run));
                        
                        % only fill the ones of the right repetition
                        conditionVec(D.FoSEx==exe,:)  = D.seqNumb(D.FoSEx==exe);
                        partition(D.FoSEx==exe,:)     = D.run(D.FoSEx==exe);
                        
                        % load their searchlight definitions and SPM file
                        searchBGDir=fullfile(BGfuncDir,'native',subj_name{s});
                        L=load(fullfile(searchBGDir,sprintf('%s_striatum_vol_searchlight_%dvox.mat',subj_name{s},voxNum)));
                        
                        name = sprintf('%s_sess%d_exe%d_striatum',subj_name{s},ss,exe);
                        
                        load SPM;
                        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                        % run the searchlight
                        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name);
                        fprintf('Searchlight done: %s sess%d exe%d',subj_name{s},ss,exe);
                    end
                end; % repInd
            end % session
        end % subject
        cd(cwd);
    case 'SEARCH_striatum_runCorrDist'      % run searchlights - correlation distance
        % run correlation distance as searchlight
        sn      = [4:9,11:25];
        sessN   = [4];
        voxNum  = 100;
        vararginoptions(varargin,{'sn','sessN'});
        
        block   = 5e7;
        cwd     = pwd;                                                        % copy current directory (to return to later)
        runs    = 1:numruns_task_sess;
        nConds  = length(num_seq);
        % make index vectors
        conditionVec  = kron(ones(numel(runs),1),[1:nConds]');      % 12 sequences
        partitionVec  = kron(runs',ones(nConds,1));
        for s=sn
            for ss=sessN
                % go to subject's glm directory
                cd(fullfile(glmSessDir{ss},subj_name{s}));
                
                if exist(sprintf('%s_sess%d_striatum_corrDist.nii',subj_name{s},ss));
                    fprintf('Searchlight map: %s sess%d already exists - skipping\n',subj_name{s},ss);
                else
                    % load their searchlight definitions and SPM file
                    searchBGDir = fullfile(BGfuncDir,'native',subj_name{s});
                    L = load(fullfile(searchBGDir,sprintf('%s_striatum_vol_searchlight_%dvox.mat',subj_name{s},voxNum)));
                    load SPM;
                    SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                    inFiles = SPM.xY.VY;
                    % define the outfiles
                    nfiles  = nConds * (nConds-1)/2;
                    name = sprintf('%s_sess%d_striatum_corrDist',subj_name{s},ss);
                    for i = 1:nfiles
                        outFiles{i} = fullfile(cd,sprintf('%s.nii,%d',name,i));
                    end;
                    % run the searchlight
                    rsa.runSearchlight(L,inFiles,outFiles,@calcCorrSearch,'optionalParams',{SPM,conditionVec,partitionVec},'idealBlock',block);
                    fprintf('Searchlight map corrDist: %s sess%d done\n\n\n',subj_name{s},ss);
                end 
            end; % sessN
        end; % sn
        cd(cwd);
    case 'SEARCH_striatum_runFoSExCorrDist' % run searhclights - 1st / 2nd exe - corrDist
        sessN=4;
        sn=[4:9,11:25];
        voxNum=100;
        repInd=[1:2]; % do searchlight for first or second execution
        vararginoptions(varargin,{'sn','sessN','repInd'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        for s=sn
            for ss=sessN
                % go to subject's glm directory
                cd(fullfile(glmFoSExDir{ss},subj_name{s}));
                for exe=repInd
                    % check first if the searchlight map already exists
                    if exist(sprintf('%s_sess%d_exe%d_striatum_corrDist.nii',subj_name{s},ss,exe));
                        fprintf('Searchlight map: %s sess%d exe%d already exists - skipping\n',subj_name{s},ss,exe);
                    else
                        D=load('SPM_info.mat');
                        % initialise vectors
                        conditionVec  = zeros(size(D.run));    % 12 sequences
                        partitionVec  = zeros(size(D.run));
                        nConds        = length(num_seq);
                        % only fill the ones of the right repetition
                        conditionVec(D.FoSEx==exe,:)  = D.seqNumb(D.FoSEx==exe);
                        partitionVec(D.FoSEx==exe,:)  = D.run(D.FoSEx==exe);
                        
                        % load their searchlight definitions and SPM file
                        searchBGDir = fullfile(BGfuncDir,'native',subj_name{s});
                        L = load(fullfile(searchBGDir,sprintf('%s_striatum_vol_searchlight_%dvox.mat',subj_name{s},voxNum)));
                        
                        load SPM;
                        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));
                        inFiles = SPM.xY.VY;
                        % define the outfiles
                        nfiles  = nConds * (nConds-1)/2;
                        name = sprintf('%s_sess%d_exe%d_striatum_corrDist',subj_name{s},ss,exe);
                        for i = 1:nfiles
                            outFiles{i} = fullfile(cd,sprintf('%s.nii,%d',name,i));
                        end;
                        % run the searchlight
                        rsa.runSearchlight(L,inFiles,outFiles,@calcCorrSearch,'optionalParams',{SPM,conditionVec,partitionVec},'idealBlock',block);
                        fprintf('Searchlight done corrDist: %s sess%d exe%d\n\n',subj_name{s},ss,exe);
                    end
                end; %repInd
            end; % sessionsml1
        end; % subject
        cd(cwd);
    case 'SEARCH_striatum_dist_map'       % make a nifti map - trained / untrained / cross dist
        sn       = 4;
        sessN    = 4;
        distType = 'LDC'; %LDC or corrDist
        vararginoptions(varargin,{'sn','sessN','distType'});

        cWD = cd;
        for s = sn
            for ss=sessN
                % Load subject surface searchlight results (1 vol per paired conds)
                LDC_file            = fullfile(glmSessDir{ss},subj_name{s},sprintf('%s_sess%d_striatum_%s.nii',subj_name{s},ss,distType)); % searchlight nifti
                [subjDir,fname,ext] = fileparts(LDC_file);
                % copy the file to BG searchlight directory
                destDir=fullfile(BGfuncDir,'native',subj_name{s});
                copyfile(LDC_file,destDir);
                cd(destDir);
                
                vol     = spm_vol([fname ext]);
                vdat    = spm_read_vols(vol); % is searchlight data
                
                % average across all paired dists
                Y.LDC   = nanmean(vdat,4);
                % prep output file
                Y.dim   = vol(1).dim;
                Y.dt    = vol(1).dt;
                Y.mat   = vol(1).mat;
                switch distType
                    case 'LDC'
                        Y.fname   = sprintf('%s_sess%d_striatum_dist.nii',subj_name{s},ss);
                    case 'corrDist'
                        Y.fname   = sprintf('%s_sess%d_striatum_corrDist_mean.nii',subj_name{s},ss);
                end
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
    
                switch distType
                    case 'LDC'
                        T.fname = sprintf('%s_sess%d_striatum_dist_trained.nii',subj_name{s},ss);
                        U.fname = sprintf('%s_sess%d_striatum_dist_untrained.nii',subj_name{s},ss);
                        Z.fname = sprintf('%s_sess%d_striatum_dist_cross.nii',subj_name{s},ss);
                    case 'corrDist'
                        T.fname = sprintf('%s_sess%d_striatum_corrDist_trained.nii',subj_name{s},ss);
                        U.fname = sprintf('%s_sess%d_striatum_corrDist_untrained.nii',subj_name{s},ss);
                        Z.fname = sprintf('%s_sess%d_striatum_corrDist_cross.nii',subj_name{s},ss);
                end
                % save outputs
                spm_write_vol(Y,Y.LDC);
                spm_write_vol(T,T.LDC);
                spm_write_vol(U,U.LDC);
                spm_write_vol(Z,Z.LDC);
                fprintf('Done striatum searchlight maps %s: %s sess%d\n',distType,subj_name{s},ss);
                
                clear vol vdat LDC Y T U Z
            end; % session
        end; % subject
        cd(cWD);  % return to working directory    
    case 'SEARCH_striatum_FoSEx_map'
        sn  = 11;
        sessN = 1;
        repInd=1;
        distType = 'LDC'; %LDC or corrDist
        vararginoptions(varargin,{'sn','sessN','repInd','distType'});
        
        cWD = cd;
        for r=repInd
            for s = sn
                for ss=sessN
                    % Load subject surface searchlight results (1 vol per paired conds)
                    LDC_file            = fullfile(glmFoSExDir{ss},subj_name{s},sprintf('%s_sess%d_exe%d_striatum_%s.nii',subj_name{s},ss,r,distType)); % searchlight nifti
                    [subjDir,fname,ext] = fileparts(LDC_file);
                    % copy the file to BG searchlight directory
                    destDir=fullfile(BGfuncDir,'native',subj_name{s});
                    copyfile(LDC_file,destDir);
                    cd(destDir);
                    
                    vol     = spm_vol([fname ext]);
                    vdat    = spm_read_vols(vol); % is searchlight data
                    
                    % average across all paired dists
                    Y.LDC   = nanmean(vdat,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat;
                    switch distType
                        case 'LDC'
                            Y.fname   = sprintf('%s_sess%d_exe%d_striatum_dist.nii',subj_name{s},ss,r);
                        case 'corrDist'
                            Y.fname   = sprintf('%s_sess%d_exe%d_striatum_corrDist_mean.nii',subj_name{s},ss,r);
                    end
                    
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
                    switch distType
                        case 'LDC'
                            T.fname = sprintf('%s_sess%d_striatum_dist_trained_exe%d.nii',subj_name{s},ss,r);
                            U.fname = sprintf('%s_sess%d_striatum_dist_untrained_exe%d.nii',subj_name{s},ss,r);
                            Z.fname = sprintf('%s_sess%d_striatum_dist_cross_exe%d.nii',subj_name{s},ss,r);
                        case 'corrDist'
                            T.fname = sprintf('%s_sess%d_striatum_corrDist_trained_exe%d.nii',subj_name{s},ss,r);
                            U.fname = sprintf('%s_sess%d_striatum_corrDist_untrained_exe%d.nii',subj_name{s},ss,r);
                            Z.fname = sprintf('%s_sess%d_striatum_corrDist_cross_exe%d.nii',subj_name{s},ss,r);
                    end
                    % save outputs
                    spm_write_vol(Y,Y.LDC);
                    spm_write_vol(T,T.LDC);
                    spm_write_vol(U,U.LDC);
                    spm_write_vol(Z,Z.LDC);
                    fprintf('Done striatum searchlight maps %s exe%d: %s sess%d\n',distType,r,subj_name{s},ss);
                    
                    clear vol vdat LDC Y T U Z
                end; % session
            end; % subject
        end; % repInd
        cd(cWD);  % return to working directory
    case 'SEARCH_striatum_MNI'
        % make MNI group average
        sn=[4:9,11:25];
        sessN=[1:4];
        seqType={'trained','untrained','cross'}; 
        distType = 'dist'; % dist for LDC, corrDist otherwise
        vararginoptions(varargin,{'sn','sessN','seqType','distType'});

        sml1_imana_BG_new('FUNC_copyFile_nativeToMNI','sn',sn,'sessN',sessN,'seqType',seqType,'imageType','searchlight','distType',distType);
    case 'SEARCH_striatum_MNI_FoSEx'
        % make MNI group average
        sn=[4:9,11:25];
        sessN=[1:4];
        distType='corrDist'; % dist for LDC, corrDist otherwise
        seqType={'trained_exe1','trained_exe2','untrained_exe1','untrained_exe2','cross_exe1','cross_exe2'}; 
        vararginoptions(varargin,{'sn','sessN','seqType'});
        sml1_imana_BG_new('FUNC_copyFile_nativeToMNI','sn',sn,'sessN',sessN,'seqType',seqType,'imageType','searchlight','distType',distType);
       
    case 'FUNC_copyFile_nativeToMNI'    % copy original files to BGfuncDir, make MNI transform (psc, spm, mask, searchlight)
        sn  = [4:9,11:22];
        sessN = [1:4];
        seqType={'TrainSeq','UntrainSeq'}; % {'AllSeq','TrainSeq','UntrainSeq'}
       % regType={'Putamen','Pallidum','CaudateN'};
        regType={'Putamen','CaudateN'};
        networks=[1:7];
        imageType='spmT';
        glm='glmSessDir';
        distType = 'dist'; % dist(LDC) or corrDist (used for searchlight)
      
        vararginoptions(varargin,{'sn','sessN','imageType','sessN','seqType','glm','distType'});
        
        for s=sn
            % subject specific MNI transformation
            subjMNItrans=fullfile(anatomicalDir,subj_name{s},sprintf('y_%s_anatomical.nii',subj_name{s}));
            outDirNative=fullfile(BGfuncDir,'native',subj_name{s});
            outDirMNI=fullfile(BGfuncDir,'MNI',subj_name{s});
            dircheck(outDirNative);dircheck(outDirMNI);
            % Select images to be realigned
            switch (imageType)
                case 'maskGLM'
                    % first create a clean copy of the mask
                    maskDir=regDir;
                    maskName=['mask_' subj_name{s} '.nii'];
                    maskFile=fullfile(maskDir,maskName);
                    maskFile_vol=spm_vol(maskFile);
                    V=spm_read_vols(maskFile_vol);
                    maskFile=maskFile_vol;
                    % copy original file
                    maskFile.fname=fullfile(outDirNative,sprintf('maskGLM_%s.nii',subj_name{s}));
                    spm_write_vol(maskFile,V);
                    % write in the MNI destination folder
                    maskFile.fname=fullfile(outDirMNI,sprintf('mask_GLM_%s.nii',subj_name{s}));
                    spm_write_vol(maskFile,V);
                    % prepare for destination - in MNI space
                    
                    Q = sprintf('%s,1',maskFile.fname);
                    % files to delete
                    D{1} = maskFile.fname;
                case 'funcMask_BG' 
                    indx=1;
                    for h=1:2 % hemisphere
                        for r=1:length(regType)
                            % Select images to be realigned - save them
                            % with a new name
                            imageDir    = outDirNative;
                            imageName   = sprintf('funcMask_%s_%s_%s.nii',subj_name{s},regType{r},hem{h});
                            imageFile   = fullfile(imageDir,imageName);
                            imageV=spm_vol(imageFile);
                            V=spm_read_vols(imageV);
                            image=imageV;
                            % copy original file
                            % to MNI folder, change name
                            image.fname=fullfile(outDirMNI,imageName);
                            spm_write_vol(image,V);
                            % prepare for destination - in MNI space
                            % Now select real the newly created one for realignment...
                            Q{indx}   = sprintf('%s,1',image.fname);
                            D{indx}   = image.fname;
                            indx=indx+1;
                        end; % seqType
                    end
                case 'funcMask_striatum_buckner' % delete, makes no sense
                    indx=1;
                    for r=1:length(networks)
                        % Select images to be realigned - save them
                        % with a new name
                        imageDir    = outDirNative;
                        imageName   = sprintf('%s_network-%d.nii',imageType,r);
                        imageFile   = fullfile(imageDir,imageName);
                        imageV=spm_vol(imageFile);
                        V=spm_read_vols(imageV);
                        image=imageV;
                        % copy original file
                        % to MNI folder, change name
                        image.fname=fullfile(outDirMNI,imageName);
                        spm_write_vol(image,V);
                        % prepare for destination - in MNI space
                        % Now select real the newly created one for realignment...
                        Q{indx}   = sprintf('%s,1',image.fname);
                        D{indx}   = image.fname;
                        indx=indx+1;
                    end; % region          
                case 'psc'
                    indx=1;
                    for ss=sessN
                        for st=1:length(seqType)
                            % Select images to be realigned - save them
                            % with a new name
                            dir=eval(glm);
                            imageDir    = fullfile(dir{ss},subj_name{s});
                            imageName = sprintf('psc_sess%d_%s.nii',ss,seqType{st});
                            imageFile = fullfile(imageDir,imageName);
                            imageV=spm_vol(imageFile);
                            V=spm_read_vols(imageV);
                            image=imageV;
                            % copy original file
                            image.fname=fullfile(outDirNative,imageName);
                            spm_write_vol(image,V);
                            % prepare for destination - in MNI space
                            image.fname=fullfile(outDirMNI,imageName);
                            spm_write_vol(image,V);
                            % Now select real the newly created one for realignment...
                            Q{indx}   = sprintf('%s,1',image.fname);
                            D{indx} = image.fname;
                            indx=indx+1;
                        end; % seqType
                    end; % session
                case 'spmT'
                    indx=1;
                    for ss=sessN
                        for st=1:length(seqType)
                            % Select images to be realigned - save them
                            % with a new name
                            dir=eval(glm);
                            imageDir    = fullfile(dir{ss},subj_name{s});
                            imageName   = sprintf('spmT_%s.nii',seqType{st});
                            imageFile   = fullfile(imageDir,imageName);
                            imageV=spm_vol(imageFile);
                            V=spm_read_vols(imageV);
                            image=imageV;
                            % copy original file
                            imageName   = sprintf('spmT_sess%d_%s.nii',ss,seqType{st}); % add session number
                            image.fname=fullfile(outDirNative,imageName);
                            spm_write_vol(image,V);
                            % prepare for destination - in MNI space
                            image.fname=fullfile(outDirMNI,imageName);
                            spm_write_vol(image,V);
                            % Now select real the newly created one for realignment...
                            Q{indx}   = sprintf('%s,1',image.fname);
                            D{indx} = image.fname;
                            indx=indx+1;
                        end; % seqType
                    end; % session
                case 'searchlight'
                    indx=1;
                   % seqType={'trained','untrained','cross'}; % different names given
                    for ss=sessN
                        for st=1:length(seqType)
                            % Select images to be realigned - save them
                            % with a new name
                            imageDir    = outDirNative;
                            imageName   = sprintf('%s_sess%d_striatum_%s_%s.nii',subj_name{s},ss,distType,seqType{st});
                            imageFile   = fullfile(imageDir,imageName);
                            imageV=spm_vol(imageFile);
                            V=spm_read_vols(imageV);
                            image=imageV;
                            % prepare for destination - in MNI space
                            image.fname=fullfile(outDirMNI,imageName);
                            spm_write_vol(image,V);
                            % Now select real the newly created files for realignment...
                            Q{indx}   = sprintf('%s,1',image.fname);
                            D{indx}   = image.fname;
                            indx=indx+1;
                        end
                    end
                case 'betas'
                    load(fullfile(baseDir,'betaNames.mat')); % saved in 'filebetas' - all beta names
                    indx=1;
                    for ss=sessN
                        imageDir = fullfile(glmSessDir{ss},subj_name{s});
                        outDirMNI = fullfile(BGfuncDir,'MNI_betas',subj_name{s},sprintf('sess-%d',ss));
                        dircheck(outDirMNI);
                        % loop over all images
                        for i=1:size(filebetas,2)
                            imageName   = filebetas{i};
                            imageFile   = fullfile(imageDir,imageName);
                            imageV      = spm_vol(imageFile);
                            V=spm_read_vols(imageV);
                            image=imageV;
                            % prepare for destination - in MNI space
                            image.fname=fullfile(outDirMNI,imageName);
                            spm_write_vol(image,V);
                            % Now select real the newly created files for realignment...
                            Q{indx}   = sprintf('%s,1',image.fname);
                            D{indx}   = image.fname;
                            indx=indx+1;
                        end; % all betas
                    end; %session
                case 'resMS'
                    indx=1;
                    for ss=sessN
                        imageDir = fullfile(glmSessDir{ss},subj_name{s});
                        outDirMNI = fullfile(BGfuncDir,'MNI_betas',subj_name{s},sprintf('sess-%d',ss));
                        dircheck(outDirMNI);
                        % loop over all images
                        imageName   = 'ResMS.nii';
                        imageFile   = fullfile(imageDir,imageName);
                        imageV      = spm_vol(imageFile);
                        V=spm_read_vols(imageV);
                        image=imageV;
                        % prepare for destination - in MNI space
                        image.fname=fullfile(outDirMNI,imageName);
                        spm_write_vol(image,V);
                        % Now select real the newly created files for realignment...
                        Q{indx}   = sprintf('%s,1',image.fname);
                        D{indx}   = image.fname;
                        indx=indx+1;
                    end; %session
                case 'betas_RS'
                    indx=1;
                    for ss=sessN
                        imageDir = fullfile(glmFoSExDir{ss},subj_name{s});
                        outDirMNI = fullfile(BGfuncDir,'MNI_betas_RS',subj_name{s},sprintf('sess-%d',ss));
                        dircheck(outDirMNI);
                        % loop over all images
                        for i=1:200 % always 200 betas
                            imageName   = sprintf('beta_%.4d.nii',i);
                            imageFile   = fullfile(imageDir,imageName);
                            imageV      = spm_vol(imageFile);
                            V=spm_read_vols(imageV);
                            image=imageV;
                            % prepare for destination - in MNI space
                            image.fname=fullfile(outDirMNI,imageName);
                            spm_write_vol(image,V);
                            % Now select real the newly created files for realignment...
                            Q{indx}   = sprintf('%s,1',image.fname);
                            D{indx}   = image.fname;
                            indx=indx+1;
                        end; % all betas
                    end; %session
                case 'resMS_RS'
                    indx=1;
                    for ss=sessN
                        imageDir = fullfile(glmFoSExDir{ss},subj_name{s});
                        outDirMNI = fullfile(BGfuncDir,'MNI_betas_RS',subj_name{s},sprintf('sess-%d',ss));
                        dircheck(outDirMNI);
                        % loop over all images
                        imageName   = 'ResMS.nii';
                        imageFile   = fullfile(imageDir,imageName);
                        imageV      = spm_vol(imageFile);
                        V=spm_read_vols(imageV);
                        image=imageV;
                        % prepare for destination - in MNI space
                        image.fname=fullfile(outDirMNI,imageName);
                        spm_write_vol(image,V);
                        % Now select real the newly created files for realignment...
                        Q{indx}   = sprintf('%s,1',image.fname);
                        D{indx}   = image.fname;
                        indx=indx+1;
                    end; %session
            end
            
            % definitions for realignment
            J.subj.def = {subjMNItrans};
            J.woptions.bb = [-78 -112 -70; 78 76 85];
            J.woptions.vox = [2 2 2];
            J.woptions.interp = 4;
            % run realignment
            if ~strcmp(imageType,'maskGLM') % for all other images
                for i=1:size(Q,2)
                    J.subj.resample = {Q{i}}; % image to realign
                    matlabbatch{1}.spm.spatial.normalise.write=J;
                    spm_jobman('run',matlabbatch);
                end
            else % for mask
            J.subj.resample = {Q}; 
            matlabbatch{1}.spm.spatial.normalise.write=J;
            spm_jobman('run',matlabbatch);
            end
            fprintf('\n Coregistration %s image for %s done\n',imageType,subj_name{s});    
            % now delete the original copies from the MNI folder
            for d=1:size(D,2)
                delete(D{d});
            end
        end
    case 'FUNC_subjAvrg'                % make subject average maps - psc, Tmaps, searchlight (in MNI space)
        sn  = [5:9,11:31];
        sessN = [1:4];
        seqType={'TrainSeq','UntrainSeq'};
        %seqType={'trained','untrained','cross'}; % different names given
        % for repsup cases
        %seqType={'TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
       % seqType={'trained_exe1','trained_exe2','untrained_exe1','untrained_exe2','cross_exe1','cross_exe2'};
        imageType='psc'; % options: psc, spmT,searchlight
        distType = 'dist'; % dist or corrDist (used for searchlight - LDC or corrDist)
      
        vararginoptions(varargin,{'sn','sessN','imageType','sessN','seqType','distType'});
    
        outDir = fullfile(BGfuncDir,'MNI','avrg');
        dircheck(outDir);
        for ss=sessN
            if strcmp(imageType,'psc') || strcmp(imageType,'spmT')                
                for st=1:length(seqType)
                    for  s=1:numel(sn)
                        IN{s}=fullfile(BGfuncDir,'MNI',subj_name{sn(s)},sprintf('w%s_sess%d_%s.nii',imageType,ss,seqType{st}));
                    end; % subject
                    outName = sprintf('subjAvrg_%s_sess%d_%s.nii',imageType,ss,seqType{st});
                    OUT=fullfile(outDir,outName);
                    spmj_imcalc_mtx(IN,OUT,'nanmean(X)');
                end; % seqType                
            elseif strcmp(imageType,'searchlight')
                dircheck(outDir);
                for st=1:length(seqType)
                    for  s=1:numel(sn)
                        IN{s}=fullfile(BGfuncDir,'MNI',subj_name{sn(s)},sprintf('w%s_sess%d_striatum_%s_%s.nii',subj_name{sn(s)},ss,distType,seqType{st}));
                    end; % subject
                    outName = sprintf('subjAvrg_%s_striatum_sess%d_%s_%s.nii',imageType,ss,distType,seqType{st});
                    OUT=fullfile(outDir,outName);
                    spmj_imcalc_mtx(IN,OUT,'nanmean(X)');
                end; % seqType 
            end
            fprintf('Done subjAvrg: %s sess-%d \n',imageType,ss);
        end % session
    case 'FUNC_MASK_ROI_native'         % Convert ROI def (.mat) into multiple .nii files
        % make a mask from ROI definitions using the glm mask (to have the
        % same resolution of voxels during functional runs)
        % used for searchlight
        % NATIVE SPACE
        sn=4;
        nativeFuncDir=fullfile(BGfuncDir,'native');
        parcelType='striatum_buckner'; % striatum_buckner or BG-striatum
        vararginoptions(varargin,{'sn','parcelType'});
        
        for s=sn
            % load ROI definition
            load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},parcelType)));
            % loop over rois
            for roi = 1:size(R,2)
                % functional mask from GLM - whole brain, containing same voxel resolution
                mask = fullfile(nativeFuncDir,subj_name{s},sprintf('maskGLM_%s.nii',subj_name{s}));
                cd(fullfile(nativeFuncDir,subj_name{s}));
                % Save region file as nifti
                region_saveasimg(R{roi},mask,'name',sprintf('funcMask_%s.nii',R{roi}.name));
            end
            fprintf('Done %s funcMask %s definition \n',subj_name{s},parcelType);
        end
    case 'FUNC_MASK_ROI_mni'            % Convert func masks (per roi) into subject MNI space
        sn=[4:9,11:22];
        parcelType = 'BG'; % BG or striatum_buckner
        vararginoptions(varargin,{'sn','parcelType'});
        parcel = sprintf('funcMask_%s',parcelType);
        % send over to FUNC_copyFile_nativeToMNI
        sml1_imana_BG_new('FUNC_copyFile_nativeToMNI','imageType',parcel,'sn',sn);
    case 'FUNC_mni_BG_subjAvrg'         % make subject average roi files (from funcMask)
        sn=[4:9,11:25];
        parcelType = 'BG'; % BG, striatum_buckner (makes no sense)
        networks = [1:7];
        vararginoptions(varargin,{'sn','parcelType'});
        
        outDir = fullfile(BGfuncDir,'MNI','avrg');
        dircheck(outDir);
        for h=1:2
            if strcmp(parcelType,'BG')
            for i=1:numregions_BG
                for s = 1:numel(sn)                  
                    IN{s} = fullfile(BGfuncDir,'MNI',subj_name{sn(s)},...
                        ['wfuncMask_',subj_name{sn(s)},'_',regname_BG{i},'_',hem{h},'.nii']);   
                end
                OUT = fullfile(outDir,...
                    ['subjAvrg_',regname_BG{i},'_',hem{h},'.nii']);
                spmj_imcalc_mtx(IN,OUT,'mean(X)');
            end
            else
                for i=1:length(networks)
                    for s = 1:numel(sn)
                        IN{s} = fullfile(BGfuncDir,'MNI',subj_name{sn(s)},...
                            ['wfuncMask_',parcelType,'_network-',num2str(i),'.nii']);
                    end
                    OUT = fullfile(outDir,...
                        ['subjAvrg_',parcelType,'_network-',num2str(i),'.nii']);
                    spmj_imcalc_mtx(IN,OUT,'mean(X)');
                end
            end
        end  
    case 'FUNC_MASK_make_subjAvrg_roi'  % maks group average roi masks (in MNI, funcMASK)
        sn=[4:9,11:25];
        vararginoptions(varargin,{'sn'});

        rawDir=fullfile(BGfuncDir,'MNI');
        for h=1:2
            for i=1:numregions_BG
                for s = 1:numel(sn)
                    % accummulate regions for all subjects  
                    subjFile = fullfile(rawDir,subj_name{sn(s)},...
                        ['wfuncMask_',subj_name{sn(s)},'_',regname_BG{i},'_',hem{h},'.nii']);                   
                    Vol1=spm_vol(subjFile);
                    Vol2=spm_read_vols(Vol1);
                    numVox(s)=length(find(Vol2));
                end
                % calculate on average how many voxels per subject
                avrgVoxNum=round(mean(numVox));
                % load in group average 
                groupAvrg = fullfile(rawDir,'avrg',sprintf('subjAvrg_%s_%s.nii',regname_BG{i},hem{h}));
                % destination file
                dest = fullfile(rawDir,'avrg',sprintf('subjAvrg_funcMASK_%s_%s.nii',regname_BG{i},hem{h}));
                makeGroupMask(avrgVoxNum,groupAvrg,dest);
                
                fprintf('Done: group funcMASK %s %s \n',regname_BG{i},hemName{h}); 
            end
        end
    case 'FUNC_MASK_combine_subjAvrg'   % combine group roi funcMASKS into BG / striatal masks
        maskType='striatum'; % whole BG or striatum
        maskNum = 2; % 1 - overall, 2 - per hemisphere
        vararginoptions(varargin,{'maskType','maskNum'});
        
        switch(maskType)
            case 'BGall'
                reg=[1:3]; % striatum, pallidum, putamen
            case 'striatum'
                reg=[1,3]; % striatum and putamen only
        end
        % makes a mask for all of the BG regions (both hemisphere)
        indx=1;
        for h=1:2
            for i=reg
                mask{indx}=fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_funcMASK_%s_%s.nii',regname_BG{i},hem{h}));
                indx=indx+1;
            end
            if maskNum==2
                M{h}=combineMask(mask);
                indx=1; % start again
            end
        end
        M{maskNum}=combineMask(mask);
        % write the mask(s)
        for m = 1:maskNum
            BG_writeMask=spm_vol(mask{1});
            if maskNum==1
                BG_writeMask.fname=fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_funcMASK_%s.nii',maskType));
            else
                BG_writeMask.fname=fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_funcMASK_%s_%s.nii',maskType,hem{m}));
            end
            BG_writeMask.descip='avrg_wholeBG_mask';
            spm_write_vol(BG_writeMask,M{m});
        end
    case 'FUNC_funcMask_group_overall'        % calculating overall functional mask (used for masking, surfaces...)
        % combine subjAvrg mask and the buckner mask - funcMask_final
        % take 1) the 'loose Buckner mask':
        % 2) the striatum_lh / rh mask subjAvrg_funcMASK_striatum_lh
        % compute an intersection of the two
        maskType = 'striatum';
        vararginoptions(varargin,{'maskType'});
        
        M{2} = fullfile(BGanatDir,'striatum_buckner','avrg','striatum_7networks_loose.nii');
        calcFunc = ('i1&i2>0');
        for h=1:2
          M{1} = fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_funcMASK_%s_%s.nii',maskType,hem{h}));  
          outFile = fullfile(BGfuncDir,'MNI','avrg',sprintf('funcMask_overall_%s.nii',hem{h}));
          spm_imcalc_ui(M,outFile,calcFunc);
        end
        fprintf('Calculated an overlap of %s subjAvrg and Buckner maps.\n',maskType);
    case 'FUNC_createBuckner_mask'
        % mask the buckner striatal mask so it only contains voxels where
        % functional subjAvrg mask is defined; taking the loose Buckner map
        % NO NEED to do this anymore - it's done
        hemi=0; % 0 - overall or 1 - per hemi
        vararginoptions(varargin,{'hemi'});
        
        calcFunc = ('i1.*i2'); 
        M{2} = fullfile(BGanatDir,'striatum_buckner','avrg','striatum_7networks_loose.nii');
        switch hemi
            case 0
                M{1} = fullfile(BGfuncDir,'MNI','avrg','subjAvrg_funcMASK_striatum.nii');
                outFile = fullfile(BGanatDir,'striatum_buckner','avrg','striatum_buckner.nii');
                spm_imcalc_ui(M,outFile,calcFunc);
            case 1
                for h=1:2
                    M{1} = fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_funcMASK_striatum_%s.nii',hem{h}));
                    outFile = fullfile(BGanatDir,'striatum_buckner','avrg',sprintf('striatum_buckner_%s.nii',hem{h}));
                    spm_imcalc_ui(M,outFile,calcFunc);
                end
        end        
        fprintf('Masked the Buckner map to fit functional mask.\n');
    case 'FUNC_MASK_individSubj'    %    create striatal masks for each subject
        maskType='striatum'; % whole BG or striatum
        maskNum = 2; % 1 - overall, 2 - per hemisphere
        sn=[4:9,11:25];
        space = 'MNI'; % MNI or native
        vararginoptions(varargin,{'maskType','maskNum','sn'});
        
        switch(maskType)
            case 'BGall'
                reg=[1:3]; % striatum, pallidum, putamen
            case 'striatum'
                reg=[1,3]; % striatum and putamen only
        end
        % makes a mask for all of the BG regions (both hemisphere)
        for s=sn
            indx=1;
            for h=1:2
                for i=reg
                    switch space
                        case 'MNI'     
                            mask{indx} = fullfile(BGfuncDir,space,subj_name{s},sprintf('wfuncMask_%s_%s_%s.nii',subj_name{s},regname_BG{i},hem{h}));
                        case 'native'
                            mask{indx} = fullfile(BGfuncDir,space,subj_name{s},sprintf('funcMask_%s_%s_%s.nii',subj_name{s},regname_BG{i},hem{h}));
                    end
                    indx=indx+1;
                end
                if maskNum==2
                    M{h}=combineMask(mask);
                    indx=1; % start again
                end
            end
            M{maskNum}=combineMask(mask);
            % write the mask(s)
            for m = 1:maskNum
                BG_writeMask=spm_vol(mask{1});
                if maskNum==1
                    BG_writeMask.fname=fullfile(BGfuncDir,space,subj_name{s},sprintf('funcMASK_%s.nii',maskType));
                else
                    BG_writeMask.fname=fullfile(BGfuncDir,space,subj_name{s},sprintf('funcMASK_%s_%s.nii',maskType,hem{m}));
                end
                BG_writeMask.descip='avrg_wholeBG_mask';
                spm_write_vol(BG_writeMask,M{m});
            end
            fprintf('Done %s mask for %s\n',maskType,subj_name{s});
        end    
    case 'FUNC_maskImages_group'              % applies the MNI group mask to whole brain images (psc, spmT)
        sessN=[1:4];
        seqType = {'corrDist_trained','corrDist_untrained','corrDist_cross'}; 
        % TrainSeq, UntrainSeq or TrainSeq_1st / _2nd, UntrainSeq_1st / 2nd; 
        %for searhclight: dist_trained / _exe1 or corrDist_trained / _exe1
        imageType = 'searchlight_striatum'; % psc, spmT or searchlight_striatum etc.
        maskType = 'striatum_rh'; % striatum or BGall; or per hemisphere - striatum_rh/_lh
        vararginoptions(varargin,{'sn','maskType','imageType','sessN','seqType'});
        
        groupDir=fullfile(BGfuncDir,'MNI','avrg');
        for ss=sessN
            for st=1:length(seqType)
                % Select images to be realigned - save them
                % with a new name
                IN  = fullfile(groupDir,sprintf('subjAvrg_%s_sess%d_%s.nii',imageType,ss,seqType{st}));
                image=spm_vol(IN); image2=spm_read_vols(image);
               
                switch maskType
                    case 'striatum'
                        mask = fullfile(groupDir,'funcMask_overall.nii'); % striatal mask - overlap of subjAvrg striatal and Buckner mask
                    case 'striatum_rh'
                        mask = fullfile(groupDir,'funcMask_overall_rh.nii'); % striatal mask - overlap of subjAvrg striatal and Buckner mask
                    case 'striatum_lh'
                        mask = fullfile(groupDir,'funcMask_overall_lh.nii'); % striatal mask - overlap of subjAvrg striatal and Buckner mask
                end
                maskRead=spm_vol(mask); maskRead2=spm_read_vols(maskRead);
               
                newImage=image2.*maskRead2;
                BG_write=image;
                BG_write.fname=fullfile(groupDir, sprintf('subjAvrg_%s_sess%d_%s_MASK_%s.nii',imageType,ss,seqType{st},maskType));
                BG_write.descip=('masked_funcImage');
                spm_write_vol(BG_write,newImage);    
                fprintf('Done sess-%d: %s-%s\n',ss,imageType,seqType{st});
            end; % seqType
        end % sess
    case 'FUNC_maskImages_individSubj' % mask whole brain images in each subject
        sessN=[1:4];
        space = 'MNI'; % MNI or native
        seqType = {'corrDist_trained','corrDist_untrained','corrDist_cross'}; 
        % TrainSeq, UntrainSeq or TrainSeq_1st / _2nd, UntrainSeq_1st / 2nd; 
        %for searhclight: dist_trained / _exe1 or corrDist_trained / _exe1
        imageType = 'corrDist'; % corrDist,dist,psc
        maskType = {'striatum_lh','striatum_rh'}; % striatum or BGall; or per hemisphere - striatum_rh/_lh
        sn=[4:9,11:25];
        vararginoptions(varargin,{'sn','maskType','imageType','sessN','seqType','sn'});
        
        for m=1:length(maskType)
            for s=sn
                subjDir=fullfile(BGfuncDir,space,subj_name{s});
                cd(subjDir);
                for ss=sessN
                    for st=1:length(seqType)
                        % Select images to be realigned - save them
                        % with a new name
                        switch imageType
                            case 'dist'
                                IN  = fullfile(sprintf('w%s_sess%d_striatum_%s.nii',subj_name{s},ss,seqType{st}));
                            case 'corrDist'
                                IN  = fullfile(sprintf('w%s_sess%d_striatum_%s.nii',subj_name{s},ss,seqType{st}));
                            case 'psc'
                                IN  = fullfile(sprintf('wpsc_sess%d_%s.nii',ss,seqType{st}));
                        end
                        image=spm_vol(IN); image2=spm_read_vols(image);
                        mask = fullfile(subjDir,sprintf('funcMask_%s.nii',maskType{m}));
                        maskRead=spm_vol(mask); maskRead2=spm_read_vols(maskRead);
                        
                        newImage=image2.*maskRead2;
                        BG_write=image;
                        if strcmp(imageType,'psc')
                            BG_write.fname=fullfile(subjDir, sprintf('%s_sess%d_%s_%s_MASK_%s.nii',subj_name{s},ss,imageType,seqType{st},maskType{m}));
                        else
                            BG_write.fname=fullfile(subjDir, sprintf('%s_sess%d_%s_MASK_%s.nii',subj_name{s},ss,seqType{st},maskType{m}));
                        end
                        BG_write.descip=('masked_funcImage');
                        spm_write_vol(BG_write,newImage);
                        fprintf('Done %s sess-%d: %s-%s\n',subj_name{s},ss,imageType,seqType{st});
                    end; % seqType
                end % sess
            end; % subject
        end; % mask
    case 'FUNC_train_vs_untrainImage'   % calculate images of train-untrain (psc, spmT, searchlight)
        sessN=[1:4];
        seqType={'TrainSeq','UntrainSeq'};
        imageType='psc'; % psc or spmT
        maskType='striatum'; % striatum or BGall
        vararginoptions(varargin,{'sn','maskType','imageType','sessN'});
        
        groupDir=fullfile(BGfuncDir,'MNI','avrg');

        for ss=sessN
            if strcmp(imageType,'searchlight')
                tr  = fullfile(groupDir, sprintf('subjAvrg_%s_%s_sess%d_dist_trained.nii',imageType,maskType,ss));
                utr = fullfile(groupDir, sprintf('subjAvrg_%s_%s_sess%d_dist_untrained.nii',imageType,maskType,ss));
            else
                tr  = fullfile(groupDir, sprintf('subjAvrg_%s_sess%d_%s_MASK_%s.nii',imageType,ss,seqType{1},maskType));
                utr = fullfile(groupDir, sprintf('subjAvrg_%s_sess%d_%s_MASK_%s.nii',imageType,ss,seqType{2},maskType));
            end
            imageT=spm_vol(tr); imageT2=spm_read_vols(imageT);
            imageU=spm_vol(utr); imageU2=spm_read_vols(imageU);
            
            newImage=imageT2-imageU2;
            diff_write=imageT;
            diff_write.fname=fullfile(groupDir, sprintf('subjAvrg_%s_DIFF_sess%d_MASK_%s.nii',imageType,ss,maskType));
            diff_write.descrip=('diff_train-untrain');
            spm_write_vol(diff_write,newImage);
        end
    case 'FUNC_calcAcrossSess_image'
        sessTr=[1,2];
        seqType={'TrainSeq','UntrainSeq'};
        imageType='psc'; % psc or spmT
        maskType='striatum'; % striatum or BGall
        seqIndx=1; % trained
        vararginoptions(varargin,{'sn','maskType','imageType','sessTr','seqIndx'});
        
        groupDir=fullfile(BGfuncDir,'MNI','avrg');
        
        
        if strcmp(imageType,'searchlight')
            seqType={'trained','untrained'};
            t1  = fullfile(groupDir, sprintf('subjAvrg_%s_%s_sess%d_dist_%s.nii',imageType,maskType,sessTr(1),seqType{seqIndx}));
            t2  = fullfile(groupDir, sprintf('subjAvrg_%s_%s_sess%d_dist_%s.nii',imageType,maskType,sessTr(2),seqType{seqIndx}));
        else
            t1  = fullfile(groupDir, sprintf('subjAvrg_%s_sess%d_%s_MASK_%s.nii',imageType,sessTr(1),seqType{seqIndx},maskType));
            t2  = fullfile(groupDir, sprintf('subjAvrg_%s_sess%d_%s_MASK_%s.nii',imageType,sessTr(2),seqType{seqIndx},maskType));
        end
        imageT1=spm_vol(t1); imageS1=spm_read_vols(imageT1);
        imageT2=spm_vol(t2); imageS2=spm_read_vols(imageT2);
        
        newImage=imageS2-imageS1;
        diff_write=imageT1;
        diff_write.fname=fullfile(groupDir, sprintf('subjAvrg_%s_%s_diffSess_sess%d-%d_MASK_%s.nii',imageType,seqType{seqIndx},sessTr(1),sessTr(2),maskType));
        diff_write.descrip=(sprintf('diff_sess%d-%d',sessTr(1),sessTr(2)));
        spm_write_vol(diff_write,newImage);
    case 'FUNC_calcAcrossSess_trainVsUntrain_image'
        sessTr=[1,2];
        seqType={'TrainSeq','UntrainSeq'};
        imageType='psc'; % psc or spmT
        maskType='striatum'; % striatum or BGall
        seqIndx=1; % trained
        vararginoptions(varargin,{'sn','maskType','imageType','sessTr','seqIndx'});
        
        groupDir=fullfile(BGfuncDir,'MNI','avrg');
        
        
        if strcmp(imageType,'searchlight')
            seqType={'trained','untrained'};
        end
            t1  = fullfile(groupDir, sprintf('subjAvrg_%s_%s_diffSess_sess%d-%d_MASK_%s.nii',imageType,seqType{1},sessTr(1),sessTr(2),maskType));
            t2  = fullfile(groupDir, sprintf('subjAvrg_%s_%s_diffSess_sess%d-%d_MASK_%s.nii',imageType,seqType{2},sessTr(1),sessTr(2),maskType));

        imageT1=spm_vol(t1); imageS1=spm_read_vols(imageT1);
        imageT2=spm_vol(t2); imageS2=spm_read_vols(imageT2);
        
        newImage=imageS1-imageS2;
        diff_write=imageT1;
        diff_write.fname=fullfile(groupDir, sprintf('subjAvrg_%s_diffSeqType_diffSess_sess%d-%d_MASK_%s.nii',imageType,sessTr(1),sessTr(2),maskType));
        diff_write.descrip=(sprintf('diff_trainVSuntrain_sess%d-%d',sessTr(1),sessTr(2)));
        spm_write_vol(diff_write,newImage);
        
    case 'ROI_network_MNItoNative'
        % transfer the MNI parcelation funcStriatum to subject native space
        sn=4;
        parcelType = 'striatum_buckner'; % striatum_buckner (_rh / _lh) or funcStriatum
        parcelDir = 'striatum_buckner'; % only because striatum_buckner_lh / _rh are in the folder of striatum_buckner
        vararginoptions(varargin,{'sn','parcelType'});
          
        spm('defaults','fmri');
        spm_jobman('initcfg');
        
        MNIfile=fullfile(BGanatDir,parcelDir,'avrg',sprintf('%s.nii',parcelType));
        
        for s=sn
            % inverse file: MNI -> native space
            subjAnatInv = fullfile(anatomicalDir,subj_name{s},sprintf('iy_%s_anatomical.nii',subj_name{s})); % inverse transformation
            
            outDir=fullfile(BGanatDir,parcelDir,subj_name{s});
            dircheck(outDir);
            
            % read MNI file
            MNIFile_vol=spm_vol(MNIfile);
            V=spm_read_vols(MNIFile_vol);
            nativeFile=MNIFile_vol;
            % copy original file
            nativeFile.fname=fullfile(outDir,sprintf('%s_%s.nii',subj_name{s},parcelType));
            spm_write_vol(nativeFile,V);
            % prepare for transformation
            Q = sprintf('%s,1',nativeFile.fname);
            J.subj.def = {subjAnatInv};
            J.woptions.bb = [-78 -112 -70; 78 76 85];
            J.woptions.vox = [2 2 2];
            J.woptions.interp = 4;
            
            J.subj.resample = {Q}; % image to realign
            matlabbatch{1}.spm.spatial.normalise.write=J;
            spm_jobman('run',matlabbatch);
           
            % now delete the copied file (in MNI space), leave only native
            % space version
            delete(fullfile(outDir,sprintf('%s_%s.nii',subj_name{s},parcelType)));
            fprintf('\n %s for %s native space done\n',parcelType,subj_name{s});    
        end     
    case 'ROI_define_BG' 
        sn=4;
        regType='striatum_buckner_hemi'; % BG, striatum, networks, striatum_buckner (change the other networks name)
        vararginoptions(varargin,{'sn','regType'});
        regFold = 'striatum_buckner'; % useful for striatum_buckner_hemi
        switch(regType)
            case 'BG'
                regNum=[1,2,3];
            case 'BG-striatum'
                regNum=[1,3];
            case 'BG-striatum-networks'
                regNum=[1:7];
            case 'striatum_buckner'
                regNum=[1:7];
            case 'striatum_buckner_hemi'
                regNum=[1:7];
        end      
        for s=sn
            R=[];
            if ~strcmp(regType,'BG-striatum-networks') & ~strcmp(regType,'striatum_buckner') % for BG, striatum - do per hemisphere
                indx=1;
                for h=1:2
                    for j=regNum
                        % Get basal ganglia
                        if strcmp(regType,'striatum_buckner_hemi')
                            INdir = fullfile(BGanatDir,regFold,subj_name{s});
                            fileBG = fullfile(INdir,sprintf('w%s_%s_%s.nii',subj_name{s},regFold,hem{h})); % functional parcellation of striatum from Carlos or Buckner
                            R{indx}.name = sprintf('%s_network-%d_%s',regFold,j,hem{h});
                            R{indx}.value = j;
                        else
                            fileBG = fullfile(BGanatDir,'native',subj_name{s},sprintf('%s_%s_%s.nii', subj_name{s}, regname_BG{j}, hem{h}));
                            R{indx}.name = [subj_name{s} '_' regname_BG{j} '_' hem{h}];
                            R{indx}.value = 1;
                        end
                        R{indx}.type = 'roi_image';
                        R{indx}.file= fileBG;
                        indx=indx+1;
                    end
                end
                R=region_calcregions(R);
            else % for networks - across hemispheres
                INdir = fullfile(BGanatDir,regType,subj_name{s});
                fileBG = fullfile(INdir,sprintf('w%s_%s.nii',subj_name{s},regType)); % functional parcellation of striatum from Carlos or Buckner
                for j=regNum
                    R{j}.type = 'roi_image';
                    R{j}.file = fileBG;
                    R{j}.name = sprintf('%s_network-%d',regType,j);
                    R{j}.value = j;
                end
               % subjToMNI=fullfile(anatomicalDir,subj_name{s},sprintf('iy_%s_anatomical.nii',subj_name{s}));
               % R=region_calcregions(R,'voxelspace',subjToMNI);
               R=region_calcregions(R);
            end   
            cd(regDir);
            save([subj_name{s} sprintf('_%s_regions.mat',regType)],'R');
            
            fprintf('\nROIs %s have been defined for %s \n',regType,subj_name{s});
        end
    case 'ROI_define_BG_group'
        % define ROI in MNI space for the group
        regType='striatum_buckner'; % BG, striatum, networks
        vararginoptions(varargin,{'sn','regType'});
        regNum = [1:7];
        R=[];
        
        indx=1;
        for h=1:2
            fileBG = fullfile(BGanatDir,regType,'avrg',sprintf('%s_%s.nii',regType,hem{h}));
            for j=regNum
                R{indx}.type = 'roi_image';
                R{indx}.file = fileBG;
                R{indx}.name = sprintf('group_%s_network-%d_%s',regType,j,hem{h});
                R{indx}.value = j;
                indx=indx+1;
            end
        end
        R=region_calcregions(R);
        
        dircheck(regDir);
        cd(regDir);
        save(sprintf('group_MNI_%s_regions.mat',regType),'R');
        
        fprintf('\nROIs %s have been defined for group \n',regType);
    case 'ROI_define_thalamus'
        sn=4;
        regType='thalamus'; % BG, striatum, networks, striatum_buckner (change the other networks name)
        vararginoptions(varargin,{'sn','regType'});       
        for s=sn
            R=[];
            for h=1:2
                % get thalamus
                fileT = fullfile(BGanatDir,'native',subj_name{s},sprintf('%s_Thalamus_%s.nii', subj_name{s},  hem{h}));
                R{h}.name = [subj_name{s} '_Thalamus_' hem{h}];
                R{h}.value = 1;
                R{h}.type = 'roi_image';
                R{h}.file= fileT;
            end
            R=region_calcregions(R);          
            cd(regDir);
            save([subj_name{s} sprintf('_%s_regions.mat',regType)],'R');
            
            fprintf('\nROIs %s defined for: %s \n',regType,subj_name{s});
        end
        
    case 'ROI_make_nii'                                                     % OPTIONAL   :  Convert ROI def (.mat) into multiple .nii files (to check!)
        parcelType   = 'striatum_buckner';
        parcelFolder = 'striatum_buckner';
        vararginoptions(varargin,{'sn','parcelType','parcelFolder'});
        
        for s=sn
            % load ROI definition
            load(fullfile(regDir,sprintf('%s_%s_regions.mat',subj_name{s},parcelType)));
            
            % loop over rois
            for roi = 1:size(R,2)
                % mask volume
                mask = fullfile(regDir,sprintf('mask_%s.nii',subj_name{s}));
                % Save region file as nifti
                cd(fullfile(BGanatDir,parcelFolder,subj_name{s}));
                region_saveasimg(R{roi},mask);
            end   
        end
    case 'ROI_betas' % -------- save betas per glm (per session) - for psc / dist
        sessN = 1;
        sn  = [4:9,11:28];
        parcelType='BG-striatum'; % options: 1) BG, 2) BG-striatum, 3) BG-striatum-networks
        vararginoptions(varargin,{'sn','sessN','roi','type','parcelType'});
        
        for ss=sessN
            % harvest
            for s=sn % for each subj
                T=[];
                fprintf('\nSubject: %d\n',s) % output to user
                cd (fullfile(glmSessDir{ss},subj_name{s}));
                
                % load files
                load(fullfile(glmSessDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)])); % load subject's region parcellation (R)
                roi=1:size(R,2);
                
                % Add a few extra images
                %----task against rest
                O{1}=sprintf('psc_sess%d_TrainSeq.nii',ss); %psc trained
                O{2}=sprintf('psc_sess%d_UntrainSeq.nii',ss); %psc untrained
                O{3}='spmT_TrainSeq.nii';
                O{4}='spmT_UntrainSeq.nii';
                % searchlights
             %   searchDir=fullfile(BGfuncDir,'native',subj_name{s});
             %   O{5}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_trained.nii',subj_name{s},ss));
             %   O{6}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_untrained.nii',subj_name{s},ss));
             %   O{7}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_cross.nii',subj_name{s},ss));
                
                P=SPM.Vbeta(SPM.xX.iC);
                V = SPM.xY.VY;
                oP=spm_vol(char(O));
                
                for r = roi % for each region
                    % get raw data for voxels in region
                    Y = region_getdata(V,R{r});  % Data Y is N x P
                    data = region_getdata(oP,R{r}); % from added images
                    % voxel position
                    S.volcoord = {R{r}.data'};
                    % estimate region betas
                    [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                    S.betaW                   = {betaW};                             % multivariate pw
                    S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))}; % univariate pw
                    S.betaRAW                 = {beta};
                    S.resMS                   = {resMS};
                    
                    % info from maps for psc
                    S.psc_train    = {data(1,:)};
                    S.psc_untrain  = {data(2,:)};
                    S.Tmap_train   = {data(3,:)};
                    S.Tmap_untrain = {data(4,:)};
              %      S.dist_train   = {data(5,:)};
              %      S.dist_untrain = {data(6,:)};
              %      S.dist_cross   = {data(7,:)};
                    
                    S.SN                      = s;
                    S.region                  = r;
                    T = addstruct(T,S);
                    fprintf('%d.',r)
                end
                dircheck(fullfile(betaDir,subj_name{s}));
                save(fullfile(betaDir,subj_name{s},sprintf('betas_%s_%s_sess%d.mat',parcelType,subj_name{s},ss)),'-struct','T');
                fprintf('\nDone beta extraction for sess%d-%s\n',ss,subj_name{s});
            end
            % save T
            save(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss)),'-struct','T');
            fprintf('\n');
        end
    case 'ROI_stats'
        sessN = [1:4];
        sn  = [4:9,11:25];
        parcelType = 'BG-striatum'; % BG, BG-striatum, BG-striatum-networks
        betaChoice = 'uni'; % uni, multi or raw
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        vararginoptions(varargin,{'sn','sessN','roi','betaChoice','type','parcelType'});
        
        for ss=sessN;
            T = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            
            % output structures / or load from before
            switch(type)
                case 'new'
                    To=[];
                case 'add'
                    To=load(fullfile(BGstatsDir,sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
            end
            
            % do stats
            for s = sn % for each subject
                D = load(fullfile(glmSessDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                fprintf('\nSubject: %d session: %d\n',s,ss)
                num_run = numruns_task_sess;
                
                roi=unique(T.region)';  % for each region
                switch(parcelType) % for determining region side and type
                    case 'BG'
                        regSide=[1 1 1 2 2 2];
                        regType=[1 2 3 1 2 3];
                    case 'BG-striatum'
                        regSide=[1 1 2 2];
                        regType=[1 2 1 2];
                end
                for r = roi
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
                    [G,Sig]     = pcm_estGCrossval(betaW(1:(12*num_run),:),D.run,D.seqNumb);   
                    So.IPM      = rsa_vectorizeIPM(G);
                    So.Sig      = rsa_vectorizeIPM(Sig);
                    % squared distances
                    So.RDM      = rsa.distanceLDC(betaW,D.run,D.seqNumb);
                    
                    % calculate G and dist
                    Do = calcDist(D,betaW,G);
                    
                    % stats from additional images
                    So.psc_train        = nanmean(S.psc_train{:});
                    So.psc_untrain      = nanmean(S.psc_untrain{:});
                    So.Tmap_train       = nanmean(S.Tmap_train{:});
                    So.Tmap_untrain     = nanmean(S.Tmap_untrain{:});
                    So.search_train     = nanmean(S.dist_train{:});
                    So.search_untrain   = nanmean(S.dist_untrain{:});
                    So.search_cross     = nanmean(S.dist_cross{:});
                    
                    % indexing fields
                    So.SN       = s;
                    So.region   = r;
                    if ~strcmp(parcelType,'BG-striatum-networks');
                        So.regSide  = regSide(r);
                        So.regType  = regType(r);
                    end
                    % data structure
                    To          = addstruct(To,So); % indexing fields, other images
                    To          = addstruct(To,Do); % distances
                    
                end; % each region
            end; % each subject
            
            % % save - stats data and simulations
            save(fullfile(BGstatsDir,sprintf('stats_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)),'-struct','To');
            
            fprintf('\nDone.\n')
        end
    case 'ROI_patternConsist_old'
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 

        % enter sn, region,  beta: 0=betaW, 1=betaU, 2=raw betas
        % (1) Set parameters
        sn  = [4:9,11:22];
        roi = [1:3,5:7];
        sessN=1;
        removeMean = 'yes'; % are we removing pattern means for patternconsistency?
        betaChoice='multiPW'; % multiPW, uniPW, raw
        parcelType='BG-striatum';
        vararginoptions(varargin,{'sn','roi','removeMean','sessN','betaChoice','parcelType'});
        
        numRun=numruns_task_sess;
        numSeq=numel(num_seq);
        
        if strcmp(removeMean,'yes')
             rmean = 1; % we are removing the mean
        else rmean = 0; % we are keeping the mean (yields higher consistencies but these are biased)
        end

        RR=[];
        switch(parcelType) % for determining region side and type
            case 'BG'
                regSide=[1 1 1 2 2 2];
                regType=[1 2 3 1 2 3];
            case 'BG-striatum'
                regSide=[1 1 2 2];
                regType=[1 2 1 2];
        end
        %========%
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_%s_sess%d',parcelType,ss))); % loads in struct 'T'
            roi=unique(T.region)';
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
                    [R.consist_crossval R.consist_voxel] = rsa_patternConsistency_crossval(beta,partVec,condVec,'removeMean',rmean,'voxSpecific',1);
                    R.sn =s;
                    R.sessN=ss;
                    R.roi=r;
                     if ~strcmp(parcelType,'BG-striatum-networks');
                        R.regSide  = regSide(r);
                        R.regType  = regType(r);
                    end
                    
                    RR=addstruct(RR,R);
                end
            end
        end
               
        % save consistency
        statsDir=fullfile(BGDir,'stats');
        dircheck(statsDir);
        save(fullfile(statsDir,sprintf('patternConsist_%s_%sbetas',parcelType,betaChoice)),'-struct','RR');
        
        keyboard;
        figure
        subplot(2,3,[1:3]);
        plt.bar(RR.regType,RR.consist,'split',RR.regSide,'leg',{'Contra','Ipsi'},'leglocation','northeast');
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
    case 'ROI_patternConsist_plot'
        parcelType='BG-striatum';
        betaChoice='multiPW';
        vararginoptions(varargin,{'parcelType','betaChoice'});
        
        B=load(fullfile(BGDir,'stats',sprintf(sprintf('patternConsist_%s_%sbetas',parcelType,betaChoice))));
         
        [regLab, ~] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);
        
        for p=1:2 % crossvalidated and normal
            if p==1
                var=B.consist;
                type='normal';
            else
                var=B.consist_crossval;
                type='crossval';
            end
            
            roi=unique(B.roi)';
            figure
            for r=roi
                subplot(1,max(roi),r)
                if ~strcmp(parcelType,'BG-striatum-networks')
                    plt.hist(var,'subset',B.regType==r,'split',B.regSide,'leg',{'contra','ipsi'},'leglocation','north');
                else
                    plt.hist(var,'subset',B.roi==r);
                end
                if r==1
                    xlabel(sprintf('%s pattern consistency %s',type,betaChoice));
                end
                title(sprintf('%s',regLab{r}));
            end
            
        end;

    case 'ROI_pattern_consist'  
        % pattern consistency for specified roi
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 

        % enter sn, region,  beta: 0=betaW, 1=betaU, 2=raw betas
        % (1) Set parameters
        sn  = [4:9,11:25];
        roi = [1:4];
        sessN=[1:4];
        parcelType='BG-striatum';
        betaChoice='multiPW'; % multiPW, uniPW, raw
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice','parcelType'});
        
        numRun=numruns_task_sess;
        numCond=numel(num_seq);
        
        RR=[];
        %========%
        for ss=sessN
            %T = load(fullfile(regDir,sprintf('betas_%s_sess%d',parcelType,ss))); % loads in struct 'T'
            %T = load(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d',parcelType,ss))); % loads in struct 'T'
            T = load(fullfile(regDir,sprintf('betas_FoSEx_%s_sess%d',parcelType,ss))); % loads in struct 'T'
            for r=roi
                for s=sn
                    S = getrow(T,(T.SN==s & T.region==r));
                    TI = load(fullfile(glmFoSExDir{ss},subj_name{s},'SPM_info'));
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
                    %partVec = kron([1:numRun]',ones(numCond,1));
                    %condVec = kron(ones(numRun,1),[1:numCond]');
                    partVec = TI.run;
                    condVec = TI.seqNumb;
                    
                    % calculate the pattern consistency
                    R.cnr                           = cnr_QC(beta,res,numCond,numRun);  % contrast to noise ratio
                    R.r2_rm                         = rsa_patternConsistency(beta,partVec,condVec,'removeMean',1);  % pattern consistency
                    R.r2                            = rsa_patternConsistency(beta,partVec,condVec,'removeMean',0);
                    [R.r2_cross_rm, R.r_cross_rm]   = rsa_patternConsistency_crossval(beta,partVec,condVec,'removeMean',1); % correlation of patterns
                    [R.r2_cross, R.r_cross]         = rsa_patternConsistency_crossval(beta,partVec,condVec,'removeMean',0);
                    
                    R.sn =s;
                    R.sessN=ss;
                    R.region=r;
                    if r<=length(roi)/2
                        R.regType=r;
                        R.hemi=1;
                    else
                        R.regType=r-(length(roi)/2);
                        R.hemi=2;
                    end
                    
                    RR=addstruct(RR,R);
                end
            end
            
        end
        
        %save the structure
        save(fullfile(QCDir,sprintf('QualityControl_%s_%sBetas',parcelType,betaChoice)),'-struct','RR');
    case 'ROI_PLOT_pattern_consist' 
       var = 'cnr'; % cnr, r2_rm, r2, r2_cross, r_cross, r2_cross_rm, r_cross_rm
       betaChoice = 'uniPW';
       roi=[1:2];
       parcelType='striatum_buckner_hemi';
       vararginoptions(varargin,{'var','betaChoice','roi','parcelType'});
       
       T = load(fullfile(QCDir,sprintf('QualityControl_%s_%sBetas',parcelType,betaChoice)));
       
       % split by hemisphere
       figure
       plt.bar(T.regType,T.(var),'split',T.hemi,'leg',{'Contra','Ipsi'},'leglocation','northeast');
       ylabel(sprintf('%s %s',var,betaChoice));
       xlabel('ROIs');

       % split by session
       figure
       plt.bar(T.regType,T.(var),'subset',T.hemi==1,'split',T.sessN,'leg',{'sess1','sess2','sess3','sess4'},'leglocation','northeast');
       ylabel(sprintf('%s %s',var,betaChoice));
       xlabel('ROIs');

            
    case 'CALC_psc' % ------------------- PSC, dist (Mahalanobis, correlation) ------------
        sn = [4:9,11:25];
        sessN = [1:4];
        parcelType = 'BG-striatum'; %BG, BG-striatum or BG-striatum-networks

        vararginoptions(varargin,{'sn','roi','seq','sessN','parcelType'});

        Stats = [];
        for ss=sessN
            T = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',parcelType,ss))); % loads region data (T)
            roi=unique(T.region)';
            for s=1:numel(sn)
                for r=roi            
                    S.psc(1,:)=nanmean(T.psc_train{(T.region==r & T.SN==sn(s)),:});
                    S.psc(2,:)=nanmean(T.psc_untrain{(T.region==r & T.SN==sn(s)),:});
                    S.seqType=[1;2];
                    S.sn=ones(size(S.seqType))*sn(s);
                    S.roi=ones(size(S.seqType))*r;
                    if ~strcmp(parcelType,'BG-striatum-networks')
                        if r>numel(roi)/2
                            S.regSide=[2;2];
                        else
                            S.regSide=[1;1];
                        end
                        S.regType=rem(S.roi,numel(roi)+1);
                    end
                    S.sessN=ones(size(S.seqType))*ss; 
                    Stats=addstruct(Stats,S);
                end
            end
        end
           
        % save structure
        statsDir=fullfile(BGDir,'stats');
        dircheck(statsDir);
        save(fullfile(statsDir,sprintf('psc_%s.mat',parcelType)),'-struct','Stats');
    case 'PLOT_psc'
        sessN=[1:4];
        parcelType = 'BG-striatum'; %BG, BG-striatum or BG-striatum-networks    
        vararginoptions(varargin,{'sessN','parcelType'});
        
        T=load(fullfile(BGDir,'stats',sprintf('psc_%s.mat',parcelType)));
        
        [regLab, hemiLab] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);

        roi=unique(T.roi)'; 
        roi=[1:2];
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
            title(sprintf('%s-%s',regLab{r},hemiLab{r}));
            if r==1
                ylabel('Percent signal change');
                xlabel('Session');
            else
                ylabel('');
            end
            plt.match('y');
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
        sn = [4:9,11:25];
        sessN = [1:4];
        betaChoice='multiPW';
        parcelType='BG-striatum'; % BG-striatum, BG-striatum-networks
        vararginoptions(varargin,{'sn','roi','sessN','betaChoice','fig','parcelType'});
        
        SAll = []; STAll = [];
        for ss=sessN
            D = load(fullfile(BGstatsDir,sprintf('stats_%s_%s_sess%d.mat',parcelType,betaChoice,ss))); % loads region data (D)
            roi = unique(D.region)';
            for s=1:numel(sn)
                for r=roi
                    % get relevant rows for subject, region
                    d=getrow(D,D.SN==sn(s)&D.region==r);
                    % extract values of interest for train / untrain
                    S.dist(1,:)     = d.dist_train;
                    S.dist(2,:)     = d.dist_untrain;
                    S.search(1,:)   = d.search_train;
                    S.search(2,:)   = d.search_untrain;
                    S.seqType       = [1;2];
                    S.sn            = [sn(s);sn(s)];
                    S.roi           = [r;r];
                    if ~strcmp(parcelType,'BG-striatum-networks')
                        if r>numel(roi)/2
                            regSide=2;
                            regType=r-numel(roi)/2;
                        else
                            regSide=1;
                            regType=r;
                        end
                        S.regSide   = [regSide;regSide];
                        S.regType   = [regType;regType];
                        ST.regType  = regType;
                        ST.regSide  = regSide;
                    end           
                    S.sessN         = [ss;ss];
                    SAll            = addstruct(SAll,S);
                    
                    % extract values of interest for seqType
                    ST.dist_seqType     = d.dist_cross;
                    ST.search_seqType   = d.search_cross;
                    ST.sn               = s;
                    ST.roi              = r;
                    ST.sessN            = ss;
                    STAll               = addstruct(STAll,ST);    
                end
            end
        end
        % save structure
        save(fullfile(BGstatsDir,sprintf('dist_%s_BG.mat',parcelType)),'-struct','SAll');
        save(fullfile(BGstatsDir,sprintf('dist_seqType_%s_BG.mat',parcelType)),'-struct','STAll');
    case 'PLOT_dist'
        sessN=[1:4];
        parcelType='BG-striatum'; % BG-striatum, BG-striatum-networks
        distType='mahalanobis'; % mahalanobis or searchlight
        
        vararginoptions(varargin,{'sessN','parcelType','distType'});
        
        D=load(fullfile(BGstatsDir,sprintf('dist_%s_BG.mat',parcelType)));
        roi=unique(D.roi)';
        
        switch(distType)
            case 'mahalanobis'
                var=ssqrt(D.dist);
            case 'searchlight'
                var=D.search;
        end
        
        [regLab, hemiLab] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);
        sub=0;
        figure
        for r=roi
            sub=sub+1;
            subplot(1,numel(roi),sub)
            plt.line([D.sessN>3 D.sessN],var,'split',D.seqType,'subset',D.roi==r,'leg',{'trained','untrained'},'leglocation','north','style',stySeq);         
            
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s-%s',regLab{r},hemiLab{r}));
            if r==1
                ylabel('Distance');
                xlabel('Session');
            else
                ylabel('');
            end
        end
    case 'PLOT_dist_seqType'
        sessN=[1:4];
        parcelType='BG-striatum'; % BG-striatum, BG-striatum-networks
        distType='mahalanobis'; % mahalanobis or searchlight
        
        vararginoptions(varargin,{'roi','sessN','parcelType','distType'});
        
        D=load(fullfile(BGstatsDir,sprintf('dist_seqType_%s_BG.mat',parcelType)));
        roi=unique(D.roi)';
        
        switch(distType)
            case 'mahalanobis'
                var=ssqrt(D.dist_seqType);
            case 'searchlight'
                var=D.search_seqType;
        end
        
        [regLab, hemiLab] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);   
        figure
        sub=0;
        for r=roi
            sub=sub+1;
            subplot(1,numel(roi),sub)
            plt.line([D.sessN>3 D.sessN],var,'subset',D.roi==r,'leg',{'seqType'},'leglocation','north');
            
            drawline(0,'dir','horz');
            plt.match('y');
            title(sprintf('%s-%s',regLab{r},hemiLab{r}));
            if r==1
                ylabel('Mahalanobis distance');
                xlabel('Session');
            else
                ylabel('');
            end
        end
    
    case 'CALC_corrDist'
        reg = [1:3,5:7];
        sn  = [4:9,11:16];
        sessN = [1:4];
        regType='native'; % 1) native or group: 2)subjAvrg, 3)SPM (MNI, MNI_step2, SPM or FNIRT) - DO NATIVE OR GROUP_AVERAGE ONLY
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        vararginoptions(varargin,{'sn','reg','sessN','subtract_mean','regType'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
          %  D   = load(fullfile(regDir,sprintf('betas_BG_%s_sess%d.mat',regType,ss)));
            D   = load(fullfile(regDir,sprintf('betas_%s_sess%d.mat',regType,ss)));
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
        save(fullfile(BGstatsDir,sprintf('corrDist_%s_ROI.mat',regType)),'-struct','SAll');
        save(fullfile(BGstatsDir,sprintf('corrDist_seqType_%s_ROI.mat',regType)),'-struct','STAll');
    case 'PLOT_corrDist'
       % roi=[1:3,5:7];
        roi=[1:4];
        sessN=[1:4];
        %regType='native';
        regType='BG-striatum';
        
        vararginoptions(varargin,{'roi','sessN','regType'});
        
        T=load(fullfile(BGstatsDir,sprintf('corrDist_%s_ROI.mat',regType)));
        
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
         %   title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Correlation distance');
                xlabel('Session');
            else
                ylabel('');
            end
        indx=indx+1;    
        end
    case 'PLOT_corrDist_seqType'
       % roi=[1:3,5:7];
        roi=[1:4];
        sessN=[1:4];
        %regType='native';
        regType='BG-striatum';
        
        vararginoptions(varargin,{'roi','sessN','regType'});
        
        D=load(fullfile(BGstatsDir,sprintf('corrDist_seqType_%s_ROI.mat',regType)));
        figure
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            if numel(sessN)==4
                plt.line([D.sessN>3 D.sessN],ssqrt(D.corr_seqType),'subset',D.roi==r);
            elseif numel(sessN)==3
                plt.line(D.sessN,ssqrt(D.corr_seqType),'subset',D.roi==r&D.sessN<4,'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            elseif numel(sessN)==2
                plt.line(D.sessN,ssqrt(D.corr_seqType),'subset',D.roi==r&(D.sessN==1|D.sessN==4),'leg',{'trained','untrained'},'leglocation','north','style',stySeqType);
            end
            drawline(0,'dir','horz');
            plt.match('y');
          %  title(sprintf('%s',regname_BG{rem(r,4)}));
            if r==1
                ylabel('Correlation distance betweeen seqTypes');
                xlabel('Session');
            else
                ylabel('');
            end
        indx=indx+1;    
        end
    
    case 'ROI_betas_RepSup' % --------- save betas for each FoSEx glm - repetition suppression
        sessN = [1:4];
        sn  = [4:9,11:25];    
        type = 'new'; % new or add - if creating from scratch (no subject or adding new ones only)
        parcelType='striatum_buckner_hemi'; % options: 1) BG, 2) BG-striatum, 3) BG-striatum-networks
        vararginoptions(varargin,{'sn','sessN','roi','type','parcelType'});
        
        for ss=sessN
            switch(type)
                case 'new'
                    T=[];
                case 'add'
                    T=load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',ss)));
            end
            if exist(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss)))
                fprintf('ROI betas already defined for sess-%d - skipping\n',ss);
            else
                % harvest
                for s=sn % for each subj
                    fprintf('\nSubject: %d\n',s) % output to user
                    
                    % load files
                    load(fullfile(glmFoSExDir{ss}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
                    load(fullfile(regDir,[subj_name{s} sprintf('_%s_regions.mat',parcelType)])); % load subject's region parcellation (R)
                    roi=1:size(R,2);
                    
                    V = SPM.xY.VY;
                    glmSubjDir = fullfile(glmFoSExDir{ss},subj_name{s});
                    cd(glmSubjDir);
                    
                    % Add a few extra images
                    %----task against rest
                    O{1}=sprintf('psc_sess%d_TrainSeq_1st.nii',ss); %psc trained - 1st execution
                    O{2}=sprintf('psc_sess%d_UntrainSeq_1st.nii',ss); %psc untrained - 1st execution
                    O{3}=sprintf('psc_sess%d_TrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                    O{4}=sprintf('psc_sess%d_UntrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                    O{5}='spmT_TrainSeq_1st.nii';
                    O{6}='spmT_TrainSeq_2nd.nii';
                    O{7}='spmT_UntrainSeq_1st.nii';
                    O{8}='spmT_UntrainSeq_2nd.nii';
                    % searchlights
                    searchDir=fullfile(BGfuncDir,'native',subj_name{s});
                    O{9}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_trained_exe1.nii',subj_name{s},ss));
                    O{10}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_trained_exe2.nii',subj_name{s},ss));
                    O{11}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_untrained_exe1.nii',subj_name{s},ss));
                    O{12}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_untrained_exe2.nii',subj_name{s},ss));
                    O{13}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_cross_exe1.nii',subj_name{s},ss));
                    O{14}=fullfile(searchDir,sprintf('%s_sess%d_striatum_dist_cross_exe2.nii',subj_name{s},ss));
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
                        S.psc_train_1st         = {data(1,:)};
                        S.psc_untrain_1st       = {data(2,:)};
                        S.psc_train_2nd         = {data(3,:)};
                        S.psc_untrain_2nd       = {data(4,:)};
                        S.Tmap_train_1st        = {data(5,:)};
                        S.Tmap_train_2nd        = {data(6,:)};
                        S.Tmap_untrain_1st      = {data(7,:)};
                        S.Tmap_untrain_2nd      = {data(8,:)};
                        S.search_train_1st      = {data(9,:)};
                        S.search_train_2nd      = {data(10,:)};
                        S.search_untrain_1st    = {data(11,:)};
                        S.search_untrain_2nd    = {data(12,:)};
                        S.search_cross_1st      = {data(13,:)};
                        S.search_cross_2nd      = {data(14,:)};
                        % add searchlight
                        S.SN                      = s;
                        S.region                  = r;
                        T = addstruct(T,S);
                        fprintf('%d.',r)
                    end
                end
                % save T
                save(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss)),'-struct','T');
                fprintf('\n');
            end
        end
    case 'ROI_betas_RepSup_group_MNI'
        % extract betas from group MNI space
        sessN   = [1:4];
        sn      = [4:9,11:25];
        type    = 'new';
        parcelType = 'group_MNI_striatum_buckner';
        
        load(fullfile(regDir,sprintf('%s_regions.mat',parcelType))); % load group region parcellation (R)
        roi=1:size(R,2);
        
        for ss=sessN
            switch(type)
                case 'new'
                    T=[];
                case 'add'
                 %   T=load(fullfile(regDir,sprintf('betas_BG_FoSEx_sess%d.mat',ss)));
            end
            % harvest
            for s=sn % for each subj
                fprintf('\nSubject: %d\n',s) % output to user       
                % 1) load betas
                betaDir = fullfile(BGfuncDir,'MNI_betas_RS',subj_name{s},sprintf('sess-%d',ss));
                rfile = fullfile(betaDir,'wResMS.nii');
                numBetas = 200; 
                for n=1:numBetas
                    B{n}=fullfile(betaDir,sprintf('wbeta_%.4d.nii',n));
                end
                betas = spm_vol(char(B));
                res   = spm_vol(char(rfile));
                % 2) load contrasts
                contrastDir = fullfile(BGfuncDir,'MNI',subj_name{s});
                cd(contrastDir);
                %----task against rest
                O{1}=sprintf('wpsc_sess%d_TrainSeq_1st.nii',ss); %psc trained - 1st execution
                O{2}=sprintf('wpsc_sess%d_UntrainSeq_1st.nii',ss); %psc untrained - 1st execution
                O{3}=sprintf('wpsc_sess%d_TrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                O{4}=sprintf('wpsc_sess%d_UntrainSeq_2nd.nii',ss); %psc trained - 2nd execution
                O{5}=sprintf('wspmT_sess%d_TrainSeq_1st.nii',ss);
                O{6}=sprintf('wspmT_sess%d_TrainSeq_2nd.nii',ss);
                O{7}=sprintf('wspmT_sess%d_UntrainSeq_1st.nii',ss);
                O{8}=sprintf('wspmT_sess%d_UntrainSeq_2nd.nii',ss);
                %----searchlights
                O{9}=sprintf('w%s_sess%d_striatum_dist_trained_exe1.nii',subj_name{s},ss);
                O{10}=sprintf('w%s_sess%d_striatum_dist_trained_exe2.nii',subj_name{s},ss);
                O{11}=sprintf('w%s_sess%d_striatum_dist_untrained_exe1.nii',subj_name{s},ss);
                O{12}=sprintf('w%s_sess%d_striatum_dist_untrained_exe2.nii',subj_name{s},ss);
                O{13}=sprintf('w%s_sess%d_striatum_dist_cross_exe1.nii',subj_name{s},ss);
                O{14}=sprintf('w%s_sess%d_striatum_dist_cross_exe2.nii',subj_name{s},ss);
                oP=spm_vol(char(O));
                
                % 3) extract betas / other images in each ROI r
                for r = roi % for each region
                    % get raw data for voxels in region
                    if size(R{r}.data,1)<3 % if less than 3 voxels
                        fprintf('Skipping network-%d - not enough voxels\n',r);
                        S.betaUW                = NaN;
                        S.betaRAW               = NaN;
                        S.resMS                 = NaN;
                        S.psc_train_1st         = NaN;
                        S.psc_untrain_1st       = NaN;
                        S.psc_train_2nd         = NaN;
                        S.psc_untrain_2nd       = NaN;
                        S.Tmap_train_1st        = NaN;
                        S.Tmap_train_2nd        = NaN;
                        S.Tmap_untrain_1st      = NaN;
                        S.Tmap_untrain_2nd      = NaN;
                        S.search_train_1st      = NaN;
                        S.search_train_2nd      = NaN;
                        S.search_untrain_1st    = NaN;
                        S.search_untrain_2nd    = NaN;
                        S.search_cross_1st      = NaN;
                        S.search_cross_2nd      = NaN;
                    else
                        % if more than 3 voxels
                        Y = region_getdata(betas,R{r});     % extract betas in region r
                        resMS = region_getdata(res,R{r});   % extract residuals in region r
                        data = region_getdata(oP,R{r});     % extract data from other contrasts
                        % construct region data
                        S.betaUW                  = {bsxfun(@rdivide,Y,ssqrt(resMS))}; % univariate pw
                        S.betaRAW                 = {Y};
                        S.resMS                   = {resMS};              
                        % info from maps for surface
                        S.psc_train_1st         = {data(1,:)};
                        S.psc_untrain_1st       = {data(2,:)};
                        S.psc_train_2nd         = {data(3,:)};
                        S.psc_untrain_2nd       = {data(4,:)};
                        S.Tmap_train_1st        = {data(5,:)};
                        S.Tmap_train_2nd        = {data(6,:)};
                        S.Tmap_untrain_1st      = {data(7,:)};
                        S.Tmap_untrain_2nd      = {data(8,:)};
                        S.search_train_1st      = {data(9,:)};
                        S.search_train_2nd      = {data(10,:)};
                        S.search_untrain_1st    = {data(11,:)};
                        S.search_untrain_2nd    = {data(12,:)};
                        S.search_cross_1st      = {data(13,:)};
                        S.search_cross_2nd      = {data(14,:)}; 
                        fprintf('Done network-%d.\n',r)
                    end
                    S.SN                      = s;
                    S.region                  = r;
                    T = addstruct(T,S);
                end
            end
            % save T
            save(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss)),'-struct','T');
            fprintf('\n');
        end       
    case 'ROI_stats_RepSup'
        sessN = [1:4];
        sn  = [4:9,11:25];
        parcelType = 'BG-striatum'; % BG, BG-striatum, BG-striatum-networks, group_MNI_striatum_buckner (only uniPW)
        betaChoice = 'multi'; % uni, multi or raw
        type='new';
        vararginoptions(varargin,{'sn','sessN','betaChoice','type','parcelType','region'});
        
        for ss=sessN
            
            T = load(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss))); % loads region data (T)
            % output structures
            switch(type)
                case 'new'
                    To=[];
                case 'add'
                    To=load(fullfile(BGstatsDir,sprintf('stats_FoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)));
            end  
            % do stats
            for s = sn % for each subject
                D = load(fullfile(glmFoSExDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                fprintf('\nSubject: %d session: %d\n',s, ss)
                num_run = numruns_task_sess;
                
                 roi=unique(T.region)';  % for each region
                 switch(parcelType) % for determining region side and type
                     case 'BG'
                         regSide=[1 1 1 2 2 2];
                         regType=[1 2 3 1 2 3];
                     case 'BG-striatum'
                         regSide=[1 1 2 2];
                         regType=[1 2 1 2];
                     case 'group_MNI_striatum_buckner'
                         regSide = [ones(1,7) ones(1,7)*2];
                         regType = [1:7 1:7];
                 end
                for r = roi % for each region
                    S = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                    fprintf('%d.',r)
                    
                    for exe = 1:2   % FoSEx
                        % check if there is data
                        if size(S.betaUW{:},1)>1
                            switch (betaChoice)
                                case 'uni'
                                    betaW  = S.betaUW{1}(D.FoSEx==exe,:);
                                case 'multi'
                                    betaW  = S.betaW{1}(D.FoSEx==exe,:);
                                case 'raw'
                                    betaW  = S.betaRAW{1}(D.FoSEx==exe,:);
                            end
                            %squeeze dimension of betaW if voxels with no
                            %data
                            indx=1;
                            for v=1:size(betaW,2)
                                if sum(isnan(betaW(:,v)))==size(betaW,1)% no data
                                    %sum(betaW(:,v))==0 || sum(isnan(betaW(:,v)))==size(betaW,1)% no data
                                else
                                    tmp(:,indx)=betaW(:,v);
                                    indx=indx+1;
                                end
                            end
                            betaW=tmp; % with no NaNs
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
                            So.psc_train_1st        = nanmean(S.psc_train_1st{:});
                            So.psc_train_2nd        = nanmean(S.psc_train_2nd{:});
                            So.psc_untrain_1st      = nanmean(S.psc_untrain_1st{:});
                            So.psc_untrain_2nd      = nanmean(S.psc_untrain_2nd{:});
                            So.Tmap_train_1st       = nanmean(S.Tmap_train_1st{:});
                            So.Tmap_train_2nd       = nanmean(S.Tmap_train_2nd{:});
                            So.Tmap_untrain_1st     = nanmean(S.Tmap_untrain_1st{:});
                            So.Tmap_untrain_2nd     = nanmean(S.Tmap_untrain_2nd{:});
                            So.search_train_1st     = nanmean(S.search_train_1st{:});
                            So.search_train_2nd     = nanmean(S.search_train_2nd{:});
                            So.search_untrain_1st   = nanmean(S.search_untrain_1st{:});
                            So.search_untrain_2nd   = nanmean(S.search_untrain_2nd{:});
                            So.search_cross_1st     = nanmean(S.search_cross_1st{:});
                            So.search_cross_2nd     = nanmean(S.search_cross_2nd{:});
                            % add vector length
                            vecLength = diag(G);
                            So.act_vecLength_train = nanmean(vecLength(1:6));
                            So.act_vecLength_untrain = nanmean(vecLength(7:12));
                            % indexing fields
                            So.SN       = s;
                            So.region   = r;
                            if sum(strcmp(parcelType,{'BG-striatum-networks','striatum_buckner'}))==0
                                if strcmp(parcelType,'striatum_buckner_hemi')
                                    if r<8
                                        So.regSide = 1;
                                        So.regType = r;
                                    else
                                        So.regSide = 2;
                                        So.regType = r-7;
                                    end
                                else
                                    So.regSide  = regSide(r);
                                    So.regType  = regType(r);
                                end
                            end
                            So.FoSEx    = exe;
                            To          = addstruct(To,So);
                        end
                    end; % FoSEx
                end; % each region
            end; % each subject
            
            % % save
            save(fullfile(BGstatsDir,sprintf('stats_FoSEx_%s_%sPW_sess%d.mat',parcelType,betaChoice,ss)),'-struct','To');
            fprintf('\nDone.\n')
        end
    case 'RepSupImage_toGroupMNI'
        sessN=[1:4];
        sn=[4:9,11:22];
        seqType={'TrainSeq_1st','UntrainSeq_1st','TrainSeq_2nd','UntrainSeq_2nd'};
        imageType='psc';
        glm=glmFoSExDir;
        
        vararginoptions(varargin,{'sn','sessN','imageType'});
        
        % first send over to transform subject images to MNI
        sml1_imana_BG_new('FUNC_copyFile_nativeToMNI','sn',sn,'seqType',seqType,'sessN',sessN,'imageType',imageType,'glm',glm);
        
        % make subject average in MNI space
        sml1_imana_BG_new('FUNC_subjAvrg','sn',sn,'seqType',seqType,'sessN',sessN,'imageType',imageType);
        
        % mask just with striatal mask
        sml1_imana_BG_new('FUNC_maskImages','seqType',seqType,'sessN',sessN,'imageType',imageType);
        
    case 'CALC_repsup' % ----------- save psc / dist for FoSEx - for cortex in sml1_imana_repsup   
        sn = [4:9,11:25];
        sessN = 1:4;
        parcelType = 'striatum_buckner_hemi';
        betaChoice = 'multiPW'; % uniPW or multiPW

        vararginoptions(varargin,{'sn','roi','seq','sessN','betaChoice','parcelType'});

        Stats = [];
        
        for ss = sessN % do per session number
            D = load(fullfile(BGstatsDir,sprintf('stats_FoSEx_%s_%s_sess%d.mat',parcelType,betaChoice,ss))); % loads region data (D)
            T = load(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss))); % loads region data (T)
            roi=unique(D.region)';
            
            runs=1:numruns_task_sess;
            conditionVec = kron(ones(numel(runs),1),[1:12]');       
            indx_train = 1:6;
            indx_untrain = 7:12;
            
            for s=1:length(sn)
                SI = load(fullfile(glmFoSExDir{ss}, subj_name{sn(s)}, 'SPM_info.mat'));   % load subject's trial structure
                for r=roi
                    for exe=1:2;    % FoSEx
                       D_exe = getrow(D,D.FoSEx==exe & D.SN==sn(s) & D.region==r); 
                       T2 = getrow(T,T.SN==sn(s) & T.region==r);
                        switch (betaChoice)
                            case 'uniPW'
                                beta = T2.betaUW{:}(SI.FoSEx==exe,:);
                            case 'multiPW'
                                beta = T2.betaW{:}(SI.FoSEx==exe,:);
                            case 'raw'
                                beta = T2.betaRAW(SI.FoSEx==exe,:);
                        end
                        
                        clear C;
                        for d = 1:6 %sequences
                            C.beta_seq_train(d,:)=mean(beta(conditionVec==indx_train(d),:),1);  % beta values for each digit (avrg across blocks)
                            C.beta_seq_untrain(d,:)=mean(beta(conditionVec==indx_untrain(d),:),1);
                        end
                        
                        AllDist = ssqrt(rsa_squareRDM(D_exe.RDM));
                        SeqTrain = triu(AllDist(indx_train,indx_train));
                        SeqUntrain = triu(AllDist(indx_untrain,indx_untrain));
                        SeqCross = triu(AllDist(indx_train,indx_untrain));
                        SeqTrainAll = SeqTrain(SeqTrain~=0);
                        SeqUntrainAll = SeqUntrain(SeqUntrain~=0);
                        SeqCrossAll = SeqCross(SeqCross~=0);
                        
                        S.dist_train=mean(SeqTrainAll);
                        S.dist_untrain=mean(SeqUntrainAll);
                        S.dist_cross=mean(SeqCrossAll);
                        S.beta_train=mean(mean(C.beta_seq_train));
                        S.beta_untrain=mean(mean(C.beta_seq_untrain));
                        S.act_vecLength_train = D_exe.act_vecLength_train;
                        S.act_vecLength_untrain = D_exe.act_vecLength_untrain;
                        if exe==1
                            S.psc_train=D_exe.psc_train_1st;
                            S.psc_untrain=D_exe.psc_untrain_1st;
                        else
                            S.psc_train=D_exe.psc_train_2nd;
                            S.psc_untrain=D_exe.psc_untrain_2nd;
                        end
                        
                        S.sn=sn(s);
                        S.roi=r;
                        if ~strcmp(parcelType,{'BG-striatum-networks','group_MNI_striatum_buckner'})
                            if r>numel(roi)/2
                                S.regSide=[2];
                                S.regType=r-(numel(roi)/2);
                            else
                                S.regSide=[1];
                                S.regType=r;
                            end
                            S.regType=rem(S.roi,numel(roi)+1);
                        end
                        if any(strcmp(parcelType,{'group_MNI_striatum_buckner','striatum_buckner_hemi'}))
                            % add contrast images
                            if exe==1
                                S.search_train      = D_exe.search_train_1st;
                                S.search_untrain    = D_exe.search_untrain_1st;
                                S.search_cross      = D_exe.search_cross_1st;
                            else
                                S.search_train      = D_exe.search_train_2nd;
                                S.search_untrain    = D_exe.search_untrain_2nd;
                                S.search_cross      = D_exe.search_cross_2nd;
                            end
                           % if r>7
                           %     S.regSide=2;
                           %     S.regType=r;
                           % else
                           %     S.regSide=1;
                           %     S.regType=r-7;
                           % end
                        end
                        S.sessN=ss;
                        S.FoSEx=exe;
                        Stats=addstruct(Stats,S);
                        fprintf('Done sess-%d %s exe-%d reg-%d\n',ss,subj_name{sn(s)},exe,r);
                    end; % FoSEx
                end; % roi
            end; % sn
        end
       % % save
        save(fullfile(BGstatsDir,sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice)),'-struct','Stats');    
    case 'PLOT_repsup_psc'
        parcelType='BG-striatum'; % striatum_buckner
        betaChoice = 'multiPW';
        roi=[1:4];
        vararginoptions(varargin,{'betaChoice','parcelType','roi'});
        
        S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice)));
        [regLab, hemiLab] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);
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
                title(sprintf('%s-%s',regLab{f},hemiLab{f}));
                %  plt.match('y');
                if f==1
                    ylabel(sprintf('Psc %s',SeqType{st}));
                    xlabel('Session');
                else
                    ylabel('');
                end
               % title(sprintf('%s',regname_BG{rem(f,4)}));
            end
        end
    case 'PLOT_repsup_dist'
       parcelType='striatum_buckner_hemi'; % striatum_buckner
       betaChoice = 'uniPW';
       seqType={'trained','untrained'};
       %roi=[1:3,5:7];
       roi=[2,4:7];
       distType = 'dist'; % dist or search
       vararginoptions(varargin,{'betaChoice','parcelType','roi','distType'});
       
       S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice))); 
       [regLab, hemiLab] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);
       for st=1:2 % seqType
           if st==1
                style=styTrained_exe;
                var=sprintf('S.%s_train',distType);
            else
                style=styUntrained_exe;
                var=sprintf('S.%s_untrain',distType);
            end
           figure
           sub=0;
           for f = roi
               sub=sub+1;
               subplot(1,numel(roi),sub);
               plt.line([S.sessN>3 S.sessN],eval(var),'split',S.FoSEx,'subset',S.roi==f,'style',style,'leg',{'1st','2nd'});
               drawline(0,'dir','horz');
               %     plt.match('y');
               if f==1
                   ylabel(sprintf('Distances %s',seqType{st})); xlabel('Session');
               else
                   ylabel('');
               end
               title(sprintf('%s-%s',regLab{f},hemiLab{f}));
           end
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
            title(sprintf('%s-%s',regLab{f},hemiLab{f}));
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
       
    case 'CALC_percent_repsup' % ------- calculate amount of RS 
        parcelType='BG-striatum';
        betaChoice = 'multiPW'; % or uniPW
        roi=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','parcelType'});
        
        S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice)));
        
        A = getrow(S,S.FoSEx==1);
        
        P.psc_beta = [S.beta_train(S.FoSEx==2)./S.beta_train(S.FoSEx==1);S.beta_untrain(S.FoSEx==2)./S.beta_untrain(S.FoSEx==1)];
        P.psc_psc  = [S.psc_train(S.FoSEx==2)./S.psc_train(S.FoSEx==1);S.psc_untrain(S.FoSEx==2)./S.psc_untrain(S.FoSEx==1)];
        P.psc_dist = [S.dist_train(S.FoSEx==2)./S.dist_train(S.FoSEx==1);S.dist_untrain(S.FoSEx==2)./S.dist_untrain(S.FoSEx==1)];
        
        P.sn = [A.sn; A.sn];
        P.roi = [A.roi; A.roi];
        P.seq = [ones(size(A.sn)); ones(size(A.sn))*2];  % 1 for trained, 2 for untrained
        P.sessN = [A.sessN; A.sessN];
        
        % % save
        save(fullfile(BGstatsDir,sprintf('RepSup_%s_percent_%s.mat',parcelType,betaChoice)),'-struct','P');
    case 'PLOT_percent_repsup_psc'
        
        betaChoice = 'multiPW';
        roi=[1:4];
        parcelType='BG-striatum';
        vararginoptions(varargin,{'betaChoice','roi','parcelType'});
        
        S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_percent_%s.mat',parcelType,betaChoice)));
        
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
         %   title(sprintf('%s',regname_BG{rem(f,4)}));
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
       parcelType = 'BG-striatum';
       vararginoptions(varargin,{'betaChoice','roi','parcelType'});
       
       S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_stats_%s.mat',parcelType,betaChoice))); 
 
       A = getrow(S,S.FoSEx==1);
       
       P.subtr_beta = [S.beta_train(S.FoSEx==1)-S.beta_train(S.FoSEx==2);S.beta_untrain(S.FoSEx==1)-S.beta_untrain(S.FoSEx==2)];
       P.subtr_psc = [S.psc_train(S.FoSEx==1)-S.psc_train(S.FoSEx==2);S.psc_untrain(S.FoSEx==1)-S.psc_untrain(S.FoSEx==2)];
       P.subtr_dist = [S.dist_train(S.FoSEx==1)-S.dist_train(S.FoSEx==2);S.dist_untrain(S.FoSEx==1)-S.dist_untrain(S.FoSEx==2)];
       P.sn = [A.sn; A.sn];
       P.roi = [A.roi; A.roi];
       P.seq = [ones(size(A.sn)); ones(size(A.sn))*2];  % 1 for trained, 2 for untrained
       P.sessN = [A.sessN; A.sessN];
       
       % % save
        save(fullfile(BGstatsDir,sprintf('RepSup_%s_subtract_%s.mat',parcelType,betaChoice)),'-struct','P');    
    case 'PLOT_subtract_repsup_psc'
        betaChoice = 'multiPW';
        roi=[1:4];
        parcelType='BG-striatum';
        vararginoptions(varargin,{'betaChoice','roi','parcelType'});
        
        S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_subtract_%s.mat',parcelType,betaChoice)));
        [regLab, hemiLab] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);
        figure
        sub=0;
        for f = roi
            sub=sub+1;
            subplot(1,numel(roi),sub);
            plt.line(S.sessN,S.subtr_psc,'split',S.seq,'subset',S.roi==f,'style',stySeq,'leg',{'Trained','Untrained'});
            plt.match('y');
            title(sprintf('%s-%s',regLab{f},hemiLab{f}));
            drawline(0,'dir','horz');
            if f==1
                ylabel('Subtract 1st - 2nd activation')
            else
                ylabel('')
            end
        end
    case 'PLOT_subtract_repsup_dist'
        betaChoice = 'multiPW';
        parcelType = 'BG-striatum';
        roi=[1:4];
        vararginoptions(varargin,{'betaChoice','roi','parcelType'});
        
        S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_subtract_%s.mat',parcelType,betaChoice)));

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
       %    title(sprintf('%s',regname_BG{rem(f,4)}))
       end
    case 'PLOT_subtract_psc_sessions'
        sessN=[1:4];
        seqType='trained';
        betaChoice='multiPW';
        parcelType='BG-striatum';
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice','parcelType'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        S = load(fullfile(BGstatsDir,sprintf('RepSup_%s_subtract_%s.mat',parcelType,betaChoice)));
        T=getrow(S,ismember(S.sessN,sessN));
        
        T.hemi=T.roi;
        T.hemi(T.roi<3,:)=1;
        T.hemi(T.roi>2,:)=2;
        
        switch(seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end
        
        figure
        for h=1:2
            subplot(1,2,h)
            plt.line(T.roi,T.subtr_psc,'split',T.sessN,'subset',T.seq==st&T.hemi==h,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
            drawline(0,'dir','horz');
            ylabel('PSC repetition suppression - 1st-2nd');
            set(gca,'XTickLabel',regname_BG(1:2));
            title(sprintf('%s %s',seqType,hemName{h}));
            xlabel('ROI');
            plt.match('y');
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
        reg = [1:4];
        sn  = [4:9,11:25];
        sessN = [1:4];
        parcelType='BG-striatum';
        betaChoice = 'multiPW'; 
        subtract_mean=0; % do NOT subtract mean - it distorts the pattern
        vararginoptions(varargin,{'sn','reg','sessN','betaChoice','parcelType','subtract_mean'});
        SAll = [];
        STAll= [];
        
        for  ss = sessN
            D   = load(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss)));
            for roi = reg;
                for s = sn;
                    SI = load(fullfile(glmFoSExDir{ss}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
                    
                    for exe = 1:2
                        t   = getrow(D,D.region==roi & D.SN==s);
                        % check if there is any data
                        if size(t.SN,1)~=0
                            SI_exe = getrow(SI,SI.FoSEx==exe);
                            switch betaChoice
                                case 'multiPW'
                                    data = t.betaW{1}(SI.FoSEx==exe,:);
                                case 'uniPW'
                                    data = t.betaUW{1}(SI.FoSEx==exe,:);
                                case 'raw'
                                    data = t.betaRAW{1}(SI.FoSEx==exe,:);
                            end
                            %squeeze dimension of betaW if voxels with no data
                            indx=1;
                            for v=1:size(data,2)
                                if sum(isnan(data(:,v)))==size(data,1)% no data
                                    %sum(betaW(:,v))==0 || sum(isnan(betaW(:,v)))==size(betaW,1)% no data
                                else
                                    tmp(:,indx)=data(:,v);
                                    indx=indx+1;
                                end
                            end
                            data=tmp; % with no NaNs
                            
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
                            [G,Sig]     = pcm_estGCrossval(data(1:(12*max(numRuns)),:),SI_exe.run,SI_exe.seqNumb);
                            indM=indicatorMatrix('identity',SI_exe.seqNumb);
                            % meanP=pinv(indM)*data;
                            % G_ncv=cov(meanP');
                            % euc = rsa_squareRDM(pdist(meanP,'euclidean'));
                            dist = rsa_squareRDM(rsa.distanceLDC(data,SI_exe.run,SI_exe.seqNumb));
                            C=corr_crossval(G,'reg','minvalue');
                            C=rsa_squareRDM(C);
                            
                            % average trained dist
                            S.corrDist(1,:) = sum(sum(triu(C(1:6,1:6))))/(6*5/2);
                            % average untrained dist
                            S.corrDist(2,:) = sum(sum(triu(C(7:12,7:12))))/(6*5/2);
                            %  S.dist_ncv(1,:) = sum(sum(triu(euc(1:6,1:6))))/(6*5/2);
                            %  S.dist_ncv(2,:) = sum(sum(triu(euc(7:12,7:12))))/(6*5/2);
                            S.dist(1,:)     = sum(sum(triu(dist(1:6,1:6))))/(6*5/2);
                            S.dist(2,:)     = sum(sum(triu(dist(7:12,7:12))))/(6*5/2);
                            S.seqType=[1;2]; % trained, untrained
                            S.sessN=[ss;ss];
                            S.roi=[roi;roi];
                            S.sn=[s;s];
                            S.FoSEx=[exe;exe];
                            SAll=addstruct(SAll,S);
                            
                            ST.corr_seqType = sum(sum(triu(C(7:12,1:6))))/(6*5);
                            ST.sessN=ss;
                            ST.roi=roi;
                            ST.sn=s;
                            ST.FoSEx=exe;
                            STAll=addstruct(STAll,ST);
                        end
                    end
                end
            end
        end
        save(fullfile(BGstatsDir,sprintf('corrDist_RS_%s_%s.mat',betaChoice,parcelType)),'-struct','SAll');
        save(fullfile(BGstatsDir,sprintf('corrDist_seqType_RS_%s_%s.mat',betaChoice,parcelType)),'-struct','STAll');
    case 'PLOT_corrDist_RS'
        roi=[1:4];
        betaChoice = 'multiPW';
        parcelType='BG-striatum';
        vararginoptions(varargin,{'roi','betaChoice','parcelType'});
        
        S=load(fullfile(BGstatsDir,sprintf('corrDist_RS_%s_%s.mat',betaChoice,parcelType)));
        
        SeqType={'trained','untrained'};
        [regLab, hemiLab] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks);
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
               % plt.line([S.sessN>3 S.sessN],S.dist,'style',style,'split',S.FoSEx,'subset',S.seqType==seq&S.roi==r,'leg',{'1st','2nd'});
                plt.line([S.sessN>3 S.sessN],S.corrDist,'style',style,'split',S.FoSEx,'subset',S.seqType==seq&S.roi==r,'leg',{'1st','2nd'});
              %  plt.match('y');
                drawline(0,'dir','horz');
                title(sprintf('%s-%s',regLab{r},hemiLab{r}));
                if r==1
                    xlabel('Session'); ylabel(sprintf('Corr distance %s',SeqType{seq}));
                else
                    ylabel('');
                end
                indx=indx+1;
            end
        end   
    case 'PLOT_corrDist_RS_seqType'
        roi=[1:4];
        betaChoice = 'multi';
        parcelType = 'BG-striatum';
        vararginoptions(varargin,{'roi','betaChoice','parcelType'});
        
        S=load(fullfile(BGstatsDir,sprintf('corrDist_seqType_RS_%sPW_%s.mat',betaChoice,parcelType)));
        
        Exe={'1st exe','2nd exe'};
        
        figure;
        indx=1;
        for r=roi
            subplot(1,numel(roi),indx)
            plt.line([S.sessN>3 S.sessN],S.corr_seqType,'style',stySeqType_exe,'split',S.FoSEx,'subset',S.roi==r,'leg',{'1st','2nd'});
            drawline(0,'dir','horz');
            plt.match('y');
            if r==1
                xlabel('Session'); ylabel('Corr distance between seqType');
            else
                ylabel('');
            end
            indx=indx+1;
        end
    case 'CALC_corrDist_subtract'
       roi=[1:4];
       betaChoice='multiPW';
       parcelType='BG-striatum';
       vararginoptions(varargin,{'roi','betaChoice','parcelType'});
       
       S=load(fullfile(BGstatsDir,sprintf('corrDist_RS_%s_%s.mat',betaChoice,parcelType)));
 
       A = getrow(S,S.FoSEx==1);
       
       C.sn=A.sn;
       C.roi=A.roi;
       C.corrDist = [S.corrDist(S.FoSEx==1)-S.corrDist(S.FoSEx==2)];
       C.sessN=A.sessN;
       C.seqType=A.seqType;
       
       % % save
       save(fullfile(BGstatsDir,sprintf('RepSup_subtract_corrDist_%s_%s.mat',betaChoice,parcelType)),'-struct','C');
    case 'PLOT_corrDist_subtract_sess'
        sessN=[1:4];
        seqType='trained';
        betaChoice = 'uniPW';
        parcelType = 'BG-striatum';
        vararginoptions(varargin,{'sessN','regExcl','seqType','betaChoice'});
        sessLegend={'sess1','sess2','sess3','sess4'};
        
        T=load(fullfile(BGstatsDir,sprintf('RepSup_subtract_corrDist_%s_%s.mat',betaChoice,parcelType))); 
        T=getrow(T,ismember(T.sessN,sessN));
        
        T.hemi=T.roi;
        T.hemi(T.roi<3,:)=1;
        T.hemi(T.roi>2,:)=2;
        
        switch(seqType)
            case 'trained'
                st=1;
            case 'untrained'
                st=2;
        end       
        figure
        for h=1:2
            subplot(1,2,h)
            plt.line(T.roi,T.corrDist,'split',T.sessN,'subset',T.seqType==st&T.hemi==h,'leg',sessLegend(sessN),'leglocation','northeast','style',stySess);
            drawline(0,'dir','horz');
            ylabel('Distance repetition suppression');
            set(gca,'XTickLabel',regname_BG(1:3));
            title(sprintf('%s %s',seqType,hemName{h}));
            xlabel('ROI');
            plt.match('y');
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
    
    case 'PCM_data_repsupModels' % new - repsup models
        runEffect  = 'fixed';
        beta_choice = 'mw';
        algorithm='NR'; % minimize or NR
        reg = [1:3,5:7];
        sn=[1:9,11:22];
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
            save(fullfile(BGDir,sprintf('PCM_repsup_reliability_%s_%s_sess%d_%s.mat',models{M_type},seqType,ss,runEffect)),'-struct','AllReg');
        end % session    
    case 'PLOT_pcm_corr_allSess_NEW'
        reg = [1:3,5:7];
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
           %     title(sprintf('%s',regname{r}));
                if sub==1
                    ylabel(sprintf('%s corr',corrType{ct}));
                    xlabel('Session');
                else
                    ylabel('');
                end
                sub=sub+1;
            end     
        end
       
    case 'getBetas_MNI'
        sn=[4:9,11:25];
        sessN=[1:4];
        maskName='striatum'; % striatum, CaudateN_lh/_rh, Putamen_lh/_rh
        vararginoptions(varargin,{'sn','sessN','maskName'});
        TT=[];
        % load the mask
        Vmask=spm_vol(fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_funcMASK_%s.nii',maskName)));
        maskFile= spm_read_vols(Vmask);
        % indicate which voxels to sample
        indx = maskFile(:);
        
        for ss=sessN
            for s=sn
                betaMNIDir  = fullfile(BGfuncDir,'MNI_betas',subj_name{s},sprintf('sess-%d',ss));
                VresMS      = spm_vol(fullfile(betaMNIDir,'wResMS.nii'));
                ResMS       = spm_read_vols(VresMS);
                % resample?
                % loop over all betas
                for i=1:104
                    nam{1}  = fullfile(betaMNIDir,sprintf('wbeta_%2.4d.nii',i));
                    Vbeta   = spm_vol(nam{1});
                    beta    = spm_read_vols(Vbeta);
                    % univariate normalisation
                    tmp     = beta./ssqrt(ResMS);
                    tmp2    = tmp(:);
                    bIncl   = tmp2(find(indx));
                    B1(i,:) = bIncl';
                end
                % prepare a structure
                T.betas = {B1};
                T.sn    = s;
                T.sessN = ss;  
                T.group = rem(s,2); % subject group: 1 - odd, 0 - even
                TT=addstruct(TT,T);
                fprintf('Done betas ses-%d %s\n',ss,subj_name{s});
            end
        end  
        save(fullfile(BGfuncDir,'MNI_betas','avrg',sprintf('groupBetas_%s',maskName)),'-struct','TT');
    case 'betas_corr'
        sessN=[1:4];
        maskName='striatum'; % striatum, CaudateN_lh/_rh, Putamen_lh/_rh
        vararginoptions(varargin,{'sn','sessN','maskName'});
        
        T          = load(fullfile(BGfuncDir,'MNI_betas','avrg',sprintf('groupBetas_%s',maskName))); 
        sn         = unique(T.sn);
        runs       = 1:numruns_task_sess;
        conds      = num_seq;
        seqT       = [ones(1,6) ones(1,6)*2];
        condVec    = repmat([conds],1,length(runs))';
        partVec    = kron([runs],ones(1,length(conds)))';
        seqTVec    = repmat([seqT],1,length(runs))';
        KK=[]; II=[];
        
        for ss=sessN
            % loop twice across subjects
            for st=1:2 % seqType
                R=zeros(length(sn));
                for i=1:length(sn);
                    for j=1:length(sn);
                        B1 = getrow(T,T.sn==sn(i) & T.sessN==ss);
                        B2 = getrow(T,T.sn==sn(j) & T.sessN==ss);
                        if (i==j) % within subject
                            cond = condVec.*(seqTVec==st);
                            if st==2
                                cond = cond-6; % make it 1-6 for untrained
                            end
                            data = B1.betas{1}(1:length(condVec),:);
                            X = indicatorMatrix('identity_p',cond);
                            mPat=[];
                            split=mod(partVec,2);
                            mPat(:,:,1)     = pinv(X(split==0,:))*data(split==0,:);
                            mPat(:,:,2)     = pinv(X(split==1,:))*data(split==1,:);
                            % indx - find only voxels which don't have a nan
                            indx            = find(~isnan(mPat(1,:,1))&~isnan(mPat(1,:,2)));
                            meanPat         = mean(mPat,1);             % Mean activity pattern over sequences
                            K.corr_mean     = corr(meanPat(:,indx,1)',meanPat(:,indx,2)');
                        else % across
                            cond = [condVec.*(seqTVec==st); condVec.*(seqTVec==st)];
                            data    = [{B1.betas{1}(1:length(condVec),:)} {B2.betas{1}(1:length(condVec),:)}];
                            
                            ses   = [ones(size(partVec));ones(size(partVec))*2];
                            X     = indicatorMatrix('identity_p',cond);
                            mPat=[];
                            
                            % Do split half correlations
                            split = mod([partVec;partVec],2);
                            mPat(:,:,1) = pinv(X(ses==1 & split==0,:)) * data{1}(split(ses==1)==0,:);
                            mPat(:,:,2) = pinv(X(ses==1 & split==1,:)) * data{1}(split(ses==1)==1,:);
                            mPat(:,:,3) = pinv(X(ses==2 & split==0,:)) * data{2}(split(ses==2)==0,:);
                            mPat(:,:,4) = pinv(X(ses==2 & split==1,:)) * data{2}(split(ses==2)==1,:);
                            % indx - find only voxels which don't have a nan (across all 4)
                            indx = find(~isnan(mPat(1,:,1))&~isnan(mPat(1,:,2))&~isnan(mPat(1,:,3))&~isnan(mPat(1,:,4)));
                            meanPat  = mean(mPat,1);  % Mean activity pattern over sequences
                            COR = [];
                            COR(1) = corr(meanPat(:,indx,1)',meanPat(:,indx,3)');
                            COR(2) = corr(meanPat(:,indx,1)',meanPat(:,indx,4)');
                            COR(3) = corr(meanPat(:,indx,2)',meanPat(:,indx,3)');
                            COR(4) = corr(meanPat(:,indx,2)',meanPat(:,indx,4)');
                            K.corr_mean = fisherinv(mean(fisherz(COR)));
                        end
                        % add other info
                        K.sn1       = sn(i);
                        K.sn2       = sn(j);
                        K.sn_idx1   = i;
                        K.sn_idx2   = j;
                        K.seqType   = st;
                        K.sessN     = ss;
                        KK = addstruct(KK,K);
                        % making matrix
                        R(j,i) = K.corr_mean;
                    end
                end
                % sequence type st, session ss complete
                % extract the R - check
                I.corrMatrix = {rsa_vectorizeIPMfull(R)};
                I.seqType    = st;
                I.sessN      = ss;
                II=addstruct(II,I);
            end
            fprintf('Done sess-%d \n',ss);
        end
        
        %save
        save(fullfile(BGstatsDir,sprintf('corr_betas_%s_MNI',maskName)),'-struct','KK');
        save(fullfile(BGstatsDir,sprintf('corr_betas_%s_MNI_matrix',maskName)),'-struct','II');
        
    case 'betas_corr_plot'
        metric='corr_mean'; % corr or corr_mean
        parcelType='striatum';
        vararginoptions(varargin,{'metric','parcelType'});
        
        T = load(fullfile(BGstatsDir,sprintf('corr_betas_%s_MNI',parcelType)));
   
        figure
        plt.line([T.sessN>3 T.sessN],T.(metric),'subset',T.sn1==T.sn2,'plotfcn','nanmean','leg','within-subj');
        hold on;
        plt.line([T.sessN>3 T.sessN],T.(metric),'subset',T.sn1~=T.sn2,'plotfcn','nanmean','style',stySeq,'leg',{'across-subj'},'leglocation','northeast');
        xlabel('Session');
        ylabel('Correlation');
        title(sprintf('%s pattern correlation',parcelType));
    case 'betas_corr_matrix'
        parcelType='striatum';
        sessN=[1:4];
        seqType=[1,2];
        vararginoptions(varargin,{'parcelType','sessN'});
        
        T = load(fullfile(BGstatsDir,sprintf('corr_betas_%s_MNI_matrix',parcelType)));
   
        figure
        indx=1;
        for ss=sessN
            for st=seqType
                g=getrow(T,T.sessN==ss & T.seqType==st);
                subplot(length(seqType),length(sessN),indx);
                imagesc(rsa_squareIPMfull(g.corrMatrix{:}));
                caxis([0 0.5]);
                indx=indx+1;
            end
        end
        
    case 'SURF_makePaint_striatum' 
    % Creates ROI boundaries on the group template/atlas (fsaverage).
    % ROIs are defined using:
    %depths = [-5:0.01:5]; % depth of sampling
    depths = [-0.1:0.01:0.9]; % depth of sampling
    for h=1:2 % loop over hemispheres
        % load the surface file
        S = gifti(fullfile(BGDir,'surface_new',hemName{h},sprintf('%s.striatum.surf.gii',hem{h})));
        [normVert normFace] = compute_normal(S.vertices,S.faces);
        normPoints = S.vertices.*normVert';
        
        % load parcelation
        file = spm_vol(fullfile(BGanatDir,'striatum_buckner','avrg',sprintf('striatum_buckner_%s.nii',hem{h})));
        %  file = spm_vol(fullfile(BGanatDir,'striatum_buckner','avrg','subjAvrg_anat_striatum_buckner.nii'));
        for d=1:length(depths)
            % sample new vertices at depth d
            vertNew = S.vertices+normPoints*depths(d);
            % transform your data in file to same coordinates
            [i,j,k] = spmj_affine_transform(vertNew(:,1),vertNew(:,2),vertNew(:,3),inv(file.mat));
            i = double(i);
            j = double(j);
            k = double(k);
            % sample the volume to vertices
            X(:,d) = spm_sample_vol(file,i,j,k,0);        
        end
        X=round(X);
        for i=1:size(X,1)
            dataROI(i,:) = mode(X(i,find(X(i,:))));
        end
        dataROI(isnan(dataROI))=0;
        
        % create and save ROIpaint file
        for i=1:7
           % roi_names{i} = sprintf('network-%d',i);
            names{i} = sprintf('network-%d',i);
            areas{i}{1} = {sprintf('network-%d',i)};
        end
        colors = [120 18 134;...
                  70 130 180;...
                  0 118 14;...
                  196 58 250;...
                  220 248 164;...
                  230 148 34;...
                  205 62 78];
        
        Paint = caret_struct('paint','data',dataROI,'paintnames',names,'column_name',{'network'});
        Paint.index = [1:length(dataROI)]';
        % save file
        newPaintPath = fullfile(BGDir,'surface','group',hemName{h},'7networks.paint');
        newAreaPath = fullfile(BGDir,'surface','group',hemName{h},'7networks.areacolor');
        caret_save(newPaintPath,Paint);
        caret_combinePaint(newPaintPath,newPaintPath,newAreaPath,...
            'areas',areas,'names',names,'colors',colors);
        clear i j k X dataROI;
    end;    
    case 'SURF_makeMetric_group' % old
        %depths = [-1:0.01:5]; % depth of sampling
        depths = [-0.1:0.01:0.9]; % depth of sampling
        sessN = [1:4];
        operation = 'nanmean';
        imageType = 'psc'; % psc, dist or corrDist
        seqType = {'TrainSeq','UntrainSeq'}; 
        % psc: TrainSeq, UntrainSeq
        % dist / corrDist: trained, untrained, cross
        vararginoptions(varargin,{'depths','operation','sessN','imageType','seqType'});

        for h=1:2
            % load in the surface file
            S = gifti(fullfile(BGDir,'surface','group',hemName{h},sprintf('%s.striatum.surf.gii',hem{h})));
            [normVert normFace] = compute_normal(S.vertices,S.faces);
            normPoints = S.vertices.*normVert';
            for st=1:length(seqType)
                for ss=sessN
                    % load the region file, mask and the file of interest
                    switch imageType
                        case 'psc'
                            file = spm_vol(fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_%s_sess%d_%s_MASK_striatum_%s.nii',imageType,ss,seqType{st},hem{h})));
                        case 'dist' % LDC searchlight
                          %  file = spm_vol(fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_searchlight_striatum_sess%d_%s_%s_MASK_striatum_%s.nii',ss,imageType,seqType{st},hem{h})));
                             file = spm_vol(fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_searchlight_striatum_sess%d_dist_%s_MASK_striatum_%s.nii',ss,seqType{st},hem{h})));
                        case 'corrDist' % correlation distance searchlight
                             file = spm_vol(fullfile(BGfuncDir,'MNI','avrg',sprintf('subjAvrg_searchlight_striatum_sess%d_corrDist_%s_MASK_striatum_%s.nii',ss,seqType{st},hem{h})));
                    end
                    clear i j k X;
                    for d=1:length(depths)
                        % sample new vertices at depth d
                        vertNew = S.vertices+normPoints*depths(d);
                        % transform your data in file to same coordinates
                        [i,j,k] = spmj_affine_transform(vertNew(:,1),vertNew(:,2),vertNew(:,3),inv(file.mat));
                        i = double(i);
                        j = double(j);
                        k = double(k);
                        % sample the volume to vertices
                        X(:,d) = spm_sample_vol(file,i,j,k,1);
                    end    
                    % apply operation across all depths
                    switch(operation)
                        case 'nanmean'
                            dataNew = nanmean(X,2);
                            if strcmp(imageType,'dist') || strcmp(imageType,'corrDist')
                                dataNew = ssqrt(dataNew);
                            end
                        case 'nanmax'
                            dataNew = nanmax(X,2);
                    end

                    columnName = sprintf('%s_sess%d_%s_%s',seqType{st},ss,imageType,operation);
                    M=caret_struct('metric','data',dataNew,'column_name',{columnName});
                    % save file
                    newFilePath = fullfile(BGDir,'surface_new',hemName{h},sprintf('%s.%s_%s_sess%d_%s.metric',hem{h},seqType{st},imageType,ss,operation));
                    caret_save(newFilePath,M);
                    fprintf('Done %s sess-%d: %s-%s\n',hem{h},ss,imageType,seqType{st});
                end
            end
        end   
    case 'SURF_summaryMetric_group' %old
        % Makes a summary metric containing all trained / untrained
        % contrasts
        % for psc / searchlights: dist, corrDist - RS or overall
        sessN = [1:4];
        %seqType = {'TrainSeq_1st','TrainSeq_2nd','UntrainSeq_1st','UntrainSeq_2nd'};
        seqType = {'trained_exe1','trained_exe2','untrained_exe1','untrained_exe2','cross_exe1','cross_exe2'}; % for searchlight
        operation = 'nanmean';
        imageType = 'dist'; % psc, dist, corrDist
        vararginoptions(varargin,{'sessN','sn','seqType','imageType'});
        
        if length(seqType)>3
            outType = strcat('RS_',imageType);
        else
            outType = imageType;
        end    
        for ss=sessN
            SummaryName = sprintf('.summary_%s_sess%d.metric',outType,ss);
            for h=1:2
                clear data column_name
                surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
                cd(surfaceGroupDir);
                %----get the full directory name of the metric files
                for i=1:length(seqType);      
                    filenames=fullfile(surfaceGroupDir,sprintf('%s.%s_%s_sess%d_%s.metric',hem{h},seqType{i},imageType,ss,operation));                  
                    Data=caret_load(filenames);
                    data(:,i) = Data.data;    
                    column_name{i}=fullfile(sprintf('%s_sess%d_%s',outType,ss,seqType{i}));
                end;
                % create metric structure
                C = caret_struct('metric','data',data,'column_name',column_name);
                caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);         
            end;
            fprintf('Done sess-%d\n',ss);
        end
    case 'SURF_RS_RGB_group' %old
        % plotting RS on surface as absolute
        % subtracting 2nd - 1st activation - ABSOLUTE RS
        sessN = [1:4];
        imageType = 'RS_psc'; % RS_psc, RS_dist
        type = 'suppression'; % suppression or enhancement 
        vararginoptions(varargin,{'sessN','sessN','imageType','maxLim','type'});
        for h=1:2
            surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);      
            for ss=sessN
                clear RGBdata;
                C = caret_load([hem{h} sprintf('.summary_%s_sess%d.metric',imageType,ss)]);
                % calculate difference 1st - 2nd exe
                minLim=0;
                % maxLim
                RGBdata(:,1) = C.data(:,1)-C.data(:,2); % Red: trained 1st-2nd BOLD signal
                RGBdata(:,3) = C.data(:,3)-C.data(:,4); % Blue: untrained 1st-2nd BOLD signal
                switch imageType
                    case 'RS_psc'
                        maxLim = 0.05;
                    case 'RS_dist'
                        maxLim = 0.01;
                        %RGBdata = ssqrt(RGBdata);
                end
                switch type
                    case 'suppression'
                        RGBdata(RGBdata(:,1)<minLim,1) = 0;
                        RGBdata(RGBdata(:,3)<minLim,3) = 0;
                    case 'enhancement'
                        RGBdata(RGBdata(:,1)>minLim,1) = 0;
                        RGBdata(RGBdata(:,3)>minLim,3) = 0;
                        RGBdata = RGBdata * (-1)';
                       % sc=[-maxLim minLim;-maxLim minLim;-maxLim minLim];  % scaling
                end
                sc=[minLim maxLim;minLim maxLim;minLim maxLim];  % scaling
                name={sprintf('%s_absolute_%s_sess%d',imageType,type,ss)};
                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_%s_sess%d.RGB_paint',imageType,type,ss)],C);
            end       
        end;
    case 'SURF_RS_diff_metric'
        % plotting RS on surface as absolute
        % subtracting 2nd - 1st activation - ABSOLUTE RS
        sessN = [1:4];
        imageType = 'psc_RS'; % RS_psc, RS_dist
        seqType = {'trained','untrained'};
        vararginoptions(varargin,{'sessN','sessN','imageType','maxLim','type'});
        for h=1:2
            surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);      
            for ss=sessN
                C = caret_load([hem{h} sprintf('.summary_%s_sess%d.metric',imageType,ss)]);
                % calculate difference 1st - 2nd exe
                for st = 1:2 % trained / untrained
                    data(:,st) = C.data(:,(st-1)*2+2)-C.data(:,(st-1)*2+1); % 2nd-1st
                    colIndx{st} = sprintf('RS_diff_%s',seqType{st});
                end
                C2=caret_struct('metric','data',data,'column_name',colIndx);
                C2.index=C.index;
                C2.column_color_mapping=C.column_color_mapping(1:2,:);
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_DIFF_sess%d.metric',imageType,ss)],C2);
                clear data;
            end       
        end;
    case 'SURF_makeMetric'
         %depths = [-1:0.01:5]; % depth of sampling
        depths = [-0.1:0.01:0.9]; % depth of sampling
        sessN = [1:4];
        imageType = 'psc'; % psc, dist or corrDist
        seqType = {'TrainSeq','UntrainSeq'}; 
        sn=[4:9,11:25];
        % psc: TrainSeq, UntrainSeq
        % dist / corrDist: trained, untrained, cross
        vararginoptions(varargin,{'depths','operation','sessN','imageType','seqType','sn'});

        for h=1:2
            % load in the surface file
            S = gifti(fullfile(BGDir,'surface','group',hemName{h},sprintf('%s.striatum.surf.gii',hem{h})));
            [normVert normFace] = compute_normal(S.vertices,S.faces);
            normPoints = S.vertices.*normVert';
            for st=1:length(seqType)
                for ss=sessN
                    for s=sn
                        % load the region file, mask and the file of interest
                        file = spm_vol(fullfile(BGfuncDir,'MNI',subj_name{s},sprintf('%s_sess%d_%s_%s_MASK_striatum_%s.nii',subj_name{s},ss,imageType,seqType{st},hem{h})));              
                        clear i j k X;
                        for d=1:length(depths)
                            % sample new vertices at depth d
                            vertNew = S.vertices+normPoints*depths(d);
                            % transform your data in file to same coordinates
                            [i,j,k] = spmj_affine_transform(vertNew(:,1),vertNew(:,2),vertNew(:,3),inv(file.mat));
                            i = double(i);
                            j = double(j);
                            k = double(k);
                            % sample the volume to vertices
                            X(:,d) = spm_sample_vol(file,i,j,k,1);
                        end
                        % apply operation across all depths          
                        dataNew = nanmean(X,2);
                        if strcmp(imageType,'dist') || strcmp(imageType,'corrDist')
                            dataNew = ssqrt(dataNew);
                        end             
                        columnName = sprintf('%s_%s_sess%d_%s_nanmean',subj_name{s},seqType{st},ss,imageType);
                        M=caret_struct('metric','data',dataNew,'column_name',{columnName});
                        % save file
                        subjDir = fullfile(BGDir,'surface',subj_name{s},hemName{h});
                        dircheck(subjDir);
                        if length(seqType)>3 % RS case - add in name
                            newFilePath = fullfile(subjDir,sprintf('%s.%s_RS_%s_sess%d_nanmean.metric',hem{h},imageType,seqType{st},ss));
                        else
                            newFilePath = fullfile(subjDir,sprintf('%s.%s_%s_sess%d_nanmean.metric',hem{h},imageType,seqType{st},ss));
                        end
                        caret_save(newFilePath,M);
                        fprintf('Done sess-%d-%s:\t%s: %s-%s\n',ss,hem{h},subj_name{s},imageType,seqType{st});
                    end
                end
            end
        end
    case 'SURF_makeGroup'
        sessN       = [1:4];
        sn          = [4:9,11:25];
        imageType   = 'dist'; % dist, corrDist, psc
        OUTname     = {'trained','untrained','cross'};
        vararginoptions(varargin,{'sessN','sn','imageType','OUTname','inputcol','replaceNaN'});

        inputcol   = 1:length(OUTname)*length(sessN);
        replaceNaN = ones(size(inputcol));
        
        % Loop over hemispheres.
        for h = 1:2
            % Go to the directory where the group surface atlas resides
            surfaceGroupDir = fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);      
            % Loop over each input metric file in 'OUTname' and make a group metric file
            indx=1; % re-set for each hemisphere
            for j = 1:length(OUTname);
                for  ss=sessN
                    % Loop over subjects...
                    for s = 1:length(sn);
                        % ...and define the names of their metric files
                        if length(OUTname)>3
                            infilenames{indx}{s} = fullfile(BGDir,'surface',subj_name{sn(s)},hemName{h},sprintf('%s.%s_RS_%s_sess%d_nanmean.metric',hem{h},imageType,OUTname{j},ss));
                        else
                            infilenames{indx}{s} = fullfile(BGDir,'surface',subj_name{sn(s)},hemName{h},sprintf('%s.%s_%s_sess%d_nanmean.metric',hem{h},imageType,OUTname{j},ss));
                        end
                        outcolnames{indx}{s} = sprintf('%s_%s_%s_sess%d',subj_name{sn(s)},imageType,OUTname{j},ss);
                        % Name the output filename for this group metric file in average surface folder
                    end;
                    if length(OUTname)>3
                        outfilenames{indx} = [surfaceGroupDir filesep hem{h} '.' imageType '_RS_' OUTname{j} '_sess' num2str(ss) '.metric'];
                    else
                        outfilenames{indx} = [surfaceGroupDir filesep hem{h} '.' imageType '_' OUTname{j} '_sess' num2str(ss) '.metric'];
                    end
                    % Finally, make the group metric file for this metric type/contrast
                    caret_metricpermute(infilenames{indx},'outfilenames',outfilenames(indx),'inputcol',inputcol(indx),'replaceNaNs',replaceNaN(indx),'outcolnames',outcolnames{indx});
                    % Verbose display to user
                    fprintf('hem: %i  sess: %i \n',h,ss);
                end;
            end;
        end;   
    case 'SURF_group_cSPM'
        sessN=[1:4];
        sn=[4:9,11:25];
        SPMname={'trained','untrained','cross'};
        imageType='dist';

        vararginoptions(varargin,{'sessN','sn','sqrtTransform','SPMname','imageType'});
                
        if ~strcmp(imageType,'psc') % Should you take ssqrt before submitting? yes for distances
            sqrtTransf = zeros(1,length(SPMname));
        else
            sqrtTransf = ones(1,length(SPMname));
        end
        
        s=1:length(sn);   
        for ss=sessN
            if length(SPMname)>3
                SummaryName = sprintf('.summary_%s_RS_sess%d.metric',imageType,ss);
                imageName = sprintf('%s_RS',imageType);
            else
                SummaryName = sprintf('.summary_%s_sess%d.metric',imageType,ss);
                imageName = imageType;
            end
            
            for h=1:2
                surfaceGroupDir = fullfile(BGDir,'surface','group',hemName{h});
                cd(surfaceGroupDir);
                %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
                for i=1:length(SPMname);
                    if length(SPMname)>3
                        filenames{i}=[surfaceGroupDir filesep hem{h} '.' imageName '_' SPMname{i} '_sess' num2str(ss) '.metric'];
                    else
                        filenames{i}=[surfaceGroupDir filesep hem{h} '.' imageName '_' SPMname{i} '_sess' num2str(ss) '.metric'];
                    end
                end;
                %----loop over the metric files and calculate the cSPM of each with the non-smoothed metrics
                for i=1:length(SPMname);
                    Data=caret_load(filenames{i});
                    if sqrtTransf(i)
                        Data.data=ssqrt(Data.data);
                    end;
                    cSPM=caret_getcSPM('onesample_t','data',Data.data(:,s),'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                    caret_savecSPM([surfaceGroupDir filesep hem{h} '.' imageName '_' SPMname{i} '_stats.metric'],cSPM);
                    save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');
                    data(:,i)=cSPM.con(1).con; % mean
                    data(:,i+length(SPMname))=cSPM.con(1).Z; % T
                    column_name{i}=['mean_' SPMname{i} '_sess' num2str(ss)];
                    column_name{i+length(SPMname)}=['T_' SPMname{i} '_sess' num2str(ss)];
                end;
                C = caret_struct('metric','data',data,'column_name',column_name);
                caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);
                clear Data C cSPM data;
            end;
            fprintf('Done sess-%d %s\n',ss,imageType);
        end
    case 'SURF_zscore'
        sessN = 1:4;
        imageType ='psc';
        for h=1:2
            surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);
            for ss=sessN
                clear RGBdata;
                C = caret_load([hem{h} sprintf('.summary_%s_sess%d.metric',imageType,ss)]);
                for i=1:C.num_cols
                    data(:,i) = zscore(C.data(:,i));
                    column_name{i} = sprintf('zscore_%s',C.column_name{i});
                end
                C_new=caret_struct('metric','data',data,'column_name',column_name);
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_zscore_sess%d.metric',imageType,ss)],C_new);
                clear C C_new data;
            end
        end;
    case 'SURF_group_RGB'
        % plot trained as red, untrained as blue
        sessN = 1:4;
        imageType ='psc';
        minLim = 0.05;
        maxLim = 0.8;
        vararginoptions(varargin,{'imageType','minLim','maxLim'});
        for h=1:2
            surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);
            for ss=sessN
                clear RGBdata;
                C = caret_load([hem{h} sprintf('.summary_%s_sess%d.metric',imageType,ss)]); 
                RGBdata(:,1) = C.data(:,2); % Red: trained
                RGBdata(:,3) = C.data(:,3); % Blue: untrained         
                sc=[minLim maxLim;minLim maxLim;minLim maxLim];  % scaling
                name={sprintf('%s_sess-%d_trained_untrained',imageType,ss)};
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_sess%d.RGB_paint',imageType,ss)],C);
            end
        end;
    case 'SURF_group_sess_RGB'
        % black asunderlay (sess 1), sess 4 : trained as red, untrained as blue
        imageType ='psc';
        minLim = 0.05;
        maxLim = 0.8;
        vararginoptions(varargin,{'imageType','minLim','maxLim'});
        sc=[minLim maxLim;minLim maxLim;minLim maxLim];  % scaling
        for h=1:2
            surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);
            C = caret_load([hem{h} sprintf('.summary_%s_sess1.metric',imageType)]);
            RGBdata(:,1:3) = repmat(C.data(:,2),1,3); % trained
            name={sprintf('%s_sess1_trained',imageType)};
            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
            caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_sess1_maskTrained.RGB_paint',imageType)],C);
            clear RGBdata;
            RGBdata(:,1:3) = repmat(C.data(:,3),1,3); % Blue: untrained
            name={sprintf('%s_sess1_untrained',imageType)};
            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
            caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_sess1_maskUntrained.RGB_paint',imageType)],C);
            clear RGBdata;
            C = caret_load([hem{h} sprintf('.summary_%s_sess4.metric',imageType)]);
            RGBdata(:,1) = C.data(:,2); % trained
            RGBdata(:,3) = C.data(:,3); % untrained
            name={sprintf('%s_sess4_trained_vs_untrained',imageType)};
            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
            caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_sess4_trained_untrained.RGB_paint',imageType)],C);
            clear RGBdata;
        end;
    case 'SURF_group_sess1-4'
        % calculate difference in activation from session 1 to 4
        imageType ='psc';
        seqType = {'AllSeq','TrainSeq','UntrainSeq'};
        vararginoptions(varargin,{'imageType','seqType'});
        for h=1:2
            surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);
            C1 = caret_load([hem{h} sprintf('.summary_%s_sess1.metric',imageType)]);
            C2 = caret_load([hem{h} sprintf('.summary_%s_sess4.metric',imageType)]);
            for i=1:length(seqType)
                data(:,i)=C1.data(:,i)-C2.data(:,i);
                column_name{i} = sprintf('diff_%s',seqType{i});
            end
            C = caret_struct('metric','data',data,'column_name',column_name);
            caret_save([surfaceGroupDir  filesep hem{h} '.psc_difference_sess1-4.metric'],C);
            clear C1 C2 data;
            
        end;
    case 'SURF_group_repsup_RGB'
        % plotting RS on surface as absolute
        % subtracting 2nd - 1st activation - ABSOLUTE RS
        sessN = [1:4];
        imageType = 'psc_RS'; % psc_RS, dist_RS
        type = 'suppression'; % suppression or enhancement 
        vararginoptions(varargin,{'sessN','sessN','imageType','maxLim','type'});
        for h=1:2
            surfaceGroupDir=fullfile(BGDir,'surface','group',hemName{h});
            cd(surfaceGroupDir);      
            for ss=sessN
                clear RGBdata;
                C = caret_load([hem{h} sprintf('.summary_%s_sess%d.metric',imageType,ss)]);
                % calculate difference 1st - 2nd exe
                minLim=0;
                % maxLim
                RGBdata(:,1) = C.data(:,1)-C.data(:,2); % Red: trained 1st-2nd BOLD signal
                RGBdata(:,3) = C.data(:,3)-C.data(:,4); % Blue: untrained 1st-2nd BOLD signal
                switch imageType
                    case 'psc_RS'
                        maxLim = 0.05;
                    case 'dist_RS'
                        maxLim = 0.01;
                        %RGBdata = ssqrt(RGBdata);
                end
                switch type
                    case 'suppression'
                        RGBdata(RGBdata(:,1)<minLim,1) = 0;
                        RGBdata(RGBdata(:,3)<minLim,3) = 0;  
                    case 'enhancement'
                        RGBdata(RGBdata(:,1)>minLim,1) = 0;
                        RGBdata(RGBdata(:,3)>minLim,3) = 0; 
                        RGBdata = RGBdata * (-1)';
                end
                sc=[minLim maxLim;minLim maxLim;minLim maxLim];  % scaling
                name={sprintf('%s_absolute_%s_sess%d',imageType,type,ss)};
                
                C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc},'column_name',name);
                caret_save([surfaceGroupDir filesep hem{h} sprintf('.%s_%s_sess%d.RGB_paint',imageType,type,ss)],C);
            end       
        end;
        
    case 'quick_hack'
        sessN   = [1:4];
        sn      = [4:9,11:25];
        parcelType = 'group_MNI_striatum_buckner';
        
        for ss=sessN
            DD = [];
            T = load(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss))); % loads region data (T)
            % harvest
            for s=sn % for each subj
                fprintf('\nSubject: %d\n',s) % output to user
                roi=unique(T.region)';  % for each region
                for r=roi
                    D = getrow(T,(T.SN==s & T.region==r)); % subject's region data
                    D = getrow(D,1);
                    if sum(sum(isnan(D.betaUW{:})))~=1%~isnan(D.betaUW{:})
                        D.betaUW{1} = real(D.betaUW{:});
                    end
                    DD = addstruct(DD,D);
                end
            end
            % save T
            save(fullfile(regDir,sprintf('betas_%s_FoSEx_sess%d.mat',parcelType,ss)),'-struct','DD');
            fprintf('\n');
        end
        
    case 'run_job'
        %sml1_imana_BG_new('FUNC_subjAvrg','seqType',{'AllSeq','TrainSeq','UntrainSeq'},'imageType','spmT','sn',[5:9,11:31]);
        %sml1_imana_BG_new('FUNC_maskImages_individSubj','seqType',{'AllSeq','TrainSeq','UntrainSeq'},'imageType','psc','sn',[5:9,11:31]);
        %sml1_imana_BG_new('SURF_makeMetric','seqType',{'AllSeq','TrainSeq','UntrainSeq'},'sn',[5:9,11:31]);
        sml1_imana_BG_new('SURF_makeGroup','OUTname',{'AllSeq','TrainSeq','UntrainSeq'},'imageType','psc','sn',[5:9,11:31]);
        sml1_imana_BG_new('SURF_group_cSPM','SPMname',{'AllSeq','TrainSeq','UntrainSeq'},'imageType','psc','sn',[5:9,11:31]);
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

function C          = calcDist(D,betaW,G)
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
function M          = pcm_repsupModel_generic

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
function M          = pcm_repsupModel_specific
    
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
function T          = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
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
function C          = pcm_correlation(Data,partVec,condVec,M,runEffect,M_type)
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
function r          = calcCorr(G)
d0=diag(G);
v1 = d0(1:6)';    % Variances first exe
v2 = d0(7:12)';   % Variances 2nd exe
cv=diag(G,6);     % Covariance
r = mean(cv)/sqrt(mean(v1)*mean(v2));
end
function M          = combineMask(maskFiles)
numFiles=size(maskFiles,2);
for i=1:numFiles
    Reg_mask=spm_vol(maskFiles{i});
    Reg_mask2=spm_read_vols(Reg_mask);
    if i==1
        M=Reg_mask2;
    else
        M=M+Reg_mask2;
    end
end

end
function [reg,hemi] = namingPlots(parcelType,regname_striatum,regname_BG,regname_striatum_networks)
% provide naming for plots
switch(parcelType)
    case 'BG-striatum'
        reg=repmat(regname_striatum,1,2);
        hemi={'lh','lh','rh','rh'};
    case 'BG-striatum-networks'
        reg=regname_striatum_networks;
        hemi=repmat({'both'},1,numel(reg));
    case 'striatum_buckner'
        reg=regname_striatum_networks;
        hemi=repmat({'both'},1,numel(reg));
    case 'group_MNI_striatum_buckner'
        reg=repmat(regname_striatum_networks,1,2);
        hemi=[repmat({'lh'},1,7),repmat({'rh'},1,7)];
    case 'striatum_buckner_hemi'
        reg=[regname_striatum_networks regname_striatum_networks];
        hemi=[repmat({'lh'},1,7),repmat({'rh'},1,7)];
end
end
function makeGroupMask(meanVoxNum,groupRaw,groupMask)
% makeGroupMask(meanVoxNum,groupRaw,groupMask)
% makes group mask file from probabilistic file given average size across subjects
% meanVoxNum - mean number of voxels in a region across subjects
% groupRaw - raw group nifti file (probabilistic)
% groupMask - resultant mask file (0 and 1 value)

    VolGroup=spm_vol(groupRaw);
    VolGroup2=spm_read_vols(VolGroup);
    voxIndx=find(VolGroup2);
    % sort voxels by numbers (how many subjects in common)
    [sortVox indx]=sort(VolGroup2(voxIndx),'descend');
    % choose the same number of voxels with most overlap
    newVox=voxIndx(indx(1:meanVoxNum));
    % calculate mask - preserving the volume
    maskVol=zeros(size(VolGroup2));
    maskVol(newVox)=1;
    maskNii=VolGroup;
    maskNii.fname=groupMask;
    maskNii=rmfield(maskNii,'descrip');
    maskNii.descrip='avrg_mask';
    spm_write_vol(maskNii,maskVol);
 
end

function [R,G] = splitHalfCorr(data,partVec,condVec,type)
% function R = splitHalfCorr(data,partVec,condVec,type)
% performs split-half correlations
switch(type)
    case 'within'
        % calculate within session split-half correlation
        X = indicatorMatrix('identity_p',condVec);
        G = crossval_estG(data,X,partVec);
        mPat=[];
        split=mod(partVec,2);
        
        mPat(:,:,1)     = pinv(X(split==0,:))*data(split==0,:);
        mPat(:,:,2)     = pinv(X(split==1,:))*data(split==1,:);
        % indx - find only voxels which don't have a nan 
        indx            = find(~isnan(mPat(1,:,1))&~isnan(mPat(1,:,2)));
        COR             = corr(mPat(:,indx,1)',mPat(:,indx,2)');
        R.corr          = fisherinv(mean(fisherz(diag(COR))));
        meanPat         = mean(mPat,1);             % Mean activity pattern over fingers
        mPat            = bsxfun(@minus,mPat,meanPat);  % Subtract out mean
        COR             = corr(mPat(:,:,1)',mPat(:,:,2)');
        R.corr_noMean   = fisherinv(mean(fisherz(diag(COR))));
        R.corr_mean     = corr(meanPat(:,:,1)',meanPat(:,:,2)');
        
    case 'across'
        Y     = [data{1};data{2}];
        part  = [partVec{1};partVec{2}];
        cond  = [condVec{1};condVec{2}];
        ses   = [ones(size(partVec{1}));ones(size(partVec{2}))*2];
        X     = indicatorMatrix('identity_p',cond);
        G = crossval_estG(Y,X,ses);
        mPat=[];
        
        % Do split half correlations
        split = mod(part,2);
        mPat(:,:,1) = pinv(X(ses==1 & split==0,:)) * data{1}(split(ses==1)==0,:);
        mPat(:,:,2) = pinv(X(ses==1 & split==1,:)) * data{1}(split(ses==1)==1,:);
        mPat(:,:,3) = pinv(X(ses==2 & split==0,:)) * data{2}(split(ses==2)==0,:);
        mPat(:,:,4) = pinv(X(ses==2 & split==1,:)) * data{2}(split(ses==2)==1,:);
        % indx - find only voxels which don't have a nan (across all 4)
        indx        = find(~isnan(mPat(1,:,1))&~isnan(mPat(1,:,2))&~isnan(mPat(1,:,3))&~isnan(mPat(1,:,4)));
        COR=[];
        COR(:,:,1)  = corr(mPat(:,indx,1)',mPat(:,indx,3)');
        COR(:,:,2)  = corr(mPat(:,indx,1)',mPat(:,indx,4)');
        COR(:,:,3)  = corr(mPat(:,indx,2)',mPat(:,indx,3)');
        COR(:,:,4)  = corr(mPat(:,indx,2)',mPat(:,indx,4)');
        R.corr      = fisherinv(mean(diag(mean(fisherz(COR),3))));
        meanPat     = mean(mPat,1);             % Mean activity pattern over fingers
        mPat        = bsxfun(@minus,mPat,meanPat);  % Subtract out mean
        COR=[];
        COR(:,:,1)  = corr(mPat(:,indx,1)',mPat(:,indx,3)');
        COR(:,:,2)  = corr(mPat(:,indx,1)',mPat(:,indx,4)');
        COR(:,:,3)  = corr(mPat(:,indx,2)',mPat(:,indx,3)');
        COR(:,:,4)  = corr(mPat(:,indx,2)',mPat(:,indx,4)');
        R.corr_noMean = fisherinv(mean(diag(mean(fisherz(COR),3))));
        COR=[];
        COR(1) = corr(meanPat(:,indx,1)',meanPat(:,indx,3)');
        COR(2) = corr(meanPat(:,indx,1)',meanPat(:,indx,4)');
        COR(3) = corr(meanPat(:,indx,2)',meanPat(:,indx,3)');
        COR(4) = corr(meanPat(:,indx,2)',meanPat(:,indx,4)');
        R.corr_mean = fisherinv(mean(fisherz(COR)));
end
end

function output = calcCorrSearch(Y,SPM,conditionVec,partitionVec)
% Called from rsa_runSearchlight in 'SEARCH_calcCorr'
% calculate correlation distance
[betaW,res_MS,~,beta_hat] = rsa.spm.noiseNormalizeBeta(Y,SPM);  % get raw beta regressor weights
betaUW  = bsxfun(@rdivide,beta_hat,sqrt(res_MS));           % apply univariate whitening to beta regressors (divide by voxel's variation)
[G,Sig] = pcm_estGCrossval(betaW,partitionVec,conditionVec); % use multivariately prewhitened betas for searchlight (discuss with Joern!)
C       = corr_crossval(G,'reg','minvalue');
output  = C';
end
