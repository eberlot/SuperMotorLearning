function varargout=sml_searchlight(what,varargin)

Dir = '/Users/eberlot/Documents/Data/SuperMotorLearning/sml_searchlight'; % change
subj_name = 's02'; 
numruns=8;


%% define searchlight
mask       = fullfile('mask_s02.nii',subj_name{s});
Vmask      = spm_vol(mask);
Vmask.data = spm_read_vols(Vmask);

LcaretDir = fullfile(Dir,['x',subj_name],'LeftHem');
RcaretDir = fullfile(caretDir,['x',subj_name],'RightHem');
white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
S         = rsa_readSurf(white,pial,topo);

L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[15 120]);
save(fullfile(Dir,sprintf('%s_searchlight_120.mat',subj_name{s})),'-struct','L');


%% run LDC
block = 5e7;
cwd   = pwd;                                                        % copy current directory (to return to later)

runs = 1:numruns;
% make index vectors
conditionVec  = kron(ones(numel(runs),1),[1:12]');      % 12 sequences
partition     = kron(runs',ones(12,1));
% go to subject's glm directory
cd(fullfile(glmSessDir{sessN},subj_name{s}));
% load their searchlight definitions and SPM file
L = load(fullfile(Dir,sprintf('%s_searchlight_120.mat',subj_name)));

SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{s}));

name = sprintf('%s',subj_name{s});
% run the searchlight
rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,'analysisName',name,'idealBlock',block);


cd(cwd);


end