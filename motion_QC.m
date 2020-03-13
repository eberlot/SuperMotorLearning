function M=motion_QC(SPM,varargin)
% function M=motion_QC(SPM,varargin)
% calculates quality control parameters based on subject motion
% - absolute motion (root mean square) 
% - framewise displacements 
%
% INPUT: 
%   SPM:  SPM structure estimated for GLM with the RWLS toolbox
% VARARGIN: 
%   'subjDir': directory where motion parameters are saved
%              otherwise directory defined from SPM
%   'run_indx': run indices for estimation
%               otherwise starts from 1 until number of runs, as defined with SPM
%   'prefix':  specify if movement parameter given any prefix (e.g. 'rp_')
%              otherwise assumes none
% OUTPUT: 
%   M: Structure with motion parameters, fields:
%      - fwd: framewise displacement (in mm)
%      - rms: root mean square displacement 


% timing information

subjDir = [];
vararginoptions(varargin,{'subjDir','prefix'});

[dirAll,~]=spm_fileparts(SPM.xY.VY(1).fname);
% determine subject directory from SPM if not given
if isempty(subjDir)
    subjDir=dirAll;
end

start=[0 cumsum(SPM.nscan)]+1;
% call text files with movement parameters
for i=1:length(SPM.nscan)
    [~,filename]=spm_fileparts(SPM.xY.VY(start(i)).fname);
    filename = filename(2:end); % temporary
    if exist('prefix')
        filename = [prefix filename];
    end
    movparam_name{i}=[subjDir filesep filename '.txt'];
end;

%-----------------------------------------------------------------------
% Load movement parameter files
MOV=[];
DTS=[];

for i=1:length(movparam_name)
    if (~isempty(movparam_name{i}))
        try
            mov=dlmread(movparam_name{i});
            MOV=[MOV;mov];
            dts=diff(mov);
            dts=[zeros(1,size(dts,2)); dts];  % first element is a zero
            DTS=[DTS;dts];
        catch
            warning(['Movementparameter file ' movparam_name{i} ' not found.']);
        end;
        if (size(mov,1)~=SPM.nscan(i))
            warning('Number of scans in movement parameter file do not match information in SPM.mat');
            mov=[];
        end;
          
    end;
end;

%-----------------------------------------------------------------------
% Calculate indices
% Translation - convert degrees into motion in mm;
radius=50;
temp=MOV(:,4:6);
temp=radius*temp;
MOV(:,4:6)=temp;

DTS2 = [zeros(1,size(DTS,2));diff(MOV)];
% replace the starts of the run with zeros
nImage = size(DTS2,1)/length(SPM.nscan);
imRepl = repmat(nImage,length(SPM.nscan),1).*(0:1:(length(SPM.nscan)-1))'+1; % images to be replaced
%DTS2(imRepl,:) = 0;
M.fwd=sum(abs(DTS2),2)';
rms=sqrt(mean(MOV.^2,1));    % root mean square for each column
M.fwd_trans = sum(abs(DTS(:,1:3)),2)'; % translation
M.fwd_rot = sum(abs(DTS2(:,4:6)),2)'; % translation
M.fwd_mean = mean(M.fwd);
    
end