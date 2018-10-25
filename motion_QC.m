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
vararginoptions(varargin,{'subjDir'});

[dirAll,~]=spm_fileparts(SPM.xY.VY(1).fname);
% determine subject directory from SPM if not given
if isempty(subjDir)
    subjDir=dirAll;
end

start=[0 cumsum(SPM.nscan)]+1;
% call text files with movement parameters
for i=1:length(SPM.nscan)
    [~,filename]=spm_fileparts(SPM.xY.VY(start(i)).fname);
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


fwd=sum(abs(DTS),2);
rms=sqrt(mean(MOV.^2,1));    % root mean square for each column
M.fwd = fwd;
M.rms = rms;
M.fwd = mean(fwd);
M.rms = mean(rms);
    
end