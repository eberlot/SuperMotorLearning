function Z = pcm_buildConditionMatrix(varargin);

% function Z = pcm_buildFeatureModel(name,varargin);
% Generates a feature component model given individual feature matrices.
%==========================================================================
% INPUT:
%       varargin: input with individual condition matrices 
%                 takes in variable amount of feature matrices 
%                 dimensions: condition x features (condNum x featNum)
%                 same condNum across matrices, featNum can differ
%--------------------------------------------------------------------------
% OUTPUT:
%   Z:      Condition matrix

T=varargin;
if iscell(T)
    T=T{:};
end
modelNum=size(T,2);     % number of input models
condNum=size(T{1},1);   % number of conditions

% total number of columns across all condition vectors / matrices
featNum=0;
for f=1:modelNum
    featNum=featNum+size(T{f},2);
end

% initialise condition matrix Z
Z=zeros(condNum,featNum);

% fill the matrix with condition vectors / matrices for each model
colCount=0;
for f=1:modelNum
    Z(:,colCount+1:colCount+size(T{f},2))=T{f};
    colCount=max(colCount+size(T{f},2));
end


end