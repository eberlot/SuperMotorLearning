function [R2 R]=rsa_patternConsistency_crossval(Y,partition,conditionVec,varargin)
% function R_cross=rsa_patternConsistency(Y,partition,conditionVec);
% Caluclates a measure of pattern consistency R2 = SS_between / SS_total
% INPUT:
%  Y                : Noise-normalized activation patterns, a KR x P matrix
%  partition        : KR x 1 integer value that indicates the partition for crossvalidation (typically run number), zeros will be ignored
%  conditionVec     : KR x 1 vector of conditions, zeros will be ignored
%
% OUTPUT:
%   R2              : crossvalidated R2 (1-SSR/SST)
%   R               : crossvalidated R (correlation)

% EBerlot, May 2018

removeMean=1; % by default removeMean
vararginoptions(varargin,{'removeMean'});

if ~isnan(Y) % if no data provided
    
    part    = unique(partition)';
    part = part(part~=0)';
    numPart = numel(part);
    numCond = length(unique(conditionVec));
    
    [N,numVox]   = size(Y);
    
    A = zeros(numCond,numVox,numPart);
    X = indicatorMatrix('identity_p',conditionVec);
    
    % Estimate condition means within each run
    for i=1:numPart
        Xa = X(partition==part(i),:);
        Ya = Y(partition==part(i),:);
        A(:,:,i) = pinv(Xa)*Ya;
        if (removeMean) % per partition
            A(:,:,i)=bsxfun(@minus,A(:,:,i),sum(A(:,:,i),1)/numCond);
        end;
    end;
    
    % prep betas - leave one out manner
    for ri=1:numPart
        testRun  = i;
        trainRun = part~=i;
        R        = bsxfun(@minus,A(:,:,testRun),mean(A(:,:,trainRun),3));
        SSR      = sum(sum(sum(R.^2)));
        SST      = sum(sum(sum(A.^2)));
        r2(i)    = 1-SSR/SST;   
        tmp      = corrcoef(A(:,:,testRun),mean(A(:,:,trainRun),3));
        r(i)     = tmp(2);
    end
    R2 = mean(r2);
    R  = mean(r);
else
    warning('no data given');
    R2 = NaN;
    R = NaN;
end


