function [out] = corr_crossval(G,varargin)
% calculates crossvalidated correlations between activity patterns across
% conditions.
% Correlation distance is then the compliment of resulting corrs.
%   INPUT:
%           G:  [numConds*numConds] crossvalivadted second moment matrix 
%   VARARGIN:
%           reg: regularization options
%           abs - make the variances into absolute values
%           minvalue - if variance lower than 0.001, make it 0.001
%   OUTPUT:
%          out:  [n x numConds] vector of crossvalidated cosine distances.
%                       Order of condition pairs is in rsa_vectorizeRDM order
%                       Use rsa_squareRDM to invert it back


reg = [];

pcm_vararginoptions(varargin,{'reg'});

numCond = size(G,1);
% prep output matrix
out = zeros(numCond,numCond);
% calculate pairwise cosine distances
for n = 1:numCond
    for nN = n+1:numCond
        
        switch(reg)
            case 'abs'
                out(nN,n) = (1 - (G(n,nN)/sqrt(abs(G(n,n)*G(nN,nN)))))/2;
            case 'minvalue'
                if G(n,n)<0.0001
                    G(n,n)=0.0001;
                end
                if G(nN,nN)<0.0001
                    G(nN,nN)=0.0001;
                end
                out(nN,n) = (1 - (G(n,nN)/sqrt(G(n,n)*G(nN,nN))))/2;
            otherwise
                out(nN,n) = (1 - (G(n,nN)/sqrt(G(n,n)*G(nN,nN))))/2;
        end
    end
end
% vectorize output matrix
out = rsa_vectorizeRDM(out);

end