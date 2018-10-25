function CNR = cnr_QC(beta,res,numCond,numRun)
% function CNR = cnr_QC(beta,res,numCond,numRun)
% 
% calculates the average CNR across voxels in a region
% as beta for all conditions across runs, divided by residuals

cnr_cond = bsxfun(@rdivide,beta,res);
CNR = mean(mean(cnr_cond));


end