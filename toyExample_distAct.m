% SML - repsup
% toy example for explaining activation / distances

% define defaults

voxNum  = 15;
seqNum  = 4;
mu1     = 1;
sigma1  = 1; 
mu2     = 0.7;
sigma2  = 0.6;

E1 = normrnd(mu1,sigma1,seqNum,voxNum);
E2 = normrnd(mu2,sigma2,seqNum,voxNum);

dist1 = pdist(E1);
dist2 = pdist(E2);


% equal axis limits
figure
subplot(8,2,[1,3,5,7])
imagesc(E1);
caxis([-1 2]);
subplot(8,2,[9,11,13,15])
imagesc(E2);
caxis([-1 2]);
subplot(8,2,[6])
imagesc(mean(E1,1));
caxis([-1 2]);
subplot(8,2,[12])
imagesc(mean(E2,1));
caxis([-1 2]);

% calculate scaled version of distances
scale=mean(E1,1)./mean(E2,1);
dist_scale = pdist(bsxfun(@rdivide,E1,scale));

% equal axis limits
figure
subplot(1,3,1)
imagesc(rsa_squareRDM(dist1));
caxis([0 10]);
subplot(1,3,2)
imagesc(rsa_squareRDM(dist_scale));
caxis([0 10]);
% scaled version
subplot(1,3,3)
imagesc(rsa_squareRDM(dist2));
caxis([0 10]);