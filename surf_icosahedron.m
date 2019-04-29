function varargout = surf_icosahedron(atlasDir)
% function varargout = surf_icosahedron(atlasDir)
% makes the label file for surface LR
% INPUT:
% OUTPUT:
% VARARGIN:
% note:
% 2 - 42 vertices comb, 80 faces
% 4 - 162 vertices comb, 320 faces
% 6 - 362 vertices comb, 720 faces
% 8 - 642 vertices comb, 1280 faces
freq    = 6;
radius  = 100;
res     = 164; % using 164 or 42k FS_LR
hem = {'L','R'}; % both hemispheres
hemName = {'CortexLeft','CortexRight'};
vararginoptions(varargin{'freq','radius','res'});
cd(atlasDir);
for h=1:2
    G = gifti(sprintf('fs_LR.%dk.%s.sphere.surf.gii',res,hem{h}));
    [x,y,z,~]=make_icosahedron(freq,radius,1,0,1); % make icosahedron
    DT = delaunayTriangulation(x',y',z'); % triangulate again (Matlab framework)
    [T,Xb] = freeBoundary(DT);
    TR = triangulation(T,Xb);
    P = incenter(TR); % compute the centroids (P)
    nCentr = length(P); % number of centroids
    distMat = zeros(nCentr,size(G.vertices,1));
    for i=1:nCentr
        distMat(i,:)=sqrt(sum(bsxfun(@minus,G.vertices,P(i,:)).^2,2)); % squared Euclidean distance
    end
    % find minimal distance and assign the label of that centroid
    [~,centr]=min(distMat);
    % create a gifti
    G2 = surf_makeFuncGifti(centr','anatomicalStruct',hemName{h});
    save(fullfile(atlasDir,sprintf('fs_LR.%dk.%s.ico-%d.label.gii',res,hem{h},nCentr)));
end

end