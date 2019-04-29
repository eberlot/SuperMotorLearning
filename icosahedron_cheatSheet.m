% test case on how to make an icosahedron
% note:
% 2 - 42 vertices comb, 80 faces
% 4 - 162 vertices comb, 320 faces
% 6 - 362 vertices comb, 720 faces
% 8 - 642 vertices comb, 1280 faces
freq = 6;
radius = 100;
[x,y,z,TRI]=make_icosahedron(freq,radius,1,1,1);
% triangulate again (Matlab framework)
DT = delaunayTriangulation(x',y',z');
% find the faces
[T,Xb] = freeBoundary(DT);
TR = triangulation(T,Xb);
% compute the centroids (P) and face normals for each triangular facet
P = incenter(TR);
F = faceNormal(TR);  
% plot all
figure
trisurf(T,Xb(:,1),Xb(:,2),Xb(:,3), ...
     'FaceColor','cyan','FaceAlpha',0.8);
axis equal
hold on  
quiver3(P(:,1),P(:,2),P(:,3), ...
     F(:,1),F(:,2),F(:,3),0.5,'color','r');
 
 figure
 plot3(P(:,1),P(:,2),P(:,3),'o')
 hold on;
 quiver3(P(:,1),P(:,2),P(:,3), ...
     F(:,1),F(:,2),F(:,3),2,'color','r');

 % then you will need to find a match between surface nodes (hemisphere)
 % and icosahedron nodes (centroids or nodes? not sure) using the computed
 % normal
 
 % probably best to:
 % 1) take the sphere - find the closest node
 % 2) assign the label
 
 
 % OLD
 % 1) assign the node that intersects with the centroid
 % 2) determine other nodes on the hemisphere to belong to the closest node
 % (Euclidean distance)