function [coord,tri] = buildBall(coord,tri,ip,cavite)
% BUILDBALL : build the ball of all vertices in the mesh
% remove the elements in the cavity and 
% define the ball of Xp

bdy = bdyCavity(coor,tri,cavity)

% new triangles are created and added at the end of table tri newTri = zeros(3,size(bdy,2));
for i=1:size(bdy,2)
  newTri(1,i) = bdy(1,i);
  newTri(2,i) = bdy(2,i);
  newTri(3,i) = ip;
  %% chech the positivness of the new triangles, if not verified error
end

% add the new triangles 
tri = [ tri , newTri];
