function cavite = builCavity(Xp,coord,tri)
% BUILDCAVITY : naive function to define the cavity of point Xp 
% by travesing all the triangles neighbor of iTri

cavity = [];
for itri=1:size(tri,2)
  % coordinates of the 3 vertices of triangle itri
  X0 = coor(:,tri(1,itri));
  X1 = coor(:,tri(2,itri));
  X2 = coor(:,tri(3,itri));
  if ( inCavity(X0,X1,X2,Xp) == 1 )
     cavity = [cavity itri];
  end 
end

