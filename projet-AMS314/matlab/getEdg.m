function edg = getEdg(tri)
% GETEDG Compute edges of mesh
% from triangles array tri

if ( size(tri,1) >= 3 )
  edg = [tri([1,2],:),tri([1,3],:),tri([2,3],:)]'; 
  edg = unique(sort(edg,2),'rows')';
else
  error('Argument tri has the wrong shape !');
end