function [acoor,asol,atri,aedg, acrn] = adapMesh(coor,sol,tri,edg,crn,met,hgrad) 
% ADAPMESH : perform 2d mesh adaptation
%      outputs new adapted mesh acoor/atri arrays 
%      that are adapted from initial coor/tri mesh with the metric met 
%      under a gradation hgrad.
%
%      Make sure to set the path to feflo.a_2d to the right directory
%

if ( hgrad <= 1 )
  error('Invalid gradation parameter ');
end  

% add ref to triangle 
if ( size(tri,1) == 3 )
  tri = [tri ; ones(1,szie(tri,2))];
end
  
if ( size(tri,1) ~= 4 )
  error('Wrong shaped triangles (tri) array must be [1:3,1:size(tri,2)]');
end

if ( size(coor,1) ~= 2 )
  error('Wrong shaped vertex (coor) array must be [1:2,1:size(coor,2)]');
end

if ( size(met,1) ~= 3 )
  error('Wrong shaped metric (met) array must be [1:3,1:size(coor,2)] or [1,1:size(coor,2)] ');
end 

if ( size(met,2) ~= size(coor,2) )
  error(' Inconsistant number of vertex and solutions size(met,2) != size(coor,2)');
end

delete temp.y.mesh
delete temp.y.sol

writemesh('temp.y.mesh',coor,tri,edg,crn);
writesol('temp.y.sol',2 , 3, met);
if ( length(sol) ~= 0 )
  delete temp.u.sol 
  writesol('temp.u.sol',2,1,sol);  
  % un comment proper executable for your host 
 system('../bin/mac/feflo.a_2d  -in temp.y.mesh -met temp.y.sol -out temp.out.mesh -hgrad hgrad -noref -itp temp.u.sol');
 %! ./bin/linux/feflo.a_2d  -in temp.y.mesh -out temp.out.mesh -hgrad hgrad -noref -itp temp.u.sol
 %! ./bin/win/feflo.a_2d  -in temp.y.mesh -out temp.out.mesh -hgrad hgrad -noref -itp temp.u.sol

 [dim,acoor,atri,atet,aedg,acrn] = readmesh('temp.out.mesh'); 
 asol                   = readsol('temp.u.itp.sol');
else
 system('../bin/mac/feflo.a_2d  -in temp.y.mesh -met temp.y.sol -out temp.out.mesh -hgrad 1.5 -noref'); 
 [dim,acoor,atri,atet,aedg,acrn] = readmesh('temp.out.mesh'); 
 asol = [];
end
  
