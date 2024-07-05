function [dim,coor,tri,tet,edg,crn] = readmesh(name)
% READMESH : open the mesh file name
%            outputs the vertices and 
%            triangles arrays
%  example: [coor,tri,edg,crn] = readmesh('toto.mesh')
%   


coor = [];
tri  = [];
tet  = [];
crn  = [];
edg  = [];


fid = fopen(name,'r');
if ( fid == -1 ) 
 error(['Mesh file ' name ' does not exist']);
else
 disp([ '% ' name ' OPENED ']);
end


while ( ~feof(fid) )

  str = fgetl( fid );
  str = str(find(str~=' '));
  switch ( lower(str) )
  
  case 'dimension'
    dim = fscanf( fid, '%d', 1 );
    if ( dim ~= 2 && dim~=3 )
      error(' Invalid inout mesh ');
    end  
        
  case 'dimension3'
    dim = 3;
  
  case 'dimension2'
    dim = 2;
  case 'vertices'
    NbrVer = fscanf( fid, '%d', 1 );
    if ( dim == 2 )
      coor = fscanf( fid, '%f %f %d', [3, NbrVer] );
      coor = coor(1:2,:);
    else
      coor = fscanf( fid, '%f %f %f %d', [4, NbrVer] );
      coor = coor(1:3,:);
    end  
   
  case 'triangles'
    NbrTri = fscanf( fid, '%d', 1 );
    tri = fscanf( fid, '%d %d %d %d', [4, NbrTri] );
    
  case 'tetrahedra'
    NbrTet = fscanf( fid, '%d', 1 );
    tet = fscanf( fid, '%d %d %d %d', [5, NbrTet] );  
  
  case 'edges'
    NbrEdg = fscanf( fid, '%d', 1 );
    edg = fscanf( fid, '%d %d %d', [3, NbrEdg] );
  
  case 'corners'
    NbrCor = fscanf( fid, '%d', 1 );
    crn = fscanf( fid, '%d ', [1, NbrCor] );

  

  end
end

fclose(fid);
end

