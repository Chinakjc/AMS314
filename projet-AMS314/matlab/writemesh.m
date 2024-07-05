function ok = writemesh(name,coor,tri,edg,crn)
% WRITEMESH : write the mesh file name
%  example: ok = writemesh('test.mesh',coor,tri,edg,crn)
%           coor vertex coordinates 
%           tri  triangles array  [is1,is2,is3,ref] 
%           edg edges array [is1,is2,ref]  
%           crn corners array 
%        

  k = strfind(name,'.mesh'); 
  if ( length(k) ~= 0 )
    in  = [ name(1:(k(1)-1)) '.mesh'];
  else
    in  = [ name '.mesh'];
  end  
  
  dim = 2;
  ok = 0;

  if ( size(tri,2) ~=0 && size(tri,1) ~= 4 )
    error('Wrong shaped triangles tab must be [1:4,1:NbrTri]');
  end
  
  if ( size(coor,1) == 2 )
    dim = 2;
  elseif ( size(coor,1) == 3 )
    dim = 3
  else
     error('Wrong shaped coordinates tab must be [1:dim,:]');    
  end

  fid = fopen(in,'w+');
  if ( fid == -1 ) 
    error(['Cannot create mesh file ' in ]);
  else
    disp([ '% '  in ' CREATED ']);
  end


fprintf(fid,'MeshVersionFormatted 1\n');
fprintf(fid,'Dimension\n%d\n',dim);
fprintf(fid,'Vertices\n%d\n',size(coor,2));
if (dim == 2 ) fprintf(fid,' %f %f %d \n',[coor ; zeros(1,size(coor,2))]);    end
if ( dim == 3) fprintf(fid,' %f %f %f %d \n',[coor ; zeros(1,size(coor,2))]); end
fprintf(fid,'Triangles\n%d\n',size(tri,2));
fprintf(fid,' %d %d %d %d\n',tri);
fprintf(fid,'Edges\n%d\n',size(edg,2));
fprintf(fid,' %d %d %d \n',edg);
fprintf(fid,'Corners\n%d\n',size(crn,2));
fprintf(fid,' %d \n', (1:1:size(crn,2)));

fprintf(fid,'\nEnd\n');  
fclose(fid);