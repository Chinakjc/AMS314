function ok = writesol(name,dim,typ,sol)
% WRITESOL : create the solution file name
%  example: ok = writesol('temp.sol',1,sol)
%   typ = 1 scalar, 2 vector, 3 metric

ok = 0;

k = strfind(name,'.sol'); 
if ( length(k) ~= 0 )
  in  = [ name(1:(k(1)-1)) '.sol'];
else
  in  = [ name '.sol'];
end  

fid = fopen(in,'w+');
if ( fid == -1 ) 
 error(['Cannot create solution file ' in ]);
else
 disp([ '% '  in ' CREATED ']);
end



fprintf(fid,'MeshVersionFormatted 2\n\n');

if ( dim == 2 )

  fprintf(fid,'Dimension\n2\n\n');
  fprintf(fid,'SolAtVertices\n%d\n1 %d\n',[size(sol,2) typ]);
 
    
  if ( size(sol,1) == 1 )
    fprintf(fid, ' %f \n', sol );
  elseif ( size(sol,1) == 2 )
    fprintf(fid, ' %f %f\n', sol );
  elseif ( size(sol,1) == 3 )
    fprintf(fid, ' %f %f %f\n', sol );
  elseif ( size(sol,1) == 4 )
    fprintf(fid, ' %f %f %f %f \n', sol );
  else
    error(' Wrong solution type ');  
  end

else
  
  fprintf(fid,'Dimension\n3\n\n');
  fprintf(fid,'SolAtVertices\n%d\n1 %d\n',[size(sol,2) typ]);
 
  if ( size(sol,1) == 1 )
    fprintf(fid, ' %f \n', sol );
  elseif ( size(sol,1) == 3 )
    fprintf(fid, ' %f %f %f\n', sol );
  elseif ( size(sol,1) == 6 )
    fprintf(fid, ' %f %f %f %f %f %f\n', sol );
  else
    error(' Wrong solution type ');  
  end 
  
  
  
  
end      
  
fprintf(fid,'\nEnd\n');
fclose(fid);

end

