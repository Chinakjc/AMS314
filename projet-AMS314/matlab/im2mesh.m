function [coor,tri,edg,crn,sol] = im2mesh(file)
% IM2MESH : convert an image to couple mesh/solution

% load image
im     = imread(file);
info   = imfinfo(file);
Width  = info.Width;
Height = info.Height;
 
disp([ '% ' file ' OPENED ']);
 
% convert image to a grey-level image from rgb
% grey = 0.3*r+0.59*g+0.11*b

im = 0.3*im(:,:,1) + 0.59*im(:,:,2) + 0.11*im(:,:,3);
%imshow(im);

nrow = Height;
ncol = Width;

tri  = [];
edg  = [];
crn  = [];

coor  = zeros(2,nrow*ncol);
sol   = zeros(1,nrow*ncol);
tri1  = zeros(3,(nrow-1)*(ncol-1));
tri2  = zeros(3,(nrow-1)*(ncol-1));
edg   = zeros(3,2*(ncol-1)+2*(nrow-1));
crn   = [ 1, ncol, 1 + ncol*(nrow-1), ncol*(nrow-1) + ncol ];

idx = 1;
for i=1:nrow

  ind = idx:(idx+ncol-1);
  coor(1,ind ) = 1:ncol;
  coor(2,ind ) = nrow-i+1;
  sol(ind)     = im(i,1:ncol);
  idx = idx + ncol ;
   
% for j=1:ncol
%   coor(:,idx) = [ j , nrow-i+1 ] ;
%   sol(idx )   = im(i,j);
%   idx = idx + 1;  
% end
 
end
  
  
%tic;  
idx = 1;  
indj = 0:(ncol-2);

for i=0:(nrow-2)
        
   ind  = idx:(idx+length(indj)-1);
   in   = indj + i*ncol+1;
   
   tri1(1,ind) = in;
   tri1(3,ind) = in+1;
   tri1(2,ind) = in+ncol+1;
  
   tri2(1,ind) = in;
   tri2(3,ind) = in+ncol+1;
   tri2(2,ind) = in+ncol;
   idx = idx+length(indj) ;
     
% very slow
% for j=0:(ncol-2)
%   in = j + i*ncol + 1 ;
%   tri(:,idx) = [in;in+1;in+ncol+1];
%   idx = idx + 1; 
%   tri(:,idx) = [in;in+ncol+1;in+ncol];
%   idx = idx + 1;
% end

end
%toc; 

tri  = [tri1,tri2];
i1 = crn(1);
i2 = crn(2);
i3 = crn(3);
i4 = crn(4); 
edg(1,:)  = [i1:(i2-1), i3:(i4-1), i1:ncol:(i3-ncol), i2:ncol:(i4-ncol) ] ; 
edg(2,:)  = [(i1+1):i2, (i3+1):i4, (i1+ncol):ncol:i3, (i2+ncol):ncol:i4 ] ;