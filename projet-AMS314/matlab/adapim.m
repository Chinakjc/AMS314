% Ce script permet de faire une adaptation de maillage
% pour controller des nx de gris pour une solution analytique
% image_sol = 1  quand la solution est un nx de gris d'une image
% image_sol = 0 pour une solution analytique 

% image_sol = 1
% [coor,tri,edg,crn,sol] = im2mesh('../data/joconde.ppm');
% [coor,tri,edg,crn,sol] = im2mesh('../data/joconde.Lowres.ppm');

 image_sol = 0
 [dim,coor,tri,tet,edg,crn] = readmesh('../data/carre_4h.mesh');


nite = 10;
N    = 1000;

ite = 0;
while ( ite < nite  )

  
  % create an analytical solution when not using an image 
  if ( image_sol == 0 )
    x = coor(1,:);
    y = coor(2,:);
    sol = tanh(-100*(y-0.5-0.25.*sin(2*pi*x))) + tanh (100*(y-x));
  end
  
  %% backup current mesh and solution 
  writemesh(['adap.' num2str(ite) '.mesh'],coor,tri,edg,crn);
  writesol(['adap.' num2str(ite) '.sol'], 2, 1, sol);
  
  
  disp(['Iteration ' num2str(ite) ' complexite ' num2str(N)]);
  
  disp([' Hessian computation ']);
  
  % hessien
  edgM = getEdg(tri);
  [link, ball]   = setBall(edgM); 
  [Grec,Hu]   = getHess(coor,tri,sol,link,ball);
  
   
  Met = zeros(3,size(coor,2));
  Det = zeros(1,size(coor,2));
  
  % compute M_L2 on each vertex
  for i=1:size(coor,2)
    error(' Ecrire le code correspondant !');
   end 
   
  
  % Compute normalization factor N
  disp([' Coefficient de normalisation ']); 
  C   = 0;
  Air = 0; 
  for i=1:size(tri,2)
    error('Ecrire le code correspondant !');
  end 

  
  
  % normalize the metric field
  disp([' Normalisation de la metrique ']); 
  Met = N/C * Met;
   
  
   
  % adapt the mesh and project back the solution 
   disp([' Adaptation du maillage ']); 
  [coor,sol,tri,edg,crn] = adapMesh(coor,sol,tri,edg,crn,Met,1.3);
  
  
  
   

  ite = ite+1;
end
