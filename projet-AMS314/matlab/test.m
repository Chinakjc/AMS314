%% Checks of provided functions


% list of functions to be used or completed 
help adapMesh.m
help buildBall.m
help getBall.m 
help im2mesh.m	
help popStack.m	
help readsol.m	
help writemesh.m
help addhashTable.m
help buildCavity.m	
help getEdg.m	
help inCavity.m	
help pushStack.m 
help setBall.m	
help writesol.m
help bdyCavity.m	
help delaunay.m	
help getHess.m	
help iniHashTable.m	
help readmesh.m


% example of reading  2D/3D meshes 
[dim,coor,tri,tet,edg,crn] = readmesh('../data/carre_4h.mesh');
disp([' dimension           : ' num2str(dim) ]);
disp([' number of vertices  : ' num2str(size(coor,2)) ]);
disp([' number of triangles : ' num2str(size(tri,2)) ]);

% getting the edge of the mesh
edg = getEdg(tri(1:3,:));
disp([' number of edge in the mesh : ' num2str(size(edg,2)) ])

[link, ball] = setBall(edg);
ver = getBall(1,link,ball);
disp([' number of neigh. vertices around vertex 1 : ' num2str(ver)]);


% example of computing the hessian of x^2 + y^2
sol =  coor(1,:).^2 + coor(2,:).^2;
[gra,hes] = getHess(coor,tri,sol,link,ball);



[dim,coor,tri,tet,edg,crn] = readmesh('../data/carre_2h.mesh');
disp([' dimension           : ' num2str(dim) ]);
disp([' number of vertices  : ' num2str(size(coor,2)) ]);
disp([' number of triangles : ' num2str(size(tri,2)) ]);


[dim,coor,tri,tet,edg,crn] = readmesh('../data/ssbj.mesh');
disp([' dimension            : ' num2str(dim) ]);
disp([' number of vertices   : ' num2str(size(coor,2)) ]);
disp([' number of triangles  : ' num2str(size(tri,2)) ]);
disp([' number of tetrahedra : ' num2str(size(tet,2)) ]);


[dim,coor,tri,tet,edg,crn] = readmesh('../data/avion_nogeom.mesh');
disp([' dimension            : ' num2str(dim) ]);
disp([' number of vertices   : ' num2str(size(coor,2)) ]);
disp([' number of triangles  : ' num2str(size(tri,2)) ]);
disp([' number of tetrahedra : ' num2str(size(tet,2)) ]);


[dim,coor,tri,tet,edg,crn] = readmesh('../data/cube_interface.mesh');
disp([' dimension           : ' num2str(dim) ]);
disp([' number of vertices  : ' num2str(size(coor,2)) ]);
disp([' number of triangles : ' num2str(size(tri,2))  ]);
disp([' number of tetrahedra : ' num2str(size(tet,2)) ]);

% for  triangles/tetrahedra, tri(4,itri) is a "reference" and tet(5,itet) is volume reference.


% reading an image
[coor,tri,edg,crn,sol] = im2mesh('../data/joconde.ppm');
disp([' number of vertices  : ' num2str(size(coor,2)) ]);
disp([' number of triangles : ' num2str(size(tri,2)) ]);
disp([' size of solution    : ' num2str(size(sol,1)) ' x ' num2str(size(sol,2)) ]);



% adapMesh :  to use adapMesh, you need to edit the file 
% to set the link to the propre fefloa_2d executable (max/linux/windows)

% Coarsening of the carre_4h
[dim,coor,tri,tet,edg,crn] = readmesh('../data/carre_4h.mesh');
h = 0.2;
met =  [1/h^2 ; 0 ; 1/h^2] * ones(1,size(coor,2));
[acoor,asol,atri,aedg, acrn] = adapMesh(coor,[],tri,edg,crn,met,1.5);  


% Example of adaptation with an interpolated solution 
[dim,coor,tri,tet,edg,crn] = readmesh('../data/carre_4h.mesh');
h = 0.1;
met =  [1/h^2 ; 0 ; 1/h^2] * ones(1,size(coor,2));
sol = coor(1,:).^2 + coor(2,:).^2;
[acoor,asol,atri,aedg, acrn] = adapMesh(coor,sol,tri,edg,crn,met,1.5);  







