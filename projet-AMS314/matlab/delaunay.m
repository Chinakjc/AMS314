function [coord,tri] =  delaunay(coord,tri,ip)
% DELAUNAY : Insert point Xp (index ip) via the Delaunay algorithm

% Get Xp coordinates
Xp          = coord(:,ip); 
cavity      = buildCavity(Xp,coord,tri);
[coord,tri] = buildBall(coord,tri,ip,cavite)


% on enleve les anciens triangles (ceux de la cavite) du tableau tri
idx    = 1:size(tri,2);
idx    = setdif(idx,cavite); 
tri    = tri(:,idx) 


