function ok = inCavity(X0,X1,X2,Xp)
% INCAVITY : checks if triangle with vertices X0,X1,X2 
% belongs to the cavity of point Xp

% compute the circumcircle (radius and centre) 
% e.g. solve a linear system
rayon = 1;               % to be written 
centre = 1/3*(X0+X1+X2); % to be written

      error(' To be written');

ok = 0
if ( norm(Xp-centre,2) <= rayon )
  ok = 1; 
end



