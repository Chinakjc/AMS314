function [gra,hes] = getHess(coor,tri,sol,link,ball)
% GETHESS hessian/gradiant recovery by least square method 
% see also: GETBALL

hes = zeros(3,size(coor,2)); 
gra = zeros(2,size(coor,2));

ListPb = [];
for iVer=1:size(coor,2)
  vois = getBall(iVer,link,ball); 
  nvois = length(vois);
  if ( nvois <= 5 )
    % for critical points (with less than 5 neigh.) 
    % we add neighboors of neighboors !
    ListPb = [ListPb iVer];
    for ivoi=1:nvois
      vois = [vois getBall(vois(ivoi),link,ball) ]; 
    end
    vois = unique(vois);
    nvois = length(vois); 
  end

  one = ones(nvois,1);


% we write the Taylor expansion around point (x,y)
% f(xv,yv) = f(x,y) + a(xv-x) + b(yv-y) + 1/2(c(xv-x)^2
% + 2d(xv-x)(yv-y) + e(yv-y)^2)

% Unknows are [a,b,c,d,e]
% Hessian matrix is then H = | c,d,|
% | d,e |

  x = coor(1,iVer)*one;
  y = coor(2,iVer)*one;
  f = sol(iVer)*one;
  
  xv = coor(1,vois)'; 
  yv = coor(2,vois)'; 
  fv = sol(vois)';
   
   error('Ecrire le code correspondant !');  
 
end

disp([' Number of critical points  ' int2str(length(ListPb)) ]);
  
  
end