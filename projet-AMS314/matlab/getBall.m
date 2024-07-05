function ver = getBall(iver,link,ball)
% GETBALL: Gives neigbooring vertices of iver
% 
% see setBall.m and getEdg.m

if ( ~isscalar(iver) || ~isvector(link) || ~isvector(ball) ) 
  error('Wrong arguments type, iver is scalar, link and ball are vectors ');
end
NbrVer = length(link) - 1;
if( (iver>= 1) & (iver <= NbrVer))
  ver = ball(link(iver):(link(iver+1)-1)); 
else
   warning('Ivalid vertex index');
   ver = []; 
 end