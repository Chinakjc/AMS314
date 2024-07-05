function [link, ball] = setBall(edg)
% SETBALL Build the link tab and neigbooring
%         vertices array ball
%

if ( size(edg,1) ~= 2 )
  error('Argument edg has the wrong shape !');
end

NbrVer = max(edg(:));
NbrEdg = size(edg,2);

NbrVoi = zeros(1,NbrVer+1);
link   = zeros(1,NbrVer);
ball   = zeros(1,2*NbrEdg);

NbrEdg = NbrEdg;
NbrVer = NbrVer;

% first pass : counting
for iEdg=1:NbrEdg
  NbrVoi(edg(1,iEdg)) = NbrVoi(edg(1,iEdg)) + 1;
  NbrVoi(edg(2,iEdg)) = NbrVoi(edg(2,iEdg)) + 1;
end  


% second pass: build the link tab
link(1) = 1;
for iVer=2:NbrVer
  link(iVer)     = link(iVer-1) + NbrVoi(iVer-1);
end
link(NbrVer+1) = link(NbrVer) + NbrVoi(NbrVer);

NbrVoi(:) = 0;

% third pass: filling the ball
for iEdg=1:NbrEdg
  i1 = edg(1,iEdg);
  i2 = edg(2,iEdg);
  
  ball(link(i1) + NbrVoi(i1)) = i2;
  ball(link(i2) + NbrVoi(i2)) = i1;
  
  NbrVoi(i1) = NbrVoi(i1) + 1;
  NbrVoi(i2) = NbrVoi(i2) + 1;

end
