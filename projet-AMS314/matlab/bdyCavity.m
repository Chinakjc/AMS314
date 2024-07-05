function bord = bdyCavity(cood,tri,cavite)
% BDYCAVITY : computes the boundary of the cavity 
% by means of the table bdy
% the edges of the boundray are
% bdy(1,i) = first vertex of the boundary edge number i 
% bdy(2,i) = second vertex
% The algorithm finds the edges seen only one time
% in the set of elements in the cavity

maxHead = 16;
maxVal = 3*length(Cavity)
[head,list,nlist] = iniHashTable(maxHead,maxVal);


%  we add all the degs of the cavity 
for itri=1:length(cavity)
  ip1 = tri(2,cavity(itri));
  ip2 = tri(3,cavity(itri));
  [head,list,nlist] = addHashTable(ip1, ip2, itri, head, list, nlist)
  ip1 = tri(3,cavity(itri));
  ip2 = tri(1,cavity(itri));
  [head,list,nlist] = addHashTable(ip1, ip2, itri, head, list, nlist)
  ip1 = tri(1,cavity(itri));
  ip2 = tri(2,cavity(itri));
  [head,list,nlist] = addHashTable(ip1, ip2, itri, head, list, nlist)
end

% table bdy is made up of the edges seen only one time in the cavity 
bdy = [];
for iar=1:nlist
  if ( list(4,iar) == 0 )
    bdy = [bdy [list(1,iar);list(2,iar)]];
  end 
end
