function [head,list,nlist] = iniHashTable(maxHead,maxVal)
% INIHASTABLE : initiliases an hash table 
% 
% maxVal can be stored
% list(1,)   = first vertex of edge
% list(2,)   = second vertex of edge
% list(3,)   = first triangle having this edge
% list(4,)   = second triangle having this edge
% list(5,)   = next edge in hash tab
% head(iPos) = link to the first edge having iPos as key 
% nlist      = last edge stored in list

head    = zeros(maxHead);
list    = zeros(5,maxVal);
nlist   = 0


