function [head,list,nlist] = addHashTable(ip1, ip2, itri, head, list, nlist)
% ADDHASHTABLE : add an element to an hash table
% add edge ip1,ip2 of triangle itri in the haching table
% either add a new entry in list at index nlist + 1
% or update list(4) (the neigbor, first choice for the definition of list))
% 

code = ip1 + ip2;
iPos = mod(code,size(head,1)) + 1;


% other keys can be tested
      error('Ecrire le code correspondant !');
% check if (ip1,ip2) is already in list head(iPos) % if not, add the edge in head(iPos) either continue
      error('Ecrire le code correspondant !');
      
      
end