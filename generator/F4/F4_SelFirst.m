% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015
%
% Selection strategy for F4 Algorithm. Selects the first pair from the set.
%
% POld - set of all pairs
% PSel - set of selected pairs
% PNew - set of remaining pairs


function [PSel, PNew] = F4_SelFirst(POld)
  
  PSel = POld(1:1,1);
  PNew = POld(2:end,1);
  
end

