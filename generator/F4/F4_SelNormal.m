% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015
%
% Selection strategy for F4 Algorithm. Selects all pairs with a minimal
% total degree. Known as the normal strategy for F4.
%
% POld - set of all pairs
% PSel - set of selected pairs
% PNew - set of remaining pairs


function [PSel, PNew] = F4_SelNormal(POld)
  
  % find minimal total degree
  minDeg = sum(POld{1}.lcm);
  select = [];
  for i = 1:length(POld)
    if sum(POld{i}.lcm) < minDeg
      minDeg = sum(POld{i}.lcm);
      select = [i];
    elseif sum(POld{i}.lcm) == minDeg
      select = [select, i];
    end
  end
  
  % select those pairs
  PSel = POld(select, 1);
  PNew = POld(setdiff(1:length(POld), select), 1);
  
end

