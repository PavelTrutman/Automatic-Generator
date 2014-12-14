% calc number of monomials from degree 'fromdeg', to degree 'todeg', 
% containing 'unkcnt' unknowns
%
% by Martin Bujnak, mar2008

function [k] = GetMonomialsCnt(unkcnt, fromdeg, todeg)

    if unkcnt == 0 
        % by def
        k = 0;
        return;
    end

    k = 0;
    for i=fromdeg:todeg
        c = nchoosek(unkcnt+i-1, i);
        k = k + c;
    end
    
end