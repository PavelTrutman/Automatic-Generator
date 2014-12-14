% calc monomial order (rank) in coefficient matrix ordered using grevlex.
%
% by Martin Bujnak, mar2008
%
% function [r] = GetMonomialOrder(mon, unknowns)
%   mon - monomial which order are we interested in (either string 
%         'x^2*y*z^2', or array form [2 1 2])
%   unknowns - set of unknown variables (define ordering) (ex. {'x' 'y' 'z'})


function [r] = GetMonomialOrder(mon, unknowns)

    if isstr(mon)
        
        [mons, degs] = ExtractMonomials(sym(mon), unknowns);
        mon = degs;
        
        if isempty(mons)
            r = 1;
            return;
        end
    end
    
    unkcnt = length(mon);
    deg = sum(mon);
    
    smaller = GetMonomialsCnt(unkcnt, 0, deg-1);
    
    % trim leading zeros (reduce unknowns)
    i = 1;
    while (i <= unkcnt) && (mon(i) == 0)
        i = i + 1;
    end
    head = i;

    % count number of monomials with the same degree
    r = smaller;
    
    remdeg = deg;
    for i=unkcnt:-1:1
        
        for j=(mon(i)+1):remdeg
            
            r = r + GetMonomialsCnt(i - 1, remdeg-j, remdeg-j);
        end
        remdeg = remdeg - mon(i);
        if remdeg < 1
            break;
        end
        
    end
    
    r = r + 1;
end

