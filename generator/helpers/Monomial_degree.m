% Given string in monomial is of the form x_i1^j1*x_i2^j2*...*x_imax^jmax
% Return degree of the monomial deg = sum(ji) and vec = [j1 j2 ... jmax]. 
%
% by Martin Bujnak, oct 2007

function [deg vec] = Monomial_degree(monstr, max)

    deg = 0;
    
    if nargin < 2
        max = 8;
    end
    
    vec = [];
    
    % extract coefs...and add sentinel 'x'
    str = [char(monstr) 'x'];
    k = strfind(str, 'x');
    
    svec = zeros(1, max);
    for j=1:(size(k,2)-1)

        mon = str(k(j):k(j+1)-1);
        [a b c] = strread(mon, '%s%d%d', 'delimiter', '_^*');
        if (size(c,1) == 0)

            c = 1; 
        end
        
        svec(b) = svec(b) + c;
        deg = deg + c;
    end

    vec = [vec; svec];
end