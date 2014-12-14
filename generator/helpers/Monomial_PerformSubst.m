% switch monomial in general form/to one where monimials are
% of the form x_i ...
%
% by Martin Bujnak, mar2007

function [mon moncnt] = Monomial_PerformSubst(from, to, str)

    if (isempty(from)) 
        eval([' from = {' sprintf('''x_%d'' ',1:size(to, 2)) '};']);
    end

    if (isempty(to))
        eval([' to = {' sprintf('''x_%d'' ',1:size(from, 2)) '};']);
    end
    
    to = strrep(to, 'x', '#');
    
    for i=min(size(from, 2), size(to, 2)):-1:1
        
        f = char(from(i));
        t = char(to(i));
        
        str = strrep(str, f, t);
    end

    k = strfind(str, '#');
    moncnt = length(k);
    
    str = strrep(str, '#', 'x');
    mon = str;
end

