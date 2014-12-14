% function compares two arbitary vectors and return 
%   -1 if a < b
%    0 if a = b
%   +1 if a > b
%
% by Martin Bujnak, nov2007

function [c] = Monomial_grevlexorder(a, b)

    sa = sum(a);
    sb = sum(b);
    
    if sa < sb
        
        c = -1;
        return;
        
    elseif sa > sb
        
        c = 1;
        return
    end

    a = a - b;
    l = length(a);

    for i=l:-1:1
        
        if a(i) > 0
          
            c = -1;
            return;
            
        elseif a(i) < 0

            c = 1;
            return;
        end
    end

    c = 0;
end