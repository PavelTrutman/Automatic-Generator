% function [corrected] = Monomials_order(strin, varscnt, unkmon)
%
%  this function fixes and update order of monomials in input term
%  unknown monomials (singletons) must be of the form "x_1" ... "x_N"
%
% by Martin Bujnak

function [corrected] = Monomials_order(strin, varscnt, unkmon)

    if nargin < 2
        varscnt = 7;
    end
    
    if nargin < 3
        
        unkmon = 'x';
    end
    
    if size(unkmon, 2) > 1 & unkmon( size(unkmon, 2) - 1 ) == '_'
        
        %mutiple monomials with the same prefix...
        multi = true;
        unkmon = unkmon(1:end-2);
    else
        multi = false;
    end

    coefs = zeros(1, varscnt);
    coefscnt = size(coefs, 2);
    
    % extract coefs...and add sentinel 'x'
    str = [char(strin) unkmon];
    p__k = strfind(str, unkmon);
    
    for p__j=1:(size(p__k,2)-1)

        mon = str(p__k(p__j):p__k(p__j+1)-1);
        if (multi) 
            [p__a p__b p__c] = strread(mon, '%s%d%d', 'delimiter', '_^*');
        else
            [p__a p__c] = strread(mon, '%s%d', 'delimiter', '^*');
            p__b = 0;
        end
        p__b = p__b + 1;
        if (size(p__c,1) == 0)

            p__c = 1; 
        end
        coefs(coefscnt - p__b + 1) = p__c;
    end

    % composite
    corrected = '';
    first = true;
    for p__j=1:coefscnt
        
        if (coefs(p__j) > 0) 

            if (multi) 
            
                if (~first)
                    corrected = [corrected '*' unkmon '_' int2str(coefscnt - p__j) '^' int2str(coefs(p__j))];
                else
                    corrected = [corrected unkmon '_' int2str(coefscnt - p__j) '^' int2str(coefs(p__j))];
                    first = false;
                end
            else
                if (~first)
                    corrected = [corrected '*' unkmon '^' int2str(coefs(p__j))];
                else
                    corrected = [corrected unkmon '^' int2str(coefs(p__j))];
                    first = false;
                end
            end
        end
    end

end