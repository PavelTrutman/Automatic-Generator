% Generate and sort all monomials of given degree from a set of given unknowns
% 
% by Martin Bujnak, feb2008

function [mons degs] = GenerateMonomials(deg, unk, order)

    if (deg == 0) 
        
        mons = {};
        degs = [];

    elseif (deg == 1)
        
        mons = unk;
        degs = eye(length(unk));
    
    elseif length(unk) == 1
        
        mons = {[unk{1} '^' int2str(deg)]};
        degs = [deg];
        
    else

        mons = {};
        degs = [];
        for i=0:deg

            [monsi degsi] = GenerateMonomials(deg-i, unk(2:end));

            if i > 0
                if (i == 1) 
                    m = unk{1};
                else
                    m = [unk{1} '^' int2str(i)];
                end
                if isempty(monsi)
                    
                    monsi = {m};
                    degsi = [zeros(1, length(unk)-1)];
                else
                    for j=1:length(monsi)

                        monsi{j} = [m '*' monsi{j}];
                    end
                end

            end
            degsi  = [i*ones(size(degsi,1),1) degsi];
            
            mons = [mons monsi];
            degs = [degs; degsi];
        end
    end
    
    if nargin > 2
        
        [mons degs] = ReorderMonomials(mons, unk, order);
        
        % remove ^1
        for i=1:length(mons)
            
            mons{i} = strrep(mons{i}, '^1', '');
        end
    end
end