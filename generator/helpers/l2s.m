% integer array to string
%
% by Martin Bujnak, mar 2008

function [s] = l2s(a, separator)

    s = '';
    for i=1:length(a)
    
        if i==1
            s = int2str(a(i));
        else
            s = [s separator int2str(a(i))];
        end
    end
end