% cell array to string
%
% by Martin Bujnak, mar 2008

function [s] = c2s(a, separator)

    s = '';
    for i=1:length(a)

        if i==1
            s = char(a{i});
        else
            s = [s separator char(a{i})];
        end
    end
end