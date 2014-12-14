% create symbolic column vector
%
% [v] = gbs_Vector(format, dim, flag)
% use: 
%   format - variable name (format specifier %d optional)
%   dim - vector dimension
%   flag (optional) - 'real', 'unreal',... see 'sym' for the reference
%
% example:
%   [v] = gbs_Vector('a', 3)
%  
% v = 
%
%  a1
%  a2
%  a3
%
%   [v] = gbs_Vector('a_%d', 3)
%
% v =
%  
%  a_1
%  a_2
%  a_3
% 
% by Martin Bujnak, mar 2008

function [v] = gbs_Vector(format, dim, flag)

    if isempty(strfind(format, '%d'))
        format = [format '%d'];
    end
    v = transpose(sym(['[' sprintf([format ' '],1:dim) ']']));

    if nargin > 2
        for i=1:dim
            sym(v(i), flag);
        end
    end
    
end