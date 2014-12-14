% create symbolic 'r' by 'c' matrix
%
%
% [A] = gbs_Matrix(format, r, c, flag)
% use: 
%   format - variable name (format specifier %d optional)
%   r - number of rows
%   c - number of cols
%   flag (optional) - 'real', 'unreal',... see 'sym' for the reference
%
% example:
%   [A] = gbs_Matrix('a', 2, 3)
%
%  A =
%  
%  [ a11, a12, a13]
%  [ a21, a22, a23]
%
%   [A] = gbs_Matrix('a(%d,%d)', 2, 3)
%
%  A =
%  
%  [ a(1,1), a(1,2), a(1,3)]
%  [ a(2,1), a(2,2), a(2,3)]
% 
% by Martin Bujnak, mar 2008

function [A] = gbs_Matrix(format, r, c, flag)

    cnt = length(strfind(format, '%d'));
    switch cnt 
        case 0
            format = [format '%d%d'];
        case 1
            format = [format '%d'];
    end

    idx1 = repmat((1:r)', 1, c);
    idx2 = repmat(1:c, r, 1);
    M=[idx1(:)'; idx2(:)'];
    
    v = transpose(sym(['[' sprintf([format ' '],M(:)) ']']));
    A = reshape(v, r, c);

    if nargin > 3
        for i=1:(r*c)
            sym(A(i), flag);
        end
    end
    
end