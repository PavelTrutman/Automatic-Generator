% create substitution table of the form x^2 -> x2 for each monomial from the
% input set. Return all monomials in compact form x2, y2, y2x2, etc.
%
% by Martin Bujnak, nov2007

function [subst smonomials] = CreateSubstBlock(monomials)

    % prepare subst table
    subst = {};
    smonomials = {};
    for i=1:size(monomials, 2)

        str = monomials(1, i);
        str2 = strrep(str, '^', '');
        str2 = strrep(str2, '*', '');
        subst = [subst {[char(str) '=' char(str2)]}];
        smonomials = [smonomials {char(str2)}];
    end
    
end