% Parse and prepare algebra B basis
% (GBsolver subroutine)
% by Martin Bujnak, mar2008


function [amLT, amLTcnt, amLTall, algB, algBidx, algBcnt] = gbs_ParseAlgebraB(algB, actMvar, unknown)

    algBcnt = size(algB, 2);
    algBidx = zeros(1, algBcnt);
    tmpLt = zeros(1, algBcnt);
    mvar = sym(actMvar);
    for i=1:algBcnt

        algBidx(i) = GetMonomialOrder(algB{i}, unknown); 

        ai = sym(algB{i}) * mvar;
        tmpLt(i) = GetMonomialOrder(char(ai), unknown); 
    end

    % sort (grevlex)
    [algBidx, p] = sort(algBidx);
    tmpLt = tmpLt(p);
    algB = algB(p);

    amLT = setdiff(tmpLt, algBidx);
    amLTcnt = size(amLT, 2);
    amLTall = tmpLt;
end