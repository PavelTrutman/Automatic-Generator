% get pointers of unknown variables in algB array
% (GBsolver subroutine)
% by Martin Bujnak, dec2008

function [one unksidx] = gbs_GetVariablesIdx(algB, unkn)

    ucnt = size(unkn, 2);

    one = findmon(algB, '1');
    
    unksidx = zeros(1, ucnt);
    for i=1:ucnt
        unksidx(i) = findmon(algB, unkn{i});
    end
end

function [idx] = findmon(algB, ch)

    res = strcmp(ch, algB);
    idx = find(res == 1);
    if (isempty(idx))
        idx = 0;
    end
end