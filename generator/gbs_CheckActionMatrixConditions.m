% Check if in given matrix are all required polynomials to build an action 
% matrix
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, March 2015

function [res] = gbs_CheckActionMatrixConditions(M, amStats, isElim, prime)

    if ~isElim

        nonzero = find(sum(M) ~= 0);
        Kk = M(:, nonzero);

        gjtime = cputime;
        B = gjzpsp(Kk, prime);
        gjtime = cputime - gjtime;
        if (gjtime > 0.5)
            fprintf('\t[t:%.2fsec]\t', gjtime);
        end

        M = zeros(size(M, 1), size(M, 2));
        M(:, nonzero) = B;
    end

    cols = size(M, 2);
        
    % check conditions for all action matrices
     for i=1:length(amStats)
    
        % test required leading terms
        if max(amStats{i}.amLT) > cols
            continue;
        end
        [r_ones, ~]=find(M(:, cols-amStats{i}.amLT+1) == 1);
        if size(r_ones,1) < amStats{i}.amLTcnt
            continue;
        end
        rowssum = find(sum(M(r_ones, amStats{i}.zero_el),2) == 1);
        if size(rowssum, 1) < amStats{i}.amLTcnt
            continue;
        end
        
        % all tests passed OK
        res = i;
        return;
    end
    
    res = 0;
end