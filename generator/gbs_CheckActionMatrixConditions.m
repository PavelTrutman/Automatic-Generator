% Check if in given matrix are all required polynomials to build an action 
% matrix
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, May 2015

function [res, ordersReq] = gbs_CheckActionMatrixConditions(M, amStats, isElim, prime)

    if ~isElim

        nonzero = find(sum(M) ~= 0);
        Kk = M(:, nonzero);

        B = gjzpsp(Kk, prime);

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
        if nargout > 1
          ordersReq = [];
        end
        return;
    end
    
    % count which leading monials are not generated yet
    if nargout > 1

      % pick leading monomials
      leadingMons = [];
      for i = 1:size(M, 1)
        index = find(M(i, :), 1, 'first');
        leadingMons = [leadingMons, size(M, 2) - index + 1];
      end
      
      % pick variable with minimal required not present leading monomials
      ordersReq = setdiff(amStats{1}.amLT, leadingMons);
      for i = 2:length(amStats)
        req = setdiff(amStats{i}.amLT, leadingMons);
        if length(req) < length(ordersReq)
          ordersReq = req;
        end
      end
    end
    
    res = 0;
end