% Build action natrix for variable "actMvar" given coefficient matrix M
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, February 2015


function [amrows, amcols, gjcols, aidx, PartitioningWorkflow] = gbs_BuildActionMatrix(cfg, M, algB, amLT, amLTall, algBidx, actMvar)

    % extract action matrix for variable "actMvar"
    
    cols = size(M, 2);
    algBcnt = size(algBidx, 2);

    amrows = zeros(algBcnt, 1);
    
    % GJ
    nonzero = find(sum(M) ~= 0);
    Kk = M(:, nonzero);
    B = gjzpsp(Kk, cfg.prime);
    A = zeros(size(M));
    A(:, nonzero) = B;
    
    %
    % detect leading '1'
    [val, idx] = min(abs(B' - 1));
    idx = idx(find(val == 0));
    
    % add algB columns

    % columns that are really required in GJ
    gjcols = union(nonzero(idx), (cols+1)-algBidx);
        
    % create monomials lookup table
    aidx = (-1)*ones(1,cols);
    aidx(gjcols) = 1:size(gjcols, 2);
    
    fprintf('Extracting the action matrix for variable ''%s'' (used coef. matrix size %dx%d) \n', char(actMvar), size(Kk, 1), length(gjcols));
    
    amcols = aidx((cols+1)-algBidx);
    
    if cfg.useMatrixPart
      %use matrix partitioning
      [B, PartitioningWorkflow] = gbs_MatrixPartitioning(M(:, gjcols), amcols, cfg.prime);
      A = zeros(size(M));
      A(:, gjcols) = B;
    else
      PartitioningWorkflow.enable = 0;
    end
    
    for i=1:algBcnt

        pos = find(algBidx == amLTall(i));
        if isempty(pos)
            
            % take remainder
            [r_one]=find(A(:,cols-amLTall(i)+1) == 1);
            amrows(i) = r_one;

        else
            % in the B
            amrows(i) = -pos;
        end
    end

end