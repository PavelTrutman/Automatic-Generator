% Generate all polynomials required to build an action matrix. Polynomials
% are generated without any strategy.
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, May 2015

function [foundVar, M, trace] = gbs_GeneratePolynomials_systematic(p, eq, unknown, maxdeg, alldegs, allmonsdeg, allmons, amStats, cfg, algorithmCfg)

  prime = cfg.prime;
  GJstep = algorithmCfg.GJstep;
  ordering = cfg.ordering;

  % count rows
  rowsalloc = 0;
  for i = 1:length(eq);
    for j = p{i}.maxdeg+1:maxdeg-1
      rowsalloc = rowsalloc + nchoosek(j - p{i}.maxdeg + length(unknown) - 1, length(unknown) - 1);
    end
  end
  rowsalloc = rowsalloc + length(eq);
    
  % initial size
  cols = size(allmonsdeg, 2);
  M = zeros(rowsalloc, cols);
  Mcoefs  = zeros(rowsalloc, cols);
  
  fprintf('  Initializing matrix size %dx%d\n', rowsalloc, cols);
  
  % insert polynomials into matrices
  equations = rowsalloc - length(p) + 1:rowsalloc;
  row = 1;
  for i = 1:length(p)
    for deg = 1:maxdeg - p{i}.maxdeg - 1
      [mons, degs] = GenerateMonomials(deg, unknown, ordering);
      for j = 1:length(mons)
        for k = 1:p{i}.monscnt
          monnew = alldegs(length(allmonsdeg) - GetMonomialOrder(char(p{i}.mons(k)), unknown) + 1, :) + degs(j, :);
          order = GetMonomialOrder(monnew, unknown);
          M(row, cols - (order - 1)) =  p{i}.coefs(k);
          Mcoefs(row, cols - (order - 1)) = p{i}.coefsIDX(k);
        end
        row = row + 1;
      end
    end
  end
  for i = 1:length(p)
    for k = 1:p{i}.monscnt
      monnew = alldegs(length(allmonsdeg) - GetMonomialOrder(char(p{i}.mons(k)), unknown) + 1, :);
      order = GetMonomialOrder(monnew, unknown);
      M(row, cols - (order - 1)) =  p{i}.coefs(k);
      Mcoefs(row, cols - (order - 1)) = p{i}.coefsIDX(k);
    end
    row = row + 1;
  end
  
  iteration = 1;
  trace{iteration}.Mcoefs = Mcoefs;
  
  nonzero = find(sum(M) ~= 0);
  Kk = M(:, nonzero);
  B = gjzpsp(Kk, prime);
  MGJ = zeros(size(M, 1), size(M, 2));
  MGJ(:, nonzero) = B;
  var = gbs_CheckActionMatrixConditions(MGJ, amStats, true, prime);
  
  todeg = maxdeg;
  if GJstep == 0
    nextElim = -1;
  else
    nextElim = 1;
  end
  
  first = true;
  equationsAddedOld = 0;
  
  % generate polynomials
  while var == 0
    if first
      MoldCoefs = Mcoefs;
      Mold = M;
      first = false;
    else
      MoldCoefs = zeros(size(Mold, 1), size(Mold, 2));
      Mold = MGJ;
    end
    
    todegOld = todeg;
    
    while nextElim <= GJstep
      alldegsold = alldegs;
      allmonsdegold = allmonsdeg;
      
      %count new rows
      rows = 0;
      degsold = zeros(1, size(Mold, 1));
      rowsold = find(sum(Mold, 2) ~= 0, 1, 'last');
      for i = equations
        [~, id] = find(Mold(i, :), 1, 'first');
        degsold(i) = allmonsdegold(id);
        if degsold(i) < todeg
          rows = rows + nchoosek(todeg - degsold(i) + length(unknown) - 1, length(unknown) - 1);
        end
      end
      rows = rows + rowsold;
      
      % reallocate cols
      if todeg > maxdeg
        maxdeg = todeg;
        [mons, degs] = GenerateMonomials(todeg, unknown, ordering);
        allmons = [mons allmons];
        alldegs = [degs; alldegs];
        allmonsdeg = [todeg*ones(1, length(mons)) allmonsdeg];
        cols = size(allmonsdeg, 2);
        for a=1:length(amStats)
          amStats{a}.zero_el = setdiff(1:cols, (cols+1)-amStats{a}.algBidx);
        end
        fprintf('  Coefficient matrix reallocaction, adding %d (%d degree) monomials\n', length(mons), todeg);
      end
      fprintf('  %d equations added\n', rows - rowsold);
      
      % prepare new matrices
      M = zeros(rows, cols);
      Mcoefs = zeros(rows, cols);
      M(rows - rowsold + 1:rows, cols - length(allmonsdegold) + 1:cols) = Mold(1:rowsold, :);
      Mcoefs(rows - rowsold + 1:rows, cols - length(allmonsdegold) + 1:cols) = MoldCoefs(1:rowsold, :);
      
      % generate new polynomials
      row = 1;
      for i = equations;
        if degsold(i) < todeg
          [mons, degs] = GenerateMonomials(todeg - degsold(i), unknown, ordering);
          idx = find(Mold(i, :));
          for j = 1:length(mons)
            for k = 1:length(idx)
              monnew = alldegsold(idx(k), :) + degs(j, :);
              order = GetMonomialOrder(monnew, unknown);
              M(row, cols - (order - 1)) =  Mold(i, idx(k));
              if (nextElim == 1) && (iteration ~= 1)
                Mcoefs(row, cols - (order - 1)) = size(Mold, 1)*(idx(k) - 1) + i;
              else
                Mcoefs(row, cols - (order - 1)) = MoldCoefs(i, idx(k));
              end
            end
            row = row + 1;
          end
        end
      end
      
      nonzero = find(sum(M) ~= 0);
      Kk = M(:, nonzero);
      B = gjzpsp(Kk, prime);
      MGJ = zeros(size(M, 1), size(M, 2));
      MGJ(:, nonzero) = B;
      var = gbs_CheckActionMatrixConditions(MGJ, amStats, true, prime);
      
      if nextElim == 1
        trace{iteration}.rowsold = rowsold;
      end
      
      todeg = todeg + 1;
      
      if var ~= 0
        break;
      end;
      
      if nextElim ~= -1
        nextElim = nextElim + 1;
      end
      
      Mold = M;
      MoldCoefs = Mcoefs;
      equations = equations + row - 1;
    end
    
    
    % remove not necesary polynomials
    filter = gbs_RemoveRedundant(M, prime);
 
    M = M(filter, :);
    Mcoefs = Mcoefs(filter, :);
    nonzero = find(sum(M) ~= 0);
    Kk = M(:, nonzero);
    B = gjzpsp(Kk, prime);
    MGJ = zeros(size(M, 1), size(M, 2));
    MGJ(:, nonzero) = B;
    
    fprintf('    GJ elimination performed on matrix with size %dx%d\n', size(M, 1), length(nonzero));
    
    equationsAdded = size(find(sum(MGJ, 2) ~= 0), 1);
    if equationsAddedOld >= equationsAdded
      %revert back to previous iteration and leave todeg increased
      iteration = iteration - 1;
      Mcoefs = trace{iteration}.Mcoefs;
      MGJ = trace{iteration}.MGJ;
      discard = 1;
    else
      %try to add new equations up to total degree todeg
      todeg = todegOld;
      discard = 0;
    end
    equationsAddedOld = equationsAdded;
    
    equations = find(sum(MGJ, 2) ~= 0)';
    
    %matrix partitioning for elimination
    if strcmp(cfg.matrixPartitioning, 'all')
      [~, partitioning] = gbs_MatrixPartitioning(M(:, nonzero), [], true, cfg.prime);
      partitioning.enable = 1;
    else
      partitioning.enable = 0;
    end
    
    % save trace
    if (discard == 0) || (var ~= 0)
      trace{iteration}.Mcoefs = Mcoefs;
      trace{iteration}.MGJ = MGJ;
      trace{iteration}.nonzerocols = nonzero;
      trace{iteration}.partitioning = partitioning;
      trace{iteration}.size = size(Mcoefs);
      if iteration > 1
        trace{iteration}.rowfrom = size(Mcoefs, 1) - trace{iteration}.rowsold + 1;
        trace{iteration}.rowto = size(Mcoefs, 1);
        trace{iteration}.columnfrom = size(Mcoefs, 2) - size(trace{iteration - 1}.Mcoefs, 2) + 1;
        trace{iteration}.columnto = size(Mcoefs, 2);
      end
    end
    
    iteration = iteration + 1;
    if GJstep ~= 0
      nextElim = 1;
    end
    
  end
  
  foundVar = var;
  
  % remove not necesary polynomials
  filter = gbs_RemoveUnnecessary(M, amStats, foundVar, prime);
  
  if iteration ~= 1
    iteration = iteration - 1;
  end
  
  M = M(filter, :);
  trace{iteration}.Mcoefs = trace{iteration}.Mcoefs(filter, :);
  if foundVar ~= 0
    trace{iteration}.nonzerocols = union(find(sum(M) ~= 0), size(M, 2) - amStats{foundVar}.algBidx + 1);
  else
    trace{iteration}.nonzerocols = find(sum(M) ~= 0);
  end
  if iteration == 1
    %matrix partitioning for elimination
    if strcmp(cfg.matrixPartitioning, 'all')
      [~, partitioning] = gbs_MatrixPartitioning(M(:, nonzero), [], true, cfg.prime);
      partitioning.enable = 1;
    else
      partitioning.enable = 0;
    end
    trace{iteration}.partitioning = partitioning;
  else
    trace{iteration}.filter = filter(filter >= trace{iteration}.rowfrom) - trace{iteration}.rowfrom + 1;
    trace{iteration}.rowfrom = size(filter(filter < trace{iteration}.rowfrom), 2) + 1;
    trace{iteration}.rowto = size(filter, 2);
    trace{iteration}.size = size(trace{iteration}.Mcoefs);
  end
  
end
