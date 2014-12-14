% Generate all polynomials required to build an action matrix for
% variable "actMvar" (detect actMvar if variable was not specified)
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, December 2014


function [M, trace, symcoefs, amVar, amLT, amLTall, algBidx, algB] = gbs_PreparePolySystemCoefsMatrix(cfg, eq, known, unknown, algB, amVar)

    ordering = cfg.ordering;
    prime = cfg.prime;
    GJstep = cfg.GJstep;

    % create random instance
    if ~isfield(cfg, 'eqinstance')

        cfg.eqinstance = cfg.InstanceGenerator(cfg, eq, known, unknown);
    end
    
    if nargin < 5 || isempty(algB)
        
        % find quotion ring basis ( monomial basis of quotient ring C[x1,...xn]/I )
        [algB] = cfg.GBSolver(cfg, eq, known, unknown);
        
        if isempty(algB)
            
            M=[];
            trace=[];
            symcoefs=[];
            amVar=[];
            amLT=[];
            amLTall=[];
            algBidx=[];
            return;
        end
    end
    
    % prepare algebra B and find which monomials we need to build the
    % action matrix
    fprintf('analyzing quotient ring basis\n');
    
    if nargin < 9
        
        amVars = unknown;
    else
        amVars = amVar;
    end
    
    amStats = {};
    for i=1:length(amVars)
        
        actMvar = amVars{i};
        
        [amLT, amLTcnt, amLTall, algB, algBidx] = gbs_ParseAlgebraB(algB, actMvar, unknown);
        
        amStats{i}.actMvar = actMvar;
        amStats{i}.amLT = amLT;
        amStats{i}.amLTcnt = amLTcnt;
        amStats{i}.amLTall = amLTall;
        amStats{i}.algB = algB;
        amStats{i}.algBidx = algBidx;
    end

    %
    % extract equations into coeffients + monomial form
    fprintf('extracting coefficient & monomials\n');
    [monomials, p, ~, symcoefs, maxdeg] = gbs_ParseEquations(cfg, eq, known, unknown);
    
    
    %
    % just statistics (can be removed from this code)
    [monomials] = ReorderMonomials(monomials, unknown, ordering);
    
    fprintf('...used %d monomials : \n   ', length(monomials));
    for mons = monomials
        fprintf('%s ', char(mons));
    end
    fprintf('\n');

    fprintf('...max. poly. degree %d\n', maxdeg);

    
    if ~isempty(cfg.crashlog)

        fprintf('...crash log file detected - recovering from previous state using: %s\n', cfg.crashlog);
        [~, ~, max_deg] = gbs_ResetFromLog(cfg.crashlog);

        maxdeg = max_deg;
        fprintf('...need to generate polynomials to degree %d\n', maxdeg);
    end
    

    
    %
    % generate monomials up-to maxdeg (matrix cols)
    allmons = {'1'};
    alldegs = zeros(1, length(unknown));
    allmonsdeg = 0;
    for i=1:maxdeg
        [mons, degs] = GenerateMonomials(i, unknown, ordering);
        allmons = [mons allmons];
        alldegs = [degs; alldegs];
        allmonsdeg = [i*ones(1, length(mons)) allmonsdeg];
    end
    
    % create action var. masks
    cols = size(allmonsdeg, 2);
    for a=1:length(amStats)
        amStats{a}.zero_el = setdiff(1:cols, (cols+1)-amStats{a}.algBidx);
    end
    
    if cfg.exportEqs
        
        fprintf('...exporting parsed equations');
        
        % debug - create coefficient matrix with all initial polynomials
        cols = size(allmonsdeg, 2);
        Minit = zeros(length(eq), cols);
        for i=1:length(eq)
            for j=1:p{i}.monscnt
                pmon2 = p{i}.deg(j,:);
                order = GetMonomialOrder(pmon2, unknown); 
                Minit(i, cols-(order-1)) =  p{i}.coefsIDX(j);
            end        
        end
        
        nonzero = find(sum(Minit) ~= 0);
        Mmons = Minit(:, nonzero);
        
        gbs_DBGExport(known, unknown, eq, symcoefs, allmons, Minit, Mmons);
    end

    fprintf('Adding polynomials\n');
    fprintf('  Initializing matrix size %dx%d\n', length(eq), length(allmons));
    
    % count rows
    rowsalloc = 0;
    for i = 1:length(eq);
      for j = p{i}.maxdeg:maxdeg
        rowsalloc = rowsalloc + nchoosek(j - p{i}.maxdeg + length(unknown) - 1, length(unknown) - 1);
      end
    end
    
    % initial size
    cols = size(allmonsdeg, 2);
    M = zeros(rowsalloc, cols);
    Mcoefs  = zeros(rowsalloc, cols);
    
    % insert polynomials into matrices
    equations = zeros(1, length(eq));
    row = 1;
    for i = 1:length(p)
      equations(i) = row;
      for deg = 0:maxdeg - p{i}.maxdeg
        if deg == 0
          mons = '1';
          degs = zeros(1, length(unknown));
        else
          [mons, degs] = GenerateMonomials(deg, unknown, ordering);
        end
        
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
    
    iteration = 1;
    trace{iteration}.Mcoefs = Mcoefs;
    
    nonzero = find(sum(M) ~= 0);
    Kk = M(:, nonzero);
    B = gjzpsp(Kk, prime);
    MGJ = zeros(size(M, 1), size(M, 2));
    MGJ(:, nonzero) = B;
    var = gbs_CheckActionMatrixConditions(MGJ, amStats, true, prime);
    
    todeg = maxdeg + 1;
    if GJstep == 0
      nextElim = -1;
    else
      nextElim = 2;
    end
    
    if GJstep ~= 1
      MGJ = M;
      first = true;
    else
      first = false;
    end
    
    % generate polynomials
    while var == 0
      Mold = MGJ;
      if first
        MoldCoefs = Mcoefs;
        first = false;
      else
        MoldCoefs = zeros(size(Mold, 1), size(Mold, 2));
      end
      
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
          rows = rows + nchoosek(todeg - degsold(i) + length(unknown) - 1, length(unknown) - 1);
        end
        rows = rows + rowsold;
      
        % reallocate cols
        [mons, degs] = GenerateMonomials(todeg, unknown, ordering);
        allmons = [mons allmons];
        alldegs = [degs; alldegs];
        allmonsdeg = [todeg*ones(1, length(mons)) allmonsdeg];            
        cols = size(allmonsdeg, 2);
        for a=1:length(amStats)
          amStats{a}.zero_el = setdiff(1:cols, (cols+1)-amStats{a}.algBidx);
        end
        fprintf('  Coefficient matrix reallocaction, adding %d (%d degree) monomials and %d equations\n', length(mons), todeg, rows - rowsold);
        
        % prepare new matrices
        M = zeros(rows, cols);
        Mcoefs = zeros(rows, cols);
        M(rows - rowsold + 1:rows, cols - length(allmonsdegold) + 1:cols) = Mold(1:rowsold, :);
        Mcoefs(rows - rowsold + 1:rows, cols - length(allmonsdegold) + 1:cols) = MoldCoefs(1:rowsold, :);
        
        % generate new polynomials
        row = 1;
        for i = equations;
          [mons, degs] = GenerateMonomials(todeg - degsold(i), unknown, ordering);
          idx = find(Mold(i, :));
          for j = 1:length(mons)
            for k = 1:length(idx)
              monnew = alldegsold(idx(k), :) + degs(j, :);
              order = GetMonomialOrder(monnew, unknown);
              M(row, cols - (order - 1)) =  Mold(i, idx(k));
              if nextElim == 1
                Mcoefs(row, cols - (order - 1)) = size(Mold, 1)*(idx(k) - 1) + i;
              else
                Mcoefs(row, cols - (order - 1)) = MoldCoefs(i, idx(k));
              end
            end
            row = row + 1;
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
        
        if var ~= 0
          break;
        end;
        
        todeg = todeg + 1;
        if nextElim ~= -1
          nextElim = nextElim + 1;
        end
        
        Mold = M;
        MoldCoefs = Mcoefs;
        equations = equations + row - 1;
      end
      
      
      % remove not necesary polynomials
      toRemove = size(find(sum(MGJ, 2) == 0), 1);
      if toRemove > 0
        fprintf('    %d equations can be removed', toRemove);
        usedRows = size(find(sum(MGJ, 2) ~= 0), 1);
        removed = 0;
        step = max([floor(toRemove/4) 1]);
        up = 1;
        filter = 1:size(M, 1);
        while up <= size(M, 1)
          down = up + step - 1;
          if down > size(M, 1)
            down = size(M, 1);
            step = down - up + 1;
          end
          filterOld = filter;
          filter = setdiff(filter, up:down);
          
          Kk = M(filter, nonzero);
          B = gjzpsp(Kk, prime);
          
          if size(find(sum(B, 2) ~= 0), 1) < usedRows
            if step == 1
              up = up + 1;
            else
              step = max([floor(step/4) 1]);
            end
            filter = filterOld;
          else
            removed = removed + step;
            if removed == toRemove
              fprintf(' - all removed');
              break;
            end
            up = down + 1;
            step = min([2*step toRemove-removed]);
          end
          
        end
        fprintf('\n');
        
        M = M(filter, :);
        Mcoefs = Mcoefs(filter, :);
        nonzero = find(sum(M) ~= 0);
        Kk = M(:, nonzero);
        B = gjzpsp(Kk, prime);
        MGJ = zeros(size(M, 1), size(M, 2));
        MGJ(:, nonzero) = B;
      end
      
      
      fprintf('    GJ elimination performed\n');
      
      equations = find(sum(MGJ, 2) ~= 0)';
      
      % save trace
      trace{iteration}.Mcoefs = Mcoefs;
      trace{iteration}.nonzerocols = nonzero;
      if iteration > 1
        trace{iteration}.rowfrom = size(Mcoefs, 1) - trace{iteration}.rowsold + 1;
        trace{iteration}.rowto = size(Mcoefs, 1);
        trace{iteration}.columnfrom = size(Mcoefs, 2) - size(trace{iteration - 1}.Mcoefs, 2) + 1;
        trace{iteration}.columnto = size(Mcoefs, 2);
      end
      
      iteration = iteration + 1;
      if GJstep ~= 0
        nextElim = 1;
      end
      
    end

    foundVar = var;
    
    % remove not necesary polynomials
    fprintf('Removing not necesary polynomials\n');
    rows = size(M, 1);
    step = max([floor(rows/32) 1]);
    up = 1;
    filter = 1:rows;
    
    while up <= rows
      down = up + step - 1;
      if down > rows
        down = rows;
        step = down - up + 1;
      end
      
      if step > 1
        fprintf('  removing equation %d - %d ', up, down);
      else
        fprintf('  removing equation %d ', down);
      end
      
      filterOld = filter;
      filter = setdiff(filter, up:down);
      
      var = gbs_CheckActionMatrixConditions(M(filter, :), amStats, false, prime);
      
      if var == foundVar
        fprintf('succeeded\n');
        up = down + 1;
        step = 2*step;
      else
        fprintf('failed\n');
        if step == 1
          up = up + 1;
        else
          step = max([floor(step / 4) 1]);
        end
        filter = filterOld;
      end
      
    end
    
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
      
    else
      trace{iteration}.filter = filter(filter >= trace{iteration}.rowfrom) - trace{iteration}.rowfrom + 1;
      trace{iteration}.rowfrom = size(filter(filter < trace{iteration}.rowfrom), 2) + 1;
      trace{iteration}.rowto = size(filter, 2);
    end
    
    if ~foundVar
        
        M = [];
        trace = [];
        
        amVar = [];
        amLT = [];
        amLTall = [];
        algBidx = [];

    else
        
        amVar = amStats{foundVar}.actMvar;
        amLT = amStats{foundVar}.amLT;
        amLTall = amStats{foundVar}.amLTall;
        algBidx = amStats{foundVar}.algBidx;
    end
end

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

%{
function [M, Mcoefs] = gbs_ReallocMatrix(M, Mcoefs, torows, tocols)

    crows = size(M, 1);
    ccols = size(M, 2);
    
    if crows ~= torows
        
        M = [M; zeros(torows - crows, ccols)];
        Mcoefs = [Mcoefs; zeros(torows - crows, ccols)];
    end
    
    if ccols ~= tocols
               
        % shift coefs
        newcols = zeros(torows, tocols - ccols);
        M = [newcols M];
        Mcoefs = [newcols Mcoefs];
    end
end
%}

%{
function [c] = InvZp(x, p)

    [g,c] = gcd(x,-p);
    c = mod(c,p);
end
%}
