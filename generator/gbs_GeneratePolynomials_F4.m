% Generate all polynomials required to build an action matrix. Polynomials
% are generated with using some strategies from F4 algorithm.

% last edit by Pavel Trutman, April 2015

function [foundVar, G, trace] = gbs_GeneratePolynomials_F4(p, eq, unknown, maxdeg, alldegs, allmonsdeg, allmons, amStats, cfg, algorithmCfg)
  
  global prime;
  global maxorder;
  global allDegs;
  global maxDeg;
  global unknowns;
  global allMons;
  global ordering;
  global GRefs;
  global input;
  
  prime = cfg.prime;
  Sel = algorithmCfg.Sel;
  maxorder = length(allmons);
  allDegs = alldegs;
  maxDeg = maxdeg;
  unknowns = unknown;
  allMons = allmons;
  ordering = cfg.ordering;
  input = p;
  
  d = 0;
  G = zeros(0, maxorder);
  P = cell(0, 1);
  GRefs = zeros(0, 2);
  GCoefs = zeros(0, maxorder);
  
  % make pairs
  for i = 1:length(p)
    f = zeros(1, maxorder);
    for j = 1:p{i}.monscnt
      order = GetMonomialOrder(p{i}.deg(j, :), unknowns);
      f(1, maxorder - order + 1) = p{i}.coefs(j);
    end
    [G, P, GRefs, GCoefs] = Update(G, P, f, GRefs, GCoefs, 0, i);
  end
  
  % itarate over pairs
  while length(P) ~= 0
    d = d + 1;
    F{d} = zeros(0, 0);
    Ft{d} = zeros(0, 0);
    
    % select pairs
    [PSel{d}, P] = Sel(P);
    
    % put left a right sides of pairs together
    L{d} = cell(0, 1);
    last = 0;
    for i = 1:length(PSel{d})
      
      left = true;
      right = true;
      % check duplicity of pairs
      for j = 1:length(L{d})
        if(left && size(L{d}{j}.polynomial, 2) == size(PSel{d}{i}.left.polynomial, 2) && sum(L{d}{j}.monomial ~= PSel{d}{i}.left.monomial) == 0 && sum(L{d}{j}.polynomial ~= PSel{d}{i}.left.polynomial) == 0)
          left = false;
        end
        if(right && size(L{d}{j}.polynomial, 2) == size(PSel{d}{i}.right.polynomial, 2) && sum(L{d}{j}.monomial ~= PSel{d}{i}.right.monomial) == 0 && sum(L{d}{j}.polynomial ~= PSel{d}{i}.right.polynomial) == 0)
          right = false;
        end
        if(~left && ~right)
          break;
        end
      end
      
      if left
        last = last + 1;
        L{d}{last} = PSel{d}{i}.left;
      end
      if right
        last = last + 1;
        L{d}{last} = PSel{d}{i}.right;
      end
    end
    
    % reducto
    [Ftplus{d}, F{d}, Ft{d}, FRefs, traceRefs, traceCoefs] = Reduction(L{d}, G, F, Ft);
    
    trace{d}.refs = traceRefs;
    trace{d}.coefs = traceCoefs;
    
    % insert new pairs
    for i = 1:size(Ftplus{d}, 1)
      [G, P, GRefs, GCoefs] = Update(G, P, Ftplus{d}(i, :), GRefs, GCoefs, d, FRefs(i, 1));
    end
    
    % recompute  amStats
    for a = 1:length(amStats)
      amStats{a}.zero_el = setdiff(1:size(G, 2), (size(G, 2) + 1) - amStats{a}.algBidx);
    end
    
    % check conditions for the action matrix
    foundVar = gbs_CheckActionMatrixConditions(G, amStats, false, prime);
    
    if foundVar
      break;
    end
    
  end

  trace{d + 1}.refs = GRefs;
  trace{d + 1}.coefs = GCoefs;
  
end



% Update (Buchberger)

function [GNew, BNew, GNewRefs, GNewCoefs] = Update(GOld, BOld, h, GOldRefs, GOldCoefs, hMatrixId, hPolRow)
  
  global maxorder;
  global allDegs;
  global input;
  
  % get HM of h
  hOrder = size(h, 2) - find(h, 1, 'first') + 1;
  hDeg = allDegs(end - hOrder + 1, :);
  
  if size(GOld, 2) < maxorder
    GOld = [zeros(size(GOld, 1), maxorder - size(GOld, 2)), GOld];
    GOldCoefs = [zeros(size(GOldCoefs, 1), maxorder - size(GOldCoefs, 2)), GOldCoefs];
  end
  
  D = zeros(0, maxorder);
  C = GOld;
  DRefs = zeros(0, 2);
  CRefs = GOldRefs;
  
  % go throught C
  for i = 1:size(C, 1)
    % get HM of g
    gOrder = size(C, 2) - find(C(i, :), 1, 'first') + 1;
    gDeg = allDegs(end - gOrder + 1, :);
    
    % are HM(h) and HM(C(i)) disjoint?
    if sum(min([hDeg; gDeg], [], 1)) == 0
      % are disjoint
      D = [D; C(i, :)];
      DRefs = [DRefs; CRefs(i, :)];
    else
      lcmHG = max([hDeg; gDeg], [], 1);
      condition1 = true;
      for j = i + 1:size(C, 1)
        %get HM of g2
        g2Order = size(C, 2) - find(C(j, :), 1, 'first') + 1;
        g2Deg = allDegs(end - g2Order + 1, :);
        lcmHG2 = max([hDeg; g2Deg], [], 1);
        % check first condition of divisibility
        if sum(lcmHG2 <= lcmHG) == size(lcmHG, 2)
          condition1 = false;
          break;
        end
      end
      
      if condition1
        condition2 = true;
        for j = 1:size(D, 1)
          g2Order = size(D, 2) - find(D(j, :), 1, 'first') + 1;
          g2Deg = allDegs(end - g2Order + 1, :);
          lcmHG2 = max([hDeg; g2Deg], [], 1);
          % check second condition of divisibility
          if sum(lcmHG2 <= lcmHG) == size(lcmHG, 2)
            condition2 = false;
            break;
          end
        end
        
        if condition2
          % all conditions satisfied
          D = [D; C(i, :)];
          DRefs = [DRefs; CRefs(i, :)];
        end
      end
      
    end
  end
  
  E = cell(0, 1);
  last = 0;
  % go throught D
  for i = 1:size(D, 1)
    % get HM of g
    gOrder = size(D, 2) - find(D(i, :), 1, 'first') + 1;
    gDeg = allDegs(end - gOrder + 1, :);
    if sum(min([hDeg; gDeg], [], 1)) ~= 0
      % are not disjoint
      lcmHG = max([hDeg; gDeg], [], 1);
      % create pair
      last = last + 1;
      E{last, 1}.left.polynomial = h;
      E{last, 1}.left.monomial = lcmHG - hDeg;
      E{last, 1}.left.matrixId = hMatrixId;
      E{last, 1}.left.polRow = hPolRow;
      E{last, 1}.right.polynomial = D(i, :);
      E{last, 1}.right.monomial = lcmHG - gDeg;
      E{last, 1}.right.matrixId = DRefs(i, 1);
      E{last, 1}.right.polRow = DRefs(i, 2);
      E{last, 1}.lcm = lcmHG;
    end
  end
  
  BNew = cell(0, 1);
  last = 0;
  % go thorought BOld
  for i = 1:length(BOld)
    % get HM of g1 and g2
    g1Order = size(BOld{i}.left.polynomial, 2) - find(BOld{i}.left.polynomial, 1, 'first') + 1;
    g1Deg = allDegs(end - g1Order + 1, :);
    g2Order = size(BOld{i}.right.polynomial, 2) - find(BOld{i}.right.polynomial, 1, 'first') + 1;
    g2Deg = allDegs(end - g2Order + 1, :);

    if (sum(hDeg < BOld{i}.lcm) < size(hDeg, 2)) || (sum(max([hDeg; g1Deg], [], 1) ~= BOld{i}.lcm) == 0) || (sum(max([hDeg; g2Deg], [], 1) ~= BOld{i}.lcm) == 0)
      % add pair
      last = last + 1;
      BNew{last, 1} = BOld{i};
    end
  end
  
  % add E into BNew
  BNew = vertcat(BNew, E);
  
  GNew = zeros(0, maxorder);
  GNewRefs = zeros(0, 2);
  GNewCoefs = zeros(0, maxorder);
  last = 0;
  % go thorought GOld
  for i = 1:size(GOld, 1)
    % get HM of g
    gOrder = size(GOld, 2) - find(GOld(i, :), 1, 'first') + 1;
    gDeg = allDegs(end - gOrder + 1, :);
    
    if sum(hDeg <= gDeg) < size(hDeg, 2)
      last = last + 1;
      GNew(last, :) = GOld(i, :);
      GNewRefs(last, :) = GOldRefs(i, :);
      GNewCoefs(last, :) = GOldCoefs(i, :);
    end
  end
  
  % add h into GNew
  GNew(end + 1, :) = h;
  GNewRefs(end + 1, :) = [hMatrixId, hPolRow];
  
  last = size(GNewCoefs, 1) + 1;
  orders = size(h, 2) - find(h) + 1;
  for order = orders
    deg = allDegs(maxorder - order + 1, :);
    if hMatrixId ~= 0
      % input polynomial from some matrix F
      GNewCoefs(last, maxorder - order + 1) = order;
    else
      % input polynomial from the input set
      GNewCoefs(last, maxorder - order + 1) = input{hPolRow}.coefsIDX(find(ismember(input{hPolRow}.deg, deg, 'rows'), 1, 'first'));
    end
  end
  
end



% Reduction (F4)
  
function [Ftplus, F, Ft, FtRefs, traceRefs, traceCoefs] = Reduction(L, G, FAll, FtAll)
   
  global prime;
  global maxorder;
 
  [F, traceRefs, traceCoefs] = SymbolicPreprocessing(L, G, FAll, FtAll);
  
  nonzero = find(sum(F) ~= 0);
  Kk = F(:, nonzero);
  fprintf('  Reducing matrix %dx%d\n', size(Kk, 1), size(Kk, 2));
  B = gjzpsp(Kk, prime);
  Ft = zeros(size(F, 1), size(F, 2));
  Ft(:, nonzero) = B;
  
  % pick head monomials
  HM = zeros(0, 1);
  for i = 1:size(F, 1)
    index = find(F(i, :), 1, 'first');
    if ~isempty(index)
      HM = [HM; index];
    end
  end
  
  % pick rows with new head monomials
  last = 0;
  Ftplus = zeros(0, maxorder);
  FtRefs = zeros(0, 1);
  for i = 1:size(F, 1)
    index =find(Ft(i, :), 1, 'first');
    if ~isempty(index)
      if sum(HM == index) == 0
        last = last + 1;
        Ftplus(last, :) = Ft(i, :);
        FtRefs(last, 1) = i;
      end
    end
  end
  
end



% Symbolic Preprocessing (F4)

function [F, traceRefs, traceCoefs] = SymbolicPreprocessing(L, G, FAll, FtAll)
  
  global maxorder;
  global allDegs;
  global GRefs;

  % multiply polynomials from pairs and parse used monomials
  monomials = zeros(0, 1);
  headMonomials = zeros(0, 1);
  F = zeros(length(L), maxorder);
  traceRefs = zeros(length(L), 2);
  traceCoefs = zeros(length(L), maxorder);
  for i = 1:length(L)
    [monomial, polynomial, matrixId, polRow] = Simplify(L{i}.monomial, L{i}.polynomial, FAll, FtAll, L{i}.matrixId, L{i}.polRow);
    [f, coefs] = Multiply(monomial, polynomial, [matrixId, polRow]);
    if size(f, 2) > size(F, 2)
      F = [zeros(size(F, 1), size(f, 2) - size(F, 2)), F];
      traceCoefs = [zeros(size(traceCoefs, 1), size(coefs, 2) - size(traceCoefs, 2)), traceCoefs];
    end
    F(i, :) = f;
    traceRefs(i, :) = [matrixId, polRow];
    traceCoefs(i, :) = coefs;
    indices = find(F(i, :));
    headMonomials = unique([headMonomials; size(F, 2) - indices(1) + 1]);
    monomials = unique([monomials; size(F, 2) - indices(2:end)' + 1]);
  end
  
  % get head monomials of all polynomials from G
  headMonsGDegs = zeros(size(G, 1), size(allDegs, 2));
  for i = 1:size(G, 1)
    index = find(G(i, :), 1, 'first');
    headMonsGDegs(i, :) = allDegs(end - size(G, 2) + index, :);
  end
  
  done = headMonomials;
  remaining = setdiff(monomials, headMonomials);
  
  % go throught all monomials
  while isempty(remaining) == 0
    m = remaining(end);
    remaining = setdiff(remaining, m);
    done = [done; m];
    
    for i = 1:size(headMonsGDegs, 1)
      if sum((allDegs(end - m + 1, :) - headMonsGDegs(i, :)) < 0) == 0
        [monomial, polynomial, matrixId, polRow] = Simplify(allDegs(end - m + 1, :) - headMonsGDegs(i, :), G(i, :), FAll, FtAll, GRefs(i, 1), GRefs(i, 2));
        [f, coefs] = Multiply(monomial, polynomial, [matrixId, polRow]);
        F = [F; f];
        traceRefs = [traceRefs; matrixId, polRow];
        traceCoefs = [traceCoefs; coefs];
        indices = find(F(end, :));
        mons = size(F, 2) - indices' + 1;
        remaining = unique([remaining; setdiff(mons, done)]);
        break;
      end
    end
    
  end
  
end



% Simplify (F4)

function [monomial, polynomial, matrixIdNew, polRowNew] = Simplify(m, f, FAll, FtAll, matrixId, polRow)
  
  % get divisors of m
  u = GetDivisors(m);
  u = setdiff(u, zeros(size(m)), 'rows');
  
  [~, IX] = sort(sum(u, 2), 'ascend');
  u = u(IX, :);
  
  % go trought all divisors of m
  for i = 1:size(u, 1)
    uf = Multiply(u(i, :), f);
    
    % go thorught FAll
    for j = 1:length(FAll) - 1
      %if exists u*f in FAll{j}
      if sum(sum(([zeros(size(FAll{j}, 1), size(uf, 2)-size(FAll{j}, 2)) FAll{j}] - ones(size(FAll{j}, 1), 1)*uf) ~= 0, 2) == 0)
        HMuf = size(uf, 2) - find(uf, 1, 'first') + 1;
        
        % find HM of FtAll{j} such equals to HM(u*f)
        for k = 1:size(FtAll{j}, 1)
          if (size(FtAll{j}, 2) - find(FtAll{j}(k, :), 1, 'first') + 1) == HMuf
            
            if sum(u(i, :) ~= m) ~= 0
              [monomial, polynomial, matrixIdNew, polRowNew] = Simplify(m - u(i, :), FtAll{j}(k, :), FAll, FtAll, j, k);
              return;
            else
              monomial = zeros(size(m));
              polynomial = FtAll{j}(k, :);
              matrixIdNew = j;
              polRowNew = k;
              return;
            end
            
          end
        end
      end
    end
  end
  
  matrixIdNew = matrixId;
  polRowNew = polRow;
  monomial = m;
  polynomial = f;
  
end



% Multiply
% Multiplies monomial m and polynomial f

function [polynomial, coefs] = Multiply(m, f, source)
  
  global maxorder;
  global maxDeg;
  global unknowns;
  global allDegs;
  global allMons;
  global ordering;
  global input;
  
  if nargin < 3
    source = [-1 -1];
  end
  
  polynomial = zeros(1, maxorder);
  coefs = zeros(1, maxorder);
  % multiply each monomial of f
  orders = size(f, 2) - find(f) + 1;
  for order = orders
    deg = allDegs(maxorder - order + 1, :);
    newOrder = GetMonomialOrder(deg + m, unknowns);
    
    if  newOrder > maxorder
      % enlarge matrix dimensions
      while newOrder > maxorder
        maxDeg = maxDeg + 1;
        [mons, degs] = GenerateMonomials(maxDeg, unknowns, ordering);
        allDegs = [degs; allDegs];
        allMons = [mons, allMons];
        maxorder = size(allMons, 2);
      end
      polynomial = [zeros(1, maxorder - size(polynomial, 2)), polynomial];
      coefs = [zeros(1, maxorder - size(coefs, 2)), coefs];
    end
    
    % numbers
    polynomial(1, maxorder - newOrder + 1) = f(1, size(f, 2) - order + 1);
    
    % coefs
    if source(1, 1) ~= 0
      % from other matrix F
      coefs(1, maxorder - newOrder + 1) = order;
    else
      % from input polynomial
      coefs(1, maxorder - newOrder + 1) = input{source(1, 2)}.coefsIDX(find(ismember(input{source(1, 2)}.deg, deg, 'rows'), 1, 'first'));
    end
  end
  
end
