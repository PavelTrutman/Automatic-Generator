% Generate all polynomials required to build an action matrix. Polynomials
% are generated with using some strategies from F4 algorithm.
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, March 2015

function [foundVar, M, trace] = gbs_GeneratePolynomials_F4(p, eq, unknown, maxdeg, alldegs, allmonsdeg, allmons, amStats, cfg, algorithmCfg)
  
  global prime;
  global maxorder;
  global allDegs;
  
  prime = cfg.prime;
  Sel = algorithmCfg.Sel;
  maxorder = length(allmons);
  allDegs = alldegs;
  
  d = 0;
  G = zeros(0, maxorder);
  P = cell(0, 1);
  
  % make pairs
  for i = 1:length(p)
    f = zeros(1, maxorder);
    for j = 1:p{i}.monscnt
      order = GetMonomialOrder(p{i}.deg(j, :), unknown);
      f(1, maxorder - order + 1) = p{i}.coefs(j);
      [G, P] = Update(G, P, f);
    end
  end
  
  % itarate over pairs
  while lenght(P) ~= 0
    d = d + 1;
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
        if(left && L{d}{j}.monomial == PSel{d}{i}.left.monomial && L{d}{j}.polynomial == PSel{d}{j}.left.polynomial)
          left = false;
        end
        if(right && L{d}{j}.monomial == PSel{d}{i}.right.monomial && L{d}{j}.polynomial == PSel{d}{j}.right.polynomial)
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
    [Ftplus{d}, F{d}, Ft{d}] = Reduction(L{d}, G, F, Ft);
    
    % insert new pairs
    for i = 1:size(Ftplus{d}, 1)
      [G, P] = Update(G, P, Ftplus{d}(i, :));
    end
    
    foundVar = gbs_CheckActionMatrixConditions(G, amStats, false, prime);
    
  end
  
end

  
  
% Reduction (F4)
  
function [Ftplus, F, Ft] = Reduction(L, G, FAll, FtAll)
   
  global prime;
  global maxorder;
 
  F = SymbolicPreprocessing(L, G, FAll, FtAll);
  Ft = gjzpsp(F, prime);
  
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
  for i = 1:size(F, 1)
    index =find(Ft(i, :), 1, 'first');
    if ~isempty(index)
      if sum(HM == index) == 0
        last = last + 1;
        Ftplus(last, :) = Ft(i, :);
      end
    end
  end
  
end



% Symbolic Preprocessing (F4)

function [F] = SymbolicPreprocessing(L, G, FAll, FtAll)
  
  global maxorder;
  global allDegs;

  % multiply polynomials from pairs and parse used monomials
  monomials = zeros(0, 1);
  headMonomials = zeros(0, 1);
  F = zeros(length(L), maxorder);
  for i = 1:length(L)
    F(i, :) = Multiply(Simplify(L{i}.monomial, L{i}.polynomial, FAll, FtAll));
    indices = find(F(i, :));
    headMonomials = unique([headMonomials; size(F, 2) - indices(1) + 1]);
    monomials = unique([monomials; size(F, 2) - indices(2:end)' + 1]);
  end
  
  % get head monomials of all polynomials from G
  headMonsGDegs = zeros(size(G, 1), size(allDegs, 2));
  for i = 1:size(G, 1)
    index = find(G(i, :), 1, 'first');
    headMonsGDegs(i, :) = allDegs(size(G, 2) - index + 1, :);
  end
  
  done = headMonomials;
  remaining = set1diff(monomials, headMonomials);
  
  % go throught all monomials
  while isempty(remaining) ~= 0
    m = remaining(end);
    remaining = setdiff(remaining, m);
    done = [done; m];
    
    for i = 1:size(headMonsGDegs, 1)
      if sum((allDegs(m, :) - headMonsGDegs(i, :)) < 0) == 0
        F = [F; Multiply(Simplify(allDegs(m, :) - headMonsGDegs(i, :), G(i, :), FAll. FtAll))];
        indices = find(F(end, :));
        mons = size(F, 2) - indices' + 1;
        remaining = unique([remaining; setdiff(mons, done)]);
        break;
      end
    end
    
  end
  
end
