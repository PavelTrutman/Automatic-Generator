% Generate all polynomials required to build an action matrix. Polynomials
% are generated with using some strategies from F4 algorithm.
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, March 2015

function [foundVar, M, trace] = gbs_GeneratePolynomials_F4(p, eq, unknown, maxdeg, alldegs, allmonsdeg, allmons, amStats, cfg, algorithmCfg)

  prime = cfg.prime;
  Sel = algorithmCfg.Sel;
  
  maxorder = length(allmons);
  
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
        if(left && L{d}{j}.t == PSel{d}{i}.leftT && L{d}{j}.f == PSel{d}{j}.leftF)
          left = false;
        end
        if(right && L{d}{j}.t == PSel{d}{i}.rightT && L{d}{j}.f == PSel{d}{j}.rightF)
          right = false;
        end
        if(~left && ~right)
          break;
        end
      end
      
      if left
        last = last + 1;
        L{d}{last}.t = PSel{d}{i}.leftT;
        L{d}{last}.f = PSel{d}{i}.leftF;
      end
      if right
        last = last + 1;
        L{d}{last}.t = PSel{d}{i}.rightT;
        L{d}{last}.f = PSel{d}{i}.rightF;
      end
    end
    
    % reducto
    [Ftplus{d}, F{d}] = Reduction(L{d}, G, F);
    
    % insert new pairs
    for i = 1:size(Ftplus{d}, 1)
      [G, P] = Update(G, P, Ftplus{d}(i, :));
    end
    
    foundVar = gbs_CheckActionMatrixConditions(G, amStats, false, prime);
    
  end
  
end