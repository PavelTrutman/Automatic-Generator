% Returns list of all divisors of given monomial
% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015
%
% [divisors] = GetDivisors(monomial)
% monomial - array with degrees of monomial
% divisors - list of all divisors of given monomial

function [divisors] = GetDivisors(monomial)
  
  N = size(monomial, 2);
  possibillities = cell(1, N); 
  for i = 1:N
    possibillities{i} = [0:monomial(1, i)];
  end
  v = cell(N, 1);
  [v{:}] = ndgrid(possibillities{:});
  divisors = reshape(cat(N + 1, v{:}), [], N);
  
end
