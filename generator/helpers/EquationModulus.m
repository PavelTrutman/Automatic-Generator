% Perform modulus on coefficients of given equation
% Pavel Trutman, trutmpav@fel.cvut.cz, aug2014
%
% [eq] = EquationModulus(eq, unknown, prime)
% Extracts coefficients of equation, perform modulus by prime and restore 
% the equation

function [eq] = EquationModulus(eq, unknown, prime)

  % use normal mod function, if only number given
  if isnumeric(eq)
    eq = mod(eq, prime);
    
  % if equation given
  else
    % use symbolic variables
    for u = unknown
      eval(['syms ' char(u) ';']);
    end
    % extract coefficients a monomials
    [coeff, monomials] = coeffs(eq, unknown);
    % mod the coefficients
    coeff = mod(coeff, prime);
    % restore the equation
    eq = sum(coeff.*monomials);
  end
end