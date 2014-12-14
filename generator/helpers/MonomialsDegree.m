% Gets the degrees of monomials
% Pavel Trutman, pavel.trutman@fel.cvut.cz, oct2014
%
% degree = MonomialsDegree(monomial, unknown)
% Degree is an array, rows for monomials, colums for unknowns

function degree = MonomialsDegree(monomial, unknown)
  degree = zeros(length(monomial), length(unknown));
  
  for i = 1:length(monomial)
    for j = 1:length(unknown)
      
      % get the exponent
      deg = regexp(monomial{i}, [unknown{j} '\^(\d*)'], 'tokens');
      if isempty(deg) == 1
        deg = regexp(monomial{i}, unknown{j}, 'tokens');
        if isempty(deg) == 1
          % unknown{j} is not in this monomial
          degree(i, j) = 0;
        else
          % unknown{j} only with exponent 1
          degree(i, j) = 1;
        end
      else
        % higher exponents of unknown{j}
        degree(i, j) = str2num(char(deg{1}));
      end
      
    end
  end
  
end