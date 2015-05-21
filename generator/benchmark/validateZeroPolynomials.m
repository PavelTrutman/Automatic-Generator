% Pavel Trutman, pavel.trutman@fel.cvut.cz, May 2015
%
% Validation function used by benchmark utility.
% This function substitutes computed solution into polynomials generated 
% by solver. This should be zero value, but because of bad stability we 
% get non zero value.
%
% inputData - a set of instanced known parameters
% solution - computed solutions by solver
% eq - input equations
% unknown - list of unknowns
% known - list of knowns

function [err] = validateZeroPolynomials(inputData, ~, solution, eq, unknown, known)
  
  err = ones(size(solution, 1), 0)*NaN;
  
  % make symbolic variables of known and unknown
  for i = 1:length(unknown)
    eval(['syms ', unknown{i}]);
  end
  for i = 1:length(known)
    eval(['syms ', known{i}]);
  end
  
  % prepare all syms into one vector
  symbolic = [known, unknown];
  
  reverseStr = '';
  % for each inputData
  for i = 1:size(inputData, 1)
    
    % prepare vector of knowns
    substitute = cell(size(symbolic));
    substitute(1, 1:length(known)) = num2cell(reshape(inputData{i}', 1, length(known)));
    
    % for each solution
    for j = 1:size(solution{i}, 2);
      
      % prepare vector of unknowns
      substitute(1, length(known)+1:end) = num2cell(solution{i}(:, j)');
      
      
      % enlarge matrix if needed
      if j > size(err, 2)
        oldErr = err;
        err = ones(size(solution, 1), j)*NaN;
        err(:, 1:size(oldErr, 2)) = oldErr;
      end
      
      err(i, j) = sum(double(subs(eq, symbolic, substitute)));
      
    end
    
    msg = sprintf('  %2.0f %%%% done', i/size(inputData, 1)*100);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg) - 1);
  end
  fprintf(reverseStr);
  
end