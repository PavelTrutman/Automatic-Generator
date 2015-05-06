% Pavel Trutman, pavel.trutman@fel.cvut.cz, May 2015
%
% Remove polynomials that are not necessary for building of the action
% matrix.
%
% M - a matrix from which we are removing polynomials
% amStats - inforamtion about what is required for buiding of the action
%   matrix
% foundVar - variable with respect to will be the action matrix build
% prime - we are computing modulo this prime
% filter - a list of rows that have to remain in the matrix M

function [filter] = gbs_RemoveUnnecessary(M, amStats, foundVar, prime)

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
        step = max([floor(step/4) 1]);
      end
      filter = filterOld;
    end
    
  end
  
end
