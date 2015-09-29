% Pavel Trutman, pavel.trutman@fel.cvut.cz, May 2015
%
% Remove polynomials that are not necessary for building of the action
% matrix with fixed solution for one variable.
%
% M - a matrix from which we are removing polynomials
% amStats - inforamtion about what is required for buiding of the action
%   matrix
% foundVar - variable with respect to will be the action matrix build
% prime - we are computing modulo this prime
% filter - a list of rows that have to remain in the matrix M

function [filter, foundVar] = gbs_RemoveUnnecessary_fixed(M, amStats, foundVar, prime)

  % remove not necesary polynomials
  fprintf('Removing not necessary polynomials:');
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
    
    filterOld = filter;
    filter = setdiff(filter, up:down);
    
    var = gbs_CheckActionMatrixConditions(M(filter, :), amStats, false, prime);
    
    if var == foundVar
      if step > 1
        fprintf([' ', int2str(up), '-', int2str(down)]);
      else
        fprintf([' ', int2str(down)]);
      end
      foundVar = var;
      up = down + 1;
      step = 2*step;
    else
      if step == 1
        up = up + 1;
      else
        step = max([floor(step/4) 1]);
      end
      filter = filterOld;
    end
    
  end
  fprintf('\n');
  
end
