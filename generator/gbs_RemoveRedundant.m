% Pavel Trutman, pavel.trutman@fel.cvut.cz, May 2015
%
% Remove polynomials that are redundant in the matrix.
%
% M - a matrix from which we are removing polynomials
% prime - we are computing modulo this prime
% filter - a list of rows that have to remain in the matrix M


function [filter] = gbs_RemoveRedundant(M, prime)
  
  nonzero = find(sum(M) ~= 0);
  Kk = M(:, nonzero);
  B = gjzpsp(Kk, prime);
  
  toRemove = size(find(sum(B, 2) == 0), 1);
  filter = 1:size(M, 1);
  
  reverseStr = '';
  if toRemove > 0
    fprintf('    removing %d equations', toRemove);
    usedRows = size(find(sum(B, 2) ~= 0), 1);
    removed = 0;
    step = max([floor(toRemove/4) 1]);
    up = 1;
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
          fprintf(reverseStr);
          fprintf(' - all removed');
          break;
        end
        up = down + 1;
        step = min([2*step toRemove-removed]);
      end
      
      % print progress
      msg = sprintf(' (%2.0f %%%% removed, remaining %d equations to check)', floor(removed/toRemove*100), size(M, 1) - up);
      fprintf([reverseStr msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg) - 1);
      
    end
    fprintf('\n');
 
  end
  
end

