% Template to generate MATLAB code of solvers using F4 method of 
% generating polynomials.

% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015

function [ ] = CodeGen_Matlab_F4(trace, lastElim, gjcols, rrefPart, coefscode, fid)
  
  % elimination part
  for i = 1:length(trace) - 1
    fprintf(fid, ['\tF', int2str(i), ' = zeros(', int2str(size(trace{i}.coefs, 1)), ', ', int2str(size(trace{i}.nonzero, 2)), ');\n']);
    for j = 1:size(trace{i}.refs, 1)
      
      [ofs] = find(trace{i}.coefs(j, trace{i}.nonzero));
      for k = ofs
        idx = sub2ind([size(trace{i}.refs, 1), length(trace{i}.nonzero)], j, k);
        if trace{i}.refs(j, 1) ~= 0
          % copy polynomial from other matrix F
          idxFrom = sub2ind([size(trace{trace{i}.refs(j, 1)}.refs, 1), length(trace{trace{i}.refs(j, 1)}.nonzero)], trace{i}.refs(j, 2), find(trace{trace{i}.refs(j, 1)}.nonzero == (size(trace{trace{i}.refs(j, 1)}.coefs, 2) - trace{i}.coefs(j, trace{i}.nonzero(k)) + 1), 1, 'first'));
          fprintf(fid, ['\tF', int2str(i), '(', int2str(idx), ') = F', int2str(trace{i}.refs(j, 1)), '(', int2str(idxFrom), ');\n']);
        else
          % copy polynomial from input
          fprintf(fid, ['\tF', int2str(i), '(', int2str(idx), ') = c(', int2str(trace{i}.coefs(j, trace{i}.nonzero(k))), ');\n']);
        end
      end
    end
    fprintf(fid, ['\tF', int2str(i), ' = rref(F', int2str(i), ');\n\n']);
  end
  
  %last elimination
  fprintf(fid, ['\tM = zeros(', int2str(size(trace{end}.refs, 1)), ', ', int2str(size(gjcols, 2)), ');\n']);
  for i = 1:size(trace{end}.refs, 1)
    for j = 1:size(gjcols, 2)
      if trace{end}.coefs(i, gjcols(j)) ~= 0
        idx = sub2ind([size(trace{end}.refs, 1), size(gjcols, 2)], i, j);
        if trace{end}.refs(i, 1) ~= 0
          % copy polynomial from other matrix F
          idxFrom = sub2ind([size(trace{trace{end}.refs(i, 1)}.refs, 1), length(trace{trace{end}.refs(i, 1)}.nonzero)], trace{end}.refs(i, 2), find(trace{trace{end}.refs(i, 1)}.nonzero == (size(trace{trace{end}.refs(i, 1)}.coefs, 2) - trace{end}.coefs(i, gjcols(j)) + 1), 1, 'first'));
          fprintf(fid, ['\tM(', int2str(idx), ') = F', int2str(trace{end}.refs(i, 1)), '(', int2str(idxFrom), ');\n']);
        else
          % copy polynomial from input
          fprintf(fid, ['\tM(', int2str(idx), ') = c(', int2str(trace{end}.coefs(i, gjcols(j))), ');\n']);
        end
      end
    end
  end
  
  if lastElim.enable
    % last elimination with partitioning
    rrefPart(lastElim, 'M', 1);
  else
    % without partitioning
    fprintf(fid, '\tM = rref(M);\n\n');
  end
  
end
