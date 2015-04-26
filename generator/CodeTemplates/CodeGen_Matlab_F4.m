% Template to generate MATLAB code of solvers using F4 method of 
% generating polynomials.

% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015

function [ ] = CodeGen_Matlab_F4(trace, lastElim, gjcols, rrefPart, coefscode, fid)
  
  % elimination part
  for i = 1:length(trace) - 1
    fprintf(fid, ['\tF', int2str(i), ' = zeros(', int2str(size(trace{i}.coefs, 1)), ', ', int2str(size(trace{i}.coefs, 2)), ');\n']);
    for j = 1:size(trace{i}.refs, 1)
      
      [ofs] = find(trace{i}.coefs(j, :));
      for k = ofs
        if trace{i}.refs(j, 1) ~= 0
          % copy polynomial from other matrix F
          fprintf(fid, ['\tF', int2str(i), '(', int2str(j), ', ',  int2str(k), ') = F', int2str(trace{i}.refs(j, 1)), '(', int2str(trace{i}.refs(j, 2)), ', ', int2str(size(trace{trace{i}.refs(j, 1)}.coefs, 2) - trace{i}.coefs(j, k) + 1), ');\n']);
        else
          % copy polynomial from input
          fprintf(fid, ['\tF', int2str(i), '(', int2str(j), ', ', int2str(k), ') = c(', int2str(trace{i}.coefs(j, k)), ');\n']);
        end
      end
    end
    fprintf(fid, ['\tF', int2str(i), ' = rref(F', int2str(i), ');\n\n']);
  end
  
  %last elimination
  fprintf(fid, ['\tM = zeros(', int2str(size(trace{end}.refs, 1)), ', ', int2str(size(gjcols, 2)), ');\n']);
  for i = 1:size(trace{end}.refs, 1)
    if trace{end}.refs(i, 1) ~= 0
      % copy polynomial from other matrix F
      for j = 1:size(gjcols, 2)
        if trace{end}.coefs(i, gjcols(j)) ~= 0
          fprintf(fid, ['\tM(', int2str(i), ', ', int2str(j), ') = F', int2str(trace{end}.refs(i, 1)), '(', int2str(trace{end}.refs(i, 2)), ', ', int2str(size(trace{trace{end}.refs(i, 1)}.coefs, 2) - trace{end}.coefs(i, gjcols(j)) + 1), ');\n']);
        end
      end
    else
      % copy polynomial from input
      for j = 1:size(gjcols, 2)
        if trace{end}.coefs(i, gjcols(j)) ~= 0
          fprintf(fid, ['\tM(', int2str(i), ', ', int2str(j), ') = c(', int2str(trace{end}.coefs(i, gjcols(j))), ');\n']);
        end
      end
    end
  end
  fprintf(fid, '\tM = rref(M);\n\n');
  
end
