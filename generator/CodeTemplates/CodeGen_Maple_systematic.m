% Template to generate Maple code of solvers using systematic method of
% generating polynomials.

% Pavel Trutman, pavel.trutman@fel.cvut.cz, May 2015

function [ ] = CodeGen_Maple_systematic(trace, lastElim, gjcols, rrefPart, coefscode, fid)
  
  % coefs matrix
  trace{end}.Mcoefs = trace{end}.Mcoefs(:, gjcols);

  if length(trace) == 1
    matrixName = 'M';
  else
    matrixName = 'M1';
    trace{1}.Mcoefs = trace{1}.Mcoefs(:, trace{1}.nonzerocols);
  end
  
  fprintf(fid, ['> \t', matrixName, ' := Matrix(' int2str(size(trace{1}.Mcoefs, 1)) ', ' int2str(size(trace{1}.Mcoefs, 2)) ', 0):\n']);
  for i=1:length(coefscode)

    [ofss] = find(trace{1}.Mcoefs == i)';
    for ofs = ofss
      fprintf(fid, ['> \t', matrixName, '(' int2str(ofs) ') := c[' int2str(i) ']:\n']);
    end
  end
  fprintf(fid, '>  \n');

  % elimination part
  if trace{1}.partitioning.enable || lastElim.enable
    %eliminate with partitioning
    if length(trace) == 1
      %last elimination
      rrefPart(lastElim, matrixName, 1);
    elseif trace{1}.partitioning.enable
      rrefPart(trace{1}.partitioning, matrixName, 0);
    else
      fprintf(fid, ['> \t', matrixName, ' := ReducedRowEchelonForm(', matrixName, '):\n']);
    end
  else
    fprintf(fid, ['> \t', matrixName, ' := ReducedRowEchelonForm(', matrixName, '):\n']);
  end
  fprintf(fid, '> \n');

  for i = 2:length(trace)
    if i == length(trace)
      matrixName = 'M';
    else
      matrixName = ['M', int2str(i)];
      trace{i}.Mcoefs = trace{i}.Mcoefs(:, trace{i}.nonzerocols);
    end
    %if trace{i - 1}.partitioning.enable
    %  fprintf(fid, ['> \tMold := M:\n']);
    %else
    %  fprintf(fid, ['> \tMold := Matrix(' int2str(size(trace{1}.Mcoefs, 1)) ', ' int2str(size(trace{1}.Mcoefs, 2)) ', 0):\n']);
    %  fprintf(fid, ['> \tMold[1..-1, [' l2s(trace{i - 1}.nonzerocols, ', ') ']] := M:\n']);
    %end
    fprintf(fid, ['> \t', matrixName, ' := Matrix(' int2str(size(trace{i}.Mcoefs, 1)) ', ' int2str(size(trace{i}.Mcoefs, 2)) ', 0):\n']);
    if isfield(trace{i}, 'filter')
      oldColumns = gjcols - trace{i}.columnfrom + 1;
      oldColumns = intersect(oldColumns(oldColumns > 0), trace{i-1}.nonzerocols);
      [~, selectCols] = ismember(oldColumns, trace{i-1}.nonzerocols);
      [~, newCols] = ismember(oldColumns + trace{i}.size(2) - trace{i-1}.size(2), gjcols);
      fprintf(fid, ['> \t', matrixName, '[' int2str(trace{i}.rowfrom) '..' int2str(trace{i}.rowto) ', [', l2s(newCols, ', '), ']] := M', int2str(i-1), '[[' l2s(trace{i}.filter, ', ') '], [' l2s(selectCols, ', ') ']]:\n']);
    else
      [~, newCols] = ismember(trace{i-1}.nonzerocols + trace{i}.size(2) - trace{i-1}.size(2), trace{i}.nonzerocols);
      fprintf(fid, ['> \t', matrixName, '[' int2str(trace{i}.rowfrom) '..' int2str(trace{i}.rowto) ', [', l2s(newCols, ', '), ']] := M', int2str(i-1), '[1..', int2str(trace{i}.rowsold), ', 1..-1]:\n']);
    end

    [ofs] = find(trace{i}.Mcoefs);
    for j = ofs'
      [k, l] = ind2sub(trace{i-1}.size, trace{i}.Mcoefs(j));
      l = find(trace{i-1}.nonzerocols == l, 1, 'first');
      idx = sub2ind(trace{i-1}.size, k, l);
      fprintf(fid, ['> \t', matrixName, '(' int2str(j) ') := M', int2str(i-1), '(' int2str(idx) '):\n']);
    end

    fprintf(fid, '> \n');
        % elimination part
    if trace{i}.partitioning.enable || lastElim.enable
      %eliminate with partitioning
      if i == length(trace)
        %last elimination
        rrefPart(lastElim, matrixName, 1);
      elseif trace{i}.partitioning.enable
        rrefPart(trace{i}.partitioning, matrixName, 0);
      else
        fprintf(fid, ['\n> \t', matrixName, ' := ReducedRowEchelonForm(', matrixName, '):\n']);
      end
    else
      fprintf(fid, ['\n> \t', matrixName, ' := ReducedRowEchelonForm(', matrixName, '):\n']);
    end
    fprintf(fid, '> \n');
  end
  
end
