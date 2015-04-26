% Template to generate Maple code of solvers using primitive method of
% generating polynomials.

% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015

function [ ] = CodeGen_Maple_Primitive(trace, lastElim, gjcols, rrefPart, coefscode, fid)
  
  % coefs matrix
  trace{end}.Mcoefs = trace{end}.Mcoefs(:, gjcols);

  fprintf(fid, ['> \tM := Matrix(' int2str(size(trace{1}.Mcoefs, 1)) ', ' int2str(size(trace{1}.Mcoefs, 2)) ', 0):\n']);
  for i=1:length(coefscode)

    [ofss] = find(trace{1}.Mcoefs == i)';
    for ofs = ofss
      fprintf(fid, ['> \tM(' int2str(ofs) ') := c[' int2str(i) ']:\n']);
    end
  end
  fprintf(fid, '>  \n');

  % elimination part
  if length(trace) == 1
    %last elimination
    if lastElim.enable
      %eliminate with partitioning
      rrefPart(lastElim, 1);
    else
      fprintf(fid, '> \tM := ReducedRowEchelonForm(M):\n');
    end
  else
    if trace{1}.partitioning.enable
      %eliminate with partitionig
      rrefPart(trace{1}.partitioning, 0);
    else
      fprintf(fid, ['> \tM := ReducedRowEchelonForm(M[1..-1, [' l2s(trace{1}.nonzerocols, ', ') ']]):\n']);
    end
  end
  fprintf(fid, '> \n');

  for i = 2:length(trace)
    if trace{i - 1}.partitioning.enable
      fprintf(fid, ['> \tMold := M:\n']);
    else
      fprintf(fid, ['> \tMold := Matrix(' int2str(size(trace{1}.Mcoefs, 1)) ', ' int2str(size(trace{1}.Mcoefs, 2)) ', 0):\n']);
      fprintf(fid, ['> \tMold[1..-1, [' l2s(trace{i - 1}.nonzerocols, ', ') ']] := M:\n']);
    end
    fprintf(fid, ['> \tM := Matrix(' int2str(size(trace{i}.Mcoefs, 1)) ', ' int2str(size(trace{i}.Mcoefs, 2)) ', 0):\n']);
    if isfield(trace{i}, 'filter')
      oldColumns = gjcols - trace{i}.columnfrom + 1;
      fprintf(fid, ['> \tM[' int2str(trace{i}.rowfrom) '..' int2str(trace{i}.rowto) ', ' int2str(size(gjcols, 2) - sum(gjcols >= trace{i}.columnfrom) + 1)  '..' int2str(size(gjcols, 2)) '] := Mold[[' l2s(trace{i}.filter, ', ') '], [' l2s(oldColumns(oldColumns > 0), ', ') ']]:\n']);
    else
      fprintf(fid, ['> \tM[' int2str(trace{i}.rowfrom) '..' int2str(trace{i}.rowto) ', ' int2str(trace{i}.columnfrom) '..' int2str(trace{i}.columnto) '] := Mold[1..', int2str(trace{i}.rowsold), ', 1..-1]:\n']);
    end

    [ofs] = find(trace{i}.Mcoefs);
    for j = ofs'
      fprintf(fid, ['> \tM(' int2str(j) ') := Mold(' int2str(trace{i}.Mcoefs(j)) '):\n']);
    end

    fprintf(fid, '> \n');
    if i == length(trace)
      %last elimination
      if lastElim.enable
        %use partitioning
        rrefPart(lastElim, 1);
      else
        fprintf(fid, ['\n> \tM := ReducedRowEchelonForm(M):\n']);
      end
    else
      if trace{i}.partitioning.enable
        %use partitioning
        rrefPart(trace{i}.partitioning, 0);
      else
        fprintf(fid, ['\n> \tM := ReducedRowEchelonForm(M[1..-1, [' l2s(trace{i}.nonzerocols, ', ') ']]):\n']);
      end
    end
  end

  fprintf(fid, '> \n');
  
end
