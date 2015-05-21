% Generate Matlab Code for given action matrix and coefficient matrices
% (GBsolver subroutine)
% by Martin Bujnak, sep2008
% last edit by Pavel Trutman, May 2015


function [res] = gbs_ExportMapleCode(filename, M, trace, coefscode, known, knowngroups, unknown, algB, actMvar, amrows, amcols, gjcols, aidx, lastElim, cfg)

  [p, probname, e] = fileparts(filename);
  if isempty(e)
    filename = [filename '.txt'];
  end;

  if (~isdir(p))
    mkdir(p);
  end

  if isempty(knowngroups)
    knowngroups = 1:length(known);
  end

  % generate coefs calculation code
  knvars=[];
  knvarnames=[];
  knvarcnt=[];
  for i=1:length(known)

    if length(knvars) >= knowngroups(i)
      knvars{knowngroups(i)} = [knvars{knowngroups(i)} sym(known{i})];
      knvarcnt(knowngroups(i)) = knvarcnt(knowngroups(i)) + 1;
    else
      knvars{knowngroups(i)} = sym(known{i});
      knvarcnt(knowngroups(i)) = 1;
    end

    if length(knvarnames) < knowngroups(i) || isempty(knvarnames(knowngroups(i)))
      name=known(i);
      knvarnames{knowngroups(i)} = name{1};
    end
  end

  fid = fopen(filename, 'w');

  fprintf(fid, '# Generated using GBSolver generator Copyright Martin Bujnak,\n');
  fprintf(fid, '# Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.\n# \n');
  fprintf(fid, '# Please refer to the following paper, when using this code :\n');
  fprintf(fid, '#      Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,\n');
  fprintf(fid, '#      ECCV 2008, Marseille, France, October 12-18, 2008\n# \n');

  fprintf(fid, '> restart:\n');
  fprintf(fid, '> with(LinearAlgebra):\n');
  fprintf(fid, '> interface(rtablesize = 210):\n');
  fprintf(fid, '> Digits:=100:\n\n');

  fprintf(fid, '# #\n');
  fprintf(fid, '# #\n');
  fprintf(fid, '# # Solver \n');
  fprintf(fid, '# #\n');
  fprintf(fid, '>  \n');
  fprintf(fid, ['> ' probname ':=proc(' c2s(knvarnames, ', ') ')' '\n> \n']);

  % local variables
  fprintf(fid, ['> \tlocal c, M, Mold, amcols, A, D1, V1, i, ' c2s((unknown), ', ') ' , mat1, mat2']);
  if strcmp(cfg.PolynomialsGenerator, 'F4')
    for i = 1:length(trace) - 1;
      fprintf(fid, [', F', int2str(i)]);
    end
  end
  fprintf(fid, ';\n> \n');
  
  % coeffs
  fprintf(fid, '> \t# precalculate polynomial equations coefficients\n');
  for i=1:length(coefscode)

    % rename coefficients according to "knowngroups"
    coefcode = char(coefscode(i));

    for j=1:length(knvars)
      if length(knvars{j}) > 1
        for k=1:length(knvars{j})
          coefcode = regexprep([coefcode, ' '], [char(knvars{j}(k)), '([^0-9\(])'], [knvarnames{j} '(' int2str(k) ')$1']);
        end
      end
    end

    % replace (,) with [,]
    coefcode = strrep(coefcode, '(', '[');
    coefcode = strrep(coefcode, ')', ']');

    fprintf(fid, ['> \tc[' int2str(i) '] := ' regexprep(coefcode, '\s*$', '') ':\n']);
  end
  fprintf(fid, '> \n');


  % Generate code using different methods
  CodeGenerator = str2func(['CodeGen_Maple_', cfg.PolynomialsGenerator]);
  CodeGenerator(trace, lastElim, gjcols, @rrefPart, coefscode, fid);
  

  % action matrix
  fprintf(fid, ['> \tA := Matrix(' int2str(length(amrows)) ', ' int2str(length(amrows)) ', 0):\n']);
  fprintf(fid, ['> \tamcols := [' l2s(amcols, ', ') ']:\n']);

  tgcols = ['1..' int2str(length(amrows))];

  for i=1:length(amrows)

    if amrows(i) < 0
      fprintf(fid, ['> \tA[' int2str(i) ', ' int2str(-amrows(i)) '] := 1:\n']);
    else
      fprintf(fid, ['> \tA[' int2str(i) ', ' tgcols '] := -M[' int2str(amrows(i)) ', amcols]:\n']);
    end
  end
  fprintf(fid, '> \n');

  % solution extraction

  fprintf(fid, '> \t(D1, V1) := Eigenvectors(evalf(A)):\n');
  fprintf(fid, '>\n');

  [oneidx, unksidx] = gbs_GetVariablesIdx(algB, unknown);
  varsinvec = find(unksidx > 0);

  ucnt = length(unknown);
  for i=1:ucnt

    fprintf(fid, ['> \t' unknown{ucnt - i + 1} ' := Vector(' int2str(length(amrows)) ', 0): \n']);
  end

  if (sum(unksidx == 0)) > 0

    idx = find(unksidx == 0);
    if (length(idx) > 1)

      fprintf(fid, '\t\tWARNING: cannot extract all unknowns at once. A back-substitution required (not implemented/automatized)\n');
    end
  end

  fprintf(fid, ['> \tfor i from 1 to ' int2str(length(amrows)) ' do  \n']);

  ucnt = length(unknown);
  for i=1:ucnt

    if unksidx(i) == 0
      fprintf(fid, ['> \t\t' unknown{i} '[i] := evalf(D1[ i, i]) \n']);
    else
      fprintf(fid, ['> \t\t' unknown{i} '[i] := evalf(V1[' int2str(unksidx(i)) ', i]) / evalf(V1[' int2str(oneidx) ', i]): \n']);
    end
  end

  fprintf(fid, '> \tend do;  \n');

  % outputs
  fprintf(fid, '> \n');
  fprintf(fid, ['> \t(' c2s((unknown), ', ') ');\n']);
  fprintf(fid, '> \n');
  fprintf(fid, '> end proc:\n');

  fclose(fid);
  
  
  function [] = rrefPart(workflow, matrixName, last)
    %matrix elimination with partitioning
    fprintf(fid, '> \n> \t# GJ elimination with partitioning\n');
    
    %first part of matrix
    mat1Cols = [workflow.noAmCols(:, [workflow.ACols1; workflow.BCols]) workflow.amCols];
    fprintf(fid, ['> \tmat1 := ', matrixName, '[[', l2s(workflow.PRows1, ', '), '], [', l2s(mat1Cols, ', '), ']]:\n']);
    fprintf(fid, ['> \tmat1[1..-1, [', l2s(workflow.mat1NonzeroCols, ', '), ']] := ReducedRowEchelonForm(mat1[1..-1, [', l2s(workflow.mat1NonzeroCols, ', '), ']]):\n']);
    
    %second part of matrix
    mat2Cols = [workflow.noAmCols(:, [workflow.ACols2; workflow.BCols]) workflow.amCols];
    fprintf(fid, ['> \tmat2 := ', matrixName, '[[', l2s(workflow.PRows2, ', '), '], [', l2s(mat2Cols, ', '), ']]:\n']);
    fprintf(fid, ['> \tmat2[1..-1, [', l2s(workflow.mat2NonzeroCols, ', '), ']] := ReducedRowEchelonForm(mat2[1..-1, [', l2s(workflow.mat2NonzeroCols, ', '), ']]):\n> \n']);
    
    %assemble both parts together
    fprintf(fid, ['> \t', matrixName, ' := Matrix(', l2s(size(workflow.res), ', '), ', 0):\n']);
    if size(workflow.mat1TopRows, 2) ~= 0
      resMat1TopRows = workflow.resMat1TopRows;
      if ~last
        resMat1TopRows = workflow.permutationRows(resMat1TopRows);
      end
      fprintf(fid, ['> \t', matrixName, '[[', l2s(resMat1TopRows, ', '), '], [', l2s(mat1Cols, ', '), ']] := mat1[[', l2s(workflow.mat1TopRows, ', '), '], 1..-1]:\n']);
    end
    if size(workflow.mat2TopRows, 2) ~= 0
      resMat2TopRows = workflow.resMat2TopRows;
      if ~last
        resMat2TopRows = workflow.permutationRows(resMat2TopRows);
      end
      fprintf(fid, ['> \t', matrixName, '[[', l2s(resMat2TopRows, ', '), '], [', l2s(mat2Cols, ', '), ']] := mat2[[', l2s(workflow.mat2TopRows, ', '), '], 1..-1]:\n']);
    end
    if size(workflow.mat1BottRows, 2) ~= 0
      resMat1BottRows = workflow.resMat1BottRows;
      if ~last
        resMat1BottRows = workflow.permutationRows(resMat1BottRows);
      end
      fprintf(fid, ['> \t', matrixName, '[[', l2s(resMat1BottRows, ', '), '], [', l2s(mat1Cols, ', '), ']] := mat1[[', l2s(workflow.mat1BottRows, ' '), '], 1..-1]:\n']);
    end
    if size(workflow.mat2BottRows, 2) ~= 0
      resMat2BottRows = workflow.resMat2BottRows;
      if ~last
        resMat2BottRows = workflow.permutationRows(resMat2BottRows);
      end
      fprintf(fid, ['> \t', matrixName, '[[', l2s(resMat2BottRows, ', '), '], [', l2s(mat2Cols, ', '), ']] := mat2[[', l2s(workflow.mat2BottRows, ' '), '], 1..-1]:\n']);
    end
    
    %eliminate bottom rows of the matrix
    if size(workflow.bottomRows, 2) ~= 0
      bottomRows = workflow.bottomRows;
      if ~last
        bottomRows = workflow.permutationRows(bottomRows);
      end
      fprintf(fid, ['> \t', matrixName, '[[', l2s(bottomRows, ', '), '], [', l2s(workflow.resBottNonzeroCols, ', '), ']] := ReducedRowEchelonForm(', matrixName, '[[', l2s(bottomRows, ', '), '], [', l2s(workflow.resBottNonzeroCols, ', '), ']]):\n']);
    end
    
    fprintf(fid, '> \n');
    
    if last
      %eliminate amrows
      elimRows = amrows(amrows > 0);
      elimRows = setdiff(elimRows, workflow.bottomRows);
      for col = workflow.noAmCols(workflow.BCols)'
        [pivotRow, ~] = find(workflow.res(workflow.bottomRows, col) == 1);
        if size(pivotRow, 1) ~= 0
          pivotRow = workflow.bottomRows(pivotRow);
          for row = elimRows'
            if row > 0
              if workflow.res(row, col) ~= 0
                fprintf(fid, ['> \t', matrixName, '[', int2str(row), ', 1..-1] := ', matrixName, '[', int2str(row), ', 1..-1] - ', matrixName, '[', int2str(row), ', ', int2str(col), ']*', matrixName, '[', int2str(pivotRow), ', 1..-1]:\n']);
              end
            end
          end
        end
      end
      fprintf(fid, '> \n');
    else
      %eliminate all
      if isfield(workflow, 'elim')
        for elim = workflow.elim
          elim = elim{1};
          if strcmp(elim.type, 'divide')
            fprintf(fid, ['> \t', matrixName, '[', int2str(elim.row), ', 1..-1] := ', matrixName, '[', int2str(elim.row), ', 1..-1]/', matrixName, '[', int2str(elim.row), ', ', int2str(elim.col), ']:\n']);
          elseif strcmp(elim.type, 'switch')
            fprintf(fid, ['> \t', matrixName, '[[', l2s(elim.rows, ', '), '], 1..-1] := ', matrixName, '[[', l2s(elim.rows(end:-1:1), ', '), '], 1..-1]:\n']);
          elseif strcmp(elim.type, 'eliminate')
            fprintf(fid, ['> \t', matrixName, '[', int2str(elim.row), ', 1..-1] := ', matrixName, '[', int2str(elim.row), ', 1..-1] - ', matrixName, '[', int2str(elim.row), ', ', int2str(elim.col), ']*', matrixName, '[', int2str(elim.pivotRow), ', 1..-1]:\n']);
          end
        end
      end
    end
  end

end
