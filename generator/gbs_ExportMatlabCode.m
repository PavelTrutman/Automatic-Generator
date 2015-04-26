% Generate Matlab Code for given action matrix and coefficient matrices
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, April 2015


function [res] = gbs_ExportMatlabCode(filename, M, trace, coefscode, known, knowngroups, unknown, algB, actMvar, amrows, amcols, gjcols, aidx, lastElim, cfg)

  [p, probname, e] = fileparts(filename);
  if isempty(e)
    filename = [filename '.m'];
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

  fprintf(fid, '%% Generated using GBSolver generator Copyright Martin Bujnak,\n');
  fprintf(fid, '%% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.\n%% \n');
  fprintf(fid, '%% Please refer to the following paper, when using this code :\n');
  fprintf(fid, '%%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,\n');
  fprintf(fid, '%%     ECCV 2008, Marseille, France, October 12-18, 2008\n');
  fprintf(fid, '\n');
  if cfg.benchmark
    fprintf(fid, ['function [unknowns, benchData] = ' probname '(args)\n\n']);
  else
    fprintf(fid, ['function [' c2s((unknown), ', ') '] = ' probname '(' c2s(knvarnames, ', ') ')\n\n']);
  end

  if cfg.benchmark
    fprintf(fid, '\t%% This is a benchmark solver!\n\n');
    for i = 1:length(knvarnames)
      fprintf(fid, ['\t', c2s(knvarnames(i)), ' = args(', int2str(i), ', :);\n']);
    end
    fprintf(fid, '\n');
  end
  
  % coeffs
  fprintf(fid, '\t%% precalculate polynomial equations coefficients\n');
  for i=1:length(coefscode)

    % rename coefficients according to "knowngroups"
    coefcode = char(coefscode(i));

    for j=1:length(knvars)
      if length(knvars{j}) > 1
        for k=1:length(knvars{j})
          coefcode = regexprep([coefcode, ' '], ['(', char(knvars{j}(k)), ')[^(0-9]'], [knvarnames{j} '(' int2str(k) ')']);
        end
      end
    end

    fprintf(fid, ['\tc(' int2str(i) ') = ' regexprep(coefcode, '\s*$', '') ';\n']);
  end
  fprintf(fid, '\n');
  
  
  % Generate code using different methods
  CodeGenerator = str2func(['CodeGen_Matlab_', cfg.PolynomialsGenerator]);
  CodeGenerator(trace, lastElim, gjcols, @rrefPart, coefscode, fid);

  
  % action matrix
  fprintf(fid, ['\tA = zeros(' int2str(length(amrows)) ');\n']);
  fprintf(fid, ['\tamcols = [' l2s(amcols, ' ') '];\n']);
  for i=1:length(amrows)

    if amrows(i) < 0
      fprintf(fid, ['\tA(' int2str(i) ', ' int2str(-amrows(i)) ') = 1;\n']);
    else
      fprintf(fid, ['\tA(' int2str(i) ', :) = -M(' int2str(amrows(i)) ', amcols);\n']);
    end
  end
  fprintf(fid, '\n');

  % solution extraction

  fprintf(fid, '\t[V D] = eig(A);\n');

  [oneidx, unksidx] = gbs_GetVariablesIdx(algB, unknown);
  varsinvec = find(unksidx > 0);

  fprintf(fid, ['\tsol =  V([' l2s(unksidx(varsinvec), ', ') '],:)./(ones(' int2str(length(varsinvec)) ', ' int2str(oneidx) ')*V(' int2str(oneidx) ',:));\n']);
  fprintf(fid, '\n');

  fprintf(fid, '\tif (find(isnan(sol(:))) > 0)\n\t\t\n');

  % division by zero filter
  ucnt = length(unknown);
  for i=1:ucnt

    fprintf(fid, ['\t\t' unknown{i} ' = zeros(1, 0);\n']);
  end

  fprintf(fid, '\telse\n\t\t\n');

  % extract variables
  if (sum(unksidx == 0)) > 0

    idx = find(unksidx == 0);
    if (length(idx) > 1)

      fprintf(fid, '\t\tWARNING: cannot extract all unknowns at once. A back-substitution required (not implemented/automatized)\n');
    end

    fprintf(fid, ['\t\tev  = diag(D);\n']);
    fprintf(fid, ['\t\tI = find(not(imag( sol(1,:) )) & not(imag( ev )));\n']);
  else
    fprintf(fid, ['\t\tI = find(not(imag( sol(1,:) )));\n']);
  end

  ui = 1;
  for i=1:ucnt

    if unksidx(i) == 0

      % action variable ?
      if strcmp(actMvar, unknown{i})

        fprintf(fid, ['\t\t' unknown{i} ' = ev(I);\n']);
      else

        fprintf(fid, '\t\tWARNING: one or more unknowns could not be extracted.\n');
      end

    else
      fprintf(fid, ['\t\t' unknown{i} ' = sol(' int2str(ui) ',I);\n']);
      ui = ui+1;
    end
  end

  fprintf(fid, '\tend\n');
  
  if cfg.benchmark
    fprintf(fid, '\n');
    fprintf(fid, '\t%%save benchmark data\n');
    for i = 1:length(unknown)
      fprintf(fid, ['\tunknowns(', int2str(i), ', :) = ', c2s(unknown(i)), ';\n']);
    end
    fprintf(fid, '\n');
    fprintf(fid, ['\tbenchData.polynomials = zeros([', l2s(size(M), ' '), ']);\n']);
    fprintf(fid, ['\tbenchData.polynomials(:, [', l2s(gjcols, ' '), ']) = M;\n']);
    fprintf(fid, '\n');
  end
  
  fprintf(fid, 'end\n');

  fclose(fid);
  
  
  function [] = rrefPart(workflow, last)
    %matrix elimination with partitioning
    fprintf(fid, '\n\t%%GJ elimination with partitioning\n');
    
    %first part of matrix
    mat1Cols = [workflow.noAmCols(:, [workflow.ACols1; workflow.BCols]) workflow.amCols];
    fprintf(fid, ['\tmat1 = M([', l2s(workflow.PRows1, ' '), '], [', l2s(mat1Cols, ' '), ']);\n']);
    fprintf(fid, ['\tmat1(:, [', l2s(workflow.mat1NonzeroCols, ' '), ']) = rref(mat1(:, [', l2s(workflow.mat1NonzeroCols, ' '), ']));\n']);
    
    %second part of matrix
    mat2Cols = [workflow.noAmCols(:, [workflow.ACols2; workflow.BCols]) workflow.amCols];
    fprintf(fid, ['\tmat2 = M([', l2s(workflow.PRows2, ' '), '], [', l2s(mat2Cols, ' '), ']);\n']);
    fprintf(fid, ['\tmat2(:, [', l2s(workflow.mat2NonzeroCols, ' '), ']) = rref(mat2(:, [', l2s(workflow.mat2NonzeroCols, ' '), ']));\n\n']);
    
    %assemble both parts together
    fprintf(fid, ['\tM = zeros([', l2s(size(workflow.res), ' '), ']);\n']);
    if size(workflow.mat1TopRows, 2) ~= 0
      resMat1TopRows = workflow.resMat1TopRows;
      if ~last
        resMat1TopRows = workflow.permutationRows(resMat1TopRows);
      end
      fprintf(fid, ['\tM([', l2s(resMat1TopRows, ' '), '], [', l2s(mat1Cols, ' '), ']) = mat1([', l2s(workflow.mat1TopRows, ' '), '], :);\n']);
    end
    if size(workflow.mat2TopRows, 2) ~= 0
      resMat2TopRows = workflow.resMat2TopRows;
      if ~last
        resMat2TopRows = workflow.permutationRows(resMat2TopRows);
      end
      fprintf(fid, ['\tM([', l2s(resMat2TopRows, ' '), '], [', l2s(mat2Cols, ' '), ']) = mat2([', l2s(workflow.mat2TopRows, ' '), '], :);\n']);
    end
    if size(workflow.mat1BottRows, 2) ~= 0
      resMat1BottRows = workflow.resMat1BottRows;
      if ~last
        resMat1BottRows = workflow.permutationRows(resMat1BottRows);
      end
      fprintf(fid, ['\tM([', l2s(resMat1BottRows, ' '), '], [', l2s(mat1Cols, ' '), ']) = mat1([', l2s(workflow.mat1BottRows, ' '), '], :);\n']);
    end
    if size(workflow.mat2BottRows, 2) ~= 0
      resMat2BottRows = workflow.resMat2BottRows;
      if ~last
        resMat2BottRows = workflow.permutationRows(resMat2BottRows);
      end
      fprintf(fid, ['\tM([', l2s(resMat2BottRows, ' '), '], [', l2s(mat2Cols, ' '), ']) = mat2([', l2s(workflow.mat2BottRows, ' '), '], :);\n']);
    end
    
    %eliminate bottom rows of the matrix
    if size(workflow.bottomRows, 2) ~= 0
      bottomRows = workflow.bottomRows;
      if ~last
        bottomRows = workflow.permutationRows(bottomRows);
      end
      fprintf(fid, ['\tM([', l2s(bottomRows, ' '), '], [', l2s(workflow.resBottNonzeroCols, ' '), ']) = rref(M([', l2s(bottomRows, ' '), '], [', l2s(workflow.resBottNonzeroCols, ' '), ']));\n']);
    end
    
    fprintf(fid, '\n');
    
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
                fprintf(fid, ['\tM(', int2str(row), ', :) = M(', int2str(row), ', :) - M(', int2str(row), ', ', int2str(col), ')*M(', int2str(pivotRow), ', :);\n']);
              end
            end
          end
        end
      end
      fprintf(fid, '\n');
    else
      %eliminate all
      if isfield(workflow, 'elim')
        for elim = workflow.elim
          elim = elim{1};
          if strcmp(elim.type, 'divide')
            fprintf(fid, ['\tM(', int2str(elim.row), ', :) = M(', int2str(elim.row), ', :)/M(', int2str(elim.row), ', ', int2str(elim.col), ');\n']);
          elseif strcmp(elim.type, 'switch')
            fprintf(fid, ['\tM([', l2s(elim.rows, ' '), '], :) = M([', l2s(elim.rows(end:1), ' '), '], :);\n']);
          elseif strcmp(elim.type, 'eliminate')
            fprintf(fid, ['\tM(', int2str(elim.row), ', :) = M(', int2str(elim.row), ', :) - M(', int2str(elim.row), ', ', int2str(elim.col), ')*M(', int2str(elim.pivotRow), ', :);\n']);
          end
        end
      end
    end
  end
  
end
