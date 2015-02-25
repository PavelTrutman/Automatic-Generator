% Generate Matlab Code for given action matrix and coefficient matrices
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, February 2015


function [res] = gbs_ExportMCode(filename, M, trace, coefscode, known, knowngroups, unknown, algB, actMvar, amrows, amcols, gjcols, aidx, PaToH)

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
            
            %name=regexp(known{i},'\D*', 'match');
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
    fprintf(fid, ['function [' c2s((unknown), ', ') '] = ' probname '(' c2s(knvarnames, ', ') ')\n\n']);
    
    % coeffs
    fprintf(fid, '\t%% precalculate polynomial equations coefficients\n');
    for i=1:length(coefscode)
        
        % rename coefficients according to "knowngroups"
        coefcode = char(coefscode(i));
        
        for j=1:length(knvars)
            if length(knvars{j}) > 1
                for k=1:length(knvars{j})
                    coefcode = strrep(coefcode, char(knvars{j}(k)), [knvarnames{j} '(' int2str(k) ')']);
                end
            end
        end
        
        fprintf(fid, ['\tc(' int2str(i) ') = ' coefcode ';\n']);
    end
    fprintf(fid, '\n');
    
    % coefs matrix
    trace{end}.Mcoefs = trace{end}.Mcoefs(:, gjcols);
    
    fprintf(fid, ['\tM = zeros(' int2str(size(trace{1}.Mcoefs, 1)) ', ' int2str(size(trace{1}.Mcoefs, 2)) ');\n']);
    for i=1:length(coefscode)
    
        [ofs] = find(trace{1}.Mcoefs == i);

        if length(ofs) == 1
            
            fprintf(fid, ['\tM(' int2str(ofs) ') = c(' int2str(i) ');\n']);
        else
            fprintf(fid, ['\tci = [']);
            fprintf(fid, l2s(ofs, ', '));
            fprintf(fid, '];\n');

            fprintf(fid, ['\tM(ci) = c(' int2str(i) ');\n\n']);
        end
    end
    fprintf(fid, '\n');
   
    % elimination part
    if length(trace) == 1
      %last elimination
      if PaToH.enable
        %use PaToH
        fprintf(fid, '\tM = PaToH_rref(M);\n');
      else
        fprintf(fid, '\tM = rref(M);\n');
      end
    else
      fprintf(fid, ['\tM = rref(M(:, [' l2s(trace{1}.nonzerocols, ' ') ']));\n']);
    end
    fprintf(fid, '\n');
    
    for i = 2:length(trace)
      fprintf(fid, ['\tMold = zeros(' int2str(size(trace{i - 1}.Mcoefs, 1)) ', ' int2str(size(trace{i - 1}.Mcoefs, 2)) ');\n']);
      fprintf(fid, ['\tMold(:, [' l2s(trace{i - 1}.nonzerocols, ' ') ']) = M;\n']);
      fprintf(fid, ['\tM = zeros(' int2str(size(trace{i}.Mcoefs, 1)) ', ' int2str(size(trace{i}.Mcoefs, 2)) ');\n']);
      if isfield(trace{i}, 'filter')
        oldColumns = gjcols - trace{i}.columnfrom + 1;
        fprintf(fid, ['\tM(' int2str(trace{i}.rowfrom) ':' int2str(trace{i}.rowto) ', ' int2str(size(gjcols, 2) - sum(gjcols >= trace{i}.columnfrom) + 1)  ':' int2str(size(gjcols, 2)) ') = Mold([' l2s(trace{i}.filter, ' ') '], [', l2s(oldColumns(oldColumns > 0), ' '), ']);\n']);
      else
        fprintf(fid, ['\tM(' int2str(trace{i}.rowfrom) ':' int2str(trace{i}.rowto) ', ' int2str(trace{i}.columnfrom)  ':' int2str(trace{i}.columnto) ') = Mold(1:', int2str(trace{i}.rowsold), ', :);\n']);
      end
      
      [ofs] = find(trace{i}.Mcoefs);
      for j = ofs'
        fprintf(fid, ['\tM(' int2str(j) ') = Mold(' int2str(trace{i}.Mcoefs(j)) ');\n']);
      end
      if i == length(trace)
        %last elimination
        if PaToH.enable
          %use PaToH
          fprintf(fid, '\n\tM = PaToH_rref(M);\n\n');
        else
          fprintf(fid, '\n\tM = rref(M);\n\n');
        end
      else
        fprintf(fid, ['\n\tM = rref(M(:, [' l2s(trace{i}.nonzerocols, ' ') ']));\n\n']);
      end
    end
    
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
        
        fprintf(fid, ['\t\t' unknown{i} ' = [];\n']);
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
    fprintf(fid, 'end\n');
    
    
    %PaToH elimination function
    if PaToH.enable
      fprintf(fid, '\nfunction [res] = PaToH_rref(matrix)\n\n');
      
      %first part of matrix
      mat1Cols = [PaToH.noAmCols(:, [PaToH.ACols1; PaToH.BCols]) PaToH.amCols];
      fprintf(fid, ['\tmat1 = matrix([', l2s(PaToH.PRows1, ' '), '], [', l2s(mat1Cols, ' '), ']);\n']);
      fprintf(fid, '\tmat1 = rref(mat1);\n');
      
      %second part of matrix
      mat2Cols = [PaToH.noAmCols(:, [PaToH.ACols2; PaToH.BCols]) PaToH.amCols];
      fprintf(fid, ['\tmat2 = matrix([', l2s(PaToH.PRows2, ' '), '], [', l2s(mat2Cols, ' '), ']);\n']);
      fprintf(fid, '\tmat2 = rref(mat2);\n\n');
      
      %assemble both parts together
      fprintf(fid, ['\tres = zeros([', l2s(size(trace{end}.Mcoefs), ' '), ']);\n']);
      if size(PaToH.mat1TopRows, 2) ~= 0
        fprintf(fid, ['\tres([', l2s(PaToH.resMat1TopRows, ' '), '], [', l2s(mat1Cols, ' '), ']) = mat1([', l2s(PaToH.mat1TopRows, ' '), '], :);\n']);
      end
      if size(PaToH.mat2TopRows, 2) ~= 0
        fprintf(fid, ['\tres([', l2s(PaToH.resMat2TopRows, ' '), '], [', l2s(mat2Cols, ' '), ']) = mat2([', l2s(PaToH.mat2TopRows, ' '), '], :);\n']);
      end
      if size(PaToH.mat1BottRows, 2) ~= 0
        fprintf(fid, ['\tres([', l2s(PaToH.resMat1BottRows, ' '), '], [', l2s(mat1Cols, ' '), ']) = mat1([', l2s(PaToH.mat1BottRows, ' '), '], :);\n']);
      end
      if size(PaToH.mat2BottRows, 2) ~= 0
        fprintf(fid, ['\tres([', l2s(PaToH.resMat2BottRows, ' '), '], [', l2s(mat2Cols, ' '), ']) = mat2([', l2s(PaToH.mat2BottRows, ' '), '], :);\n']);
      end
      
      %eliminate bottom rows of the matrix
      if size(PaToH.bottomRows, 2) ~= 0
        fprintf(fid, ['\tres([', l2s(PaToH.bottomRows, ' '), '], [', l2s(PaToH.resBottNonzeroCols, ' '), ']) = rref(res([', l2s(PaToH.bottomRows, ' '), '], [', l2s(PaToH.resBottNonzeroCols, ' '), ']));\n']);
      end
      
      fprintf(fid, '\n');
      
      %eliminate amrows
      elimRows = amrows(amrows > 0);
      elimRows = setdiff(elimRows, PaToH.bottomRows);
      for col = PaToH.BCols'
        [pivotRow, ~] = find(PaToH.res(PaToH.bottomRows, col) == 1);
        if size(pivotRow, 1) ~= 0
          pivotRow = PaToH.bottomRows(pivotRow);
          for row = elimRows'
            if row > 0    
              if PaToH.res(row, col) ~= 0
                fprintf(fid, ['\tres(', int2str(row), ', :) = res(', int2str(row), ', :) - res(', int2str(row), ', ', int2str(col), ')*res(', int2str(pivotRow), ', :);\n']);
              end
            end
          end
        end
      end
      
      fprintf(fid, 'end\n');
    end
    
    
    fclose(fid);
end
