% Generate Matlab Code for given action matrix and coefficient matrices
% (GBsolver subroutine)
% by Martin Bujnak, sep2008
% last edit by Pavel Trutman, February 2015


function [res] = gbs_ExportMapleCode(filename, M, trace, coefscode, known, knowngroups, unknown, algB, actMvar, amrows, amcols, gjcols, aidx, PaToH)

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

    fprintf(fid, ['> \tlocal c, M, Mold, amcols, A, D1, V1, i, ' c2s((unknown), ', ') ';\n> \n']);
    
    % coeffs
    fprintf(fid, '> \t# precalculate polynomial equations coefficients\n');
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
        
        % replace (,) with [,]
        coefcode = strrep(coefcode, '(', '[');
        coefcode = strrep(coefcode, ')', ']');
        
        fprintf(fid, ['> \tc[' int2str(i) '] := ' coefcode ':\n']);
    end
    fprintf(fid, '> \n');
    
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
      fprintf(fid, '> \tM := ReducedRowEchelonForm(M):\n');
    else
  	  fprintf(fid, ['> \tM := ReducedRowEchelonForm(M[1..-1, [' l2s(trace{1}.nonzerocols, ', ') ']]):\n']);
    end
    fprintf(fid, '> \n');
    
    for i = 2:length(trace)
      fprintf(fid, ['> \tMold := Matrix(' int2str(size(trace{1}.Mcoefs, 1)) ', ' int2str(size(trace{1}.Mcoefs, 2)) ', 0):\n']);
      fprintf(fid, ['> \tMold[1..-1, [' l2s(trace{i - 1}.nonzerocols, ', ') ']] := M:\n']);
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
        fprintf(fid, ['\n> \tM := ReducedRowEchelonForm(M):\n']);
      else
        fprintf(fid, ['\n> \tM := ReducedRowEchelonForm(M[1..-1, [' l2s(trace{i}.nonzerocols, ', ') ']]):\n']);
      end
    end
    
    fprintf(fid, '> \n');
    
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
end
