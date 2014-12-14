% Generate Matlab Code for given action matrix and coefficient matrices
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, oct 2014


function [res] = gbs_ExportMCode(filename, M, trace, coefscode, known, knowngroups, unknown, algB, actMvar, amrows, amcols, gjcols, aidx)

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
      fprintf(fid, '\tM = rref(M);\n');
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
        fprintf(fid, '\n\tM = rref(M);\n\n');  
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
    fclose(fid);
end
