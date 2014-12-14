% find algebra B basis (macaulay2 GB solver)
% (GBsolver subroutine)
%
% by Martin Bujnak, mar 2008

function [algB res] = gbs_findAlgB_macaulay(cfg, eq, known, unknown)

    algB = [];
    res = 0;
    
    mcExecCode = 'gbsMacaulay\calc.bat';
    mcCode = 'gbsMacaulay\code.m2';
    mcCodeTemplate = 'gbsMacaulay\code_template.m2';

    fprintf('calculating Grobner basis using Macaulay2 \n');

    if ~isfield(cfg, 'eqinstance')

        eqs = cfg.InstanceGenerator(cfg, eq, known, unknown);
    else
        eqs = cfg.eqinstance;
    end

    if ~isempty(cfg.GBrestrictEq)
        eqs = eqs(cfg.GBrestrictEq);
    end
    
    % subst unknowns to x_0 ... x_N
    unkcnt = length(unknown);
    eqcnt = length(eqs);
    mcideal = 'f = (';
    mceqs = '';
    for i = 1:eqcnt

        eq = Monomial_PerformSubst(unknown, {}, char(eqs(i)));
        
        if eq == '0' 
            continue;
        end
        
        % 
        mceqs = [mceqs 'f' int2str(i) '=' eq ';\r\n'];
        mcideal = [mcideal 'f' int2str(i)];
        if i < eqcnt
            mcideal = [mcideal ' || '];
        end
    end
    mcideal = [mcideal ' );'];
    
    
    % prepare macaylay code
    
    repm = [{'$PRIME$'} {int2str(cfg.prime)}];
    repm = [repm; {'$UNKCNT$'} {int2str(unkcnt)} ];
    repm = [repm; {'$UNKNOWN$'} {c2s(unknown, ', ')} ];
    repm = [repm; {'$EQUATIONS$'} {mceqs} ];
    repm = [repm; {'$EQUATIONSVEC$'} {mcideal} ];

    % create output file
    fid = fopen(mcCode, 'w');

    templ = ReadTemplate(mcCodeTemplate);
    for re = repm'
        
        templ = strrep(templ, re{1}, re{2});
    end
    
    fprintf(fid, templ);
    fclose(fid);    

    % call macaulay
    
    [res, val] = dos(mcExecCode);

    dimpos = strfind(val, 'dim:');
    degpos = strfind(val, 'deg:');
    if isempty(dimpos) || isempty(degpos)

        fprintf('Error executing Macaulay2. Check paths.\n');
        return;
    end
    
    % parse dimension 
    dim = sscanf(val(dimpos:dimpos+10), 'dim:\r\n%d');
    
    % parse number of solutions
    deg = sscanf(val(degpos:degpos+10), 'deg:\r\n%d');
    
    if dim < 0
        
        fprintf('WARNING ! non-zero-dimensional variety (select a subset of equations and try again).\n');
        
    elseif dim > 0
        
        fprintf('WARNING ! expected module to be a finite dimensional module.\n');

    else
    
        basispos = strfind(val, '@@@@@@@@@@@@@@');
        basisstr = val((basispos(1)+14):basispos(2)-1);
        sep = strfind(basisstr, '|');
        basisstr = basisstr(sep(1)+1:sep(2)-1);        

        % add * between terms
        xs = strfind(basisstr, 'x');
        res = '';
        prev = 1;
        for xidx = xs
            
            star = '';
            if (xidx > 0) && (basisstr(xidx-1) ~= ' ')
                star = '*';
            end
            res = [res basisstr(prev:(xidx-1)) star];
            prev = xidx;
        end
        res = [res basisstr(prev:end)];

        % rename x_i to original names
        basisstr = Monomial_PerformSubst({}, unknown, res);
        
        % create string list...
        algB = {};
        fprintf('...algebra B basis :\n   ');

        rem = basisstr;
        i = 1;
        while true
            
            [t, rem] = strtok(rem, ' ');
            if isempty(t),  break;  end
            
            algB{i} = t;
            fprintf('%s ', t);
            i = i + 1;
        end
        
        fprintf('\n');
        
        fprintf('...system is expected to have %d solutions\n', deg);
        
    end
end
