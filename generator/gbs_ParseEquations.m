% Parse input equations
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, oct 2014


function [monomials, p, p_currdeg, symcoefs, maxdeg] = gbs_ParseEquations(cfg, eq, known, unknown)

    if ~isfield(cfg, 'eqinstance')

        eqs = cfg.InstanceGenerator(cfg, eq, known, unknown);
    else
        eqs = cfg.eqinstance;
    end
    
    prime = cfg.prime;
    
    monomials=[];
    eqcnt = length(eq);
    usedcoefidx = 1;
    clear symcoefs;
    clear p;
    p = cell([1 eqcnt]);
    for i=1:eqcnt

        eq(i) = expand(eq(i));
        
        if eq(i) == 0
            
            warning('empty/zero equation (%d), remove it and run solver again', i);
        end

        p{i}.eq = eq(i);

        % extract monomials & coefficents
        [coefs, monoms] = coeffs(eq(i), unknown);
        [coefs_rc, monoms] = coeffs(eqs(i), unknown);
        [monoms, sortArray] = sort(monoms, 'descend');
        p{i}.coefs = coefs(sortArray);
        p{i}.coefs_rc = coefs_rc(sortArray);
        p{i}.mons = arrayfun(@(x) char(x), monoms, 'UniformOutput', 0);
        p{i}.deg = MonomialsDegree(p{i}.mons, unknown);
        p{i}.maxdeg = max( sum(p{i}.deg, 2) );

        p_currdeg(i) = p{i}.maxdeg;
        p{i}.monsidx = 0;
        p{i}.monscnt = length(p{i}.mons);
        monomials = union(monomials, p{i}.mons);

        % convert (fractions) coeffs to Zp
        zpcoefs = [];
        coefidx = [];
        for j=1:length(p{i}.coefs)

           cmpl = p{i}.coefs_rc(j);
           [a]=sscanf(char(cmpl),'%d/%d');
           if length(a) > 1
               zpcoefs(j) = double(Zp(a(1), prime) / Zp(a(2), prime));
           else
               zpcoefs(j) = double(Zp(a(1), prime));
           end
           
           symcoefs(usedcoefidx) = p{i}.coefs(j);
           coefidx(j) = usedcoefidx;
           usedcoefidx = usedcoefidx+1;
        end

        p{i}.coefs = zpcoefs;
        p{i}.coefsIDX = coefidx;

    end 
    maxdeg = max(p_currdeg);
    p_currdeg = p_currdeg - min(p_currdeg);
end
