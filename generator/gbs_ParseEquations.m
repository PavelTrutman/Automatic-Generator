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
        [coefs_used, monoms_used] = coeffs(eq(i), unknown);
        [coefs_rc_used, monoms_rc_used] = coeffs(eqs(i), unknown);

        monoms = union(monoms_used, monoms_rc_used);
        monoms = sort(monoms, 'descend');
        coefs = sym(zeros(size(monoms)));
        coefs_rc = sym(zeros(size(monoms)));
        for j = 1:size(monoms, 2)
          for k = 1:size(monoms_used, 2)
            if isequal(monoms(j), monoms_used(k))
              coefs(j) = coefs_used(k);
            end
          end
          for k = 1:size(monoms_rc_used, 2)
            if isequal(monoms(j), monoms_rc_used(k))
              coefs_rc(j) = coefs_rc_used(k);
            end
          end
        end

        p{i}.coefs = coefs;
        p{i}.coefs_rc = coefs_rc;
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
