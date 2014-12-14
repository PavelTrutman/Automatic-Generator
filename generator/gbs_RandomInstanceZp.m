% create random instance of the equation in Zp
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, oct 2014

function [eqi] = gbs_RandomInstanceZp(cfg, eq, known, unknown)

    vars = [unknown known];
    for mon = vars
        eval(['syms ' char(mon) ';']);
    end

    for mon = known
        eval([char(mon) '=' num2str(floor(rand()*cfg.ZpGeneratorMaxNumber)) ';']);
    end

    fprintf('Creating random instance in Z_%d for %d equations\n', cfg.prime, length(eq));
    
    t = clock();
    for i=1:length(eq)

        st = char(expand(eq(i)));
        
        % split to + and - and evaluate separately
        ep = strfind(st, '-');
        em = strfind(st, '+');
        els = [ep em length(st)+1];
        els = sort(els);
        
        p1 = 0;
        prev = 1;
        for j=1:size(els, 2)

            if prev < (els(j)-1)
                if st(prev) == '+'
                    prev = prev+1;
                end
                ee = eval( st( prev:(els(j)-1)) );
                p1 = EquationModulus(p1 + ee, unknown, cfg.prime);
            end
            
            prev = els(j);
            
            dt = etime(clock, t);
            if (dt > 2)

                fprintf(' %d of %d equations instanced [subterm %d / %d]\n', i, length(eq), j, size(els, 2));
                t = clock();
            end            
        end
        
        eqi(i) = p1;
        %eqi2(i) = eval(eq(i));
    end
end