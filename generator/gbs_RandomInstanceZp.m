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
    
    reverseStr = '';
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

            % print status
            msg = sprintf('  instancing %d equation of %d (%2.0f %%%%)', i, length(eq), j/ size(els, 2)*100);
            fprintf([reverseStr msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg) - 1);
          
        end
        
        eqi(i) = p1;
        %eqi2(i) = eval(eq(i));
    end
    fprintf(reverseStr);
end