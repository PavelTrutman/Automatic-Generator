% create random instance of the equation in Zp
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, oct 2014

function [eqi] = gbs_RandomInstanceZp(cfg, eq, known, unknown)

    vars = [unknown known];
    for mon = vars
        eval(['syms ' char(mon) ';']);
    end

    random = sym(randi([1 cfg.ZpGeneratorMaxNumber - 1], 1, length(known)));

    fprintf('Creating random instance in Z_%d for %d equations\n', cfg.prime, length(eq));

    reverseStr = '';
    for i=1:length(eq)
        msg = sprintf('  instancing %d equation of %d', i, length(eq));
        fprintf([reverseStr msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        [C, M] = coeffs(expand(subs(eq(i), known, random)), unknown);
        C = mod(C, cfg.prime);
        eqi(i) = sum(C.*M);
    end
    fprintf(reverseStr);
end
