% create random instance of the equation in Zp
% (GBsolver subroutine)
% by Martin Bujnak, mar2008

function [eqi] = gbs_RandomInstance(cfg, eq, known, unknown)

    vars = [unknown known];
    for mon = vars
        eval(['syms ' char(mon) ';']);
    end

    for mon = known
        eval([char(mon) '=' num2str(floor(rand()*cfg.ZpGeneratorMaxNumber)) ';']);
    end

    fprintf('Creating random instance for %d equations\n', length(eq));
    
    for i=1:length(eq)
        eqi(i) = eval(eq(i));
    end
end