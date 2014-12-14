% find algebra B basis (maple GB solver)
% (GBsolver subroutine)
%
% by Martin Bujnak, mar 2008

function [algB res] = gbs_findAlgB_maple(cfg, eq, known, unknown)

    fprintf('calculating Grobner basis\n');
    
%    if ~isfield(cfg, 'eqinstance')
%
%        eqs = cfg.InstanceGenerator(cfg, eq, known, unknown);
%    else
%        eqs = cfg.eqinstance;
%    end

    eqs = gbs_RandomInstance(cfg, eq, known, unknown);
     
    if ~isempty(cfg.GBrestrictEq)
        eqs = eqs(cfg.GBrestrictEq);
    end
    
    maple('with(Groebner)');
    ord = maple(['tdeg(' c2s(unknown, ', ') ')']);
    
    % load equations in maple
    eqstr = '[';
    for i=1:length(eqs)
        if i > 1
            eqstr = [eqstr ',' char(eqs(i))];
        else
            eqstr = [eqstr char(eqs(i))];
        end
    end
    eqstr = [eqstr ']'];
    
    F = maple(eqstr);
    G = maple('gbasis', F, ord);

    str = strrep(char(G), 'matrix([', '');
    str = strrep(str, '])', '');

    GG = maple(str);
    [ns, res] = maple('SetBasis', GG , ord);
    
    if res > 0
        
        fprintf('WARNING ! non-zero-dimensional variety (select a subset of equations and try again).\n');
        algB = [];
    else
        reg = regexp(char(ns), '\[.*?\]', 'match');
        ns = reg{1};
        ns = sym(ns);

        fprintf('...algebra B basis :\n   ');
        
        % create string list...
        algB = {};
        for i=1:length(ns)
        
            algB{i} = char(ns(i));
            fprintf('%s ', algB{i});
        end
        
        fprintf('\n');
    end
end