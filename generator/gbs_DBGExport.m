% debug export (equations, variables, instance, polynomials,...)
% (GBsolver subroutine)
% by Martin Bujnak, oct2008

function [] = gbs_DBGExport(known, unknown, eq, symcoefs, allmons, Minit, Mmons)

    % export code
    fid = fopen('coefs.txt', 'w');

    fprintf(fid, '%% Generated using GBSolver generator Copyright Martin Bujnak,\n'); 
    fprintf(fid, '%% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.\n%% \n');
    fprintf(fid, '%% Please refer to the following paper, when using this code :\n'); 
    fprintf(fid, '%%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,\n');
    fprintf(fid, '%%     ECCV 2008, Marseille, France, October 12-18, 2008\n');
    fprintf(fid, '%%\n');

    fprintf(fid, '%% unknowns\n\t%s', c2s(unknown, ' '));

    fprintf(fid, '\n\n');

    fprintf(fid, '%% known\n\t%s', c2s(known, ' '));

    fprintf(fid, '\n\n');        
    % export equations in raw form
    fprintf(fid, '%% raw equations \n');
    for i=1:length(eq)
        fprintf(fid, '\teq(%d) = %s;\n', i, char(eq(i)));
    end        
    fprintf(fid, '\n\n');

    % export coefs vars
    fprintf(fid, '%% coeficients \n');

    for i=1:length(symcoefs)
        fprintf(fid, '\tc(%d) = %s;\n', i, char(symcoefs(i)));
    end        
    fprintf(fid, '\n\n');

    % export monimials vector
    fprintf(fid, '%% monomials [%d] \n   ', length(allmons));
    fprintf(fid, '\tmons=[...\n');
    for mons = allmons
        fprintf(fid, '%s ', char(mons));
    end
    fprintf(fid, '];\n');

    fprintf(fid, '\n\n');

    % export coefs matrix
    fprintf(fid, '%% coefficients matrix \n');
    fprintf(fid, '\tM=[...\n');
    for i=1:size(Minit,1)
        fprintf(fid, '\t\t');
        for j=1:size(Minit,2)
            fprintf(fid, '%d ', Minit(i,j));
        end
        fprintf(fid, '...\n');
    end
    fprintf(fid, '\t\t];\n');


    % export used monimials vector
    fprintf(fid, '%% used monomials [%d] \n   ', length(monomials));
    fprintf(fid, '\tmons=[...\n');
    for mons = monomials
        fprintf(fid, '%s ', char(mons));
    end
    fprintf(fid, '1];\n');

    fprintf(fid, '\n\n');

    % export coefs matrix
    fprintf(fid, '%% used coefficients matrix \n');
    fprintf(fid, '\tM=[...\n');
    for i=1:size(Mmons,1)
        fprintf(fid, '\t\t');
        for j=1:size(Mmons,2)
            fprintf(fid, '%d ', Mmons(i,j));
        end
        fprintf(fid, '...\n');
    end
    fprintf(fid, '\t\t];\n');        

    fclose(fid);
end