% log recovery parser
% (GBsolver subroutine)
% by Martin Bujnak, sep2008

function [used last max_deg] = gbs_ResetFromLog(filename)

    fid = fopen(filename, 'r');

    used = [];
    max_deg = 0;
    last = 9999999;
    while (true)
        tline = fgets(fid);
        if isnumeric(tline) && tline == -1
            break;
        end

        [eq cnt err] = sscanf(tline, '...removing equation %d failed *\n');

        if (~isempty(err))
            [eq cnt err] = sscanf(tline, '...removing equation %d \t[t:%fsec]\tfailed *\n');
            if isempty(err)
                err = [];
                eq = eq(1);
            end
        end

        if (isempty(err))

            used = [used;eq];
            last = min(last, eq);
        end

        [eq cnt err] = sscanf(tline, '...removing equation %d succeeded\n');
        if (~isempty(err))
            [eq cnt err] = sscanf(tline, '...removing equation %d \t[t:%fsec]\tsucceeded\n');
            if isempty(err)
                err = [];
                eq = eq(1);
            end
        end
        if (isempty(err))

            last = min(last, eq);
        end

        [deg cnt err] = sscanf(tline, '...max. poly. degree %d\n');
        if (isempty(err))

            max_deg = deg;
        end

        pos = strfind(tline, '...coefficient matrix ');
        if ~isempty(pos)
        
            [deg cnt err] = sscanf(tline(pos:end), '...coefficient matrix reallocaction, adding %d(%d degree) monomials\n');
            if (isempty(err))

                max_deg = deg(2);
            end
        end
    end
    
    fclose(fid);
end