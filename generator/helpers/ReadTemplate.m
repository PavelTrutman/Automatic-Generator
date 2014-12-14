% read text file to string
%
% by Martin Bujnak, oct 2005

function [str] = ReadTemplate(filename);

    fid = fopen(filename, 'r');
    
    str = [];
    while (fid >= 0)
       
        tline = fgets(fid);
        if isnumeric(tline) && tline == -1
            break;
        end
        str = [str tline];
    end
    
    fclose(fid);
end