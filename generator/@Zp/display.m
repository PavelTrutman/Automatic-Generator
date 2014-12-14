
function display(p)
% Zp/DISPLAY Command window display of a Zp number
    disp(' ');
    disp([inputname(1),' = ']);
    disp(' ');
    disp(['   ' int2str(p.c) '(mod ' int2str(p.p) ')']);
    disp(' ');