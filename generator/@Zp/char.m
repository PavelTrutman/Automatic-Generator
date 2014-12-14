
function s = char(p) 

	% Zp/CHAR   
	% CHAR(p) is the string representation of p.c

	s=[int2str(p.c) '(mod ' int2str(p.p) ')'];
end