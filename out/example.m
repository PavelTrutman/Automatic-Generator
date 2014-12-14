% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [x, y] = example(a0, a1, a2, a3, a4, b0, b1, b2, b3, b4)

	% precalculate polynomial equations coefficients
	c(1) = a2;
	c(2) = a0;
	c(3) = a3;
	c(4) = a1;
	c(5) = a4;
	c(6) = -b2;
	c(7) = b0;
	c(8) = b3;
	c(9) = b1;
	c(10) = b4;

	M = zeros(6, 12);
	ci = [31, 51];
	M(ci) = c(1);

	ci = [1, 15];
	M(ci) = c(2);

	ci = [49, 63];
	M(ci) = c(3);

	ci = [43, 57];
	M(ci) = c(4);

	ci = [61, 69];
	M(ci) = c(5);

	ci = [8, 29, 36, 52];
	M(ci) = c(6);

	ci = [2, 17, 24, 40];
	M(ci) = c(7);

	ci = [26, 47, 54, 64];
	M(ci) = c(8);

	ci = [20, 41, 48, 58];
	M(ci) = c(9);

	ci = [44, 59, 66, 70];
	M(ci) = c(10);


	M = rref(M);

	A = zeros(6);
	amcols = [12 11 10 9 8 6];
	A(1, 3) = 1;
	A(2, 5) = 1;
	A(3, :) = -M(6, amcols);
	A(4, :) = -M(5, amcols);
	A(5, :) = -M(4, amcols);
	A(6, :) = -M(2, amcols);

	[V D] = eig(A);
	sol =  V([3, 2],:)./(ones(2, 1)*V(1,:));

	if (find(isnan(sol(:))) > 0)
		
		x = [];
		y = [];
	else
		
		I = find(not(imag( sol(1,:) )));
		x = sol(1,I);
		y = sol(2,I);
	end
end
