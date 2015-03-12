function [codename, eq, known, unknown, kngroups, cfg, algB] = sw6pt()
%%
% 6point focal length problem
% http://cmp.felk.cvut.cz/minimal/6_pt_relative.php
codename = 'sw6pt';

%%
setpaths;

%%
% create symbolic matrices 
% these matrices represent basis of null space of linearized equation y'*E*x = 0
% matrices F1, F2, F3 will be treated as known
F1 = gbs_Matrix('F1%d%d', 3 ,3);
F2 = gbs_Matrix('F2%d%d', 3 ,3);
F3 = gbs_Matrix('F3%d%d', 3 ,3);

% create fundamental matrix as a linear combination of the null space,
% w = 1/focal^2
% (symbolicaly - variables x, y, w are unknown)
syms x y w real;
F = x*F1 + y*F2 + F3;

% build equations
%  1. fundamental matrix is singular
eq(1) = det(F);

Q = diag([1 1 w]);

% calibrate fundamental matrix using unknown focal length and use "trace"
% constrain for obtained essential matrix
Ft = transpose(F);
te = 2*(F*Q*Ft*Q)*F - trace(F*Q*Ft*Q)*F;
eq(2:10) = te(:);

%%
% collect known & unknown variables
unknown = {'x' 'y' 'w'};
vars = transpose([F1(:); F2(:); F3(:)]);
known = {};
for var = vars
    known = [known {char(var)}];
end

%%

% define variable groups (optional)
kngroups = ones(9,1)*[1 2 3];

% define configuration
cfg = gbs_InitConfig();

% no algB yet computed
algB = [];

end