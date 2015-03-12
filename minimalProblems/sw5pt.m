function [codename, eq, known, unknown, kngroups, cfg, algB] = sw5pt()
%%
% 5point relative pose problem
% http://cmp.felk.cvut.cz/minimal/5_pt_relative.php
codename = 'sw5pt';

%%
% register the generator
setpaths;

%% 
% create symbolic matrices 
% these matrices represent basis of null space of linearized equation y'*E*x = 0
% matrices E1, E2, E3, E4 will be treated as known
E1 = gbs_Matrix('E1%d%d', 3 ,3);
E2 = gbs_Matrix('E2%d%d', 3 ,3);
E3 = gbs_Matrix('E3%d%d', 3 ,3);
E4 = gbs_Matrix('E4%d%d', 3 ,3);

% create essential matrix as a   linear combination of the null space
% (symbolicaly - variables x, y, z are unknown)
syms x y z;
E = x*E1 + y*E2 + z*E3 +  E4;

% build equations
%  1. essential matrix is singular
eq(1) = det(E);

%  2. essential matrix "trace" constrain - two singular values are equal
%  and one zero (note that this constrain ensures essential matrix
%  sinularity too. However, adding the equation comming from determinant 
%  causes that generator does not need to generate new polynomials in order
%  to find solution.
% 
Et = transpose(E);
te = 2*(E*Et)*E - trace(E*Et)*E;
eq(2:10) = te(:);

%%
% collect known & unknown variables
unknown = {'x' 'y' 'z'};
vars = transpose([E1(:); E2(:); E3(:); E4(:)]);
known = {};
for var = vars
    known = [known {char(var)}];
end

%%
% specify which "known" variables should be grouped into a single input argument (as
% a vector). "kngroups" is a vector where kngroups(k) = l says that k-th
% known variable should be grouped into l-th vector.
kngroups = ones(9,1)*[1 2 3 4];

% optionaly you can call the code generator with 
% different configuration (cfg) or custom basis of algebra A (algB).
cfg = gbs_InitConfig();
algB = [];

end