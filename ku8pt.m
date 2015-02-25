g1 = transpose(gbs_Vector('g1', 7));
g2 = transpose(gbs_Vector('g2', 7));
g3 = transpose(gbs_Vector('g3', 7));
g4 = transpose(gbs_Vector('g4', 7));
g5 = transpose(gbs_Vector('g5', 7));
g6 = transpose(gbs_Vector('g6', 7));
g7 = transpose(gbs_Vector('g7', 7));
g8 = transpose(gbs_Vector('g8', 7));
syms f31 f32 k;
mon = [k*f31,k*f32,k^2,f31,f32,k,1].';

% parametrization of the fundamental matrix with three unknowns
f11 = -g1*mon;
f12 = -g2*mon;
f13 = -g6*mon;
f21 = -g3*mon;
f22 = -g4*mon;
f23 = -g8*mon;

% fundamental matrix
F = [f11 f12 f13; f21 f22 f23; f31 f32 1];

% three polynomial equations
eq(1) = -(k*g6)*mon + g5*mon;
eq(2) = -(k*g8)*mon + g7*mon;
eq(3) = det(F);

% three unknowns
unknown = {'f31' 'f32' 'k'};

% known parameters
vars = transpose([g1(:); g2(:); g3(:); g4(:); g5(:); g6(:); g7(:); g8(:)]);
known = {};
for var = vars
  known = [known {char(var)}];
end
% call code generator
[res export] = gbs_CreateCode('ku8pt', eq, known, unknown);