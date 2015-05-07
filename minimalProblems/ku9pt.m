function [eq, known, unknown, kngroups, cfg, algB] = ku9pt()

g1 = transpose(gbs_Vector('g1', 7));
g2 = transpose(gbs_Vector('g2', 7));
g3 = transpose(gbs_Vector('g3', 7));
g4 = transpose(gbs_Vector('g4', 7));
g5 = transpose(gbs_Vector('g5', 7));
g6 = transpose(gbs_Vector('g6', 7));
g7 = transpose(gbs_Vector('g7', 7));
g8 = transpose(gbs_Vector('g8', 7));
g9 = transpose(gbs_Vector('g9', 7));

syms f31 f32 k1 k2;

mon = [k1*f32,k1*k2,f31,f32,k1,k2,1].';

% parametrization of the fundamental matrix with four unknowns
f11 = -g1*mon;
f12 = -g2*mon;
f13 = -g6*mon;
f21 = -g3*mon;
f22 = -g4*mon;
f23 = -g8*mon;

% fundamental matrix
F = [f11 f12 f13; f21 f22 f23; f31 f32 1];

% four polynomial equations
eq(1) = -(k2*g6)*mon + g5*mon;
eq(2) = -(k2*g8)*mon + g7*mon;
eq(3) = det(F);
eq(4) = k1*f31+g9*mon;

% four unknowns
unknown = {'f31' 'f32' 'k1' 'k2'};

% known parameters
vars = transpose([g1(:); g2(:); g3(:); g4(:); g5(:); g6(:); g7(:); g8(:); g9(:)]);
known = {};
for var = vars
  known = [known {char(var)}];
end

kngroups = [];
cfg = gbs_InitConfig();
algB = [];

end