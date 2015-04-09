function [eq, known, unknown, kngroups, cfg, algB] = three_quadratic()

% symbolic variables
syms v1 v2 v3;

g1 = transpose(gbs_Vector('g1', 10));
g2 = transpose(gbs_Vector('g2', 10));
g3 = transpose(gbs_Vector('g3', 10));

mon = [v1^2 v2^2 v3^2 v1*v2 v1*v3 v2*v3 v1 v2  v3 1];

% equations
eq(1) =   g1*transpose(mon);
eq(2) =   g2*transpose(mon);
eq(3) =   g3*transpose(mon);

% unknowns
unknown = {'v1' 'v2' 'v3'};

% known parameters
vars = transpose([g1(:); g2(:); g3(:)]);
known = {};

for var = vars
    known = [known {char(var)}];
end    

% create symbolic vars
for mon = unknown
    eval(['syms ' char(mon) ';']);
end

% kngroups
kngroups = ones(10,1)*[1,2,3];

%config file
cfg = gbs_InitConfig();

% we do not have precomputed algB
algB = [];

end