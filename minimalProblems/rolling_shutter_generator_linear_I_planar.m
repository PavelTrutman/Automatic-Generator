% 
% P6P rolling shutter (eliminated)
function [eq, known, unknown, kngroups, cfg, algB] = minimal_rolling_shutter_generator_linear_I_planar()

g1 = transpose(gbs_Vector('g1', 8));
g2 = transpose(gbs_Vector('g2', 8));
g3 = transpose(gbs_Vector('g3', 8));
g4 = transpose(gbs_Vector('g4', 8));
g5 = transpose(gbs_Vector('g5', 8));
g6 = transpose(gbs_Vector('g6', 8));



syms v_1 v_2 v_3 w_1 w_2 w_3;


mon = [v_3*w_3 v_1 v_2 v_3,  w_1, w_2, w_3, 1];

cfg = gbs_InitConfig();

clear eq;


% g1(1) = 1
% g2(1)=0
% g3(1)=0
% g4(1)=1
% g5(1)=0
% g6(1)=0
% 
% g1(5)=0
% g2(5)=0
% g3(5)=0
% g4(5)=0
% g5(5)=0
% g6(5) = 1;
% 
% g1(6)=0
% g2(6)=0
% g3(6)=0
% g4(6)=0
% g5(6)=-1
% g6(6)=0
% 
% g1(7) = 0
% g2(7) = -1
% g3(7) = 1
% g4(7) = 0
% g5(7)=0
% g6(7)=0


eq(1) = v_1*w_1 + g1*transpose(mon);
eq(2) = v_1*w_2 + g2*transpose(mon);
eq(3) = v_2*w_1 + g3*transpose(mon);
eq(4) = v_2*w_2 + g4*transpose(mon);
eq(5) = v_3*w_1 + g5*transpose(mon);
eq(6) = v_3*w_2 + g6*transpose(mon);


% g1 = transpose(gbs_Vector('g1', 8));
% g2 = transpose(gbs_Vector('g2', 8));
% g3 = transpose(gbs_Vector('g3', 8));
% g4 = transpose(gbs_Vector('g4', 8));
% g5 = transpose(gbs_Vector('g5', 8));
% g6 = transpose(gbs_Vector('g6', 8));

unknown = {'v_1' 'v_2' 'v_3' 'w_1' 'w_2' 'w_3'};
vars = transpose([g1(:); g2(:); g3(:); g4(:); g5(:); g6(:)]);
known = {};
for var = vars
    known = [known {char(var)}];
end    

% create symbolic vars
for mon = unknown
    eval(['syms ' char(mon) ';']);
end

%cfg.eqinstance = R6P_planar_inst(cfg);

% call code generator
kngroups = ones(8,1)*[1 2 3 4 5 6];
%[res export] = gbs_CreateCode('p6p_rs_lin_I_planar', eq, known, unknown, kngroups);
algB = [];

%[A symcoefs] = rsSolver('imu3pr_peieg.m', eq, 'tan', unknown, known, kngroups);
end