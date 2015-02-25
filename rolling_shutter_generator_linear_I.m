% 
% P6P rolling shutter (eliminated)
%clear all;




g1 = transpose(gbs_Vector('g1', 10));
g2 = transpose(gbs_Vector('g2', 10));
g3 = transpose(gbs_Vector('g3', 10));
g4 = transpose(gbs_Vector('g4', 10));
g5 = transpose(gbs_Vector('g5', 10));
g6 = transpose(gbs_Vector('g6', 10));





syms v_1 v_2 v_3 w_1 w_2 w_3;


mon = [  v_3*w_1 v_3*w_2 v_3*w_3 v_1 v_2 v_3,  w_1, w_2, w_3, 1];


clear eq;

% g1(9)=0
% g2(9)=-1
% g3(9)=0
% g4(9)=1
% g5(9)=0
% g6(9)=0


eq(1) = v_1*w_1 + g1*transpose(mon);
eq(2) = v_1*w_2 + g2*transpose(mon);
eq(3) = v_1*w_3 + g3*transpose(mon);
eq(4) = v_2*w_1 + g4*transpose(mon);
eq(5) = v_2*w_2 + g5*transpose(mon);
eq(6) = v_2*w_3 + g6*transpose(mon);


g1 = transpose(gbs_Vector('g1', 10));
g2 = transpose(gbs_Vector('g2', 10));
g3 = transpose(gbs_Vector('g3', 10));
g4 = transpose(gbs_Vector('g4', 10));
g5 = transpose(gbs_Vector('g5', 10));
g6 = transpose(gbs_Vector('g6', 10));

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


% call code generator
kngroups = [];
[res export] = gbs_CreateCode('p6p_rs_lin_I_bezstlpca', eq, known, unknown, kngroups);

%[A symcoefs] = rsSolver('imu3pr_peieg.m', eq, 'tan', unknown, known, kngroups);
