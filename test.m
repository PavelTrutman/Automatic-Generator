% symbolic variables
syms a0 a1 a2 a3 a4 b0 b1 b2 b3 b4;
syms x y;
% two equations representing an ellipse and a hyperbola
eq(1) = a0*x^3 + a1*x + a2*y^2 + a3*y + a4;
eq(2) = b0*x^2 + b1*x - b2*y^2 + b3*y + b4;
% known parameters
known = {'a0' 'a1' 'a2' 'a3' 'a4' 'b0' 'b1' 'b2' 'b3' 'b4'};
%two unknowns
unknown = {'x' 'y'};
% call code generator
[res, export] = gbs_CreateCode('example', eq, known, unknown);