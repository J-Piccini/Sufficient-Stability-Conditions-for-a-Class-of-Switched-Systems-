close all;
clear;
clc;
fsolve_opt = optimset('Display', 'off');
solver_opt.solver = 'sedumi';

A = [-1, -1; 1, -1];
b1 = [1; 1]; 
b2 = [-1; 1];
b3 = [1; -1];

pvar x1 x2 real
f1 = A * [x1; x2] + b1;
f2 = A * [x1; x2] + b2;
f3 = A * [x1; x2] + b3;

x = [x1; x2];
prog = sosprogram(x);

l = 2;
vec_V = monomials(x, 2 * l);
vec_p = monomials(x, 2 : 2 : 2);

[prog, V] = sospolymatrixvar(prog, vec_V,[1 1]);
[prog, p1] = sospolymatrixvar(prog, vec_p,[1 1]);
[prog, p2] = sospolymatrixvar(prog, vec_p,[1 1]);
[prog, p3] = sospolymatrixvar(prog, vec_p,[1 1]);

dV = jacobian(V, x);

F1 = dV * f1;
F2 = dV * f2;
F3 = dV * f3;

delta = 1;

prog = sosineq(prog, p1);
prog = sosineq(prog, p2);
prog = sosineq(prog, p3);

c = 2;
prog = sosineq(prog, -F1 - delta * x'.^l * x.^l - p1 * (x' * x - c));
prog = sosineq(prog, -F2 - delta * x'.^l * x.^l - p2 * (x' * x - c));
prog = sosineq(prog, -F3 - delta * x'.^l * x.^l - p3 * (x' * x - c));

prog = sosineq(prog, V - delta * x'.^l * x.^l);
[prog, info] = sossolve(prog, solver_opt);
V = sosgetsol(prog,V);

prog = sosprogram(x);
vec_2 = monomials(x, 2 : 2 : 12);
[prog, p4] = sospolymatrixvar(prog, vec_2, [1 1]);
% g = 2.21e2;
g = 38.43;
prog = sosineq(prog, p4);
prog = sosineq(prog, -(V - g) + p4 * (x' * x - c));
[prog, info1] = sossolve(prog,solver_opt);

fprintf('\nFirst SOS:  \n     feasratio: %4.2f\n     pinf %d \n     numerr %d',[info.feasratio, info.pinf, info.numerr])
fprintf('\nSecond SOS: \n     feasratio: %4.2f\n     pinf %d \n     numerr %d',[info1.feasratio, info1.pinf, info1.numerr])
