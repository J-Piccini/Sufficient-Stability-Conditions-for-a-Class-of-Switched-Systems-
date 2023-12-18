close all;
clear;
clc;
fsolve_opt = optimset('Display', 'off');

solver_opt.solver = 'sedumi';
pvar x1 x2 real

A1 = [0, 1;  -1, -(x1^2 - 1)];
A2 = [0, 1; -6, -2];

x = [x1; x2];
prog = sosprogram(x);

l = 1;
vec = monomials(x, 2 : 2 : 8);
vec1 = monomials(x, 2 : 2 : 2);

[prog, V] = sospolymatrixvar(prog, vec, [1 1]);
[prog, p1] = sospolymatrixvar(prog, vec1, [1 1]);
[prog, p2] = sospolymatrixvar(prog, vec1, [1 1]);

dV = jacobian(V, x);
delta = 1e-4;
F1 = dV * A1 * x; % f(x) * \nabla V_u(x)
% F1 = dV * A1;tt_times
F2 = dV * A2 * x;

c = 1.4e1;

prog = sosineq(prog, p1);
prog = sosineq(prog, p2);
prog = sosineq(prog, -F1 - delta * x'.^l * x.^l - p1 * (x' * x - c)); % (16) in paper
prog = sosineq(prog, -F2 - delta * x'.^l * x.^l - p2 * (x' * x - c));

prog = sosineq(prog, V - delta * x'.^l * x.^l); % Radially unbounded (16), with \delta = 1)

[prog, info1] = sossolve(prog, solver_opt);

V = sosgetsol(prog, V);
% return; % Safe so far
prog = sosprogram(x);

vec12 = monomials(x, 2 : 2 : 8);
[prog, p3] = sospolymatrixvar(prog, vec12, [1 1]);
g = 4.65e2;
prog = sosineq(prog, p3);
prog = sosineq(prog, -(V - g) + p3 * (x' * x - c));
[prog,info] = sossolve(prog, solver_opt);
close all;

fprintf('\nFirst SOS:  \n     feasratio: %4.2f\n     pinf %d \n     numerr %d',[info1.feasratio, info1.pinf, info1.numerr])
fprintf('\nSecond SOS: \n     feasratio: %4.2f\n     pinf %d \n     numerr %d',[info.feasratio, info.pinf, info.numerr])
