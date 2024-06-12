close all;
clear;
clc;
%%
fsolve_opt = optimset('Display', 'off');
solver_opt.solver = 'sedumi';

A = [-1, -1; 1, -1];

% Switching signal-dependent disturbances
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

[prog, V] = sospolymatrixvar(prog, vec_V,[1 1]); % Creating a SOS Lyapunov function with polynomials defined in vec_V
[prog, p1] = sospolymatrixvar(prog, vec_p,[1 1]); % Creating SOS p_i polynomials with degrees specified in vecP
[prog, p2] = sospolymatrixvar(prog, vec_p,[1 1]);
[prog, p3] = sospolymatrixvar(prog, vec_p,[1 1]);

dV = jacobian(V, x); % \nabla V

% \dot{V} for the different modes
F1 = dV * f1;
F2 = dV * f2;
F3 = dV * f3;

delta = 1;

prog = sosineq(prog, p1);
prog = sosineq(prog, p2);
prog = sosineq(prog, p3);

c = 2;
% Setting the constraint: -\dot{V}_i - p_i (||x||^2_2 - \beta) - \delta||x||^{2\ell}_{2\ell} Eqn. 2 in the paper
prog = sosineq(prog, -F1 - delta * x'.^l * x.^l - p1 * (x' * x - c));
prog = sosineq(prog, -F2 - delta * x'.^l * x.^l - p2 * (x' * x - c));
prog = sosineq(prog, -F3 - delta * x'.^l * x.^l - p3 * (x' * x - c));

% Requiring V to be radially unbounded
prog = sosineq(prog, V - delta * x'.^l * x.^l);
[prog, info] = sossolve(prog, solver_opt);
V = sosgetsol(prog,V);

% The switched system is bounded, now we minimise the absorbing set
prog = sosprogram(x);
vec_2 = monomials(x, 2 : 2 : 12);
[prog, p4] = sospolymatrixvar(prog, vec_2, [1 1]);
% g = 2.21e2;
g = 38.43;
prog = sosineq(prog, p4);
% Now imposing the constraint: -(V - \gamma) + q(||x||^2_2 - \beta) Eqn. 6 in the paper
prog = sosineq(prog, -(V - g) + p4 * (x' * x - c));
[prog, info1] = sossolve(prog,solver_opt);

fprintf('\nFirst SOS:  \n     feasratio: %4.2f\n     pinf %d \n     numerr %d',[info.feasratio, info.pinf, info.numerr])
fprintf('\nSecond SOS: \n     feasratio: %4.2f\n     pinf %d \n     numerr %d',[info1.feasratio, info1.pinf, info1.numerr])
