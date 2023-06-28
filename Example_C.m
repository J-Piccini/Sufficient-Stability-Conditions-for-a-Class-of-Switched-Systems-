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
% return
%%
close all;
% g = 1e2;
i = 0;
r = c^.5;
X = -r * 2 : 0.01 : r * 2;
X = -1.41 : 0.01 : 1.41;
v = zeros(length(X), 1);
w = zeros(length(X), 1);
syms x2
for x1 = X
    i = i + 1;

    M = @(x2) ...
        8.7957*x1^4 + 1.8977*x1^3*x2 + 17.4811*x1^2*x2^2 ...
        - 1.5706*x1*x2^3 + 9.3477*x2^4 ...
        -g;
    v(i) = fsolve(M, 10, fsolve_opt);
    w(i) = fsolve(M, -10, fsolve_opt);
end

syms z1 z2
M = 8.7957*z1^4 + 1.8977*z1^3*z2 + 17.4811*z1^2*z2^2 ...
    - 1.5706*z1*z2^3 + 9.3477*z2^4 ...
    -g;

[V1, W1] = meshgrid(linspace(-3.15, 3.15), linspace(-2.25, 2.5));
MM = 8.7957*V1.^4 + 1.8977*V1.^3.*W1 + 17.4811*V1.^2 .* W1.^2 ...
    - 1.5706*V1.*W1.^3 + 9.3477*W1.^4 ...
    -g;

figure;
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
contour(V1, W1, MM, [0, 0], 'LineWidth', 1)

f1_sim = @(t, x) [-x(1) - x(2) + 1; x(1) - x(2) + 1];
f2_sim = @(t, x) [-x(1) - x(2) - 1; x(1) - x(2) + 1];
f3_sim = @(t, x) [-x(1) - x(2) + 1; x(1) - x(2) - 1];
tspan = 10;

for i = 1 : 30
    i_select = ceil(10 + (size(v) - 10) * rand);
    i_select = i_select(1);
    X_select = X(i_select);
    v_select = v(i_select);
    w_select = w(i_select);
    x1 = 2* rand - 1;
    s = sign(randn);
    X1 = [X_select, w_select + .2 * rand];
    X1b = [X_select, v_select - .2 * rand];
    [~, x] = ode45(f1_sim, [0 50], X1);
    plot(x(:, 1), x(:, 2), 'r--')
    [~, x] = ode45(f1_sim, [0 50], X1b);
    plot(x(:, 1), x(:, 2), 'r--')
    [~, x] = ode45(f2_sim, [0 50], X1);
    plot(x(:, 1), x(:, 2), 'k--')
    [~, x] = ode45(f2_sim, [0 50], X1b);
    plot(x(:, 1), x(:, 2), 'k--')
    [~, x] = ode45(f3_sim, [0 50], X1);
    plot(x(:, 1), x(:, 2), '--', 'color', [0 0.4470 0.7410])
    [~, x] = ode45(f3_sim, [0 50], X1b);
    plot(x(:, 1), x(:, 2), '--', 'color', [0 0.4470 0.7410])
    clearvars x
end
% axis tight

[V, W] = meshgrid(linspace(-3, 3), linspace(-2.5, 2.5));
Z1 = V.^2 + (W - 1).^2;
Z2 = (V + 1).^2 + W.^2;
Z3 = (V - 1).^2 + W.^2;

contour(V, W, Z1, [2, 2], '-.', 'LineWidth', 1)
contour(V, W, Z2, [4, 4], '-.', 'LineWidth', 1)
contour(V, W, Z3, [4, 4], '-.', 'LineWidth', 1)
return
