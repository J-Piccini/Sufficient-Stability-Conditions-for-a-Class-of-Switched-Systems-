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
% V

[V, W] = meshgrid(linspace(-3.8, 3.8), linspace(-9.5, 9.5));

Z = 0.0090473*V.^8 + 2.7791e-13*V.^7.*W + 0.0060313*V.^6.*W.^2 + 0.00084595*V.^5.*W.^3 ...
    - 0.00071884*V.^4.*W.^4 - 9.8373e-06*V.^3.*W.^5 + 0.00048711*V.^2.*W.^6 - 7.7045e-05*V.*W.^7 ...
    + 7.5971e-06*W.^8 + 0.062842*V.^6 + 0.045846*V.^5.*W + 0.075806*V.^4.*W.^2 ...
    - 0.0053731*V.^3.*W.^3 + 0.0077372*V.^2.*W.^4 - 0.005043*V.*W.^5 + 0.0010251*W.^6 ...
    - 0.59898*V.^4 + 0.065958*V.^3.*W - 0.38174*V.^2.*W.^2 + 0.16454*V.*W.^3 ...
    - 0.047215*W.^4 + 3.4372*V.^2 - 4.3946*V.*W + 2.026*W.^2 ...
    -g;

figure;
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
contour(V, W, Z, [0, 0], 'LineWidth', 1)

X = linspace(-2.8, 2.8);
v = zeros(length(X), 1);
w = zeros(length(X), 1);
i = 0;
for x1 = X
    i = i + 1;

    M = @(x2) ...
        0.0090473*x1^8 + 2.7791e-13*x1^7*x2 + 0.0060313*x1^6*x2^2 + 0.00084595*x1^5*x2^3 ...
        - 0.00071884*x1^4*x2^4 - 9.8373e-06*x1^3*x2^5 + 0.00048711*x1^2*x2^6 ...
        - 7.7045e-05*x1*x2^7 + 7.5971e-06*x2^8 + 0.062842*x1^6 + 0.045846*x1^5*x2 ...
        + 0.075806*x1^4*x2^2 - 0.0053731*x1^3*x2^3 + 0.0077372*x1^2*x2^4 - 0.005043*x1*x2^5 ...
        + 0.0010251*x2^6 - 0.59898*x1^4 + 0.065958*x1^3*x2 - 0.38174*x1^2*x2^2 ...
        + 0.16454*x1*x2^3 - 0.047215*x2^4 + 3.4372*x1^2 - 4.3946*x1*x2 + 2.026*x2^2 ...
        -g; 
    v(i) = fsolve(M, 10, fsolve_opt);
    w(i) = fsolve(M, -10, fsolve_opt);
end

f1_p = @(t, x) [x(2); -x(1) - x(2) * (x(1)^2 - 1)];
f2_p = @(t, x) [0, 1; -6, -2] * [x(1); x(2)];

for i = 1 : 50
    i_select = ceil(10 + (size(v) - 10) * rand);
    i_select = i_select(1);
    X_select = X(i_select);
    v_select = v(i_select);
    w_select = w(i_select);

    X1 = [X_select, w_select + .2 * rand];
    X1b = [X_select, v_select - .2 * rand];
    [~, x] = ode45(f1_p, [0 150], X1);
    plot(x(:, 1), x(:, 2), 'r--')
    [~, x] = ode45(f1_p, [0 150], X1b);
    plot(x(:, 1), x(:, 2), 'r--')
    [~, x] = ode45(f2_p, [0 150], X1);
    plot(x(:, 1), x(:, 2), 'k--')
    [~, x] = ode45(f2_p, [0 150], X1b);
    plot(x(:, 1), x(:, 2), 'k--')
end
