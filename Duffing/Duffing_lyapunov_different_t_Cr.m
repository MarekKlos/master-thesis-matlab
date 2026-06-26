% =========================================================================
% DUFFING FERRORESONANCE: LYAPUNOV SWEEP OVER C  (Circuit C)
% =========================================================================
% Rows:    C values from Table 5.15 of the thesis
% Columns: lambda_max at t_final = 2 s and t_final = 10 s
%
% Purpose: identify transient-chaos cases (positive at 2 s but decaying by
%          10 s) across the parallel-capacitance sweep.
% =========================================================================
clear; clc; close all;
global par;

% -------------------------------------------------------------------------
% 1. FIXED PARAMETERS  (Circuit C baseline, Case 2 of the thesis)
% -------------------------------------------------------------------------
Im  = 1.5;          % Sinusoidal current source amplitude [A]
I0  = 0.5;          % DC component of the current source [A]
R   = 19.38E+6;     % Core-loss (parallel) resistance [Ohm]
m1  = 0.32;         % Linear coefficient of magnetization curve
m3  = 0.9;          % Cubic coefficient of magnetization curve
omg = 100*pi;       % Angular frequency [rad/s]

y0  = [0; 0; 0];    % Initial conditions [Psi; Uc; phase]

% -------------------------------------------------------------------------
% 2. SWEEP RANGE  (matches Table 5.15 in the thesis)
% -------------------------------------------------------------------------
C_vec = [0.25, 0.5, 0.75, 1.0, 5, 10, 15, 20, 25, 30, ...
         35, 40, 45, 50, 60, 70, 80, 90, 100] * 1e-6;     % [F]
nC = length(C_vec);

% -------------------------------------------------------------------------
% 3. LYAPUNOV SETTINGS
% -------------------------------------------------------------------------
stept     = 0.005;
ss_window = 0.20;

% -------------------------------------------------------------------------
% 4. RUN SWEEP
% -------------------------------------------------------------------------
lambda_2s  = zeros(nC, 1);
lambda_10s = zeros(nC, 1);

fprintf('Running %d C points  (2 Lyapunov runs per point)...\n', nC);
tic;
for i = 1:nC
    C = C_vec(i);
    par = [Im; I0; C; R; m1; m3; omg];

    % ---- t_final = 2 s ------------------------------------------------
    [~, L2]  = lyapunov(3, @duffing_ext, @ode45, 0, stept, 2.0,  y0, 0);
    idx2     = round((1-ss_window)*size(L2,1)) : size(L2,1);
    lambda_2s(i) = max( mean(L2(idx2,:), 1) );

    % ---- t_final = 10 s -----------------------------------------------
    [~, L10] = lyapunov(3, @duffing_ext, @ode45, 0, stept, 10.0, y0, 0);
    idx10    = round((1-ss_window)*size(L10,1)) : size(L10,1);
    lambda_10s(i) = max( mean(L10(idx10,:), 1) );

    fprintf('  [%2d/%2d]  C = %8.3f uF  -->  lambda(2s) = %8.3f   lambda(10s) = %8.3f\n', ...
            i, nC, C*1e6, lambda_2s(i), lambda_10s(i));
end
fprintf('Sweep finished in %.1f s.\n', toc);

% -------------------------------------------------------------------------
% 5. BUILD TABLE
% -------------------------------------------------------------------------
C_uF     = C_vec(:) * 1e6;
T_lambda = table(C_uF, lambda_2s, lambda_10s, ...
                 'VariableNames', {'C_uF', 'lambda_max_2s', 'lambda_max_10s'});

disp(' ');
disp('==============================================================');
disp(' Maximum Lyapunov exponent vs. parallel capacitance  (Circuit C)');
disp('==============================================================');
disp(T_lambda);

% -------------------------------------------------------------------------
% 6. FLAG TRANSIENT CHAOS
% -------------------------------------------------------------------------
eps_tol = 0.5;
is_transient  = (lambda_2s  >  eps_tol) & (lambda_10s <=  eps_tol);
is_persistent = (lambda_2s  >  eps_tol) & (lambda_10s >   eps_tol);

disp('Cases flagged as TRANSIENT chaos (positive at 2 s, decays by 10 s):');
disp(C_uF(is_transient));

disp('Cases flagged as PERSISTENT chaos (positive at both 2 s and 10 s):');
disp(C_uF(is_persistent));

% -------------------------------------------------------------------------
% 7. EXPORT CSV
% -------------------------------------------------------------------------
writetable(T_lambda, 'Circuit_C_lambda_vs_C.csv');
fprintf('Results saved to Circuit_C_lambda_vs_C.csv\n');

% -------------------------------------------------------------------------
% 8. PLOT
% -------------------------------------------------------------------------
figure('Color', 'w', 'Position', [120 120 700 450]);
semilogx(C_uF, lambda_2s,  '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto'); hold on;
semilogx(C_uF, lambda_10s, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
yline(0, '--k');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Parallel capacitance C [\muF]', 'FontName', 'Times New Roman', 'FontSize', 13);
ylabel('\lambda_{max}',                 'FontName', 'Times New Roman', 'FontSize', 13);
title('Maximum Lyapunov Exponent vs. Parallel Capacitance (Circuit C)', ...
      'FontName', 'Times New Roman', 'FontSize', 14);
legend({'t_{final} = 2 s', 't_{final} = 10 s', '\lambda = 0'}, ...
       'Location', 'best');

% =========================================================================
% LOCAL FUNCTION:  extended ODE for Lyapunov calculation
% (identical to duffing_ext in your main analysis script)
% =========================================================================
function f = duffing_ext(t, X)
    global par
    Im=par(1); I0=par(2); C=par(3); R=par(4);
    m1=par(5); m3=par(6); omg=par(7);

    n  = 3;
    x1 = X(1); x2 = X(2); x3 = X(3);
    Y  = reshape(X(n+1:end), n, n);

    f    = zeros(n + n^2, 1);
    f(1) = x2;
    f(2) = (1/C)*(Im*cos(x3) + I0) - (1/(R*C))*x2 ...
           - (m3/C)*(x1^3) - (m1/C)*x1;
    f(3) = omg;

    Jac        = zeros(n,n);
    Jac(1,2)   = 1;
    Jac(2,1)   = -(3*m3/C)*x1^2 - (m1/C);
    Jac(2,2)   = -1/(R*C);
    Jac(2,3)   = -(Im/C)*sin(x3);

    dY         = Jac * Y;
    f(n+1:end) = dY(:);
end