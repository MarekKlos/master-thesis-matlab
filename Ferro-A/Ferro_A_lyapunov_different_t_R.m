% =========================================================================
% FERRO_A: LYAPUNOV SWEEP OVER SERIES RESISTANCE Rs  (Circuit A)
% =========================================================================
% Rows:    Rs values from Table 5.4 of the thesis
% Columns: lambda_max at t_final = 2 s and t_final = 10 s
%
% Purpose: identify transient-chaos cases (lambda_max > 0 at 2 s,
%          but lambda_max <= 0 at 10 s) across the damping sweep.
% =========================================================================
clear; clc; close all;
global par;

% -------------------------------------------------------------------------
% 1. FIXED PARAMETERS  (Circuit A baseline, Case 3 of the thesis)
% -------------------------------------------------------------------------
U   = 90000;        % Source voltage amplitude [V] (chaotic operating point)
alf = 0.0;
Rm  = 50000;
C   = 0.78E-6;      % Capacitance [F]
q   = 11;
m1  = 4.079845e-04;
mq  = 2.108275e-27;
omg = 100*pi;

y0  = [280; 1530; 0];           % initial conditions [Psi; Uc; phi]

% -------------------------------------------------------------------------
% 2. SWEEP RANGE  (matches Table 5.4 in the thesis)
% -------------------------------------------------------------------------
Rs_vec = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, ...
          1.8, 1.9, 2.0, 3.0, 4.0, 5.0, 7.5, 10, 12.5, 15, ...
          18, 20, 25, 30, 35, 40, 45, 50, 60, 70, ...
          80, 90, 100, 150, 200, 250, 300, 350, 400, 450, ...
          500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, ...
          3500, 4000, 4500, 5000, 10000];
nR = length(Rs_vec);

% -------------------------------------------------------------------------
% 3. LYAPUNOV SETTINGS
% -------------------------------------------------------------------------
stept     = 0.005;              % Gram-Schmidt renormalization step
ss_window = 0.20;               % fraction of trajectory averaged at the end

% -------------------------------------------------------------------------
% 4. RUN SWEEP
% -------------------------------------------------------------------------
lambda_2s  = zeros(nR, 1);
lambda_10s = zeros(nR, 1);

fprintf('Running %d resistance points  (2 Lyapunov runs per point)...\n', nR);
tic;
for i = 1:nR
    Rs = Rs_vec(i);

    % Derived parameters (must be recomputed every time Rs changes)
    p2 = 1/(C*(Rm+Rs)); p1 = U*p2;
    p3 = p2*m1*Rm;      p4 = p2*mq*Rm;
    b2 = Rm/(Rm+Rs);    b1 = U*b2;
    b3 = b2*m1*Rs;      b4 = b2*mq*Rs;

    par = [b1; b2; b3; b4; p1; p2; p3; p4; alf; omg; q];

    % ---- t_final = 2 s ------------------------------------------------
    [~, L2]  = lyapunov(3, @ferro_A_ext, @ode45, 0, stept, 2.0,  y0, 0);
    idx2     = round((1-ss_window)*size(L2,1)) : size(L2,1);
    lambda_2s(i) = max( mean(L2(idx2,:), 1) );

    % ---- t_final = 10 s -----------------------------------------------
    [~, L10] = lyapunov(3, @ferro_A_ext, @ode45, 0, stept, 10.0, y0, 0);
    idx10    = round((1-ss_window)*size(L10,1)) : size(L10,1);
    lambda_10s(i) = max( mean(L10(idx10,:), 1) );

    fprintf('  [%2d/%2d]  Rs = %8.2f Ohm  -->  lambda(2s) = %8.3f   lambda(10s) = %8.3f\n', ...
            i, nR, Rs, lambda_2s(i), lambda_10s(i));
end
fprintf('Sweep finished in %.1f s.\n', toc);

% -------------------------------------------------------------------------
% 5. BUILD TABLE
% -------------------------------------------------------------------------
T_lambda = table(Rs_vec(:), lambda_2s, lambda_10s, ...
                 'VariableNames', {'Rs_Ohm', 'lambda_max_2s', 'lambda_max_10s'});

disp(' ');
disp('==============================================================');
disp(' Maximum Lyapunov exponent vs. series resistance  (Circuit A)');
disp('==============================================================');
disp(T_lambda);

% -------------------------------------------------------------------------
% 6. FLAG TRANSIENT CHAOS
% -------------------------------------------------------------------------
eps_tol = 0.5;     % small positive threshold (numerical noise margin)
is_transient  = (lambda_2s  >  eps_tol) & (lambda_10s <=  eps_tol);
is_persistent = (lambda_2s  >  eps_tol) & (lambda_10s >   eps_tol);

disp('Cases flagged as TRANSIENT chaos (positive at 2 s, decays by 10 s):');
disp(Rs_vec(is_transient).');

disp('Cases flagged as PERSISTENT chaos (positive at both 2 s and 10 s):');
disp(Rs_vec(is_persistent).');

% -------------------------------------------------------------------------
% 7. EXPORT CSV
% -------------------------------------------------------------------------
writetable(T_lambda, 'Circuit_A_lambda_vs_Rs.csv');
fprintf('Results saved to Circuit_A_lambda_vs_Rs.csv\n');

% -------------------------------------------------------------------------
% 8. PLOT
% -------------------------------------------------------------------------
figure('Color', 'w', 'Position', [120 120 700 450]);
semilogx(Rs_vec, lambda_2s,  '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto'); hold on;
semilogx(Rs_vec, lambda_10s, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
yline(0, '--k');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Series resistance R_s [\Omega]', 'FontName', 'Times New Roman', 'FontSize', 13);
ylabel('\lambda_{max}',                  'FontName', 'Times New Roman', 'FontSize', 13);
title('Maximum Lyapunov Exponent vs. Series Resistance (Circuit A)', ...
      'FontName', 'Times New Roman', 'FontSize', 14);
legend({'t_{final} = 2 s', 't_{final} = 10 s', '\lambda = 0'}, ...
       'Location', 'best');

% =========================================================================
% LOCAL FUNCTION:  extended ODE for Lyapunov calculation
% (identical to ferro_A_ext in your main analysis script)
% =========================================================================
function f = ferro_A_ext(t, X)
    global par
    b1=par(1); b2=par(2); b3=par(3); b4=par(4);
    p1=par(5); p2=par(6); p3=par(7); p4=par(8);
    alf=par(9); omg=par(10);

    x1 = X(1); x2 = X(2); x3 = X(3);
    Y  = [X(4), X(7), X(10);
          X(5), X(8), X(11);
          X(6), X(9), X(12)];

    f    = zeros(12,1);
    f(1) = b1*cos(x3+alf) - b2*x2 - b3*x1 - b4*x1^par(11);
    f(2) = p1*cos(x3+alf) - p2*x2 + p3*x1 + p4*x1^par(11);
    f(3) = omg;

    Jac = [ -b3 - par(11)*b4*x1^(par(11)-1),  -b2,  -b1*sin(x3+alf);
             p3 + par(11)*p4*x1^(par(11)-1),  -p2,  -p1*sin(x3+alf);
             0,                                  0,   0              ];

    dY      = Jac * Y;
    f(4:12) = [dY(1,1); dY(2,1); dY(3,1); ...
               dY(1,2); dY(2,2); dY(3,2); ...
               dY(1,3); dY(2,3); dY(3,3)];
end