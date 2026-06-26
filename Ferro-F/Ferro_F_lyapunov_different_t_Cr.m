% =========================================================================
% FERRO_F: LYAPUNOV SWEEP OVER PARALLEL CAPACITANCE Cr  (Circuit B, 4th order)
% =========================================================================
% Rows:    Cr values from Table 5.11 of the thesis
% Columns: lambda_max at t_final = 2 s
%          lambda_max at t_final = 10 s
%          max steady-state |Uc| at t_final = 10 s   (in p.u.)
%
% Purpose: identify transient-chaos cases AND verify the long-term
%          overvoltage levels across the parallel-capacitance sweep.
% =========================================================================
clear; clc; close all;
global par;

% -------------------------------------------------------------------------
% 1. FIXED PARAMETERS  (Circuit B baseline, Case 3 of the thesis)
% -------------------------------------------------------------------------
U   = 33550;        % Source voltage amplitude [V] (chaotic operating point)
alf = 0.0;
Rs  = 0.2;
RN  = 1.25;
Cs  = 3.1522E-6;    % Series capacitance [F] (fixed for this sweep)
q   = 11;
m1  = 0.04320728;
mq  = 1.444531e-18;

% Core-loss polynomial coefficients
h1  = 3.396059e-05;
h3  = -6.906720e-14;
h5  = 1.593581e-22;

% Saturation switch for the core-loss model
p13 = 1.2 * 2.122891e+04;
Ro  = 1 / (h1 + p13*p13*(h3 + h5*p13*p13));

omg = 100*pi;

% Initial conditions [Psi; Uc; UL; phase]
y0  = [80; 25300; 16000; 0];

% Derived parameters that DO NOT depend on Cr
Rz  = Rs + RN;
p2  = 1/(Cs*Rz);
p1  = U*p2;

% -------------------------------------------------------------------------
% 2. SWEEP RANGE  (matches Table 5.11 in the thesis)
% -------------------------------------------------------------------------
Cr_vec = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, ...
          100, 500, 1000] * 1e-6;     % [F]
nCr = length(Cr_vec);

% -------------------------------------------------------------------------
% 3. SETTINGS
% -------------------------------------------------------------------------
stept     = 0.005;
ss_window = 0.20;
options   = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% -------------------------------------------------------------------------
% 4. RUN SWEEP
% -------------------------------------------------------------------------
lambda_2s   = zeros(nCr, 1);
lambda_10s  = zeros(nCr, 1);
Uss_max_10s = zeros(nCr, 1);   % max steady-state |Uc| at t=10 s [p.u.]

fprintf('Running %d Cr points  (2 Lyapunov runs + 1 ODE run per point, 4D system)...\n', nCr);
fprintf('This will take a while — the 4D variational system has 20 equations.\n');
tic;
for i = 1:nCr
    Cr = Cr_vec(i);

    % Parameters that DEPEND on Cr
    b2  = 1/(Rz*Cr);    b1  = U*b2;
    b3  = h1/Cr;        b4  = h3/Cr;       b5 = h5/Cr;
    b6  = m1/Cr;        b7  = mq/Cr;
    p14 = 1.0/(Ro*Cr);

    par = [p1; p2; b1; b2; b3; b4; b5; b6; b7; alf; omg; q; p13; p14];

    % ---- Lyapunov at t_final = 2 s ------------------------------------
    [~, L2]  = lyapunov(4, @ferro_F_ext, @ode15s, 0, stept, 2.0,  y0, 0);
    idx2     = round((1-ss_window)*size(L2,1)) : size(L2,1);
    lambda_2s(i) = max( mean(L2(idx2,:), 1) );

    % ---- Lyapunov at t_final = 10 s -----------------------------------
    [~, L10] = lyapunov(4, @ferro_F_ext, @ode15s, 0, stept, 10.0, y0, 0);
    idx10    = round((1-ss_window)*size(L10,1)) : size(L10,1);
    lambda_10s(i) = max( mean(L10(idx10,:), 1) );

    % ---- Long ODE run for steady-state overvoltage at t = 10 s --------
    sol_long = ode15s(@ferro_F_ode_local, [0, 10.0], y0, options);
    T_long   = sol_long.x(:);
    Y_long   = sol_long.y.';
    [T_long, uidx] = unique(T_long);
    Y_long = Y_long(uidx, :);

    idx_ss = T_long >= 0.8 * T_long(end);
    Uss_max_10s(i) = max(abs(Y_long(idx_ss, 2))) / U;   % [p.u.]

    fprintf('  [%2d/%2d]  Cr = %6.1f uF  -->  lambda(2s) = %7.3f   lambda(10s) = %7.3f   Uss(10s) = %.3f p.u.\n', ...
            i, nCr, Cr*1e6, lambda_2s(i), lambda_10s(i), Uss_max_10s(i));
end
fprintf('Sweep finished in %.1f s.\n', toc);

% -------------------------------------------------------------------------
% 5. BUILD TABLE
% -------------------------------------------------------------------------
Cr_uF    = Cr_vec(:) * 1e6;
T_lambda = table(Cr_uF, lambda_2s, lambda_10s, Uss_max_10s, ...
                 'VariableNames', ...
                 {'Cr_uF', 'lambda_max_2s', 'lambda_max_10s', 'Uss_max_10s_pu'});

disp(' ');
disp('==============================================================');
disp(' Circuit B: lambda_max and steady-state Uc vs. parallel capacitance');
disp('==============================================================');
disp(T_lambda);

% -------------------------------------------------------------------------
% 6. FLAG TRANSIENT CHAOS
% -------------------------------------------------------------------------
eps_tol = 0.5;
is_transient  = (lambda_2s  >  eps_tol) & (lambda_10s <=  eps_tol);
is_persistent = (lambda_2s  >  eps_tol) & (lambda_10s >   eps_tol);

disp('Cases flagged as TRANSIENT chaos (positive at 2 s, decays by 10 s):');
disp(Cr_uF(is_transient));

disp('Cases flagged as PERSISTENT chaos (positive at both 2 s and 10 s):');
disp(Cr_uF(is_persistent));

% -------------------------------------------------------------------------
% 7. EXPORT CSV
% -------------------------------------------------------------------------
writetable(T_lambda, 'Circuit_B_lambda_vs_Cr.csv');
fprintf('Results saved to Circuit_B_lambda_vs_Cr.csv\n');

% -------------------------------------------------------------------------
% 8. PLOT
% -------------------------------------------------------------------------
figure('Color', 'w', 'Position', [120 120 800 700]);

subplot(2,1,1);
semilogx(Cr_uF, lambda_2s,  '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto'); hold on;
semilogx(Cr_uF, lambda_10s, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
yline(0, '--k');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Parallel capacitance C_r [\muF]', 'FontName', 'Times New Roman', 'FontSize', 13);
ylabel('\lambda_{max}',                   'FontName', 'Times New Roman', 'FontSize', 13);
title('Maximum Lyapunov Exponent vs. Parallel Capacitance (Circuit B)', ...
      'FontName', 'Times New Roman', 'FontSize', 14);
legend({'t_{final} = 2 s', 't_{final} = 10 s', '\lambda = 0'}, ...
       'Location', 'best');

subplot(2,1,2);
semilogx(Cr_uF, Uss_max_10s, '-d', 'LineWidth', 1.5, ...
         'MarkerFaceColor', 'auto', 'Color', '#D95319');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Parallel capacitance C_r [\muF]', 'FontName', 'Times New Roman', 'FontSize', 13);
ylabel('U_{ss,max} [p.u.]',               'FontName', 'Times New Roman', 'FontSize', 13);
title('Steady-state capacitor overvoltage at t = 10 s (Circuit B)', ...
      'FontName', 'Times New Roman', 'FontSize', 14);

% =========================================================================
% LOCAL FUNCTIONS  (identical to those in Ferro_F_Analysis.m)
% =========================================================================
function dydt = ferro_F_ode_local(t, y)
    global par
    if abs(y(3)) > par(13)
        iR = par(14) * y(3);
    else
        p  = y(3) * y(3);
        iR = y(3) * (par(5) + p*(par(6) + par(7)*p));
    end
    dydt    = zeros(4,1);
    dydt(1) = y(3);
    dydt(2) = par(1)*cos(y(4)+par(10)) - par(2)*y(2) - par(2)*y(3);
    dydt(3) = par(3)*cos(y(4)+par(10)) - par(4)*(y(2)+y(3)) - iR ...
              - par(8)*y(1) - par(9)*y(1)^par(12);
    dydt(4) = par(11);
end

function f = ferro_F_ext(t, X)
    global par
    p1=par(1);  p2=par(2);  b1=par(3);  b2=par(4);
    b3=par(5);  b4=par(6);  b5=par(7);  b6=par(8);  b7=par(9);
    alf=par(10); omg=par(11); q=par(12); p13=par(13); p14=par(14);

    x1 = X(1); x2 = X(2); x3 = X(3); x4 = X(4);
    n  = 4;
    Y  = reshape(X(n+1:end), n, n);

    if abs(x3) > p13
        iR      = p14 * x3;
        diR_dx3 = p14;
    else
        p       = x3 * x3;
        iR      = x3 * (b3 + p*(b4 + b5*p));
        diR_dx3 = b3 + 3*b4*x3^2 + 5*b5*x3^4;
    end

    f    = zeros(n + n^2, 1);
    f(1) = x3;
    f(2) = p1*cos(x4+alf) - p2*x2 - p2*x3;
    f(3) = b1*cos(x4+alf) - b2*(x2+x3) - iR - b6*x1 - b7*x1^q;
    f(4) = omg;

    Jac        = zeros(4,4);
    Jac(1,3)   = 1;
    Jac(2,2)   = -p2;
    Jac(2,3)   = -p2;
    Jac(2,4)   = -p1*sin(x4+alf);
    Jac(3,1)   = -b6 - q*b7*x1^(q-1);
    Jac(3,2)   = -b2;
    Jac(3,3)   = -b2 - diR_dx3;
    Jac(3,4)   = -b1*sin(x4+alf);

    dY         = Jac * Y;
    f(n+1:end) = dY(:);
end