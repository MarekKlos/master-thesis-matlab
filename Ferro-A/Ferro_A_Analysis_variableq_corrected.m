% =========================================================================
% FERRO_A: MAIN ANALYSIS
% =========================================================================
% Circuit A - Ferroresonance Analysis
% Computes: time-domain waveforms, phase portrait, 3D trajectory,
%           FFT spectrum, stroboscopic diagram, Lyapunov exponents.
%
% Reference for Lyapunov algorithm:
%   A. Wolf, J. B. Swift, H. L. Swinney, and J. A. Vastano,
%   "Determining Lyapunov Exponents from a Time Series,"
%   Physica D, Vol. 16, pp. 285-317, 1985.
% =========================================================================
clear; clc; close all;
global par;

disp('Starting Ferro_A Analysis...');

% -------------------------------------------------------------------------
% 1. SYSTEM PARAMETERS
% -------------------------------------------------------------------------
U   = 50000;        % Source voltage amplitude [V]
alf = 0.0;          % Initial phase angle [rad]
%Rs  = 1.5;          % Series resistance [Ohm]
Rs=1.5;
Rm  = 50000;        % Core-loss (magnetizing) resistance [Ohm]
C   = 0.78E-6;      % Capacitance [F]
%C   = 1085E-6;
q   = 11;           % Degree of nonlinearity in flux-current characteristic
m1  = 4.079845e-04; % Linear coefficient of magnetization curve
mq  = 2.108275e-27; % Nonlinear coefficient of magnetization curve
omg = 100*pi;       % Angular frequency [rad/s]  (50 Hz)

% Derived parameters
p2 = 1/(C*(Rm+Rs)); p1 = U*p2;
p3 = p2*m1*Rm;      p4 = p2*mq*Rm;
b2 = Rm/(Rm+Rs);    b1 = U*b2;
b3 = b2*m1*Rs;      b4 = b2*mq*Rs;

%   Index: 1   2   3   4   5   6   7   8    9    10   11
par = [b1; b2; b3; b4; p1; p2; p3; p4; alf; omg; q];

% -------------------------------------------------------------------------
% 2. SIMULATION SETTINGS
% -------------------------------------------------------------------------
t_final = 2.0;          % Total simulation time [s]  (100 cycles at 50 Hz)
tspan   = [0, t_final];
y0      = [280; 1530; 0]; % Initial conditions: [flux Psi; voltage Uc; phase]
%y0      = [200; 300; 0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% -------------------------------------------------------------------------
% 3. RUN MAIN SIMULATION
% -------------------------------------------------------------------------
disp('Running ODE simulation...');
tic;
% Capture the full solution structure so the stroboscopic section can use
% deval (the solver's own consistent interpolant) instead of interp1.
sol = ode45(@ferro_A_ode, tspan, y0, options);
elapsed = toc;
fprintf('ODE simulation completed in %.2f s.\n', elapsed);

% Reconstruct T and Y so all downstream charts work unchanged.
T = sol.x(:);
Y = sol.y.';

x1 = Y(:,1);   % Flux linkage Psi [Wb]
x2 = Y(:,2);   % Capacitor voltage Uc [V]

% Extract last 20% of simulation as steady state
t_cutoff  = 0.8 * T(end);
idx_steady = T >= t_cutoff;
T_steady  = T(idx_steady);
x1_steady = x1(idx_steady);
x2_steady = x2(idx_steady);

T(end)
%T(idx_steady)
% -------------------------------------------------------------------------
% 4. OVERVOLTAGE REPORT
% -------------------------------------------------------------------------
max_transient_x2 = max(abs(x2));
max_steady_x2    = max(abs(x2_steady));
max_transient_x2pu = max(abs(x2))/U;
max_steady_x2pu    = max(abs(x2_steady))/U;

disp('===================================================');
disp('   Capacitor Overvoltages (x2 = Uc)');
disp('===================================================');
fprintf('> Max overvoltage (full transient):  %.2f p.u.\n', max_transient_x2pu);
fprintf('> Max overvoltage (steady state):    %.2f p.u.\n', max_steady_x2pu);
disp('===================================================');

% -------------------------------------------------------------------------
% 5. LYAPUNOV EXPONENTS
% -------------------------------------------------------------------------
disp('Calculating Lyapunov exponents...');
[T_lyap, L_exp] = lyapunov(3, @ferro_A_ext, @ode45, 0, 0.005, t_final, y0, 0);

% Average over the last 20% to suppress fluctuations
idx_lyap_ss = round(0.8 * length(T_lyap)) : length(T_lyap);
lyap_final  = mean(L_exp(idx_lyap_ss, :), 1);
lambda_max  = max(lyap_final);

disp('===================================================');
disp('        Maximum Lyapunov Exponent (MLE)');
disp('===================================================');
fprintf('> lambda_max = %.6f\n', lambda_max);
if lambda_max > 0
    disp('> Result: CHAOTIC behaviour confirmed (lambda_max > 0)');
else
    disp('> Result: Regular (non-chaotic) behaviour (lambda_max <= 0)');
end
disp('===================================================');

% =========================================================================
% 6. FIGURE SETTINGS  (edit here for consistent thesis formatting)
% =========================================================================
FIG_W  = 560;   % Figure width  [px]
FIG_H  = 420;   % Figure height [px]
FSIZE  = 14;    % Axis tick/label font size
FTITLE = 12;    % Title font size
FNAME  = 'Times New Roman';
LW     = 1.5;   % Default line width

COL1   = '#0072BD';   % Blue
COL2   = '#D95319';   % Orange
COL3   = '#77AC30';   % Green

% -------------------------------------------------------------------------
% CHART 1: Time-Domain Waveforms (Full Simulation)
% -------------------------------------------------------------------------
fig1 = figure('Name', '1. Time Domain', 'Color', 'w', ...
              'Position', [100, 100, FIG_W, FIG_H]);

subplot(2,1,1);
plot(T, x2, 'Color', COL1, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Capacitor Voltage U_c', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]',  'FontName', FNAME);
ylabel('U_c [V]',     'FontName', FNAME);
xlim([0, t_final]);
%xlim([0, 1]);

subplot(2,1,2);
plot(T, x1, 'Color', COL2, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Flux Linkage \Psi', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]',        'FontName', FNAME);
ylabel('\Psi [Wb]',         'FontName', FNAME);
xlim([0, t_final]);
%xlim([0, 1]);
% -------------------------------------------------------------------------
% CHART 2: Phase Portrait (Steady State — last 20%)
% -------------------------------------------------------------------------
fig2 = figure('Name', '2. Phase Portrait', 'Color', 'w', ...
              'Position', [120, 120, FIG_W, FIG_H]);
plot(x1_steady, x2_steady, 'Color', COL1, 'LineWidth', 0.8);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Phase Portrait', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME, 'Interpreter', 'tex');
ylabel('U_c [V]',   'FontName', FNAME, 'Interpreter', 'tex');

% -------------------------------------------------------------------------
% CHART 3: 3D Extended Phase Space Trajectory (Full Simulation)
% -------------------------------------------------------------------------
fig3 = figure('Name', '3. 3D Phase Space', 'Color', 'w', ...
              'Position', [140, 140, FIG_W, FIG_H]);
plot3(x1, x2, T, 'Color', COL1, 'LineWidth', 0.8);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Phase Space Trajectory', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]',   'FontName', FNAME);
ylabel('U_c [V]',     'FontName', FNAME);
zlabel('t [s]',  'FontName', FNAME);
view(45, 30);
rotate3d on;

% -------------------------------------------------------------------------
% CHART 4: FFT Spectrum (Steady State — last 20%)
% -------------------------------------------------------------------------
dt_fft    = mean(diff(T_steady));
t_fft     = T_steady(1) : dt_fft : T_steady(end);
x2_interp = interp1(T_steady, x2_steady, t_fft, 'spline');

L    = length(x2_interp);
Fs   = 1 / dt_fft;
Y_fft = fft(x2_interp);
P2   = abs(Y_fft / L);
P1   = P2(1:floor(L/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
f    = Fs * (0:(L/2)) / L;

fig4 = figure('Name', '4. FFT Spectrum', 'Color', 'w', ...
              'Position', [160, 160, FIG_W, FIG_H]);
plot(f, P1, 'Color', COL1, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('FFT Spectrum', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('f [Hz]',  'FontName', FNAME);
ylabel('|U_c| [V]',         'FontName', FNAME);
xlim([0, 350]);

% Mark fundamental frequency and its harmonics
f0 = omg / (2*pi);  % 50 Hz
hold on;
for k = 1:7
    xline(k*f0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
hold off;

% -------------------------------------------------------------------------
% CHART 5: Stroboscopic Diagram / Poincare Section (Steady State)
% -------------------------------------------------------------------------
% Sample at t = k*T0 using deval (the solver's own consistent interpolant)
% instead of interp1. interp1 fits an external spline through the solver's
% irregular output grid, and its small fitting error varies from cycle to
% cycle — so a true period-1 orbit appears as a tiny cluster of points
% rather than a single point. deval evaluates the same continuous extension
% the solver used during integration, removing that artifact.
T0      = 2*pi / omg;
t_strob = 0 : T0 : T(end);
t_strob = t_strob(t_strob >= t_cutoff);  % Keep only steady-state samples

XY       = deval(sol, t_strob);   % rows: [Psi; Uc; phase], cols: sample times
x1_strob = XY(1, :).';
x2_strob = XY(2, :).';

fig5 = figure('Name', '5. Stroboscopic Diagram', 'Color', 'w', ...
              'Position', [180, 180, FIG_W, FIG_H]);
plot(x1_strob, x2_strob, '.k', 'MarkerSize', 12);
%ylim([1.5e05, 1.6e05])
%xlim([134.29, 134.31])
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Stroboscopic Diagram', ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME);
%xlim([34, 35]);
%ylim([1.5e05, 1.6e05]);
ylabel('U_c [V]',   'FontName', FNAME);

% -------------------------------------------------------------------------
% CHART 6: Lyapunov Exponents (Full Simulation)
% -------------------------------------------------------------------------
fig6 = figure('Name', '6. Lyapunov Exponents', 'Color', 'w', ...
              'Position', [200, 200, FIG_W, FIG_H]);
plot(T_lyap, L_exp(:,1), 'Color', 'r', 'LineWidth', LW); hold on;
plot(T_lyap, L_exp(:,2), 'Color', COL3, 'LineWidth', LW);
plot(T_lyap, L_exp(:,3), 'Color', COL1, 'LineWidth', LW);
yline(0, '--k', 'LineWidth', 1.0, 'Label', '\lambda = 0', ...
      'FontName', FNAME, 'FontSize', FSIZE-1);
hold off;
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Lyapunov Exponents', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]',   'FontName', FNAME);
ylabel('\lambda_i',    'FontName', FNAME);
legend({'\lambda_1', '\lambda_2', '\lambda_3'}, ...
       'FontName', FNAME, 'FontSize', FSIZE, 'Location', 'best');

disp('===================================================');
disp('Analysis complete.');
disp('===================================================');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function dydt = ferro_A_ode(t, y)
% Right-hand side of the ferroresonant circuit ODE (Circuit A).
% States: y(1) = Psi [Wb], y(2) = Uc [V], y(3) = phase [rad]
    global par
    dydt    = zeros(3,1);
    dydt(1) = par(1)*cos(y(3)+par(9)) - par(2)*y(2) - par(3)*y(1) - par(4)*y(1)^par(11);
    dydt(2) = par(5)*cos(y(3)+par(9)) - par(6)*y(2) + par(7)*y(1) + par(8)*y(1)^par(11);
    dydt(3) = par(10);
end

function f = ferro_A_ext(t, X)
% Extended ODE system for Lyapunov exponent calculation.
% X(1:3)   - state vector [Psi; Uc; phase]
% X(4:12)  - columns of the variational matrix Y (3x3), stored column-major
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

    % Jacobian of the vector field
    Jac = [ -b3 - par(11)*b4*x1^(par(11)-1),  -b2,  -b1*sin(x3+alf);
             p3 + par(11)*p4*x1^(par(11)-1),  -p2,  -p1*sin(x3+alf);
             0,                                  0,   0              ];

    dY       = Jac * Y;
    f(4:12)  = [dY(1,1); dY(2,1); dY(3,1); ...
                dY(1,2); dY(2,2); dY(3,2); ...
                dY(1,3); dY(2,3); dY(3,3)];
end
