% =========================================================================
% DUFFING FERRORESONANCE: MAIN ANALYSIS  (Circuit C)
% =========================================================================
% Circuit C - Duffing-type ferroresonant circuit driven by a current source
% with a sinusoidal component and a DC offset. The DC component I0 is
% responsible for the appearance of a subharmonic at f0/2.
%
% Computes: time-domain waveforms, phase portrait, 3D trajectory,
%           FFT spectrum, stroboscopic diagram, Lyapunov exponents.
%
% States: y(1) = Psi  [Wb]  - magnetic flux linkage
%         y(2) = Uc   [V]   - capacitor voltage
%         y(3) = phi  [rad] - excitation phase
%
% Reference for Lyapunov algorithm:
%   A. Wolf, J. B. Swift, H. L. Swinney, and J. A. Vastano,
%   "Determining Lyapunov Exponents from a Time Series,"
%   Physica D, Vol. 16, pp. 285-317, 1985.
% =========================================================================
clear; clc; close all;
global par;

disp('Starting Duffing Ferroresonance Analysis...');

% -------------------------------------------------------------------------
% 1. SYSTEM PARAMETERS
% -------------------------------------------------------------------------
Im  = 1.5;          % Sinusoidal current source amplitude [A]
I0  = 1.75;          % DC component of the current source [A]
C   = 40.53E-6;     % Parallel capacitance [F]
%C   = 100E-6;
R   = 19.38E+6;     % Core-loss (parallel) resistance [Ohm]
%R   = 30E+6;
m1  = 0.32;         % Linear coefficient of magnetization curve
m3  = 0.9;          % Cubic coefficient of magnetization curve
omg = 100*pi;       % Angular frequency [rad/s]  (50 Hz)

% The degree of nonlinearity q = 3 is hardcoded as y(1)^3 in the ODE
% (cubic Duffing-type nonlinearity).

%   Index: 1   2   3   4    5    6    7
par = [Im; I0; C; R; m1; m3; omg];

% -------------------------------------------------------------------------
% 2. SIMULATION SETTINGS
% -------------------------------------------------------------------------
t_final = 2.0;     % Total simulation time [s]
tspan   = [0, t_final];
y0      = [0; 0; 0];   % Initial conditions: [Psi; Uc; phase]

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% -------------------------------------------------------------------------
% 3. RUN MAIN SIMULATION
% -------------------------------------------------------------------------
disp('Running ODE simulation...');
tic;
% Capture the solution structure so the stroboscopic section can use
% deval (the solver's own consistent interpolant).
sol = ode45(@duffing_ode, tspan, y0, options);
elapsed = toc;
fprintf('ODE simulation completed in %.2f s.\n', elapsed);

T = sol.x(:);
Y = sol.y.';

x1_psi = Y(:,1);   % Magnetic flux linkage Psi [Wb]
x2_uc  = Y(:,2);   % Capacitor voltage Uc      [V]

% Extract last 20% of simulation as steady state
t_cutoff   = 0.8 * T(end);
idx_steady = T >= t_cutoff;
T_steady   = T(idx_steady);
x1_steady  = x1_psi(idx_steady);
x2_steady  = x2_uc(idx_steady);

% -------------------------------------------------------------------------
% 4. OVERVOLTAGE REPORT
% -------------------------------------------------------------------------
max_transient_x2 = max(abs(x2_uc));
max_steady_x2    = max(abs(x2_steady));

disp('===================================================');
disp('   Capacitor Overvoltages (x2 = Uc)');
disp('===================================================');
fprintf('> Max overvoltage (full transient):  %.2f V\n', max_transient_x2);
fprintf('> Max overvoltage (steady state):    %.2f V\n', max_steady_x2);
disp('===================================================');

% -------------------------------------------------------------------------
% 5. LYAPUNOV EXPONENTS
% -------------------------------------------------------------------------
disp('Calculating Lyapunov exponents...');
[T_lyap, L_exp] = lyapunov(3, @duffing_ext, @ode45, 0, 0.005, t_final, y0, 0);

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
FSIZE  = 12;    % Axis tick/label font size
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
plot(T, x2_uc, 'Color', COL1, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Capacitor Voltage U_c', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('U_c [V]',    'FontName', FNAME);
xlim([0, t_final]);

subplot(2,1,2);
plot(T, x1_psi, 'Color', COL2, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Flux Linkage \Psi', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('\Psi [Wb]',  'FontName', FNAME);
xlim([0, t_final]);

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
plot3(x1_psi, x2_uc, T, 'Color', COL1, 'LineWidth', 0.8);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Phase Space Trajectory', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]',  'FontName', FNAME);
ylabel('U_c [V]',    'FontName', FNAME);
zlabel('Time t [s]', 'FontName', FNAME);
view(45, 30);
rotate3d on;

% -------------------------------------------------------------------------
% CHART 4: FFT Spectrum (Steady State — last 20%)
% -------------------------------------------------------------------------
dt_fft    = mean(diff(T_steady));
t_fft     = T_steady(1) : dt_fft : T_steady(end);
x2_interp = interp1(T_steady, x2_steady, t_fft, 'spline');

L     = length(x2_interp);
Fs    = 1 / dt_fft;
Y_fft = fft(x2_interp);
P2    = abs(Y_fft / L);
P1    = P2(1:floor(L/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
f     = Fs * (0:(L/2)) / L;

fig4 = figure('Name', '4. FFT Spectrum', 'Color', 'w', ...
              'Position', [160, 160, FIG_W, FIG_H]);
plot(f, P1, 'Color', COL1, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('FFT Spectrum', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('f [Hz]', 'FontName', FNAME);
ylabel('|U_c| [V]',        'FontName', FNAME);
xlim([0, 350]);

% Mark fundamental frequency and its harmonics; mark f0/2 separately as
% the subharmonic produced by the DC offset I0 is the diagnostic feature.
f0 = omg / (2*pi);  % 50 Hz
hold on;
%xline(f0/2, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0, ...
%      'Label', 'f_0/2', 'FontName', FNAME, 'FontSize', FSIZE-1, ...
%      'Interpreter', 'tex');
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
% rather than a single point. deval removes that artifact.
T0      = 2*pi / omg;
t_strob = 0 : T0 : T(end);
t_strob = t_strob(t_strob >= t_cutoff);  % Keep only steady-state samples

XY       = deval(sol, t_strob);   % rows: [Psi; Uc; phase], cols: sample times
x1_strob = XY(1, :).';
x2_strob = XY(2, :).';

fig5 = figure('Name', '5. Stroboscopic Diagram', 'Color', 'w', ...
              'Position', [180, 180, FIG_W, FIG_H]);
plot(x1_strob, x2_strob, '.k', 'MarkerSize', 12);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Stroboscopic Diagram', ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME);
ylabel('U_c [V]',   'FontName', FNAME);

% -------------------------------------------------------------------------
% CHART 6: Lyapunov Exponents (Full Simulation)
% -------------------------------------------------------------------------
fig6 = figure('Name', '6. Lyapunov Exponents', 'Color', 'w', ...
              'Position', [200, 200, FIG_W, FIG_H]);
plot(T_lyap, L_exp(:,1), 'Color', 'r',  'LineWidth', LW); hold on;
plot(T_lyap, L_exp(:,2), 'Color', COL3, 'LineWidth', LW);
plot(T_lyap, L_exp(:,3), 'Color', COL1, 'LineWidth', LW);
yline(0, '--k', 'LineWidth', 1.0, 'Label', '\lambda = 0', ...
      'FontName', FNAME, 'FontSize', FSIZE-1);
hold off;
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Lyapunov Exponents', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('Time t [s]', 'FontName', FNAME);
ylabel('\lambda_i',  'FontName', FNAME);
legend({'\lambda_1', '\lambda_2', '\lambda_3'}, ...
       'FontName', FNAME, 'FontSize', FSIZE, 'Location', 'best');

disp('===================================================');
disp('Analysis complete.');
disp('===================================================');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function dydt = duffing_ode(t, y)
% Right-hand side of the Duffing-type ferroresonant circuit ODE (Circuit C).
% States: y(1) = Psi [Wb], y(2) = Uc [V], y(3) = phase [rad]
    global par
    Im=par(1); I0=par(2); C=par(3); R=par(4);
    m1=par(5); m3=par(6); omg=par(7);

    dydt    = zeros(3,1);
    dydt(1) = y(2);
    dydt(2) = (1/C)*(Im*cos(y(3)) + I0) - (1/(R*C))*y(2) ...
              - (m3/C)*(y(1)^3) - (m1/C)*y(1);
    dydt(3) = omg;
end

function f = duffing_ext(t, X)
% Extended ODE system for Lyapunov exponent calculation.
% X(1:3)   - state vector [Psi; Uc; phase]
% X(4:12)  - columns of the variational matrix Y (3x3), stored column-major
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

    % Jacobian of the vector field
    Jac        = zeros(n,n);
    Jac(1,2)   = 1;
    Jac(2,1)   = -(3*m3/C)*x1^2 - (m1/C);
    Jac(2,2)   = -1/(R*C);
    Jac(2,3)   = -(Im/C)*sin(x3);

    dY         = Jac * Y;
    f(n+1:end) = dY(:);
end
