% =========================================================================
% FERRO_F: MAIN ANALYSIS
% =========================================================================
% Circuit B (4th-order ferroresonant circuit) - Ferroresonance Analysis
% Computes: time-domain waveforms, phase portraits, 3D trajectory,
%           FFT spectrum, stroboscopic diagram, Lyapunov exponents.
%
% States: y(1) = Psi  [Wb]  - magnetic flux linkage
%         y(2) = Uc   [V]   - capacitor voltage
%         y(3) = UL   [V]   - inductor voltage
%         y(4) = phi  [rad] - excitation phase
%
% Reference for Lyapunov algorithm:
%   A. Wolf, J. B. Swift, H. L. Swinney, and J. A. Vastano,
%   "Determining Lyapunov Exponents from a Time Series,"
%   Physica D, Vol. 16, pp. 285-317, 1985.
% =========================================================================
clear; clc; close all;
global par;

disp('Starting Ferro_F Analysis...');

% -------------------------------------------------------------------------
% 1. SYSTEM PARAMETERS
% -------------------------------------------------------------------------
U = 35000;
%U   = 27500;        % Source voltage amplitude [V]  (chaotic operating point)
alf = 0.0;          % Initial phase angle [rad]
Rs  = 0.2;          % Series resistance [Ohm]
%Rs  = 5.0; 
RN  = 1.25;         % Neutral-point grounding resistance [Ohm]
Cs  = 3.1522E-6;    % Series capacitance [F]
%Cs  = 10E-6;
Cr  = 4.6607E-6;    % Parallel capacitance [F]
%Cr  = 40E-6;
q   = 11;           % Degree of nonlinearity in flux-current characteristic
m1  = 0.04320728;   % Linear coefficient of magnetization curve
mq  = 1.444531e-18; % Nonlinear coefficient of magnetization curve

% Core-loss polynomial coefficients (i_R = h1*UL + h3*UL^3 + h5*UL^5)
h1  = 3.396059e-05;
h3  = -6.906720e-14;
h5  = 1.593581e-22;

% Saturation switch for the core-loss model
p13 = 1.2 * 2.122891e+04;                          % Switching voltage [V]
Ro  = 1 / (h1 + p13*p13*(h3 + h5*p13*p13));        % Linear conductance above p13
p14 = 1.0 / (Ro*Cr);

omg = 100*pi;       % Angular frequency [rad/s]  (50 Hz)

% Derived parameters
Rz = Rs + RN;
p2 = 1/(Cs*Rz);     p1 = U*p2;
b2 = 1/(Rz*Cr);     b1 = U*b2;
b3 = h1/Cr;         b4 = h3/Cr;       b5 = h5/Cr;
b6 = m1/Cr;         b7 = mq/Cr;

%   Index: 1  2  3  4  5  6  7  8  9  10   11   12  13   14
par = [p1; p2; b1; b2; b3; b4; b5; b6; b7; alf; omg; q; p13; p14];

% -------------------------------------------------------------------------
% 2. SIMULATION SETTINGS
% -------------------------------------------------------------------------
t_final = 2.0;     % Total simulation time [s]
tspan   = [0, t_final];
y0      = [80; 25300; 16000; 0];   % [Psi; Uc; UL; phase]
%y0      = [0 ; 0; 0; 0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% -------------------------------------------------------------------------
% 3. RUN MAIN SIMULATION
% -------------------------------------------------------------------------
disp('Running ODE simulation...');
tic;
% Capture the solution structure so the stroboscopic section can use
% deval (the solver's own consistent interpolant).
sol = ode15s(@ferro_F_ode, tspan, y0, options);
elapsed = toc;
fprintf('ODE simulation completed in %.2f s.\n', elapsed);

% Reconstruct T and Y from sol; remove any duplicate time points that
% ode15s may emit at saturation switches.
T = sol.x(:);
Y = sol.y.';
[T, unique_idx] = unique(T);
Y = Y(unique_idx, :); %for removing duplicate points...

x1_flux = Y(:,1);   % Magnetic flux linkage Psi [Wb]
x2_uc   = Y(:,2);   % Capacitor voltage Uc      [V]
x3_ul   = Y(:,3);   % Inductor voltage UL       [V]

% Extract last 20% of simulation as steady state
t_cutoff   = 0.8 * T(end);
idx_steady = T >= t_cutoff;
T_steady   = T(idx_steady);
x1_steady  = x1_flux(idx_steady);  %only these elements which have true value
x2_steady  = x2_uc(idx_steady);
x3_steady  = x3_ul(idx_steady);

% -------------------------------------------------------------------------
% 4. OVERVOLTAGE REPORT
% -------------------------------------------------------------------------
max_transient_x2 = max(abs(x2_uc));
max_steady_x2    = max(abs(x2_steady));
max_transient_x2pu = max(abs(x2_uc))/U;
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
disp('Calculating Lyapunov exponents (4D system, this may take a minute)...');
[T_lyap, L_exp] = lyapunov(4, @ferro_F_ext, @ode15s, 0, 0.005, t_final, y0, 0);

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
COL1   = '#0072BD';   % Blue   (capacitor voltage)
COL2   = '#D95319';   % Orange (flux)
COL3   = '#77AC30';   % Green  (auxiliary)
COL4   = '#7E2F8E';   % Purple (inductor voltage)

% -------------------------------------------------------------------------
% CHART 1: Time-Domain Waveforms (Full Simulation)
% -------------------------------------------------------------------------
fig1 = figure('Name', '1. Time Domain', 'Color', 'w', ...
              'Position', [100, 50, FIG_W, 700]);

subplot(3,1,1);
plot(T, x2_uc, 'Color', COL1, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Capacitor Voltage U_c', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('U_c [V]',    'FontName', FNAME);
%xlim([0, t_final]);
xlim([0, 1]);

subplot(3,1,2);
plot(T, x3_ul, 'Color', COL4, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Inductor Voltage U_L', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('U_L [V]',    'FontName', FNAME);
%xlim([0, t_final]);
xlim([0, 1]);

subplot(3,1,3);
plot(T, x1_flux, 'Color', COL2, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Flux Linkage \Psi', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('\Psi [Wb]',  'FontName', FNAME);
%xlim([0, t_final]);
xlim([0, 1]);

% -------------------------------------------------------------------------
% CHART 2a: Phase Portrait Psi vs Uc (Steady State — last 20%)
% -------------------------------------------------------------------------
fig2a = figure('Name', '2a. Phase Portrait (Psi vs Uc)', 'Color', 'w', ...
               'Position', [120, 120, FIG_W, FIG_H]);
plot(x1_steady, x2_steady, 'Color', COL1, 'LineWidth', 0.8);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Phase Portrait', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME, 'Interpreter', 'tex');
ylabel('U_c [V]',   'FontName', FNAME, 'Interpreter', 'tex');

% -------------------------------------------------------------------------
% CHART 2b: Phase Portrait Psi vs UL (Steady State — last 20%)
% -------------------------------------------------------------------------
fig2b = figure('Name', '2b. Phase Portrait (Psi vs UL)', 'Color', 'w', ...
               'Position', [140, 140, FIG_W, FIG_H]);
plot(x1_steady, x3_steady, 'Color', COL1, 'LineWidth', 0.8);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Phase Portrait', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME, 'Interpreter', 'tex');
ylabel('U_L [V]',   'FontName', FNAME, 'Interpreter', 'tex');

% -------------------------------------------------------------------------
% CHART 3: 3D Extended Phase Space Trajectory (Full Simulation)
% -------------------------------------------------------------------------
fig3 = figure('Name', '3. 3D Phase Space', 'Color', 'w', ...
              'Position', [160, 160, FIG_W, FIG_H]);
plot3(x1_flux, x2_uc, T, 'Color', COL1, 'LineWidth', 0.8);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('3D Extended Phase Space Trajectory', 'FontName', FNAME, 'FontSize', FTITLE);
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
              'Position', [180, 180, FIG_W, FIG_H]);
plot(f, P1, 'Color', COL1, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('FFT Spectrum', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('f [Hz]', 'FontName', FNAME);
ylabel('|U_c| [V]',        'FontName', FNAME);
xlim([0, 500]);

% Mark fundamental frequency and its harmonics
f0 = omg / (2*pi);  % 50 Hz
hold on;
for k = 1:10
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

XY       = deval(sol, t_strob);   % rows: [Psi; Uc; UL; phase], cols: sample times
x1_strob = XY(1, :).';
x2_strob = XY(2, :).';

fig5 = figure('Name', '5. Stroboscopic Diagram', 'Color', 'w', ...
              'Position', [200, 200, FIG_W, FIG_H]);
plot(x1_strob, x2_strob, '.k', 'MarkerSize', 12);
%ylim([4e04,5e04])
%xlim([0,3.5])
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
              'Position', [220, 220, FIG_W, FIG_H]);
plot(T_lyap, L_exp(:,1), 'Color', 'r',  'LineWidth', LW); hold on;
plot(T_lyap, L_exp(:,2), 'Color', COL3, 'LineWidth', LW);
plot(T_lyap, L_exp(:,3), 'Color', COL1, 'LineWidth', LW);
plot(T_lyap, L_exp(:,4), 'Color', 'm',  'LineWidth', LW);
yline(0, '--k', 'LineWidth', 1.0, 'Label', '\lambda = 0', ...
      'FontName', FNAME, 'FontSize', FSIZE-1);
hold off;
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Lyapunov Exponents', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('\lambda_i',  'FontName', FNAME);
legend({'\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4'}, ...
       'FontName', FNAME, 'FontSize', FSIZE, 'Location', 'best');

disp('===================================================');
disp('Analysis complete.');
disp('===================================================');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function dydt = ferro_F_ode(t, y)
% Right-hand side of the ferroresonant circuit ODE (Circuit B, 4th order).
% States: y(1) = Psi [Wb], y(2) = Uc [V], y(3) = UL [V], y(4) = phase [rad]
%
% Core-loss current i_R is a 5th-order polynomial of UL inside the
% saturation band |UL| < p13, and switches to a linear branch outside it.
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
% Extended ODE system for Lyapunov exponent calculation.
% X(1:4)   - state vector [Psi; Uc; UL; phase]
% X(5:20)  - columns of the variational matrix Y (4x4), stored column-major
    global par
    p1=par(1);  p2=par(2);  b1=par(3);  b2=par(4);
    b3=par(5);  b4=par(6);  b5=par(7);  b6=par(8);  b7=par(9);
    alf=par(10); omg=par(11); q=par(12); p13=par(13); p14=par(14);

    x1 = X(1); x2 = X(2); x3 = X(3); x4 = X(4);
    n  = 4;
    Y  = reshape(X(n+1:end), n, n);

    % Piecewise core-loss current and its derivative wrt UL
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

    % Jacobian of the vector field
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
