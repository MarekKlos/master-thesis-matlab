% =========================================================================
% FERRO_A WITH MOV ARRESTER  —  Abbasi et al. 2010 model
% =========================================================================
% Circuit A modified with a Metal Oxide Varistor (MOV) connected in
% parallel with the transformer magnetising branch.  The MOV is modelled
% following the formulation of:
%
%   A. Abbasi, M. Rostami, S. H. Fathi, H. R. Abbasi, H. Abdollahi,
%   "Effect of Metal Oxide Arrester on Chaotic Behavior of Power
%    Transformers," Energy and Power Engineering, Vol. 2, pp. 254-261,
%    2010.
%
% MOV characteristic, paper equation (10):
%       i_MOV = (|V_m| / k)^alpha * sign(V_m)
%
% with typical parameters alpha = 25 and k chosen relative to the
% transformer rated voltage.
% =========================================================================
clear; clc; close all;
global par;

disp('Starting Ferro_A with MOV Analysis (Abbasi et al. 2010 model)...');

% -------------------------------------------------------------------------
% 1. CIRCUIT PARAMETERS  (Circuit A base values — known chaotic point)
% -------------------------------------------------------------------------
U   = 90000;        % Source voltage amplitude [V]
alf = 0.0;
Rs  = 1.5;          % Series resistance [Ohm]
Rm  = 50000;        % Core-loss resistance [Ohm]
C   = 0.78E-6;      % Capacitance [F]
q   = 11;
m1  = 4.079845e-04;
mq  = 2.108275e-27;
omg = 100*pi;

% -------------------------------------------------------------------------
% Reset text interpreters to MATLAB's default 'tex' (in case a previous
% script set them to 'latex' globally, which would reject these labels).
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot, 'defaultTextInterpreter', 'tex');
set(groot, 'defaultLegendInterpreter', 'tex');

% 2. MOV PARAMETERS  (Abbasi et al. 2010, equation 10)
% -------------------------------------------------------------------------
MOV_ENABLED = true;     % Toggle MOV on/off
alpha = 15;             % Nonlinearity exponent (paper uses 25)

% k_MOV is the voltage at which the MOV conducts 1 A — its effective
% clamping level.  IMPORTANT: it must sit ABOVE the normal fundamental-mode
% voltage (~160 kV here) but BELOW the chaotic overvoltage (~284 kV here),
% so that the MOV stays off during normal operation and clamps the chaotic
% overvoltage.  If k is set too low, the MOV is driven far past its knee
% and the (V/k)^25 term explodes, causing the solver to fail.
%
% NOTE on the paper's value: Abbasi et al. give k = 2.5101 in PER-UNIT
% (their base voltage is 635.1 kV).  That value cannot be copied directly
% into this SI-unit model; what matters is the physical placement of k
% relative to the circuit's own voltages.
k_MOV = 2.8 * U;        % [V]  ~216 kV — between normal and chaotic levels

% -------------------------------------------------------------------------
% 3. DERIVED CIRCUIT COEFFICIENTS
% -------------------------------------------------------------------------
p2 = 1/(C*(Rm+Rs)); p1 = U*p2;
p3 = p2*m1*Rm;      p4 = p2*mq*Rm;
b2 = Rm/(Rm+Rs);    b1 = U*b2;
b3 = b2*m1*Rs;      b4 = b2*mq*Rs;

%   Index:  1  2  3  4  5  6  7  8   9   10  11  12          13     14
par = [b1; b2; b3; b4; p1; p2; p3; p4; alf; omg; q;
       double(MOV_ENABLED); k_MOV; alpha];

% -------------------------------------------------------------------------
% 4. SIMULATION SETTINGS
% -------------------------------------------------------------------------
t_final = 2.0;
tspan   = [0, t_final];
y0      = [280; 1530; 0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% -------------------------------------------------------------------------
% 5. RUN PAIRED SIMULATIONS  (with and without MOV)
% -------------------------------------------------------------------------
% The MOV term makes the system stiff (alpha = 25 produces very sharp
% nonlinearity), so a stiff solver (ode15s) is used instead of ode45.
disp('Running simulation WITHOUT MOV...');
par(12) = 0;
[T_off, Y_off] = ode15s(@ferro_mov_ode, tspan, y0, options);

disp('Running simulation WITH MOV...');
par(12) = 1;
[T_on, Y_on] = ode15s(@ferro_mov_ode, tspan, y0, options);

% Sanity check: confirm both solvers reached the end of tspan
if T_off(end) < 0.99*t_final
    warning('WITHOUT-MOV simulation stopped early at t = %.3f s.', T_off(end));
end
if T_on(end) < 0.99*t_final
    warning(['WITH-MOV simulation stopped early at t = %.3f s. ', ...
             'The MOV term is likely too stiff — try lowering alpha or ', ...
             'raising k_MOV.'], T_on(end));
end

Uc_off = Y_off(:,2);  Psi_off = Y_off(:,1);
Uc_on  = Y_on(:,2);   Psi_on  = Y_on(:,1);

ss_off = T_off >= 0.8*t_final;
ss_on  = T_on  >= 0.8*t_final;

% -------------------------------------------------------------------------
% 6. OVERVOLTAGE REPORT
% -------------------------------------------------------------------------
fprintf('\n===================================================\n');
fprintf('   OVERVOLTAGE COMPARISON  (U = %.0f kV)\n', U/1e3);
fprintf('   MOV parameters: alpha = %d, k = %.0f kV\n', alpha, k_MOV/1e3);
fprintf('===================================================\n');
fprintf('Without MOV:\n');
fprintf('   Max transient    U_c = %.2f kV\n', max(abs(Uc_off))/1e3);
fprintf('   Max steady-state U_c = %.2f kV\n', max(abs(Uc_off(ss_off)))/1e3);
fprintf('With MOV:\n');
fprintf('   Max transient    U_c = %.2f kV\n', max(abs(Uc_on))/1e3);
fprintf('   Max steady-state U_c = %.2f kV\n', max(abs(Uc_on(ss_on)))/1e3);
reduction = (1 - max(abs(Uc_on(ss_on)))/max(abs(Uc_off(ss_off)))) * 100;
fprintf('   ==> Steady-state overvoltage reduction: %.1f%%\n', reduction);
fprintf('===================================================\n');

% =========================================================================
% 7. FIGURE FORMATTING
% =========================================================================
%FIG_W  = 800;  FIG_H  = 500;
FIG_W  = 800;  FIG_H  = 400;
FSIZE  = 10;   FTITLE = 12;
FNAME  = 'Times New Roman';
LW     = 1.5;
COL_OFF = '#D95319';   % Orange — without MOV
COL_ON  = '#0072BD';   % Blue   — with MOV

% -------------------------------------------------------------------------
% CHART 1: Time-domain comparison
% -------------------------------------------------------------------------
figure('Name', '1. U_c Time-Domain Comparison', 'Color', 'w', ...
       'Position', [100, 100, FIG_W, FIG_H]);

subplot(2,1,1);
plot(T_off, Uc_off/1e3, 'Color', COL_OFF, 'LineWidth', LW); grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Capacitor Voltage U_c — Without MOV', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('U_c [kV]', 'FontName', FNAME);
xlim([0, t_final]);

subplot(2,1,2);
plot(T_on, Uc_on/1e3, 'Color', COL_ON, 'LineWidth', LW); grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title(sprintf('Capacitor Voltage U_c — With MOV (k = %.0f kV, \\alpha = %d)', ...
              k_MOV/1e3, alpha), 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('t [s]', 'FontName', FNAME);
ylabel('U_c [kV]', 'FontName', FNAME);
xlim([0, t_final]);

% -------------------------------------------------------------------------
% CHART 2: Phase portrait comparison
% -------------------------------------------------------------------------
figure('Name', '2. Phase Portrait Comparison', 'Color', 'w', ...
       'Position', [120, 120, FIG_W, FIG_H]);

subplot(1,2,1);
plot(Psi_off(ss_off), Uc_off(ss_off)/1e3, 'Color', COL_OFF, 'LineWidth', 0.6);
grid on; axis tight;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Without MOV', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME);
ylabel('U_c [kV]', 'FontName', FNAME);

subplot(1,2,2);
plot(Psi_on(ss_on), Uc_on(ss_on)/1e3, 'Color', COL_ON, 'LineWidth', 0.8);
grid on; axis tight;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('With MOV', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME);
ylabel('U_c [kV]', 'FontName', FNAME);

% -------------------------------------------------------------------------
% CHART 3: MOV i-v characteristic (Abbasi et al. 2010, eq. 10)
% -------------------------------------------------------------------------
figure('Name', '3. MOV Characteristic', 'Color', 'w', ...
       'Position', [140, 140, FIG_W, FIG_H]);
v_range = linspace(-2.5*k_MOV, 2.5*k_MOV, 1000);
i_curve = (abs(v_range)/k_MOV).^alpha .* sign(v_range);

subplot(1,2,1);
plot(v_range/1e3, i_curve, 'Color', COL_ON, 'LineWidth', LW); grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('MOV i-u char. (linear scale)', ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlabel('U_m [kV]', 'FontName', FNAME);
ylabel('i_{MOV} [A]', 'FontName', FNAME);
xline(k_MOV/1e3,  '--k', 'LineWidth', 1.0, 'Label', 'k');
xline(-k_MOV/1e3, '--k', 'LineWidth', 1.0);
ylim([-5, 5]);

subplot(1,2,2);
semilogy(v_range/1e3, abs(i_curve) + 1e-15, 'Color', COL_ON, 'LineWidth', LW);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('MOV i-u char. (log scale)', ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlabel('U_m [kV]', 'FontName', FNAME);
ylabel('|i_{MOV}| [A]', 'FontName', FNAME);
ylim([1e-10, 1e6]);

% -------------------------------------------------------------------------
% CHART 4: FFT comparison
% -------------------------------------------------------------------------
figure('Name', '4. FFT Comparison', 'Color', 'w', ...
       'Position', [160, 160, FIG_W, FIG_H]);

[f1, P_off] = compute_fft(T_off(ss_off), Uc_off(ss_off));
[f2, P_on]  = compute_fft(T_on(ss_on),  Uc_on(ss_on));

subplot(2,1,1);
plot(f1, P_off/1e3, 'Color', COL_OFF, 'LineWidth', LW); grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('FFT Spectrum — Without MOV', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('f [Hz]', 'FontName', FNAME);
ylabel('|U_c| [kV]', 'FontName', FNAME);
xlim([0, 350]);

subplot(2,1,2);
plot(f2, P_on/1e3, 'Color', COL_ON, 'LineWidth', LW); grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('FFT Spectrum — With MOV', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('f [Hz]', 'FontName', FNAME);
ylabel('|U_c| [kV]', 'FontName', FNAME);
xlim([0, 350]);

% -------------------------------------------------------------------------
% CHART 5: Stroboscopic diagram comparison
% -------------------------------------------------------------------------
figure('Name', '5. Stroboscopic Comparison', 'Color', 'w', ...
       'Position', [180, 180, FIG_W, FIG_H]);

T0 = 2*pi/omg;
[Psi_s_off, Uc_s_off] = stroboscopic_sample(T_off, Psi_off, Uc_off, T0, 0.8*t_final);
[Psi_s_on,  Uc_s_on ] = stroboscopic_sample(T_on,  Psi_on,  Uc_on,  T0, 0.8*t_final);

subplot(1,2,1);
plot(Psi_s_off, Uc_s_off/1e3, '.', 'Color', COL_OFF, 'MarkerSize', 10);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Stroboscopic — Without MOV', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME);
ylabel('U_c [kV]', 'FontName', FNAME);

subplot(1,2,2);
plot(Psi_s_on, Uc_s_on/1e3, '.', 'Color', COL_ON, 'MarkerSize', 14);
grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
title('Stroboscopic — With MOV', 'FontName', FNAME, 'FontSize', FTITLE);
xlabel('\Psi [Wb]', 'FontName', FNAME);
ylabel('U_c [kV]', 'FontName', FNAME);
ylim([150, 160]);
xlim([30, 40]);
disp('Analysis complete.');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function dydt = ferro_mov_ode(t, y)
% Ferroresonant ODE with MOV in parallel with the transformer.
% MOV current uses Abbasi et al. 2010, equation (10):
%       i_MOV = (|V_m| / k)^alpha * sign(V_m)
%
% Including the MOV adds an extra current branch at the transformer node,
% modifying the dU_c/dt equation accordingly.
    global par

    b1  = par(1);  b2  = par(2);  b3  = par(3);  b4  = par(4);
    p1  = par(5);  p2  = par(6);  p3  = par(7);  p4  = par(8);
    alf = par(9);  omg = par(10); q   = par(11);
    use_mov = par(12);
    k_MOV   = par(13);
    alpha   = par(14);

    Psi = y(1);  Uc = y(2);  phi = y(3);

    % Reconstruct transformer (magnetising-branch) voltage u_m
    % from the algebraic relation u_m = u - U_c - R_s*i  where i is the
    % series current through the capacitor and resistor.
    u_source = (b1/b2) * cos(phi + alf);
    u_m      = u_source - Uc - (b3/b2)*Psi - (b4/b2)*Psi^q;

    % MOV current  — Abbasi et al. 2010, equation (10).
    % i_mov is a physical current [A] through the arrester, which is
    % connected in parallel with the capacitor branch.  At the node, this
    % current diverts away from the capacitor, so it enters dU_c/dt
    % divided by C (the same way any node current affects capacitor voltage).
    if use_mov > 0.5
        % Guard against u_m == 0 to avoid 0^alpha numerical issues
        if abs(u_m) < 1e-12
            i_mov = 0;
        else
            i_mov = (abs(u_m)/k_MOV)^alpha * sign(u_m);
        end
        % Physical current limit: a real arrester has finite energy-handling
        % capacity and cannot pass unbounded current.  Limiting i_mov keeps
        % the model physically realistic and prevents the (V/k)^25 term from
        % causing solver failure during the initial transient.
        I_MAX = 500;   % [A] — maximum arrester current
        i_mov = max(min(i_mov, I_MAX), -I_MAX);
    else
        i_mov = 0;
    end

    % Capacitance value (must match the C used to build the par vector)
    C_val    = 0.78e-6;
    mov_term = i_mov / C_val;

    dydt    = zeros(3,1);
    dydt(1) = b1*cos(phi+alf) - b2*Uc - b3*Psi - b4*Psi^q;
    dydt(2) = p1*cos(phi+alf) - p2*Uc + p3*Psi + p4*Psi^q + mov_term;
    dydt(3) = omg;
end

function [f, P1] = compute_fft(T, signal)
    % Guard against empty or too-short input (e.g. if the ODE solver
    % stopped early and the steady-state mask selected no samples).
    if numel(T) < 4
        warning('compute_fft received %d sample(s) — returning empty spectrum. The simulation likely stopped early.', numel(T));
        f = []; P1 = [];
        return;
    end
    dt = mean(diff(T));
    t_unif = T(1):dt:T(end);
    sig_interp = interp1(T, signal, t_unif, 'spline');
    L  = length(sig_interp);
    Fs = 1/dt;
    Y_fft = fft(sig_interp);
    P2 = abs(Y_fft / L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    f = Fs * (0:(L/2)) / L;
end

function [Psi_s, Uc_s] = stroboscopic_sample(T, Psi, Uc, T0, t_start)
    t_sample = t_start : T0 : T(end);
    Psi_s = interp1(T, Psi, t_sample, 'spline');
    Uc_s  = interp1(T, Uc,  t_sample, 'spline');
end
