% =========================================================================
% FERRO_F: BIFURCATION DIAGRAM (VOLTAGE SWEEP)
% =========================================================================
% Generates a bifurcation diagram for Circuit B (4th-order ferroresonant
% circuit) by sweeping the source voltage amplitude U and sampling the
% capacitor voltage Uc at each period of the fundamental frequency
% (Poincare section / stroboscopic map).
%
% Method: continuation — the final state of each simulation step is used
%         as the initial condition for the next step, which accelerates
%         convergence to steady state and preserves solution branches.
%
% Sampling note: the capacitor voltage is sampled at exactly k*T0 by
%         placing the sample instants directly in the solver's output
%         vector (tspan). ode15s then returns the state at those instants
%         using its own internally consistent interpolant. This avoids the
%         cycle-to-cycle scatter produced by an external interp1 spline,
%         so a true period-1 orbit collapses to a single point.
% =========================================================================
clear; clc; close all;
global par;

disp('Starting Bifurcation Diagram - Ferro_F (Voltage Sweep)...');

% -------------------------------------------------------------------------
% 1. FIXED SYSTEM PARAMETERS
% -------------------------------------------------------------------------
alf = 0.0;          % Initial phase angle [rad]
Rs  = 0.2;          % Series resistance [Ohm]
RN  = 1.25;         % Neutral-point grounding resistance [Ohm]
Cs  = 3.1522E-6;    % Series capacitance [F]
Cr  = 4.6607E-6;    % Parallel capacitance [F]
q   = 11;           % Degree of nonlinearity
m1  = 0.04320728;
mq  = 1.444531e-18;
h1  = 3.396059e-05;
h3  = -6.906720e-14;
h5  = 1.593581e-22;
p13 = 1.2 * 2.122891e+04;
Ro  = 1 / (h1 + p13*p13*(h3 + h5*p13*p13));
p14 = 1.0 / (Ro*Cr);
omg = 100*pi;       % Angular frequency [rad/s]  (50 Hz)
T0  = 2*pi / omg;   % Fundamental period [s]    (0.02 s)

% -------------------------------------------------------------------------
% 2. BIFURCATION SWEEP SETTINGS
% -------------------------------------------------------------------------
U_start   = 30000;      % Start voltage [V]
U_end     = 35000;      % End voltage [V]
num_steps = 250;        % Number of voltage steps
U_vector  = linspace(U_start, U_end, num_steps);

% Simulation length per step (cycle-based for consistency)
num_cycles     = 50;    % Total cycles to simulate per step (~1.0 s)
discard_cycles = 35;    % Transient cycles to discard before sampling

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
y0      = [80; 25300; 16000; 0];   % [Psi; Uc; UL; phase]
%y0      = [200 ; 0; 0; 0];
% Sample instants k*T0 (held constant across the sweep). These are placed
% directly in tspan so the solver reports Uc at exactly these times.
sample_times = (discard_cycles : num_cycles) * T0;

% -------------------------------------------------------------------------
% 3. FIGURE SETUP  (thesis formatting)
% -------------------------------------------------------------------------
FNAME  = 'Times New Roman';
FSIZE  = 12;
FTITLE = 14;

fig = figure('Name', 'Bifurcation Diagram: Voltage Sweep (Ferro_F)', ...
             'Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
xlabel('Source Voltage Amplitude U [V]',    'FontName', FNAME);
ylabel('Sampled Capacitor Voltage U_c [V]', 'FontName', FNAME);
title('Bifurcation Diagram: Source Voltage Sweep', ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlim([U_start, U_end]);

% -------------------------------------------------------------------------
% 4. MAIN CONTINUATION LOOP
% -------------------------------------------------------------------------
tic;
for i = 1:num_steps

    U = U_vector(i);

    % Recalculate parameters dependent on U
    Rz = Rs + RN;
    p2 = 1/(Cs*Rz);   p1 = U*p2;
    b2 = 1/(Rz*Cr);   b1 = U*b2;
    b3 = h1/Cr;       b4 = h3/Cr;     b5 = h5/Cr;
    b6 = m1/Cr;       b7 = mq/Cr;

    %   Index: 1  2  3  4  5  6  7  8  9  10   11   12  13   14
    par = [p1; p2; b1; b2; b3; b4; b5; b6; b7; alf; omg; q; p13; p14];

    % Integrate. The sample instants are placed directly in tspan, so the
    % returned Y holds the state at exactly [0, k*T0] — no interpolation.
    [T, Y] = ode15s(@ferro_F_ode, [0, sample_times], y0, options);

    % Continuation: pass final state as next initial condition
    y0    = Y(end, :)';
    y0(4) = mod(y0(4), 2*pi);   % Wrap phase to [0, 2*pi]

    % Sampled Uc values (skip the prepended t = 0 row)
    sampled_points = Y(2:end, 2);

    % Plot sampled Uc values as black dots
    plot(repmat(U, length(sampled_points), 1), sampled_points, ...
         '.k', 'MarkerSize', 2);
    drawnow;

    % Progress indicator every 50 steps
    if mod(i, 50) == 0
        fprintf('Progress: %d / %d steps (U = %.0f V)\n', i, num_steps, U);
    end

end
elapsed = toc;
fprintf('Bifurcation diagram completed in %.1f s.\n', elapsed);

disp('Done.');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function dydt = ferro_F_ode(t, y)
% Right-hand side of the ferroresonant circuit ODE (Circuit B, 4th order).
% States: y(1) = Psi [Wb], y(2) = Uc [V], y(3) = UL [V], y(4) = phase [rad]
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
