% =========================================================================
% DUFFING FERRORESONANCE: BIFURCATION DIAGRAM (DAMPING RESISTANCE SWEEP)
% =========================================================================
% Generates a bifurcation diagram for Circuit C by sweeping the parallel
% (core-loss) resistance R at a fixed chaotic forcing amplitude Im, and
% sampling the capacitor voltage Uc at each period of the fundamental
% frequency (Poincare section / stroboscopic map).
%
% Note on damping: in this parallel topology a LARGER R means LESS
%       damping (less dissipation through the resistor). Chaos is
%       therefore expected at high R, regular response at low R.
%
% Method: continuation. Sampling done with the [0, sample_times] trick
%         (see Duffing_Bifurcation.m for the explanation).
% =========================================================================
clear; clc; close all;

disp('Starting Bifurcation Diagram - Duffing (Damping Resistance Sweep)...');

% -------------------------------------------------------------------------
% 1. FIXED SYSTEM PARAMETERS
% -------------------------------------------------------------------------
Im  = 1.5;          % Sinusoidal current source amplitude [A]  (chaotic point)
I0  = 0.5;          % DC component of the current source [A]
C   = 40.53E-6;     % Parallel capacitance [F]
m1  = 0.32;         % Linear coefficient of magnetization curve
m3  = 0.9;          % Cubic coefficient of magnetization curve
omg = 100*pi;       % Angular frequency [rad/s]  (50 Hz)
T0  = 2*pi / omg;   % Fundamental period [s]    (0.02 s)

% -------------------------------------------------------------------------
% 2. BIFURCATION SWEEP SETTINGS
% -------------------------------------------------------------------------
R_start   = 1E+6;       % Start resistance [Ohm]   (heavy damping)
R_end     = 40E+6;      % End resistance [Ohm]     (light damping → chaos)
num_steps = 400;        % Number of sweep points
R_vector  = linspace(R_start, R_end, num_steps);

% Simulation length per step (cycle-based for consistency)
num_cycles     = 40;    % Total cycles to simulate per step
discard_cycles = 25;    % Transient cycles to discard before sampling

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
y0      = [0; 0; 0];    % Initial condition [Psi; Uc; phase]

% Sample instants k*T0 (held constant across the sweep).
sample_times = (discard_cycles : num_cycles) * T0;

% -------------------------------------------------------------------------
% 3. FIGURE SETUP  (thesis formatting)
% -------------------------------------------------------------------------
FNAME  = 'Times New Roman';
FSIZE  = 12;
FTITLE = 14;

fig = figure('Name', 'Bifurcation Diagram: Damping Sweep (Duffing)', ...
             'Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
xlabel('Damping Resistance R [\Omega]',     'FontName', FNAME);
ylabel('Sampled Capacitor Voltage U_c [V]', 'FontName', FNAME);
title(sprintf('Bifurcation Diagram: Damping Sweep (I_m = %.2f A)', Im), ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlim([R_start, R_end]);

% -------------------------------------------------------------------------
% 4. MAIN CONTINUATION LOOP
% -------------------------------------------------------------------------
tic;
for i = 1:num_steps

    R = R_vector(i);

    duffing_ode = @(t, y) [
        y(2);
        (1/C)*(Im*cos(y(3)) + I0) - (1/(R*C))*y(2) ...
              - (m3/C)*(y(1)^3) - (m1/C)*y(1);
        omg
    ];

    % Integrate. Sample instants placed directly in tspan — no interp1.
    [T, Y] = ode45(duffing_ode, [0, sample_times], y0, options);

    % Continuation: pass final state as next initial condition
    y0    = Y(end, :)';
    y0(3) = mod(y0(3), 2*pi);

    % Sampled Uc values (skip the prepended t = 0 row)
    sampled_points = Y(2:end, 2);

    % Plot sampled Uc values as black dots
    plot(repmat(R, length(sampled_points), 1), sampled_points, ...
         '.k', 'MarkerSize', 2);
    drawnow;

    % Progress indicator every 50 steps
    if mod(i, 50) == 0
        fprintf('Progress: %d / %d steps (R = %.2e Ohm)\n', i, num_steps, R);
    end

end
elapsed = toc;
fprintf('Bifurcation diagram completed in %.1f s.\n', elapsed);

disp('Done.');
