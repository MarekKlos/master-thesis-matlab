% =========================================================================
% DUFFING FERRORESONANCE: BIFURCATION DIAGRAM (FORCING-AMPLITUDE SWEEP)
% =========================================================================
% Generates a bifurcation diagram for Circuit C by sweeping the
% sinusoidal current-source amplitude Im, and sampling the capacitor
% voltage Uc at each period of the fundamental frequency (Poincare
% section / stroboscopic map).
%
% Method: continuation — the final state of each simulation step is used
%         as the initial condition for the next step.
%
% Sampling note: Uc is sampled at exactly k*T0 by placing the sample
%         instants directly in the solver's output vector (tspan). ode45
%         then returns the state at those instants using its own
%         internally consistent interpolant — no external interp1 spline,
%         so a true period-1 orbit collapses to a single point.
% =========================================================================
clear; clc; close all;

disp('Starting Bifurcation Diagram - Duffing (Forcing-Amplitude Sweep)...');

% -------------------------------------------------------------------------
% 1. FIXED SYSTEM PARAMETERS
% -------------------------------------------------------------------------
I0  = 0.5;          % DC component of the current source [A]
C   = 40.53E-6;     % Parallel capacitance [F]
R   = 19.38E+6;     % Core-loss (parallel) resistance [Ohm]
m1  = 0.32;         % Linear coefficient of magnetization curve
m3  = 0.9;          % Cubic coefficient of magnetization curve
omg = 100*pi;       % Angular frequency [rad/s]  (50 Hz)
T0  = 2*pi / omg;   % Fundamental period [s]    (0.02 s)

% -------------------------------------------------------------------------
% 2. BIFURCATION SWEEP SETTINGS
% -------------------------------------------------------------------------
Im_start  = 0.0;        % Start forcing amplitude [A]
Im_end    = 5.0;        % End forcing amplitude [A]
num_steps = 400;        % Number of sweep points
Im_vector = linspace(Im_start, Im_end, num_steps);

% Simulation length per step (cycle-based for consistency)
num_cycles     = 40;    % Total cycles to simulate per step
discard_cycles = 25;    % Transient cycles to discard before sampling

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
y0      = [0; 0; 0];    % Initial condition [Psi; Uc; phase]

% Sample instants k*T0 (held constant across the sweep). These are placed
% directly in tspan so the solver reports Uc at exactly these times.
sample_times = (discard_cycles : num_cycles) * T0;

% -------------------------------------------------------------------------
% 3. FIGURE SETUP  (thesis formatting)
% -------------------------------------------------------------------------
FNAME  = 'Times New Roman';
FSIZE  = 12;
FTITLE = 14;

fig = figure('Name', 'Bifurcation Diagram: Forcing-Amplitude Sweep', ...
             'Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
xlabel('I_m [A]',         'FontName', FNAME);
ylabel('U_c [V]', 'FontName', FNAME);
title('Bifurcation Diagram', ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlim([Im_start, Im_end]);

% -------------------------------------------------------------------------
% 4. MAIN CONTINUATION LOOP
% -------------------------------------------------------------------------
tic;
for i = 1:num_steps

    Im = Im_vector(i);

    % Anonymous RHS (no globals needed for this smooth 3D system).
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
    y0(3) = mod(y0(3), 2*pi);   % Wrap phase to [0, 2*pi]

    % Sampled Uc values (skip the prepended t = 0 row)
    sampled_points = Y(2:end, 2);

    % Plot sampled Uc values as black dots
    plot(repmat(Im, length(sampled_points), 1), sampled_points, ...
         '.k', 'MarkerSize', 2);
    drawnow;

    % Progress indicator every 50 steps
    if mod(i, 50) == 0
        fprintf('Progress: %d / %d steps (Im = %.2f A)\n', i, num_steps, Im);
    end

end
elapsed = toc;
fprintf('Bifurcation diagram completed in %.1f s.\n', elapsed);

disp('Done.');
