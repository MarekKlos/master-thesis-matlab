% =========================================================================
% FERRO_A: BIFURCATION DIAGRAM (VOLTAGE SWEEP)
% =========================================================================
% Generates a bifurcation diagram for Circuit A by sweeping the source
% voltage amplitude U and sampling the capacitor voltage Uc at each
% period of the fundamental frequency (Poincare section / stroboscopic map).
%
% Method: continuation — the final state of each simulation step is used
%         as the initial condition for the next step, which accelerates
%         convergence to steady state and preserves solution branches.
% =========================================================================
clear; clc; close all;

disp('Starting Bifurcation Diagram (Voltage Sweep)...');

% -------------------------------------------------------------------------
% 1. FIXED SYSTEM PARAMETERS
% -------------------------------------------------------------------------
alf = 0.0;          % Initial phase angle [rad]
Rs  = 1.5;          % Series resistance [Ohm]
Rm  = 50000;        % Core-loss (magnetizing) resistance [Ohm]
C   = 0.78E-6;      % Capacitance [F]
m1  = 4.079845e-04; % Linear coefficient of magnetization curve
mq  = 2.108275e-27; % Nonlinear coefficient of magnetization curve
omg = 100*pi;       % Angular frequency [rad/s]  (50 Hz)
T0  = 2*pi / omg;   % Fundamental period [s]  (0.02 s)

% -------------------------------------------------------------------------
% 2. BIFURCATION SWEEP SETTINGS
% -------------------------------------------------------------------------
U_start   = 10000;      % Start voltage [V]
U_end     = 200000;     % End voltage [V]
num_steps = 400;        % Number of voltage steps (higher = smoother diagram)
U_vector  = linspace(U_start, U_end, num_steps);

% Simulation length per step (cycle-based for consistency)
num_cycles     = 40;    % Total cycles to simulate per step
discard_cycles = 25;    % Transient cycles to discard before sampling
t_final        = num_cycles * T0;   % ~0.8 s per step

tspan   = [0, t_final];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
y0      = [280; 1530; 0];   % Initial condition [Psi; Uc; phase]

% -------------------------------------------------------------------------
% 3. FIGURE SETUP  (thesis formatting)
% -------------------------------------------------------------------------
FNAME  = 'Times New Roman';
FSIZE  = 12;
FTITLE = 14;

fig = figure('Name', 'Bifurcation Diagram: Voltage Sweep', ...
             'Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
set(gca, 'FontName', FNAME, 'FontSize', FSIZE);
xlabel('U [V]', 'FontName', FNAME);
ylabel('U_c [V]', 'FontName', FNAME);
title('Bifurcation Diagram', ...
      'FontName', FNAME, 'FontSize', FTITLE);
xlim([U_start, U_end]);

% -------------------------------------------------------------------------
% 4. MAIN CONTINUATION LOOP
% -------------------------------------------------------------------------
tic;
for i = 1:num_steps

    U = U_vector(i);

    % Recalculate parameters dependent on U
    p2 = 1/(C*(Rm+Rs)); p1 = U*p2;
    p3 = p2*m1*Rm;      p4 = p2*mq*Rm;
    b2 = Rm/(Rm+Rs);    b1 = U*b2;
    b3 = b2*m1*Rs;      b4 = b2*mq*Rs;

    ferro_ode = @(t, y) [
        b1*cos(y(3)+alf) - b2*y(2) - b3*y(1) - b4*y(1)^11;
        p1*cos(y(3)+alf) - p2*y(2) + p3*y(1) + p4*y(1)^11;
        omg
    ];

    % Integrate
    [T, Y] = ode15s(ferro_ode, tspan, y0, options);

    % Continuation: pass final state as next initial condition
    y0    = Y(end, :)';
    y0(3) = mod(y0(3), 2*pi);  % Wrap phase to [0, 2*pi]

    % Stroboscopic sampling at t = k*T0, discarding first transient cycles
    sample_times   = (discard_cycles : num_cycles) * T0;
    sampled_points = interp1(T, Y(:,2), sample_times, 'spline');

    % Plot sampled Uc values as black dots
    plot(repmat(U, 1, length(sampled_points)), sampled_points, ...
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
