% =========================================================================
% DUFFING FERRORESONANCE: PARAMETER SWEEP (Im x I0)  -- Circuit C
% =========================================================================
% Rows:    Im in {0.50, 0.75, ..., 5.00}  (step 0.25)  -> 19 values
% Columns: I0 in {0.00, 0.50, ..., 2.50}  (step 0.50)  ->  6 values
% Cell value: max |Uc| over the last 20% of the simulation (steady state)
% =========================================================================
clear; clc; close all;
global par;

% -------------------------------------------------------------------------
% 1. FIXED PARAMETERS  (match the Circuit C baseline from Chapter 5)
% -------------------------------------------------------------------------
C   = 40.53E-6;     % Parallel capacitance [F]   <-- change if needed
R   = 19.38E+6;     % Core-loss resistance [Ohm]
m1  = 0.32;         % Linear coefficient of magnetization curve
m3  = 0.9;          % Cubic  coefficient of magnetization curve
omg = 100*pi;       % Angular frequency [rad/s]

% par order: [Im; I0; C; R; m1; m3; omg]

% -------------------------------------------------------------------------
% 2. SWEEP RANGES
% -------------------------------------------------------------------------
Im_vec = 0.50 : 0.25 : 5.00;    % rows
I0_vec = 0.00 : 0.50 : 2.50;    % columns
nIm = length(Im_vec);
nI0 = length(I0_vec);

% -------------------------------------------------------------------------
% 3. SIMULATION SETTINGS
% -------------------------------------------------------------------------
t_final = 2.0;
tspan   = [0, t_final];
y0      = [0; 0; 0];                       % [Psi0; Uc0; phi0]
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

Uss_max = zeros(nIm, nI0);                 % max steady-state |Uc|
Utr_max = zeros(nIm, nI0);                 % max transient   |Uc| (extra)

% -------------------------------------------------------------------------
% 4. SWEEP
% -------------------------------------------------------------------------
total   = nIm * nI0;
counter = 0;
fprintf('Running %d simulations...\n', total);
tic;

for i = 1:nIm
    for j = 1:nI0
        counter = counter + 1;
        Im = Im_vec(i);
        I0 = I0_vec(j);
        par = [Im; I0; C; R; m1; m3; omg];

        [T, Y] = ode45(@duffing_ode, tspan, y0, options);

        x2_uc      = Y(:,2);
        idx_steady = T >= 0.8 * T(end);

        Utr_max(i,j) = max(abs(x2_uc));
        Uss_max(i,j) = max(abs(x2_uc(idx_steady)));

        fprintf('  [%3d/%3d]  Im=%.2f  I0=%.2f  -->  Uss,max = %8.2f V\n', ...
                counter, total, Im, I0, Uss_max(i,j));
    end
end
fprintf('Sweep finished in %.1f s.\n', toc);

% -------------------------------------------------------------------------
% 5. DISPLAY AS A MATLAB TABLE
% -------------------------------------------------------------------------
colNames = arrayfun(@(x) sprintf('I0_%.1f', x), I0_vec, 'UniformOutput', false);
rowNames = arrayfun(@(x) sprintf('Im_%.2f', x), Im_vec, 'UniformOutput', false);

T_uss = array2table(round(Uss_max, 1), 'VariableNames', colNames, ...
                                       'RowNames',      rowNames);

disp(' ');
disp('==============================================================');
disp(' Max steady-state |Uc| [V]   (rows = Im [A], columns = I0 [A])');
disp('==============================================================');
disp(T_uss);

% -------------------------------------------------------------------------
% 6. EXPORT TO CSV  (paste straight into Word/Excel)
% -------------------------------------------------------------------------
writetable(T_uss, 'Circuit_C_Im_I0_sweep.csv', 'WriteRowNames', true);
fprintf('Results saved to Circuit_C_Im_I0_sweep.csv\n');

% -------------------------------------------------------------------------
% 7. HEATMAP (optional, but very useful for spotting structure)
% -------------------------------------------------------------------------
figure('Name', 'Im-I0 Sweep', 'Color', 'w', 'Position', [100 100 700 700]);
imagesc(I0_vec, Im_vec, Uss_max);
set(gca, 'YDir', 'normal', 'FontName', 'Times New Roman', 'FontSize', 12);
colormap('hot');
cb = colorbar;
cb.Label.String   = 'Max steady-state |U_c| [V]';
cb.Label.FontName = 'Times New Roman';
xlabel('I_0 [A]',  'FontName', 'Times New Roman', 'FontSize', 13);
ylabel('I_m [A]',  'FontName', 'Times New Roman', 'FontSize', 13);
title('Max steady-state capacitor voltage  (Circuit C)', ...
      'FontName', 'Times New Roman', 'FontSize', 14);
xticks(I0_vec);
yticks(Im_vec(1:2:end));    % every second tick to avoid clutter

% =========================================================================
% LOCAL FUNCTION  (Duffing-type RHS, same as the main analysis script)
% =========================================================================
function dydt = duffing_ode(~, y)
    global par
    Im=par(1); I0=par(2); C=par(3); R=par(4);
    m1=par(5); m3=par(6); omg=par(7);

    dydt    = zeros(3,1);
    dydt(1) = y(2);
    dydt(2) = (1/C)*(Im*cos(y(3)) + I0) - (1/(R*C))*y(2) ...
              - (m3/C)*(y(1)^3) - (m1/C)*y(1);
    dydt(3) = omg;
end