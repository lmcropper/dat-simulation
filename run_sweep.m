% Sweep Script: Run multiple simulations with varying parameters
% This script automates parameter sweeps and saves results as a 3D matrix

clear all;
close all;

% ========== SWEEP PARAMETERS ==========
% Define which parameter to sweep and the range
sweep_param = 'mask_radius_min';  % Parameter to sweep: 'mask_radius_min' or 'mask_radius_max'
sweep_start = 0.0;     % Start value (0%)
sweep_end = 0.1;       % End value (70%)
sweep_increment = 0.01; % Increment (1%)

% Fixed parameters
fixed_mask_radius_min = 0.0;   % Inner radius (when not sweeping)
fixed_mask_radius_max = 1.0;   % Outer radius (when not sweeping)

% ========== SIMULATION PARAMETERS ==========
lambda = 532e-9;
a = 12e-3;
w0 = 1.1e-3;
f = 25.4e-3;
k = 2*pi/lambda;

z_scale = f^2 / (k*a^2);
r_scale = f / (k*a);

z_min = -6000e-6;
z_max = -2000e-6;
r_min = -30e-6;
r_max = 30e-6;

u_values = linspace(z_min/z_scale, z_max/z_scale, 1000);
v_values = linspace(r_min/r_scale, r_max/r_scale, 400);

% ========== CREATE SWEEP ARRAY ==========
sweep_values = sweep_start : sweep_increment : sweep_end;
num_sweep_points = length(sweep_values);

fprintf('Sweep configuration:\n');
fprintf('  Parameter: %s\n', sweep_param);
fprintf('  Range: %.2f to %.2f\n', sweep_start, sweep_end);
fprintf('  Increment: %.3f\n', sweep_increment);
fprintf('  Total simulations: %d\n\n', num_sweep_points);

% ========== PRE-ALLOCATE RESULTS MATRIX ==========
% 3D matrix: (num_v, num_u, num_sweep_points)
% Note: run_simulation creates 400x400 grids internally
intensity_sweep = zeros(400, 400, num_sweep_points);
sweep_param_values = zeros(1, num_sweep_points);

% ========== RUN SIMULATIONS ==========
for sweep_idx = 1:num_sweep_points
    sweep_value = sweep_values(sweep_idx);
    
    % Set parameters based on what we're sweeping
    if strcmp(sweep_param, 'mask_radius_min')
        mask_radius_min = sweep_value;
        mask_radius_max = fixed_mask_radius_max;
    elseif strcmp(sweep_param, 'mask_radius_max')
        mask_radius_min = fixed_mask_radius_min;
        mask_radius_max = sweep_value;
    else
        error('Unknown sweep parameter: %s', sweep_param);
    end
    
    % Log progress
    fprintf('Running simulation %d/%d: %s = %.3f\n', sweep_idx, num_sweep_points, sweep_param, sweep_value);
    
    % Run simulation with verbose output
    intensity_flat = run_simulation(mask_radius_min, mask_radius_max, true);
    
    % Reshape to 2D grid (400x400 as created by run_simulation)
    intensity_2d = reshape(intensity_flat, 400, 400);
    
    % Store in 3D matrix
    intensity_sweep(:, :, sweep_idx) = intensity_2d;
    sweep_param_values(sweep_idx) = sweep_value;
end

fprintf('\nAll simulations complete!\n');

% ========== SAVE RESULTS ==========
% Create filename with timestamp
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
filename = sprintf('sweep_results_%s_%s.mat', sweep_param, char(timestamp));

save(filename, 'intensity_sweep', 'sweep_param_values', 'sweep_param', ...
    'u_values', 'v_values', 'z_scale', 'r_scale', ...
    'lambda', 'a', 'w0', 'f', 'k', ...
    'z_min', 'z_max', 'r_min', 'r_max', ...
    'sweep_start', 'sweep_end', 'sweep_increment');

fprintf('Results saved to: %s\n', filename);
fprintf('Matrix size: %d x %d x %d\n', size(intensity_sweep, 1), size(intensity_sweep, 2), size(intensity_sweep, 3));
