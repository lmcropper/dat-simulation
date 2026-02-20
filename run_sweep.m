% Sweep Script: Run multiple simulations with varying parameters
% This script automates parameter sweeps and saves results as a 3D matrix

clear all;
close all;

% ========== SIMULATION PARAMETERS ==========
% Create parameter struct using helper function
params = create_simulation_params();

% Override specific parameters if needed
params.z_min = -50e-6;
params.z_max = 200e-6;
params.r_max = 50e-6;
params.r_min = -50e-6;
params.grid_points_z = 1100;
params.grid_points_r = 800;
params.Nrho = 2048;     % Keep this high, reducing it causes spatial aliasing in the simulation
params.lambda = 405e-9;  % Wavelength (m)
params.a = 15e-3;        % Aperture radius (m)
params.w0 = 1.1e-3;      % Beam waist (m)
params.f = 40e-3;      % Focal length (m)
params.p = -1;        % Position factor of system
params.q = -1;        % Shape factor of lens

% ========== SWEEP PARAMETERS ==========
% Define which parameter to sweep and the range
sweep_param = 'mask_radius_min';  % Options: 'mask_radius_min', 'mask_radius_max', 'annulus_fixed_width'
sweep_start = 0e-3 / params.a;     % Start value (0%)
sweep_end = params.a / params.a;       % End value (70%)
sweep_increment = 0.1e-3 / params.a; % Increment (1%)

% Fixed parameters
fixed_mask_radius_min = 0.0;   % Inner radius (when not sweeping)
fixed_mask_radius_max = 11e-3 / params.a;   % Outer radius (when not sweeping)
annulus_width = 1e-3 / params.a;           % Fixed annulus width (for 'annulus_fixed_width' mode)

% Extract for convenience
lambda = params.lambda;
a = params.a;
w0 = params.w0;
f = params.f;
k = params.k;
z_min = params.z_min;
z_max = params.z_max;
r_min = params.r_min;
r_max = params.r_max;
grid_points_z = params.grid_points_z;
grid_points_r = params.grid_points_r;

z_scale = f^2 / (k*a^2);
r_scale = f / (k*a);

u_values = linspace(z_min/z_scale, z_max/z_scale, grid_points_z);
v_values = linspace(r_min/r_scale, r_max/r_scale, grid_points_r);

% ========== CREATE SWEEP ARRAY ==========
sweep_values = sweep_start : sweep_increment : sweep_end;
num_sweep_points = length(sweep_values);

fprintf('Sweep configuration:\n');
fprintf('  Parameter: %s\n', sweep_param);
fprintf('  Range: %.2f to %.2f\n', sweep_start, sweep_end);
fprintf('  Increment: %.3f\n', sweep_increment);
fprintf('  Total simulations: %d\n\n', num_sweep_points);

% ========== PRE-ALLOCATE RESULTS MATRIX ==========
% 3D matrix: (num_r, num_z, num_sweep_points)
intensity_sweep = zeros(grid_points_r, grid_points_z, num_sweep_points);
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
    elseif strcmp(sweep_param, 'annulus_fixed_width')
        % Sweep outer radius with fixed annulus width
        mask_radius_max = sweep_value;
        mask_radius_min = max(0.0, sweep_value - annulus_width);  % Ensure non-negative
    else
        error('Unknown sweep parameter: %s', sweep_param);
    end
    
    % Log progress
    fprintf('Running simulation %d/%d: %s = %.3f', sweep_idx, num_sweep_points, sweep_param, sweep_value);
    if strcmp(sweep_param, 'annulus_fixed_width')
        fprintf(' (inner = %.3f, width = %.3f)', mask_radius_min, mask_radius_max - mask_radius_min);
    end
    fprintf('\n');
    
    % Run simulation with verbose output
    intensity_flat = run_simulation_gpu_ultra(params, mask_radius_min, mask_radius_max, false);
    
    % Reshape to 2D grid (grid_points_r x grid_points_z as created by simulation)
    intensity_2d = reshape(intensity_flat, grid_points_r, grid_points_z);
    
    % Store in 3D matrix
    intensity_sweep(:, :, sweep_idx) = intensity_2d;
    sweep_param_values(sweep_idx) = sweep_value;
end

fprintf('\nAll simulations complete!\n');

% ========== SAVE RESULTS ==========
% Create filename with timestamp
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
filename = sprintf('Outputs/sweep_results_%s_%s.mat', sweep_param, char(timestamp));

% Save with -v7.3 flag for HDF5 format (enables partial loading in sweep_viewer_fast)
fprintf('Saving results (HDF5 format for fast loading)...\n');
save(filename, 'intensity_sweep', 'sweep_param_values', 'sweep_param', ...
    'u_values', 'v_values', 'z_scale', 'r_scale', ...
    'lambda', 'a', 'w0', 'f', 'k', ...
    'z_min', 'z_max', 'r_min', 'r_max', ...
    'sweep_start', 'sweep_end', 'sweep_increment', 'annulus_width', '-v7.3');

fprintf('Results saved to: %s\n', filename);
fprintf('Matrix size: %d x %d x %d\n', size(intensity_sweep, 1), size(intensity_sweep, 2), size(intensity_sweep, 3));
