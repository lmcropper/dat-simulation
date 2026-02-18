% Example: How to run simulations with custom parameters
clear; clc;

fprintf('=== Example: Running simulations with custom parameters ===\n\n');

%% Example 1: Using default parameters
fprintf('Example 1: Using default parameters\n');
params = create_simulation_params();

% Run a simulation
mask_min = 0.0;
mask_max = 1.0;
intensity = run_simulation_fast(params, mask_min, mask_max, true);

% Reshape and plot
intensity_2d = reshape(intensity, params.grid_points_r, params.grid_points_z);

% Physical coordinates for plotting
z_scale = params.f^2 / (params.k * params.a^2);
r_scale = params.f / (params.k * params.a);
z_axis = linspace(params.z_min, params.z_max, params.grid_points_z) * 1e6; % micrometers
r_axis = linspace(params.r_min, params.r_max, params.grid_points_r) * 1e6; % micrometers

figure('Color', 'w');
imagesc(z_axis, r_axis, intensity_2d);
axis image; colormap hot; colorbar;
title('Example 1: Default Parameters');
xlabel('z (\mum)'); ylabel('r (\mum)');

%% Example 2: Custom wavelength (633nm instead of 532nm)
fprintf('\n\nExample 2: Custom wavelength (633nm red laser)\n');
params_red = create_simulation_params('lambda', 633e-9);

intensity_red = run_simulation_fast(params_red, mask_min, mask_max, true);
intensity_red_2d = reshape(intensity_red, params_red.grid_points_r, params_red.grid_points_z);

figure('Color', 'w');
imagesc(z_axis, r_axis, intensity_red_2d);
axis image; colormap hot; colorbar;
title('Example 2: 633nm Wavelength');
xlabel('z (\mum)'); ylabel('r (\mum)');

%% Example 3: Different focal length and aperture
fprintf('\n\nExample 3: Longer focal length (50mm) and larger aperture (20mm)\n');
params_long = create_simulation_params('f', 50e-3, 'a', 10e-3);

intensity_long = run_simulation_fast(params_long, mask_min, mask_max, true);
intensity_long_2d = reshape(intensity_long, params_long.grid_points_r, params_long.grid_points_z);

z_axis_long = linspace(params_long.z_min, params_long.z_max, params_long.grid_points_z) * 1e6;
r_axis_long = linspace(params_long.r_min, params_long.r_max, params_long.grid_points_r) * 1e6;

figure('Color', 'w');
imagesc(z_axis_long, r_axis_long, intensity_long_2d);
axis image; colormap hot; colorbar;
title('Example 3: f=50mm, D=20mm');
xlabel('z (\mum)'); ylabel('r (\mum)');

%% Example 4: Annular aperture (mask inner 30%)
fprintf('\n\nExample 4: Annular aperture (blocking inner 30%%)\n');
mask_min_annular = 0.3;
mask_max_annular = 1.0;

intensity_annular = run_simulation_fast(params, mask_min_annular, mask_max_annular, true);
intensity_annular_2d = reshape(intensity_annular, params.grid_points_r, params.grid_points_z);

figure('Color', 'w');
imagesc(z_axis, r_axis, intensity_annular_2d);
axis image; colormap hot; colorbar;
title('Example 4: Annular Aperture (inner 30% blocked)');
xlabel('z (\mum)'); ylabel('r (\mum)');

%% Example 5: Higher resolution grid
fprintf('\n\nExample 5: Higher resolution (800x800 grid)\n');
params_hires = create_simulation_params('grid_points_z', 800, 'grid_points_r', 800);

intensity_hires = run_simulation_fast(params_hires, mask_min, mask_max, true);
intensity_hires_2d = reshape(intensity_hires, params_hires.grid_points_r, params_hires.grid_points_z);

z_axis_hires = linspace(params_hires.z_min, params_hires.z_max, params_hires.grid_points_z) * 1e6;
r_axis_hires = linspace(params_hires.r_min, params_hires.r_max, params_hires.grid_points_r) * 1e6;

figure('Color', 'w');
imagesc(z_axis_hires, r_axis_hires, intensity_hires_2d);
axis image; colormap hot; colorbar;
title('Example 5: High Resolution (800x800)');
xlabel('z (\mum)'); ylabel('r (\mum)');

fprintf('\n\nAll examples complete!\n');
