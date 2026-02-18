% SIMULATION PARAMETERS GUIDE
% ===========================
% All simulation functions now accept a 'params' struct as the first argument.
% This allows all parameters to be controlled from the calling function.

%% FUNCTION SIGNATURES (ALL METHODS):
% intensity = run_simulation(params, mask_radius_min, mask_radius_max, verbose)
% intensity = run_simulation_optimized(params, mask_radius_min, mask_radius_max, verbose)
% intensity = run_simulation_fast(params, mask_radius_min, mask_radius_max, verbose)
% intensity = run_simulation_gpu(params, mask_radius_min, mask_radius_max, verbose)
% intensity = run_simulation_gpu_ultra(params, mask_radius_min, mask_radius_max, verbose)

%% CREATING PARAMETERS:
% Use the helper function:
params = create_simulation_params();

% Or create manually:
params = struct();
params.lambda = 532e-9;           % Wavelength (m)
params.a = 12e-3;                 % Aperture radius (m)
params.w0 = 1.1e-3;               % Beam waist (m)
params.f = 25.4e-3;               % Focal length (m)
params.k = 2*pi/params.lambda;    % Wave number
params.n = 1.5;                   % Refractive index
params.p = -1;                    % Position factor
params.q = -1;                    % Shape factor
params.z_min = -3000e-6;          % Min z coordinate (m)
params.z_max = 1000e-6;           % Max z coordinate (m)
params.r_min = -30e-6;            % Min r coordinate (m)
params.r_max = 30e-6;             % Max r coordinate (m)
params.grid_points_z = 400;       % Number of z grid points
params.grid_points_r = 400;       % Number of r grid points
params.Nrho = 500;                % Radial integration samples

%% CUSTOMIZING PARAMETERS:
% Option 1: Use helper with name-value pairs
params = create_simulation_params('lambda', 633e-9, 'f', 50e-3);

% Option 2: Modify after creation
params = create_simulation_params();
params.lambda = 633e-9;
params.k = 2*pi/params.lambda;  % Update k when changing lambda
params.f = 50e-3;

%% PARAMETER DESCRIPTIONS:

% PHYSICAL PARAMETERS:
% lambda     - Wavelength in meters (default: 532e-9 = 532nm green)
% a          - Aperture radius in meters (default: 12e-3 = 12mm)
% w0         - Gaussian beam waist in meters (default: 1.1e-3 = 1.1mm)
% f          - Focal length in meters (default: 25.4e-3 = 25.4mm)
% n          - Refractive index of lens material (default: 1.5)
% p          - Position factor for aberration (default: -1)
% q          - Shape factor for aberration (default: -1)
% k          - Wave number = 2*pi/lambda (auto-calculated)

% GRID PARAMETERS:
% z_min, z_max       - Axial extent of observation region (m)
% r_min, r_max       - Radial extent of observation region (m)
% grid_points_z      - Number of samples along z axis
% grid_points_r      - Number of samples along r axis
% Nrho               - Number of integration points in radial direction
%                      (higher = more accurate but slower)

%% EXAMPLES OF USE:

% Example 1: Default parameters
params = create_simulation_params();
intensity = run_simulation_fast(params, 0.0, 1.0, true);

% Example 2: Different wavelength
params = create_simulation_params('lambda', 633e-9);
intensity = run_simulation_fast(params, 0.0, 1.0, true);

% Example 3: Higher resolution
params = create_simulation_params('grid_points_z', 800, 'grid_points_r', 800);
intensity = run_simulation_gpu_ultra(params, 0.0, 1.0, true);

% Example 4: Smaller observation region (zoom in near focus)
params = create_simulation_params('z_min', -500e-6, 'z_max', 500e-6, ...
                                  'r_min', -10e-6, 'r_max', 10e-6);
intensity = run_simulation_fast(params, 0.0, 1.0, true);

% Example 5: Run parameter sweep (see run_sweep.m for full example)
params = create_simulation_params();
params.z_min = -6000e-6;
params.z_max = -2000e-6;
for mask_val = 0:0.1:0.5
    intensity = run_simulation_gpu_ultra(params, mask_val, 1.0, true);
    % Process results...
end

%% FILES IN THIS PROJECT:

% SIMULATION FUNCTIONS (choose one based on speed/hardware):
% - run_simulation.m              - Original (slow, uses integral())
% - run_simulation_optimized.m    - Vectorized (20-30x faster)
% - run_simulation_fast.m         - Vectorized + parallel (80-150x faster)
% - run_simulation_gpu.m          - GPU accelerated (200-500x faster)
% - run_simulation_gpu_ultra.m    - GPU ultra (500-1000x faster, needs more VRAM)

% UTILITY SCRIPTS:
% - create_simulation_params.m    - Helper to create parameter struct
% - run_sweep.m                   - Automated parameter sweep
% - benchmark_simulations.m       - Compare speed of all methods
% - example_custom_params.m       - Examples of using custom parameters

%% TYPICAL WORKFLOW:

% 1. Create parameters
params = create_simulation_params();

% 2. Customize if needed
params.lambda = 633e-9;  % Change wavelength
params.k = 2*pi/params.lambda;  % Update k

% 3. Choose simulation method based on hardware
if gpuDeviceCount > 0
    intensity = run_simulation_gpu_ultra(params, 0, 1, true);
else
    intensity = run_simulation_fast(params, 0, 1, true);
end

% 4. Reshape and plot
intensity_2d = reshape(intensity, params.grid_points_r, params.grid_points_z);
imagesc(intensity_2d); colormap hot; colorbar;
