% Optical Intensity Pattern Calculator
% This script runs a single simulation and displays the results
% It calls the run_simulation function for the core calculation

% Run simulation with current mask parameters
mask_radius_min = 0.0; % Inner radius (normalized, 0 to 1)
mask_radius_max = 1.0; % Outer radius (normalized, 0 to 1)

% Parameters
lambda = 532e-9; % Wavelength in meters
a = 12e-3; % Aperture radius in meters
w0 = 1.1e-3; % Beam waist in meters
f = 25.4e-3; % Focal length in meters
P = 1; % Beam power in arbitrary units
k = 2*pi/lambda; % Wave number
n = 1.5;

p = -1; % Position Factor equal to 1-(2f/Si). For a collimated input beam, equal to 1.
q = -1; % Shape factor. Equal to 1 for a plano-convex lens with the plane side facing the source.

% Define the grid for visualization
z_min = -6000e-6;
z_max = -2000e-6;
r_min = -30e-6;
r_max = 30e-6;

z_scale = f^2 / (k*a^2);
r_scale = f / (k*a);

u_values = linspace(z_min/z_scale, z_max/z_scale, 1000);
v_values = linspace(r_min/r_scale, r_max/r_scale, 400);

% Run simulation
intensity_flat = run_simulation(mask_radius_min, mask_radius_max);

% Reshape back to grid
intensity = reshape(intensity_flat, length(v_values), length(u_values));

% Normalize for display with gamma correction
gamma = 0.5;
intensity_gamma = intensity.^gamma;

% Plot the intensity pattern
figure;
imagesc(u_values*z_scale*1e6, v_values*r_scale*1e6, intensity_gamma);
colormap(hot);
colorbar;
xlabel('Axial Distance z (μm)');
ylabel('Radial Distance r (μm)');
title('Optical Intensity Pattern with Spherical Aberration');

% Plot intensity along the optical axis (r = 0)
figure;
[~, v_zero_idx] = min(abs(v_values)); % Find index where v is closest to 0
intensity_axis = intensity(v_zero_idx, :); % Extract intensity at r = 0
plot(u_values*z_scale*1e6, intensity_axis, 'LineWidth', 2);
xlabel('Axial Distance z (μm)');
ylabel('Normalized Intensity');
title('Intensity Distribution Along Optical Axis');
grid on;

