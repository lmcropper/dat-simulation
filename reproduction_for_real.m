% Parameters
lambda = 633e-9; % Wavelength in meters
a = 10e-3; % Aperture radius in meters
w0 = 1.4e-3; % Beam waist in meters
f = 30e-3; % Focal length in meters
P = 1; % Beam power in arbitrary units
k = 2*pi/lambda; % Wave number

% Define the grid
u_values = linspace(-1200, 200, 200); % u coordinate
v_values = linspace(-100, 100, 200); % v coordinate
[u, v] = meshgrid(u_values, v_values);

% Spherical aberration phase term
phi = @(rho) -(a*rho).^4 * (1.5)^2 / (8*f^3*(1.5 - 1)^2); % Modify for lens index

% Intensity calculation
intensity = zeros(size(u));
for i = 1:numel(u)
    % Parameters for current u and v
    current_u = u(i);
    current_v = v(i);
    
    % Define the integrand
    integrand = @(rho) exp(-rho.^2 / (w0/a)^2) .* ...
        exp(-1i * (k * phi(rho) - current_u * rho.^2 / 2)) .* ...
        besselj(0, current_v * rho) .* rho;
    
    % Numerical integration over rho (0 to 1)
    intensity(i) = abs(integral(integrand, 0, 1))^2;
end

% Normalize intensity
intensity = intensity / max(intensity(:));

gamma = 0.5;
intensity_gamma = intensity.^gamma;

map = [0 0 0
    0.5 0 0
    0.7 0 0
    0.8 0 0
    0.9 0.1 0.1
    0.95 0.2 0.2
    1 1 1];

% Plot the intensity pattern
figure;
imagesc(u_values, v_values, intensity_gamma);
colormap(hot);
colorbar;
xlabel('u');
ylabel('v');
title('Optical Intensity Pattern with Spherical Aberration');
% set(gca, 'color', 'none');
% set(gca,'XColor', 'none','YColor','none');
% grid off;

