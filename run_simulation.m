% Function to run a single simulation
function intensity = run_simulation(mask_radius_min, mask_radius_max, verbose)
    % Optional verbose parameter for progress reporting
    if nargin < 3
        verbose = false;
    end
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

    % Define the grid
    % Physical plot range in real units
    z_min = -3000e-6; % Axial distance in meters (micrometers shown as -2000)
    z_max = 1000e-6;   % Axial distance in meters
    r_min = -30e-6;  % Radial distance in meters
    r_max = 30e-6;   % Radial distance in meters

    % Conversion factors from normalized (u,v) to real units (z,r)
    % u is scaled axial distance, v is scaled radial distance
    % Typical Fresnel scaling: u ~ z*k*a^2/(2*f), v ~ k*a*r/f
    z_scale = f^2 / (k*a^2);
    r_scale = f / (k*a);       % Conversion factor: r = v * r_scale

    % Create normalized coordinate grids from real unit ranges
    grid_points_z = 400;
    grid_points_r = 400;
    u_values = linspace(z_min/z_scale, z_max/z_scale, grid_points_z);
    v_values = linspace(r_min/r_scale, r_max/r_scale, grid_points_r);
    [u, v] = meshgrid(u_values, v_values);

    % Spherical aberration phase term (Shvedov's Default)
    %phi = @(rho) -(a*rho).^4 * (1.5)^2 / (8*f^3*(1.5 - 1)^2); % Modify for lens index

    % Generalized phase term (From Seidel 3rd-order aberration structural polynomial)
    phi = @(rho) -((a*rho).^4 / (32*n*f^3*(n-1))) * (((n+2)/(n-1))*q^2 + (4*(n+1)*p*q) + ((3*n+2)*(n-1)*p^2) + ((n^3)/(n-1)));

    % Mask function: 1 in annulus region (between min and max radius), 0 outside
    mask = @(rho) (rho >= mask_radius_min) & (rho <= mask_radius_max);

    % Intensity calculation
    intensity = zeros(size(u));
    total_points = numel(u);
    
    if verbose
        fprintf('  Calculating intensity for %d points...\n', total_points);
    end

    for i = 1:numel(u)
        % Progress update every 10% or at the start
        if verbose && (mod(i, ceil(total_points/10)) == 0 || i == 1)
            fprintf('    Progress: %d/%d points (%.0f%%)\n', i, total_points, 100*i/total_points);
        end
        
        % Parameters for current u and v
        current_u = u(i);
        current_v = v(i);
        
        % Define the integrand
        integrand = @(rho) mask(rho) .* exp(-rho.^2 / (w0/a)^2) .* ...
            exp(-1i * (k * phi(rho) - current_u * rho.^2 / 2)) .* ...
            besselj(0, current_v * rho) .* rho;
        
        % Numerical integration over rho (0 to 1)
        intensity(i) = abs(integral(integrand, 0, 1))^2;
    end

    % Normalize intensity
    intensity = intensity / max(intensity(:));
    
    % Flatten the 2D array for return
    intensity = intensity(:);
end
