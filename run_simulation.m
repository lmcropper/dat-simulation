% Function to run a single simulation
function intensity = run_simulation(params, mask_radius_min, mask_radius_max, verbose)
    % Optional verbose parameter for progress reporting
    if nargin < 4
        verbose = false;
    end
    
    % Extract parameters from struct
    lambda = params.lambda;
    a = params.a;
    w0 = params.w0;
    f = params.f;
    k = params.k;
    n = params.n;
    p = params.p;
    q = params.q;
    z_min = params.z_min;
    z_max = params.z_max;
    r_min = params.r_min;
    r_max = params.r_max;
    grid_points_z = params.grid_points_z;
    grid_points_r = params.grid_points_r;
    
    % Conversion factors
    z_scale = f^2 / (k*a^2);
    r_scale = f / (k*a);
    
    % Create normalized coordinate grids
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
