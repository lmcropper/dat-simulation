% Maximum speed: vectorized + parallel
function intensity = run_simulation_fast(params, mask_radius_min, mask_radius_max, verbose)
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
    Nrho = params.Nrho;

    % Grid
    z_scale = f^2 / (k*a^2);
    r_scale = f / (k*a);
    
    u_values = linspace(z_min/z_scale, z_max/z_scale, grid_points_z);
    v_values = linspace(r_min/r_scale, r_max/r_scale, grid_points_r);
    rho = linspace(0, 1, Nrho);
    
    % Pre-compute 1D arrays
    mask_vals = double((rho >= mask_radius_min) & (rho <= mask_radius_max));
    phi_vals = -((a*rho).^4 / (32*n*f^3*(n-1))) * ...
               (((n+2)/(n-1))*q^2 + (4*(n+1)*p*q) + ((3*n+2)*(n-1)*p^2) + ((n^3)/(n-1)));
    gauss_vals = exp(-rho.^2 / (w0/a)^2);
    base_vals = mask_vals .* gauss_vals .* rho;
    
    if verbose
        fprintf('  Computing %d x %d field (parallel)...\n', grid_points_r, grid_points_z);
    end
    
    intensity = zeros(grid_points_r, grid_points_z);
    
    % Parallel loop over z (faster, less memory than parfor over all points)
    parfor iu = 1:grid_points_z
        u_val = u_values(iu);
        phase_u = exp(-1i * (k * phi_vals - u_val * rho.^2 / 2));
        common = base_vals .* phase_u;
        
        col_intensity = zeros(grid_points_r, 1);
        for iv = 1:grid_points_r
            v_val = v_values(iv);
            J0_vals = besselj(0, v_val * rho);
            field = trapz(rho, common .* J0_vals);
            col_intensity(iv) = abs(field)^2;
        end
        intensity(:, iu) = col_intensity;
    end
    
    intensity = intensity / max(intensity(:));
    intensity = intensity(:);
end
