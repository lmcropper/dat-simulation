% Ultra-fast GPU version with maximum parallelization
function intensity = run_simulation_gpu_ultra(params, mask_radius_min, mask_radius_max, verbose)
    if nargin < 4
        verbose = false;
    end
    
    % Check GPU
    if gpuDeviceCount == 0
        error('No GPU device found.');
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
    
    if verbose
        gpu = gpuDevice();
        fprintf('  GPU: %s\n', gpu.Name);
        fprintf('  Computing full 3D array on GPU...\n');
    end
    
    % Move all to GPU
    rho_gpu = gpuArray(single(rho));  % Use single precision for speed
    u_gpu = gpuArray(single(u_values));
    v_gpu = gpuArray(single(v_values));
    
    % Pre-compute 1D arrays
    mask_vals = gpuArray(single((rho >= mask_radius_min) & (rho <= mask_radius_max)));
    
    %PHI FOR SPHERICAL ABERRATION
    phi_const = single(-((a^4) / (32*n*f^3*(n-1))) * ...
               (((n+2)/(n-1))*q^2 + (4*(n+1)*p*q) + ((3*n+2)*(n-1)*p^2) + ((n^3)/(n-1))));
    phi_vals = phi_const * (rho_gpu.^4);

    %PHI FOR NO ABERRATION
    %phi_vals = gpuArray(single(zeros(size(rho_gpu))));  % No aberration
    
    % gauss_vals = exp(-rho_gpu.^2 / (w0/a)^2);
    gauss_vals = 1;
    base_vals = mask_vals .* gauss_vals .* rho_gpu;
    
    % Create 3D grids: (grid_points_r x grid_points_z x Nrho)
    % This is memory-intensive but maximally parallel
    
    % Reshape for 3D broadcasting
    v_3d = reshape(v_gpu, [grid_points_r, 1, 1]);      % r dimension
    u_3d = reshape(u_gpu, [1, grid_points_z, 1]);      % z dimension  
    rho_3d = reshape(rho_gpu, [1, 1, Nrho]);           % integration dimension
    base_3d = reshape(base_vals, [1, 1, Nrho]);
    phi_3d = reshape(phi_vals, [1, 1, Nrho]);
    
    if verbose
        fprintf('    Computing phase array...\n');
    end
    
    % Compute full 3D phase: (1 x grid_points_z x Nrho)
    rho_sq = rho_3d.^2;
    phase_3d = exp(-1i * (k * phi_3d - u_3d .* rho_sq / 2));
    
    % Multiply with base: (1 x grid_points_z x Nrho)
    common_3d = phase_3d .* base_3d;
    
    if verbose
        fprintf('    Computing Bessel functions...\n');
    end
    
    % Bessel arguments: (grid_points_r x 1 x Nrho)
    bessel_args = v_3d .* rho_3d;
    J0_3d = besselj(0, bessel_args);
    
    if verbose
        fprintf('    Computing full integrand and integrating...\n');
    end
    
    % Full integrand: (grid_points_r x grid_points_z x Nrho)
    integrand_3d = J0_3d .* common_3d;  % Broadcasting
    
    % Integrate over rho (dimension 3) using trapezoidal rule
    % Manual trapz for GPU efficiency
    drho = rho_gpu(2) - rho_gpu(1);
    intensity_gpu = drho * (sum(integrand_3d(:,:,2:end-1), 3) + ...
                            0.5 * (integrand_3d(:,:,1) + integrand_3d(:,:,end)));
    
    % Compute intensity
    intensity_gpu = abs(intensity_gpu).^2;
    
    if verbose
        fprintf('    Gathering results...\n');
    end
    
    % Transfer back
    intensity = gather(intensity_gpu);
    intensity = double(intensity);
    % intensity = intensity / max(intensity(:));
    intensity = intensity(:);
    
    if verbose
        fprintf('    Done!\n');
    end
end
