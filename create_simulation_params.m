% Helper function to create simulation parameters struct
function params = create_simulation_params(varargin)
    % Create default parameters
    params = struct();
    
    % Physical parameters
    params.lambda = 532e-9;       % Wavelength (m)
    params.a = 12e-3;             % Aperture radius (m)
    params.w0 = 1.1e-3;           % Beam waist (m)
    params.f = 25.4e-3;           % Focal length (m)
    params.n = 1.5;               % Refractive index
    params.p = -1;                % Position factor
    params.q = -1;                % Shape factor
    
    % Derived
    params.k = 2*pi/params.lambda;
    
    % Grid parameters (observation space)
    params.z_min = -3000e-6;      % Min axial distance (m)
    params.z_max = 1000e-6;       % Max axial distance (m)
    params.r_min = -30e-6;        % Min radial distance (m)
    params.r_max = 30e-6;         % Max radial distance (m)
    params.grid_points_z = 400;   % Number of z samples
    params.grid_points_r = 400;   % Number of r samples
    
    % Integration parameters
    params.Nrho = 500;            % Radial integration samples
    
    % Parse optional name-value pairs to override defaults
    for i = 1:2:length(varargin)
        if isfield(params, varargin{i})
            params.(varargin{i}) = varargin{i+1};
        else
            warning('Unknown parameter: %s', varargin{i});
        end
    end
    
    % Recalculate k if lambda was changed
    if any(strcmp(varargin, 'lambda'))
        params.k = 2*pi/params.lambda;
    end
end
