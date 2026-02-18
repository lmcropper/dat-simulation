# Simulation Code Update Summary

## Changes Made

All simulation functions have been updated to accept parameters via a `params` struct instead of having them hardcoded internally.

## Updated Files

### Simulation Functions (ALL updated):
1. **run_simulation.m** - Original method
2. **run_simulation_optimized.m** - Vectorized method
3. **run_simulation_fast.m** - Vectorized + parallel method
4. **run_simulation_gpu.m** - GPU method
5. **run_simulation_gpu_ultra.m** - GPU ultra method

**New signature:** `intensity = run_simulation_*(params, mask_radius_min, mask_radius_max, verbose)`

### Support Scripts (updated):
- **run_sweep.m** - Now uses `create_simulation_params()` and passes params to simulation
- **benchmark_simulations.m** - Now uses `create_simulation_params()` and passes params

### New Files Created:
- **create_simulation_params.m** - Helper function to create params struct with defaults
- **example_custom_params.m** - Examples showing how to use custom parameters
- **PARAMETERS_GUIDE.m** - Comprehensive guide to all parameters and usage

## Benefits

1. **Centralized Control**: All parameters now controlled from calling function (e.g., sweep script)
2. **Flexibility**: Easy to modify parameters without editing simulation functions
3. **Consistency**: All simulation methods use the same parameter interface
4. **Convenience**: Helper function provides sensible defaults
5. **Documentation**: New guide explains all parameters

## Migration Guide

### Before:
```matlab
intensity = run_simulation_fast(0.0, 1.0, true);
```

### After:
```matlab
params = create_simulation_params();
intensity = run_simulation_fast(params, 0.0, 1.0, true);
```

### With Custom Parameters:
```matlab
% Option 1: Use helper with name-value pairs
params = create_simulation_params('lambda', 633e-9, 'f', 50e-3);
intensity = run_simulation_fast(params, 0.0, 1.0, true);

% Option 2: Modify after creation
params = create_simulation_params();
params.lambda = 633e-9;
params.k = 2*pi/params.lambda;  % Remember to update k!
intensity = run_simulation_fast(params, 0.0, 1.0, true);
```

## Parameters Available

### Physical:
- `lambda` - Wavelength (m)
- `a` - Aperture radius (m)
- `w0` - Beam waist (m)
- `f` - Focal length (m)
- `n` - Refractive index
- `p` - Position factor
- `q` - Shape factor
- `k` - Wave number (auto-calculated)

### Grid:
- `z_min`, `z_max` - Axial range (m)
- `r_min`, `r_max` - Radial range (m)
- `grid_points_z` - Z axis samples
- `grid_points_r` - R axis samples
- `Nrho` - Integration samples

## Example Workflows

See **example_custom_params.m** for working examples including:
- Default parameters
- Custom wavelength
- Custom focal length/aperture
- Annular apertures
- High resolution grids

## Testing

Run **benchmark_simulations.m** to verify all methods work correctly with new parameter system.
