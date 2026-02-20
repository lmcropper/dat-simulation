# Simulation Workflow Guide

This project simulates the intensity pattern of a laser beam as it passes through a lens and aperture. It calculates what you would see in the focal region of the beam (the bright spot where the lens focuses the light). 

The simulation supports multiple speed options:
- **CPU versions**: For computers without a graphics card
- **GPU versions**: For computers with NVIDIA GPUs (much faster)

This guide will show you how to set up parameters, run simulations, compare different conditions, and visualize the results.

## What This Simulation Does

When light passes through a circular opening (aperture) and then a lens, it creates a diffraction pattern in the focal region. This simulation calculates the intensity (brightness) at each point in that focal region by summing contributions from every part of the aperture. The result is a 2D intensity map that shows the bright spot and any rings around it.

You can modify:
- **The light**: wavelength (color)
- **The aperture**: size and shape (circular or ring-shaped)
- **The lens**: focal length, material, and quality
- **The observation region**: how far from focus and how wide an area to look at

## Requirements

- MATLAB R2016b or later
- For faster parallel simulations: Parallel Computing Toolbox
- For GPU acceleration (much faster): Parallel Computing Toolbox + compatible NVIDIA GPU

## Quick Start: Running Your First Simulation

This section will walk you through the complete process of running a single simulation with default settings, then viewing your results.

### Step 1: Create Default Parameters

Open the MATLAB command window and type:

```matlab
params = create_simulation_params();
```

This loads a set of sensible default parameters from the file `create_simulation_params.m`. These defaults include:
- Wavelength: 633 nm (red laser light)
- Aperture radius: 1 mm
- Beam waist: 1 mm
- Focal length: 50 mm
- And many others...

**Where to find the defaults**: Open the file `create_simulation_params.m` to see what the default values are. Look for lines like `params.lambda = ...` and `params.f = ...`

### Step 2: (Optional) Change One Parameter

If you want to try different wavelengths or lens properties BEFORE running your first simulation, you can modify parameters:

```matlab
params.lambda = 532e-9;    % Change to 532 nm (green)
params.k = 2*pi/params.lambda;  % Important: update this when changing lambda
```

Or change the focal length (how strong the lens is):

```matlab
params.f = 100e-3;  % Change focal length to 100 mm
```

⚠️ **Important**: If you change the wavelength (`lambda`), you MUST also update `k` with the line shown above. The variable `k` is the wave number and depends on the wavelength.

**Where to find all possible parameters**: Open `PARAMETERS_GUIDE.m` to see descriptions of every parameter you can modify, what it means physically, and what values are reasonable to try.

### Step 3: Choose and Run a Simulation Method

Depending on whether you have a GPU, choose one of these:

**If you have an NVIDIA GPU** (fastest):
```matlab
intensity = run_simulation_gpu_ultra(params, 0.0, 1.0, true);
```

**If you only have a CPU** (slower but still reasonable):
```matlab
intensity = run_simulation_fast(params, 0.0, 1.0, true);
```

The four numbers mean:
- `params`: your parameters (from Step 1)
- `0.0`: minimum radius of aperture mask (0.0 = full circle)
- `1.0`: maximum radius of aperture mask (1.0 = full aperture)
- `true`: normalize the intensity to brightness 0-1 (if `false`, you get raw values)

**What does this actually do?** This runs the diffraction integral calculation for every point in your observation grid. It may take a few seconds to a few minutes depending on which method you choose and your computer's speed.

### Step 4: Display Your Result

Once the simulation finishes, you have the intensity pattern as a long vector. Reshape it into a 2D image and display it:

```matlab
% Reshape into 2D grid (rows = radius, columns = axial distance)
intensity_2d = reshape(intensity, params.grid_points_r, params.grid_points_z);

% Display the result
figure;
imagesc(intensity_2d);
colormap hot;
colorbar;
xlabel('Axial distance (z)');
ylabel('Radial distance (r)');
title('Beam Intensity Profile');
```

This shows you a colormap image where red is bright (high intensity) and dark is dim (low intensity). You should see a bright central spot with possible rings around it.

## Choosing a Simulation Method

All methods compute the same result; they just differ in speed.

- **`run_simulation.m`**
  - Very slow (reference method)
  - Use only for checking if results are correct
  - Recommended: Don't use unless you're debugging

- **`run_simulation_fast.m`**
  - Vectorized CPU version with parallel processing
  - Recommended for CPU-only systems
  - Takes 10-60 seconds depending on grid resolution

- **`run_simulation_gpu_ultra.m`**
  - Full GPU parallelization
  - Recommended for systems with NVIDIA GPUs
  - Takes 1-5 seconds for the same problem
  - Requires VRAM; reduce `grid_points_z` or `grid_points_r` if you run out of memory

**Quick decision tree**:
```matlab
% Check if you have a GPU
if gpuDeviceCount > 0
    % You have a GPU, use the fastest method
    intensity = run_simulation_gpu_ultra(params, 0.0, 1.0, true);
else
    % No GPU, use parallel CPU version (still quite fast)
    intensity = run_simulation_fast(params, 0.0, 1.0, true);
end
```

## Understanding and Changing Parameters

Parameters control everything about your simulation: what light you're using, the size and shape of your aperture, what lens you're using, and which region you want to look at.

### Where to Find Parameter Descriptions

**File to open**: `PARAMETERS_GUIDE.m`

This file contains detailed descriptions of every parameter you can change, including:
- What each parameter means physically
- What units it uses
- What range of values is reasonable
- How changing it affects the simulation results

### Most Important Parameters to Understand

When you first start working with this code, these are the parameters you'll most likely want to change:

| Parameter | What it controls | Typical values | File to change |
|-----------|------------------|-----------------|-----------------|
| `lambda` | Color of light (wavelength) | 405e-9 to 1550e-9 (UV to infrared) | `create_simulation_params.m` |
| `a` | Size of the hole/aperture | 0.5e-3 to 5e-3 (0.5 mm to 5 mm) | `create_simulation_params.m` |
| `f` | Lens strength (focal length) | 10e-3 to 300e-3 (10 mm to 300 mm) | `create_simulation_params.m` |
| `w0` | Beam width at the aperture | Usually same as `a` or smaller | `create_simulation_params.m` |
| `n` | Lens material (refractive index) | 1.5 for glass, 1.33 for water | `create_simulation_params.m` |
| `grid_points_z` | Resolution along the focus axis | 100 to 2000 (more = sharper image) | `create_simulation_params.m` |
| `grid_points_r` | Resolution in radial direction | 100 to 2000 (more = sharper image) | `create_simulation_params.m` |

### Step-by-Step: How to Change a Parameter

**Option A: Quick change in MATLAB command window**

```matlab
% Load defaults
params = create_simulation_params();

% Change one parameter
params.lambda = 532e-9;  % Try green light instead of red
params.k = 2*pi / params.lambda;  % Update the wave number

% Run simulation with new parameters
intensity = run_simulation_fast(params, 0.0, 1.0, true);
```

**Option B: Edit the defaults file (makes change permanent)**

1. Open with text editor: `create_simulation_params.m`
2. Find the line that says something like: `params.lambda = 633e-9;`
3. Change the number to what you want (e.g., `params.lambda = 532e-9;`)
4. Save the file
5. Now any future time you run `create_simulation_params()`, it will use your new default

### Changing Grid Resolution (for faster/sharper simulation)

The grid resolution determines how detailed your result image is:
- **Higher resolution** (e.g., `grid_points_z = 2000`) = sharper, more detailed image, but much slower
- **Lower resolution** (e.g., `grid_points_z = 200`) = faster, but blurrier image

To change resolution:

1. Open `create_simulation_params.m`
2. Find the lines:
   ```matlab
   params.grid_points_z = 512;
   params.grid_points_r = 256;
   ```
3. Try smaller values to speed up: `grid_points_z = 256; grid_points_r = 128;`
4. Or try larger values for sharper images: `grid_points_z = 1024; grid_points_r = 512;`

### Important: Keep `k` in Sync with `lambda`

If you change the wavelength, you MUST update the wave number in the same edit:

```matlab
params.lambda = 532e-9;  % Change wavelength
params.k = 2*pi / params.lambda;  % Always do this when changing lambda
```

The variable `k` is the wave number and it's used in the diffraction integral calculation. If you forget to update it, your results will be wrong.

## Working with Masks (Ring Apertures)

By default, the aperture is a full circle. You can make it a ring (circular aperture with the center blocked) by using mask parameters.

The mask parameters are passed to the simulation functions:

```matlab
% Full circle (default)
intensity = run_simulation_fast(params, 0.0, 1.0, true);

% Ring: block the center (inner radius 0.3, outer radius 1.0)
intensity = run_simulation_fast(params, 0.3, 1.0, true);

% Smaller ring: block center and edges
intensity = run_simulation_fast(params, 0.2, 0.8, true);
```

The values are in normalized units (0 to 1), where:
- 0.0 = center of aperture
- 1.0 = edge of aperture

For a working example, look at `example_custom_params.m` which shows different mask configurations.

## Running a Parameter Sweep (Comparing Many Simulations)

A "sweep" means running many simulations while changing one parameter each time. For example, you might want to see how the focal spot changes as you vary the wavelength from 400 nm to 700 nm in small steps.

### Step 1: Open and Configure the Sweep File

Open the file: `run_sweep.m`

In the beginning of the file, you'll see a section labeled "SWEEP CONFIGURATION". This is where you control what parameter to vary:

```matlab
%% SWEEP CONFIGURATION
sweep_param = 'lambda';        % Change this to sweep a different parameter
sweep_start = 400e-9;          % Start wavelength: 400 nm (UV/blue light)
sweep_end = 700e-9;            % End wavelength: 700 nm (red light)
sweep_increment = 50e-9;       % Step size: 50 nm
```

This example sweeps the wavelength from blue (400 nm) to red (700 nm) in 50 nm steps.

**What other parameters can you sweep?** Any parameter in the `params` struct, such as:
- `'lambda'`: Change the wavelength
- `'f'`: Change the focal length
- `'a'`: Change the aperture size
- `'w0'`: Change the beam width
- `'n'`: Change the lens material
- `'grid_points_z'`: Change the resolution (though this is less physically interesting)

**Example 1: Sweep focal length**
```matlab
sweep_param = 'f';
sweep_start = 10e-3;      % Start at 10 mm
sweep_end = 100e-3;       % End at 100 mm
sweep_increment = 10e-3;  % Step by 10 mm each time
```

**Example 2: Sweep aperture size**
```matlab
sweep_param = 'a';
sweep_start = 0.5e-3;     % Start at 0.5 mm
sweep_end = 5e-3;         % End at 5 mm
sweep_increment = 0.5e-3; % Step by 0.5 mm
```

### Step 2: Configure Simulation Settings

Below the sweep configuration, there's a "SIMULATION SETTINGS" section. You can change:
- Which simulation method to use (`run_simulation_fast` or `run_simulation_gpu_ultra`)
- Grid resolution (for faster/slower sweeps)
- Mask configuration

### Step 3: Run the Sweep

From the MATLAB command window, type:

```matlab
run_sweep
```

This will run many simulations (one for each step in your sweep range) and save the results to a file in the `Outputs/` folder. This can take a few minutes to a few hours depending on:
- How many steps are in your sweep
- Which simulation method you're using
- Your computer's speed
- How high the grid resolution is

**Console output**: You'll see the progress printed to the console, like "Computing wavelength 450 nm..." as it works.

### Step 4: View the Results with the Sweep Viewer

After the sweep completes, you can visualize all the results:

```matlab
sweep_viewer
```

This will:
1. Ask you to select which sweep file you want to view
2. Load all the simulation results from that file
3. Show you an interactive viewer where you can:
   - Use a slider to move through different parameter values
   - Automatically play through all values
   - See plots of how the intensity changes with the parameter
   - See 1D cuts through the data (cross-sections)

This is a great way to understand how your results change as you vary a parameter.

## Benchmarking: Compare Speed of Different Methods

Want to know which simulation method is fastest on YOUR computer? Run the benchmarking script:

```matlab
benchmark_simulations
```

This will:
1. Run the same simulation using different methods (CPU, GPU, etc.)
2. Measure how long each method takes
3. Tell you which is fastest and by how much

This is useful when deciding which method to use for your own simulations. Typically:
- GPU methods are 10-100x faster than CPU methods
- But this depends on your specific computer and GPU

## Working Examples

Open the file `example_custom_params.m` to see working code examples of:
- Using default parameters
- Changing the wavelength
- Changing focal length and aperture
- Creating ring-shaped apertures (annular masks)
- Using high-resolution grids for publication-quality images

Copy and paste the examples into your MATLAB command window, or modify them to try your own combinations.

## Output Files and Storage

When you run a parameter sweep, results are saved to the `Outputs/` folder as `.mat` files (MATLAB data files). These files contain:
- All simulation results from the sweep
- The parameters used for each simulation
- Metadata about when the sweep was run

These files are ignored by git (see `.gitignore`), so they won't be backed up if you use git. Keep local backups if you want to preserve important results.

## Troubleshooting Guide

### "Out of Memory" Error

**Problem**: You get an error like "CUDA out of memory" or "out of memory" when running a simulation.

**Solutions**:
1. **Reduce grid resolution** - Open `create_simulation_params.m` and reduce these values:
   ```matlab
   params.grid_points_z = 256;  % was 512, now smaller
   params.grid_points_r = 128;  % was 256, now smaller
   ```
   This makes simulations run faster and use less memory, but the result image will be less detailed.

2. **Reduce integration points** - Also in `create_simulation_params.m`:
   ```matlab
   params.Nrho = 512;  % was 1024, integration points reduced
   ```

3. **Use CPU version instead of GPU** - If using `run_simulation_gpu_ultra`:
   ```matlab
   intensity = run_simulation_fast(params, 0.0, 1.0, true);  % Use CPU instead
   ```

4. **Restart MATLAB** - Sometimes MATLAB's memory gets fragmented. Close and reopen MATLAB.

### Simulation Results Look Wrong or Unexpected

**Problem**: The intensity pattern looks strange, doesn't match expectations, or changes unexpectedly when you change a parameter.

**Most common cause**: You changed `lambda` but forgot to update `k`.

**Solution**: After changing wavelength, ALWAYS do this:
```matlab
params.lambda = 532e-9;               % Your new wavelength
params.k = 2*pi / params.lambda;     % UPDATE K! This is essential!
```

**Other things to check**:
1. Did you change `grid_points_z` or `grid_points_r` very small? Small grid = blurry result.
2. Did you use unusual mask values (like mask_min = 0.8)? This creates unusual ring patterns, which is expected.
3. Are you looking at the right part of the output? Try adjusting the plot limits.

### Simulation is Very Slow

**Problem**: The simulation is taking too long (more than a few minutes for a single run).

**Solution 1: Use a faster method**
```matlab
% If you have a GPU, this is MUCH faster
intensity = run_simulation_gpu_ultra(params, 0.0, 1.0, true);

% If no GPU, at least use the parallel version
intensity = run_simulation_fast(params, 0.0, 1.0, true);
```

**Solution 2: Reduce resolution**
```matlab
params.grid_points_z = 256;   % Reduce from 512
params.grid_points_r = 128;   % Reduce from 256
```

**Solution 3: Reduce integration accuracy**
```matlab
params.Nrho = 512;  % Reduce from 1024
```
This computes the integral using fewer sample points. Results may be slightly less accurate but will run much faster.

### Can't Find a File

**Problem**: You get an error like "File not found" or "Undefined function".

**Solution**: Make sure all the project files are in the same folder. Check that you have:
- `create_simulation_params.m`
- `run_simulation_fast.m` (or the GPU version you're using)
- `sweep_viewer.m`
- And others from the file list at the top of the README

All files should be in the same directory when you run them.

### Parameter Changes Don't Have Any Effect

**Problem**: You change a parameter but the result looks identical.

**Possible causes**:
1. **Grid is too coarse**: You're changing a parameter that has subtle effects, but your grid is too coarse to see them. Try increasing `grid_points_z` and `grid_points_r`.
2. **It's actually the same effect, just hard to see**: Try making a more dramatic change. Instead of `lambda = 630e-9` to `632e-9`, try `600e-9` to `700e-9` (much bigger difference).
3. **The result is normalized**: If you normalized the result (`true` in the simulation call), it's hard to see small intensity changes. Try running with `false`:
   ```matlab
   intensity = run_simulation_fast(params, 0.0, 1.0, false);  % Raw values, not normalized
   ```

### Sweep Viewer Not Responding

**Problem**: You run `sweep_viewer` and nothing happens, or it's very slow.

**Solutions**:
1. **Check if the file exists**: Make sure you actually ran and completed a sweep before trying to view it.
2. **Reduce file size**: In your next sweep, use lower resolution to create a smaller result file:
   ```matlab
   params.grid_points_z = 256;
   params.grid_points_r = 128;
   ```
3. **Restart MATLAB**: Sometimes the viewer gets stuck. Close and reopen MATLAB.

## For More Details

- **Understanding the math**: See `MATH_GUIDE.md`
- **All parameter descriptions**: See `PARAMETERS_GUIDE.m`
- **What each variable does in the code**: Look at comments in `create_simulation_params.m`
