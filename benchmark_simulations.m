% Benchmark all simulation methods
clear; clc;

fprintf('=== Simulation Speed Benchmark ===\n\n');

% Create simulation parameters using helper function
params = create_simulation_params();

% Can customize if needed:
% params.grid_points_z = 200;  % Reduce for faster benchmark
% params.grid_points_r = 200;

% Test parameters
mask_min = 0.0;
mask_max = 1.0;

% Check GPU availability
has_gpu = gpuDeviceCount > 0;
if has_gpu
    gpu = gpuDevice();
    fprintf('GPU detected: %s (%.1f GB available)\n\n', gpu.Name, gpu.AvailableMemory/1e9);
else
    fprintf('No GPU detected. GPU tests will be skipped.\n\n');
end

% Test original (if you want - typically very slow)
run_original = false;
if run_original
    fprintf('1. Testing ORIGINAL (integral-based, single-threaded)...\n');
    tic;
    I1 = run_simulation(params, mask_min, mask_max, false);
    t1 = toc;
    fprintf('   Time: %.2f seconds\n\n', t1);
else
    fprintf('1. ORIGINAL version skipped (too slow)\n\n');
    t1 = nan;
end

% Test optimized (vectorized)
fprintf('2. Testing OPTIMIZED (vectorized integration)...\n');
tic;
I2 = run_simulation_optimized(params, mask_min, mask_max, true);
t2 = toc;
fprintf('   Time: %.2f seconds\n\n', t2);

% Test fast (vectorized + parallel)
fprintf('3. Testing FAST (vectorized + parallel CPU)...\n');
tic;
I3 = run_simulation_fast(params, mask_min, mask_max, true);
t3 = toc;
fprintf('   Time: %.2f seconds\n\n', t3);

% Test GPU versions
if has_gpu
    fprintf('4. Testing GPU (chunked processing)...\n');
    try
        tic;
        I4 = run_simulation_gpu(params, mask_min, mask_max, true);
        t4 = toc;
        fprintf('   Time: %.2f seconds\n\n', t4);
    catch ME
        fprintf('   ERROR: %s\n\n', ME.message);
        t4 = nan;
    end
    
    fprintf('5. Testing GPU ULTRA (full 3D parallelization)...\n');
    try
        tic;
        I5 = run_simulation_gpu_ultra(params, mask_min, mask_max, true);
        t5 = toc;
        fprintf('   Time: %.2f seconds\n\n', t5);
    catch ME
        fprintf('   ERROR: %s\n   (Likely out of GPU memory - try reducing grid size)\n\n', ME.message);
        t5 = nan;
    end
else
    t4 = nan;
    t5 = nan;
end

% Summary
fprintf('=== SUMMARY ===\n');
if run_original
    fprintf('Original:    %.2f s  (1.0x)\n', t1);
end
fprintf('Optimized:   %.2f s  ', t2);
if ~isnan(t1)
    fprintf('(%.1fx faster)', t1/t2);
end
fprintf('\n');

fprintf('Fast (CPU):  %.2f s  ', t3);
if ~isnan(t1)
    fprintf('(%.1fx faster)', t1/t3);
end
fprintf('\n');

if has_gpu && ~isnan(t4)
    fprintf('GPU:         %.2f s  ', t4);
    if ~isnan(t1)
        fprintf('(%.1fx faster)', t1/t4);
    end
    fprintf('\n');
end

if has_gpu && ~isnan(t5)
    fprintf('GPU Ultra:   %.2f s  ', t5);
    if ~isnan(t1)
        fprintf('(%.1fx faster)', t1/t5);
    end
    fprintf('\n');
end

fprintf('\nRecommendation: ');
if has_gpu && ~isnan(t5) && t5 < t3
    fprintf('Use run_simulation_gpu_ultra() for maximum speed!\n');
elseif has_gpu && ~isnan(t4) && t4 < t3
    fprintf('Use run_simulation_gpu() for best speed!\n');
else
    fprintf('Use run_simulation_fast() for best speed!\n');
end
