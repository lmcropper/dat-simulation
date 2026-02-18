function sweep_viewer()
% Fast Sweep Viewer - Rebuilt using proven fast architecture from debug version
% Key difference: Pre-create all UI elements once, just update data with set()

close all;

% ========== FILE SELECTION & LOADING ==========
[filename, pathname] = uigetfile('*.mat', 'Select sweep results file');
if isequal(filename, 0)
    return;
end

file_path = fullfile(pathname, filename);
file_info = dir(file_path);
file_size_mb = file_info.bytes / 1e6;

% Check format
try
    h5info(file_path);
    is_v73 = true;
catch
    is_v73 = false;
end

% Load strategy - ALWAYS FULL-LOAD for single viewer (fastest for interactive use)
use_lazy_loading = false;
fprintf('Loading entire dataset (%.1f MB) into memory...\n', file_size_mb);
tic;
data = load(file_path);
t_load = toc;

intensity_sweep_full = data.intensity_sweep;
fprintf('  -> Loaded in %.2f sec, array size: %.1f MB\n', t_load, numel(intensity_sweep_full)*8/1e6);
fprintf('  -> FULL-LOAD MODE: All slices in RAM, ~1-2ms per frame access\n\n');

sweep_param_values = data.sweep_param_values;
sweep_param = data.sweep_param;
u_values = data.u_values;
v_values = data.v_values;
z_scale = data.z_scale;
r_scale = data.r_scale;
matData = [];

% Get dimensions
num_r = size(intensity_sweep_full, 1);
num_u = size(intensity_sweep_full, 2);
num_slices = size(intensity_sweep_full, 3);

fprintf('Array: %d × %d × %d\n', num_r, num_u, num_slices);
fprintf('Parameter: %s (%.3f to %.3f)\n\n', sweep_param, sweep_param_values(1), sweep_param_values(end));

% Compute axis vectors
z_axis = u_values * z_scale * 1e6;
r_axis = v_values * r_scale * 1e6;
[~, v_zero_idx] = min(abs(v_values));
[~, u_zero_idx] = min(abs(u_values));

% Pre-compute gamma LUT with high resolution for smooth display
gamma_lut = struct();
gamma_lut.gamma_val = 0.5;
gamma_lut.lut = linspace(0, 1, 4096).^0.5;  % 4096 levels for smooth gamma

% ========== BUILD UI (Pre-create all axes once) ==========
fig = figure('Name', sprintf('Sweep Viewer (Fast): %s', sweep_param), ...
    'NumberTitle', 'off', ...
    'Position', [100 100 1600 600], ...
    'Renderer', 'painters', ...
    'DoubleBuffer', 'on');

% Main image - left side, takes upper portion
ax_main = axes('Position', [0.08 0.40 0.55 0.5]);
img_handle = imagesc(ax_main, z_axis, r_axis, zeros(num_r, num_u));
colormap(ax_main, hot);
cb_handle = colorbar(ax_main);
xlabel(ax_main, 'Axial Distance z (μm)');
ylabel(ax_main, 'Radial Distance r (μm)');
title_main = title(ax_main, '');
hold(ax_main, 'on');
line_z_indicator = line(ax_main, [z_axis(u_zero_idx) z_axis(u_zero_idx)], [r_axis(1) r_axis(end)], ...
    'Color', [0 0.7 1], 'LineWidth', 2, 'LineStyle', '--');
hold(ax_main, 'off');

% Secondary plots - right side, stacked vertically
% Axial plot - top right
ax_axial = axes('Position', [0.67 0.65 0.3 0.25], 'Visible', 'off');
line_axial = plot(ax_axial, z_axis, zeros(1, num_u), 'LineWidth', 2, 'Color', [1 0.5 0]);
xlabel(ax_axial, 'Axial Distance z (μm)');
ylabel(ax_axial, 'Normalized Intensity');
title(ax_axial, 'Intensity Along Optical Axis (r=0)');
grid(ax_axial, 'on');

% Radial plot - bottom right
ax_radial = axes('Position', [0.67 0.3 0.3 0.25], 'Visible', 'off');
line_radial = plot(ax_radial, r_axis, zeros(num_r, 1), 'LineWidth', 2, 'Color', [0 0.7 1]);
xlabel(ax_radial, 'Radial Distance r (μm)');
ylabel(ax_radial, 'Normalized Intensity');
title_radial = title(ax_radial, 'Radial Profile');
grid(ax_radial, 'on');

% Slider for radial plot z-position
slider_radial_step = min(1, [1/(num_u-1), 10/(num_u-1)]);
slider_radial_h = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', num_u, 'Value', u_zero_idx, ...
    'SliderStep', slider_radial_step, ...
    'Position', [1100 100 300 20], ...
    'Callback', @update_plot);

% Label for radial slider
uicontrol('Style', 'text', 'String', 'Radial z-position:', ...
    'Position', [1100 50 100 20], 'HorizontalAlignment', 'right');

text_radial_z = uicontrol('Style', 'text', 'String', '', ...
    'Position', [1200 50 80 20]);

% Slider - bottom, full width
slider_step = min(1, [1/(num_slices-1), 10/(num_slices-1)]);
slider_h = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', num_slices, 'Value', 1, ...
    'SliderStep', slider_step, ...
    'Position', [100 165 650 20], ...
    'Callback', @update_plot);

% Navigation buttons - left side
btn_prev = uicontrol('Style', 'pushbutton', 'String', '< Previous', ...
    'Position', [100 125 80 30], ...
    'Callback', @(src,evt) move_slider(-1));

btn_next = uicontrol('Style', 'pushbutton', 'String', 'Next >', ...
    'Position', [190 125 80 30], ...
    'Callback', @(src,evt) move_slider(1));

btn_play = uicontrol('Style', 'pushbutton', 'String', 'Play', ...
    'Position', [300 125 60 30], ...
    'Callback', @play_animation);

% Gamma control - center
uicontrol('Style', 'text', 'String', 'Gamma:', 'Position', [450 130 40 20]);
gamma_edit = uicontrol('Style', 'edit', 'String', '0.5', ...
    'Position', [400 130 40 20], 'Callback', @update_plot);

% Display checkboxes - right side
cb_axial = uicontrol('Style', 'checkbox', 'String', 'Show axial plot', ...
    'Position', [800 180 150 20], 'Value', 1, ...
    'Callback', @update_plot);

cb_radial = uicontrol('Style', 'checkbox', 'String', 'Show radial plot', ...
    'Position', [800 160 150 20], 'Value', 1, ...
    'Callback', @update_plot);

cb_normalize = uicontrol('Style', 'checkbox', 'String', 'Normalize intensity', ...
    'Position', [800 140 150 20], 'Value', 1, ...
    'Callback', @update_plot);

cb_log_scale = uicontrol('Style', 'checkbox', 'String', 'Log scale', ...
    'Position', [800 120 100 20], 'Value', 0, ...
    'Callback', @update_plot);

% Info display - bottom
info_text = uicontrol('Style', 'text', ...
    'Position', [100 20 900 90], ...
    'FontSize', 10,'HorizontalAlignment', 'left');

is_playing = false;

% ========== CALLBACKS ==========
    function intensity = load_slice(idx)
        % Direct memory access (all data loaded upfront)
        intensity = intensity_sweep_full(:, :, idx);
    end

    function update_plot(src, evt)
        tic_main = tic;
        
        idx = round(get(slider_h, 'Value'));
        idx = max(1, min(idx, num_slices));
        set(slider_h, 'Value', idx);
        
        % Load slice
        tic; intensity_current = load_slice(idx); t_load = toc;
        param_value = sweep_param_values(idx);
        
        % Gamma correction
        try
            gamma = str2double(get(gamma_edit, 'String'));
            if isnan(gamma) || gamma <= 0
                gamma = 0.5;
                set(gamma_edit, 'String', '0.5');
            end
        catch
            gamma = 0.5;
        end
        
        % Update LUT if gamma changed
        if gamma ~= gamma_lut.gamma_val
            gamma_lut.gamma_val = gamma;
            gamma_lut.lut = linspace(0, 1, 4096).^gamma;
        end
        
        % Apply log scale and gamma via LUT (with optional normalization)
        tic;
        if get(cb_normalize, 'Value')
            % Normalize to [0, 1] before gamma correction
            max_val = max(intensity_current(:));
            intensity_norm = intensity_current / max_val;
        else
            % Use raw intensity values (assumes already in reasonable range)
            max_val = max(intensity_current(:));
            intensity_norm = intensity_current / max(max_val, 1e-10);  % Prevent division by zero
        end
        
        % Apply log scale if enabled
        if get(cb_log_scale, 'Value')
            % Log scale with small offset to avoid log(0)
            eps = 1e-6;
            intensity_norm = log10(intensity_norm + eps);
            % Normalize log values to [0, 1] range
            min_log = log10(eps);
            max_log = log10(1 + eps);
            intensity_norm = (intensity_norm - min_log) / (max_log - min_log);
            intensity_norm = max(0, min(intensity_norm, 1));  % Clamp to [0, 1]
        end
        
        intensity_idx = max(1, min(4096, round(intensity_norm * 4095) + 1));
        intensity_gamma = reshape(gamma_lut.lut(intensity_idx), size(intensity_current));
        t_gamma = toc;
        
        % Update main image
        tic; set(img_handle, 'CData', intensity_gamma); t_img = toc;
        set(title_main, 'String', sprintf('%s = %.3f', sweep_param, param_value));
        
        % Update colorbar label and ticks based on scale
        if get(cb_log_scale, 'Value')
            cb_handle.Label.String = 'Log Scale';
            % Set custom ticks showing orders of magnitude
            cb_handle.Ticks = [0 0.2 0.4 0.6 0.8 1];
            % Map normalized colorbar positions back to log values
            % 0->-6, 1->0 (approximately)
            eps = 1e-6;
            min_log = log10(eps);
            max_log = log10(1 + eps);
            tick_vals = cb_handle.Ticks * (max_log - min_log) + min_log;
            tick_labels = arrayfun(@(x) sprintf('10^%.0f', x), tick_vals, 'UniformOutput', false);
            cb_handle.TickLabels = tick_labels;
        else
            cb_handle.Label.String = 'Linear Scale';
            % Reset to automatic ticks
            cb_handle.Ticks = linspace(0, 1, 6);
            cb_handle.TickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(0, 1, 6), 'UniformOutput', false);
        end
        
        % Update secondary plots visibility and data
        tic;
        if get(cb_axial, 'Value')
            set(ax_axial, 'Visible', 'on');
            intensity_axis = intensity_gamma(v_zero_idx, :);
            set(line_axial, 'YData', intensity_axis);
        else
            set(ax_axial, 'Visible', 'off');
        end
        
        if get(cb_radial, 'Value')
            set(ax_radial, 'Visible', 'on');
            u_idx_radial = round(get(slider_radial_h, 'Value'));
            u_idx_radial = max(1, min(u_idx_radial, num_u));
            intensity_rad = intensity_gamma(:, u_idx_radial);
            z_value = z_axis(u_idx_radial);
            set(line_radial, 'YData', intensity_rad);
            set(title_radial, 'String', sprintf('Radial Profile at z = %.2f μm', z_value));
            set(text_radial_z, 'String', sprintf('%.2f μm', z_value));
            % Update vertical line indicator on main plot
            set(line_z_indicator, 'XData', [z_value z_value]);
        else
            set(ax_radial, 'Visible', 'off');
        end
        t_secondary = toc;
        
        % Update info text
        tic;
        norm_status = '';
        if get(cb_normalize, 'Value')
            norm_status = 'normalized';
        else
            norm_status = 'raw';
        end
        set(info_text, 'String', ...
            sprintf(['Slice: %d / %d\n' ...
                     '%s = %.4f\n' ...
                     'Grid: %d × %d | Intensity: %s\n' ...
                     'Load: %.1f ms | Gamma: %.1f ms | Image: %.1f ms | Secondary: %.1f ms'], ...
                idx, num_slices, ...
                sweep_param, param_value, ...
                num_r, num_u, norm_status, ...
                t_load*1000, t_gamma*1000, t_img*1000, t_secondary*1000));
        t_info = toc;
        
        % Optimized rendering
        tic; drawnow limitrate; t_draw = toc;
        
        total = toc(tic_main);
        fprintf('Total: %.0f ms (load:%.0f + gamma:%.0f + img:%.0f + second:%.0f + info:%.0f + draw:%.0f)\n', ...
            total*1000, t_load*1000, t_gamma*1000, t_img*1000, t_secondary*1000, t_info*1000, t_draw*1000);
    end

    function move_slider(direction)
        current_val = get(slider_h, 'Value');
        new_val = current_val + direction;
        new_val = max(1, min(new_val, num_slices));
        set(slider_h, 'Value', new_val);
        update_plot();
    end

    function play_animation(src, evt)
        is_playing = ~is_playing;
        
        if is_playing
            set(btn_play, 'String', 'Pause');
            
            while is_playing && get(slider_h, 'Value') < num_slices
                move_slider(1);
                pause(0.05);
            end
            
            is_playing = false;
            set(btn_play, 'String', 'Play');
        else
            set(btn_play, 'String', 'Play');
        end
    end

% ========== INITIALIZE ==========
fprintf('Viewer ready!\n');
fprintf('Full-load mode: all slices in memory for instant navigation (~1-2ms per frame).\n');
fprintf('Tip: Keep secondary plot checkboxes OFF for maximum speed.\n\n');

update_plot();
end
