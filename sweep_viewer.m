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
line_r_indicator = line(ax_main, [z_axis(1) z_axis(end)], [r_axis(v_zero_idx) r_axis(v_zero_idx)], ...
    'Color', [1 0.5 0], 'LineWidth', 2, 'LineStyle', '--');
hold(ax_main, 'off');

% Secondary plots - right side, stacked vertically
% Axial plot - top right
ax_axial = axes('Position', [0.67 0.65 0.3 0.25], 'Visible', 'off');
line_axial = plot(ax_axial, z_axis, zeros(1, num_u), 'LineWidth', 2, 'Color', [1 0.5 0]);
xlabel(ax_axial, 'Axial Distance z (μm)');
ylabel(ax_axial, 'Normalized Intensity');
title_axial = title(ax_axial, 'Intensity Along Optical Axis (r=0)');
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
uicontrol('Style', 'text', 'String', 'Radial plot z-pos:', ...
    'Position', [1100 78 100 20], 'HorizontalAlignment', 'right');

text_radial_z = uicontrol('Style', 'text', 'String', '', ...
    'Position', [1200 78 80 20]);

% Slider for axial plot r-position
slider_axial_step = min(1, [1/(num_r-1), 10/(num_r-1)]);
slider_axial_h = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', num_r, 'Value', v_zero_idx, ...
    'SliderStep', slider_axial_step, ...
    'Position', [1100 50 300 20], ...
    'Callback', @update_plot);

% Label for axial slider
uicontrol('Style', 'text', 'String', 'Axial plot r-pos:', ...
    'Position', [1100 28 100 20], 'HorizontalAlignment', 'right');

text_axial_r = uicontrol('Style', 'text', 'String', '', ...
    'Position', [1200 28 80 20]);

% View mode dropdown (extensible for future analysis methods)
uicontrol('Style', 'text', 'String', 'View:', ...
    'Position', [960 180 40 20], 'HorizontalAlignment', 'right');
view_mode_popup = uicontrol('Style', 'popupmenu', ...
    'String', {'Standard plots', '2D autocorrelation', '2D FFT'}, ...
    'Position', [1005 180 170 20], ...
    'Callback', @on_view_mode_change);

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

% Autocorrelation plots area setup (will be created on demand)
autocorr_axes_created = false;
autocorr_handles = struct();
autocorr_handles.cb_2d = [];
autocorr_handles.cb_fft_2d = [];
intensity_gamma = [];  % Store current gamma-corrected slice for autocorr use

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
        if get(view_mode_popup, 'Value') == 1
            set(cb_handle, 'Visible', 'on');
        else
            set(cb_handle, 'Visible', 'off');
        end
        
        % Update secondary plots visibility and data
        tic;
        
        % Always read slider values and update text readouts (regardless of view mode)
        v_idx_axial = round(get(slider_axial_h, 'Value'));
        v_idx_axial = max(1, min(v_idx_axial, num_r));
        r_value = r_axis(v_idx_axial);
        set(text_axial_r, 'String', sprintf('%.2f μm', r_value));
        
        u_idx_radial = round(get(slider_radial_h, 'Value'));
        u_idx_radial = max(1, min(u_idx_radial, num_u));
        z_value = z_axis(u_idx_radial);
        set(text_radial_z, 'String', sprintf('%.2f μm', z_value));
        
        % Update plots based on view mode
        if get(view_mode_popup, 'Value') == 1
            if get(cb_axial, 'Value')
                set(ax_axial, 'Visible', 'on');
                intensity_axis = intensity_gamma(v_idx_axial, :);
                set(line_axial, 'YData', intensity_axis);
                set(title_axial, 'String', sprintf('Axial Profile at r = %.2f μm', r_value));
                % Update horizontal line indicator on main plot
                set(line_r_indicator, 'YData', [r_value r_value]);
            else
                set(ax_axial, 'Visible', 'off');
            end
            
            if get(cb_radial, 'Value')
                set(ax_radial, 'Visible', 'on');
                intensity_rad = intensity_gamma(:, u_idx_radial);
                set(line_radial, 'YData', intensity_rad);
                set(title_radial, 'String', sprintf('Radial Profile at z = %.2f μm', z_value));
                % Update vertical line indicator on main plot
                set(line_z_indicator, 'XData', [z_value z_value]);
            else
                set(ax_radial, 'Visible', 'off');
            end
        else
            set(ax_axial, 'Visible', 'off');
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
        
        % Update analysis plots only when the view mode is selected
        if get(view_mode_popup, 'Value') == 2
            update_autocorr_plots(intensity_gamma);
        elseif get(view_mode_popup, 'Value') == 3
            update_fft_plots(intensity_gamma);
        end
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

    function on_view_mode_change(src, evt)
        mode_idx = get(view_mode_popup, 'Value');
        apply_view_mode(mode_idx);
    end

    function apply_view_mode(mode_idx)
        if mode_idx == 1
            set(ax_main, 'Visible', 'on');
            set(get(ax_main, 'Children'), 'Visible', 'on');
            set(cb_handle, 'Visible', 'on');
            dcm = datacursormode(fig);
            set(dcm, 'UpdateFcn', []);
            if get(cb_axial, 'Value')
                set(ax_axial, 'Visible', 'on');
                set(get(ax_axial, 'Children'), 'Visible', 'on');
            else
                set(ax_axial, 'Visible', 'off');
                set(get(ax_axial, 'Children'), 'Visible', 'off');
            end
            if get(cb_radial, 'Value')
                set(ax_radial, 'Visible', 'on');
                set(get(ax_radial, 'Children'), 'Visible', 'on');
            else
                set(ax_radial, 'Visible', 'off');
                set(get(ax_radial, 'Children'), 'Visible', 'off');
            end
            if autocorr_axes_created
                set(autocorr_handles.ax_2d, 'Visible', 'off');
                set(autocorr_handles.ax_z, 'Visible', 'off');
                set(autocorr_handles.ax_r, 'Visible', 'off');
                set(get(autocorr_handles.ax_2d, 'Children'), 'Visible', 'off');
                set(get(autocorr_handles.ax_z, 'Children'), 'Visible', 'off');
                set(get(autocorr_handles.ax_r, 'Children'), 'Visible', 'off');
                if ~isempty(autocorr_handles.cb_2d) && isvalid(autocorr_handles.cb_2d)
                    set(autocorr_handles.cb_2d, 'Visible', 'off');
                end
                if ~isempty(autocorr_handles.cb_fft_2d) && isvalid(autocorr_handles.cb_fft_2d)
                    set(autocorr_handles.cb_fft_2d, 'Visible', 'off');
                end
            end
            drawnow();
        elseif mode_idx == 2
            ensure_autocorr_axes();
            set(ax_main, 'Visible', 'off');
            set(get(ax_main, 'Children'), 'Visible', 'off');
            set(cb_handle, 'Visible', 'off');
            set(ax_axial, 'Visible', 'off');
            set(get(ax_axial, 'Children'), 'Visible', 'off');
            set(ax_radial, 'Visible', 'off');
            set(get(ax_radial, 'Children'), 'Visible', 'off');
            set(autocorr_handles.ax_2d, 'Visible', 'on');
            set(autocorr_handles.ax_z, 'Visible', 'on');
            set(autocorr_handles.ax_r, 'Visible', 'on');
            if ~isempty(autocorr_handles.cb_2d) && isvalid(autocorr_handles.cb_2d)
                set(autocorr_handles.cb_2d, 'Visible', 'on');
            end
            if ~isempty(autocorr_handles.cb_fft_2d) && isvalid(autocorr_handles.cb_fft_2d)
                set(autocorr_handles.cb_fft_2d, 'Visible', 'off');
            end
            uistack(autocorr_handles.ax_2d, 'top');
            uistack(autocorr_handles.ax_z, 'top');
            uistack(autocorr_handles.ax_r, 'top');
            update_autocorr_plots(intensity_gamma);
            dcm = datacursormode(fig);
            set(dcm, 'UpdateFcn', []);
        else
            ensure_autocorr_axes();
            set(ax_main, 'Visible', 'off');
            set(get(ax_main, 'Children'), 'Visible', 'off');
            set(cb_handle, 'Visible', 'off');
            set(ax_axial, 'Visible', 'off');
            set(get(ax_axial, 'Children'), 'Visible', 'off');
            set(ax_radial, 'Visible', 'off');
            set(get(ax_radial, 'Children'), 'Visible', 'off');
            set(autocorr_handles.ax_2d, 'Visible', 'on');
            set(autocorr_handles.ax_z, 'Visible', 'on');
            set(autocorr_handles.ax_r, 'Visible', 'on');
            if ~isempty(autocorr_handles.cb_2d) && isvalid(autocorr_handles.cb_2d)
                set(autocorr_handles.cb_2d, 'Visible', 'off');
            end
            if ~isempty(autocorr_handles.cb_fft_2d) && isvalid(autocorr_handles.cb_fft_2d)
                set(autocorr_handles.cb_fft_2d, 'Visible', 'on');
            end
            uistack(autocorr_handles.ax_2d, 'top');
            uistack(autocorr_handles.ax_z, 'top');
            uistack(autocorr_handles.ax_r, 'top');
            update_fft_plots(intensity_gamma);
            dcm = datacursormode(fig);
            set(dcm, 'UpdateFcn', @fft_datatip_txt);
        end
    end

    function ensure_autocorr_axes()
        if autocorr_axes_created
            return;
        end
        autocorr_handles.ax_2d = axes('Parent', fig, 'Position', [0.08 0.40 0.55 0.5]);
        set(autocorr_handles.ax_2d, 'XDir', 'normal', 'YDir', 'normal', 'Layer', 'top');
        autocorr_handles.ax_z = axes('Parent', fig, 'Position', [0.67 0.65 0.3 0.25]);
        set(autocorr_handles.ax_z, 'XDir', 'normal', 'YDir', 'normal', 'Layer', 'top');
        autocorr_handles.ax_r = axes('Parent', fig, 'Position', [0.67 0.3 0.3 0.25]);
        set(autocorr_handles.ax_r, 'XDir', 'normal', 'YDir', 'normal', 'Layer', 'top');
        autocorr_axes_created = true;
    end

    function update_autocorr_plots(intensity_data)
        if ~autocorr_axes_created
            return;
        end
        
        visibility = get(autocorr_handles.ax_2d, 'Visible');
        if ~strcmp(visibility, 'on')
            return;
        end
        
        try
            % Use provided data (already gamma-corrected from update_plot)
            if isempty(intensity_data) || ~isnumeric(intensity_data)
                % Fallback: create synthetic periodic data
                ny = 50;
                nz = 50;
                [yy, zz] = meshgrid(1:nz, 1:ny);
                temp = sin(2*pi*yy/20) .* cos(2*pi*zz/20) + 0.5;
                intensity_data = temp;
            end
            
            [ny, nz] = size(intensity_data);
            
            % Handle edge cases
            if ny < 2 || nz < 2
                return;
            end
            
            % Normalize intensity data to [0, 1]
            intensity_min = min(intensity_data(:));
            intensity_max = max(intensity_data(:));
            if intensity_max > intensity_min
                intensity_norm = (intensity_data - intensity_min) / (intensity_max - intensity_min);
            else
                intensity_norm = ones(ny, nz) * 0.5;  % Constant data fallback
            end
            
            % Compute 2D autocorrelation using FFT with padding
            pad_y = ny;
            pad_z = nz;
            intensity_padded = pad_array(intensity_norm, [pad_y, pad_z], 0, 'post');
            
            % FFT-based autocorrelation computation
            fft_result = fft2(intensity_padded);
            power_spectrum = abs(fft_result) .^ 2;
            autocorr_full = ifft2(power_spectrum);
            autocorr_full = real(autocorr_full);
            
            % Extract relevant portion
            autocorr_2d = autocorr_full(1:ny, 1:nz);
            
            % Normalize by the zero-lag autocorrelation value
            max_val = max(autocorr_2d(:));
            if max_val > 1e-10
                autocorr_2d = autocorr_2d / max_val;
            end
            
            % Create centered view by mirroring
            if ny > 1 && nz > 1
                autocorr_2d_centered = [fliplr(autocorr_2d(:, 2:end)), autocorr_2d];
                autocorr_2d_centered = [flipud(autocorr_2d_centered(2:end, :)); autocorr_2d_centered];
            else
                autocorr_2d_centered = autocorr_2d;
            end
            
            % Compute physical shift axes
            if length(r_axis) > 1
                dr = r_axis(2) - r_axis(1);
            else
                dr = 1;
            end
            
            if length(z_axis) > 1
                dz = z_axis(2) - z_axis(1);
            else
                dz = 1;
            end
            
            n_r_centered = size(autocorr_2d_centered, 1);
            n_z_centered = size(autocorr_2d_centered, 2);
            
            % Center indices
            center_r_idx = (n_r_centered + 1) / 2;
            center_z_idx = (n_z_centered + 1) / 2;
            
            % Shift axes in physical units
            r_shifts = ((1:n_r_centered) - center_r_idx) * dr;
            z_shifts = ((1:n_z_centered) - center_z_idx) * dz;
            
            % ===== PLOT 1: 2D Autocorrelation =====
            cla(autocorr_handles.ax_2d);
            set(autocorr_handles.ax_2d, 'NextPlot', 'replacechildren');
            
            % Use imagesc with proper coordinate specification
            im_handle = imagesc(autocorr_handles.ax_2d, z_shifts, r_shifts, autocorr_2d_centered);
            
            % Set proper axis properties
            set(autocorr_handles.ax_2d, 'YDir', 'normal');
            axis(autocorr_handles.ax_2d, 'auto');
            
            % Add colormap and colorbar
            colormap(autocorr_handles.ax_2d, hot);
            if isempty(autocorr_handles.cb_2d) || ~isvalid(autocorr_handles.cb_2d)
                autocorr_handles.cb_2d = colorbar(autocorr_handles.ax_2d);
            end
            set(autocorr_handles.cb_2d, 'Visible', 'on');
            autocorr_handles.cb_2d.Label.String = 'Normalized Autocorr';
            
            % Labels and title
            xlabel(autocorr_handles.ax_2d, 'Delta z (μm)', 'FontSize', 10);
            ylabel(autocorr_handles.ax_2d, 'Delta r (μm)', 'FontSize', 10);
            title(autocorr_handles.ax_2d, '2D Autocorrelation', 'FontSize', 11, 'FontWeight', 'bold');
            grid(autocorr_handles.ax_2d, 'off');
            
            % ===== PLOT 2 & 3: Extract cross-sections =====
            center_r_int = round(center_r_idx);
            center_z_int = round(center_z_idx);
            
            % Ensure indices are valid
            center_r_int = max(1, min(center_r_int, n_r_centered));
            center_z_int = max(1, min(center_z_int, n_z_centered));
            
            % Use selected radial position for 1D autocorr along z (axial)
            v_idx_axial = round(get(slider_axial_h, 'Value'));
            v_idx_axial = max(1, min(v_idx_axial, num_r));
            z_slice = intensity_norm(v_idx_axial, :);
            z_pad = [z_slice, zeros(1, nz)];
            fft_z = fft(z_pad);
            ac_full_1d_z = ifft(abs(fft_z).^2);
            ac_1d_z = real(ac_full_1d_z(1:nz));
            if ac_1d_z(1) > 1e-10
                ac_1d_z = ac_1d_z / ac_1d_z(1);
            end
            autocorr_z_line = [fliplr(ac_1d_z(2:end)), ac_1d_z];
            z_shifts_1d = ((1:(2*nz-1)) - nz) * dz;
            
            % Use selected radial slice for 1D autocorr along r
            u_idx_radial = round(get(slider_radial_h, 'Value'));
            u_idx_radial = max(1, min(u_idx_radial, num_u));
            r_slice = intensity_norm(:, u_idx_radial);
            r_pad = [r_slice; zeros(ny, 1)];
            fft_r = fft(r_pad);
            ac_full_1d = ifft(abs(fft_r).^2);
            ac_1d = real(ac_full_1d(1:ny));
            if ac_1d(1) > 1e-10
                ac_1d = ac_1d / ac_1d(1);
            end
            autocorr_r_line = [flipud(ac_1d(2:end)); ac_1d];
            r_shifts_1d = ((1:(2*ny-1)) - ny) * dr;
            
            % ===== PLOT 2: Z-axis Cross-section =====
            cla(autocorr_handles.ax_z);
            set(autocorr_handles.ax_z, 'NextPlot', 'replacechildren');
            plot(autocorr_handles.ax_z, z_shifts_1d, autocorr_z_line, 'Color', [1 0.5 0], 'LineWidth', 2.5);
            hold(autocorr_handles.ax_z, 'off');
            
            set(autocorr_handles.ax_z, 'XGrid', 'on', 'YGrid', 'on');
            xlabel(autocorr_handles.ax_z, 'Δz (μm)', 'FontSize', 9);
            ylabel(autocorr_handles.ax_z, 'Correlation', 'FontSize', 9);
            title(autocorr_handles.ax_z, 'Z-Axis', 'FontSize', 10, 'FontWeight', 'bold');
            
            ylim(autocorr_handles.ax_z, [min(autocorr_z_line) - 0.05, max(autocorr_z_line) + 0.1]);
            xlim(autocorr_handles.ax_z, [min(z_shifts_1d) max(z_shifts_1d)]);
            
            % ===== PLOT 3: R-axis Cross-section =====
            cla(autocorr_handles.ax_r);
            set(autocorr_handles.ax_r, 'NextPlot', 'replacechildren');
            plot(autocorr_handles.ax_r, r_shifts_1d, autocorr_r_line, 'Color', [0 0.7 1], 'LineWidth', 2.5);
            hold(autocorr_handles.ax_r, 'off');
            
            set(autocorr_handles.ax_r, 'XGrid', 'on', 'YGrid', 'on');
            xlabel(autocorr_handles.ax_r, 'Δr (μm)', 'FontSize', 9);
            ylabel(autocorr_handles.ax_r, 'Correlation', 'FontSize', 9);
            title(autocorr_handles.ax_r, 'R-Axis', 'FontSize', 10, 'FontWeight', 'bold');
            
            ylim(autocorr_handles.ax_r, [min(autocorr_r_line) - 0.05, max(autocorr_r_line) + 0.1]);
            xlim(autocorr_handles.ax_r, [min(r_shifts_1d) max(r_shifts_1d)]);
            
            % Force redraw
            drawnow();
            
        catch
        end
    end

    function update_fft_plots(intensity_data)
        if ~autocorr_axes_created
            return;
        end
        
        visibility = get(autocorr_handles.ax_2d, 'Visible');
        if ~strcmp(visibility, 'on')
            return;
        end
        
        try
            if isempty(intensity_data) || ~isnumeric(intensity_data)
                ny = 50;
                nz = 50;
                [yy, zz] = meshgrid(1:nz, 1:ny);
                intensity_data = sin(2*pi*yy/20) .* cos(2*pi*zz/20) + 0.5;
            end
            
            [ny, nz] = size(intensity_data);
            if ny < 2 || nz < 2
                return;
            end
            
            intensity_min = min(intensity_data(:));
            intensity_max = max(intensity_data(:));
            if intensity_max > intensity_min
                intensity_norm = (intensity_data - intensity_min) / (intensity_max - intensity_min);
            else
                intensity_norm = ones(ny, nz) * 0.5;
            end
            
            fft2d = fftshift(fft2(intensity_norm));
            fft_mag = abs(fft2d);
            fft_log = log10(1 + fft_mag);
            
            if length(r_axis) > 1
                dr = r_axis(2) - r_axis(1);
            else
                dr = 1;
            end
            if length(z_axis) > 1
                dz = z_axis(2) - z_axis(1);
            else
                dz = 1;
            end
            
            r_freq = ((-floor(ny/2)):(ceil(ny/2)-1)) / (ny * dr);
            z_freq = ((-floor(nz/2)):(ceil(nz/2)-1)) / (nz * dz);

            r_pos = r_freq > 0;
            z_pos = z_freq > 0;
            if ~any(r_pos) || ~any(z_pos)
                return;
            end

            r_period = 1 ./ r_freq(r_pos);
            z_period = 1 ./ z_freq(z_pos);
            [r_period_sorted, r_idx] = sort(r_period, 'ascend');
            [z_period_sorted, z_idx] = sort(z_period, 'ascend');

            fft_log_pos = fft_log(r_pos, z_pos);
            fft_log_pos = fft_log_pos(r_idx, z_idx);

            % Use evenly spaced indices for display, with long periods on the left
            r_period_desc = flip(r_period_sorted);
            z_period_desc = flip(z_period_sorted);
            fft_log_pos = flipud(fliplr(fft_log_pos));
            r_idx_display = 1:numel(r_period_desc);
            z_idx_display = 1:numel(z_period_desc);

            cla(autocorr_handles.ax_2d);
            set(autocorr_handles.ax_2d, 'NextPlot', 'replacechildren');
            imagesc(autocorr_handles.ax_2d, z_idx_display, r_idx_display, fft_log_pos);
            set(autocorr_handles.ax_2d, 'YDir', 'normal');
            axis(autocorr_handles.ax_2d, 'auto');

            colormap(autocorr_handles.ax_2d, hot);
            if isempty(autocorr_handles.cb_fft_2d) || ~isvalid(autocorr_handles.cb_fft_2d)
                autocorr_handles.cb_fft_2d = colorbar(autocorr_handles.ax_2d);
            end
            set(autocorr_handles.cb_fft_2d, 'Visible', 'on');
            autocorr_handles.cb_fft_2d.Label.String = 'Log FFT Magnitude';

            xlabel(autocorr_handles.ax_2d, 'Period z (μm)', 'FontSize', 10);
            ylabel(autocorr_handles.ax_2d, 'Period r (μm)', 'FontSize', 10);
            title(autocorr_handles.ax_2d, '2D FFT (Period Domain)', 'FontSize', 11, 'FontWeight', 'bold');
            grid(autocorr_handles.ax_2d, 'off');
            set(autocorr_handles.ax_2d, 'UserData', struct('z_period', z_period_desc, 'r_period', r_period_desc));

            x_tick_idx = unique(round(linspace(1, numel(z_period_desc), min(6, numel(z_period_desc)))));
            y_tick_idx = unique(round(linspace(1, numel(r_period_desc), min(6, numel(r_period_desc)))));
            set(autocorr_handles.ax_2d, 'XTick', x_tick_idx);
            set(autocorr_handles.ax_2d, 'YTick', y_tick_idx);
            set(autocorr_handles.ax_2d, 'XTickLabel', arrayfun(@(v) sprintf('%.3g', v), z_period_desc(x_tick_idx), 'UniformOutput', false));
            set(autocorr_handles.ax_2d, 'YTickLabel', arrayfun(@(v) sprintf('%.3g', v), r_period_desc(y_tick_idx), 'UniformOutput', false));

            % Use selected radial position for 1D FFT along z (axial)
            v_idx_axial = round(get(slider_axial_h, 'Value'));
            v_idx_axial = max(1, min(v_idx_axial, num_r));
            z_slice = intensity_norm(v_idx_axial, :);
            fft_z_1d = fftshift(fft(z_slice));
            fft_z_mag = abs(fft_z_1d);
            fft_z_log = log10(1 + fft_z_mag);
            z_freq_1d = ((-floor(nz/2)):(ceil(nz/2)-1)) / (nz * dz);
            z_pos_1d = z_freq_1d > 0;
            z_period_1d = 1 ./ z_freq_1d(z_pos_1d);
            [z_period_1d_sorted, z_idx_1d] = sort(z_period_1d, 'ascend');
            fft_z_line = fft_z_log(z_pos_1d);
            fft_z_line = fft_z_line(z_idx_1d);

            z_period_1d_desc = flip(z_period_1d_sorted);
            fft_z_line = fliplr(fft_z_line);
            z_idx_display_1d = 1:numel(z_period_1d_desc);

            cla(autocorr_handles.ax_z);
            set(autocorr_handles.ax_z, 'NextPlot', 'replacechildren');
            plot(autocorr_handles.ax_z, z_idx_display_1d, fft_z_line, 'Color', [1 0.5 0], 'LineWidth', 2.5);
            hold(autocorr_handles.ax_z, 'off');
            set(autocorr_handles.ax_z, 'XGrid', 'on', 'YGrid', 'on');
            xlabel(autocorr_handles.ax_z, 'Period z (μm)', 'FontSize', 9);
            ylabel(autocorr_handles.ax_z, 'Log |FFT|', 'FontSize', 9);
            title(autocorr_handles.ax_z, 'Z-Axis FFT', 'FontSize', 10, 'FontWeight', 'bold');
            xlim(autocorr_handles.ax_z, [1 numel(z_period_1d_desc)]);
            ylim(autocorr_handles.ax_z, [min(fft_z_line) - 0.05, max(fft_z_line) + 0.1]);
            z_tick_idx_1d = unique(round(linspace(1, numel(z_period_1d_desc), min(6, numel(z_period_1d_desc)))));
            set(autocorr_handles.ax_z, 'XTick', z_tick_idx_1d);
            set(autocorr_handles.ax_z, 'XTickLabel', arrayfun(@(v) sprintf('%.3g', v), z_period_1d_desc(z_tick_idx_1d), 'UniformOutput', false));
            set(autocorr_handles.ax_z, 'UserData', struct('z_period', z_period_1d_desc));

            % Use selected radial slice for 1D FFT along r
            u_idx_radial = round(get(slider_radial_h, 'Value'));
            u_idx_radial = max(1, min(u_idx_radial, num_u));
            r_slice = intensity_norm(:, u_idx_radial);
            fft_r_1d = fftshift(fft(r_slice));
            fft_r_mag = abs(fft_r_1d);
            fft_r_log = log10(1 + fft_r_mag);
            r_freq_1d = ((-floor(ny/2)):(ceil(ny/2)-1)) / (ny * dr);
            r_pos_1d = r_freq_1d > 0;
            r_period_1d = 1 ./ r_freq_1d(r_pos_1d);
            [r_period_1d_sorted, r_idx_1d] = sort(r_period_1d, 'ascend');
            fft_r_line = fft_r_log(r_pos_1d);
            fft_r_line = fft_r_line(r_idx_1d);

            r_period_1d_desc = flip(r_period_1d_sorted);
            fft_r_line = flipud(fft_r_line);
            r_idx_display_1d = 1:numel(r_period_1d_desc);

            cla(autocorr_handles.ax_r);
            set(autocorr_handles.ax_r, 'NextPlot', 'replacechildren');
            plot(autocorr_handles.ax_r, r_idx_display_1d, fft_r_line, 'Color', [0 0.7 1], 'LineWidth', 2.5);
            hold(autocorr_handles.ax_r, 'off');
            set(autocorr_handles.ax_r, 'XGrid', 'on', 'YGrid', 'on');
            xlabel(autocorr_handles.ax_r, 'Period r (μm)', 'FontSize', 9);
            ylabel(autocorr_handles.ax_r, 'Log |FFT|', 'FontSize', 9);
            title(autocorr_handles.ax_r, 'R-Axis FFT', 'FontSize', 10, 'FontWeight', 'bold');
            xlim(autocorr_handles.ax_r, [1 numel(r_period_1d_desc)]);
            ylim(autocorr_handles.ax_r, [min(fft_r_line) - 0.05, max(fft_r_line) + 0.1]);
            r_tick_idx = unique(round(linspace(1, numel(r_period_1d_desc), min(6, numel(r_period_1d_desc)))));
            set(autocorr_handles.ax_r, 'XTick', r_tick_idx);
            set(autocorr_handles.ax_r, 'XTickLabel', arrayfun(@(v) sprintf('%.3g', v), r_period_1d_desc(r_tick_idx), 'UniformOutput', false));
            set(autocorr_handles.ax_r, 'UserData', struct('r_period', r_period_1d_desc));
            
            drawnow();
            
        catch
        end
    end

    function txt = fft_datatip_txt(~, event_obj)
        pos = get(event_obj, 'Position');
        target = get(event_obj, 'Target');
        ax = ancestor(target, 'axes');
        txt = {
            sprintf('X: %.3f', pos(1))
            sprintf('Y: %.3f', pos(2))
        };
        if isempty(ax) || ~isvalid(ax)
            return;
        end
        if ax == autocorr_handles.ax_2d
            ud = get(ax, 'UserData');
            if isstruct(ud) && isfield(ud, 'z_period') && isfield(ud, 'r_period')
                xi = round(pos(1));
                yi = round(pos(2));
                xi = max(1, min(xi, numel(ud.z_period)));
                yi = max(1, min(yi, numel(ud.r_period)));
                txt = {
                    sprintf('Period z (um): %.4g', ud.z_period(xi))
                    sprintf('Period r (um): %.4g', ud.r_period(yi))
                    sprintf('Log |FFT|: %.4g', pos(3))
                };
            end
        elseif ax == autocorr_handles.ax_z
            ud = get(ax, 'UserData');
            if isstruct(ud) && isfield(ud, 'z_period')
                xi = round(pos(1));
                xi = max(1, min(xi, numel(ud.z_period)));
                txt = {
                    sprintf('Period z (um): %.4g', ud.z_period(xi))
                    sprintf('Log |FFT|: %.4g', pos(2))
                };
            end
        elseif ax == autocorr_handles.ax_r
            ud = get(ax, 'UserData');
            if isstruct(ud) && isfield(ud, 'r_period')
                xi = round(pos(1));
                xi = max(1, min(xi, numel(ud.r_period)));
                txt = {
                    sprintf('Period r (um): %.4g', ud.r_period(xi))
                    sprintf('Log |FFT|: %.4g', pos(2))
                };
            end
        end
    end


    function B = pad_array(A, padsize, padval, direction)
        % Custom padding function without Image Processing Toolbox dependency
        % Syntax: B = pad_array(A, padsize, padval, direction)
        % direction: 'post' (default), 'pre', or 'both'
        
        if nargin < 4, direction = 'post'; end
        if nargin < 3, padval = 0; end
        
        [m, n] = size(A);
        if length(padsize) == 1
            padsize = [padsize, padsize];
        end
        
        if strcmp(direction, 'post')
            B = padval * ones(m + padsize(1), n + padsize(2));
            B(1:m, 1:n) = A;
        elseif strcmp(direction, 'pre')
            B = padval * ones(m + padsize(1), n + padsize(2));
            B(padsize(1)+1:end, padsize(2)+1:end) = A;
        else
            B = padval * ones(m + 2*padsize(1), n + 2*padsize(2));
            B(padsize(1)+1:padsize(1)+m, padsize(2)+1:padsize(2)+n) = A;
        end
    end

fprintf('Viewer ready!\n');
fprintf('Full-load mode: all slices in memory for instant navigation (~1-2ms per frame).\n');
fprintf('Tip: Keep secondary plot checkboxes OFF for maximum speed.\n\n');

update_plot();
end
