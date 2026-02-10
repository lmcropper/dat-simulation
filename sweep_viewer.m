function sweep_viewer()
% Sweep Viewer: Load and browse sweep results interactively
% This script loads saved sweep results and provides an interactive viewer

close all;

% ========== LOAD SWEEP RESULTS ==========
[filename, pathname] = uigetfile('*.mat', 'Select sweep results file');
if isequal(filename, 0)
    disp('No file selected.');
    return;
end

file_path = fullfile(pathname, filename);
data = load(file_path);

required_fields = {'intensity_sweep', 'sweep_param_values', 'sweep_param', ...
    'u_values', 'v_values', 'z_scale', 'r_scale'};
missing_fields = required_fields(~isfield(data, required_fields));
if ~isempty(missing_fields)
    error('Missing required fields in file: %s', strjoin(missing_fields, ', '));
end

intensity_sweep = data.intensity_sweep;
sweep_param_values = data.sweep_param_values;
sweep_param = data.sweep_param;
u_values = data.u_values;
v_values = data.v_values;
z_scale = data.z_scale;
r_scale = data.r_scale;

% Ensure axis vectors match the intensity grid dimensions
num_r = size(intensity_sweep, 1);
num_u = size(intensity_sweep, 2);
if length(u_values) ~= num_u || length(v_values) ~= num_r
    if isfield(data, 'z_min') && isfield(data, 'z_max') && isfield(data, 'r_min') && isfield(data, 'r_max')
        u_values = linspace(data.z_min / z_scale, data.z_max / z_scale, num_u);
        v_values = linspace(data.r_min / r_scale, data.r_max / r_scale, num_r);
    else
        u_values = linspace(0, 1, num_u);
        v_values = linspace(0, 1, num_r);
    end
end

fprintf('Loaded: %s\n', filename);
fprintf('Sweep parameter: %s\n', sweep_param);
fprintf('Matrix size: %d x %d x %d\n', size(intensity_sweep, 1), size(intensity_sweep, 2), size(intensity_sweep, 3));
fprintf('Parameter range: %.3f to %.3f\n\n', sweep_param_values(1), sweep_param_values(end));

% ========== CREATE INTERACTIVE VIEWER ==========
% Create figure with controls
fig = figure('Name', sprintf('Sweep Viewer: %s', sweep_param), ...
    'NumberTitle', 'off', ...
    'Position', [100 100 1200 900]);

% Create main axes
ax_main = axes('Position', [0.1 0.35 0.8 0.6]);

% Create slider
slider_h = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', length(sweep_param_values), 'Value', 1, ...
    'Position', [100 280 500 20], ...
    'Callback', @update_plot);

% Create text label for parameter value
text_info = uicontrol('Style', 'text', ...
    'Position', [100 300 500 40], ...
    'FontSize', 12, ...
    'BackgroundColor', [0.94 0.94 0.94], ...
    'HorizontalAlignment', 'left');

% Create previous/next buttons
btn_prev = uicontrol('Style', 'pushbutton', 'String', '< Previous', ...
    'Position', [100 240 80 30], ...
    'Callback', @(src, evt) move_slider(-1));

btn_next = uicontrol('Style', 'pushbutton', 'String', 'Next >', ...
    'Position', [190 240 80 30], ...
    'Callback', @(src, evt) move_slider(1));

% Create step size spinner
uicontrol('Style', 'text', 'String', 'Step size:', ...
    'Position', [290 245 60 20], 'HorizontalAlignment', 'right');

spinner_step = uicontrol('Style', 'edit', 'String', '1', ...
    'Position', [360 245 50 20], ...
    'BackgroundColor', [1 1 1]);

% Create play/pause button
btn_play = uicontrol('Style', 'pushbutton', 'String', 'Play', ...
    'Position', [420 240 60 30], ...
    'Callback', @play_animation);

is_playing = false;

% Create radial line plot checkbox
cb_axial = uicontrol('Style', 'checkbox', 'String', 'Show axial intensity plot', ...
    'Position', [100 120 200 20], 'Value', 1, ...
    'Callback', @update_plot);

cb_radial = uicontrol('Style', 'checkbox', 'String', 'Show radial intensity plot', ...
    'Position', [100 90 200 20], 'Value', 0, ...
    'Callback', @update_plot);

% Create info text
info_text = uicontrol('Style', 'text', ...
    'Position', [100 10 400 70], ...
    'FontSize', 10, ...
    'BackgroundColor', [0.94 0.94 0.94], ...
    'HorizontalAlignment', 'left');

% ========== CALLBACK FUNCTIONS ==========
    function update_plot(src, evt)
        idx = round(get(slider_h, 'Value'));
        
        % Clamp index
        idx = max(1, min(idx, length(sweep_param_values)));
        set(slider_h, 'Value', idx);
        
        % Get current data
        intensity_current = intensity_sweep(:, :, idx);
        param_value = sweep_param_values(idx);
        
        % Plot 2D intensity
        cla(ax_main);
        hold(ax_main, 'on');
        
        % Normalize for display
        gamma = 0.5;
        intensity_gamma = intensity_current.^gamma;
        
        imagesc(ax_main, u_values*z_scale*1e6, v_values*r_scale*1e6, intensity_gamma);
        colormap(ax_main, hot);
        colorbar(ax_main);
        xlabel(ax_main, 'Axial Distance z (μm)');
        ylabel(ax_main, 'Radial Distance r (μm)');
        title(ax_main, sprintf('Intensity Pattern - %s = %.3f', sweep_param, param_value));
        
        % Update info text
        set(text_info, 'String', ...
            sprintf('%s = %.3f (%d / %d)', sweep_param, param_value, idx, length(sweep_param_values)));
        
        % Add secondary plots if requested
        show_axial = get(cb_axial, 'Value');
        show_radial = get(cb_radial, 'Value');
        
        if show_axial || show_radial
            if show_axial && show_radial
                num_subplots = 2;
                pos_axial = [0.1 0.05 0.35 0.25];
                pos_radial = [0.55 0.05 0.35 0.25];
            elseif show_axial
                num_subplots = 1;
                pos_axial = [0.1 0.05 0.8 0.25];
            else
                num_subplots = 1;
                pos_radial = [0.1 0.05 0.8 0.25];
            end
            
            if show_axial
                % Plot intensity along optical axis
                ax_axial = axes('Position', pos_axial);
                [~, v_zero_idx] = min(abs(v_values));
                intensity_axis = intensity_current(v_zero_idx, :);
                plot(ax_axial, u_values*z_scale*1e6, intensity_axis, 'LineWidth', 2, 'Color', [1 0.5 0]);
                xlabel(ax_axial, 'Axial Distance z (μm)');
                ylabel(ax_axial, 'Normalized Intensity');
                title(ax_axial, 'Intensity Along Optical Axis (r=0)');
                grid(ax_axial, 'on');
            end
            
            if show_radial
                % Plot intensity along radial axis
                ax_radial = axes('Position', pos_radial);
                [~, u_zero_idx] = min(abs(u_values));
                intensity_rad = intensity_current(:, u_zero_idx);
                r_vals = v_values * r_scale * 1e6;
                plot(ax_radial, r_vals, intensity_rad, 'LineWidth', 2, 'Color', [0 0.7 1]);
                xlabel(ax_radial, 'Radial Distance r (μm)');
                ylabel(ax_radial, 'Normalized Intensity');
                title(ax_radial, 'Intensity Along Radial Axis (z=z_min)');
                grid(ax_radial, 'on');
            end
        end
        
        % Update info text
        set(info_text, 'String', ...
            sprintf(['Sweep Index: %d / %d\n', ...
                     '%s = %.4f\n', ...
                     'Matrix size: %d × %d\n', ...
                     'Max intensity: %.4f\n', ...
                     'File: %s'], ...
                idx, length(sweep_param_values), ...
                sweep_param, param_value, ...
                size(intensity_current, 1), size(intensity_current, 2), ...
                max(intensity_current(:)), ...
                filename));
    end

    function move_slider(direction)
        try
            step_size = str2double(get(spinner_step, 'String'));
        catch
            step_size = 1;
            set(spinner_step, 'String', '1');
        end
        
        current_val = get(slider_h, 'Value');
        new_val = current_val + direction * step_size;
        new_val = max(1, min(new_val, length(sweep_param_values)));
        set(slider_h, 'Value', new_val);
        update_plot();
    end

    function play_animation(src, evt)
        is_playing = ~is_playing;
        
        if is_playing
            set(btn_play, 'String', 'Pause');
            
            while is_playing && get(slider_h, 'Value') < length(sweep_param_values)
                move_slider(1);
                pause(0.2);  % Delay between frames
            end
            
            % Stop at end
            is_playing = false;
            set(btn_play, 'String', 'Play');
        else
            set(btn_play, 'String', 'Play');
        end
    end

% ========== INITIAL PLOT ==========
update_plot();
end
