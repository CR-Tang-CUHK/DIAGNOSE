clear all;
clc;
%% Load unit vector data and original CT data
% Assuming previously saved files are in current directory
vx_file = 'vector_components_x.xlsx';
vy_file = 'vector_components_y.xlsx';
vz_file = 'vector_components_z.xlsx';
filename = '3D_greyscale_voxel_table_of_fiber_masks_top_right_3.4um_cut a section.xls'; % Original CT data

% Get sheet information
[~, sheets] = xlsfinfo(vx_file);
num_sheets = length(sheets);

% Ensure sheets is a cell array
if isstring(sheets)
    sheets = cellstr(sheets);
end

% Check sheet names, convert to string if numeric
sheet_names = cell(1, num_sheets);
for i = 1:num_sheets
    if isnumeric(sheets{i})
        sheet_names{i} = num2str(sheets{i});
    else
        sheet_names{i} = sheets{i};
    end
end

% Read first sheet to get matrix dimensions
first_sheet = xlsread(vx_file, sheet_names{1});
[m, n] = size(first_sheet);

% Initialize vector component matrices
vx_3D = zeros(m, n, num_sheets);
vy_3D = zeros(m, n, num_sheets);
vz_3D = zeros(m, n, num_sheets);

% Fill vector component matrices
for k = 1:num_sheets
    vx_3D(:,:,k) = xlsread(vx_file, sheet_names{k});
    vy_3D(:,:,k) = xlsread(vy_file, sheet_names{k});
    vz_3D(:,:,k) = xlsread(vz_file, sheet_names{k});
end

% Load original CT data
ct_data_3D = zeros(m, n, num_sheets);
for k = 1:num_sheets
    % Try reading using both numeric and string names
    try
        ct_data_3D(:,:,k) = xlsread(filename, sheet_names{k});
    catch
        % If name fails, try using index
        try
            ct_data_3D(:,:,k) = xlsread(filename, k);
        catch
            % If index fails, try using default name
            try
                ct_data_3D(:,:,k) = xlsread(filename, 'Sheet1');
            catch
                % If all methods fail, display error message
                error('Cannot read sheet %d, please check Excel file format', k);
            end
        end
    end
end

%% Create fiber region mask
% According to provided grayscale information, fiber parts have grayscale 100 or 200, background is 0
fiber_mask = ct_data_3D > 50; % Threshold set to 50 to ensure capturing fibers with grayscale 100 and 200

% Use morphological operations to improve mask
se = strel('sphere', 1);
fiber_mask = imopen(fiber_mask, se);
fiber_mask = imclose(fiber_mask, se);

%% Region averaging (3x3x3 cube) and apply fiber mask
% Set region size
block_size = 3;
half_block = floor(block_size/2);

% Calculate new matrix dimensions
new_m = floor(m / block_size);
new_n = floor(n / block_size);
new_k = floor(num_sheets / block_size);

% Initialize averaged vector matrices and mask
vx_avg = zeros(new_m, new_n, new_k);
vy_avg = zeros(new_m, new_n, new_k);
vz_avg = zeros(new_m, new_n, new_k);
mask_avg = false(new_m, new_n, new_k); % Averaged mask

% Calculate average vector and mask for each region
for i = 1:new_m
    for j = 1:new_n
        for k = 1:new_k
            % Calculate current region indices
            i_start = (i-1)*block_size + 1;
            i_end = i*block_size;
            j_start = (j-1)*block_size + 1;
            j_end = j*block_size;
            k_start = (k-1)*block_size + 1;
            k_end = k*block_size;
            
            % Extract current region vector components and mask
            vx_block = vx_3D(i_start:i_end, j_start:j_end, k_start:k_end);
            vy_block = vy_3D(i_start:i_end, j_start:j_end, k_start:k_end);
            vz_block = vz_3D(i_start:i_end, j_start:j_end, k_start:k_end);
            mask_block = fiber_mask(i_start:i_end, j_start:j_end, k_start:k_end);
            
            % Calculate average vector
            avg_vx = mean(vx_block(:));
            avg_vy = mean(vy_block(:));
            avg_vz = mean(vz_block(:));
            
            % Calculate fiber voxel ratio in region
            fiber_ratio = sum(mask_block(:)) / numel(mask_block);
            
            % Only keep region if fiber ratio exceeds threshold
            if fiber_ratio > 0.3 % Threshold set to 30%
                % Normalize to unit vector
                norm_avg = norm([avg_vx, avg_vy, avg_vz]);
                if norm_avg > 0.1 % Only keep meaningful vectors
                    vx_avg(i,j,k) = avg_vx / norm_avg;
                    vy_avg(i,j,k) = avg_vy / norm_avg;
                    vz_avg(i,j,k) = avg_vz / norm_avg;
                    mask_avg(i,j,k) = true;
                end
            end
        end
    end
end

%% Create coordinate grid (only containing fiber regions)
[X, Y, Z] = meshgrid(1:new_n, 1:new_m, 1:new_k);

% Only keep fiber regions
valid_indices = find(mask_avg);
X_valid = X(valid_indices);
Y_valid = Y(valid_indices);
Z_valid = Z(valid_indices);
vx_valid = vx_avg(valid_indices);
vy_valid = vy_avg(valid_indices);
vz_valid = vz_avg(valid_indices);

%% Calculate color-coded angles (only for valid regions)
% XY plane angle with X-axis (0-180 degree range)
angle_xy = atan2(vy_valid, vx_valid); % radians
angle_xy_deg = mod(rad2deg(angle_xy), 180); % Convert to 0-180 degree range

% ZX plane angle with X-axis (0-180 degree range)
angle_zx = atan2(vz_valid, vx_valid); % radians
angle_zx_deg = mod(rad2deg(angle_zx), 180); % Convert to 0-180 degree range

% YZ plane angle with Y-axis (0-180 degree range)
angle_yz = atan2(vz_valid, vy_valid); % radians
angle_yz_deg = mod(rad2deg(angle_yz), 180); % Convert to 0-180 degree range

%% Create 3D plot - colored by XY plane angle with X-axis
figure;
set(gcf, 'Position', [100, 100, 1000, 800]);
hold on;

% Use quiver3 to draw arrows, colored by XY plane angle
for i = 1:length(X_valid)
    % Get current point color (based on XY plane angle)
    color_val = angle_xy_deg(i) / 180; % Normalized to 0-1 range
    color = hsv2rgb([color_val, 1, 1]); % Use HSV color space, hue determined by angle
    
    % Draw arrow
    quiver3(X_valid(i), Y_valid(i), Z_valid(i), ...
            vx_valid(i), vy_valid(i), vz_valid(i), ...
            0.5, 'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
end

% Set graph properties
axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Unit Vector 3D Plot - Colored by XY Plane Angle with X-axis (Fiber Regions Only)');
view(3);
rotate3d on;

% Add color bar
colormap(hsv);
caxis([0 180]); % Set color range corresponding to 0-180 degrees
colorbar('Ticks', 0:30:180, 'TickLabels', {'0°', '30°', '60°', '90°', '120°', '150°', '180°'});

hold off;

%% Create 3D plot - colored by ZX plane angle with X-axis
figure;
set(gcf, 'Position', [100, 100, 1000, 800]);
hold on;

% Use quiver3 to draw arrows, colored by ZX plane angle
for i = 1:length(X_valid)
    % Get current point color (based on ZX plane angle)
    color_val = angle_zx_deg(i) / 180; % Normalized to 0-1 range
    color = hsv2rgb([color_val, 1, 1]); % Use HSV color space, hue determined by angle
    
    % Draw arrow
    quiver3(X_valid(i), Y_valid(i), Z_valid(i), ...
            vx_valid(i), vy_valid(i), vz_valid(i), ...
            0.5, 'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
end

% Set graph properties
axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Unit Vector 3D Plot - Colored by ZX Plane Angle with X-axis (Fiber Regions Only)');
view(3);
rotate3d on;

% Add color bar
colormap(hsv);
caxis([0 180]); % Set color range corresponding to 0-180 degrees
colorbar('Ticks', 0:30:180, 'TickLabels', {'0°', '30°', '60°', '90°', '120°', '150°', '180°'});

hold off;

%% Create 3D plot - colored by YZ plane angle with Y-axis
figure;
set(gcf, 'Position', [100, 100, 1000, 800]);
hold on;

% Use quiver3 to draw arrows, colored by YZ plane angle
for i = 1:length(X_valid)
    % Get current point color (based on YZ plane angle)
    color_val = angle_yz_deg(i) / 180; % Normalized to 0-1 range
    color = hsv2rgb([color_val, 1, 1]); % Use HSV color space, hue determined by angle
    
    % Draw arrow
    quiver3(X_valid(i), Y_valid(i), Z_valid(i), ...
            vx_valid(i), vy_valid(i), vz_valid(i), ...
            0.5, 'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
end

% Set graph properties
axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Unit Vector 3D Plot - Colored by YZ Plane Angle with Y-axis (Fiber Regions Only)');
view(3);
rotate3d on;

% Add color bar
colormap(hsv);
caxis([0 180]); % Set color range corresponding to 0-180 degrees
colorbar('Ticks', 0:30:180, 'TickLabels', {'0°', '30°', '60°', '90°', '120°', '150°', '180°'});

hold off;

%% Create composite plot - three subplots showing different angle color coding
figure;
set(gcf, 'Position', [100, 100, 1400, 500]);

% Subplot 1: XY plane angle with X-axis
subplot(1,3,1);
scatter3(X_valid, Y_valid, Z_valid, 20, angle_xy_deg, 'filled');
hold on;
quiver3(X_valid, Y_valid, Z_valid, vx_valid, vy_valid, vz_valid, 0.5, 'k', 'LineWidth', 0.5);
hold off;
axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('XY Plane Angle with X-axis (Fiber Regions Only)');
colormap(hsv);
caxis([0 180]);
colorbar('Ticks', 0:30:180, 'TickLabels', {'0°', '30°', '60°', '90°', '120°', '150°', '180°'});
view(3);

% Subplot 2: ZX plane angle with X-axis
subplot(1,3,2);
scatter3(X_valid, Y_valid, Z_valid, 20, angle_zx_deg, 'filled');
hold on;
quiver3(X_valid, Y_valid, Z_valid, vx_valid, vy_valid, vz_valid, 0.5, 'k', 'LineWidth', 0.5);
hold off;
axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('ZX Plane Angle with X-axis (Fiber Regions Only)');
colormap(hsv);
caxis([0 180]);
colorbar('Ticks', 0:30:180, 'TickLabels', {'0°', '30°', '60°', '90°', '120°', '150°', '180°'});
view(3);

% Subplot 3: YZ plane angle with Y-axis
subplot(1,3,3);
scatter3(X_valid, Y_valid, Z_valid, 20, angle_yz_deg, 'filled');
hold on;
quiver3(X_valid, Y_valid, Z_valid, vx_valid, vy_valid, vz_valid, 0.5, 'k', 'LineWidth', 0.5);
hold off;
axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('YZ Plane Angle with Y-axis (Fiber Regions Only)');
colormap(hsv);
caxis([0 180]);
colorbar('Ticks', 0:30:180, 'TickLabels', {'0°', '30°', '60°', '90°', '120°', '150°', '180°'});
view(3);

%% Display statistical information
fprintf('Statistical Information:\n');
fprintf('Original data dimensions: %d x %d x %d\n', m, n, num_sheets);
fprintf('After region averaging: %d x %d x %d\n', new_m, new_n, new_k);
fprintf('Number of fiber regions: %d (out of %d total regions)\n', length(X_valid), new_m*new_n*new_k);
fprintf('Fiber region ratio: %.2f%%\n', 100*length(X_valid)/(new_m*new_n*new_k));