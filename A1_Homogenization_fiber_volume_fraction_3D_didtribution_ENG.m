clear all;
clc;
% Read Excel file and convert to 3D matrix
filename = '3D_greyscale_voxel.xls'; % Replace with your Excel filename
[~, sheets] = xlsfinfo(filename);
num_sheets = numel(sheets);

% Read first sheet to get data dimensions
first_sheet = xlsread(filename, sheets{1});
[rows, cols] = size(first_sheet);
data_3d = zeros(rows, cols, num_sheets); % Initialize 3D matrix

% Read all sheets in order
for k = 1:num_sheets
    sheet_name = sheets{k};
    % Try to convert sheet name to number
    if ~isnan(str2double(sheet_name))
        current_data = xlsread(filename, str2double(sheet_name));
    else
        current_data = xlsread(filename, sheet_name);
    end
    
    % Ensure data dimensions are consistent
    if all(size(current_data) == [rows, cols])
        data_3d(:, :, k) = current_data;
    else
        error(['Sheet ' sheet_name ' has inconsistent dimensions with other sheets']);
    end
end

% Set block size
block_size = 15; % Block size for each dimension

% Calculate number of blocks
num_blocks_x = ceil(cols / block_size);
num_blocks_y = ceil(rows / block_size);
num_blocks_z = ceil(num_sheets / block_size);

% Initialize block statistics matrices
block_ratios = zeros(num_blocks_y, num_blocks_x, num_blocks_z);
block_centers_x = zeros(num_blocks_y, num_blocks_x, num_blocks_z);
block_centers_y = zeros(num_blocks_y, num_blocks_x, num_blocks_z);
block_centers_z = zeros(num_blocks_y, num_blocks_x, num_blocks_z);

% Calculate non-zero element ratio for each block
% Calculation the volume fraction of every block
for z_block = 1:num_blocks_z
    z_start = (z_block - 1) * block_size + 1;
    z_end = min(z_block * block_size, num_sheets);
    
    for y_block = 1:num_blocks_y
        y_start = (y_block - 1) * block_size + 1;
        y_end = min(y_block * block_size, rows);
        
        for x_block = 1:num_blocks_x
            x_start = (x_block - 1) * block_size + 1;
            x_end = min(x_block * block_size, cols);
            
            % Extract current block data
            block_data = data_3d(y_start:y_end, x_start:x_end, z_start:z_end);
            
            % Calculate non-zero element ratio
            non_zero_count = nnz(block_data > 0);
            total_elements = numel(block_data);
            ratio = non_zero_count / total_elements;
            
            % Store ratio value
            block_ratios(y_block, x_block, z_block) = ratio;
            
            % Calculate block center position (apply scaling factor)
            block_centers_x(y_block, x_block, z_block) = mean([x_start, x_end]) * 3.4;
            block_centers_y(y_block, x_block, z_block) = mean([y_start, y_end]) * 3.4;
            block_centers_z(y_block, x_block, z_block) = mean([z_start, z_end]) * 3.4;
        end
    end
end

% Create 3D visualization
figure;
hold on;
grid on;
view(3); % 3D view
xlabel('X/um');
ylabel('Y/um');
zlabel('Z/um');
title('Visualization of fiber volume fraction') %[Homogenization degree: 5 voxels]');

% Get minimum and maximum ratio values
min_ratio = min(block_ratios(:));
max_ratio = max(block_ratios(:));

% Create color map (from blue to red)
color_map = jet(256); % Using jet colormap, can also use hot, parula, etc.

% Plot each block
for z = 1:num_blocks_z
    for y = 1:num_blocks_y
        for x = 1:num_blocks_x
            ratio_val = block_ratios(y, x, z);
            
            % Only plot blocks with non-zero ratio
            if ratio_val > 0
                % Calculate color index
                color_idx = round(1 + (ratio_val - min_ratio) / (max_ratio - min_ratio) * 255);
                color_idx = max(1, min(256, color_idx)); % Ensure within 1-256 range
                block_color = color_map(color_idx, :);
                
                % Draw cube representing the block
                plot_cube(block_centers_x(y, x, z), ...
                          block_centers_y(y, x, z), ...
                          block_centers_z(y, x, z), ...
                          block_size * 3.4, block_size * 3.4, block_size * 3.4, ...
                          block_color);
            end
        end
    end
end

% Add color bar
colormap(color_map);
caxis([min_ratio max_ratio]);
colorbar;
title(colorbar, 'Fiber volume fraction');

hold off;

% Cube plotting function
function plot_cube(x, y, z, dx, dy, dz, color)
    % Define 8 vertices of the cube
    vertices = [x-dx/2, y-dy/2, z-dz/2;
                x+dx/2, y-dy/2, z-dz/2;
                x+dx/2, y+dy/2, z-dz/2;
                x-dx/2, y+dy/2, z-dz/2;
                x-dx/2, y-dy/2, z+dz/2;
                x+dx/2, y-dy/2, z+dz/2;
                x+dx/2, y+dy/2, z+dz/2;
                x-dx/2, y+dy/2, z+dz/2];
    
    % Define 6 faces of the cube
    faces = [1,2,3,4;   % Bottom face
             5,6,7,8;   % Top face
             1,2,6,5;   % Front face
             3,4,8,7;   % Back face
             1,4,8,5;   % Left face
             2,3,7,6];  % Right face
    
    % Draw semi-transparent cube
    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', color, 'EdgeColor', 'k', ...
          'FaceAlpha', 0.7, 'EdgeAlpha', 0.5);
end