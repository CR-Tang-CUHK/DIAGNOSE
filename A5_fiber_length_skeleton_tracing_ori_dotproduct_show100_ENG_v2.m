%% Close/clear all
close all
clear variables
clc

%% Variables
filename = "3D_greyscale_voxel"; % Filename for the voxel file

voxel_size = 3.4; % Voxel size in micrometers
diameter = 15e-6; % Fiber diameter in meters (used for analysis based FVF calculation)

% Skeleonization parameters
algorithm = 'bwulterode&bwskel'; % Use bwskel, bwulterode or bwulterode&bwskel algorithm
skel_size = 4; % Size for bwulterode

% Orientation parameters
neighborhood_size = 3; % Use larger neighborhood for more stable direction estimation

% Tracking parameters
min_fiber_length = 30; % Minimum fiber length in µm
max_step_size = 3; % Maximum step size
max_angle_threshold = 45; % Maximum allowed angle change (degrees)

%% Some initialization
min_fiber_length = ceil(min_fiber_length/voxel_size); % convert min_fiber_length to a number of voxels

%% Convert data if not done yet
if exist(filename + ".mat","file")
    load(filename + ".mat")
else
    fprintf("Converting 3D voxel file\n");
    tic;
    % Read first sheet to get data dimensions
    sheets = sheetnames(filename + ".xls");
    data = readmatrix(filename + ".xls","Sheet",sheets(1));
    N_x = size(data,1);
    N_y = size(data,2);
    N_z = length(sheets);
    data_3d = zeros(N_x,N_y,N_z, "uint16");

    % Read all sheets in order
    for i=1:length(sheets)
        data_3d(:,:,i) = readmatrix(filename + ".xls","Sheet",sheets(i));
    end
    volume_data = zeros(size(data_3d));
    volume_data(data_3d>0) = 1;
    save(filename + ".mat","data_3d","volume_data","N_x","N_y","N_z");
    fprintf("done with converting 3D voxel file\n");
    toc;
end

%% 3. Skeletonization processing - get fiber centerlines
fprintf('Performing skeletonization processing...\n');

% Use MATLAB Image Processing Toolbox for 3D skeletonization
if license('test', 'image_toolbox')
    switch algorithm
        case 'bwskel'
            % Use bwskel function for skeletonization
            skeleton_data = bwskel(volume_data > 0);
            fprintf('Skeletonization completed using bwskel\n');
        case 'bwulterode'
            % Use bwulterode function for skeletonization
            skeleton_data = bwulterode(volume_data>0,'euclidean',skel_size);
            fprintf('Skeletonization completed using bwulterode\n');
        case 'bwulterode&bwskel'
            % Use bwulterode function for skeletonization
            skeleton_data = bwskel(bwulterode(volume_data>0,'euclidean',skel_size));
            fprintf('Skeletonization completed using bwulterode & bwskel\n');
    end
else
    % If Image Processing Toolbox is not available, use simplified skeletonization method
    fprintf('Image Processing Toolbox not detected, using simplified skeletonization method\n');
    skeleton_data = simple_skeleton_3d(volume_data);
end

skeleton_data = bwmorph3(skeleton_data,'clean'); % Clean data (removes lonely voxels which saves some time later on)

% Check skeletonization results
skeleton_voxels = sum(skeleton_data(:));
fprintf('Skeletonization completed, number of skeleton voxels: %d\n', skeleton_voxels);

if skeleton_voxels == 0
    error('Skeletonization failed, no skeleton voxels detected.');
end

% Plot the skeleton voxels
figure(1); t = tiledlayout(1,2,"TileSpacing","compact","Padding","compact"); 
nexttile(1); hold on; grid on;
[x,y,z] = ndgrid(voxel_size.*(1:N_x),voxel_size.*(1:N_y),voxel_size.*(1:N_z));
scatter3(x(skeleton_data),y(skeleton_data),z(skeleton_data),[],double(skeleton_data(skeleton_data)),'.k');
axis equal tight; view(3);
title("Skeletonized voxel image");
drawnow;
    
%% 4. Calculate structure tensor and direction vector field (on skeleton)
fprintf('Calculating structure tensor and direction vector field on skeleton...\n');

% Calculate direction vectors for skeleton voxels
vector_field = zeros(N_x, N_y, N_z, 3);

% Get coordinates of all skeleton voxels
[skel_x, skel_y, skel_z] = ind2sub(size(skeleton_data), find(skeleton_data > 0));
num_skel_voxels = length(skel_x);

for i = 1:num_skel_voxels
    x = skel_x(i);
    y = skel_y(i);
    z = skel_z(i);

    % Define neighborhood range (on original volume data)
    x_min = max(1, x-neighborhood_size);
    x_max = min(N_x, x+neighborhood_size);
    y_min = max(1, y-neighborhood_size);
    y_max = min(N_y, y+neighborhood_size);
    z_min = max(1, z-neighborhood_size);
    z_max = min(N_z, z+neighborhood_size);

    % Extract neighborhood (using original fiber data, not skeleton data)
    neighborhood = volume_data(x_min:x_max, y_min:y_max, z_min:z_max);

    if sum(neighborhood(:)) > 5 % Enough fiber points
        % Calculate gradient
        [gx, gy, gz] = gradient(double(neighborhood));

        % Calculate structure tensor
        J11 = sum(gx.^2,"all")/numel(gx)^3;
        J12 = sum(gx.*gy,"all")/numel(gx)^3;
        J13 = sum(gx.*gz,"all")/numel(gx)^3;
        J22 = sum(gy.^2,"all")/numel(gx)^3;
        J23 = sum(gy.*gz,"all")/numel(gx)^3;
        J33 = sum(gz.^2,"all")/numel(gx)^3;

        structure_tensor = [J11, J12, J13;
            J12, J22, J23;
            J13, J23, J33];

        % Calculate eigenvalues and eigenvectors
        [eigenvectors, eigenvalues] = eig(structure_tensor);
        eigenvalues = diag(eigenvalues);

        % Find eigenvector corresponding to minimum eigenvalue
        [~, min_idx] = min(eigenvalues);
        direction_vector = eigenvectors(:, min_idx);

        % Normalize and store
        if norm(direction_vector) > 0
            direction_vector = direction_vector / norm(direction_vector);
            % Ensure direction consistency
            if direction_vector(1) < 0
                direction_vector = -direction_vector;
            end
        end

        vector_field(x, y, z, :) = direction_vector;
    end
end

fprintf('Direction vector field calculation completed\n');
    
%% 5. Fiber tracking based on skeleton and direction vector field
fprintf('Performing fiber tracking based on skeleton...\n');
tic;

% Mark visited skeleton voxels
visited = false(size(skeleton_data));

fibers = {}; % Store skeleton voxel coordinates of each fiber
fiber_paths = {}; % Store physical coordinate paths of each fiber
fiber_fitted_paths = {}; % Store fitted curve paths of each fiber

% Get coordinates of all skeleton voxels
[skel_x, skel_y, skel_z] = ind2sub(size(skeleton_data), find(skeleton_data > 0));

for i = 1:length(skel_x)
    start_x = skel_x(i);
    start_y = skel_y(i);
    start_z = skel_z(i);

    % If this skeleton voxel has been visited, skip
    if visited(start_x, start_y, start_z)
        continue;
    end

    % Start bidirectional tracking from current skeleton voxel
    [fiber_voxels_forward,visited] = track_fiber_on_skeleton(start_x, start_y, start_z, ...
        vector_field, visited, skeleton_data, max_step_size, 'forward', max_angle_threshold);

    [fiber_voxels_backward,visited] = track_fiber_on_skeleton(start_x, start_y, start_z, ...
        vector_field, visited, skeleton_data, max_step_size, 'backward', max_angle_threshold);

    % Merge tracking results from both directions
    if ~isempty(fiber_voxels_backward)
        fiber_voxels = [fiber_voxels_backward; [start_x, start_y, start_z]; fiber_voxels_forward];
    else
        fiber_voxels = [[start_x, start_y, start_z]; fiber_voxels_forward];
    end

    % If fiber is long enough, save it
    if size(fiber_voxels, 1) >= min_fiber_length
        fibers{end+1} = fiber_voxels;

        % Convert to physical coordinates
        physical_path = fiber_voxels * voxel_size;
        fiber_paths{end+1} = physical_path;

        % Perform curve fitting on path
        fitted_path = fit_3d_curve(physical_path);
        fiber_fitted_paths{end+1} = fitted_path;

        % Mark all skeleton voxels as visited
        for j = 1:size(fiber_voxels, 1)
            visited(fiber_voxels(j, 1), fiber_voxels(j, 2), fiber_voxels(j, 3)) = true;
        end
    end
end

fprintf('Fiber tracking completed, found %d fibers\n', length(fibers));
toc;
    
%% 6. Calculate fiber length (using fitted curve paths)
fprintf('Calculating fiber path lengths...\n');

fiber_lengths = [];
fiber_fitted_lengths = []; % Length calculated using fitted curves
fiber_directions = [];

for i = 1:length(fiber_paths)
    path = fiber_paths{i};
    fitted_path = fiber_fitted_paths{i};

    % Calculate original path length - accumulate distances between adjacent points
    path_length = 0;
    for j = 1:size(path, 1)-1
        segment_length = norm(path(j+1, :) - path(j, :));
        path_length = path_length + segment_length;
    end

    % Calculate fitted curve length - more accurate length estimation
    fitted_length = 0;
    for j = 1:size(fitted_path, 1)-1
        segment_length = norm(fitted_path(j+1, :) - fitted_path(j, :));
        fitted_length = fitted_length + segment_length;
    end

    fiber_lengths = [fiber_lengths; path_length];
    fiber_fitted_lengths = [fiber_fitted_lengths; fitted_length];

    % Calculate fiber direction (unit vector from start to end point)
    direction = path(end, :) - path(1, :);
    if norm(direction) > 0
        direction = direction / norm(direction);
    end
    fiber_directions = [fiber_directions; direction];
end
    
%% 7. Visualize results

% Figure 1: 3D skeleton path display of first 100 fibers (original paths and fitted curves)
figure('Position', [100, 100, 1400, 600]);

subplot(1, 2, 1);
hold on;

% Determine number of fibers to display (maximum 100)
num_fibers_to_display = min(100, length(fiber_paths));

% Generate colors for each fiber
colors = lines(num_fibers_to_display);

for i = 1:num_fibers_to_display
    path = fiber_paths{i};
    fitted_path = fiber_fitted_paths{i};
    color = colors(i, :);

    % Draw original skeleton path (dashed line)
    plot3(path(:, 1), path(:, 2), path(:, 3), '--', ...
        'Color', color*0.7, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 8);

    % Draw fitted curve path (solid line)
    plot3(fitted_path(:, 1), fitted_path(:, 2), fitted_path(:, 3), ...
        'Color', color, 'LineWidth', 2);

    % Mark endpoints
    scatter3(path(1, 1), path(1, 2), path(1, 3), ...
        80, color, 'filled', 'MarkerEdgeColor', 'k', 'Marker', 's');
    scatter3(path(end, 1), path(end, 2), path(end, 3), ...
        80, color, 'filled', 'MarkerEdgeColor', 'k', 'Marker', 'o');

    % Add length label at path midpoint (using fitted length)
    if size(path, 1) > 5 && fiber_fitted_lengths(i) > 30
        mid_point = round(size(fitted_path, 1)/2);
        text(fitted_path(mid_point, 1), fitted_path(mid_point, 2), fitted_path(mid_point, 3), ...
            sprintf('%.0fμm', fiber_fitted_lengths(i)), 'FontSize', 8, 'Color', color, ...
            'BackgroundColor', 'white', 'Margin', 1);
    end
end

xlabel('X (μm)'); ylabel('Y (μm)'); zlabel('Z (μm)');
title(sprintf('Fiber Path Reconstruction Based on Skeleton and Curve Fitting (Showing First %d Fibers)', num_fibers_to_display));
grid on;
axis equal;

% Add legend explanation
legend({'Original Skeleton', 'Fitted Curve', 'Start Point', 'End Point'}, 'Location', 'best');

% Figure 2: Fiber length distribution comparison (original path length vs fitted curve length)
subplot(1, 2, 2);

if ~isempty(fiber_lengths)
    % Remove outliers
    mean_length = mean(fiber_fitted_lengths);
    std_length = std(fiber_fitted_lengths);
    valid_fitted_lengths = fiber_fitted_lengths(fiber_fitted_lengths < mean_length + 3*std_length);
    valid_raw_lengths = fiber_lengths(fiber_fitted_lengths < mean_length + 3*std_length);

    % Draw comparison of two length distributions
    hold on;
    histogram(valid_raw_lengths, 20, 'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'black', 'FaceAlpha', 0.7);
    histogram(valid_fitted_lengths, 20, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'black', 'FaceAlpha', 0.7);
    xlabel('Fiber Length (μm)');
    ylabel('Number of Fibers');
    title(sprintf('Fiber Length Distribution Comparison (Total %d Fibers)', length(fiber_lengths)));
    grid on;
    legend({'Original Path Length', 'Fitted Curve Length'}, 'Location', 'best');

    % Add statistical information
    text(0.7, 0.9, sprintf('Original Path Average Length: %.2f μm\nFitted Curve Average Length: %.2f μm\nLength Difference: %.2f μm (%.1f%%)', ...
        mean(valid_raw_lengths), mean(valid_fitted_lengths), ...
        mean(valid_raw_lengths - valid_fitted_lengths), ...
        100*mean((valid_raw_lengths - valid_fitted_lengths)./valid_raw_lengths)), ...
        'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
else
    text(0.5, 0.5, 'No fibers detected', 'HorizontalAlignment', 'center', 'Units', 'normalized');
end
    
%% 8. Length difference analysis
if ~isempty(fiber_lengths)
    figure('Position', [100, 600, 1000, 400]);

    % Calculate length differences
    length_differences = fiber_lengths - fiber_fitted_lengths;
    relative_differences = length_differences ./ fiber_lengths * 100;

    subplot(1, 2, 1);
    scatter(fiber_lengths, length_differences, 40, 'filled');
    xlabel('Original Path Length (μm)');
    ylabel('Length Difference (μm)');
    title('Relationship Between Original Path Length and Length Difference');
    grid on;

    subplot(1, 2, 2);
    histogram(relative_differences, 20, 'FaceColor', [0.5, 0.8, 0.3], 'EdgeColor', 'black');
    xlabel('Length Difference Percentage (%)');
    ylabel('Number of Fibers');
    title('Length Difference Percentage Distribution');
    grid on;

    % Add statistical information
    fprintf('\nLength Difference Analysis:\n');
    fprintf('  Average Original Path Length: %.2f μm\n', mean(fiber_lengths));
    fprintf('  Average Fitted Curve Length: %.2f μm\n', mean(fiber_fitted_lengths));
    fprintf('  Average Length Difference: %.2f μm\n', mean(length_differences));
    fprintf('  Average Length Difference Percentage: %.2f%%\n', mean(relative_differences));
    fprintf('  Maximum Length Difference: %.2f μm\n', max(length_differences));
end
    
%% 9. Save results
if ~isempty(fibers)
    fprintf('\n=== Analysis Results ===\n');
    fprintf('Number of fibers detected: %d\n', length(fibers));
    fprintf('Fiber Length Statistics (based on fitted curve length):\n');
    fprintf('  Average: %.2f μm\n', mean(fiber_fitted_lengths));
    fprintf('  Standard Deviation: %.2f μm\n', std(fiber_fitted_lengths));
    fprintf('  Minimum: %.2f μm\n', min(fiber_fitted_lengths));
    fprintf('  Maximum: %.2f μm\n', max(fiber_fitted_lengths));
    fprintf('  Median: %.2f μm\n', median(fiber_fitted_lengths));
    fprintf('  FVF (voxel based): %.2f %%\n',100*(sum(data_3d>0,"all")/numel(data_3d)));
    fprintf('  FVF (analysis based): %.2f %%\n',100*((sum(fiber_lengths)*1e-6*pi*(diameter/2)^2)/(N_x*N_y*N_z*(voxel_size*1e-6)^3)));

    % Save fiber data
    fiber_data = struct();
    for i = 1:length(fibers)
        fiber_data(i).skeleton_voxels = fibers{i};
        fiber_data(i).path = fiber_paths{i};
        fiber_data(i).fitted_path = fiber_fitted_paths{i};
        fiber_data(i).raw_length_um = fiber_lengths(i);
        fiber_data(i).fitted_length_um = fiber_fitted_lengths(i);
        fiber_data(i).direction = fiber_directions(i, :);
    end

    save('fiber_skeleton_curve_results.mat', 'fiber_data', 'fiber_lengths', 'fiber_fitted_lengths', 'fiber_directions', 'skeleton_data');

    % Save as CSV
    results_table = table((1:length(fiber_lengths))', fiber_lengths, fiber_fitted_lengths, ...
        fiber_directions(:, 1), fiber_directions(:, 2), fiber_directions(:, 3), ...
        'VariableNames', {'FiberID', 'RawLength_um', 'FittedLength_um', 'Direction_X', 'Direction_Y', 'Direction_Z'});
    writetable(results_table, 'fiber_skeleton_curve_analysis.csv');

    fprintf('Results saved to fiber_skeleton_curve_results.mat and fiber_skeleton_curve_analysis.csv\n');
else
    fprintf('No fibers detected.\n');
end

%% plot fitted fibers
figure(1); 
nexttile(2); hold on; grid on;
for i=1:min([length(fiber_data) 5000]) % Plot fibers up to a maximum of 5000
    plot3(fiber_data(i).path(:,1),fiber_data(i).path(:,2),fiber_data(i).path(:,3),'k'); plot3(fiber_data(i).fitted_path(:,1),fiber_data(i).fitted_path(:,2),fiber_data(i).fitted_path(:,3),'r--');
end
axis equal tight; view(3);
xlim([0 voxel_size*N_x]); ylim([0 voxel_size*N_y]); zlim([0 voxel_size*N_z]);
title("Traced fiber paths");
drawnow;

%% Functions
function fitted_path = fit_3d_curve(path)
    % Perform curve fitting on 3D path
    
    n_points = size(path, 1);
    
    % If too few points, return original path directly
    if n_points < 4
        fitted_path = path;
        return;
    end
    
    % Parameterize path - use cumulative chord length as parameter
    t = zeros(n_points, 1);
    for i = 2:n_points
        t(i) = t(i-1) + norm(path(i, :) - path(i-1, :));
    end
    
    % Normalize parameter to [0, 1] range
    if t(end) > 0
        t = t / t(end);
    else
        fitted_path = path;
        return;
    end
    
    % Perform spline interpolation for each coordinate component
    try
        % Use spline interpolation
        pp_x = spline(t, path(:, 1));
        pp_y = spline(t, path(:, 2));
        pp_z = spline(t, path(:, 3));
        
        % Generate denser interpolation points
        t_fine = linspace(0, 1, max(50, n_points*3)); % At least 50 points, or 3 times original points
        
        % Calculate points on fitted curve
        x_fine = ppval(pp_x, t_fine);
        y_fine = ppval(pp_y, t_fine);
        z_fine = ppval(pp_z, t_fine);
        
        fitted_path = [x_fine', y_fine', z_fine'];
    catch
        % If spline interpolation fails, use smoothing spline or return original path
        fprintf('Spline interpolation failed, using smoothing\n');
        fitted_path = smooth_3d_path(path);
    end
end

function smoothed_path = smooth_3d_path(path)
    % Use moving average to smooth 3D path
    
    n_points = size(path, 1);
    window_size = min(3, floor(n_points/2)); % Window size
    
    if window_size < 1
        smoothed_path = path;
        return;
    end
    
    smoothed_path = zeros(size(path));
    
    for i = 1:n_points
        % Determine window range
        start_idx = max(1, i - window_size);
        end_idx = min(n_points, i + window_size);
        
        % Calculate average within window
        smoothed_path(i, :) = mean(path(start_idx:end_idx, :), 1);
    end
end

function [fiber_voxels,visited] = track_fiber_on_skeleton(start_x, start_y, start_z, vector_field, visited, skeleton_data, max_step_size, direction, max_angle_threshold)
    % Track fiber along skeleton
    
    fiber_voxels = [];
    current_x = start_x;
    current_y = start_y;
    current_z = start_z;
    
    max_steps = 1000; % Prevent infinite loop
    step_count = 0;
    
    while step_count < max_steps
        visited(current_x, current_y, current_z) = 1;
        % Get direction vector of current point
        current_vector = squeeze(vector_field(current_x, current_y, current_z, :))';
        
        % Find next skeleton voxel
        [next_x, next_y, next_z, found] = find_next_skeleton_voxel(current_x, current_y, current_z, ...
                                                                  current_vector, vector_field, visited, skeleton_data, max_step_size, direction, max_angle_threshold);
        
        if ~found
            break;
        end
        
        % Check direction continuity
        next_vector = squeeze(vector_field(next_x, next_y, next_z, :))';
        if norm(next_vector) > 0.1 && norm(current_vector) > 0.1
            dot_product = dot(current_vector, next_vector);
            if dot_product < 0.3 % Direction sudden change threshold
                break;
            end
        end
        
        % Add to fiber path
        if strcmp(direction, 'forward')
            fiber_voxels = [fiber_voxels; next_x, next_y, next_z];
        else
            fiber_voxels = [next_x, next_y, next_z; fiber_voxels];
        end
        
        % Update current position
        current_x = next_x;
        current_y = next_y;
        current_z = next_z;
        
        step_count = step_count + 1;
    end
end

function [next_x, next_y, next_z, found] = find_next_skeleton_voxel(current_x, current_y, current_z, ...
                                                                  current_vector, vector_field, visited, skeleton_data, max_step_size, track_direction, max_angle_threshold)
    % Find next voxel in skeleton neighborhood
    
    [rows, cols, depths] = size(skeleton_data);
    found = false;
    next_x = current_x;
    next_y = current_y;
    next_z = current_z;
    
    best_score = -inf;
    
    % Search neighborhoods
    for dx = -max_step_size:max_step_size
        for dy = -max_step_size:max_step_size
            for dz = -max_step_size:max_step_size
                % Skip current point
                if dx == 0 && dy == 0 && dz == 0
                    continue;
                end
                
                % Limit step size
                if sqrt(dx^2 + dy^2 + dz^2) > max_step_size
                    continue;
                end
                
                x_candidate = current_x + dx;
                y_candidate = current_y + dy;
                z_candidate = current_z + dz;
                
                % Check boundaries
                if x_candidate < 1 || x_candidate > cols || ...
                   y_candidate < 1 || y_candidate > rows || ...
                   z_candidate < 1 || z_candidate > depths
                    continue;
                end
                
                % Check if it's a skeleton voxel and not visited
                if skeleton_data(x_candidate, y_candidate, z_candidate) > 0 && ...
                   ~visited(x_candidate, y_candidate, z_candidate)
                    
                    % Calculate direction vector of candidate point
                    candidate_vector = squeeze(vector_field(x_candidate, y_candidate, z_candidate, :))';
                    
                    % Check if direction vector is valid
                    if norm(candidate_vector) > 0.1
                        % Calculate angle between current vector and candidate vector
                        dot_product = dot(current_vector, candidate_vector);
                        angle_rad = acos(min(max(dot_product, -1), 1)); % Limit to [-1,1] range
                        angle_deg = rad2deg(angle_rad);
                        
                        % If angle change exceeds threshold, skip this candidate
                        if angle_deg > max_angle_threshold
                            continue;
                        end
                        
                        % Calculate direction consistency score
                        displacement = [dx, dy, dz];
                        displacement = displacement / norm(displacement);
                        
                        % Adjust scoring based on tracking direction
                        if strcmp(track_direction, 'forward')
                            direction_consistency = dot(current_vector, displacement);
                        else
                            direction_consistency = dot(-current_vector, displacement);
                        end
                        
                        vector_consistency = dot(current_vector, candidate_vector);
                        
                        % Use angle as penalty term - smaller angle gets higher score
                        angle_penalty = 1 - (angle_deg / max_angle_threshold);
                        
                        score = direction_consistency + 0.5 * vector_consistency + 0.3 * angle_penalty;
                        
                        if score > best_score
                            best_score = score;
                            next_x = x_candidate;
                            next_y = y_candidate;
                            next_z = z_candidate;
                            found = true;
                        end
                    end
                end
            end
        end
    end
end

function skeleton = simple_skeleton_3d(binary_volume)
    % Simplified 3D skeletonization algorithm
    % Use this simplified version if Image Processing Toolbox is not available
    
    fprintf('Using simplified skeletonization algorithm...\n');
    
    [rows, cols, depths] = size(binary_volume);
    skeleton = false(size(binary_volume));
    
    % Simple distance transform + centerline extraction method
    for z = 1:depths
        slice = binary_volume(:,:,z);
        if any(slice(:))
            % Use 2D skeletonization
            if license('test', 'image_toolbox')
                skeleton_slice = bwmorph(slice, 'skel', Inf);
            else
                % If Image Processing Toolbox is not available, use simplified 2D skeletonization
                skeleton_slice = simple_skeleton_2d(slice);
            end
            skeleton(:,:,z) = skeleton_slice;
        end
    end
    
    % Simple 3D connectivity processing
    for iter = 1:3
        skeleton = clean_skeleton_3d(skeleton);
    end
end

function skeleton_slice = simple_skeleton_2d(binary_slice)
    % Simplified 2D skeletonization algorithm
    [rows, cols] = size(binary_slice);
    skeleton_slice = false(rows, cols);
    
    % Use distance transform to find centerline
    dt_slice = bwdist(~binary_slice);
    max_dt = max(dt_slice(:));
    
    if max_dt > 0
        % Extract distance transform ridge as skeleton
        for i = 2:rows-1
            for j = 2:cols-1
                if binary_slice(i,j)
                    % Check if it's a local maximum
                    if dt_slice(i,j) >= dt_slice(i-1,j) && ...
                       dt_slice(i,j) >= dt_slice(i+1,j) && ...
                       dt_slice(i,j) >= dt_slice(i,j-1) && ...
                       dt_slice(i,j) >= dt_slice(i,j+1)
                        skeleton_slice(i,j) = true;
                    end
                end
            end
        end
    end
end

function cleaned_skeleton = clean_skeleton_3d(skeleton)
    % Clean 3D skeleton, remove isolated points
    cleaned_skeleton = skeleton;
    [rows, cols, depths] = size(skeleton);
    
    for z = 2:depths-1
        for y = 2:rows-1
            for x = 2:cols-1
                if skeleton(y,x,z)
                    % Calculate number of skeleton points in 26-neighborhood
                    neighborhood = skeleton(y-1:y+1, x-1:x+1, z-1:z+1);
                    neighbor_count = sum(neighborhood(:)) - 1; % Subtract self
                    
                    % If too few neighbors, it might be noise
                    if neighbor_count < 2
                        cleaned_skeleton(y,x,z) = false;
                    end
                end
            end
        end
    end
end
