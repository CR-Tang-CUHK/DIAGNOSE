clear all;
clc;
%% Read Excel file and build 3D matrix
filename = '3D_greyscale_voxel.xls'; % Input Excel filename
[~, sheets] = xlsfinfo(filename);
num_sheets = length(sheets);

% Read first sheet to get matrix dimensions
first_sheet = xlsread(filename, sheets{1});
[m, n] = size(first_sheet);
data_3D = zeros(m, n, num_sheets); % Initialize 3D matrix

% Fill 3D matrix
for k = 1:num_sheets
    sheet_data = xlsread(filename, sheets{k});
    data_3D(:,:,k) = sheet_data;
end

%% Data normalization
min_val = min(data_3D(:));
max_val = max(data_3D(:));

if abs(max_val - min_val) < eps
    data_3D_normalized = 0.5 * ones(size(data_3D));
else
    data_3D_normalized = (data_3D - min_val) / (max_val - min_val);
end

%% Calculate 3D gradient (using central difference) - using normalized data
[gx, gy, gz] = gradient(data_3D_normalized);

%% Initialize result matrices
vx_matrix = zeros(m, n, num_sheets); % Unit vector x-component
vy_matrix = zeros(m, n, num_sheets); % Unit vector y-component
vz_matrix = zeros(m, n, num_sheets); % Unit vector z-component

angle_xy = nan(m, n, num_sheets); % Angle between XY plane projection and X-axis (radians), initialized as NaN
angle_xz = nan(m, n, num_sheets); % Angle between XZ plane projection and X-axis (radians)
angle_yz = nan(m, n, num_sheets); % Angle between YZ plane projection and Y-axis (radians)

%% Mark computed points (initialized as false)
computed_mask = false(m, n, num_sheets);

%% Calculate structure tensor, eigenvectors and angles
window_size = 3; % Neighborhood window size
half_win = floor(window_size/2);

% Iterate through each non-boundary point
for i = (1+half_win):(m-half_win)
    for j = (1+half_win):(n-half_win)
        for k = (1+half_win):(num_sheets-half_win)
            % Only process points with non-zero normalized data
            if data_3D_normalized(i,j,k) == 0
                continue;
            end
            
            % Mark this point for computation
            computed_mask(i,j,k) = true;
            
            % Extract gradients in neighborhood
            Ix = gx(i-half_win:i+half_win, j-half_win:j+half_win, k-half_win:k+half_win);
            Iy = gy(i-half_win:i+half_win, j-half_win:j+half_win, k-half_win:k+half_win);
            Iz = gz(i-half_win:i+half_win, j-half_win:j+half_win, k-half_win:k+half_win);
            
            % Calculate structure tensor
            J11 = mean(Ix(:).^2);
            J12 = mean(Ix(:).*Iy(:));
            J13 = mean(Ix(:).*Iz(:));
            J22 = mean(Iy(:).^2);
            J23 = mean(Iy(:).*Iz(:));
            J33 = mean(Iz(:).^2);
            
            J = [J11, J12, J13;
                 J12, J22, J23;
                 J13, J23, J33];
            
            % Calculate eigenvalues and eigenvectors
            [V, D] = eig(J);
            eigenvalues = diag(D);
            [~, min_idx] = min(eigenvalues); % Minimum eigenvalue index
            eigenvec = V(:, min_idx); % Corresponding eigenvector
            
            % Normalize to unit vector
            norm_vec = norm(eigenvec);
            if norm_vec > eps
                eigenvec = eigenvec / norm_vec;
            else
                eigenvec = [0; 0; 0]; % Avoid zero vector
            end
            
            % Store vector components
            vx_matrix(i,j,k) = eigenvec(1);
            vy_matrix(i,j,k) = eigenvec(2);
            vz_matrix(i,j,k) = eigenvec(3);
            
            % Calculate angles (only when vector is non-zero)
            if norm_vec > eps
                % Angle between XY plane projection and X-axis
                if abs(eigenvec(1)) > eps || abs(eigenvec(2)) > eps
                    angle_xy(i,j,k) = atan2(eigenvec(2), eigenvec(1));
                end
                
                % Angle between XZ plane projection and X-axis
                if abs(eigenvec(1)) > eps || abs(eigenvec(3)) > eps
                    angle_xz(i,j,k) = atan2(eigenvec(3), eigenvec(1));
                end
                
                % Angle between YZ plane projection and Y-axis
                if abs(eigenvec(2)) > eps || abs(eigenvec(3)) > eps
                    angle_yz(i,j,k) = atan2(eigenvec(3), eigenvec(2));
                end
            end
        end
    end
end

%% Process angle data and create polar histograms
% Extract valid angles (only computed and non-NaN values)
valid_angles_xy = angle_xy(computed_mask & ~isnan(angle_xy));
valid_angles_xz = angle_xz(computed_mask & ~isnan(angle_xz));
valid_angles_yz = angle_yz(computed_mask & ~isnan(angle_yz));

% Convert to degrees and map to 0-360 degree range
valid_angles_xy_deg = mod(rad2deg(valid_angles_xy), 360);
valid_angles_xz_deg = mod(rad2deg(valid_angles_xz), 360);
valid_angles_yz_deg = mod(rad2deg(valid_angles_yz), 360);

% Merge symmetric angles (e.g., 90 degrees and 270 degrees)
% Map angles to 0-180 degree range and record counts
function [theta_bins, counts] = merge_symmetric_angles(angles_deg, num_bins)
    % Map angles to 0-180 degree range
    mapped_angles = mod(angles_deg, 180);
    
    % Create histogram
    [counts, edges] = histcounts(mapped_angles, num_bins, 'BinLimits', [0, 180]);
    
    % Calculate center angle for each bin
    bin_width = edges(2) - edges(1);
    theta_bins = edges(1:end-1) + bin_width/2;
end

% Create merged histograms for each angle set
num_bins = 36; % One bin per 5 degrees (180/36=5)
[theta_bins_xy, counts_xy] = merge_symmetric_angles(valid_angles_xy_deg, num_bins);
[theta_bins_xz, counts_xz] = merge_symmetric_angles(valid_angles_xz_deg, num_bins);
[theta_bins_yz, counts_yz] = merge_symmetric_angles(valid_angles_yz_deg, num_bins);

% Create full angle data (0-360 degrees)
function [full_angles, full_counts] = create_full_polar_data(theta_bins, counts)
    % Create full 0-360 degree angle array
    full_angles = [theta_bins, theta_bins + 180];
    
    % Create full count array
    full_counts = [counts, counts];
end

% Create full data for each angle set
[full_angles_xy, full_counts_xy] = create_full_polar_data(theta_bins_xy, counts_xy);
[full_angles_xz, full_counts_xz] = create_full_polar_data(theta_bins_xz, counts_xz);
[full_angles_yz, full_counts_yz] = create_full_polar_data(theta_bins_yz, counts_yz);

% Plot polar histograms
figure;

% Angle between XY plane projection and X-axis
subplot(1,3,1);
% Create angle data array, repeating each angle according to its count
angles_data_xy = [];
for i = 1:length(full_angles_xy)
    angles_data_xy = [angles_data_xy, repmat(full_angles_xy(i), 1, round(full_counts_xy(i)))];
end
polarhistogram(deg2rad(angles_data_xy), 72, 'FaceColor', [0.2 0.6 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('Fiber orientation on XY plane');
set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');

% Angle between XZ plane projection and X-axis
subplot(1,3,2);
% Create angle data array, repeating each angle according to its count
angles_data_xz = [];
for i = 1:length(full_angles_xz)
    angles_data_xz = [angles_data_xz, repmat(full_angles_xz(i), 1, round(full_counts_xz(i)))];
end
polarhistogram(deg2rad(angles_data_xz), 72, 'FaceColor', [0.9 0.5 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('Fiber orientation on XZ plane');
set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');

% Angle between YZ plane projection and Y-axis
subplot(1,3,3);
% Create angle data array, repeating each angle according to its count
angles_data_yz = [];
for i = 1:length(full_angles_yz)
    angles_data_yz = [angles_data_yz, repmat(full_angles_yz(i), 1, round(full_counts_yz(i)))];
end
polarhistogram(deg2rad(angles_data_yz), 72, 'FaceColor', [0.4 0.8 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('Fiber orientation on YZ plane');
set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');

% Adjust layout
set(gcf, 'Position', [100, 100, 800, 900]);

%% Display statistical information
fprintf('Statistical Information:\n');
fprintf('Original data range: [%.4f, %.4f]\n', min_val, max_val);
fprintf('Normalized data range: [%.4f, %.4f]\n', min(data_3D_normalized(:)), max(data_3D_normalized(:)));
fprintf('Number of valid angles in XY plane: %d\n', numel(valid_angles_xy_deg));
fprintf('Number of valid angles in XZ plane: %d\n', numel(valid_angles_xz_deg));
fprintf('Number of valid angles in YZ plane: %d\n', numel(valid_angles_yz_deg));