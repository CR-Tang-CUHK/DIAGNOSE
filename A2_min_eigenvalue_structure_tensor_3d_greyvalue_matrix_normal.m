%% Read Excel files and construct 3D matrices
filename = '3D_greyscale_voxel_table_of_fiber_masks_top_right_3.4um_cut a section.xls'; % Input Excel filename
[~, sheets] = xlsfinfo(filename);
num_sheets = length(sheets);

% Read the first worksheet to obtain the matrix size
first_sheet = xlsread(filename, sheets{1});
[m, n] = size(first_sheet);
data_3D = zeros(m, n, num_sheets); % Initialize 3D matrix

% Fill 3D matrix
for k = 1:num_sheets
    sheet_data = xlsread(filename, sheets{k});
    data_3D(:,:,k) = sheet_data;
end

%% Data normalization
% Find the minimum and maximum values of the entire 3D matrix
min_val = min(data_3D(:));
max_val = max(data_3D(:));

% Avoid dividing by zero
if abs(max_val - min_val) < eps
    % If all values are the same, normalize to 0.5
    data_3D_normalized = 0.5 * ones(size(data_3D));
else
    %Linear normalization to the range of [0,1]
    data_3D_normalized = (data_3D - min_val) / (max_val - min_val);
end

%% Calculate 3D gradient (using center difference) - using normalized data
[gx, gy, gz] = gradient(data_3D_normalized);

%% Initialization result matrix
vx_matrix = zeros(m, n, num_sheets); % Unit vector x component
vy_matrix = zeros(m, n, num_sheets); % Unit vector y component
vz_matrix = zeros(m, n, num_sheets); % Unit vector z component

angle_xy = nan(m, n, num_sheets); % Angle between XY plane projection and X-axis (radians), initialized to NaN
angle_xz = nan(m, n, num_sheets); % Angle between XZ plane projection and X-axis (radians)
angle_yz = nan(m, n, num_sheets); % Angle between YZ plane projection and Y-axis (radians)

%% Mark calculated points (initialized to false)
computed_mask = false(m, n, num_sheets);

%% Calculate structural tensors, eigenvectors, and angles
window_size = 3; % Kernel size
half_win = floor(window_size/2);

% Traverse each non-boundary point
for i = (1+half_win):(m-half_win)
    for j = (1+half_win):(n-half_win)
        for k = (1+half_win):(num_sheets-half_win)
            % Only calculate non-zero points in normalized data
            if data_3D_normalized(i,j,k) == 0
                continue;
            end
            
            % Mark this point to be calculated
            computed_mask(i,j,k) = true;
            
            % Extract gradients within the neighborhood
            Ix = gx(i-half_win:i+half_win, j-half_win:j+half_win, k-half_win:k+half_win);
            Iy = gy(i-half_win:i+half_win, j-half_win:j+half_win, k-half_win:k+half_win);
            Iz = gz(i-half_win:i+half_win, j-half_win:j+half_win, k-half_win:k+half_win);
            
            % Calculate the structural tensor
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
            eigenvec = V(:, min_idx); % Corresponding feature vector
            
            % Normalize to unit vector
            norm_vec = norm(eigenvec);
            if norm_vec > eps
                eigenvec = eigenvec / norm_vec;
            else
                eigenvec = [0; 0; 0]; % Avoid zero vectors
            end
            
            % Store vector components
            vx_matrix(i,j,k) = eigenvec(1);
            vy_matrix(i,j,k) = eigenvec(2);
            vz_matrix(i,j,k) = eigenvec(3);
            
            % Calculate angle (only when the vector is non-zero)
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

%% Save vector components to Excel
output_vector_x = 'vector_components_x.xlsx';
output_vector_y = 'vector_components_y.xlsx';
output_vector_z = 'vector_components_z.xlsx';

for k = 1:num_sheets
    xlswrite(output_vector_x, vx_matrix(:,:,k), k);
    xlswrite(output_vector_y, vy_matrix(:,:,k), k);
    xlswrite(output_vector_z, vz_matrix(:,:,k), k);
end

%% Save angle to Excel
output_angle_xy = 'angles_xy.xlsx';
output_angle_xz = 'angles_xz.xlsx';
output_angle_yz = 'angles_yz.xlsx';

for k = 1:num_sheets
    % Replace NaN with 0 for Excel to save
    angle_xy_sheet = angle_xy(:,:,k);
    angle_xy_sheet(isnan(angle_xy_sheet)) = 0;
    
    angle_xz_sheet = angle_xz(:,:,k);
    angle_xz_sheet(isnan(angle_xz_sheet)) = 0;
    
    angle_yz_sheet = angle_yz(:,:,k);
    angle_yz_sheet(isnan(angle_yz_sheet)) = 0;
    
    xlswrite(output_angle_xy, angle_xy_sheet, k);
    xlswrite(output_angle_xz, angle_xz_sheet, k);
    xlswrite(output_angle_yz, angle_yz_sheet, k);
end

%% Draw an angle distribution histogram (only counting the points that have been actually calculated)
% Extract effective angles (only including calculated and non NaN values)
valid_angles_xy = angle_xy(computed_mask & ~isnan(angle_xy));
valid_angles_xz = angle_xz(computed_mask & ~isnan(angle_xz));
valid_angles_yz = angle_yz(computed_mask & ~isnan(angle_yz));

% Convert radius to degrees 
valid_angles_xy_deg = rad2deg(valid_angles_xy);
valid_angles_xz_deg = rad2deg(valid_angles_xz);
valid_angles_yz_deg = rad2deg(valid_angles_yz);

% Draw angle distribution histogram 
figure;
subplot(3,1,1);
histogram(valid_angles_xy_deg, 50, 'BinLimits', [-180, 180], 'FaceColor', [0.2 0.6 0.9]);
title('Orientation projection on XY plane');
xlabel('Angle (deg)');
ylabel('Voxel number');
grid on;

subplot(3,1,2);
histogram(valid_angles_xz_deg, 50, 'BinLimits', [-180, 180], 'FaceColor', [0.9 0.5 0.2]);
title('Orientation projection on XZ plane');
xlabel('Angle (deg)');
ylabel('Voxel number');
grid on;

subplot(3,1,3);
histogram(valid_angles_yz_deg, 50, 'BinLimits', [-180, 180], 'FaceColor', [0.4 0.8 0.4]);
title('Orientation projection on XZ plane');
xlabel('Angle (deg)');
ylabel('Voxel number');
grid on;

% Adjust histogram layout
set(gcf, 'Position', [100, 100, 800, 900]);

%% print messages
fprintf('统计信息：\n');
fprintf('原始数据范围: [%.4f, %.4f]\n', min_val, max_val);
fprintf('归一化后数据范围: [%.4f, %.4f]\n', min(data_3D_normalized(:)), max(data_3D_normalized(:)));
fprintf('XY平面有效角度数量: %d\n', numel(valid_angles_xy_deg));
fprintf('XZ平面有效角度数量: %d\n', numel(valid_angles_xz_deg));
fprintf('YZ平面有效角度数量: %d\n', numel(valid_angles_yz_deg));