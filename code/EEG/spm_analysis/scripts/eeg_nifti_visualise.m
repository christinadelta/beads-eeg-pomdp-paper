% EEG visualisation plots 
% BEADS task June 2024

%% set figure-docking as default 

clear all
clc

% set(0,'DefaultFigureWindowStyle','docked')

%% Add SPM to MATLAB path and load nifti files

tmp_path        = pwd; % use parent directory for now
addpath(tmp_path);
 

% Load the spmT_001.nii file
V_tstat         = spm_vol(fullfile(tmp_path,'spmT_0001.nii'));
spmT_data       = spm_read_vols(V_tstat);

% Load the con_0001.nii file to get the contrast data
V_con           = spm_vol(fullfile(tmp_path,'con_0001.nii'));
con_data        = spm_read_vols(V_con);

% Load mask volume 
V_mask          = spm_vol('mask.nii');  % Adjust to your file path
mask_data       = spm_read_vols(V_mask);

% Ensure the mask is binary (in case it's not)
mask_data       = mask_data > 0;

%% Identify statisticaly significant tps in the spmT data

% Threshold for statistical significance
threshold               = 3.31;

% Initialize array to store significant time points
significant_time_points = [];

% Loop through each time slice in the 3rd dimension
for t = 1:size(spmT_data, 3)

    % Extract the T-map slice for the current time point
    t_slice             = spmT_data(:, :, t);
    
    % Check if there are any values above the threshold
    if any(t_slice(:) > threshold)

        % Add this time point to the list
        significant_time_points = [significant_time_points, t];
    end
end

% Display significant time points
disp('Significant Time Points:');
disp(significant_time_points);

%% Define the epoch parameters

epoch_start                 = -500;  % Start of the epoch in ms
epoch_end                   = 800;
% time_resolution             = 1000 / 256;  % Time resolution in ms (3.90625 ms/sample)
time_resolution             = 1000 / 128;

% Define significant indices (already calculated earlier)
significant_indices         = significant_time_points; 

% Calculate the time points for the significant indices
significant_time_points_ms  = epoch_start + (significant_indices - 1) * time_resolution;

% Display the results
disp('Significant Time Points (in ms):');
disp(significant_time_points_ms);

%% Volume-Based Visualization 

% Initialize the masked T-stat data matrix
masked_spmT_data = zeros(size(spmT_data));

% Loop through each slice to apply the mask
for slice = 1:size(spmT_data, 3)
    masked_spmT_data(:, :, slice) = spmT_data(:, :, slice) .* mask_data(:, :, slice);
end

% Calculate global max absolute value for consistent color bar limits
global_max = max(abs(masked_spmT_data(:)));

% Define significant clusters
sig_clusters = masked_spmT_data > threshold;

% Plot a single slice with significant regions
slice_idx = significant_indices(47);  % Example: Pick a significant time point
figure;
imagesc(masked_spmT_data(:, :, slice_idx));  % Plot the masked data slice
caxis([-global_max, global_max]);  % Consistent color bar limits
hold on;
contour(sig_clusters(:, :, slice_idx), 1, 'r', 'LineWidth', 2);  % Overlay significant regions
colorbar;
title(['Significant Clusters at Time Point: ', num2str(significant_time_points_ms(47)), ' ms']);

% Adjust figure layout
set(gca, 'FontSize', 12);

%%
% Visualize the mask itself
figure;
imagesc(mask_data(:, :, 8));  % Plot a slice of the mask
title(['Mask Slice at Index ', num2str(slice_idx)]);
colorbar;

%%
% Calculate the difference between consecutive slices in the mask
mask_diff = diff(mask_data, 1, 3);  % Compute difference along the 3rd dimension

% Visualize the difference for a few slices
figure;
for i = 1:5  % Adjust the range as needed
    subplot(1, 5, i);
    imagesc(mask_diff(:, :, i));  % Plot the difference between slices
    title(['Slice ', num2str(i)]);
    colorbar;
end

% Check if the mask is identical across all slices
is_same = all(mask_diff(:) == 0);
if is_same
    disp('The mask is the same across all slices.');
else
    disp('The mask varies across slices.');
end

%% Region-of-Interest (ROI) Visualisation

% Focus on an ROI (e.g., central grid region)
roi_x = 10:20;  % Example: X-coordinates for the ROI
roi_y = 10:20;  % Example: Y-coordinates for the ROI

% Extract ROI data
roi_data = spmT_data(roi_x, roi_y, slice_idx);

% Visualize ROI
figure;
imagesc(roi_data);
colorbar;
title('ROI Data');

%%

% Create an animation of significant regions over time
figure;
for t = significant_indices
    imagesc(con_data(:, :, t));
    hold on;
    contour(sig_clusters(:, :, t), 1, 'r', 'LineWidth', 2);  % Overlay significant regions
    colorbar;
    title(['Time: ', num2str(-500 + (t-1) * (1000/128)), ' ms']);
    pause(0.5);  % Pause for animation
end