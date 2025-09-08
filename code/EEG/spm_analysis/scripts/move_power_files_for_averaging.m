
%%
clc
clear all

%%

% Define directories
source_dir = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/output/';  % Base directory containing participant folders
% dest_slow = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/averages/TFRs_beta_slow/'; % Destination for slow beta files
% dest_fast = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/averages/TFRs_beta_fast/'; % Destination for fast beta files
dest_alpha = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/averages/alpha_averaged';
dest_gamma = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/averages/gamma_averaged';

% Ensure destination folders exist
if ~exist(dest_alpha, 'dir')
    mkdir(dest_alpha);
end
if ~exist(dest_gamma, 'dir')
    mkdir(dest_gamma);
end

%%
% Get list of all participant subfolders
sub_folders = dir(fullfile(source_dir, 'sub-*')); % Adjust pattern if needed

% Loop through each participant folder
for i = 1:length(sub_folders)
    sub_folder = fullfile(source_dir, sub_folders(i).name); % Full path to participant folder
    
    % Define expected file names
%     slow_beta_file = sprintf('mPslow_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', i);
%     fast_beta_file = sprintf('mPfast_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', i);

    alpha_file = sprintf('wud_mPalpha_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', i);
    gamma_file = sprintf('wud_mPgamma_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', i);
    
    % Check if slow beta file exists and move it
    if exist(fullfile(sub_folder, alpha_file), 'file')
        movefile(fullfile(sub_folder, alpha_file), fullfile(dest_alpha, alpha_file));
        fprintf('Moved %s to %s\n', alpha_file, dest_alpha);
    else
        fprintf('File not found: %s\n', alpha_file);
    end
    
    % Check if fast beta file exists and move it
    if exist(fullfile(sub_folder, gamma_file), 'file')
        movefile(fullfile(sub_folder, gamma_file), fullfile(dest_gamma, gamma_file));
        fprintf('Moved %s to %s\n', gamma_file, dest_gamma);
    else
        fprintf('File not found: %s\n', gamma_file);
    end
end

fprintf('File moving complete.\n');
