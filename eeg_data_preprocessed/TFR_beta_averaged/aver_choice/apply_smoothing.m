% ========================================================================
% Apply Gaussian Smoothing to Time-Frequency Data in an SPM MEEG Object
% ========================================================================
% This script:
% 1. Loads the SPM MEEG object containing time-frequency data.
% 2. Extracts the power data.
% 3. Applies Gaussian smoothing to each channel.
% 4. Saves the smoothed data back into the MEEG object properly.

% ----------------
% User Parameters
% ----------------
% input_file = 'aver_tfr_choices_frequencies.mat'; % 
input_file = 'choices_aver_tfr_freq.mat';
output_prefix = 'smoothed_';  % Name for the smoothed output file

%%
% ----------------
% Step 1: Load SPM MEEG Time-Frequency Data
% ----------------
D = spm_eeg_load(input_file); % Load the SPM EEG/MEG object
tf_data = D(:,:,:,:); % Extract TF power data (Channels x Frequencies x Time)

%%
% ----------------
% Step 2: Define Gaussian Smoothing Kernel
% ----------------
smooth_kernel = fspecial('gaussian',[5 5], 1.5); % Adjust kernel size & sigma

%%
% ----------------
% Step 3: Apply Smoothing to Each Channel
% ----------------
for ch = 1:size(tf_data, 1)  % Loop through channels
    tf_data(ch, :, :) = imfilter(squeeze(tf_data(ch, :, :)), smooth_kernel, 'replicate');
end

%%
% ----------------
% Step 4: Create a New SPM MEEG Object (Avoid Overwriting)
% ----------------
Dnew = clone(D, [output_prefix D.fname]); % Create a new copy with modified name
Dnew(:,:,:,:) = tf_data; % Store smoothed TF data in the new object

%%
% ----------------
% Step 5: Save the New Smoothed MEEG Object
% ----------------
save(Dnew); % Save the new file

%%
% ----------------
% Completion Message
% ----------------
fprintf('Smoothing complete! Smoothed TF file saved as: %s\n', Dnew.fname);