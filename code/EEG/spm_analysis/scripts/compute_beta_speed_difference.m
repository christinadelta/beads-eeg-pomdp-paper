
% Define directories
fast_beta_dir = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/tfr_contr_smoothed_beta_fast/';  % Update this path
slow_beta_dir = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/tfr_contr_smoothed_beta_slow/';  % Update this path
output_dir = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/tfr_contrasts_speed_diff/'; % Output directory

% Ensure output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Number of subjects
nsubs = 40; 

%%
% Loop over subjects
for sub = 1:nsubs

    % Define subject ID with leading zeros
    sub_id = sprintf('%02d', sub); % Ensures sub_01, sub_02, ..., sub_40
    
    % Define input file paths
    fast_beta_img = fullfile(fast_beta_dir, sprintf('ssub_%s_condition_urndraw.nii', sub_id));
    slow_beta_img = fullfile(slow_beta_dir, sprintf('ssub_%s_condition_urndraw.nii', sub_id));
    
    % Define output file path
    output_img = fullfile(output_dir, sprintf('sub_%s_condition_urndraw_diff.nii', sub_id));

    % Check if both input images exist
    if exist(fast_beta_img, 'file') && exist(slow_beta_img, 'file')
        
        % Set up SPM ImCalc
        matlabbatch = {};
        matlabbatch{1}.spm.util.imcalc.input = {fast_beta_img; slow_beta_img};
        matlabbatch{1}.spm.util.imcalc.output = output_img;
        matlabbatch{1}.spm.util.imcalc.outdir = {output_dir};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1 - i2'; % Fast - Slow
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16; % 16-bit float

        % Run batch
        spm_jobman('run', matlabbatch);
        fprintf('Computed difference image for sub-%s\n', sub_id);
    else
        fprintf('Skipping sub-%s: missing input files\n', sub_id);
    end
end
