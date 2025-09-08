%% Create required directories and define paths %%

% Add SPM12 to the matlab path if needed:
% addpath /Users/christinadelta/software/neuroscience/spm/spm12 % change this to your
% own SPM12 path
% savepath

% create spm directories that will be used for storing the output MEEG
% objects mat files and images. MEEG objects are stored in the output dir,
% .mat files in the jobs dir and images are stored in directories that are
% created automatically during conversion of the MEEG objects to .nii
% files.

% clear workspace
clear all
clc

%% 

basedir     = pwd;
% spmdir      = fullfile(basedir, 'spmDir');
% addpath(genpath(basedir));

% if output and jobs directories do not exist, create them:
% outDir      = fullfile(basedir, 'output');
% jobsDir     = fullfile(basedir, 'jobs');
% cropDir     = fullfile(basedir, 'cropped');
outDir      = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/output';
jobsDir     = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/jobs';
cropDir     = '/Volumes/beadsData/optimal_stopping_data/beads/spm_analysis/cropped_v3';

if ~exist(outDir, 'dir') && ~exist(jobsDir, 'dir')
%     mkdir(outDir)
%     mkdir(jobsDir)
    mkdir(cropDir)
end

% add directories to the path
addpath(genpath(jobsDir));
addpath(genpath(outDir));
addpath(genpath(fullfile(basedir, 'scripts')));
addpath(genpath(fullfile(basedir, 'utilities')));

% set data path
datadir         = '/Volumes/beadsData/optimal_stopping_data/data/beads/eeg/';  % change to your own data path 
subs            = dir(fullfile(datadir, '*sub*'));                      % how many subject folders?
nsubs           = length(subs);                                         % how many subjects?

taskname        = 'beads';
blocks          = 4;
% conditions      = 2;
choices         = 2;
analysest       = 2; % analyses types (erp & tf)
nconstrasts     = 5;
contrastpref    = {'wud_' 'wde_' 'wi_' 'wu_' 'wd_'};
condtypes       = {'easydraw' 'easyurn' 'diffdraw' 'diffurn'};

%% Run TF analysis for slow and fast beta

% these steps are used to preprocess the raw .bdf files, for each
% participant. Given that we will run evoked analysis and time-frequency
% analysis, we will low-pass filter the data in two ways, one for evoked
% analysis and one for TF analysis. Both MEEG objects will be preprocessed
% up until merging, however, only the ERPs MEEG object will undergo
% artefact rejection (the last step of preprocessing).

% init spm 
spm('defaults', 'eeg');

%% loop over subejcts and compute power for slow and fast beta

for sub = 1:nsubs

    % get subject's output path:
    subout          = fullfile(outDir, sprintf('sub-%02d', sub));

    % before averaging conditions, we will need to average the TF
    % object over slow beta frequency (as this is the focus of the analysis)

    S                       = [];
    S.D                     = fullfile(subout, sprintf('rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
    S.freqwin               = 13:20; % average over beta
    S.prefix                = 'Pslow_beta_';
    D                       = spm_eeg_avgfreq(S);

    % now average power of fast beta
    S                       = [];
    S.D                     = fullfile(subout, sprintf('rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
    S.freqwin               = 20:30; % average over beta
    S.prefix                = 'Pfast_beta_';
    D                       = spm_eeg_avgfreq(S);


    %% average over conditions

    % average slow beta power over time 
    S                       = [];
    S.D                     = fullfile(subout, sprintf('Pslow_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
    S.robust.ks             = 3;
    S.robust.bycondition    = true;
    S.robust.savew          = false;
    S.robust.removebad      = false;
    S.circularise           = false;
    S.prefix                = 'm';
    D                       = spm_eeg_average(S);

    % average fast beta power over time
    S                       = [];
    S.D                     = fullfile(subout, sprintf('Pfast_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
    S.robust.ks             = 3;
    S.robust.bycondition    = true;
    S.robust.savew          = false;
    S.robust.removebad      = false;
    S.circularise           = false;
    S.prefix                = 'm';
    D                       = spm_eeg_average(S);

    %% compute contrasts

    % compute contrasts for slow beta
    S                       = [];
    S.D                     = fullfile(subout, sprintf('mPslow_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
    S.object                = spm_eeg_load(S.D);
    tmp                     = S.object.conditions;
    
    % estimate TF power contrasts
    beadsEstimateContrast(S, tmp, sub)

    % compute contrasts for fast beta
    S                       = [];
    S.D                     = fullfile(subout, sprintf('mPfast_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
    S.object                = spm_eeg_load(S.D);
    tmp                     = S.object.conditions;
    
    % estimate TF power contrasts
    beadsEstimateContrast(S, tmp, sub)

    %% convert contrasts yo 3D volumes

    % convert the 5 slow beta contrast MEEG objects into .nii files
    for c = 1:nconstrasts

        S               = [];
        S.D             = fullfile(subout, sprintf('%smPslow_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', contrastpref{c}, sub));
        S.mode          = 'scalp x time';
        S.conditions    = {};
        S.channels      = 'EEG';
        S.timewin       = [-Inf Inf];
        S.freqwin       = [-Inf Inf];
        S.prefix        = '';

        D               = spm_eeg_convert2images(S);

    end

    % convert the 5 fast beta contrast MEEG objects into .nii files
    for c = 1:nconstrasts

        S               = [];
        S.D             = fullfile(subout, sprintf('%smPfast_beta_rtf_cetfrfdfMspmeeg_sub_%02d_beads_block_01.mat', contrastpref{c}, sub));
        S.mode          = 'scalp x time';
        S.conditions    = {};
        S.channels      = 'EEG';
        S.timewin       = [-Inf Inf];
        S.freqwin       = [-Inf Inf];
        S.prefix        = '';

        D               = spm_eeg_convert2images(S);

    end


end % end of subejcts loop


