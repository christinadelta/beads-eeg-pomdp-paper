
% The tests looks whether there is a relationship between the difference in
% the number of draws for humans - CNM and nll

% the hypothesis is that: if the model fails to predict human sampling
% (with either higher or lower draws), the nll is weak 

clear all 
clc 

%% Set Figure-Docking as Default
set(0, 'DefaultFigureWindowStyle', 'docked')

%% Load and Organise Data 

% load human draws, cnm draws and nll values for correlations
load('human_cnm_draws_nll_data.mat')


% re-orginise the data
human_draws_08  = easy_avdraws;
human_draws_06  = diff_avdraws;

cnm_draws_08    = fit_samples_v2(:,1);
cnm_draws_06    = fit_samples_v2(:,2);

nll_08          = allsub_ll_v2(:,1);
nll_06          = allsub_ll_v2(:,2);

% Number of participants
nsubs           = length(human_draws_08);  % Assuming 40 participants

%% Compute Absolute Differences in Sampling Rates

diff_08         = abs(human_draws_08 - cnm_draws_08);
diff_06         = abs(human_draws_06 - cnm_draws_06);


%% Statistical Analysis: Correlations
% Pearson correlation for 0.8 condition
[corr_08, p_08] = corr(diff_08, nll_08);

% Pearson correlation for 0.6 condition
[corr_06, p_06] = corr(diff_06, nll_06);

% Display correlations
disp('Correlation between Sampling Difference and NLL:');
disp(['0.8 Condition: r = ', num2str(corr_08), ', p = ', num2str(p_08)]);
disp(['0.6 Condition: r = ', num2str(corr_06), ', p = ', num2str(p_06)]);

%% Visualization: Scatter Plots with Line of Best Fit (Supplementary figure 6)
figure;
subplot(1,2,1);
scatter(diff_08, nll_08, 'filled');
title('0.8 Condition: Sampling Difference vs NLL');
xlabel('Absolute Difference in Draws (Human - CNM)');
ylabel('NLL');
grid on;
lsline;  % Adds line of best fit

subplot(1,2,2);
scatter(diff_06, nll_06, 'filled');
title('0.6 Condition: Sampling Difference vs NLL');
xlabel('Absolute Difference in Draws (Human - CNM)');
ylabel('NLL');
grid on;
lsline;  % Adds line of best fit
fontsize(gcf, 20, "points");