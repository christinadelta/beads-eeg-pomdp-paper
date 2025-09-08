% parameter recovery of the beads task using Bruno's adapted code
% Version 1 June 2025

% Parameter recovery method: Random Search 
%% Housekeeping Commands
clear all
clc

%% Set Figure-Docking as Default
set(0, 'DefaultFigureWindowStyle', 'docked')

%% Define Initial Variables for Simulations
totaltrials         = 52;
conditions          = 2;
simvars.ntrials     = totaltrials;
simvars.maxDraws    = 10;
simvars.qvals       = [0.8 0.6];
simvars.correct     = 10;
simvars.error       = -10;
simvars.difference  = -20;
simvars.conditions  = conditions;
simvars.contrials   = totaltrials / conditions;

%% Define total samples for random search

total_samples   = 2000;
num_starts      = 5; % based on computational resources

%% Preallocate arrays
true_cs         = zeros(total_samples, 1);
true_beta       = zeros(total_samples, 1);
fit_cs          = zeros(total_samples, conditions);
fit_beta        = zeros(total_samples, conditions);
fit_bestNLL     = zeros(total_samples, conditions);
fit_NLL         = zeros(total_samples, conditions);
sim_draws       = zeros(total_samples, conditions);
fit_draws       = zeros(total_samples, conditions);

%% Generate true parameters and perform recovery

for sample = 1:total_samples

    % Sample true parameters from full ranges
    true_cs(sample)     = unifrnd(-2, 0);
    true_beta(sample)   = unifrnd(0, 10);
    
    param               = [true_cs(sample), true_beta(sample)];
    
    for cond = 1:conditions

        % Condition-wise simulation
        simvars.cond                = cond;
        simvars.thisq               = simvars.qvals(cond);
        
        % Simulate data
        simoutput                   = sim_POMDP_Beads(simvars, param);
        
        % store simulation results
        sim_draws(sample, cond)     = simoutput.avsamples;
        sim_sequences{cond}         = simoutput.simsequences;
        sim_choices{cond}           = simoutput.simchoicevec;
        simvars.urntype             = simoutput.simurns;
        
        % Define objective function
        obFunc          = @(x) mdpBeads_paramRec([x(1) x(2)], simvars, simoutput.simsequences, simoutput.simchoicevec);
        
        % Define full bounds for fitting
        lb              = [-2, 0];
        ub              = [0, 10];
        best_NegLL      = inf;
        best_Xfit       = [];
        options         = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIter', 5000, 'TolFun', 1e-6);


        % Optimisation with random initial conditions
        for start = 1:num_starts
            start_vals      = [unifrnd(-2, 0), unifrnd(0,10)]; % Random starting point within full bounds
            [Xfit, NegLL]   = fmincon(obFunc, start_vals, [], [], [], [], lb, ub, [], options);
            
            if NegLL < best_NegLL
                best_NegLL  = NegLL;
                best_Xfit   = Xfit;
            end
        end

        fit_cs(sample, cond)        = best_Xfit(1);
        fit_beta(sample, cond)      = best_Xfit(2);
        fit_bestNLL(sample, cond)   = best_NegLL;
        fit_NLL(sample, cond)       = NegLL;
        
        % Run MDP model with fitted parameters
        fMDP_output                 = fit_MDPBeads_paramRec(Xfit, simvars, simoutput.simsequences, simoutput.simchoicevec);
        fit_draws(sample, cond)     = fMDP_output.avsamples;
    end
end

%% evaluate parameter recovery

% pearsons correlations
corr_cs_08          = corr(true_cs, fit_cs(:, 1));
corr_beta_08        = corr(true_beta, fit_beta(:, 1));
corr_cs_06          = corr(true_cs, fit_cs(:, 2));
corr_beta_06        = corr(true_beta, fit_beta(:, 2));

% bias
bias_cs_08          = mean(fit_cs(:, 1) - true_cs);
bias_beta_08        = mean(fit_beta(:, 1) - true_beta);
bias_cs_06          = mean(fit_cs(:, 2) - true_cs);
bias_beta_06        = mean(fit_beta(:, 2) - true_beta);

% RMSE
rmse_cs_08          = sqrt(mean((fit_cs(:, 1) - true_cs).^2));
rmse_beta_08        = sqrt(mean((fit_beta(:, 1) - true_beta).^2));
rmse_cs_06          = sqrt(mean((fit_cs(:, 2) - true_cs).^2));
rmse_beta_06        = sqrt(mean((fit_beta(:, 2) - true_beta).^2));

% standard deviation 
std_cs_08       = std(fit_cs(:, 1) - true_cs);
std_beta_08     = std(fit_beta(:, 1) - true_beta);
std_cs_06       = std(fit_cs(:, 2) - true_cs);
std_beta_06     = std(fit_beta(:, 2) - true_beta);

%% display results

fprintf('Parameter Recovery Results:\n');
fprintf('0.8 Condition:\n');
fprintf('Correlation (Cs): %.3f, (Beta): %.3f\n', corr_cs_08, corr_beta_08);
fprintf('Bias (Cs): %.3f, (Beta): %.3f\n', bias_cs_08, bias_beta_08);
fprintf('RMSE (Cs): %.3f, (Beta): %.3f\n', rmse_cs_08, rmse_beta_08);
fprintf('Standard Deviation (Cs): %.3f, (Beta): %.3f\n', std_cs_08, std_beta_08);
fprintf('0.6 Condition:\n');
fprintf('Correlation (Cs): %.3f, (Beta): %.3f\n', corr_cs_06, corr_beta_06);
fprintf('Bias (Cs): %.3f, (Beta): %.3f\n', bias_cs_06, bias_beta_06);
fprintf('RMSE (Cs): %.3f, (Beta): %.3f\n', rmse_cs_06, rmse_beta_06);
fprintf('Standard Deviation (Cs): %.3f, (Beta): %.3f\n', std_cs_06, std_beta_06);

%% visualisation (Scatterplots)

figure;
subplot(2, 2, 1);
scatter(true_cs, fit_cs(:, 1), 'filled');
xlabel('True C_s'); ylabel('Recovered C_s (0.8)'); title('C_s Recovery (0.8)');
line([-2 0], [-2 0], 'Color', 'k', 'LineStyle', '--');
grid on;

subplot(2, 2, 2);
scatter(true_beta, fit_beta(:, 1), 'filled');
xlabel('True \beta'); ylabel('Recovered \beta (0.8)'); title('\beta Recovery (0.8)');
line([0.5 10], [0.5 10], 'Color', 'k', 'LineStyle', '--');
grid on;

subplot(2, 2, 3);
scatter(true_cs, fit_cs(:, 2), 'filled');
xlabel('True C_s'); ylabel('Recovered C_s (0.6)'); title('C_s Recovery (0.6)');
line([-2 0], [-2 0], 'Color', 'k', 'LineStyle', '--');
grid on;

subplot(2, 2, 4);
scatter(true_beta, fit_beta(:, 2), 'filled');
xlabel('True \beta'); ylabel('Recovered \beta (0.6)'); title('\beta Recovery (0.6)');
line([0.5 10], [0.5 10], 'Color', 'k', 'LineStyle', '--');
grid on;

%% evaluate recovery of sampling rates

% correlations
corr_draws_08 = corr(sim_draws(:, 1), fit_draws(:, 1));
corr_draws_06 = corr(sim_draws(:, 2), fit_draws(:, 2));

% bias
bias_draws_08 = mean(fit_draws(:, 1) - sim_draws(:, 1));
bias_draws_06 = mean(fit_draws(:, 2) - sim_draws(:, 2));

% RMSE
rmse_draws_08 = sqrt(mean((fit_draws(:, 1) - sim_draws(:, 1)).^2));
rmse_draws_06 = sqrt(mean((fit_draws(:, 2) - sim_draws(:, 2)).^2));

% standard deviation 
std_draws_08    = std(fit_draws(:, 1) - sim_draws(:, 1));
std_draws_06    = std(fit_draws(:, 2) - sim_draws(:, 2));

%% display results

fprintf('Sampling Rate Recovery Results:\n');
fprintf('0.8 Condition:\n');
fprintf('Correlation: %.3f\n', corr_draws_08);
fprintf('Bias: %.3f\n', bias_draws_08);
fprintf('RMSE: %.3f\n', rmse_draws_08);
fprintf('Standard Deviation: %.3f\n', std_draws_08);
fprintf('0.6 Condition:\n');
fprintf('Correlation: %.3f\n', corr_draws_06);
fprintf('Bias: %.3f\n', bias_draws_06);
fprintf('RMSE: %.3f\n', rmse_draws_06);
fprintf('Standard Deviation: %.3f\n', std_draws_06);

%% Visualisation (Scatterplots)

figure;
% 0.8 Condition
subplot(1, 2, 1);
scatter(sim_draws(:, 1), fit_draws(:, 1), 'filled');
xlabel('True Sampling Rate (0.8)'); 
ylabel('Recovered Sampling Rate (0.8)'); 
title('Sampling Rate Recovery (0.8)');
% Add diagonal line for perfect recovery
max_draws = max([max(sim_draws(:, 1)), max(fit_draws(:, 1))]);
min_draws = min([min(sim_draws(:, 1)), min(fit_draws(:, 1))]);
line([min_draws max_draws], [min_draws max_draws], 'Color', 'k', 'LineStyle', '--');
grid on;

% 0.6 Condition
subplot(1, 2, 2);
scatter(sim_draws(:, 2), fit_draws(:, 2), 'filled');
xlabel('True Sampling Rate (0.6)'); 
ylabel('Recovered Sampling Rate (0.6)'); 
title('Sampling Rate Recovery (0.6)');
% Add diagonal line for perfect recovery
max_draws = max([max(sim_draws(:, 2)), max(fit_draws(:, 2))]);
min_draws = min([min(sim_draws(:, 2)), min(fit_draws(:, 2))]);
line([min_draws max_draws], [min_draws max_draws], 'Color', 'k', 'LineStyle', '--');
grid on;

%% thresholds for "well-recovered"
% parameter recovery thresholds (based on RMSE from 0.8 condition)
threshold_cs_08     = rmse_cs_08;
threshold_beta_08   = rmse_beta_08;
threshold_cs_06     = rmse_cs_06;
threshold_beta_06   = rmse_beta_06;

% sampling rate recovery thresholds (based on RMSE from 0.8 condition)
threshold_draws_08 = rmse_draws_08; % RMSE for sampling rates in 0.8 condition
threshold_draws_06 = rmse_draws_06; % RMSE for sampling rates in 0.8 condition

%% compute absolute errors
% parameter recovery errors
error_cs_08     = abs(fit_cs(:, 1) - true_cs);
error_beta_08   = abs(fit_beta(:, 1) - true_beta);
error_cs_06     = abs(fit_cs(:, 2) - true_cs);
error_beta_06   = abs(fit_beta(:, 2) - true_beta);

% sampling rate recovery errors
error_draws_08  = abs(fit_draws(:, 1) - sim_draws(:, 1));
error_draws_06  = abs(fit_draws(:, 2) - sim_draws(:, 2));

%% visualise absolute error

% Visualize absolute error histograms in two subplots with overlaid conditions
figure;

% Cs Errors
subplot(1, 3, 1);
histogram(error_cs_08, 20, 'FaceColor', [1 0 0], 'FaceAlpha', 0.5, 'DisplayName', '0.8 Condition');
hold on;
histogram(error_cs_06, 20, 'FaceColor', [0 0 1], 'FaceAlpha', 0.5, 'DisplayName', '0.6 Condition');
xlabel('Absolute Error for C_s');
ylabel('Frequency');
title('Histogram of C_s Recovery Errors');
legend('show');
grid on;

% Beta Errors
subplot(1, 3, 2);
histogram(error_beta_08, 20, 'FaceColor', [1 0 0], 'FaceAlpha', 0.5, 'DisplayName', '0.8 Condition');
hold on;
histogram(error_beta_06, 20, 'FaceColor', [0 0 1], 'FaceAlpha', 0.5, 'DisplayName', '0.6 Condition');
xlabel('Absolute Error for \beta');
ylabel('Frequency');
title('Histogram of \beta Recovery Errors');
legend('show');
grid on;

% sapling rates error
subplot(1, 3, 3);
histogram(error_draws_08, 20, 'FaceColor', [1 0 0], 'FaceAlpha', 0.5, 'DisplayName', '0.8 Condition');
hold on;
histogram(error_draws_06, 20, 'FaceColor', [0 0 1], 'FaceAlpha', 0.5, 'DisplayName', '0.6 Condition');
xlabel('Absolute Error for sampling rates');
ylabel('Frequency');
title('Histogram of Samples Recovery Errors');
legend('show');
grid on;

fontsize(gcf, 20, "points");



