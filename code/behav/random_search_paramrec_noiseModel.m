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
true_cs             = -0.25;

%% Define total samples for random search

total_samples   = 2000;
num_starts      = 5; % based on computational resources

%% Preallocate arrays

true_beta       = zeros(total_samples, 1);
fit_beta        = zeros(total_samples, conditions);
fit_bestNLL     = zeros(total_samples, conditions);
fit_NLL         = zeros(total_samples, conditions);
sim_draws       = zeros(total_samples, conditions);
fit_draws       = zeros(total_samples, conditions);

%% Generate true parameters and perform recovery

for sample = 1:total_samples

    % Sample true parameters from full ranges
    true_beta(sample)   = unifrnd(0, 10);
    
    param               = [true_cs, true_beta(sample)];
    
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
        simvars.cs                  = true_cs;
        
        % Define objective function
        obFunc          = @(x) mdpBeadsNoise_paramRec(x(1), simvars, simoutput.simsequences, simoutput.simchoicevec);
        
        % Define full bounds for fitting
        lb              = 0;
        ub              = 10;
        best_NegLL      = inf;
        best_Xfit       = [];
        options         = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIter', 5000, 'TolFun', 1e-6);


        % Optimisation with random initial conditions
        for start = 1:num_starts
            start_vals      = unifrnd(0, 10); % Random starting point within full bounds
            [Xfit, NegLL]   = fmincon(obFunc, start_vals, [], [], [], [], lb, ub, [], options);
            
            if NegLL < best_NegLL
                best_NegLL  = NegLL;
                best_Xfit   = Xfit;
            end
        end

        fit_beta(sample, cond)      = best_Xfit;
        fit_bestNLL(sample, cond)   = best_NegLL;
        fit_NLL(sample, cond)       = NegLL;
        
        % Run MDP model with fitted parameters
        fMDP_output                 = fit_MDPBeadsNoise_paramRec(Xfit, simvars, simoutput.simsequences, simoutput.simchoicevec);
        fit_draws(sample, cond)     = fMDP_output.avsamples;
    end
end

%% evaluate parameter recovery

% % pearsons correlations
% corr_beta_08        = corr(true_beta, fit_beta(:, 1));
% corr_beta_06        = corr(true_beta, fit_beta(:, 2));
% 
% % bias
% bias_beta_08        = mean(fit_beta(:, 1) - true_beta);
% bias_beta_06        = mean(fit_beta(:, 2) - true_beta);
% 
% % RMSE
% rmse_beta_08        = sqrt(mean((fit_beta(:, 1) - true_beta).^2));
% rmse_beta_06        = sqrt(mean((fit_beta(:, 2) - true_beta).^2));
% 
% % standard deviation 
% std_beta_08     = std(fit_beta(:, 1) - true_beta);
% std_beta_06     = std(fit_beta(:, 2) - true_beta);
% 
% %% display results
% fprintf('Parameter Recovery Results:\n');
% fprintf('0.8 Condition:\n');
% fprintf('Correlation (Beta): %.3f\n', corr_beta_08);
% fprintf('Bias (Beta): %.3f\n', bias_beta_08);
% fprintf('RMSE (Beta): %.3f\n', rmse_beta_08);
% fprintf('Standard Deviation (Beta): %.3f\n', std_beta_08);
% fprintf('0.6 Condition:\n');
% fprintf('Correlation (Beta): %.3f\n', corr_beta_06);
% fprintf('Bias (Beta): %.3f\n', bias_beta_06);
% fprintf('RMSE (Beta): %.3f\n', rmse_beta_06);
% fprintf('Standard Deviation (Beta): %.3f\n', std_beta_06);

%% visualisation (Scatterplots; Supplementary Figure 5A)
figure;

subplot(1, 2, 1);
scatter(true_beta, fit_beta(:, 1), 'filled');
xlabel('True \beta'); ylabel('Recovered \beta (0.8)'); title('\beta Recovery (0.8)');
line([0 10], [0 10], 'Color', 'k', 'LineStyle', '--');
grid on;

subplot(1, 2, 1);
scatter(true_beta, fit_beta(:, 2), 'filled');
xlabel('True \beta'); ylabel('Recovered \beta (0.6)'); title('\beta Recovery (0.6)');
line([0 10], [0 10], 'Color', 'k', 'LineStyle', '--');
grid on;

%% evaluate recovery of sampling rates

% % correlations
% corr_draws_08 = corr(sim_draws(:, 1), fit_draws(:, 1));
% corr_draws_06 = corr(sim_draws(:, 2), fit_draws(:, 2));
% 
% % bias
% bias_draws_08 = mean(fit_draws(:, 1) - sim_draws(:, 1));
% bias_draws_06 = mean(fit_draws(:, 2) - sim_draws(:, 2));
% 
% % RMSE
% rmse_draws_08 = sqrt(mean((fit_draws(:, 1) - sim_draws(:, 1)).^2));
% rmse_draws_06 = sqrt(mean((fit_draws(:, 2) - sim_draws(:, 2)).^2));
% 
% % standard deviation 
% std_draws_08    = std(fit_draws(:, 1) - sim_draws(:, 1));
% std_draws_06    = std(fit_draws(:, 2) - sim_draws(:, 2));
% 
% %% display results
% 
% fprintf('Sampling Rate Recovery Results:\n');
% fprintf('0.8 Condition:\n');
% fprintf('Correlation: %.3f\n', corr_draws_08);
% fprintf('Bias: %.3f\n', bias_draws_08);
% fprintf('RMSE: %.3f\n', rmse_draws_08);
% fprintf('Standard Deviation: %.3f\n', std_draws_08);
% fprintf('0.6 Condition:\n');
% fprintf('Correlation: %.3f\n', corr_draws_06);
% fprintf('Bias: %.3f\n', bias_draws_06);
% fprintf('RMSE: %.3f\n', rmse_draws_06);
% fprintf('Standard Deviation: %.3f\n', std_draws_06);

%% Visualisation (Scatterplots; Supplementary Figure 5B)

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
% % parameter recovery thresholds (based on RMSE from 0.8 condition)
% threshold_beta_08   = rmse_beta_08;
% threshold_beta_06   = rmse_beta_06;
% 
% % sampling rate recovery thresholds (based on RMSE from 0.8 condition)
% threshold_draws_08 = rmse_draws_08; % RMSE for sampling rates in 0.8 condition
% threshold_draws_06 = rmse_draws_06; % RMSE for sampling rates in 0.8 condition
% 
% %% compute absolute errors
% % parameter recovery errors
% error_beta_08 = abs(fit_beta(:, 1) - true_beta);
% error_beta_06 = abs(fit_beta(:, 2) - true_beta);
% 
% % sampling rate recovery errors
% error_draws_08 = abs(fit_draws(:, 1) - sim_draws(:, 1));
% error_draws_06 = abs(fit_draws(:, 2) - sim_draws(:, 2));
% 
% %% visualise absolute error
% 
% % Visualize absolute error histograms in two subplots with overlaid conditions
% figure;
% 
% % Beta Errors
% subplot(1, 2, 1);
% histogram(error_beta_08, 20, 'FaceColor', [1 0 0], 'FaceAlpha', 0.5, 'DisplayName', '0.8 Condition');
% hold on;
% histogram(error_beta_06, 20, 'FaceColor', [0 0 1], 'FaceAlpha', 0.4, 'DisplayName', '0.6 Condition');
% xlabel('Absolute Error for \beta');
% ylabel('Frequency');
% title('Histogram of \beta Recovery Errors');
% legend('show');
% grid on;
% 
% % sapling rates error
% subplot(1, 2, 2);
% histogram(error_draws_08, 20, 'FaceColor', [1 0 0], 'FaceAlpha', 0.5, 'DisplayName', '0.8 Condition');
% hold on;
% histogram(error_draws_06, 20, 'FaceColor', [0 0 1], 'FaceAlpha', 0.4, 'DisplayName', '0.6 Condition');
% xlabel('Absolute Error for sampling rates');
% ylabel('Frequency');
% title('Histogram of Samples Recovery Errors');
% legend('show');
% grid on;
% 
% fontsize(gcf, 20, "points");



%%

