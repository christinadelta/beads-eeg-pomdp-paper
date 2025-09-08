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

% Define parameters for the ideal observer
R.error             = -10;
R.correct           = 10;
R.diff              = -20;
R.q                 = [0.8 0.6];

% Model recovery parameters
num_simulations     = 2000; % Number of simulations per model per condition
num_starts          = 5; % Number of starting points for optimization
nmodels             = 3; % Number of models

% Preallocate confusion and inverse confusion matrices for BIC
confusion_08_BIC        = zeros(nmodels, nmodels); % For 0.8 condition, BIC
confusion_06_BIC        = zeros(nmodels, nmodels); % For 0.6 condition, BIC
inv_confusion_08_BIC    = zeros(nmodels, nmodels); % Inverse for 0.8 condition, BIC
inv_confusion_06_BIC    = zeros(nmodels, nmodels); % Inverse for 0.6 condition, BIC

% Preallocate for average BIC values (true_model x fitted_model x simulations)
BIC_sims_08 = zeros(nmodels, nmodels, num_simulations); % For 0.8 condition
BIC_sims_06 = zeros(nmodels, nmodels, num_simulations); % For 0.6 condition

%% Model Recovery Loop

for cond = 1:conditions

    simvars.cond    = cond;
    simvars.thisq   = simvars.qvals(cond);
    R.thisq         = R.q(cond);

    for true_model = 1:nmodels

        for sim = 1:num_simulations

            % Simulate data from the true model
            if true_model == 1 % Model 1: Cs only
                true_Cs     = unifrnd(-2, 0);
                true_beta   = 2; % Fixed
                param       = [true_Cs];
                simoutput   = sim_POMDP_Beads(simvars, [true_Cs, true_beta]);

            elseif true_model == 2 % Model 2: Cs + beta 
               
                true_Cs     = unifrnd(-2, 0);
                true_beta   = unifrnd(0, 10);
                param       = [true_Cs, true_beta];
                simoutput   = sim_POMDP_Beads(simvars, param);

            elseif true_model == 3 % Model 3: beta only

                true_Cs     = -0.25; % Fixed
                true_beta   = unifrnd(0, 10);
                param       = [true_beta];
                simoutput   = sim_POMDP_Beads(simvars, [true_Cs, true_beta]);
                
            end

            % Store simulation results
            sim_sequences   = simoutput.simsequences;
            sim_choices     = simoutput.simchoicevec;
            sim_urns        = simoutput.simurns;
            R.urntype       = sim_urns;
            R.beta          = 2;
            D               = sum(simoutput.simdraws); % Total draws for BIC/AIC
            logD            = log(D);

            % Fit all three models
            options         = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIter', 5000, 'TolFun', 1e-6, 'TolX', 1e-6);

            % Fit Model 1 (Cs only)
            obFunc1         = @(x) mdpBeadsCost_paramRec([x(1)], R, sim_sequences, sim_choices);
            best_NegLL1     = inf;
            best_Xfit1      = [];

            for start = 1:num_starts
                start_vals      = unifrnd(-2, 0);
                [Xfit, NegLL]   = fmincon(obFunc1, start_vals, [], [], [], [], -2, 0, [], options);
                if NegLL < best_NegLL1
                    best_NegLL1 = NegLL;
                    best_Xfit1  = Xfit;
                end
            end

            BIC1 = 1 * logD + 2 * best_NegLL1; % k = 1 for Model 1

            % Fit Model 2 (Cs + beta)
            clear R.cs start_vals
            
            obFunc2             = @(x) mdpBeads_paramRec([x(1), x(2)], R, sim_sequences, sim_choices);
            best_NegLL2         = inf;
            best_Xfit2          = [];

            for start = 1:num_starts
                start_vals      = [unifrnd(-2, 0), unifrnd(0, 10)];
                [Xfit, NegLL]   = fmincon(obFunc2, start_vals, [], [], [], [], [-2, 0], [0, 10], [], options);
                if NegLL < best_NegLL2
                    best_NegLL2 = NegLL;
                    best_Xfit2  = Xfit;
                end
            end

            BIC2 = 2 * logD + 2 * best_NegLL2; % k = 2 for Model 2

            % Fit Model 3 (beta only)
            clear R.beta start_vals
            
            R.cs                = -0.25;
            obFunc3             = @(x) mdpBeadsNoise_paramRec([x(1)], R, sim_sequences, sim_choices);
            best_NegLL3         = inf;
            best_Xfit3          = [];

            for start = 1:num_starts
                
                start_vals      = unifrnd(0, 10);
                [Xfit, NegLL]   = fmincon(obFunc3, start_vals, [], [], [], [], 0, 10, [], options);
                
                if NegLL < best_NegLL3
                    best_NegLL3 = NegLL;
                    best_Xfit3  = Xfit;
                end
            end

            BIC3 = 1 * logD + 2 * best_NegLL3; % k = 1 for Model 3

            % Find the model with the smallest BIC
            BIC_values          = [BIC1, BIC2, BIC3];
            [~, best_model_BIC] = min(BIC_values);

            % Update confusion matrices
            if cond == 1
                confusion_08_BIC(true_model, best_model_BIC) = confusion_08_BIC(true_model, best_model_BIC) + 1;
                BIC_sims_08(true_model, 1, sim) = BIC1;
                BIC_sims_08(true_model, 2, sim) = BIC2;
                BIC_sims_08(true_model, 3, sim) = BIC3;
            else
                confusion_06_BIC(true_model, best_model_BIC) = confusion_06_BIC(true_model, best_model_BIC) + 1;
                BIC_sims_06(true_model, 1, sim) = BIC1;
                BIC_sims_06(true_model, 2, sim) = BIC2;
                BIC_sims_06(true_model, 3, sim) = BIC3;
            end
           
        end % end of simulations 
    end % end of model loop
end % end of condition loop

%% Compute Inverse Confusion Matrices

% Normalize confusion matrices
confusion_08_BIC = confusion_08_BIC / num_simulations;
confusion_06_BIC = confusion_06_BIC / num_simulations;

% Compute inverse confusion matrices using Bayes' rule
% P(simulated model X | fit model Y) = P(fit model Y | simulated model X) * P(simulated model X) / P(fit model Y)
% Assuming equal priors: P(simulated model X) = 1/nmodels
prior       = 1/nmodels;

for fit_model = 1:nmodels
    % For 0.8 condition, BIC
    p_fit_model_08_BIC = sum(confusion_08_BIC(:, fit_model)) / nmodels; % P(fit model Y)
    if p_fit_model_08_BIC > 0
        inv_confusion_08_BIC(:, fit_model) = (confusion_08_BIC(:, fit_model) * prior) / p_fit_model_08_BIC;
    end
    % For 0.6 condition, BIC
    p_fit_model_06_BIC = sum(confusion_06_BIC(:, fit_model)) / nmodels; % P(fit model Y)
    if p_fit_model_06_BIC > 0
        inv_confusion_06_BIC(:, fit_model) = (confusion_06_BIC(:, fit_model) * prior) / p_fit_model_06_BIC;
    end
end

%% Compute Average BIC Matrices
mean_BIC_08 = mean(BIC_sims_08, 3, 'omitnan'); % Average over simulations for 0.8
mean_BIC_06 = mean(BIC_sims_06, 3, 'omitnan'); % Average over simulations for 0.6

%% Plot the Model Recovery Results as Heatmaps

% Display confusion and inverse confusion matrices
disp('Confusion Matrix for 0.8 Condition (BIC):');
disp(['  Cost Model  Cost-Noise Model Noise Model']);
disp(confusion_08_BIC);
disp('Inverse Confusion Matrix for 0.8 Condition (BIC):');
disp(['  Cost Model  Cost-Noise Model Noise Model']);
disp(inv_confusion_08_BIC);
disp('Confusion Matrix for 0.6 Condition (BIC):');
disp(['  Cost Model  Cost-Noise Model Noise Model ']);
disp(confusion_06_BIC);
disp('Inverse Confusion Matrix for 0.6 Condition (BIC):');
disp(['  Cost Model  Cost-Noise Model Noise Model']);
disp(inv_confusion_06_BIC);

% Display average BIC matrices
disp('Average BIC Matrix for 0.8 Condition:');
disp(['  Cost Model  Cost-Noise Model Noise Model']);
disp(mean_BIC_08);
disp('Average BIC Matrix for 0.6 Condition:');
disp(['  Cost Model  Cost-Noise Model Noise Model']);
disp(mean_BIC_06);

%% Visualize confusion and inverse confusion matrices

figure;

% BIC Confusion Matrices
subplot(2, 2, 1);
imagesc(confusion_08_BIC);
colorbar;
title('Confusion Matrix - 0.8 Condition (BIC)');
xlabel('Best-Fitting Model');
ylabel('True Model');
set(gca, 'XTick', 1:3, 'XTickLabel', {'CM', 'CNM', 'NM'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'CM', 'CNM', 'NM'});
for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.2f', confusion_08_BIC(i,j)), 'Color', 'white', 'HorizontalAlignment', 'center');
    end
end

subplot(2, 2, 2);
imagesc(confusion_06_BIC);
colorbar;
title('Confusion Matrix - 0.6 Condition (BIC)');
xlabel('Best-Fitting Model');
ylabel('True Model');
set(gca, 'XTick', 1:3, 'XTickLabel', {'CM', 'CNM', 'NM'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'CM', 'CNM', 'NM'});
for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.2f', confusion_06_BIC(i,j)), 'Color', 'white', 'HorizontalAlignment', 'center');
    end
end

% BIC Inverse Confusion Matrices
subplot(2, 2, 3);
imagesc(inv_confusion_08_BIC);
colorbar;
title('Inverse Confusion Matrix - 0.8 Condition (BIC)');
xlabel('Best-Fitting Model');
ylabel('True Model');
set(gca, 'XTick', 1:3, 'XTickLabel', {'CM', 'CNM', 'NM'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'CM', 'CNM', 'NM'});
for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.2f', inv_confusion_08_BIC(i,j)), 'Color', 'white', 'HorizontalAlignment', 'center');
    end
end

subplot(2, 2, 4);
imagesc(inv_confusion_06_BIC);
colorbar;
title('Inverse Confusion Matrix - 0.6 Condition (BIC)');
xlabel('Best-Fitting Model');
ylabel('True Model');
set(gca, 'XTick', 1:3, 'XTickLabel', {'CM', 'CNM', 'NM'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'CM', 'CNM', 'NM'});
for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.2f', inv_confusion_06_BIC(i,j)), 'Color', 'white', 'HorizontalAlignment', 'center');
    end
end

fontsize(gcf, 20, "points");

%% Visualize average BIC matrices as heatmaps

figure;

% Average BIC Matrices
subplot(1, 2, 1);
imagesc(mean_BIC_08);
colorbar;
title('Average BIC - 0.8 Condition');
xlabel('Fitted Model');
ylabel('True Model');
set(gca, 'XTick', 1:3, 'XTickLabel', {'CM', 'CNM', 'NM'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'CM', 'CNM', 'NM'});
for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.2f', mean_BIC_08(i,j)), 'Color', 'white', 'HorizontalAlignment', 'center');
    end
end

subplot(1, 2, 2);
imagesc(mean_BIC_06);
colorbar;
title('Average BIC - 0.6 Condition');
xlabel('Fitted Model');
ylabel('True Model');
set(gca, 'XTick', 1:3, 'XTickLabel', {'CM', 'CNM', 'NM'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'CM', 'CNM', 'NM'});
for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.2f', mean_BIC_06(i,j)), 'Color', 'white', 'HorizontalAlignment', 'center');
    end
end

fontsize(gcf, 20, "points");

%% Visualize BIC Distributions as Side-by-Side Boxplots (9 per condition)

% Labels for true-fitted pairs
pair_labels = {'CM-CM', 'CM-CNM', 'CM-NM', 'CNM-CM', 'CNM-CNM', 'CNM-NM', 'NM-CM', 'NM-CNM', 'NM-NM'};

% For 0.8 Condition: Prepare data as matrix (num_simulations rows x 9 columns)
BIC_08_matrix = zeros(num_simulations, 9);
col_idx = 1;
for true_idx = 1:3
    for fit_idx = 1:3
        BIC_08_matrix(:, col_idx) = squeeze(BIC_sims_08(true_idx, fit_idx, :));
        col_idx = col_idx + 1;
    end
end

% Plot for 0.8 Condition using boxplot
figure('Name', 'BIC Distributions - 0.8 Condition');
boxplot(BIC_08_matrix, 'Notch', 'off', 'Labels', pair_labels, ...
        'LabelOrientation', 'horizontal', 'Widths', 0.5);
ylabel('BIC');
title('BIC Distributions per True-Fitted Pair (0.8 Condition)');
grid on;
set(gca, 'XTickLabelRotation', 45);  % Rotate labels for readability
fontsize(gcf, 20, "points");

% Customise colours (blue for CM fitted, red for CNM, green for NM)
h = findobj(gca, 'Tag', 'Box');
colors = [0 0 1; 1 0 0; 0 1 0; 0 0 1; 1 0 0; 0 1 0; 0 0 1; 1 0 0; 0 1 0];  % Repeating per fitted model
for j = 1:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(j, :), 'FaceAlpha', 0.5);
end
legend({'Fitted: CM', 'Fitted: CNM', 'Fitted: NM'}, 'Location', 'northeast');

% For 0.6 Condition: Prepare data similarly
BIC_06_matrix = zeros(num_simulations, 9);
col_idx = 1;
for true_idx = 1:3
    for fit_idx = 1:3
        BIC_06_matrix(:, col_idx) = squeeze(BIC_sims_06(true_idx, fit_idx, :));
        col_idx = col_idx + 1;
    end
end

% Plot for 0.6 Condition using boxplot
figure('Name', 'BIC Distributions - 0.6 Condition');
boxplot(BIC_06_matrix, 'Notch', 'off', 'Labels', pair_labels, ...
        'LabelOrientation', 'horizontal', 'Widths', 0.5);
ylabel('BIC');
title('BIC Distributions per True-Fitted Pair (0.6 Condition)');
grid on;
set(gca, 'XTickLabelRotation', 45);
fontsize(gcf, 20, "points");

% Customise colours (same as above)
h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(j, :), 'FaceAlpha', 0.5);
end
legend({'Fitted: CM', 'Fitted: CNM', 'Fitted: NM'}, 'Location', 'northeast');

%%


