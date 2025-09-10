%% Housekeeping Commands
clear all
clc

%% Set Figure-Docking as Default
set(0, 'DefaultFigureWindowStyle', 'docked')

 %% Run Model Comparison for Each Condition

% load model fitting results
load('fitting_final.mat')

% Preallocate BIC and AIC for each model per subject per condition
BIC_model1 = nan(nsubs, conditions);
AIC_model1 = nan(nsubs, conditions);
BIC_model2 = nan(nsubs, conditions);
AIC_model2 = nan(nsubs, conditions);
BIC_model3 = nan(nsubs, conditions);
AIC_model3 = nan(nsubs, conditions);

% Loop over subjects and conditions to calculate BIC and AIC
for sub = 1:nsubs
    for cond = 1:conditions

        % Get total draws for this subject and condition (sample size D)
        % D = sum_cond_draws(sub, cond);  % From your fitting code, sum_cond_draws(sub,cond)
        D = 52;
        logD = log(D);

        % Model 1: free Cs (k = 1)
        NegLL1 = allsub_ll(sub, cond);
        BIC_model1(sub, cond) = 1 * logD + 2 * NegLL1;
        AIC_model1(sub, cond) = 2 * 1 + 2 * NegLL1;

        % Model 2: free Cs + beta (k = 2)
        NegLL2 = allsub_ll_v2(sub, cond);
        BIC_model2(sub, cond) = 2 * logD + 2 * NegLL2;
        AIC_model2(sub, cond) = 2 * 2 + 2 * NegLL2;

        % Model 3: free beta (k = 1)
        NegLL3 = allsub_ll_v3(sub, cond);
        BIC_model3(sub, cond) = 1 * logD + 2 * NegLL3;
        AIC_model3(sub, cond) = 2 * 1 + 2 * NegLL3;
    end
end

% Preallocate best model counts for BIC and AIC per condition
best_model_counts_BIC = zeros(conditions, 3);  % Rows: conditions, Columns: models
best_model_counts_AIC = zeros(conditions, 3);

%%
% Loop over subjects and conditions to determine best model
for sub = 1:nsubs
    for cond = 1:conditions
        
        % BIC comparison
        BIC_values = [BIC_model1(sub, cond), BIC_model2(sub, cond), BIC_model3(sub, cond)];
        [~, best_model_BIC] = min(BIC_values);
        best_model_counts_BIC(cond, best_model_BIC) = best_model_counts_BIC(cond, best_model_BIC) + 1;

        % AIC comparison
        AIC_values = [AIC_model1(sub, cond), AIC_model2(sub, cond), AIC_model3(sub, cond)];
        [~, best_model_AIC] = min(AIC_values);
        best_model_counts_AIC(cond, best_model_AIC) = best_model_counts_AIC(cond, best_model_AIC) + 1;
    end
end

%%
% Normalize to get proportions
prop_best_BIC = best_model_counts_BIC / nsubs;
prop_best_AIC = best_model_counts_AIC / nsubs;

prop_best_AIC = prop_best_AIC'; prop_best_BIC=prop_best_BIC';

% Display results
disp('Proportion of Subjects Best Fit by Each Model (BIC):');
disp('Condition | Model 1 (Cs) | Model 2 (Cs+beta) | Model 3 (beta)');
disp([ (1:conditions)', prop_best_BIC ]);

disp('Proportion of Subjects Best Fit by Each Model (AIC):');
disp('Condition | Model 1 (Cs) | Model 2 (Cs+beta) | Model 3 (beta)');
disp([ (1:conditions)', prop_best_AIC ]);

%%
%% Visualise Model Comparison with Bar Plots (Proportions; Figure 3A)
figure;
subplot(1,2,1);
barHandles = bar(prop_best_BIC'); % Rows: conditions, columns: models
% Set FaceAlpha for each bar group
for i = 1:length(barHandles)
    barHandles(i).FaceAlpha = 0.4; % Set bar transparency to 0.4
end
ylim([0 1])
title('Model Comparison - BIC (Proportions)');
xlabel('Model');
ylabel('Proportion of Subjects');
set(gca, 'XTickLabel', {'CM', 'CNM', 'NM'});
legend({'0.8 condition', '0.6 condition'}, 'Location', 'Best');
grid on;
subplot(1,2,2);
barHandles = bar(prop_best_AIC');
% Set FaceAlpha for each bar group
for i = 1:length(barHandles)
    barHandles(i).FaceAlpha = 0.4; % Set bar transparency to 0.4
end
ylim([0 1])
title('Model Comparison - AIC (Proportions)');
xlabel('Model');
ylabel('Proportion of Subjects');
set(gca, 'XTickLabel', {'CM', 'CNM', 'NM'});
legend({'0.8 condition', '0.6 condition'}, 'Location', 'Best');
grid on;
fontsize(gcf, 20, "points");

%% Group-Level Statistics on BIC Averages or Differences
% Rather than pairwise t-tests, perform repeated-measures ANOVA per condition
% with Model as the factor (three levels), then post-hoc multiple comparisons
% if ANOVA is significant. This handles overall differences and corrects for
% multiple comparisons in post-hoc tests.

% Requires Statistics and Machine Learning Toolbox (fitrm, ranova, multcompare)

for cond = 1:conditions
    
    % Extract BIC data for this condition (nsubs x 3 models)
    BIC_data = [BIC_model1(:,cond), BIC_model2(:,cond), BIC_model3(:,cond)];
    
    % Compute mean BIC per model for reporting
    mean_BIC(cond,:) = mean(BIC_data, 1, 'omitnan');
    fprintf('Mean BIC (Model1,2,3) for Cond %d: %.2f, %.2f, %.2f\n', cond, mean_BIC);
    
    % Create table for repeated-measures ANOVA
    % Columns: Subject ID, Model1 BIC, Model2 BIC, Model3 BIC
    t = table((1:nsubs)', BIC_data(:,1), BIC_data(:,2), BIC_data(:,3), ...
              'VariableNames', {'Subject', 'Model1', 'Model2', 'Model3'});
    
    % Define repeated-measures model: BIC ~ 1 + (1|Subject)
    % Within-subject design: Models as repeated measures
    rm = fitrm(t, 'Model1-Model3 ~ 1', 'WithinDesign', [1 2 3]);
    
    % Run repeated-measures ANOVA
    anova_table = ranova(rm);
    disp(['Repeated-Measures ANOVA for Cond ' num2str(cond) ':']);
    disp(anova_table);
    
    % Check if overall Model effect is significant
    p_anova = anova_table.pValue(1);  % p-value for the Model factor (intercept term is the main effect)
    if p_anova < 0.05
        fprintf('ANOVA significant for Cond %d (p=%.4f). Performing post-hoc comparisons.\n', cond, p_anova);
        
        % Post-hoc multiple comparisons (Tukey's HSD, corrects for multiples)
        mc = multcompare(rm, 'Time', 'ComparisonType', 'bonferroni');  % 'Time' is the default within-subject factor
        disp(['Post-Hoc Comparisons (Bonferroni) for Cond ' num2str(cond) ':']);
        disp(mc);
    else
        fprintf('ANOVA not significant for Cond %d (p=%.4f). No post-hoc needed.\n', cond, p_anova);
    end
end

%% visualise Figure 3B

% Now, create a single grouped boxplot for all models across both conditions
% Prepare data in long format for boxchart with grouping

% Prepare data: Ensure ydata is a vector (stacked BIC values)
% Order: For 0.8: CM, CNM, NM; For 0.6: CM, CNM, NM
BIC_vector = [BIC_model1(:,1); BIC_model2(:,1); BIC_model3(:,1); ...
              BIC_model1(:,2); BIC_model2(:,2); BIC_model3(:,2)];

% Categorical for Conditions (x-groups): Match length of BIC_vector
% Two main groups: 0.8 and 0.6
Condition = categorical([repmat({'0.8'}, nsubs*3, 1); repmat({'0.6'}, nsubs*3, 1)]);

% Categorical for Models (colours): Match length, repeating for each condition
Model = categorical([repmat({'CM'}, nsubs, 1); repmat({'CNM'}, nsubs, 1); repmat({'NM'}, nsubs, 1); ...
                     repmat({'CM'}, nsubs, 1); repmat({'CNM'}, nsubs, 1); repmat({'NM'}, nsubs, 1)]);

% Create grouped boxchart: Group on x by Condition, color by Model
figure;
b = boxchart(Condition, BIC_vector, 'GroupByColor', Model, 'Notch', 'on');

% Customise: Assign colors for models (e.g., blue for CM, red for CNM, green for NM)
% b(1).BoxFaceColor = [0 0 1];  % Blue for CM
% b(2).BoxFaceColor = [1 0 0];  % Red for CNM
% b(3).BoxFaceColor = [0 1 0];  % Green for NM
% b(1).BoxFaceAlpha = 1;        % Filled
% b(2).BoxFaceAlpha = 1;
% b(3).BoxFaceAlpha = 1;

ylabel('BIC');
title('BIC per Condition and Model');
grid on;
legend({'CM', 'CNM', 'NM'}, 'Location', 'northeast');
fontsize(gcf, 20, "points");


%% clear all to run model recovery

clear all 

%% \
% Either run model recovery (takes ~ 5-7 hours to run) or go to line 413 to
% load the model recovery results and visualise the heatmaps and box plots

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

%% Visualize confusion and inverse confusion matrices (Figure 3 D,E)

% load the model recovery results because it takes forever to run (2000
% simulations!!)
load('model_recovery_final.mat')

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

%%