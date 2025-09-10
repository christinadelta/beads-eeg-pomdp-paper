% fit MDP beads model (see Furl and Averbeck, 2011; Averbeck, 2015) using
% the fmincon tool

% Model 1: fits a simple Cs model while all other parameters are fixed 
% Model 2: fits a Cs + Beta model 
% Model 3: fits a simple beta model while all other parameters are fixed 


% I will try inluding other parameters as I know that 

%% housekeeping commands

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% INIT LOAD DATA %%
% GET PATHS & DEFINE VARIABLES
% The four next lines (paths) should be changed to your paths 
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
datapath        = '/Users/christinadelta/gitrepos/SocialDeMa/experiments/results';
resultspath     = fullfile(datapath, 'beads');
croppedpath     = fullfile(startpath, 'analysis', 'beads', 'behav', 'cropped');

task            = 'beads';
subs            = dir(fullfile(resultspath, '*sub*'));
nsubs           = length(subs);
%nsubs           = 40;

totaltrials     = 52; 
conditions      = 2;
condtrials      = totaltrials/conditions;
nmodels         = 3;

% init required variables
avdraws                 = nan(nsubs,1);
avacc                   = nan(nsubs,1);
easy_avdraws            = nan(nsubs,1);
diff_avdraws            = nan(nsubs,1);
easy_avacc              = nan(nsubs,1);
diff_avacc              = nan(nsubs,1);

% Define number of starts
num_starts              = 1;

%% fit cost model to participant data

for sub = 1:nsubs

    fprintf('loading beads block data\n') 
    subject = subs(sub).name;
    subdir  = fullfile(resultspath,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % extract blocktrial data
    [subsequences,subchoiceVec,all_data, draws_index]    = extract_blockdata(subdir,sub,task);
    
    % store in cell for each participant
    allsub_alldata{1,sub}                  = all_data;
    allsub_sequences{1,sub}                = subsequences;
    allsub_choiceVec{1,sub}                = subchoiceVec;
    allsub_drawinfo{1,sub}                 = draws_index; % col1: draw, col2: trial

    % SPLIT SUB DATA INTO CONDITIONS
    
    for cond = 1:conditions % loop over conditions
        
        tmp                                 = find(all_data(:,8) == cond);
        cond_data{1,sub}{cond}              = all_data((tmp),:);
        clear tmp
        
    end % end of condition loop

    % AVERAGE PARTICIPANT DRAWS & ACCURACY & CALCULATE POINTS 
    % create a nx1 vector (n=number of participants) with the averaged number
    % of draws for each participant.
    % This vector will be used as a covariate for the individual differences
    % analysis in SPM12.
    
    sub_draws           = all_data(:,5);
    sub_acc             = all_data(:,7);
    avdraws(sub,1)      = mean(sub_draws-1); % -1 because the last draw is the urn choice
    avacc(sub,1)        = mean(sub_acc);
    allsub_draws(:,sub) = all_data(:,5);
    allsub_conds(:,sub) = all_data(:,8);
    
    clear sub_draws sub_acc
    
    % average draws and acc for each condition
    for cond = 1:conditions
        
        tmp_cond                = cond_data{1,sub}{1,cond};
        cond_draws              = tmp_cond(:,5);
        cond_acc                = tmp_cond(:,7);
        
        if cond == 1
            easy_avdraws(sub,1) = mean(cond_draws-1);
            easy_avacc(sub,1)   = mean(cond_acc);
        else
            diff_avdraws(sub,1) = mean(cond_draws-1);
            diff_avacc(sub,1)   = mean(cond_acc);
            
        end

        % sum draws for testing 
        sum_cond_draws(sub,cond) = sum(cond_draws);

        draws_for_testing{cond} = cond_draws;

        
        % we will now calculate participant points. we can use this to se how each participant performed and to compare
        % each participant with their corresponding model instant 
    
        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allsub_points(sub, cond) = (sum(cond_acc==1)*10) + (sum(cond_acc==0)*-10) + (sum(cond_draws)*-0.25);

    end 
    
    % sum participant draws for all conditions
    sum_draws(sub,1) = sum(sum_cond_draws(sub,:));

    
    % fit MDP beads model with only one free parameter (cost-sample)
    
    % define parameters of the ideal observer
    R.beta              = 2;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
    R.error             = -10;          % cost for being wrong
    R.correct           = 10;           % reward for being correct
    R.diff              = -20; 
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    % R.sample            = -0.25;        % the cost to sample

    
    for cond = 1:conditions

        init_sample     = -0.25; % initial value for sampling 
        cond_urns       = cond_data{1,sub}{1,cond}(:,4);
        cond_sequences  = subsequences{1,cond};
        cond_choiceData = subchoiceVec{1,cond};
        R.urns          = cond_urns';
        R.thisq         = R.q(cond);
        
        % testing the number of draws
        R.subdraws          = draws_for_testing{cond};

        % define the objective function and get optimal parameters
        obFunc                          = @(x) mdpBeads([x(1)], R, cond_sequences, cond_choiceData);
        lb                              = -2;
        ub                              = 0;
        options                         = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIter', 5000, 'TolFun', 1e-6);

        % Initialize best fit variables
        best_NegLL  = inf;
        best_Xfit   = [];
    
        % Run optimization with multiple starting points
        for start = 1:num_starts
            % Generate random starting points within bounds
            % Note: Adjust bounds if necessary (e.g., ub(1) should be 0 if Cs <= 0)
            starting_points = [unifrnd(lb(1), ub(1))];
    
            % Call fmincon
            [Xfit, NegLL]   = fmincon(obFunc, starting_points, [], [], [], [], lb, ub, [], options);
    
            % Update best fit if current is better
            if NegLL < best_NegLL
                best_NegLL  = NegLL;
                best_Xfit   = Xfit;
            end
        end
        allsub_fitted_sample(sub, cond) = best_Xfit;
        allsub_ll(sub, cond)            = best_NegLL;

        % quickly test the obejctive function outside the optimisation
        % toolbox
        % ll = mdpBeads(init_sample,R,cond_sequences,cond_choiceData);

        % run the MDP model with the Xfit parameter values and compute
        % % model sampling rate
        fMDP_output = fit_subMDPBeads_v1(best_Xfit,R,cond_sequences,cond_choiceData);
        % 
        % % extract output
        fit_correct_v1(sub,cond)        = fMDP_output.performance;
        fit_samples_v1(sub,cond)        = fMDP_output.avsamples;
        fit_AQs_v1{sub,cond}            = fMDP_output.actionVals;
        allsub_fit_model_v1{sub,cond}   = fMDP_output;

        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
       % allCsModel_points(sub, cond) = (sum(fMDP_output.choice==1)*10) + (sum(fMDP_output.choice==0)*-10) + (sum(fMDP_output.samples)*Xfit);

        clear Xfit NegLL cond_sequences cond_choiceData R.urns R.thisq cond_urns fMDP_output

    end % end of conditions loop
end % end of subjects loop

%% plot parameter values and nll (Supplementary figure S1 B)

% Define Data for Two Conditions
Cs_0_8          = allsub_fitted_sample(:,1); 
Cs_0_6          = allsub_fitted_sample(:,2);  

% Calculate Means and SEM
mean_Cs = [mean(Cs_0_8), mean(Cs_0_6)];
sem_Cs = [std(Cs_0_8) / sqrt(length(Cs_0_8)), std(Cs_0_6) / sqrt(length(Cs_0_6))];

% Create Bar Plot
figure;
bar(1:2, mean_Cs, 0.5, 'FaceColor', [0.3 0.6 0.9]); % Solid light blue bars
hold on;

% Add Error Bars
errorbar(1:2, mean_Cs, sem_Cs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Customize Plot
xticks([1 2]);
xticklabels({'0.8 Condition', '0.6 Condition'});
ylabel('Cs Parameter');
xlabel('Condition');
title('Differences in Cs Parameter Between Conditions');
ylim([-2.2, 0.5]); % Adjust Y-axis range as needed

% Add Individual Data Points with Jitter
jitter = 0.05; % Amount of random noise
scatter(1 + (rand(size(Cs_0_8)) - 0.5) * jitter, Cs_0_8, 'filled', 'r'); % Jittered 0.8 condition
scatter(2 + (rand(size(Cs_0_6)) - 0.5) * jitter, Cs_0_6, 'filled', 'r'); % Jittered 0.6 condition

% Add Legend
% legend({'Mean ± SEM', '0.8 Condition Data', '0.6 Condition Data'}, 'Location', 'Best');
fontsize(gcf, 20, "points");
hold off;

% plot nll 

% Create Density Plot
figure;
[f1, xi1] = ksdensity(allsub_ll(:, 1)); % 0.8 Condition
[f2, xi2] = ksdensity(allsub_ll(:, 2)); % 0.6 Condition

% Plot filled density curves
fill(xi1, f1, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2); % Filled density for 0.8
hold on;
fill(xi2, f2, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'LineWidth', 2); % Filled density for 0.6

% Customise Plot
legend({'0.8 Condition', '0.6 Condition'}, 'Location', 'Best');
xlabel('Negative Log-Likelihood (NLL)');
ylabel('Density');
title('Density Plot of NLL Values');
fontsize(gcf, 20, "points");
hold off;


%% fit cost-noise model to participant data

% loop over subjects and extarct behavioural data 

for sub = 1:nsubs

    fprintf('loading beads block data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(resultspath,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % extract blocktrial data
    [subsequences,subchoiceVec,all_data, draws_index]    = extract_blockdata(subdir,sub,task);
    
    % store in cell for each participant
    allsub_alldata{1,sub}                  = all_data;
    allsub_sequences{1,sub}                = subsequences;
    allsub_choiceVec{1,sub}                = subchoiceVec;
    allsub_drawinfo{1,sub}                 = draws_index; % col1: draw, col2: trial

    % SPLIT SUB DATA INTO CONDITIONS
    
    for cond = 1:conditions % loop over conditions
        
        tmp                                 = find(all_data(:,8) == cond);
        cond_data{1,sub}{cond}                = all_data((tmp),:);
        clear tmp
        
    end % end of condition loop

    % AVERAGE PARTICIPANT DRAWS & ACCURACY & CALCULATE POINTS 
    % create a nx1 vector (n=number of participants) with the averaged number
    % of draws for each participant.
    % This vector will be used as a covariate for the individual differences
    % analysis in SPM12.
    
    sub_draws           = all_data(:,5);
    sub_acc             = all_data(:,7);
    avdraws(sub,1)      = mean(sub_draws-1); % -1 because the last draw is the urn choice
    avacc(sub,1)        = mean(sub_acc);
    allsub_draws(:,sub) = all_data(:,5);
    allsub_conds(:,sub) = all_data(:,8);
    
    clear sub_draws sub_acc
    
    % average draws and acc for each condition
    for cond = 1:conditions
        
        tmp_cond                = cond_data{1,sub}{1,cond};
        cond_draws              = tmp_cond(:,5);
        cond_acc                = tmp_cond(:,7);
        
        if cond == 1
            easy_avdraws(sub,1) = mean(cond_draws-1);
            easy_avacc(sub,1)   = mean(cond_acc);
        else
            diff_avdraws(sub,1) = mean(cond_draws-1);
            diff_avacc(sub,1)   = mean(cond_acc);
            
        end

        % sum draws for testing 
        sum_cond_draws(sub,cond) = sum(cond_draws);

        
        % we will now calculate participant points. we can use this to se how each participant performed and to compare
        % each participant with their corresponding model instant 
    
        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allsub_points(sub, cond) = (sum(cond_acc==1)*10) + (sum(cond_acc==0)*-10) + (sum(cond_draws)*-0.25);

    end 
    
    % sum participant draws for all conditions
    sum_draws(sub,1) = sum(sum_cond_draws(sub,:));

    % fit MDP beads model with two free parameters (cost-sample)
    
    % define parameters of the ideal observer
    R.error             = -10;          % cost for being wrong
    R.correct           = 10;           % reward for being correct
    R.diff              = -20; 
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    % R.sample            = -0.25;        % the cost to sample

    for cond = 1:conditions

        init_beta       = 2; % initial value for inverse temperature
        init_sample     = -0.25; % initial value for sampling 
        % starting_points = [init_sample init_beta];
        cond_urns       = cond_data{1,sub}{1,cond}(:,4);
        cond_sequences  = subsequences{1,cond};
        cond_choiceData = subchoiceVec{1,cond};
        R.urns          = cond_urns';
        R.thisq         = R.q(cond);

        % % define the objective function and get optimal parameters
        lb = [-2, 0]; % Lower bounds: cost-sample, beta
        ub = [0, 10];   % Upper bounds: cost-sample, beta
        
        % Set optimisation options for fmincon
        options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIter', 5000, 'TolFun', 1e-6);
        
        % Define the objective function
        obFunc = @(x) mdpBeads_model2([x(1), x(2)], R, cond_sequences, cond_choiceData);
        
        % Initialize best fit variables
        best_NegLL = inf;
        best_Xfit   = [];
    
        % Run optimization with multiple starting points
        for start = 1:num_starts
            % Generate random starting points within bounds
            % Note: Adjust bounds if necessary (e.g., ub(1) should be 0 if Cs <= 0)
            starting_points = [unifrnd(lb(1), ub(1)), unifrnd(lb(2), ub(2))];
    
            % Call fmincon
            [Xfit, NegLL] = fmincon(obFunc, starting_points, [], [], [], [], lb, ub, [], options);
    
            % Update best fit if current is better
            if NegLL < best_NegLL
                best_NegLL  = NegLL;
                best_Xfit   = Xfit;
            end
        end
        
        % Store results
        allsub_fitted_sample_v2(sub,cond) = best_Xfit(1);
        allsub_fitted_beta_v2(sub,cond) = best_Xfit(2);
        allsub_ll_v2(sub,cond) = best_NegLL;

        % ll = mdpBeads_model2(starting_points,R,cond_sequences,cond_choiceData);

        % run the MDP model with the Xfit parameter values and compute
        % model sampling rate
        fMDP_output = fit_subMDPBeads_v2(best_Xfit,R,cond_sequences,cond_choiceData);

        % extract output
        fit_correct_v2(sub,cond)        = fMDP_output.performance;
        fit_samples_v2(sub,cond)        = fMDP_output.avsamples;
        fit_AQs_v2{sub,cond}            = fMDP_output.actionVals;
        allsub_fit_model_v2{sub,cond}   = fMDP_output;

        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allCsBetaModel_points(sub, cond) = (sum(fMDP_output.choice==1)*10) + (sum(fMDP_output.choice==0)*-10) + (sum(fMDP_output.samples)*Xfit(1));


        clear best_Xfit best_NegLL cond_sequences cond_choiceData R.urns R.thisq cond_urns fMDP_output starting_points
    end % end of conditions loop

end % end of subjects loop

%% plot parameter values and nll (Figure 4 B and C)

% Define Data for Two Conditions
Cs_0_8_model2          = allsub_fitted_sample_v2(:,1);
Cs_0_6_model2          = allsub_fitted_sample_v2(:,2);

beta_0_8_model2          = allsub_fitted_beta_v2(:,1);
beta_0_6_model2          = allsub_fitted_beta_v2(:,2);

% Calculate Means and SEM
mean_Cs = [mean(Cs_0_8_model2), mean(Cs_0_6_model2)];
sem_Cs = [std(Cs_0_8_model2) / sqrt(length(Cs_0_8_model2)), std(Cs_0_6_model2) / sqrt(length(Cs_0_6_model2))];

% Create Bar Plot
figure;
bar(1:2, mean_Cs, 0.5, 'FaceColor', [0.3 0.6 0.9]); % Solid light blue bars
hold on;

% Add Error Bars
errorbar(1:2, mean_Cs, sem_Cs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Customize Plot
xticks([1 2]);
xticklabels({'0.8 Condition', '0.6 Condition'});
ylabel('Cs Parameter');
xlabel('Condition');
title('Differences in Cs Parameter Between Conditions');
ylim([-2.2, 0.5]); % Adjust Y-axis range as needed

% Add Individual Data Points with Jitter
jitter = 0.05; % Amount of random noise
scatter(1 + (rand(size(Cs_0_8_model2)) - 0.5) * jitter, Cs_0_8_model2, 'filled', 'r'); % Jittered 0.8 condition
scatter(2 + (rand(size(Cs_0_6_model2)) - 0.5) * jitter, Cs_0_6_model2, 'filled', 'r'); % Jittered 0.6 condition

% Add Legend
% legend({'Mean ± SEM', '0.8 Condition Data', '0.6 Condition Data'}, 'Location', 'Best');
fontsize(gcf, 20, "points");
hold off;

%%% plot beta values
% Calculate Means and SEM
mean_beta = [mean(beta_0_8_model2), mean(beta_0_6_model2)];
sem_beta = [std(beta_0_8_model2) / sqrt(length(beta_0_8_model2)), std(beta_0_6_model2) / sqrt(length(beta_0_6_model2))];

% Create Bar Plot
figure;
bar(1:2, mean_beta, 0.5, 'FaceColor', [0.3 0.6 0.9]); % Solid light blue bars
hold on;

% Add Error Bars
errorbar(1:2, mean_beta, sem_beta, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Customize Plot
xticks([1 2]);
xticklabels({'0.8 Condition', '0.6 Condition'});
ylabel('Beta Parameter');
xlabel('Condition');
title('Differences in Beta Parameter Between Conditions');
ylim([0 12]); % Adjust Y-axis range as needed

% Add Individual Data Points with Jitter
jitter = 0.05; % Amount of random noise
scatter(1 + (rand(size(beta_0_8_model2)) - 0.5) * jitter, beta_0_8_model2, 'filled', 'r'); % Jittered 0.8 condition
scatter(2 + (rand(size(beta_0_6_model2)) - 0.5) * jitter, beta_0_6_model2, 'filled', 'r'); % Jittered 0.6 condition

% Add Legend
% legend({'Mean ± SEM', '0.8 Condition Data', '0.6 Condition Data'}, 'Location', 'Best');
fontsize(gcf, 20, "points");
hold off;

% plot nll 

% Create Density Plot
figure;
[f1, xi1] = ksdensity(allsub_ll_v2(:, 1)); % 0.8 Condition
[f2, xi2] = ksdensity(allsub_ll_v2(:, 2)); % 0.6 Condition

% Plot filled density curves
fill(xi1, f1, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2); % Filled density for 0.8
hold on;
fill(xi2, f2, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'LineWidth', 2); % Filled density for 0.6

% Customise Plot
legend({'0.8 Condition', '0.6 Condition'}, 'Location', 'Best');
xlabel('Negative Log-Likelihood (NLL)');
ylabel('Density');
title('Density Plot of NLL Values');
fontsize(gcf, 20, "points");
hold off;


%% fit noise model to participant data 

% loop over subjects and extarct behavioural data 
for sub = 1:nsubs

    fprintf('loading beads block data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(resultspath,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % extract blocktrial data
    [subsequences,subchoiceVec,all_data, draws_index]    = extract_blockdata(subdir,sub,task);
    
    % store in cell for each participant
    allsub_alldata{1,sub}                  = all_data;
    allsub_sequences{1,sub}                = subsequences;
    allsub_choiceVec{1,sub}                = subchoiceVec;
    allsub_drawinfo{1,sub}                 = draws_index; % col1: draw, col2: trial

    % SPLIT SUB DATA INTO CONDITIONS
    
    for cond = 1:conditions % loop over conditions
        
        tmp                                 = find(all_data(:,8) == cond);
        cond_data{1,sub}{cond}                = all_data((tmp),:);
        clear tmp
        
    end % end of condition loop

    % AVERAGE PARTICIPANT DRAWS & ACCURACY & CALCULATE POINTS 
    % create a nx1 vector (n=number of participants) with the averaged number
    % of draws for each participant.
    % This vector will be used as a covariate for the individual differences
    % analysis in SPM12.
    
    sub_draws           = all_data(:,5);
    sub_acc             = all_data(:,7);
    avdraws(sub,1)      = mean(sub_draws-1); % -1 because the last draw is the urn choice
    avacc(sub,1)        = mean(sub_acc);
    allsub_draws(:,sub) = all_data(:,5);
    allsub_conds(:,sub) = all_data(:,8);
    
    clear sub_draws sub_acc
    
    % average draws and acc for each condition
    for cond = 1:conditions
        
        tmp_cond                = cond_data{1,sub}{1,cond};
        cond_draws              = tmp_cond(:,5);
        cond_acc                = tmp_cond(:,7);
        
        if cond == 1
            easy_avdraws(sub,1) = mean(cond_draws-1);
            easy_avacc(sub,1)   = mean(cond_acc);
        else
            diff_avdraws(sub,1) = mean(cond_draws-1);
            diff_avacc(sub,1)   = mean(cond_acc);
            
        end

        % sum draws for testing 
        sum_cond_draws(sub,cond) = sum(cond_draws);

        
        % we will now calculate participant points. we can use this to se how each participant performed and to compare
        % each participant with their corresponding model instant 
    
        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allsub_points(sub, cond) = (sum(cond_acc==1)*10) + (sum(cond_acc==0)*-10) + (sum(cond_draws)*-0.25);

    end 
    
    % sum participant draws for all conditions
    sum_draws(sub,1) = sum(sum_cond_draws(sub,:));

    % fit MDP beads model with two free parameters (cost-sample)
    
    % define parameters of the ideal observer
    R.error             = -10;          % cost for being wrong
    R.correct           = 10;           % reward for being correct
    R.diff              = -20; 
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.sample            = -0.25;        % the cost to sample

    for cond = 1:conditions

        init_beta       = 2; % initial value for inverse temperature
        cond_urns       = cond_data{1,sub}{1,cond}(:,4);
        cond_sequences  = subsequences{1,cond};
        cond_choiceData = subchoiceVec{1,cond};
        R.urns          = cond_urns';
        R.thisq         = R.q(cond);

        % define the objective function and get optimal parameters
        obFunc                          = @(x) mdpBeads_model3([x(1)], R, cond_sequences, cond_choiceData);
       
        lb                              = 0;
        ub                              = 10;
        options                         = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIter', 5000, 'TolFun', 1e-6);

        % Initialize best fit variables
        best_NegLL  = inf;
        best_Xfit   = [];
    
        % Run optimization with multiple starting points
        for start = 1:num_starts
            % Generate random starting points within bounds
            % Note: Adjust bounds if necessary (e.g., ub(1) should be 0 if Cs <= 0)
            starting_points = [unifrnd(lb(1), ub(1))];
    
            % Call fmincon
            [Xfit, NegLL]   = fmincon(obFunc, starting_points, [], [], [], [], lb, ub, [], options);
    
            % Update best fit if current is better
            if NegLL < best_NegLL
                best_NegLL  = NegLL;
                best_Xfit   = Xfit;
            end
        end
        
        allsub_fitted_beta_v3(sub,cond)     = best_Xfit(1);
        allsub_ll_v3(sub,cond)              = best_NegLL;

        % % run the model to participant sequences with fitted beta values
        fMDP_output = fit_subMDPBeads_v3(best_Xfit,R,cond_sequences,cond_choiceData);

        % extract output
        fit_correct_v3(sub,cond)        = fMDP_output.performance;
        fit_samples_v3(sub,cond)        = fMDP_output.avsamples;
        fit_AQs_v3{sub,cond}            = fMDP_output.actionVals;
        allsub_fit_model_v3{sub,cond}   = fMDP_output;

        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allBetaModel_points(sub, cond)  = (sum(fMDP_output.choice==1)*10) + (sum(fMDP_output.choice==0)*-10) + (sum(fMDP_output.samples)*-0.25);


        clear Xfit NegLL cond_sequences cond_choiceData R.urns R.thisq cond_urns fMDP_output

    end % end of conditions loop

end % end of subjects loop


%% plot beta parameter values and nll (Supplementary figure S1 C)

beta_0_8_model2          = allsub_fitted_beta_v3(:,1);
beta_0_6_model2          = allsub_fitted_beta_v3(:,2);


%%% plot beta values
% Calculate Means and SEM
mean_beta = [mean(beta_0_8_model2), mean(beta_0_6_model2)];
sem_beta = [std(beta_0_8_model2) / sqrt(length(beta_0_8_model2)), std(beta_0_6_model2) / sqrt(length(beta_0_6_model2))];

% Create Bar Plot
figure;
bar(1:2, mean_beta, 0.5, 'FaceColor', [0.3 0.6 0.9]); % Solid light blue bars
hold on;

% Add Error Bars
errorbar(1:2, mean_beta, sem_beta, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Customize Plot
xticks([1 2]);
xticklabels({'0.8 Condition', '0.6 Condition'});
ylabel('Beta Parameter');
xlabel('Condition');
title('Differences in Beta Parameter Between Conditions');
ylim([0 12]); % Adjust Y-axis range as needed

% Add Individual Data Points with Jitter
jitter = 0.05; % Amount of random noise
scatter(1 + (rand(size(beta_0_8_model2)) - 0.5) * jitter, beta_0_8_model2, 'filled', 'r'); % Jittered 0.8 condition
scatter(2 + (rand(size(beta_0_6_model2)) - 0.5) * jitter, beta_0_6_model2, 'filled', 'r'); % Jittered 0.6 condition

% Add Legend
% legend({'Mean ± SEM', '0.8 Condition Data', '0.6 Condition Data'}, 'Location', 'Best');
fontsize(gcf, 20, "points");
hold off;

% plot nll 

% Create Density Plot
figure;
[f1, xi1] = ksdensity(allsub_ll_v3(:, 1)); % 0.8 Condition
[f2, xi2] = ksdensity(allsub_ll_v3(:, 2)); % 0.6 Condition

% Plot filled density curves
fill(xi1, f1, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2); % Filled density for 0.8
hold on;
fill(xi2, f2, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'LineWidth', 2); % Filled density for 0.6

% Customise Plot
legend({'0.8 Condition', '0.6 Condition'}, 'Location', 'Best');
xlabel('Negative Log-Likelihood (NLL)');
ylabel('Density');
title('Density Plot of NLL Values');
fontsize(gcf, 20, "points");
hold off;

%% run ideal observer 

% loop over subjects and extarct behavioural data 
for sub = 1:nsubs

    fprintf('loading beads block data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(resultspath,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % extract blocktrial data
    [subsequences,subchoiceVec,all_data, draws_index]    = extract_blockdata(subdir,sub,task);
    
    % store in cell for each participant
    allsub_alldata{1,sub}                  = all_data;
    allsub_sequences{1,sub}                = subsequences;
    allsub_choiceVec{1,sub}                = subchoiceVec;
    allsub_drawinfo{1,sub}                 = draws_index; % col1: draw, col2: trial

    % SPLIT SUB DATA INTO CONDITIONS
    
    for cond = 1:conditions % loop over conditions
        
        tmp                                 = find(all_data(:,8) == cond);
        cond_data{1,sub}{cond}              = all_data((tmp),:);
        clear tmp
        
    end % end of condition loop

    % AVERAGE PARTICIPANT DRAWS & ACCURACY & CALCULATE POINTS 
    % create a nx1 vector (n=number of participants) with the averaged number
    % of draws for each participant.
    % This vector will be used as a covariate for the individual differences
    % analysis in SPM12.
    
    sub_draws           = all_data(:,5);
    sub_acc             = all_data(:,7);
    avdraws(sub,1)      = mean(sub_draws-1); % -1 because the last draw is the urn choice
    avacc(sub,1)        = mean(sub_acc);
    allsub_draws(:,sub) = all_data(:,5);
    allsub_conds(:,sub) = all_data(:,8);
    
    clear sub_draws sub_acc
    
    % average draws and acc for each condition
    for cond = 1:conditions
        
        tmp_cond                = cond_data{1,sub}{1,cond};
        cond_draws              = tmp_cond(:,5);
        cond_acc                = tmp_cond(:,7);
        
        if cond == 1
            easy_avdraws(sub,1) = mean(cond_draws-1);
            easy_avacc(sub,1)   = mean(cond_acc);
        else
            diff_avdraws(sub,1) = mean(cond_draws-1);
            diff_avacc(sub,1)   = mean(cond_acc);
            
        end

        % sum draws for testing 
        sum_cond_draws(sub,cond) = sum(cond_draws);

        draws_for_testing{cond} = cond_draws;

        
        % we will now calculate participant points. we can use this to se how each participant performed and to compare
        % each participant with their corresponding model instant 
    
        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allsub_points(sub, cond) = (sum(cond_acc==1)*10) + (sum(cond_acc==0)*-10) + (sum(cond_draws)*-0.25);

    end 
    
    % sum participant draws for all conditions
    sum_draws(sub,1) = sum(sum_cond_draws(sub,:));

    % % Run ideal observer 
    % % define parameters of the ideal observer
    R.beta              = 2;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
    R.error             = -10;          % cost for being wrong
    R.correct           = 10;           % reward for being correct
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.sample            = -0.25;        % the cost to sample

    % loop over conditions
    for cond = 1:conditions

        thiscond_data                   = cond_data{1,sub}{1,cond};
        thiscond_seq                    = subsequences{1,cond};
        R.cond                          = cond;
        
        io_output                       = run_MDP_BeadsIo(R,thiscond_seq,thiscond_data);

        % extract output
        allsub_io_output{1,sub}{1,cond} = io_output;
        all_ioacc(sub,cond)             = io_output.accuracy;
        all_iodraws(sub,cond)           = io_output.draws;  
        all_iopoints(sub,cond)          = io_output.points;

        clear io_output
    end 

 
end % end of subjects loop

% remove the final choice from io's number of ddraws
allmean_iodraws     =  all_iodraws-1;

%% run anova  for agents: humans, io, cnm model 

% add draws and performance in one vec
all_acc             = [easy_avacc diff_avacc];
all_draws           = [easy_avdraws diff_avdraws];
all_model_draws     = fit_samples_v2; % beta + cs model
all_model_acc       = fit_correct_v2;

anova_struct        = struct('draws_humans', all_draws, 'correct_humans', all_acc, 'correct_io', all_ioacc,...
    'draws_io', allmean_iodraws, 'draws_cnm', all_model_draws, 'correct_cnm', all_model_acc);

output_structure    = run_stats(nsubs,anova_struct);

%% exatract multicompare results for ploting with significance stars (Figure 4A)

pc_tables_draws     = output_structure.pc_tables.draws;
pc_tables_acc       = output_structure.pc_tables.acc;

% get means and SD
mean_draws          = [mean(all_draws(:,1)), mean(allmean_iodraws(:,1)), mean(fit_samples_v2(:,1));
              mean(all_draws(:,2)), mean(allmean_iodraws(:,2)), mean(fit_samples_v2(:,2))];

std_error_draws     = [std(all_draws(:,1))/sqrt(size(all_draws,1)), std(allmean_iodraws(:,1))/sqrt(size(allmean_iodraws,1)), std(fit_samples_v2(:,1))/sqrt(size(fit_samples_v2,1));
                   std(all_draws(:,2))/sqrt(size(all_draws,1)), std(allmean_iodraws(:,2))/sqrt(size(allmean_iodraws,1)), std(fit_samples_v2(:,2))/sqrt(size(fit_samples_v2,1))];

mean_acc            = [mean(all_acc(:,1)), mean(all_ioacc(:,1)), mean(fit_correct_v2(:,1));
            mean(all_acc(:,2)), mean(all_ioacc(:,2)), mean(fit_correct_v2(:,2))];

std_error_acc       = [std(all_acc(:,1))/sqrt(size(all_acc,1)), std(all_ioacc(:,1))/sqrt(size(all_ioacc,1)), std(fit_correct_v2(:,2))/sqrt(size(fit_correct_v2,1));
                 std(all_acc(:,2))/sqrt(size(all_acc,1)), std(all_ioacc(:,2))/sqrt(size(all_ioacc,1)), std(fit_correct_v2(:,2))/sqrt(size(fit_correct_v2,1))];

conditionLabels     = {'0.8 Condition', '0.6 Condition'};
legendLabels        = {'Humans', 'Ideal Observer', 'Main Model'};

plotGroupedBars_v2(mean_draws, std_error_draws, mean_acc, std_error_acc,...
    conditionLabels, legendLabels, pc_tables_draws, pc_tables_acc)


%% 
%% plot all sampling rates and accuracy for all agents (Supplementary figure S1 A)

% Calculate means and SEMs for draws
nsubs = length(easy_avdraws); % Assuming nsubs is the number of subjects

% Draws: 0.8 condition
means_draws_08 = [mean(easy_avdraws), mean(allmean_iodraws(:,1)), mean(fit_samples_v1(:,1)), mean(fit_samples_v2(:,1)), mean(fit_samples_v3(:,1))];
sems_draws_08 = [std(easy_avdraws)/sqrt(nsubs), std(allmean_iodraws(:,1))/sqrt(nsubs), std(fit_samples_v1(:,1))/sqrt(nsubs), std(fit_samples_v2(:,1))/sqrt(nsubs), std(fit_samples_v3(:,1))/sqrt(nsubs)];

% Draws: 0.6 condition
means_draws_06 = [mean(diff_avdraws), mean(allmean_iodraws(:,2)), mean(fit_samples_v1(:,2)), mean(fit_samples_v2(:,2)), mean(fit_samples_v3(:,2))];
sems_draws_06 = [std(diff_avdraws)/sqrt(nsubs), std(allmean_iodraws(:,2))/sqrt(nsubs), std(fit_samples_v1(:,2))/sqrt(nsubs), std(fit_samples_v2(:,2))/sqrt(nsubs), std(fit_samples_v3(:,2))/sqrt(nsubs)];

% Calculate means and SEMs for accuracy
% Accuracy: 0.8 condition
means_acc_08 = [mean(easy_avacc), mean(all_ioacc(:,1)), mean(fit_correct_v1(:,1)), mean(fit_correct_v2(:,1)), mean(fit_correct_v3(:,1))];
sems_acc_08 = [std(easy_avacc)/sqrt(nsubs), std(all_ioacc(:,1))/sqrt(nsubs), std(fit_correct_v1(:,1))/sqrt(nsubs), std(fit_correct_v2(:,1))/sqrt(nsubs), std(fit_correct_v3(:,1))/sqrt(nsubs)];

% Accuracy: 0.6 condition
means_acc_06 = [mean(diff_avacc), mean(all_ioacc(:,2)), mean(fit_correct_v1(:,2)), mean(fit_correct_v2(:,2)), mean(fit_correct_v3(:,2))];
sems_acc_06 = [std(diff_avacc)/sqrt(nsubs), std(all_ioacc(:,2))/sqrt(nsubs), std(fit_correct_v1(:,2))/sqrt(nsubs), std(fit_correct_v2(:,2))/sqrt(nsubs), std(fit_correct_v3(:,2))/sqrt(nsubs)];

% Define color matrix for the five groups
colors = [1 0 0; 0 0 1; 0 1 0; 1 1 0; 1 0 1]; % Red, Blue, Green, Yellow, Magenta

% Plot Draws with different colors for each group
figure;
subplot(1,2,1);
bar_handle = bar(1:5, means_draws_08); % Replace 'means_draws_08' with your data
set(bar_handle, 'FaceColor', 'flat', 'CData', colors);
hold on;
errorbar(1:5, means_draws_08, sems_draws_08, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTick', 1:5, 'XTickLabel', {'Humans', 'IO', 'Model 1', 'Model 2', 'Model 3'});
title('Draws - 0.8 Condition');
ylabel('Average Draws');
ylim([0, max([means_draws_08 + sems_draws_08, means_draws_06 + sems_draws_06])]);

subplot(1,2,2);
bar_handle = bar(1:5, means_draws_06); % Replace 'means_draws_06' with your data
set(bar_handle, 'FaceColor', 'flat', 'CData', colors);
hold on;
errorbar(1:5, means_draws_06, sems_draws_06, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTick', 1:5, 'XTickLabel', {'Humans', 'IO', 'Model 1', 'Model 2', 'Model 3'});
title('Draws - 0.6 Condition');
ylabel('Average Draws');
ylim([0, max([means_draws_08 + sems_draws_08, means_draws_06 + sems_draws_06])]);
fontsize(gcf, 17, "points");

% Plot Accuracy with different colors for each group
figure;
subplot(1,2,1);
bar_handle = bar(1:5, means_acc_08); % Replace 'means_acc_08' with your data
set(bar_handle, 'FaceColor', 'flat', 'CData', colors);
hold on;
errorbar(1:5, means_acc_08, sems_acc_08, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTick', 1:5, 'XTickLabel', {'Humans', 'IO', 'Model 1', 'Model 2', 'Model 3'});
title('Accuracy - 0.8 Condition');
ylabel('Accuracy');
ylim([0, 1]);

subplot(1,2,2);
bar_handle = bar(1:5, means_acc_06); % Replace 'means_acc_06' with your data
set(bar_handle, 'FaceColor', 'flat', 'CData', colors);
hold on;
errorbar(1:5, means_acc_06, sems_acc_06, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTick', 1:5, 'XTickLabel', {'Humans', 'IO', 'Model 1', 'Model 2', 'Model 3'});
title('Accuracy - 0.6 Condition');
ylabel('Accuracy');
ylim([0, 1]);

fontsize(gcf, 17, "points");




