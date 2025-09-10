function plotGroupedBars_v2(mean_draws, std_error_draws, mean_acc, std_error_acc, conditionLabels, legendLabels, pc_tables_draws, pc_tables_acc)
% plotGroupedBars plots grouped bar plots for draws and accuracy with significance annotations.
%
% @christinadelta June 2025
%
%   plotGroupedBars(mean_draws, std_error_draws, mean_acc, std_error_acc, 
%                   conditionLabels, legendLabels, pc_tables_draws, pc_tables_acc)
%
% Inputs:
%   - mean_draws:       Matrix of mean values for number of draws. 
%                       Each row represents a condition (e.g., 0.8, 0.6) and each
%                       column a different group (e.g., Humans, Ideal Observer, Model, etc.).
%
%   - std_error_draws:  Matrix of standard errors corresponding to mean_draws.
%
%   - mean_acc:         Matrix of mean values for accuracy.
%
%   - std_error_acc:    Matrix of standard errors corresponding to mean_acc.
%
%   - conditionLabels:  Cell array of strings for each condition (x-axis ticks).
%
%   - legendLabels:     Cell array of strings for each group (legend labels).
%
%   - pc_tables_draws:  Table of pairwise comparison results for draws.
%
%   - pc_tables_acc:    Table of pairwise comparison results for accuracy.
%
% This function creates a figure with two subplots:
%   Left: Grouped bar plot for number of draws with significance annotations.
%   Right: Grouped bar plot for accuracy with significance annotations.

%% Check that dimensions match between means and errors
if any(size(mean_draws) ~= size(std_error_draws)) || any(size(mean_acc) ~= size(std_error_acc))
    error('Mean and standard error matrices must be the same size.');
end

% Determine number of conditions and groups
nConditions = size(mean_draws, 1);
nGroups = size(mean_draws, 2);
% Generate a set of colors automatically (using the "lines" colormap)
colors = lines(nGroups);
% Create a figure with two subplots
figure;
%% Plot for Number of Draws
subplot(1,2,1);
barHandles = bar(mean_draws, 'grouped');
hold on;
% Assign colors and FaceAlpha to each group
for i = 1:nGroups
    barHandles(i).FaceColor = colors(i,:);
    barHandles(i).FaceAlpha = 0.4; % Set bar transparency to 0.4
end

% Calculate group width and add error bars
groupWidth = min(0.8, nGroups/(nGroups + 1.5));
for i = 1:nGroups
    x = (1:nConditions) - groupWidth/2 + (2*i-1) * groupWidth/(2*nGroups);
    errorbar(x, mean_draws(:,i), std_error_draws(:,i), 'k', 'linestyle', 'none');
end
title('Number of Draws');
set(gca, 'XTick', 1:nConditions, 'XTickLabel', conditionLabels);
ylabel('Mean Number of Draws');
% legend(legendLabels, 'Location', 'northeast');
% Get bar positions for significance annotations
barPositions = zeros(nConditions, nGroups);
for i = 1:nGroups
    x = (1:nConditions) - groupWidth/2 + (2*i-1) * groupWidth/(2*nGroups);
    barPositions(:,i) = x;
end

% Filter significant comparisons for draws
sig_comparisons = pc_tables_draws(pc_tables_draws.("P-value") < 0.05,:);

% Initialize cell arrays for sigstar inputs
pairs_within_condition = {};
pairs_between_conditions = {};
p_values_within = [];
p_values_between = [];
% Define group labels mapping
% agentvec: 1=humans, 2=IO, 3=CNM
% probvec: 1=0.8, 2=0.6
group_mapping = {
    'agentvec=1,probvec=1', [1,1]; % Humans, 0.8
    'agentvec=1,probvec=2', [1,2]; % Humans, 0.6
    'agentvec=2,probvec=1', [2,1]; % IO, 0.8
    'agentvec=2,probvec=2', [2,2]; % IO, 0.6
    'agentvec=3,probvec=1', [3,1]; % CNM, 0.8
    'agentvec=3,probvec=2', [3,2]; % CNM, 0.6
};

% Create a map for quick lookup
group_map = containers.Map;
for i = 1:size(group_mapping,1)
    group_map(group_mapping{i,1}) = group_mapping{i,2};
end

% Filter for between-agent within-condition and within-agent between-condition comparisons
for r = 1:height(sig_comparisons)
    groupA = string(sig_comparisons.("Group A")(r)); % Ensure string type
    groupB = string(sig_comparisons.("Group B")(r)); % Ensure string type
    p_val = sig_comparisons.("P-value")(r);
    % Debug: Print to check keys
    if ~isKey(group_map, groupA)
        disp(['Key not found for Group A: ', groupA]);
        continue;
    end
    if ~isKey(group_map, groupB)
        disp(['Key not found for Group B: ', groupB]);
        continue;
    end
    idxA = group_map(groupA);
    idxB = group_map(groupB);
    agentA = idxA(1);
    condA = idxA(2);
    agentB = idxB(1);
    condB = idxB(2);
    % Between agents within condition
    if condA == condB && agentA ~= agentB
        pairs_within_condition{end+1} = [barPositions(condA, agentA) barPositions(condB, agentB)];
        p_values_within(end+1) = p_val;
    end
    % Within agent between conditions
    if agentA == agentB && condA ~= condB
        pairs_between_conditions{end+1} = [barPositions(condA, agentA) barPositions(condB, agentB)];
        p_values_between(end+1) = p_val;
    end
end

% Add significance stars for within-condition comparisons
if ~isempty(pairs_within_condition)
    sigstar(pairs_within_condition, p_values_within); % Assumes sigstar uses solid lines or needs configuration
end

% Add significance stars and horizontal lines for between-conditions comparisons
if ~isempty(pairs_between_conditions)
    y_base = max(mean_draws(:)) + 0.2 * range(mean_draws(:)); % Higher y-position for between-condition lines
    for i = 1:length(pairs_between_conditions)
        x1 = pairs_between_conditions{i}(1);
        x2 = pairs_between_conditions{i}(2);
        p = p_values_between(i);
        y = y_base + (i-1) * 0.05 * range(mean_draws(:)); % Stack vertically for multiple comparisons
        % Draw horizontal line (solid, not dashed)
        line([x1 x2], [y y], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        % Place asterisk above the line
        x_mid = mean([x1 x2]);
        y_text = y + 0.02 * range(mean_draws(:));
        if p < 0.001
            text(x_mid, y_text, '***', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');
        elseif p < 0.01
            text(x_mid, y_text, '**', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');
        elseif p < 0.05
            text(x_mid, y_text, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');
        end
    end
end
hold off;

%% Plot for Accuracy

subplot(1,2,2);
barHandles = bar(mean_acc, 'grouped');
hold on;
% Use the same colors and FaceAlpha as for draws
for i = 1:nGroups
    barHandles(i).FaceColor = colors(i,:);
    barHandles(i).FaceAlpha = 0.4; % Set bar transparency to 0.4
end

% Add error bars for accuracy
for i = 1:nGroups
    x = (1:nConditions) - groupWidth/2 + (2*i-1) * groupWidth/(2*nGroups);
    errorbar(x, mean_acc(:,i), std_error_acc(:,i), 'k', 'linestyle', 'none');
end

title('Accuracy');
set(gca, 'XTick', 1:nConditions, 'XTickLabel', conditionLabels);
ylabel('Mean Accuracy');
% legend(legendLabels, 'Location', 'northeast');
% Get bar positions for accuracy plot (same as draws for consistency)
barPositions = zeros(nConditions, nGroups);
for i = 1:nGroups
    x = (1:nConditions) - groupWidth/2 + (2*i-1) * groupWidth/(2*nGroups);
    barPositions(:,i) = x;
end

% Filter significant comparisons for accuracy
sig_comparisons = pc_tables_acc(pc_tables_acc.("P-value") < 0.05,:);
% Initialize cell arrays for sigstar inputs
pairs_within_condition = {};
pairs_between_conditions = {};
p_values_within = [];
p_values_between = [];
% Use the same group mapping as for draws
% Filter for between-agent within-condition and within-agent between-condition comparisons
for r = 1:height(sig_comparisons)
    groupA = string(sig_comparisons.("Group A")(r)); % Ensure string type
    groupB = string(sig_comparisons.("Group B")(r)); % Ensure string type
    p_val = sig_comparisons.("P-value")(r);
    % Debug: Print to check keys
    if ~isKey(group_map, groupA)
        disp(['Key not found for Group A: ', groupA]);
        continue;
    end
    if ~isKey(group_map, groupB)
        disp(['Key not found for Group B: ', groupB]);
        continue;
    end
    idxA = group_map(groupA);
    idxB = group_map(groupB);
    agentA = idxA(1);
    condA = idxA(2);
    agentB = idxB(1);
    condB = idxB(2);
    % Between agents within condition
    if condA == condB && agentA ~= agentB
        pairs_within_condition{end+1} = [barPositions(condA, agentA) barPositions(condB, agentB)];
        p_values_within(end+1) = p_val;
    end
    % Within agent between conditions
    if agentA == agentB && condA ~= condB
        pairs_between_conditions{end+1} = [barPositions(condA, agentA) barPositions(condB, agentB)];
        p_values_between(end+1) = p_val;
    end
end
% Add significance stars for within-condition comparisons
if ~isempty(pairs_within_condition)
    sigstar(pairs_within_condition, p_values_within); % Assumes sigstar uses solid lines or needs configuration
end

% Add significance stars and horizontal lines for between-conditions comparisons
if ~isempty(pairs_between_conditions)
    y_base = max(mean_acc(:)) + 0.2 * range(mean_acc(:)); % Higher y-position for between-condition lines
    for i = 1:length(pairs_between_conditions)
        x1 = pairs_between_conditions{i}(1);
        x2 = pairs_between_conditions{i}(2);
        p = p_values_between(i);
        y = y_base + (i-1) * 0.05 * range(mean_acc(:)); % Stack vertically for multiple comparisons
        % Draw horizontal line (solid, not dashed)
        line([x1 x2], [y y], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        % Place asterisk above the line
        x_mid = mean([x1 x2]);
        y_text = y + 0.02 * range(mean_acc(:));
        if p < 0.001
            text(x_mid, y_text, '***', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');
        elseif p < 0.01
            text(x_mid, y_text, '**', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');
        elseif p < 0.05
            text(x_mid, y_text, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k');
        end
    end
end

hold off;
% Adjust figure size if desired
set(gcf, 'Position', [100, 100, 1000, 500]);
fontsize(gcf, 20, "points");

end