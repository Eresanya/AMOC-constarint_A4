clear; clc;

%% Observational correlations
r1_obs = 0.57;
r2_obs = 0.52;

% Sample size and significance threshold
n = 21;%no of sampled years
alpha = 0.05;
t_crit = tinv(1 - alpha/2, n - 2);  % two-tailed
r_sig = t_crit / sqrt(t_crit^2 + n - 2);  % minimum r for significance

% Model correlations and names
r1 = [-0.03, 0.10, 0.13, 0.17, 0.33, 0.82, 0.61, 0.41, 0.62, 0.07, -0.36, 0.54, ...
      0.32, -0.12, 0.04, 0.73, 0.13, -0.29, 0.17];
r2 = [0.21, 0.05, 0.51, 0.13, 0.85, 0.71, 0.62, 0.91, -0.37, -0.04, 0.59, 0.22, ...
      -0.17, 0.15, 0.53, 0.77, 0.02, -0.25, 0.16];

models_name = {...
    'ACCESS-CM2', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5', ...
    'CESM2', 'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV2', ...
    'CMCC-CM2-HR4', 'E3SM-1-0', 'E3SM-1-1', 'FGOALS-f3-L', ...
    'GISS-E2-1-G', 'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', ...
    'MIROC6', 'NESM3', 'NorESM2-LM'};

% Distance to observation
dist = sqrt((r1 - r1_obs).^2 + (r2 - r2_obs).^2);
[sortedDist, sortIdx] = sort(dist);
sortedModels = models_name(sortIdx);
topN = 6;
radius = sortedDist(topN);

% Setup figure
figure('Color', 'w', 'Position', [100, 100, 1400, 600]);

% Colors and markers
colors = lines(length(r1));
defaultMarkers = {'o','s','^','d','v','>','<','p','h','x','+','*','.','o','s','^','d','v','>'};
markers = defaultMarkers(1:length(r1));

% === PANEL 1: Scatter Plot ===
subplot(1, 2, 1); hold on;

% Highlight region: r > 0.5 and significant
patch([0.5, 1, 1, 0.5], [0.5, 0.5, 1, 1], [0.85 1 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Draw observation
plot(r1_obs, r2_obs, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
text(r1_obs + 0.02, r2_obs, 'Observation', 'FontSize', 11, 'FontWeight', 'bold');

% Circle around top-N
theta = linspace(0, 2*pi, 300);
plot(r1_obs + radius*cos(theta), r2_obs + radius*sin(theta), '--k', 'LineWidth', 2);

% Threshold lines
xline(0.5, ':k', 'LineWidth', 1.5);
yline(0.5, ':k', 'LineWidth', 1.5);
xline(r_sig, '--b', 'LineWidth', 1); yline(r_sig, '--b', 'LineWidth', 1);
xline(-r_sig, '--b', 'LineWidth', 1); yline(-r_sig, '--b', 'LineWidth', 1);
plot([-1 1], [-1 1], 'k--', 'LineWidth', 1);  % Diagonal

% Plot each model
for i = 1:length(r1)
    isTop = ismember(i, sortIdx(1:topN));
    isHighCorr = (r1(i) > 0.5) && (r2(i) > 0.5);
    isSignificant = (abs(r1(i)) > r_sig) && (abs(r2(i)) > r_sig);

    faceColor = colors(i,:);
    edgeColor = 'k';
    lineW = 1.5;

    if isTop
        faceColor = [1 0 0];  % Red fill
    end
    if isHighCorr && isSignificant
        edgeColor = [0 0.6 0];  % Green border
        lineW = 2.5;
    end

    plot(r1(i), r2(i), markers{i}, ...
        'MarkerSize', 10, ...
        'MarkerEdgeColor', edgeColor, ...
        'MarkerFaceColor', faceColor, ...
        'LineWidth', lineW);

    text(r1(i)+0.015, r2(i), models_name{i}, 'FontSize', 9, 'Interpreter', 'none');
end

xlabel('r_1: AMOC vs ZOS Coast Index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('r_2: AMOC vs ZOS Index', 'FontSize', 13, 'FontWeight', 'bold');
title('(a) Model Correlations vs Observation', 'FontSize', 14, 'FontWeight', 'bold');
axis equal; xlim([-1 1]); ylim([-1 1]);
box on;

% Annotation for r_sig
annotation('textbox', [0.15, 0.87, 0.2, 0.05], 'String', ...
    sprintf('95%% significance: |r| > %.2f (n = %d)', r_sig, n), ...
    'EdgeColor', 'none', 'FontSize', 10, 'Color', 'b');

% === PANEL 2: GROUPED BAR PLOT OF r1 and r2 ===
subplot(1, 2, 2);

% Sort data by distance
r1_sorted = r1(sortIdx);
r2_sorted = r2(sortIdx);
model_labels = models_name(sortIdx);

% Bar data
bar_data = [r1_sorted(:), r2_sorted(:)];

% Bar chart
hb = bar(bar_data, 'grouped');
hb(1).FaceColor = [0.2 0.4 0.8];  % r1
hb(2).FaceColor = [0.9 0.4 0.2];  % r2

% X labels
xticks(1:length(r1));
xticklabels(model_labels);
xtickangle(45);

% Threshold lines (and save handles for legend)
hold on;
h1 = yline(r_sig, '--b', 'LineWidth', 1); % 95% significance
h2 = yline(0.5, ':k', 'LineWidth', 1.5);  % r > 0.5

% Dummy handles for legend
h3 = plot(NaN, NaN, '--b', 'LineWidth', 1);   % 95% line legend
h4 = plot(NaN, NaN, ':k', 'LineWidth', 1.5);  % r > 0.5 line legend

ylabel('Correlation (r)', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) Correlation Strength per Model', 'FontSize', 14, 'FontWeight', 'bold');
ylim([-1 1]);
box on;

% Updated legend to include lines
legend([hb(1), hb(2), h3, h4], ...
    {'r_1: Coast Index', 'r_2: ZOS Index', ...
     '95% significance (|r| > 0.43)', 'r > 0.5'}, ...
    'Location', 'northwest', 'FontSize', 9);
