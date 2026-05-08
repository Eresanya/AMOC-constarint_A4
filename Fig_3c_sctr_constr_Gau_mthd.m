%% ========================================================================
% OVERVIEW OF ANALYSIS PIPELINE
%
% This script quantifies the consistency between CMIP6 model-simulated
% and observed relationships between the Atlantic Meridional Overturning
% Circulation (AMOC) and coastal sea-level variability (ZOS).
%
% The workflow consists of two main components:
%
% (1) Estimation of observational uncertainty:
%     A block bootstrap method is applied to observed AMOC and ZOS
%     time series to account for temporal autocorrelation. For each
%     bootstrap realization, lagged AMOC–ZOS correlations are recomputed
%     exactly as in the original analysis, and the maximum lagged
%     correlation is extracted. The standard deviation of the resulting
%     distribution provides an estimate of observational uncertainty (σ).
%
% (2) Observational constraint on CMIP6 models:
%     The observed maximum lagged correlation is compared to model-
%     simulated correlations from CMIP6 models. A Gaussian likelihood
%     function, centered on the observed value and scaled by the bootstrap
%     uncertainty (σ), is used to assign a weight to each model reflecting
%     its consistency with observations.
%
%     These weights are used to:
%       - Compute a weighted ensemble mean correlation
%       - Rank models by observational consistency
%       - Identify well-constrained models relative to a threshold
%
% Finally, results are saved in both MATLAB and CSV formats, and a
% publication-ready visualization is generated showing model performance,
% observational constraints, and relative agreement with observations.
%
% ========================================================================
clc; clear

%% =========================================================
%% 1. LOAD OBSERVATIONAL DATA
%% =========================================================
load('obs_AMOC_ZOS_all.mat')
load('AMOC_ZOS_scalar_corr_1850_2014.mat');  % contains best_r_CI

%% =========================================================
%% 2. BOOTSTRAP ESTIMATE OF OBSERVATIONAL UNCERTAINTY (SIGMA)
%% =========================================================
nboot = 1000;
boot_r = nan(nboot,1);

N = length(amoc_s);

for b = 1:nboot
    
    % --- Block bootstrap (preserve autocorrelation) ---
    block_size = 5;  % matches smoothing scale
    idx = [];
    
    while length(idx) < N
        start = randi(N - block_size + 1);
        idx = [idx, start:(start+block_size-1)];
    end
    
    idx = idx(1:N);
    
    am_boot = amoc_s(idx);
    zos_boot = zos_s(idx);
    
    % --- Recompute lagged correlation EXACTLY like original ---
    corrVals_boot = nan(size(lags));
    
    for k = 1:length(lags)
        lag = lags(k);
        am = am_boot(1:end-lag);
        zo = zos_boot(1+lag:end);
        
        R = corrcoef(am, zo, 'rows', 'complete');
        corrVals_boot(k) = R(1,2);
    end
    
    % Store max correlation
    boot_r(b) = max(corrVals_boot);
end

% --- Observational uncertainty ---
sigma = std(boot_r, 'omitnan');
fprintf('Estimated sigma (bootstrap) = %.3f\n', sigma);

%% =========================================================
%% 3. OBSERVATIONAL + MODEL DATA
%% =========================================================
r_obs = maxCorr;        % Observed maximum lagged correlation (scalar)
r_model = best_r_CI;    % Model correlations (vector)

model_name = { ...
    'BCC-CSM2-MR','CAMS-CSM1-0','CanESM5','CESM2','CESM2-WACCM','CNRM-CM6-1', ...
    'CNRM-CM6-1-HR','E3SM-1-1','FIO-ESM-2-0','HadGEM3-GC31-LL', ...
    'HadGEM3-GC31-MM','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR', ...
    'MPI-ESM1-2-HR','MRI-ESM2-0','NESM3','NorESM2-LM', ...
    'SAM0-UNICON','UKESM1-0-LL'};

%% =========================================================
%% 4. OBSERVATIONAL CONSTRAINT WEIGHTING (USING BOOTSTRAP SIGMA)
%% =========================================================
weights = exp(-((r_obs - r_model).^2) / (2 * sigma^2));

weighted_avg_r = sum(weights .* r_model) / sum(weights);

fprintf('Weighted mean correlation (Gaussian constraint): %.3f\n', weighted_avg_r);

%% =========================================================
%% 5. MODEL SELECTION THRESHOLD
%% =========================================================
weight_threshold = 0.3 * max(weights);

%% =========================================================
%% 6. SORT MODELS
%% =========================================================
[r_sorted, sortIdx] = sort(r_model, 'descend');
weights_sorted = weights(sortIdx);
model_sorted = model_name(sortIdx);

%% =========================================================
%% 7. SAVE RESULTS
%% =========================================================
model_weights.model  = model_sorted(:)';
model_weights.r_model = r_sorted(:)';
model_weights.weight  = weights_sorted(:)';

save('model_exponential_weights.mat', 'model_weights');

T = table(model_sorted(:), r_sorted(:), weights_sorted(:), ...
    'VariableNames', {'Model','Correlation_r','GaussianWeight'});

writetable(T, 'model_exponential_weights.csv');

fprintf('Model weights saved (MAT + CSV format)\n');

%% =========================================================
%% 8. VISUALIZATION
%% =========================================================
col_constraint   = [0 0.45 0.74];
col_unconstraint = [0.85 0.33 0.10];

figure('Color','w','Position',[100,100,1500,600]);
hold on; box on;

max_bar_width = max(r_sorted) * 1.05;

%% =========================================================
%% 9. BAR PLOT
%% =========================================================
for i = 1:length(r_sorted)

    if weights_sorted(i) > weight_threshold
        col = col_constraint * weights_sorted(i) + [0.8 0.8 0.8] * (1 - weights_sorted(i));
    else
        col = col_unconstraint;
    end

    barh(i, r_sorted(i), ...
        'FaceColor', col, ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.2);
end

%% =========================================================
%% 10. OBSERVATIONAL LINE
%% =========================================================
xline(r_obs, 'r--', 'LineWidth', 2);

text(r_obs + 0.01 * max_bar_width, length(r_sorted) + 0.5, ...
    sprintf('Observed r = %.2f', r_obs), ...
    'Color', 'r', 'FontWeight', 'bold');

%% =========================================================
%% 11. AXES
%% =========================================================
yticks(1:length(r_sorted));
yticklabels(model_sorted);
ytickangle(0);

xlabel('AMOC–ZOS correlation (r)', 'FontWeight','bold');
ylabel('CMIP6 Models', 'FontWeight','bold');

title(sprintf('(b) Models Ranked by Observational Consistency (Weighted Mean r = %.3f)', ...
    weighted_avg_r), ...
    'FontWeight','bold');

grid on;

ax = gca;
ax.GridColor = [0.8 0.8 0.8];
ax.GridAlpha = 0.5;

xlim([-0.05, max_bar_width * 1.1]);

set(findall(gcf,'Type','axes'),'FontSize',12,'LineWidth',1.2);

%% =========================================================
%% 12. COLORBAR
%% =========================================================
colormap(gca, [linspace(0.8, col_constraint(1), 256)' ...
               linspace(0.8, col_constraint(2), 256)' ...
               linspace(0.8, col_constraint(3), 256)']);

c = colorbar('eastoutside');
c.Label.String = 'Gaussian weight (consistency with observations)';
c.Label.FontWeight = 'bold';
c.FontSize = 11;

c.Ticks = [0 1];
c.TickLabels = {'Low','High'};

%% =========================================================
%% 13. LEGEND
%% =========================================================
h1 = barh(NaN, NaN, 'FaceColor', col_constraint, 'EdgeColor','k');
h2 = barh(NaN, NaN, 'FaceColor', col_unconstraint, 'EdgeColor','k');
h3 = xline(NaN,'r--','LineWidth',2);

legend([h1 h2 h3], {'Constrained', 'Others', 'Observation'}, ...
    'Location','eastoutside', ...
    'FontSize',10);