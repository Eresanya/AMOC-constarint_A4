%% ============================================================
%  AMOC (26.5°N) – CMIP6 SSP585 & SSP126 + TG Constraint
%  1850–2100 | 5-year running mean | Baseline: 1850–2014
% =============================================================
clc; clear;
%% -------------------- USER PARAMETERS -----------------------
ncfile585 = 'CMIP6_AMOC_ssp585_26N.nc';
ncfile126 = 'CMIP6_AMOC_ssp126_26N.nc';
%%
% load('TG_adjstd_SSH_updated.mat');% TG detrended 
load('TG_adjstd_SSH.mat'); %Both SSH and TG detrended
%%
years_full = 1850:2100;
hist_ref_idx = years_full <= 2014;
fut_idx = years_full > 2014;
win = 5;   % running mean

models_target = { ...
    'BCC-CSM2-MR','CAMS-CSM1-0','CanESM5','CESM2', ...
    'CESM2-WACCM','E3SM-1-1','FGOALS-f3-L','GISS-E2-1-G', ...
    'INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC6','NESM3','NorESM2-LM'};

% model_constraint = { ...
%     'CanESM5','CESM2-WACCM','CNRM-CM6-1-HR','FIO-ESM-2-0','HadGEM3-GC31-LL',...
%     'BCC-CSM2-MR','CAMS-CSM1-0','CESM2','CNRM-CM6-1','HadGEM3-GC31-MM', ...
%     'MPI-ESM1-2-HR','SAM0-UNICON','MRI-ESM2-0','NESM3','NorESM2-LM','UKESM1-0-LL'};

model_constraint = { ...
    'CanESM5','CESM2-WACCM','CNRM-CM6-1-HR','FIO-ESM-2-0','HadGEM3-GC31-LL',...
    'BCC-CSM2-MR','CAMS-CSM1-0','CESM2','CNRM-CM6-1','HadGEM3-GC31-MM', ...
    'MPI-ESM1-2-HR','SAM0-UNICON','NESM3','UKESM1-0-LL'};
%
% %% -------------------- SCENARIO PROCESSOR --------------------
processScenario = @(ncfile) ...
    (function_process(ncfile, models_target, model_constraint, ...
                      years_full, hist_ref_idx, win));
[ens585_mean, ens585_std, trend585_hist, trend585_fut] = ...
    processScenario(ncfile585);

[ens126_mean, ens126_std, trend126_hist, trend126_fut] = ...
    processScenario(ncfile126);

%% -------------------- TG PROCESSING -------------------------
tg_years = 1880:2025;
tg_mask = tg_years >= 1850 & tg_years <= 2100;


% Smooth full series FIRST
% tg_amoc_smooth = movmean(amoc_obs_from_tg_full, win, 'omitnan');
tg_amoc_smooth = movmean(tg_grad_adj_full, win, 'omitnan');
tg_amoc_smooth = tg_amoc_smooth/100;
tg_amoc_smooth = tg_amoc_smooth - mean(tg_amoc_smooth);
% Define analysis period
common_years = 1940:2014;

% Indices
idx_tg_common = ismember(tg_years, common_years);
idx_cmip = ismember(years_full, common_years);

% Extract aligned data
tg_hist_overlap = tg_amoc_smooth(idx_tg_common);
cmip6_hist_overlap = ens585_mean(idx_cmip);

% Ensure column vectors
tg_hist_overlap    = tg_hist_overlap(:);
cmip6_hist_overlap = cmip6_hist_overlap(:);

% b = regress(cmip6_hist_overlap, tg_hist_overlap);
% tg_amoc_smooth = tg_amoc_smooth * b;
%-------------standardize--------------%
% tg_amoc_std = (tg_hist_overlap - mean(tg_hist_overlap,'omitnan')) ./ ...
%                std(tg_hist_overlap,'omitnan');
% 
% cmip6_hist_std = (cmip6_hist_overlap - mean(cmip6_hist_overlap,'omitnan')) ./ ...
               % std(cmip6_hist_overlap,'omitnan');

% Detrend
% tg_detrend   = detrend(tg_amoc_std);
% cmip_detrend = detrend(cmip6_hist_std );

% Correlation
% RAW correlation (trend + variability)
[r_hist, p_hist] = corr(cmip6_hist_overlap, tg_hist_overlap,'Rows','complete');

% DETRENDED correlation (variability only)
[r_detr, p_detr] = corr(detrend(cmip6_hist_overlap), ...
                        detrend(tg_hist_overlap), ...
                        'Rows','complete');
% [r_hist, p_hist] = corr(cmip6_hist_std, tg_amoc_std, 'Rows','complete');
% [r_hist, p_hist] = corr(cmip6_hist_overlap, tg_hist_overlap, 'Rows','complete');
% [r_hist, p_hist] = corr(cmip_detrend, tg_detrend, 'Rows','complete');
%% ================= LAG CORRELATION =================

max_lag = 10; % years 
[lags, ~] = meshgrid(-max_lag:max_lag, 1); % just for reference

% [r_lag, lag_vals] = xcorr(cmip_detrend, tg_detrend, max_lag, 'coeff');
% [r_lag, lag_vals] = xcorr(cmip6_hist_std, tg_amoc_std, max_lag, 'coeff');
[r_lag, lag_vals] = xcorr(detrend(cmip6_hist_overlap), ...
                         detrend(tg_hist_overlap), ...
                         max_lag, 'coeff');
% Find maximum absolute correlation
[~, i_max] = max(abs(r_lag));
best_r   = r_lag(i_max);
best_lag = lag_vals(i_max);

if best_lag > 0
    disp('CMIP6 leads TG')
elseif best_lag < 0
    disp('TG leads CMIP6')
else
    disp('No lag (simultaneous)')
end
%% ================= LAG CORRELATION PLOT =================

figure('Color','w');
plot(lag_vals, r_lag, 'k', 'LineWidth',1.8); hold on
yline(0,'k--')
xline(0,'k--')

% Highlight best lag
plot(best_lag, best_r, 'ro', 'MarkerFaceColor','r')

xlabel('Lag (years)')
ylabel('Correlation (r)')
title('Lag correlation: CMIP6 vs TG (detrended)')
grid on

text(best_lag, best_r, ...
     sprintf('  r = %.2f at lag = %d', best_r, best_lag), ...
     'VerticalAlignment','bottom');

fprintf('\n=====================================\n')
fprintf(' Lag Correlation (Detrended)\n')
fprintf('=====================================\n')
fprintf(' Max |r| = %.3f at lag = %d years\n', best_r, best_lag)
fprintf('=====================================\n\n')
%% ================= HISTORICAL CORRELATION PRINT =================

fprintf('\n=====================================\n')
fprintf(' TG–CMIP6 Historical Correlation\n')
fprintf(' Period: 1940–2014\n')
fprintf('=====================================\n')
fprintf(' r = %.3f\n', r_hist)
fprintf(' p = %.4f\n', p_hist)
fprintf('Detrended r = %.3f\n', r_detr)
fprintf('Detrended p = %.4f\n', p_detr)
fprintf('=====================================\n\n')


%% ================= FUTURE OVERLAP CORRELATION =================
% Period: 1990–2025 (TG vs SSP585 & SSP126)

future_years =1990:2025;

% Indices
idx_tg_fut   = ismember(tg_years, future_years);
idx_cmip_fut = ismember(years_full, future_years);

% Extract data
tg_fut   = tg_amoc_smooth(idx_tg_fut);
ssp585_fut = ens585_mean(idx_cmip_fut);
ssp126_fut = ens126_mean(idx_cmip_fut);

% Ensure column vectors
tg_fut     = tg_fut(:);
ssp585_fut = ssp585_fut(:);
ssp126_fut = ssp126_fut(:);

% ================= RAW CORRELATION =================
[r585_fut, p585_fut] = corr(ssp585_fut, tg_fut, 'Rows','complete');
[r126_fut, p126_fut] = corr(ssp126_fut, tg_fut, 'Rows','complete');

% ================= DETRENDED CORRELATION =================
[r585_fut_d, p585_fut_d] = corr(detrend(ssp585_fut), detrend(tg_fut), 'Rows','complete');
[r126_fut_d, p126_fut_d] = corr(detrend(ssp126_fut), detrend(tg_fut), 'Rows','complete');

%% ================= LAG CORRELATION (2014–2025) =================

max_lag = 10;

% SSP585
[r_lag_585, lag_585] = xcorr(detrend(ssp585_fut), detrend(tg_fut), max_lag, 'coeff');
[~, i585] = max(abs(r_lag_585));
best_r_585   = r_lag_585(i585);
best_lag_585 = lag_585(i585);

% SSP126
[r_lag_126, lag_126] = xcorr(detrend(ssp126_fut), detrend(tg_fut), max_lag, 'coeff');
[~, i126] = max(abs(r_lag_126));
best_r_126   = r_lag_126(i126);
best_lag_126 = lag_126(i126);

%% ================= PRINT RESULTS =================
fprintf('\n=====================================\n')
fprintf(' TG vs SSP Correlation (1990–2025)\n')
fprintf('=====================================\n')

fprintf('\nSSP585:\n')
fprintf(' Raw r = %.3f (p = %.4f)\n', r585_fut, p585_fut)
fprintf(' Detrended r = %.3f (p = %.4f)\n', r585_fut_d, p585_fut_d)
fprintf(' Max lag r = %.3f at lag = %d yr\n', best_r_585, best_lag_585)

fprintf('\nSSP126:\n')
fprintf(' Raw r = %.3f (p = %.4f)\n', r126_fut, p126_fut)
fprintf(' Detrended r = %.3f (p = %.4f)\n', r126_fut_d, p126_fut_d)
fprintf(' Max lag r = %.3f at lag = %d yr\n', best_r_126, best_lag_126)

fprintf('=====================================\n\n')
%% ===================== PUBLICATION FIGURE ====================

figure('Color','w','Units','inches','Position',[1 1 9 5]);
hold on;box on;grid on;

% --- Historical shading (SSP585 only) ---
hHistSpread = fill([years_full(hist_ref_idx) fliplr(years_full(hist_ref_idx))], ...
     [ens585_mean(hist_ref_idx)+ens585_std(hist_ref_idx) ...
     fliplr(ens585_mean(hist_ref_idx)-ens585_std(hist_ref_idx))], ...
     [0.75 0.85 1], 'EdgeColor','none','FaceAlpha',0.25);

% --- Future shading SSP585 ---
h585Spread = fill([years_full(fut_idx) fliplr(years_full(fut_idx))], ...
     [ens585_mean(fut_idx)+ens585_std(fut_idx) ...
      fliplr(ens585_mean(fut_idx)-ens585_std(fut_idx))], ...
     [1 0.7 0.7], 'EdgeColor','none','FaceAlpha',0.25);

% --- Future shading SSP126 ---
h126Spread = fill([years_full(fut_idx) fliplr(years_full(fut_idx))], ...
     [ens126_mean(fut_idx)+ens126_std(fut_idx) ...
      fliplr(ens126_mean(fut_idx)-ens126_std(fut_idx))], ...
     [0.7 0.85 1], 'EdgeColor','none','FaceAlpha',0.25);

% ================= ENSEMBLE MEANS =================

% ---- Historical (shared baseline) ----
hHist = plot(years_full(hist_ref_idx), ...
             ens585_mean(hist_ref_idx), ...
             'k','LineWidth',2.2);

% ---- SSP585 Future (red) ----
h585 = plot(years_full(fut_idx), ...
            ens585_mean(fut_idx), ...
            'Color',[0.8 0 0], ...
            'LineWidth',2.2);

% ---- SSP126 Future (blue) ----
h126 = plot(years_full(fut_idx), ...
            ens126_mean(fut_idx), ...
            'Color',[0 0.25 0.8], ...
            'LineWidth',2.2);
% --- TG overlay (1940–2025 only) ---
tg_plot_years = tg_years(tg_mask);

% idx_tg_plot = tg_plot_years >= 1940 & tg_plot_years <= 2025;
idx_tg_plot = tg_plot_years >= 1940 & tg_plot_years <= 2014;
% hTG = plot(tg_plot_years(idx_tg_plot), ...
%            tg_anom(idx_tg_plot), ...
hTG = plot(tg_years(idx_tg_plot), ...
           tg_amoc_smooth(idx_tg_plot), ...
           'Color',[0.1 0.1 0.1], ...
           'LineStyle','--','LineWidth',1.8);

% --- Historical/Future divider ---
xline(2014,'k--','LineWidth',1);

% --- Axis styling (Nature style) ---
xlim([1850 2100])
ylim([-15 5])
set(gca,'FontName','Helvetica',...
        'FontSize',10,...
        'LineWidth',1,...
        'TickDir','in',...
        'Box','off');

xlabel('Year')
ylabel('AMOC anomaly (Sv)')

legend([hHist, h585, h126, hTG], ...
       {'Historical (mean ± spread)', ...
        'SSP585 (mean ± spread)', ...
        'SSP126 (mean ± spread)', ...
        'Tide-gauge AMOC'}, ...
       'Location','northeast','Box','off');
% ================= TREND ANNOTATION =================

% Position for text (bottom-left clean placement)
x_text = 1860;
y_text_585 = -13;
y_text_126 = -14.2;

% SSP585 future trend (red)
text(x_text, y_text_585, ...
     sprintf('SSP585 (2015–2100): %.2f Sv decade^{-1}', trend585_fut), ...
     'Color',[0.8 0 0], ...
     'FontSize',10, ...
     'FontName','Helvetica');

% SSP126 future trend (blue)
text(x_text, y_text_126, ...
     sprintf('SSP126 (2015–2100): %.2f Sv decade^{-1}', trend126_fut), ...
     'Color',[0 0.25 0.8], ...
     'FontSize',10, ...
     'FontName','Helvetica');

% ================= CORRELATION ANNOTATION =================

text(1860, -11.5, ...
     sprintf('Historical r = %.2f (p = %.4f)', r_hist, p_hist), ...
     'Color','k', ...
     'FontSize',10, ...
     'FontName','Helvetica');
title('AMOC evolution (26.5°N), 1850–2100', ...
      'FontWeight','normal')

text(1860, -10.5, ...
     sprintf('Max r = %.2f (lag = +%d yr)', best_r, best_lag), ...
     'Color','k', ...
     'FontSize',10, ...
     'FontName','Helvetica');

%% ================= AMOC STATISTICS =================

% Historical statistics (1850–2014)
hist_mean = mean(ens585_mean(hist_ref_idx),'omitnan');
hist_std  = std(ens585_mean(hist_ref_idx),'omitnan');

% Future statistics
ssp585_mean = mean(ens585_mean(fut_idx),'omitnan');
ssp585_std  = std(ens585_mean(fut_idx),'omitnan');

ssp126_mean = mean(ens126_mean(fut_idx),'omitnan');
ssp126_std  = std(ens126_mean(fut_idx),'omitnan');

% Percentage changes relative to historical
% pct_mean_585  = (ssp585_mean - hist_mean) / abs(hist_mean) * 100;
% pct_mean_126  = (ssp126_mean - hist_mean) / abs(hist_mean) * 100;

pct_trend_585 = (trend585_fut - trend585_hist) / abs(trend585_hist) * 100;
pct_trend_126 = (trend126_fut - trend126_hist) / abs(trend126_hist) * 100;

fprintf('\n==============================\n')
fprintf(' AMOC Statistics (1850–2100)\n')
fprintf('==============================\n')

fprintf('Historical Period (1850–2014):\n')
fprintf(' Constrained: mean = %.3f Sv, std = %.3f Sv, trend = %.3f Sv/decade\n\n', ...
        hist_mean, hist_std, trend585_hist)

fprintf('Future Period (2015–2100):\n')
fprintf(' SSP585: mean = %.3f Sv, std = %.3f Sv, trend = %.3f Sv/decade\n', ...
        ssp585_mean, ssp585_std, trend585_fut)

fprintf(' SSP126: mean = %.3f Sv, std = %.3f Sv, trend = %.3f Sv/decade\n\n', ...
        ssp126_mean, ssp126_std, trend126_fut)

% % fprintf('Percentage Change (Future vs Historical):\n')
% % fprintf(' SSP585: mean = %.1f%%, trend = %.1f%%\n', ...
% %         pct_mean_585, pct_trend_585)
% % 
% % fprintf(' SSP126: mean = %.1f%%, trend = %.1f%%\n', ...
% %         pct_mean_126, pct_trend_126)

fprintf('==============================\n\n')
function [ens_mean, ens_std, trend_hist, trend_fut] = ...
    function_process(ncfile, models_target, model_constraint, ...
                     years_full, hist_ref_idx, win)

AMOC_raw = ncread(ncfile,'AMOC_CMIP6');
AMOC = squeeze(permute(AMOC_raw,[5 3 1 2 4]));
AMOC(~isfinite(AMOC)) = NaN;
AMOC(abs(AMOC)>1e4) = NaN;

model_names = strtrim(cellstr(ncread(ncfile,'model')'));
valid_idx = ismember(model_names, models_target);
selected_models = model_names(valid_idx);

AMOC_sel = AMOC(valid_idx,:,:);
AMOC_mean = squeeze(mean(AMOC_sel,2,'omitnan'));

% Baseline anomaly
AMOC_anom = AMOC_mean - ...
    mean(AMOC_mean(:,hist_ref_idx),2,'omitnan');

AMOC_smooth = movmean(AMOC_anom,win,2,'omitnan');

idx_constr = ismember(selected_models, model_constraint);
AMOC_constr = AMOC_smooth(idx_constr,:);

ens_mean = mean(AMOC_constr,1,'omitnan');
ens_std  = std(AMOC_constr,0,1,'omitnan');

fut_idx = years_full > 2014;

p = polyfit(years_full(hist_ref_idx),ens_mean(hist_ref_idx),1);
trend_hist = p(1)*10;

p = polyfit(years_full(fut_idx),ens_mean(fut_idx),1);
trend_fut = p(1)*10;
end
