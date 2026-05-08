%% ========================================================================
% This script reproduces multi-panel time-series comparisons between
% Atlantic Meridional Overturning Circulation (AMOC) variability and
% coastal dynamic sea-level (ZOS) variability using CMIP6 climate models
% and observational datasets.
%
% The primary objective is to evaluate whether CMIP6 models reproduce the
% observed temporal relationship between AMOC strength and coastal sea-level
% variability over the North Atlantic, and to assess the consistency of
% this relationship across individual models and the multi-model ensemble.
%
% The script performs the following major steps:
%
% 1. Loads CMIP6 AMOC and coastal ZOS diagnostics
%    - Historical simulations from CMIP6 models are loaded from a pre-
%      processed results structure.
%    - Two analysis periods are available:
%         • 1850–2014 (historical simulations)
%         • 1993–2014 (satellite era)
%
% 2. Aligns model time series
%    - All model AMOC and coastal ZOS time series are truncated to a common
%      temporal length to ensure direct inter-model comparison.
%
% 3. Normalizes model variability
%    - Each model AMOC and ZOS time series is standardized by removing its
%      mean and dividing by its standard deviation.
%    - This normalization emphasizes temporal co-variability rather than
%      differences in absolute magnitude among models.
%
% 4. Computes CMIP6 ensemble statistics
%    - Multi-model ensemble means are calculated for AMOC and coastal ZOS.
%    - Ensemble spread is quantified using the inter-model standard
%      deviation.
%    - The correlation between ensemble-mean AMOC and coastal ZOS
%      variability is computed.
%
% 5. Loads and processes observational datasets
%    - Observed AMOC and coastal sea-level indices are loaded for the
%      satellite period (1993–2014).
%    - Missing values are removed and overlapping observational years are
%      retained for comparison with CMIP6 simulations.
%
% 6. Evaluates model–observation agreement
%    - Correlations are computed between:
%         • CMIP6 ensemble AMOC and ensemble coastal ZOS
%         • Observed and simulated AMOC variability
%         • Observed and simulated coastal ZOS variability
%
% 7. Produces publication-quality multi-panel figures
%
%    Ensemble panel:
%       • Displays ensemble-mean AMOC and coastal ZOS variability
%       • Shows ensemble spread as shaded uncertainty envelopes
%       • Overlays observational AMOC and coastal ZOS records
%       • Annotates ensemble and observation-model correlations
%
%    Individual model panels:
%       • Show normalized AMOC and coastal ZOS variability for each CMIP6
%         model
%       • Display the correlation between AMOC and coastal ZOS variability
%         within each model
%
% 8. Visualization features
%    - Consistent axis scaling across all panels
%    - Publication-quality formatting for Nature-style figures
%    - Shaded uncertainty ranges for ensemble variability
%    - Standardized panel labels and annotations
%
% Overall, this figure assesses the robustness of the temporal relationship
% between AMOC variability and coastal dynamic sea-level changes in CMIP6
% models and observations, providing insight into the capability of climate
% models to reproduce observed North Atlantic sea-level fingerprints of
% AMOC variability.
%% ========================================================================
clear; clc;

%% ============================================================
% Load CMIP6 results
%% ============================================================

load('AMOC_ZOS_CMIP6_Results.mat');   % loads struct "Results"

% Choose analysis period
% 1 = 1850–2014
% 2 = 1993–2013
p = 1;

R = Results(p);

%% ============================================================
% Align models to common length
%% ============================================================

minLength = min(arrayfun(@(x) numel(x.AMOC), R.Data));
years = R.period(1):(R.period(1)+minLength-1);
N = numel(R.Data);

amoc_all = arrayfun(@(x) x.AMOC(1:minLength), R.Data, 'uni', 0);
zos_all  = arrayfun(@(x) x.coastal_index(1:minLength), R.Data, 'uni', 0);

amoc_mat = cell2mat(amoc_all');
zos_mat  = cell2mat(zos_all');

%% ============================================================
% Normalize each model
%% ============================================================

for m = 1:N

    amoc_mat(m,:) = (amoc_mat(m,:) - mean(amoc_mat(m,:),'omitnan')) ...
                     / std(amoc_mat(m,:),'omitnan');

    zos_mat(m,:)  = (zos_mat(m,:) - mean(zos_mat(m,:),'omitnan')) ...
                     / std(zos_mat(m,:),'omitnan');

end

%% ============================================================
% Ensemble mean
%% ============================================================

amoc_ens = mean(amoc_mat,1,'omitnan');
zos_ens  = mean(zos_mat,1,'omitnan');
r_ens = corr(amoc_ens', zos_ens','rows','complete');

amoc_std_ens = std(amoc_mat,0,1,'omitnan');
zos_std_ens  = std(zos_mat,0,1,'omitnan');
%% ============================================================
% Load Observations
%% ============================================================

load('obs_AMOC_ZOS_all.mat');  % loads amoc_s, zos_s, yrs_trimmed

% Restrict observations to 1993–2014
yr_idx = yrs_trimmed >= 1993 & yrs_trimmed <= 2014;

obs_AMOC = amoc_s(yr_idx);
obs_ZOS  = zos_s(yr_idx);
obs_years = yrs_trimmed(yr_idx);

%% ============================================================
% Clean + align observations (CMIP-style block)
%% ============================================================

X = [obs_AMOC(:) obs_ZOS(:) obs_years(:)];

% remove rows where either variable is NaN
X = X(all(~isnan(X(:,1:2)),2),:);

obs_AMOC_valid  = X(:,1);
obs_ZOS_valid   = X(:,2);
obs_years_valid = X(:,3);

% observed correlation
r_obs = corr(obs_AMOC_valid, obs_ZOS_valid);

% % Indices of overlapping years in CMIP6
% overlap_idx = years >= 1993 & years <= 2014;
% 
% % CMIP6 ensemble for overlapping years
% amoc_ens_overlap = amoc_ens(overlap_idx);
% zos_ens_overlap  = zos_ens(overlap_idx);
% 
% % Observations already restricted
% obs_AMOC_overlap = obs_AMOC_valid;
% obs_ZOS_overlap  = obs_ZOS_valid;
% 
% % Correlations: Obs vs CMIP6 ensemble
% r_obs_amoc = corr(obs_AMOC_overlap, amoc_ens_overlap', 'rows','complete');
% r_obs_zos  = corr(obs_ZOS_overlap,  zos_ens_overlap',  'rows','complete');
%% ============================================================
% Colors
%% ============================================================

% CMIP6
col_amoc_cmip = [0.20 0.40 0.80];
col_zos_cmip  = [0.85 0.33 0.10];

col_amoc_obs = [0 0 0];       % black
col_zos_obs  = [0.55 0 0.75]; % deep purple

%% ============================================================
% Figure layout
%% ============================================================

ncol = 4;
nrow = ceil((N+1)/ncol);

figure('Units','inches','Position',[1 1 10 6],'Color','w');

tiledlayout(nrow,ncol,'TileSpacing','compact','Padding','compact')

panel_labels = arrayfun(@(x) ['(' char('a'+x-1) ')'], ...
                        1:(N+1),'UniformOutput',false);

%% ============================================================
% Ensemble panel (updated with Obs vs CMIP6 correlations)
%% ============================================================

nexttile; hold on; box on

% Panel label
text(0.01,0.98,panel_labels{1}, ...
    'Units','normalized','FontSize',11, ...
    'FontWeight','bold','VerticalAlignment','top');

% Fill CMIP6 ensemble std shading
fill([years fliplr(years)], ...
     [amoc_ens-amoc_std_ens fliplr(amoc_ens+amoc_std_ens)], ...
     col_amoc_cmip, 'FaceAlpha',0.15,'EdgeColor','none');

fill([years fliplr(years)], ...
     [zos_ens-zos_std_ens fliplr(zos_ens+zos_std_ens)], ...
     col_zos_cmip, 'FaceAlpha',0.15,'EdgeColor','none');

% Vertical line for observation period start
xline(1993,'k--','LineWidth',0.8,'HandleVisibility','off');

% Plot CMIP6 ensemble mean
plot(years, amoc_ens,'-','Color',col_amoc_cmip,'LineWidth',2);
plot(years, zos_ens,'-','Color',col_zos_cmip,'LineWidth',2);

% --- OBS AMOC halo ---
plot(obs_years_valid,obs_AMOC_valid,'-','Color','w','LineWidth',5);
% --- OBS AMOC line ---
h_obs_amoc = plot(obs_years_valid,obs_AMOC_valid,'-', ...
     'Color',col_amoc_obs,'LineWidth',2);

% --- OBS ZOS halo ---
plot(obs_years_valid,obs_ZOS_valid,'-','Color','w','LineWidth',5);
% --- OBS ZOS line ---
h_obs_zos = plot(obs_years_valid,obs_ZOS_valid,'-', ...
     'Color',col_zos_obs,'LineWidth',2);

% Legend
legend([h_obs_amoc h_obs_zos], ...
       {'Obs AMOC','Obs Zoscoast'}, ...
       'Location','southwest', ...
       'FontSize',8, ...
       'Box','off');

%% =========================
% Compute correlations
%% =========================

% Overlap period indices for CMIP6
overlap_idx = years >= 1993 & years <= 2014;

amoc_ens_overlap = amoc_ens(overlap_idx);
zos_ens_overlap  = zos_ens(overlap_idx);

obs_AMOC_overlap = obs_AMOC_valid;
obs_ZOS_overlap  = obs_ZOS_valid;

% Correlations
r_ens      = corr(amoc_ens', zos_ens','rows','complete');          % CMIP6 AMOC vs ZOS
r_obs_amoc = corr(obs_AMOC_overlap, amoc_ens_overlap','rows','complete');  % Obs AMOC vs CMIP6 AMOC
r_obs_zos  = corr(obs_ZOS_overlap,  zos_ens_overlap','rows','complete');   % Obs ZOS vs CMIP6 ZOS

%Display
text(0.98, 0.95, ...
     sprintf('r_{ens} = %.2f, r_{AMOC(Obs vs CMIP6)} = %.2f', ...
             r_ens, r_obs_amoc), ...
     'Units','normalized', ...
     'FontSize',10, ...
     'FontWeight','bold', ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','top', ...
     'Interpreter','tex');
% 
% % Second line: r_ZOS directly under r_ens
text(0.98, 0.91, ...  % slightly lower y-position
     sprintf('r_{ZOS(Obs vs CMIP6)} = %.2f', r_obs_zos), ...
     'Units','normalized', ...
     'FontSize',10, ...
     'FontWeight','bold', ...
     'HorizontalAlignment','right', ...  % align left 
     'VerticalAlignment','top', ...
     'Interpreter','tex');
%% =========================
% Axes formatting
%% =========================

xlim([1850 2014])
xticks(1850:50:2014)

ylim([-5 5])
yticks(-5:2:5)

title('Ensemble mean','FontSize',10,'FontWeight','bold')
ylabel('Normalized anomaly','FontWeight','bold')

set(gca,'FontSize',10,'TickDir','out','LineWidth',1.1)
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off')
xticklabels([]);

%% ============================================================
% Individual model panels
%% ============================================================

for m = 1:N

    nexttile; hold on; box on

    text(0.01,0.98,panel_labels{m+1}, ...
        'Units','normalized','FontSize',10, ...
        'FontWeight','bold','VerticalAlignment','top');

    plot(years,amoc_mat(m,:),'-','Color',col_amoc_cmip,'LineWidth',1.2)
    plot(years,zos_mat(m,:),'-','Color',col_zos_cmip,'LineWidth',1.2)

    r = corr(amoc_mat(m,:)',zos_mat(m,:)','rows','complete');

    text(0.97,0.92,sprintf('r = %.2f',r), ...
        'Units','normalized','FontSize',9, ...
        'FontWeight','bold','HorizontalAlignment','right');

    title(R.Data(m).model_name,'FontSize',9)

    xlim([1850 2014])
    xticks(1850:50:2014)

    ylim([-5 5])
    yticks(-5:2:5)

    set(gca,'FontSize',9,'TickDir','out','LineWidth',0.8)

    yline(0,'k:','LineWidth',0.6,'HandleVisibility','off')

    if m == N
        xlabel('Year','FontWeight','bold')
    else
        xticklabels([])
    end

    yticklabels([])

end

%% ============================================================
% Legend
%% ============================================================

legend({'CMIP6 AMOC','CMIP6 Zoscoast'}, ...
        'Location','southoutside', ...
        'Orientation','horizontal', ...
        'Box','off')
