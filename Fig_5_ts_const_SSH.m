%% ============================================================
% Dual-panel Figure: TG–ZOS & TG–AMOC (1940–2014)
% Data: CMIP6 models + Tide Gauge
% ============================================================
clear; clc;

%% === LOAD DATA =============================================
load('AMOC_ZOS_CMIP6_Results.mat');  % CMIP6 model results (1850–2014)
load('TG_adjstd_SSH.mat') %Both SSH and TG detrended
% load('TG_adjstd_SSH_updated.mat');% TG detrended

R = Results(1);   % CMIP6 results (1850–2014)
D = R.Data;
nmod = length(D);

% Time vector
time = R.period(1):R.period(2);
ntime = length(time);

% Extract model data
AMOC = nan(ntime,nmod);
ZOS  = nan(ntime,nmod);
model_names = cell(1,nmod);

for m = 1:nmod
    AMOC(:,m) = D(m).AMOC;
    ZOS(:,m)  = D(m).coastal_index;
    model_names{m} = D(m).model_name;
end

%% === STANDARDIZE CMIP6 SERIES =============================
AMOC_std = (AMOC - mean(AMOC,1)) ./ std(AMOC,[],1);
ZOS_std  = (ZOS  - mean(ZOS,1))  ./ std(ZOS,[],1);

%% === MODEL GROUPS ==========================================

unconstraint_models= {'E3SM-1-1',...
    'INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR', ...
    'MRI-ESM2-0','NorESM2-LM',...
    };

 constraint_models= { ...
     'CanESM5','CESM2-WACCM','CNRM-CM6-1-HR','FIO-ESM-2-0','HadGEM3-GC31-LL',...
    'BCC-CSM2-MR','CAMS-CSM1-0','CESM2','CNRM-CM6-1','HadGEM3-GC31-MM', ...
    'MPI-ESM1-2-HR','SAM0-UNICON','NESM3','UKESM1-0-LL'}; 

idx_con = find(ismember(model_names,constraint_models));
idx_un  = find(ismember(model_names,unconstraint_models));
idx_all = 1:nmod;

% === GROUP SPREAD (PHYSICAL SPACE) =======================
AMOC_con_spread_full = std(AMOC_std(:,idx_con), 0, 2, 'omitnan');
AMOC_un_spread_full  = std(AMOC_std(:,idx_un),  0, 2, 'omitnan');

ZOS_con_spread_full = std(ZOS_std(:,idx_con), 0, 2, 'omitnan');
ZOS_un_spread_full  = std(ZOS_std(:,idx_un),  0, 2, 'omitnan');


% Ensemble mean function
ensTS = @(X,idx) mean(X(:,idx),2,'omitnan');

AMOC_con_std = ensTS(AMOC_std,idx_con); %ensemble mean of standardized models
AMOC_un_std  = ensTS(AMOC_std,idx_un);
ZOS_con_std  = ensTS(ZOS_std,idx_con);
ZOS_un_std   = ensTS(ZOS_std,idx_un);

%% === TARGET YEARS & TG SERIES =============================
target_years = 1940:2014;

% Time masks
model_mask = ismember(time,target_years);
tg_years   = 1880:2025;   % adjust if needed
tg_mask    = ismember(tg_years,target_years);

% --- Extract spread (same mask as means)
AMOC_con_spread = AMOC_con_spread_full(model_mask);
AMOC_un_spread  = AMOC_un_spread_full(model_mask);

ZOS_con_spread = ZOS_con_spread_full(model_mask);
ZOS_un_spread  = ZOS_un_spread_full(model_mask);

% ---- Extract CMIP6 (already standardized earlier) ----
ZOS_con = ZOS_con_std(model_mask);
ZOS_un  = ZOS_un_std(model_mask);
AMOC_con = AMOC_con_std(model_mask);
AMOC_un  = AMOC_un_std(model_mask);

% ---- Extract TG ---- 
tg_zos_raw = tg_grad_adj_full; %TG index (TG zoscoast)
% tg_zos_raw = movmean(tg_grad_adj_full,7,'omitnan'); %TG index (TG zoscoast)
tg_zos_raw = tg_zos_raw(tg_mask); %TG index 1940-2014
% tg_zos_raw = tg_zos_raw/100;
 
% tg_amoc_raw = amoc_obs_from_tg_full;
tg_amoc_raw = tg_grad_adj_full; 
% tg_amoc_raw = movmean(amoc_obs_from_tg_full,7,'omitnan');  
tg_amoc_raw = tg_amoc_raw(tg_mask); %TG index 1940-2014
% tg_amoc_raw = tg_amoc_raw/100;
time_common = datetime(target_years,1,1);

%% === STANDARDIZE TG ========================
tg_zos_std = (tg_zos_raw - mean(tg_zos_raw,'omitnan')) ./ ...
              std(tg_zos_raw,'omitnan');

tg_amoc_std = (tg_amoc_raw - mean(tg_amoc_raw,'omitnan')) ./ ...
               std(tg_amoc_raw,'omitnan');
%% === CORRELATIONS ==================
r1_con = corr(tg_zos_std, ZOS_con,'rows','complete');
r1_un  = corr(tg_zos_std, ZOS_un,'rows','complete');

r2_con = corr(tg_amoc_std, AMOC_con,'rows','complete');
r2_un  = corr(tg_amoc_std, AMOC_un,'rows','complete');

%% === LEAD-LAG ANALYSIS: AMOC LEADING =======================

maxLag = 10;                     % test 0–10 year lead
lags = 0:maxLag;

r_lag_con = nan(length(lags),1);
r_lag_un  = nan(length(lags),1);

for i = 1:length(lags)
    L = lags(i);
    
    % AMOC leading TG-AMOC estimate
    r_lag_con(i) = corr(AMOC_con(1:end-L), ...
                        tg_amoc_std(1+L:end), ...
                        'rows','complete');
                    
    r_lag_un(i)  = corr(AMOC_un(1:end-L), ...
                        tg_amoc_std(1+L:end), ...
                        'rows','complete');
end

% Find maximum correlation and corresponding lag
[max_r_con, idx_max_con] = max(r_lag_con);
best_lag_con = lags(idx_max_con);

[max_r_un, idx_max_un] = max(r_lag_un);
best_lag_un = lags(idx_max_un);

%% === 5-YEAR RUNNING MEAN ==================================
tg_smooth = movmean(tg_zos_std,5,'omitnan');
zos_con_smooth = movmean(ZOS_con,5,'omitnan');
zos_un_smooth  = movmean(ZOS_un,5,'omitnan');

amoc_tg_smooth  = movmean(tg_amoc_std,5,'omitnan');
amoc_con_smooth = movmean(AMOC_con,5,'omitnan');
amoc_un_smooth  = movmean(AMOC_un,5,'omitnan');

% --- SMOOTH SPREAD ---------------------------------------
AMOC_con_spread_s = movmean(AMOC_con_spread,5,'omitnan');
AMOC_un_spread_s  = movmean(AMOC_un_spread,5,'omitnan');

ZOS_con_spread_s = movmean(ZOS_con_spread,5,'omitnan');
ZOS_un_spread_s  = movmean(ZOS_un_spread,5,'omitnan');
%% === SAVE TG AMOC SERIES ==================================
save('TG_AMOC_processed_1940_2014.mat', ...
     'tg_amoc_std', ...
     'amoc_tg_smooth', ...
     'time_common');
%% === COLOR PALETTE ========================================
col_tg      = [0 0 0];        % black TG
col_con     = [0.90 0.60 0.00]; % orange constrained
col_un      = [0.00 0.60 0.50]; % teal unconstrained ZOS
col_un_dark = [0.35 0.70 0.90]; % sky blue unconstrained AMOC
col_amoc    = [0.35 0.35 0.35]; % dark gray TG AMOC

%% ===============================================
% TG–AMOC: Trend Analysis & Nature-standard Visualization
% Periods: 1950–1990 and 1990–2014
% ===============================================
% clearvars -except time_common tg_amoc_std amoc_tg_smooth AMOC_con AMOC_un amoc_con_smooth amoc_un_smooth col_amoc col_con col_un_dark
clearvars -except time_common tg_amoc_std amoc_tg_smooth ...
                     AMOC_con AMOC_un amoc_con_smooth amoc_un_smooth ...
                     AMOC_con_spread_s AMOC_un_spread_s ...
                     col_amoc col_con col_un_dark

% --- Convert datetime to numeric year for plotting
time_num = year(time_common);

% --- Define periods --------------------------------
period1_mask = (time_num >= 1950) & (time_num <= 1990);
period2_mask = (time_num > 1990);

t_p1 = time_num(period1_mask);
t_p2 = time_num(period2_mask);

tg_p1 = amoc_tg_smooth(period1_mask);
tg_p2 = amoc_tg_smooth(period2_mask);

con_p1 = amoc_con_smooth(period1_mask);
con_p2 = amoc_con_smooth(period2_mask);

un_p1 = amoc_un_smooth(period1_mask);
un_p2 = amoc_un_smooth(period2_mask);

% --- Linear trends (slope per year) ----------------
tg_trend_p1  = polyfit(t_p1, tg_p1, 1);
tg_trend_p2  = polyfit(t_p2, tg_p2, 1);

con_trend_p1 = polyfit(t_p1, con_p1, 1);
con_trend_p2 = polyfit(t_p2, con_p2, 1);

un_trend_p1  = polyfit(t_p1, un_p1, 1);
un_trend_p2  = polyfit(t_p2, un_p2, 1);

% --- Trend differences & correlations -------------
delta_con_p1 = tg_trend_p1(1) - con_trend_p1(1);
delta_con_p2 = tg_trend_p2(1) - con_trend_p2(1);

delta_un_p1 = tg_trend_p1(1) - un_trend_p1(1);
delta_un_p2 = tg_trend_p2(1) - un_trend_p2(1);

r_con_p1 = corr(tg_p1, con_p1,'rows','complete');
r_con_p2 = corr(tg_p2, con_p2,'rows','complete');

r_un_p1  = corr(tg_p1, un_p1,'rows','complete');
r_un_p2  = corr(tg_p2, un_p2,'rows','complete');

% % Agreement metric
% agreement_con_p1 = abs(r_con_p1)/(1+abs(delta_con_p1));
% agreement_con_p2 = abs(r_con_p2)/(1+abs(delta_con_p2));
% agreement_un_p1  = abs(r_un_p1)/(1+abs(delta_un_p1));
% agreement_un_p2  = abs(r_un_p2)/(1+abs(delta_un_p2));
%% =========================================================
%% TAYLOR SKILL SCORE
%% =========================================================

% --- Standard deviation ratios (model / observation)
sigma_con_p1 = std(con_p1,'omitnan') / std(tg_p1,'omitnan');
sigma_con_p2 = std(con_p2,'omitnan') / std(tg_p2,'omitnan');

sigma_un_p1  = std(un_p1,'omitnan') / std(tg_p1,'omitnan');
sigma_un_p2  = std(un_p2,'omitnan') / std(tg_p2,'omitnan');

% --- Taylor Skill Score
taylorSkill = @(r,sigma) (2 * (1 + r)) ./ ((sigma + 1./sigma).^2);

skill_con_p1 = taylorSkill(r_con_p1, sigma_con_p1);
skill_con_p2 = taylorSkill(r_con_p2, sigma_con_p2);

skill_un_p1  = taylorSkill(r_un_p1, sigma_un_p1);
skill_un_p2  = taylorSkill(r_un_p2, sigma_un_p2);

% --- Summary table --------------------------------
% Metric        = {'Trend (slope/yr)'; 'Trend diff'; 'Correlation'; 'Agreement'};
% Constrained_p1 = [con_trend_p1(1); delta_con_p1; r_con_p1; agreement_con_p1];
% Constrained_p2 = [con_trend_p2(1); delta_con_p2; r_con_p2; agreement_con_p2];
% Unconstrained_p1 = [un_trend_p1(1); delta_un_p1; r_un_p1; agreement_un_p1];
% Unconstrained_p2 = [un_trend_p2(1); delta_un_p2; r_un_p2; agreement_un_p2];
Metric        = {'Trend (slope/yr)'; 'Trend diff'; 'Correlation'; 'Taylor Skill'};
Constrained_p1 = [con_trend_p1(1); delta_con_p1; r_con_p1; skill_con_p1];
Constrained_p2 = [con_trend_p2(1); delta_con_p2; r_con_p2; skill_con_p2];
Unconstrained_p1 = [un_trend_p1(1); delta_un_p1; r_un_p1; skill_un_p1];
Unconstrained_p2 = [un_trend_p2(1); delta_un_p2; r_un_p2; skill_un_p2];


T = table(Metric, Constrained_p1, Constrained_p2, Unconstrained_p1, Unconstrained_p2);
disp('TG–AMOC Agreement Summary (1950–1990 vs 1990–2014):');
disp(T);

% --- Enhanced figure (Nature style) ----------------
figure('Color','w','Position',[100 100 1000 550]); hold on;

% --- Determine y-limits for shading
yl = [min([amoc_tg_smooth; amoc_con_smooth; amoc_un_smooth])-0.2, ...
      max([amoc_tg_smooth; amoc_con_smooth; amoc_un_smooth])+0.2];

% --- Highlight periods with light shading
patch([1950 1990 1990 1950], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'FaceAlpha',0.25,'EdgeColor','none');
patch([1990 2014 2014 1990], [yl(1) yl(1) yl(2) yl(2)], [0.95 0.95 0.95], 'FaceAlpha',0.25,'EdgeColor','none');

t = time_num(:);  % force column vector

% --- Constrained envelope ---
hConSpread = fill([t; flipud(t)], ...
     [amoc_con_smooth(:) - AMOC_con_spread_s(:); ...
      flipud(amoc_con_smooth(:) + AMOC_con_spread_s(:))], ...
     col_con, 'FaceAlpha',0.12,'EdgeColor','none');

% --- Unconstrained envelope ---
hUnSpread = fill([t; flipud(t)], ...
     [amoc_un_smooth(:) - AMOC_un_spread_s(:); ...
      flipud(amoc_un_smooth(:) + AMOC_un_spread_s(:))], ...
     col_un_dark, 'FaceAlpha',0.12,'EdgeColor','none');

hConSpread.Annotation.LegendInformation.IconDisplayStyle = 'on';
hUnSpread.Annotation.LegendInformation.IconDisplayStyle  = 'on';

% --- Plot smoothed series on top ---------------------------
hTG  = plot(time_num, amoc_tg_smooth, '-', 'Color', col_amoc,'LineWidth',1.5);
hCon = plot(time_num, amoc_con_smooth, '--', 'Color', col_con,'LineWidth',1.5);
hUn  = plot(time_num, amoc_un_smooth, '-.', 'Color', col_un_dark,'LineWidth',1.5);

% --- Trend lines
plot(t_p1, polyval(tg_trend_p1, t_p1), '-', 'Color', col_amoc, 'LineWidth',1);
plot(t_p1, polyval(con_trend_p1, t_p1), '--', 'Color', col_con, 'LineWidth',1);
plot(t_p1, polyval(un_trend_p1, t_p1), '-.', 'Color', col_un_dark, 'LineWidth',1);

plot(t_p2, polyval(tg_trend_p2, t_p2), '-', 'Color', col_amoc, 'LineWidth',1, 'LineStyle', ':');
plot(t_p2, polyval(con_trend_p2, t_p2), '--', 'Color', col_con, 'LineWidth',1, 'LineStyle', ':');
plot(t_p2, polyval(un_trend_p2, t_p2), '-.', 'Color', col_un_dark, 'LineWidth',1, 'LineStyle', ':');

yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1);
xline(1950,'--','Color',[0.5 0.5 0.5],'LineWidth',1);
% --- Axis, labels, title
ylabel('Standardized AMOC','FontSize',14,'FontWeight','bold','FontName','Times New Roman');
xlabel('Year','FontSize',14,'FontWeight','bold','FontName','Times New Roman');
title('TG–CMIP6 AMOC Relationship (1940-2014)','FontSize',18,'FontWeight','bold','FontName','Times New Roman');

legend([hTG, hCon, hConSpread, hUn, hUnSpread], ...
       {'TG AMOC', ...
        'Constrained AMOC', 'Constrained Spread', ...
        'Unconstrained AMOC', 'Unconstrained Spread'}, ...
       'Location','northwest','FontSize',11,'Box','off');
% --- Grid and ticks
grid on; box on;
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.XLim = [1940 2014];
ax.YLim = [-2.5 2.5];
ax.YTick = -2.5:1:2.5;
% --- Annotate period 1 (1950–1990) trend values on the left
x_left = t_p1(1);  % first year of period 1
text(x_left, polyval(tg_trend_p1, x_left), ...
    sprintf('TG: %.2f/yr', tg_trend_p1(1)), ...
    'FontSize',10,'Color',col_amoc,'HorizontalAlignment','left','VerticalAlignment','bottom');

text(x_left, polyval(con_trend_p1, x_left), ...
    sprintf('Constr: %.2f/yr', con_trend_p1(1)), ...
    'FontSize',10,'Color',col_con,'HorizontalAlignment','left','VerticalAlignment','bottom');

text(x_left, polyval(un_trend_p1, x_left), ...
    sprintf('Unconstr: %.2f/yr', un_trend_p1(1)), ...
    'FontSize',10,'Color',col_un_dark,'HorizontalAlignment','left','VerticalAlignment','bottom');

% --- Annotate period 2 (1990–2014) trend values on the right
x_right = t_p2(end);  % last year of period 2
text(x_right, polyval(tg_trend_p2, x_right), ...
    sprintf('TG: %.2f/yr', tg_trend_p2(1)), ...
    'FontSize',10,'Color',col_amoc,'HorizontalAlignment','right','VerticalAlignment','bottom');

text(x_right, polyval(con_trend_p2, x_right), ...
    sprintf('Constr: %.2f/yr', con_trend_p2(1)), ...
    'FontSize',10,'Color',col_con,'HorizontalAlignment','right','VerticalAlignment','bottom');

text(x_right, polyval(un_trend_p2, x_right), ...
    sprintf('Unconstr: %.2f/yr', un_trend_p2(1)), ...
    'FontSize',10,'Color',col_un_dark,'HorizontalAlignment','right','VerticalAlignment','bottom');

% --- Compact metrics summary string
% summaryStr = sprintf([ ...
%     'Metric        Constr   Unconstr\n' ...
%     '--------------------------------\n' ...
%     'Trend P1      %.2f     %.2f\n' ...
%     'Trend P2      %.2f     %.2f\n' ...
%     'Trend diff P1 %.2f     %.2f\n' ...
%     'Trend diff P2 %.2f     %.2f\n' ...
%     'Corr P1       %.2f     %.2f\n' ...
%     'Corr P2       %.2f     %.2f\n' ...
%     % 'Agree P1      %.2f     %.2f\n' ...
%     % 'Agree P2      %.2f     %.2f'], ...
%     'Skill P1      %.2f     %.2f\n'...
%     'Skill P2      %.2f     %.2f'], ...
%     con_trend_p1(1), un_trend_p1(1), ...
%     con_trend_p2(1), un_trend_p2(1), ...
%     delta_con_p1, delta_un_p1, ...
%     delta_con_p2, delta_un_p2, ...
%     r_con_p1, r_un_p1, ...
%     r_con_p2, r_un_p2, ...
%     agreement_con_p1, agreement_un_p1, ...
%     agreement_con_p2, agreement_un_p2);
summaryStr = sprintf([ ...
    'Metric        Constr   Unconstr\n' ...
    '--------------------------------\n' ...
    'Trend P1      %.2f     %.2f\n' ...
    'Trend P2      %.2f     %.2f\n' ...
    'Trend diff P1 %.2f     %.2f\n' ...
    'Trend diff P2 %.2f     %.2f\n' ...
    'Corr P1       %.2f     %.2f\n' ...
    'Corr P2       %.2f     %.2f\n' ...
    'Skill P1      %.2f     %.2f\n' ...
    'Skill P2      %.2f     %.2f'], ...
    con_trend_p1(1), un_trend_p1(1), ...
    con_trend_p2(1), un_trend_p2(1), ...
    delta_con_p1, delta_un_p1, ...
    delta_con_p2, delta_un_p2, ...
    r_con_p1, r_un_p1, ...
    r_con_p2, r_un_p2, ...
    skill_con_p1, skill_un_p1, ...
    skill_con_p2, skill_un_p2);

% --- Smaller box inside figure
text(0.97,0.93,summaryStr,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top',...
    'FontSize',9,...
    'FontName','Consolas',...
    'BackgroundColor',[1 1 1],...
    'EdgeColor',[0.7 0.7 0.7],...
    'Margin',4);

%% % ---------plot for TG and constrain models alone-----
figure('Color','w','Position',[100 100 1000 550]); hold on;

% --- Determine y-limits for shading
yl = [min([amoc_tg_smooth; amoc_con_smooth])-0.2, ...
      max([amoc_tg_smooth; amoc_con_smooth])+0.2];

% --- Highlight periods with light shading
patch([1950 1990 1990 1950], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'FaceAlpha',0.25,'EdgeColor','none');
patch([1990 2014 2014 1990], [yl(1) yl(1) yl(2) yl(2)], [0.95 0.95 0.95], 'FaceAlpha',0.25,'EdgeColor','none');

t = time_num(:);  % force column vector
% --- Constrained envelope ---
hConSpread = fill([t; flipud(t)], ...
     [amoc_con_smooth(:) - AMOC_con_spread_s(:); ...
      flipud(amoc_con_smooth(:) + AMOC_con_spread_s(:))], ...
     col_con, 'FaceAlpha',0.12,'EdgeColor','none');

hConSpread.Annotation.LegendInformation.IconDisplayStyle = 'on';
% --- Plot smoothed series on top ---------------------------
hTG  = plot(time_num, amoc_tg_smooth, '-', 'Color', col_amoc,'LineWidth',1.5);
hCon = plot(time_num, amoc_con_smooth, '--', 'Color', col_con,'LineWidth',1.5);

% --- Trend lines
plot(t_p1, polyval(tg_trend_p1, t_p1), '--', 'Color', col_amoc, 'LineWidth',1);
plot(t_p1, polyval(con_trend_p1, t_p1), '--', 'Color', col_con, 'LineWidth',1);

plot(t_p2, polyval(tg_trend_p2, t_p2), '-', 'Color', col_amoc, 'LineWidth',2, 'LineStyle', ':');
plot(t_p2, polyval(con_trend_p2, t_p2), '--', 'Color', col_con, 'LineWidth',2, 'LineStyle', ':');

yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1);
xline(1950,'--','Color',[0.5 0.5 0.5],'LineWidth',1);
% --- Axis, labels, title
ylabel('Standardized AMOC','FontSize',14,'FontWeight','bold','FontName','Times New Roman');
xlabel('Year','FontSize',14,'FontWeight','bold','FontName','Times New Roman');
title('TG–CMIP6 AMOC Relationship (1940-2014)','FontSize',18,'FontWeight','bold','FontName','Times New Roman');

legend([hTG, hCon, hConSpread], ...
       {'TG AMOC', ...
        'Constrained AMOC', 'Constrained Spread'}, ...
       'Location','northwest','FontSize',11,'Box','off');
% --- Grid and ticks
grid on; box on;
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.XLim = [1940 2014];
ax.YLim = [-2.5 2.5];
ax.YTick = -2.5:1:2.5;
% --- Annotate period 1 (1950–1990) trend values on the left
x_left = t_p1(1);  % first year of period 1
text(x_left, polyval(tg_trend_p1, x_left), ...
    sprintf('TG: %.2f/yr', tg_trend_p1(1)), ...
    'FontSize',10,'Color',col_amoc,'HorizontalAlignment','left','VerticalAlignment','bottom');

text(x_left, polyval(con_trend_p1, x_left), ...
    sprintf('Constr: %.2f/yr', con_trend_p1(1)), ...
    'FontSize',10,'Color',col_con,'HorizontalAlignment','left','VerticalAlignment','bottom');

% --- Annotate period 2 (1990–2014) trend values on the right
x_right = t_p2(end);  % last year of period 2
text(x_right, polyval(tg_trend_p2, x_right), ...
    sprintf('TG: %.2f/yr', tg_trend_p2(1)), ...
    'FontSize',10,'Color',col_amoc,'HorizontalAlignment','right','VerticalAlignment','bottom');

text(x_right, polyval(con_trend_p2, x_right), ...
    sprintf('Constr: %.2f/yr', con_trend_p2(1)), ...
    'FontSize',10,'Color',col_con,'HorizontalAlignment','right','VerticalAlignment','bottom');
