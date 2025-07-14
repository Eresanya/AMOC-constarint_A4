clear;clc;
%%
load('amoc_zos_indices.mat');%load computed obs. zosindices
load('AMOC_ZOS_GroupedData.mat');%load AMOC_ZOS_GroupedData
load('Ensemble_AMOC_ZOS_GroupedData.mat');%load cEnsemble_AMOC_ZOS_GroupedData
%
% CMIP6 Model Files and Names
model_names_zos = { ...  
    'zos_Omon_CESM2-WACCM-FV2_historical_r1i1p1f1_gn_1850-2014.nc',...
    'zos_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_1850-2014.nc',...
    'zos_Omon_CESM2-FV2_historical_r1i1p1f1_gn_1850-2014.nc',... 
    'zos_Omon_CESM2-WACCM_historical_r1i1p1f1_gr_1850-2014.nc'};

model_names_amoc = {... 
    'Atlantic_trans_None_CESM2-WACCM-FV2_historical_r1i1p1f1_gn_230122.nc',...
    'Atlantic_trans_None_IPSL-CM6A-LR_historical_r1i1p1f1_gn_230122.nc',...
    'Atlantic_trans_None_CESM2-FV2_historical_r1i1p1f1_gn_230122.nc',...
    'Atlantic_trans_None_CESM2-WACCM_historical_r1i1p1f1_gn_230122.nc'};

% Model Configuration
models = {'CESM2-WACCM-FV2','IPSL-CM6A-LR','CESM2-FV2','CESM2-WACCM',};
grid_options = {'gn', 'gr', 'gr1'};
start_year = 1850; end_year = 2014;
yrs = start_year:end_year;
time_full = datetime(start_year, 1, 1) + calyears(0:numel(yrs)-1);
target_years = 1993:2013;
t_idx = ismember(year(time_full), target_years);

% Initialize storage
n_models = numel(models);
n_years = numel(yrs);
zos_index_all = NaN(n_models, n_years);
zos_coast_index_all = NaN(n_models, n_years);
amoc_all = NaN(n_models, n_years);

% Main Loop: Process Each Model
for m = 1:n_models
    model = models{m};

    % --- Load ZOS Data ---
    zos_file = dir(['zos_Omon_' model '_historical_r1i1p1f1_*.nc']);
    if isempty(zos_file)
        warning("Missing ZOS file for %s", model);
        continue;
    end

    fpath = zos_file(1).name;
    zos = squeeze(ncread(fpath, 'zos'));
    lat = ncread(fpath, 'lat');
    lon = ncread(fpath, 'lon');
    time = ncread(fpath, 'time');

    % Convert time
    try
        time = datetime(1850,1,1) + days(time);
    catch
        time = datetime(1850,1,1) + calmonths(0:numel(time)-1);
    end

    % Monthly to Annual
    if mod(size(zos, 3), 12) ~= 0
        warning("ZOS monthly data not divisible by 12: %s", model);
        continue;
    end
    zos = zos - mean(zos, [1 2], 'omitnan');

    % Indices
    zos1 = squeeze(mean(mean(zos(lon >= 270 & lon <= 360, lat >= 0 & lat <= 40, :), 1, 'omitnan'), 2, 'omitnan'));
    zos2 = squeeze(mean(mean(zos(lon >= 270 & lon <= 360, lat >= 40 & lat <= 60, :), 1, 'omitnan'), 2, 'omitnan'));
    zos_index = zos1 - zos2;

    coast1 = squeeze(mean(mean(zos(lon >= 290 & lon <= 320, lat >= 20 & lat <= 40, :), 1, 'omitnan'), 2, 'omitnan'));
    coast2 = squeeze(mean(mean(zos(lon >= 300 & lon <= 330, lat >= 40 & lat <= 60, :), 1, 'omitnan'), 2, 'omitnan'));
    zos_coast_index = coast1 - coast2;

    zos_index_ann = mean(reshape(zos_index, 12, []), 1, 'omitnan');
    zos_coast_ann = mean(reshape(zos_coast_index, 12, []), 1, 'omitnan');

    if numel(zos_index_ann) ~= n_years
        warning("ZOS annual length mismatch: %s", model);
        continue;
    end

    % Detrend
    zos_index_all(m, :) = zos_index_ann - mean(zos_index_ann, 'omitnan');
    zos_coast_index_all(m, :) = zos_coast_ann - mean(zos_coast_ann, 'omitnan');

    % --- Load AMOC Data ---
    amoc_found = false;
    for g = 1:numel(grid_options)
        amoc_file = sprintf('Atlantic_trans_None_%s_historical_r1i1p1f1_%s_230122.nc', model, grid_options{g});
        if exist(amoc_file, 'file')
            amoc_found = true;
            break;
        end
    end
    if ~amoc_found
        warning("Missing AMOC file for %s", model);
        continue;
    end

    try
        vlat = ncread(amoc_file, 'vlat');
        moc = ncread(amoc_file, 'moc_section');

        moc_ann = mean(reshape(moc, size(moc,1), size(moc,2), 12, []), 3, 'omitnan');
        lat_idx = find(vlat >= 26 & vlat <= 26.5);
        if isempty(lat_idx)
            warning("No valid latitude for AMOC: %s", model);
            continue;
        end

        amoc_lat = max(moc_ann(lat_idx, :, :), [], 1);
        amoc_depth = squeeze(max(amoc_lat, [], 2));
        amoc_anom = amoc_depth - mean(amoc_depth, 'omitnan');
        amoc_all(m, :) = smoothdata(amoc_anom, 'movmean', 5, 'omitnan');
    catch
        warning("Error reading AMOC: %s", model);
    end
end

% Compute Ensemble Mean and Std (1993–2013)
time_sub = time_full(t_idx);
amoc_ens = mean(amoc_all(:, t_idx), 1, 'omitnan');
amoc_std = std(amoc_all(:, t_idx), 0, 1, 'omitnan');

zos_index_ens = mean(zos_index_all(:, t_idx), 1, 'omitnan');
zos_index_std = std(zos_index_all(:, t_idx), 0, 1, 'omitnan');

zos_coast_index_ens = mean(zos_coast_index_all(:, t_idx), 1, 'omitnan');
zos_coast_index_std = std(zos_coast_index_all(:, t_idx), 0, 1, 'omitnan');

% Create publication-quality figure: AMOC vs ZOS Coast Index
figure('Units', 'inches', 'Position', [1 1 6.5 4], 'PaperPositionMode', 'auto'); 
hold on;

% Define colors (consistent with GRL aesthetic: distinct but not too bright)
amoc_color = [0 0.45 0.74];       % Blue
zos_coast_color = [0.3 0.6 0.2];  % Muted green

% Plot AMOC (left axis)
yyaxis left
fill([time_sub fliplr(time_sub)], ...
     [amoc_ens + amoc_std, fliplr(amoc_ens - amoc_std)], ...
     amoc_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(time_sub, amoc_ens, '-', 'LineWidth', 1.8, 'Color', amoc_color);
plot(time_sub, ref_ssh_trimmed * 2, '--', 'LineWidth', 1.5, 'Color', [0 0 0]); % Obs: Black dashed
plot(time_sub, ens_amoc_unconstraint, ':', 'LineWidth', 1.5, 'Color', amoc_color);
ylabel('AMOC max (Sv)', 'Color', amoc_color);
ylim([-3 3]);
set(gca, 'YColor', amoc_color);

% Plot ZOS Coast (right axis)
yyaxis right
fill([time_sub fliplr(time_sub)], ...
     [zos_coast_index_ens + zos_coast_index_std, ...
      fliplr(zos_coast_index_ens - zos_coast_index_std)], ...
     zos_coast_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(time_sub, zos_coast_index_ens, '-', 'LineWidth', 1.8, 'Color', zos_coast_color);
plot(time_sub, zos_Obs_coast_index, '--', 'LineWidth', 1.5, 'Color', [0.2 0.5 0.1]); % Obs: darker green
plot(time_sub, ens_zos_coast_unconstraint, ':', 'LineWidth', 1.5, 'Color', zos_coast_color);
ylabel('ZOS coast index (m)', 'Color', zos_coast_color);
ylim([-0.1 0.1]);
set(gca, 'YColor', zos_coast_color);

% Final styling
xlabel('');
title('AMOC and ZOS Coast Index Ensemble Mean ±1 Std (1993–2013)', ...
    'FontWeight', 'normal', 'FontSize', 12);

grid on;
box on;
set(gca, 'FontSize', 10);
set(gca, 'TickDir', 'out');
set(gca, 'LineWidth', 1);
xticks(datetime(1993:2:2013, 1, 1));

% Create the legend entries
lgd = legend({'AMOC ±1 std', 'AMOC constraint', 'AMOC (observed)', 'AMOC (unconstrained)', ...
              'ZOS Coast ±1 std', 'ZOS Coast constraint', 'ZOS Coast (observed)', 'ZOS Coast (unconstrained)'}, ...
              'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9);

% Force the legend into two rows (4 items per row)
lgd.NumColumns = 4;

% Position the legend below the figure
lgd.Position = [0.1, 0.02, 0.8, 0.05];  % [left, bottom, width, height]


% Tight layout for publication
set(gcf, 'Color', 'w');