clear;clc;

%% CMIP6 Model Files and Names

model_names_zos = { ...
    'zos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_1850-2014.nc',... 
    'zos_Omon_BCC-CSM2-MR_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_CAMS-CSM1-0_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_CanESM5_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_CESM2_historical_r1i1p1f1_gr_1850-2014.nc',...  
    'zos_Omon_CESM2-FV2_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_CESM2-WACCM_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_CESM2-WACCM-FV2_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_CMCC-CM2-HR4_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_E3SM-1-0_historical_r1i1p1f1_gr_1850-2014.nc',... 
    'zos_Omon_E3SM-1-1_historical_r1i1p1f1_gr_1850-2014.nc',...  
    'zos_Omon_FGOALS-f3-L_historical_r1i1p1f1_gn_1850-2014.nc',...   
    'zos_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_INM-CM4-8_historical_r1i1p1f1_gr1_1850-2014.nc',...  
    'zos_Omon_INM-CM5-0_historical_r1i1p1f1_gr1_1850-2014.nc',... 
    'zos_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_MIROC6_historical_r1i1p1f1_gn_1850-2014.nc',... 
    'zos_Omon_NESM3_historical_r1i1p1f1_gn_1850-2014.nc',...  
    'zos_Omon_NorESM2-LM_historical_r1i1p1f1_gn_1850-2014.nc'};

model_names_amoc = {...
    'Atlantic_trans_None_ACCESS-CM2_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_BCC-CSM2-MR_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_CAMS-CSM1-0_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_CanESM5_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_CESM2_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_CESM2-FV2_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_CESM2-WACCM_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_CESM2-WACCM-FV2_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_CMCC-CM2-HR4_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_E3SM-1-0_historical_r1i1p1f1_gr_230122.nc',... 
    'Atlantic_trans_None_E3SM-1-1_historical_r1i1p1f1_gr_230122.nc',...  
    'Atlantic_trans_None_FGOALS-f3-L_historical_r1i1p1f1_gn_230122.nc',...   
    'Atlantic_trans_None_GISS-E2-1-G_historical_r1i1p1f1_gn_230122.nc',...  
    'Atlantic_trans_None_INM-CM4-8_historical_r1i1p1f1_gr1_230122.nc',...  
    'Atlantic_trans_None_INM-CM5-0_historical_r1i1p1f1_gr1_230122.nc',... 
    'Atlantic_trans_None_IPSL-CM6A-LR_historical_r1i1p1f1_gn_230122.nc',...  
    'Atlantic_trans_None_MIROC6_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_NESM3_historical_r1i1p1f1_gn_230122.nc',... 
    'Atlantic_trans_None_NorESM2-LM_historical_r1i1p1f1_gr_230122.nc'};


% CMIP6 Model Files and Names (FIO-ESM2 removed)
models_name = {...
    'ACCESS-CM2', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5', 'CESM2', ...
    'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV2', 'CMCC-CM2-HR4', ...
    'E3SM-1-0', 'E3SM-1-1', 'FGOALS-f3-L','GISS-E2-1-G', ...
    'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MIROC6', 'NESM3', 'NorESM2-LM'};

grid_options = {'gn', 'gr1', 'gr'};
yr = 1993:2013; % Updated year range to 1993-2013
lat_bands = {[25, 26]};
anom_models = cell(length(models_name), length(lat_bands));


% Define categories as cell arrays
AMOC_strong = {'CESM2', 'CESM2-WACCM', 'GISS-E2-1-G', 'NorESM2-LM'};
AMOC_weak = {'NESM3', 'MIROC6', 'IPSL-CM6A-LR', 'INM-CM5-0', 'BCC-CSM2-MR'};
Eqt_pos = {'ACCESS-CM2', 'CAMS-CSM1-0', 'CESM2', 'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV2', 'E3SM-1-0', 'E3SM-1-1', 'FGOALS-f3-L', 'INM-CM4-8', 'INM-CM5-0'};
Eqt_neg = {'CMCC-CM2-HR4', 'GISS-E2-1-G', 'MIROC6', 'NESM3'};
Eqt_bipol = {'BCC-CSM2-MR', 'CanESM5', 'IPSL-CM6A-LR', 'NorESM2-LM'};
wb_np={'ACCESS-CM2', 'BCC-CSM2-MR','CanESM5', 'CESM2','CESM2-FV2','CESM2-WACCM','CESM2-WACCM-FV2',...
  'CMCC-CM2-HR4','IPSL-CM6A-LR','MIROC6', 'NESM3', 'NorESM2-LM'};
wb_nn={'CAMS-CSM1-0','FGOALS-f3-L','GISS-E2-1-G','INM-CM4-8','INM-CM5-0'};
wb_pn={'E3SM-1-0', 'E3SM-1-1'};

% Ensure consistent model count
assert(length(model_names_zos) == length(model_names_amoc), 'Mismatch in the number of models between ZOS and AMOC.');

% Preallocate variables for storing results across models
zos_index_all = [];
zos_coast_index_all = [];
amoc_max_all = [];
time_annual_all = [];

for i = 1:length(model_names_zos)
    try
        ncfile_zos = model_names_zos{i};
        ncfile_amoc = model_names_amoc{i};

        % Read ZOS
        zos_data = ncread(ncfile_zos, 'zos');
        time = ncread(ncfile_zos, 'time');
        lat_zos = ncread(ncfile_zos, 'lat');
        lon_zos = ncread(ncfile_zos, 'lon');

        ref_date = datetime(1850, 1, 1);
        time = ref_date + days(time);
        time_annual = time(1:12:end);

        % Update the time filter for 1993-2013
        time_idx = time_annual >= datetime(1993, 1, 1) & time_annual <= datetime(2013, 12, 31);
        zos_data_a = zos_data(:, :, time_idx); % Filter ZOS data for the new time range
        zos_data = zos_data(:, :, :);
        time_annual = time_annual(time_idx);  % Update the time vector

        zos_data = squeeze(zos_data); % [lon x lat x time]
        zos_data = zos_data - mean(zos_data, 3, 'omitnan');

        lon_range_zos_1 = find(lon_zos >= 270 & lon_zos <= 360); %270°W and 360°W (the western part of the North Atlantic).
        lat_range_zos_1 = find(lat_zos >= 0 & lat_zos <= 40);    %0°N and 40°N (lower latitudes).
        lat_range_zos_2 = find(lat_zos >= 40 & lat_zos <= 60);   %40°N and 60°N (higher latitudes).

        zos_mean_1 = squeeze(mean(mean(zos_data(lon_range_zos_1, lat_range_zos_1, :), 1, 'omitnan'), 2, 'omitnan')); %zos in lower lat
        zos_mean_2 = squeeze(mean(mean(zos_data(lon_range_zos_1, lat_range_zos_2, :), 1, 'omitnan'), 2, 'omitnan')); %zos in higher lat
        zos_index = zos_mean_1 - zos_mean_2; %diff btw lower and higher lat
        zos_index = squeeze(mean(reshape(zos_index, 12, []), 1));
        zos_index = zos_index - mean(zos_index, 'omitnan');
        zos_index=zos_index(144:164);%for 1993 to 2013

        lon_range_zos_2 = find(lon_zos >= 290 & lon_zos <= 320);%70°W and 40°W (a different section of the North Atlantic).
        lon_range_zos_2a = find(lon_zos >= 300 & lon_zos <= 330);%60°W and 30°W (another coastal section).
        lat_range_zos_3 = find(lat_zos >= 20 & lat_zos <= 40); %20°N and 40°N (coastal region in the lower latitudes).
        lat_range_zos_4 = find(lat_zos >= 40 & lat_zos <= 60); %40°N and 60°N (coastal region in the higher latitudes).

        zos_mean_coast_1 = squeeze(mean(mean(zos_data(lon_range_zos_2, lat_range_zos_3, :), 1, 'omitnan'), 2, 'omitnan'));
        zos_mean_coast_2 = squeeze(mean(mean(zos_data(lon_range_zos_2a, lat_range_zos_4, :), 1, 'omitnan'), 2, 'omitnan'));
        zos_coast_index1 = zos_mean_coast_1 - zos_mean_coast_2;
        zos_coast_index = squeeze(mean(reshape(zos_coast_index1, 12, []), 1));
        zos_coast_index = zos_coast_index - mean(zos_coast_index, 'omitnan');
        zos_coast_index=zos_coast_index(144:164);%for 1993 to 2013

        zos_coast_index_all = [zos_coast_index_all; zos_coast_index];
        zos_index_all = [zos_index_all; zos_index];
    end
end

% Load Data for Each Model
for m = 1:length(models_name)
    for b = 1:length(lat_bands)
        lat_band = lat_bands{b};
        anom_models{m, b} = NaN(length(yr), 1); % Initialize with NaNs

        % Check for file existence
        found_file = false;
        for i = 1:length(grid_options)
            fn = sprintf('Atlantic_trans_None_%s_historical_r1i1p1f1_%s_230122.nc', models_name{m}, grid_options{i});
            if exist(fn, 'file')
                found_file = true; break;
            end
        end
        if ~found_file, continue; end

        try
            vlat = ncread(fn, 'vlat');
            moc_section = ncread(fn, 'moc_section');
            fil_val = ncinfo(fn).Variables(1).FillValue;
            moc_section(moc_section == fil_val) = NaN;
            moc_data = squeeze(mean(reshape(moc_section, size(moc_section,1), size(moc_section,2), 12, []), 3));

              
            lati = find(vlat >= lat_band(1)-0.5 & vlat <= lat_band(2)+0.5);
            if ~isempty(lati)
                
               
                % AMOC_full = max(moc_data(lati,:,:), [], 1); % max over lat
                AMOC_full = max(moc_data(lati,:,:), [], 1); % max over lat b
                AMOC_ts = squeeze(max(AMOC_full, [], 2));   % max over depth
                AMOC_anom = AMOC_ts - mean(AMOC_ts, 'omitnan');
                AMOC_anom=AMOC_anom(144:164);%for 1993 to 2013
                anom_models{m, b} = smoothdata(AMOC_anom, 'movmean', 5, 'omitnan');
            end
        catch
            anom_models{m, b} = NaN(length(yr), 1);
        end
    end
end

  % amoc_max_all = NaN(length(models_name), length(yr));
for m = 1:length(models_name)
    amoc_max_all(m, :) = anom_models{m, 1};  % assuming one lat_band
end

%% Plot Each Model Separately with Correlation Coefficients, Regression Line, and Significance Test
%% GRL-Style Multi-Panel Figure: AMOC vs ZOS Indices
figure;
set(gcf, 'Position', [100, 100, 1400, 1000]); % Higher resolution layout
tiledlayout(5, 4, 'TileSpacing', 'compact', 'Padding', 'tight'); 

fontSize = 11;

for i = 1:length(model_names_zos)
    nexttile;

    % Data
    x = amoc_max_all(i, :);
    y1 = zos_index_all(i, :);
    y2 = zos_coast_index_all(i, :);

    % Correlations
    [r1, p1] = corr(x', y1', 'Rows', 'pairwise');
    [r2, p2] = corr(x', y2', 'Rows', 'pairwise');

    % Scatter plots
    hold on;
    s1 = scatter(x, y1, 28, 'filled', 'MarkerFaceColor', [0 0.45 0.74], 'MarkerEdgeColor', 'k');
    s2 = scatter(x, y2, 28, 'filled', 'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerEdgeColor', 'k');

    % Regression lines
    coeffs1 = polyfit(x, y1, 1);
    yfit1 = polyval(coeffs1, x);
    plot(x, yfit1, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);

    coeffs2 = polyfit(x, y2, 1);
    yfit2 = polyval(coeffs2, x);
    plot(x, yfit2, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);

    % Correlation text
    textLocY = 0.085;
    text(min(x)+0.1, textLocY, sprintf('r_1 = %.2f%s', r1, ternary(p1<0.05, ' *', '')), ...
         'FontSize', fontSize, 'Color', [0 0.45 0.74]);
    text(min(x)+0.1, textLocY-0.025, sprintf('r_2 = %.2f%s', r2, ternary(p2<0.05, ' *', '')), ...
         'FontSize', fontSize, 'Color', [0.85 0.33 0.1]);

    % Axis labels and limits
    xlabel('AMOC max (Sv)', 'FontSize', fontSize);
    if ismember(i, [1, 5, 9, 13, 17])
        ylabel('ZOS anomaly (m)', 'FontSize', fontSize);
    end
    xlim([-5 3]);
    ylim([-0.1 0.1]);

    % Title
    title(models_name{i}, 'FontSize', fontSize, 'FontWeight', 'bold', 'Interpreter', 'none');

    % Formatting
    % grid on;
    box on;
    set(gca, 'FontSize', fontSize - 1, 'LineWidth', 1);

    if i == 1
        legend([s1 s2], {'ZOS index', 'ZOS coast index'}, 'FontSize', fontSize-2, ...
               'Location', 'southoutside', 'Orientation', 'horizontal', ...
               'Box', 'off');
    end
end

function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end

%% Create figure for ZOS Index vs AMOC Max categorization
% Create figure for ZOS Index vs AMOC Max
category_indices = {AMOC_strong, AMOC_weak, Eqt_pos, Eqt_neg, Eqt_bipol, wb_np, wb_nn, wb_pn};
category_names = {'AMOC Strong', 'AMOC Weak', 'Eqt Positive', 'Eqt Negative', 'Eqt Bipolar', 'Wb np', 'Wb nn', 'Wb pn'};

figure('Position', [100, 100, 1300, 1000], 'Color', 'w');
sgtitle('Classification of AMOC Max vs ZOS Indices', ...
    'FontSize', 16, 'FontWeight', 'bold');

for cat_idx = 1:length(category_indices)
    models_in_category = category_indices{cat_idx};
    model_indices = ismember(models_name, models_in_category);
    zos_index_cat = zos_index_all(model_indices, :);
    zos_coast_index_cat = zos_coast_index_all(model_indices, :);
    amoc_max_cat = amoc_max_all(model_indices, :);
    
    subplot(3, 3, cat_idx);
    hold on;

    % === Correlation calculations ===
    [r1, p1] = corr(amoc_max_cat(:), zos_coast_index_cat(:), 'Rows', 'pairwise');
    [r2, p2] = corr(amoc_max_cat(:), zos_index_cat(:), 'Rows', 'pairwise');

    % === Scatter plots ===
    scatter(amoc_max_cat(:), zos_index_cat(:), 45, ...
        'MarkerFaceColor', [0 0.45 0.74], 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'ZOS Index');
    scatter(amoc_max_cat(:), zos_coast_index_cat(:), 45, ...
        'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'ZOS Coast Index');

    % === Regression lines ===
    coeffs1 = polyfit(amoc_max_cat(:), zos_index_cat(:), 1);
    y_fit1 = polyval(coeffs1, amoc_max_cat(:));
    plot(amoc_max_cat(:), y_fit1, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);

    coeffs2 = polyfit(amoc_max_cat(:), zos_coast_index_cat(:), 1);
    y_fit2 = polyval(coeffs2, amoc_max_cat(:));
    plot(amoc_max_cat(:), y_fit2, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);

    % === Annotate correlation coefficients ===
    xpos = min(amoc_max_cat(:)) + 0.2;

if p1 < 0.05
    str1 = sprintf('r_1 = %.2f (p < 0.05)', r1);
else
    str1 = sprintf('r_1 = %.2f (p = %.3f)', r1, p1);
end
text(xpos, 0.1, str1, 'FontSize', 10, 'Color', [0.85 0.33 0.1]);

if p2 < 0.05
    str2 = sprintf('r_2 = %.2f (p < 0.05)', r2);
else
    str2 = sprintf('r_2 = %.2f (p = %.3f)', r2, p2);
end
text(xpos, 0.08, str2, 'FontSize', 10, 'Color', [0 0.45 0.74]);

    % === Axes labels & limits ===
    xlabel('AMOC Max (Sv)', 'FontSize', 11);
    if ismember(cat_idx, [1, 4, 7])
        ylabel('ZOS Anomaly (m)', 'FontSize', 11);
    end
    xlim([-5 4]);
    ylim([-0.15 0.15]);

    % === Title & formatting ===
    title(category_names{cat_idx}, 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 11);
    grid on;
    box on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);

    % Optional: add legend only to first subplot
    if cat_idx == 1
        legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9);
    end
end

    