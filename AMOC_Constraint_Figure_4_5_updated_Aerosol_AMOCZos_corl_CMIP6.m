clear;clc;

%% Load data
load('AMOC_change_26p5.mat')
load('ZOS_change_26p5.mat')   
load('lat.mat');
load('lat_range.mat');
load('lon.mat');
load('lon_range.mat');
load('na.dat');  % For land outline
% CMIP6 Model Files and Names
%%
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
yr = 1850:2014; % Full year range
lat_bands = {[25, 26]};
anom_models = cell(length(models_name), length(lat_bands));

% Define categories as cell arrays
AMOC_strong = {'CESM2', 'CESM2-WACCM', 'GISS-E2-1-G', 'NorESM2-LM'};
AMOC_weak = {'NESM3', 'MIROC6', 'IPSL-CM6A-LR', 'INM-CM5-0', 'BCC-CSM2-MR'};

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

        % zos_data = squeeze(zos_data); % [lon x lat x time]
        % zos_data = zos_data - mean(zos_data, 3, 'omitnan');

        zos_data = squeeze(zos_data); % [lon x lat x time]
        [nlon, nlat, ntime] = size(zos_data);
   for t = 1:ntime
    snapshot = zos_data(:, :, t);
    global_mean = mean(snapshot(:), 'omitnan');
    zos_data(:, :, t) = snapshot - global_mean;
   end

        lon_range_zos_1 = find(lon_zos >= 270 & lon_zos <= 360); %270°W and 360°W (the western part of the North Atlantic).
        lat_range_zos_1 = find(lat_zos >= 0 & lat_zos <= 40);    %0°N and 40°N (lower latitudes).
        lat_range_zos_2 = find(lat_zos >= 40 & lat_zos <= 60);   %40°N and 60°N (higher latitudes).

        zos_mean_1 = squeeze(mean(mean(zos_data(lon_range_zos_1, lat_range_zos_1, :), 1, 'omitnan'), 2, 'omitnan')); %zos in lower lat
        zos_mean_2 = squeeze(mean(mean(zos_data(lon_range_zos_1, lat_range_zos_2, :), 1, 'omitnan'), 2, 'omitnan')); %zos in higher lat
        zos_index = zos_mean_1 - zos_mean_2; %diff btw lower and higher lat
        zos_index = squeeze(mean(reshape(zos_index, 12, []), 1));
        zos_index = zos_index - mean(zos_index, 'omitnan');

        lon_range_zos_2 = find(lon_zos >= 290 & lon_zos <= 320);%70°W and 40°W (a different section of the North Atlantic).
        lon_range_zos_2a = find(lon_zos >= 300 & lon_zos <= 330);%60°W and 30°W (another coastal section).
        lat_range_zos_3 = find(lat_zos >= 20 & lat_zos <= 40); %20°N and 40°N (coastal region in the lower latitudes).
        lat_range_zos_4 = find(lat_zos >= 40 & lat_zos <= 60); %40°N and 60°N (coastal region in the higher latitudes).

        zos_mean_coast_1 = squeeze(mean(mean(zos_data(lon_range_zos_2, lat_range_zos_3, :), 1, 'omitnan'), 2, 'omitnan'));
        zos_mean_coast_2 = squeeze(mean(mean(zos_data(lon_range_zos_2a, lat_range_zos_4, :), 1, 'omitnan'), 2, 'omitnan'));
        zos_coast_index = zos_mean_coast_1 - zos_mean_coast_2;
        zos_coast_index = squeeze(mean(reshape(zos_coast_index, 12, []), 1));
        zos_coast_index = zos_coast_index - mean(zos_coast_index, 'omitnan');
        
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
               
                AMOC_full = max(moc_data(lati,:,:), [], 1); % max over lat
                AMOC_ts = squeeze(max(AMOC_full, [], 2));   % max over depth
                AMOC_anom = AMOC_ts - mean(AMOC_ts, 'omitnan');
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


%% Improved GRL-Style Figure: AMOC Strength vs ZOS Indices
category_indices = {AMOC_strong, AMOC_weak};
group_names = {'AMOC Strong', 'AMOC Weak'};
model_groups = {
    {'CESM2', 'CESM2-WACCM', 'GISS-E2-1-G', 'NorESM2-LM'}, ...
    {'NESM3', 'MIROC6', 'IPSL-CM6A-LR', 'INM-CM5-0', 'BCC-CSM2-MR'}
};

figure('Units', 'inches', 'Position', [1, 1, 12, 10], 'Color', 'w');

% --- Top row: Spatial correlation maps ---
for g = 1:length(group_names)
    subplot(2, 2, g); % Positions 1 and 2
    
    % Data gathering and correlation calculations stay the same as before
    AMOC_all = [];
    ZOS_all = [];

    for i = 1:length(model_groups{g})
        model_name = strrep(model_groups{g}{i}, '-', '_'); 
        if isfield(AMOC_data_now, model_name) && isfield(zos_avg_now, model_name)
            AMOC_all = [AMOC_all, AMOC_data_now.(model_name).AMOC_change];
            ZOS_all = cat(3, ZOS_all, zos_avg_now.(model_name).zos_change_nw);
        end
    end

    [nx, ny, nt] = size(ZOS_all);
    AMOC_tsa = NaN(nx, ny);
    p_values = NaN(nx, ny);

    for j = 1:ny
        for i = 1:nx
            [r, p] = corr(squeeze(AMOC_all(:)), squeeze(ZOS_all(i, j, :)), 'rows', 'complete');
            AMOC_tsa(i, j) = r;
            p_values(i, j) = p;
        end
    end

    sig_mask = p_values <= 0.05;

    hold on;
    pcolor(lon(lon_range), lat(lat_range), AMOC_tsa');  
    shading flat;
    title([group_names{g}, ' Correlation']);
    xlabel('Longitude'); ylabel('Latitude');

    [C, h] = contour(lon(lon_range), lat(lat_range), AMOC_tsa', 'k', 'LineWidth', 1.0);
    clabel(C, h, 'FontSize', 8);

    [sig_x, sig_y] = meshgrid(lon(lon_range), lat(lat_range));
    scatter(sig_x(sig_mask'), sig_y(sig_mask'), 5, 'k', 'filled');

    plot(na(:,1) + 360, na(:,2), 'k', 'LineWidth', 1.5);

    xlim([270 360]);
    ylim([0 60]);
    grid on; box on;
    set(gca, 'FontSize', 13);
    set(gca, 'XTick', 270:20:360, 'XTickLabel', {'90°W', '70°W', '50°W', '30°W', '10°W'});
    set(gca, 'YTick', 0:20:80, 'YTickLabel', {'0', '20°N', '40°N', '60°N', '80°N'});

    cmap = makeColorMap([0 0 1],[1 1 1],[1 0 0],10);
    colormap(gca, cmap);
    caxis([-0.5 0.5]);
end

% Add colorbar for top row spatial plots
cb = colorbar('Position', [0.92 0.55 0.02 0.35]); % Adjust position for top half
cb.Label.String = 'Correlation coefficient';
cb.FontSize = 13;

% --- Bottom row: Scatter plots ---
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];

for cat_idx = 1:length(category_indices)
    subplot(2, 2, cat_idx + 2); % Positions 3 and 4
    hold on;

    models_in_category = category_indices{cat_idx};
    model_indices = ismember(models_name, models_in_category);
    zos_index_cat = zos_index_all(model_indices, :);
    zos_coast_index_cat = zos_coast_index_all(model_indices, :);
    amoc_max_cat = amoc_max_all(model_indices, :);

    s1 = scatter(amoc_max_cat(:), zos_index_cat(:), 50, ...
        'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.8, 'DisplayName', 'ZOS Index');
    s2 = scatter(amoc_max_cat(:), zos_coast_index_cat(:), 50, ...
        'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.8, 'DisplayName', 'ZOS Coast Index');

    xvals = linspace(min(amoc_max_cat(:)), max(amoc_max_cat(:)), 100);
    coeffs1 = polyfit(amoc_max_cat(:), zos_index_cat(:), 1);
    plot(xvals, polyval(coeffs1, xvals), '-', 'Color', colors(1,:), 'LineWidth', 2);

    coeffs2 = polyfit(amoc_max_cat(:), zos_coast_index_cat(:), 1);
    plot(xvals, polyval(coeffs2, xvals), '-', 'Color', colors(2,:), 'LineWidth', 2);

    [r1, p1] = corr(amoc_max_cat(:), zos_index_cat(:), 'Rows', 'pairwise');
    [r2, p2] = corr(amoc_max_cat(:), zos_coast_index_cat(:), 'Rows', 'pairwise');

    if p1 < 0.05
        txt1 = sprintf('r_{ZOS} = %.2f (p < 0.05)', r1);
    else
        txt1 = sprintf('r_{ZOS} = %.2f (p = %.2f)', r1, p1);
    end

    if p2 < 0.05
        txt2 = sprintf('r_{ZOS coast} = %.2f (p < 0.05)', r2);
    else
        txt2 = sprintf('r_{ZOS coast} = %.2f (p = %.2f)', r2, p2);
    end

    x_text = min(xvals) + 0.5;
    ylims = ylim;
    text(x_text, ylims(2) - 0.015, txt1, 'FontSize', 12, ...
        'Color', colors(1,:), 'FontWeight', 'normal', 'FontName', 'Arial');
    text(x_text, ylims(2) - 0.045, txt2, 'FontSize', 12, ...
        'Color', colors(2,:), 'FontWeight', 'normal', 'FontName', 'Arial');

    xlabel('AMOC Maximum Anomaly (Sv)', 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Arial');
    if cat_idx == 1
        ylabel('ZOS Anomaly (m)', 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Arial');
    end
    title(group_names{cat_idx}, 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial');

    xlim([-5 4]);
    ylim([-0.15 0.15]);
    set(gca, 'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 1.2);
    xticks(-5:1:4);
    yticks(-0.15:0.05:0.15);

    box on;

    xlims = xlim;
    ylims = ylim;
    text(xlims(1) - 0.7, ylims(2) - 0.01, ['(' char('a' + cat_idx + 1) ')'], ...
        'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial');
end

% Centralized legend below scatter plots
legend([s1, s2], {'ZOS Index', 'ZOS Coast Index'}, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', ...
    'Box', 'off', 'FontSize', 12, 'FontName', 'Arial');

% Overall figure title
sgtitle({'Relationship Between AMOC Anomalies and ZOS Indices'}, ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
