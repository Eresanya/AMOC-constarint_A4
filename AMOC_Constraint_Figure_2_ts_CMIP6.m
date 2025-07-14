clear; clc;

%% CMIP6 Model Files and Names
AMOC_files = {...
    'Harmonized_Atlantic_trans_None_ACCESS-CM2_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_BCC-CSM2-MR_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_CAMS-CSM1-0_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_CanESM5_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_CESM2_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_CESM2-FV2_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_CESM2-WACCM_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_CESM2-WACCM-FV2_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_CMCC-CM2-HR4_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_E3SM-1-0_historical_r1i1p1f1_gr_230122.nc',... 
    'Harmonized_Atlantic_trans_None_E3SM-1-1_historical_r1i1p1f1_gr_230122.nc',...  
    'Harmonized_Atlantic_trans_None_FGOALS-f3-L_historical_r1i1p1f1_gn_230122.nc',...   
    'Harmonized_Atlantic_trans_None_GISS-E2-1-G_historical_r1i1p1f1_gn_230122.nc',...  
    'Harmonized_Atlantic_trans_None_INM-CM4-8_historical_r1i1p1f1_gr1_230122.nc',...  
    'Harmonized_Atlantic_trans_None_INM-CM5-0_historical_r1i1p1f1_gr1_230122.nc',... 
    'Harmonized_Atlantic_trans_None_IPSL-CM6A-LR_historical_r1i1p1f1_gn_230122.nc',...  
    'Harmonized_Atlantic_trans_None_MIROC6_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_NESM3_historical_r1i1p1f1_gn_230122.nc',... 
    'Harmonized_Atlantic_trans_None_NorESM2-LM_historical_r1i1p1f1_gr_230122.nc'};

model_names = {...
    'ACCESS-CM2', 'BCC-CSM2-MR', 'CAMS-CSM2', 'CanESM5', 'CESM2', ...
    'CESM2-FV2', 'CESM2-WACCM2', 'CESM2-WACCM-FV', 'CMCC-CM2-HR42', ...
    'E3SM-1-0', 'E3SM-1-1', 'FGOALS-f3-L2','GISS-E2-2-G2', ...
    'INM-CM4-8', 'INM-CM5', 'IPSL-CM6A-LR', 'MIROC', 'NESM3', 'NorESM2-LM'};

% Initialize ensemble mean storage
ensemble_moc_section_mean = [];

for i = 1:length(AMOC_files)
    fn = AMOC_files{i};
    ncinf = ncinfo(fn);
    lat = ncread(fn, 'lat');
    lev  = ncread(fn, 'lev');
    year_modelsa = ncread(fn, 'time');
    time = year_modelsa + datenum(1850,1,1);

    fil_val = ncinf.Variables(1).FillValue;
    lat(lat == fil_val) = NaN;
    lev(lev == fil_val) = NaN;

    moc_section = ncread(fn, 'moc_section');
    moc_section(moc_section == fil_val) = NaN;

    dims = size(moc_section);
    n_years = floor(dims(3) / 12);
    moc_section_reshaped = moc_section(:,:,1:(n_years*12));
    annual_mean = squeeze(mean(reshape(moc_section_reshaped, dims(1), dims(2), 12, n_years), 3));
    data2 = permute(annual_mean,[2 1 3]);
    moc_section_mean = squeeze(mean(data2, 3));

    if isempty(ensemble_moc_section_mean)
        ensemble_moc_section_mean = moc_section_mean;
    else
        ensemble_moc_section_mean = moc_section_mean;
    end
end

ensemble_moc_section_mean = ensemble_moc_section_mean / length(AMOC_files);
interp_ensemble_moc_section = interp1(lev, ensemble_moc_section_mean, 100:100:4000);

subplot(5, 4, 1);
hold on;
sgtitle('MOC Sections for 19 CMIP6 Models and Ensemble Mean','FontSize', 10);

% ----------------------- Anomaly Analysis -----------------------
models = {...
    'ACCESS-CM2', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5', 'CESM2', ...
    'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV2', 'CMCC-CM2-HR4', ...
    'E3SM-1-0', 'E3SM-1-1', 'FGOALS-f3-L','GISS-E2-1-G', ...
    'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MIROC6', 'NESM3', 'NorESM2-LM'};

grid_options = {'gn', 'gr1', 'gr'};
yr = 1850:2014;
lat_bands = {[26, 26.5]};
colors = [0 0 0; 1 0 0; 0 0 1];

anom_models = cell(length(models), length(lat_bands));
ensemble_mean = cell(length(lat_bands), 1);
for b = 1:length(lat_bands)
    ensemble_mean{b} = NaN(length(yr), 1);
end

for m = 1:length(models)
    for b = 1:length(lat_bands)
        lat_band = lat_bands{b};
        anom_models{m, b} = NaN(length(yr), 1);

        found_file = false;
        for i = 1:length(grid_options)
            fn = sprintf('Harmonized_Atlantic_trans_None_%s_historical_r1i1p1f1_%s_230122.nc', models{m}, grid_options{i});
            if exist(fn, 'file')
                found_file = true; break;
            end
        end
        if ~found_file, continue; end

        try
            lat = ncread(fn, 'lat');
            moc_section = ncread(fn, 'moc_section');
            fil_val = ncinfo(fn).Variables(1).FillValue;
            moc_section(moc_section == fil_val) = NaN;

            moc_data = squeeze(mean(reshape(moc_section, size(moc_section,1), size(moc_section,2), 12, []), 3));
            lati = find(lat >= lat_band(1)-0.5 & lat <= lat_band(2)+0.5);

            if ~isempty(lati)
                AMOC_full = squeeze(max(moc_data(lati,:,:), [], 1));
                AMOC_ts = squeeze(max(AMOC_full));
                AMOC_anom = AMOC_ts - mean(AMOC_ts, 'omitnan');
                anom_models{m, b} = smoothdata(AMOC_anom, 'movmean', 5, 'omitnan');
            end
        catch
            anom_models{m, b} = NaN(length(yr), 1);
        end
    end
end

for b = 1:length(lat_bands)
    all_data = NaN(length(yr), length(models));
    for m = 1:length(models)
        if ~isempty(anom_models{m, b})
            all_data(:, m) = anom_models{m, b};
        end
    end
    valid_counts = sum(~isnan(all_data), 2);
    valid_counts(valid_counts == 0) = NaN;
    ensemble_mean{b} = sum(all_data, 2, 'omitnan') ./ valid_counts;
end

% figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 9]);
set(gcf, 'Color', 'w');

% Font & visual settings
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);
set(groot, 'defaultAxesLineWidth', 1);
set(groot, 'defaultLineLineWidth', 1.2);

% Colorblind-friendly palette
colors = [0.1216, 0.4667, 0.7059;  % blue
          1.0000, 0.4980, 0.0549;  % orange
          0.1725, 0.6275, 0.1725]; % green

% Layout
t = tiledlayout(5, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, 'AMOC Anomalies in 19 CMIP6 Models and Ensemble Mean', ...
    'FontSize', 12, 'FontWeight', 'bold');

% Panel labels (a) to (t)
panel_labels = arrayfun(@(x) ['(', char(96 + x), ')'], 1:20, 'UniformOutput', false);

% === Panel 1: Ensemble Mean ===
nexttile;
hold on;
for b = 1:length(lat_bands)
    plot(yr, ensemble_mean{b}, 'Color', colors(b,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%g–%g°N', lat_bands{b}));
end
yline(0, 'k--', 'HandleVisibility', 'off');
legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
       'FontSize', 8, 'Box', 'off');
title('Ensemble Mean', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('AMOC anomaly (Sv)', 'FontWeight', 'bold');
xlabel('Year', 'FontWeight', 'bold');
set(gca, 'XTick', 1850:20:2010, 'XLim', [1850 2014], 'YLim', [-5 5]);
grid on; box on;
text(0.01, 0.90, panel_labels{1}, 'Units', 'normalized', 'FontWeight', 'bold');

% === Panels 2–20: Individual Models ===
for m = 1:length(models)
    nexttile;
    hold on;
    for b = 1:length(lat_bands)
        if ~isempty(anom_models{m, b}) && any(~isnan(anom_models{m, b}))
            plot(yr, anom_models{m, b}, 'Color', colors(b,:), 'LineWidth', 1.2);
        end
    end
    yline(0, 'k--', 'HandleVisibility', 'off');
    title(models{m}, 'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none');
    set(gca, 'XTick', 1850:40:2010, 'XLim', [1850 2014], 'YLim', [-5 5]);
    grid on; box on;

    % Label axes only where necessary
    if mod(m, 4) == 1
        ylabel('Anomaly (Sv)', 'FontWeight', 'bold');
    end
    if m >= 16
        xlabel('Year', 'FontWeight', 'bold');
    end

    % Panel label
    text(0.01, 0.90, panel_labels{m+1}, 'Units', 'normalized', 'FontWeight', 'bold');
end