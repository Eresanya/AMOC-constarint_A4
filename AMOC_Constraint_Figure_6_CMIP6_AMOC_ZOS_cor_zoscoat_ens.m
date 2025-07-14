clear; clc;

%% Define model files for AMOC and ZOS
model_filenames_amoc = { ...
    'Harmonized_Atlantic_trans_None_ACCESS-CM2_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_BCC-CSM2-MR_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_CAMS-CSM1-0_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_CanESM5_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_CESM2_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_CESM2-FV2_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_CESM2-WACCM_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_CESM2-WACCM-FV2_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_CMCC-CM2-HR4_historical_r1i1p1f1_gn_230122.nc', ...
    'Harmonized_Atlantic_trans_None_E3SM-1-0_historical_r1i1p1f1_gr_230122.nc', ...
    'Harmonized_Atlantic_trans_None_E3SM-1-1_historical_r1i1p1f1_gr_230122.nc', ...  
    'Harmonized_Atlantic_trans_None_FGOALS-f3-L_historical_r1i1p1f1_gn_230122.nc', ...   
    'Harmonized_Atlantic_trans_None_GISS-E2-1-G_historical_r1i1p1f1_gn_230122.nc', ...  
    'Harmonized_Atlantic_trans_None_INM-CM4-8_historical_r1i1p1f1_gr1_230122.nc', ...  
    'Harmonized_Atlantic_trans_None_INM-CM5-0_historical_r1i1p1f1_gr1_230122.nc', ... 
    'Harmonized_Atlantic_trans_None_IPSL-CM6A-LR_historical_r1i1p1f1_gn_230122.nc', ...  
    'Harmonized_Atlantic_trans_None_MIROC6_historical_r1i1p1f1_gn_230122.nc', ... 
    'Harmonized_Atlantic_trans_None_NESM3_historical_r1i1p1f1_gn_230122.nc', ...  
    'Harmonized_Atlantic_trans_None_NorESM2-LM_historical_r1i1p1f1_gr_230122.nc'};

model_filenames_zos = {...
    'zos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_1850-2014.nc', ... 
    'zos_Omon_BCC-CSM2-MR_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_CAMS-CSM1-0_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_CanESM5_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_CESM2_historical_r1i1p1f1_gr_1850-2014.nc', ...  
    'zos_Omon_CESM2-FV2_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_CESM2-WACCM_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_CESM2-WACCM-FV2_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_CMCC-CM2-HR4_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_E3SM-1-0_historical_r1i1p1f1_gr_1850-2014.nc', ... 
    'zos_Omon_E3SM-1-1_historical_r1i1p1f1_gr_1850-2014.nc', ...  
    'zos_Omon_FGOALS-f3-L_historical_r1i1p1f1_gn_1850-2014.nc', ...   
    'zos_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_INM-CM4-8_historical_r1i1p1f1_gr1_1850-2014.nc', ...  
    'zos_Omon_INM-CM5-0_historical_r1i1p1f1_gr1_1850-2014.nc', ... 
    'zos_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_MIROC6_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_NESM3_historical_r1i1p1f1_gn_1850-2014.nc', ...  
    'zos_Omon_NorESM2-LM_historical_r1i1p1f1_gn_1850-2014.nc'};

model_names_zoscoast = {'CESM2-WACCM','CESM2-FV2','IPSL-CM6A-LR'};
model_ensm = {...
    'ACCESS-CM2', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5', 'CESM2', ...
    'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV', 'CMCC-CM2-HR4', ...
    'E3SM-1-0', 'E3SM-1-1', 'FGOALS-f3-L','GISS-E2-1-G', ...
    'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MIROC', 'NESM3', 'NorESM2-LM'};

model_groups = {model_ensm,model_names_zoscoast};
group_titles = {'CMIP6 Ensemble','CMIP6 Constraint'};


for g = 1:2  % Loop over zoscoast and zosindex groups
    model_names = model_groups{g};
    ensemble_AMOC_change_all = [];
    ensemble_ZOS_change_nw_all = [];

    for m = 1:length(model_names)
        model_name = model_names{m};

        % Find matching index
        idx = find(contains(model_filenames_amoc, model_name));

        if isempty(idx)
            warning('Model %s not found.', model_name);
            continue;
        end

        %% Load AMOC
        fn_amoc = model_filenames_amoc{idx};
        msftmz = ncread(fn_amoc, 'moc_section');
        lat_amoc = ncread(fn_amoc, 'lat');
        lev = ncread(fn_amoc, 'lev');
        year_model = ncread(fn_amoc, 'time');
        time_model = year_model + datenum(1800, 01, 01);
        dat_model = datevec(time_model);
        yr = dat_model(:, 1);

        % lat_amoc_range = lat_amoc >= -20 & lat_amoc <= 60;
        lat_amoc_range = lat_amoc >= 26 & lat_amoc <= 26.5;
        msftmz_avg = squeeze(mean(reshape(msftmz, size(msftmz,1), size(msftmz,2), 12, []), 3));
        msftmz_avg = permute(msftmz_avg, [2 1 3]);
        AMOC_data = msftmz_avg(:, lat_amoc_range, :);
        AMOC_ = squeeze(tsnanmean(AMOC_data, 1));

        yr_amoc = 1850:2010;
        AMOC_yr1 = AMOC_(:, yr_amoc >= 1850 & yr_amoc <= 1900);
        AMOC_yr2 = AMOC_(:, yr_amoc >= 1960 & yr_amoc <= 2010);
        AMOC_change = squeeze(mean(AMOC_yr1 - AMOC_yr2, 1));

        %% Load ZOS
        fn_zos = model_filenames_zos{idx};
        zos_data = ncread(fn_zos, 'zos');
        lon = ncread(fn_zos, 'lon');
        lat = ncread(fn_zos, 'lat');
        year_zos = ncread(fn_zos, 'time');
        time_zos = year_zos + datenum(1850, 01, 01);
        dat_zos = datevec(time_zos);

        lon_range = lon >= 270 & lon <= 360;
        lat_range = lat >= 0 & lat <= 60;

        zos_avg = squeeze(mean(reshape(zos_data, size(zos_data,1), size(zos_data,2), 12, []), 3));
        zos_region = zos_avg(lon_range, lat_range, :);

        yr_zos = 1850:2010;
        zos_yr1 = zos_region(:, :, yr_zos >= 1850 & yr_zos <= 1900);
        zos_yr2 = zos_region(:, :, yr_zos >= 1960 & yr_zos <= 2010);
        zos_change = zos_yr1 - zos_yr2;

        zos_change_nw = zeros(size(zos_change));
        for i = 1:size(zos_change, 3)
            tmp = zos_change(:, :, i);
            zos_change_nw(:, :, i) = tmp - mean(tmp(:), 'omitnan');
        end

        ensemble_AMOC_change_all = cat(3, ensemble_AMOC_change_all, AMOC_change);
        ensemble_ZOS_change_nw_all = cat(4, ensemble_ZOS_change_nw_all, zos_change_nw);
    end

    ensemble_AMOC_changes{g} = mean(ensemble_AMOC_change_all, 3, 'omitnan');
    ensemble_ZOS_changes{g} = mean(ensemble_ZOS_change_nw_all, 4, 'omitnan');
end

% Visualization
figure;
for g = 1:2
    AMOC_change_ensemble = ensemble_AMOC_changes{g};
    ZOS_change_ensemble = ensemble_ZOS_changes{g};

    [nx, ny, nt] = size(ZOS_change_ensemble);
    ensemble_corr = zeros(nx, ny);
    p_values_ensemble = ones(nx, ny);

    for j = 1:ny
        for i = 1:nx
            y = squeeze(ZOS_change_ensemble(i,j,:));
            x = AMOC_change_ensemble(:);
            if all(isfinite(x)) && all(isfinite(y))
                [r, p] = corr(x, y, 'Rows', 'complete');
                ensemble_corr(i,j) = r;
                p_values_ensemble(i,j) = p;
            else
                ensemble_corr(i,j) = NaN;
                p_values_ensemble(i,j) = NaN;
            end
        end
    end

    sig_mask = p_values_ensemble <= 0.05;

    subplot(1,2,g);
    pcolor(lon(lon_range), lat(lat_range), ensemble_corr'); shading flat;
    hold on;
    contour(lon(lon_range), lat(lat_range), ensemble_corr', [-1:0.2:1], '-k', 'LineWidth', 1.0);
    load('na.dat');
    plot(na(:,1)+360, na(:,2), 'k', 'LineWidth', 1.5);
    [x_sig, y_sig] = meshgrid(lon(lon_range), lat(lat_range));
    scatter(x_sig(sig_mask'), y_sig(sig_mask'), 5, 'k', 'filled');

    xlim([270 360]); ylim([0 60]);
    set(gca, 'FontSize', 13, 'YTick', 0:20:80, 'XTick', 270:20:360, ...
             'YTickLabel', {'0', '20°N', '40°N', '60°N', '80°N'}, ...
             'XTickLabel', {'90°W', '70°W', '50°W', '30°W', '10°W'});
    title(group_titles{g}, 'FontWeight', 'bold');
    cmap = makeColorMap([0 0 1],[1 1 1],[1 0 0],10);
    colormap(cmap);
    caxis([-0.5 0.5]);
end

cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Correlation coefficient';
cb.FontSize = 13;
