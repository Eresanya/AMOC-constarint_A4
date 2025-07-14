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

model_names = {...
    'ACCESS-CM2', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CanESM5', 'CESM2', ...
    'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV', 'CMCC-CM2-HR4', ...
    'E3SM-1-0', 'E3SM-1-1', 'FGOALS-f3-L','GISS-E2-1-G', ...
    'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MIROC6', 'NESM3', 'NorESM2-LM'};

% Preallocate storage
ensemble_AMOC_change_all = [];
ensemble_ZOS_change_nw_all = [];

figure;
num_models = length(model_filenames_amoc);
t = tiledlayout(5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:num_models
    nexttile(m + 1);
    title(model_names{m}, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

    % Load and Process AMOC Data
    fn_amoc = model_filenames_amoc{m};
    msftmz = ncread(fn_amoc, 'moc_section');
    lat_amoc = ncread(fn_amoc, 'lat');
    lev = ncread(fn_amoc, 'lev');
    year_model = ncread(fn_amoc, 'time');
    time_model = year_model + datenum(1800, 01, 01);
    dat_model = datevec(time_model);
    yr = dat_model(:, 1);
    lat_amoc_range = lat_amoc >= 26 & lat_amoc <= 26.5;
    msftmz_avg = squeeze(mean(reshape(msftmz, 300, 50, 12, []), 3));
    msftmz_avg = permute(msftmz_avg, [2 1 3]);
    AMOC_data = msftmz_avg(:, lat_amoc_range, :);
    AMOC_ = squeeze(tsnanmean(AMOC_data, 1));
    yr_amoc = 1850:2010;
    AMOC_yr1 = AMOC_(:, yr_amoc >= 1850 & yr_amoc <= 1900);
    AMOC_yr2 = AMOC_(:, yr_amoc >= 1960 & yr_amoc <= 2010);
    AMOC_change = squeeze(mean(AMOC_yr1 - AMOC_yr2, 1));

    % Load and Process ZOS Data
    fn_zos = model_filenames_zos{m};
    zos_data = ncread(fn_zos, 'zos');
    lon = ncread(fn_zos, 'lon');
    lat = ncread(fn_zos, 'lat');
    year_zos = ncread(fn_zos, 'time');
    time_zos = year_zos + datenum(1850, 01, 01);
    lat_range = lat >= 0 & lat <= 60;
    lon_range = lon >= 270 & lon <= 360;
    zos_avg = squeeze(mean(reshape(zos_data, 144, 72, 12, []), 3));
    zos_region = zos_avg(lon_range, lat_range, :);
    yr_zos = 1850:2010;
    zos_yr1 = zos_region(:, :, yr_zos >= 1850 & yr_zos <= 1900);
    zos_yr2 = zos_region(:, :, yr_zos >= 1960 & yr_zos <= 2010);
    zos_change = zos_yr1 - zos_yr2;

    % Remove Global Mean from ZOS Change
    zos_change_nw = zeros(size(zos_change));
    for i = 1:size(zos_change, 3)
        tmp = zos_change(:, :, i);
        zos_change_nw(:, :, i) = tmp - mean(tmp(:), 'omitnan');
    end



    % For ensemble
    ensemble_AMOC_change_all = cat(3, ensemble_AMOC_change_all, AMOC_change);
    ensemble_ZOS_change_nw_all = cat(4, ensemble_ZOS_change_nw_all, zos_change_nw);

    % Correlation Calculation
    [nx, ny, nt] = size(zos_change_nw);
    AMOC_tsa = zeros(nx, ny);
    p_values = ones(nx, ny);
    for j = 1:ny
        for i = 1:nx
            [r, p] = corrcoef(AMOC_change', squeeze(zos_change_nw(i,j,:)));
            AMOC_tsa(i,j) = r(2);
            p_values(i,j) = p(2);
        end
    end
    sig_mask = p_values <= 0.05;

    %% Plot the Regression of Dynamic Sea Level Against AMOC
    pcolor(lon(lon_range), lat(lat_range), AMOC_tsa'); 
    shading flat; 
    hold on;

    % Add contour lines without labels
    contour(lon(lon_range), lat(lat_range), AMOC_tsa', '-k', 'LineWidth', 1.0);

    % Overlay continental boundaries
    load('na.dat'); % Make sure 'na' is loaded for boundaries
    plot(na(:,1) + 360, na(:,2), 'k', 'LineWidth', 1.5);

    % Adjust axis limits
    xlim([270 360]); % Keep within valid range
    ylim([0 60]);

    % Overlay black dots at significant locations
    [x_sig, y_sig] = meshgrid(lon(lon_range), lat(lat_range));
    scatter(x_sig(sig_mask'), y_sig(sig_mask'), 5, 'k', 'filled'); % Significant points

    % Add latitude and longitude grid lines
    grid on; 
    ax = gca;
    ax.XGrid = 'on'; % Longitude grid lines
    ax.YGrid = 'on'; % Latitude grid lines

    % Set up colormap and axis labels
    load('newcolormaps', 'bwr'); 
    colormap(ax, bwr);
    caxis([-0.5 0.5]);
  
   set(gca, 'FontSize', 13);

    % Set y-axis labels only for specific indices (1, 5, 9, 13, 17)
    if ismember(m, [4, 8, 12, 16])
        set(gca, 'YTick', 0:20:60, 'YTickLabel', {'0', '20°N', '40°N', '60°N'});
    else
        set(gca, 'YTickLabel', []);  % Remove y-axis labels for other indices
    end

    % Remove y-axis labels for specific indices (2, 3, 4, 13, 14, 15)
    if ismember(m, [1,2,3,5,6,7,9,10,11,13,14,15,17,18,19])
        set(gca, 'YTickLabel', []);  % No y-axis tick labels
    end

  % Set x-axis labels only for specific indices (17,18,19)
   if ismember(m, [16,17,18,19])
        set(gca, 'XTick', 270:20:360, 'XTickLabel', {'90°W', '70°W', '50°W', '30°W', '10°W'});
    else
        set(gca, 'XTickLabel', []);  % Remove x-axis labels for other indices
    end

    % Remove x-axis labels for specific indices (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
    if ismember(m, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
        set(gca, 'XTickLabel', []);  % No y-axis tick labels
    end   
    
    % Title and grid
    title([model_names{m}],'FontSize', 12);
    box on;
    grid on;

  end
  
% Compute ensemble means
AMOC_change_ensemble = mean(ensemble_AMOC_change_all, 3, 'omitnan');
ZOS_change_ensemble = mean(ensemble_ZOS_change_nw_all, 4, 'omitnan');

% Compute correlation
[nx, ny, nt] = size(ZOS_change_ensemble);
ensemble_corr = zeros(nx, ny);
p_values_ensemble = ones(nx, ny);


for j = 1:ny
    for i = 1:nx
        y = squeeze(ZOS_change_ensemble(i,j,:));
        x = AMOC_change_ensemble(:);  % Ensure column vector

        if all(isfinite(x)) && all(isfinite(y))
            [r, p] = corr(x, y, 'Rows', 'complete');  % Pearson by default
            ensemble_corr(i,j) = r;
            p_values_ensemble(i,j) = p;
        else
            ensemble_corr(i,j) = NaN;
            p_values_ensemble(i,j) = NaN;
        end
    end
end


% Find significant regions (p ≤ 0.05)
sig_mask_ensemble = p_values_ensemble <= 0.05;

% Plot the ensemble correlation at tile 1
nexttile(1);
pcolor(lon(lon_range), lat(lat_range), ensemble_corr'); 
shading flat;
hold on;

contour(lon(lon_range), lat(lat_range), ensemble_corr', [-1:0.2:1],'-k', 'LineWidth', 1.0);

% Overlay continental boundaries
load('na.dat'); 
plot(na(:,1) + 360, na(:,2), 'k', 'LineWidth', 1.5);
% Set axis limits
xlim([270 360]); 
ylim([0 60]);

% Overlay significant points
[x_sig, y_sig] = meshgrid(lon(lon_range), lat(lat_range));
scatter(x_sig(sig_mask_ensemble'), y_sig(sig_mask_ensemble'), 5, 'k', 'filled');

% Set up colormap
load('newcolormaps', 'bwr'); 
colormap(bwr);
caxis([-0.5 0.5]);

% Customize axes
set(gca, 'FontSize', 13, 'YTick', 0:20:80, 'YTickLabel', {'0', '20°N', '40°N', '60°N', '80°N'});
% set(gca, 'FontSize', 13, 'XTick', 270:20:360, 'XTickLabel', {'90°W', '70°W', '50°W', '30°W', '10°W'});
set(gca, 'FontSize', 13, 'XTick', []);

title('Ensemble-Mean', 'FontWeight', 'bold');

cmap = makeColorMap([0 0 1],[1 1 1],[1  0 0],10);
colorbar
caxis([-0.5 0.5]); colormap(cmap)
cb = colorbar('Location', 'eastoutside');
cb.Layout.Tile = 'east';  % Puts it beside the tiled layout
cb.Label.String = 'Correlation coefficient';
cb.FontSize = 13;

