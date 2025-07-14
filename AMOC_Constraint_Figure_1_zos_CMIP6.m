clc; clear;
%% -------------------------------------------------------------------------
% Description:
% Computes dynamic sea level change between two periods for historical runs
% across multiple CMIP6 models, performs a t-test across models, and plots
% significant ensemble mean changes in the North Atlantic.
% -------------------------------------------------------------------------

%% Load model file list
zos_files = {
    'zos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_1850-2014.nc',
    % ... add all other model files here ...
    'zos_Omon_NorESM2-LM_historical_r1i1p1f1_gn_1850-2014.nc'
};

nmodels = numel(zos_files);

% Read latitude and longitude from first file
lat_mod = ncread(zos_files{1}, 'lat');
lon_mod = ncread(zos_files{1}, 'lon');
lon_mod(lon_mod > 180) = lon_mod(lon_mod > 180) - 360;  % Convert to -180 to 180
[lon_mod, loni_mod] = sort(lon_mod);  % Sort for consistency

nlat = length(lat_mod);
nlon = length(lon_mod);

% Preallocate
zos_change_models = nan(nlon, nlat, nmodels);

% Define time periods
yr_mod = 1850:2014;
yr1_mod = [1850, 1879];  % Early period
yr2_mod = [1985, 2014];  % Late period

% Process each model
for i = 1:nmodels
    % Read ZOS data
    data = ncread(zos_files{i}, 'zos');
    data = data(loni_mod, :, :);  % Align longitudes
    
    % Reshape to (lon x lat x year) and compute annual means
    data_yr = squeeze(mean(reshape(data, size(data,1), size(data,2), 12, []), 3));

    % Calculate mean ZOS for early and late periods
    mean1 = mean(data_yr(:, :, yr_mod >= yr1_mod(1) & yr_mod <= yr1_mod(2)), 3, 'omitnan');
    mean2 = mean(data_yr(:, :, yr_mod >= yr2_mod(1) & yr_mod <= yr2_mod(2)), 3, 'omitnan');

    % Compute dynamic sea level change and remove global mean
    change = mean2 - mean1;
    change = change - mean(change, 'all', 'omitnan');
    zos_change_models(:, :, i) = change;
end

% Ensemble Mean and Significance Test
ensemble_mean = mean(zos_change_models, 3, 'omitnan');

% One-sample t-test against zero
[h, pval] = ttest(zos_change_models, 0, 'dim', 3, 'alpha', 0.05);  % h = 1 if significant

% Define North Atlantic Region
lat_inds = find(lat_mod >= 0 & lat_mod <= 60);
lon_inds = find(lon_mod >= -90 & lon_mod <= 0);

% Subset to North Atlantic
zos_na = zos_change_models(lon_inds, lat_inds, :);
ensemble_mean_na = mean(zos_na, 3, 'omitnan');

% Perform t-test in North Atlantic
[h_na, pval_na] = ttest(zos_na, 0, 'dim', 3, 'alpha', 0.05);

% Grids for plotting
lon_na = lon_mod(lon_inds);
lat_na = lat_mod(lat_inds);
[lon_grid_na, lat_grid_na] = meshgrid(lon_na, lat_na);

% Plot Results
figure;
m_proj('robinson', 'lon', [-90, 0], 'lat', [0, 60]);
m_pcolor(lon_na, lat_na, ensemble_mean_na'); shading interp;

m_coast('line', 'color', [0,0,0]);
m_grid('tickdir', 'in', 'linewi', 2);

% Colormap
cmap = makeColorMap([0 0 1], [1 1 1], [1 0 0], 10);
colormap(cmap);
clim([-0.05 0.05]);
a = colorbar;
ylabel(a, '(m)', 'FontSize', 14, 'Rotation', 270);
title('Dynamic Sea Level Change (95% Confidence)');

% Overlay statistically significant grid points
hold on;
sig_mask_na = h_na';  % transpose for plotting
step = 1;  % change to 2 or 3 for fewer points
lon_sub = lon_grid_na(1:step:end, 1:step:end);
lat_sub = lat_grid_na(1:step:end, 1:step:end);
sig_sub = sig_mask_na(1:step:end, 1:step:end);

[row, col] = find(sig_sub == 1);
lon_sig = lon_sub(sub2ind(size(sig_sub), row, col));
lat_sig = lat_sub(sub2ind(size(sig_sub), row, col));

[x, y] = m_ll2xy(lon_sig, lat_sig);
scatter(x, y, 10, 'k', 'filled');
