clear; clc;

% Define the common grid
lat_min = -20; lat_max = 60;      % Latitude range
lev_min = 1000; lev_max = 4000;   % Depth range
lat_new = linspace(lat_min, lat_max, 300);  % Uniform latitude grid (300 points)
lev_new = linspace(lev_min, lev_max, 50);   % Uniform depth grid (50 points)

% List of model filenames
model_names_amoc = { ...
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

for m = 1:length(model_names_amoc)
    filename = model_names_amoc{m};
    
    % Display progress
    disp(['Processing ', filename, '...']);
    
    % Read lat, lev, time, and variable
    lat_old = ncread(filename, 'vlat');   % Adjust if needed
    lev_old = ncread(filename, 'lev');   % Adjust if needed
    time_old = ncread(filename, 'time'); % Time remains unchanged
    data_old = ncread(filename, 'moc_section');  % Adjust variable name if needed
    
    % Extract within desired range
    lat_mask = (lat_old >= lat_min & lat_old <= lat_max);
    lev_mask = (lev_old >= lev_min & lev_old <= lev_max);
    lat_old = lat_old(lat_mask);
    lev_old = lev_old(lev_mask);
    data_old = data_old(lat_mask, lev_mask, :);  % Subset data
    
    % Initialize new regridded data array
    data_new = NaN(length(lat_new), length(lev_new), length(time_old));
    
    % Loop over time dimension to interpolate data
    for t = 1:length(time_old)
        temp_data = interp1(lat_old, squeeze(data_old(:, :, t)), lat_new, 'linear', 'extrap');
        for j = 1:length(lat_new)
            data_new(j, :, t) = interp1(lev_old, squeeze(temp_data(j, :)), lev_new, 'linear', 'extrap');
        end
    end
    
    % Define new filename based on original
    new_filename = ['Harmonized_', filename];
    
    % Create new NetCDF file with harmonized dimensions
    nccreate(new_filename, 'lat', 'Dimensions', {'lat', length(lat_new)});
    nccreate(new_filename, 'lev', 'Dimensions', {'lev', length(lev_new)});
    nccreate(new_filename, 'time', 'Dimensions', {'time', length(time_old)});
    nccreate(new_filename, 'moc_section', 'Dimensions', {'lat', length(lat_new), 'lev', length(lev_new), 'time', length(time_old)});

    % Write new data
    ncwrite(new_filename, 'lat', lat_new);
    ncwrite(new_filename, 'lev', lev_new);
    ncwrite(new_filename, 'time', time_old);
    ncwrite(new_filename, 'moc_section', data_new);

    disp(['Saved harmonized data: ', new_filename]);
end

disp('Harmonization complete for all models.');
