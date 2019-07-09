clear; clc

% The lake detection algorithm for ICESat-2 (first iteration)
% Currently, only the center strong beam (gt2r) is considered
% Subroutines used: is2_class_merge, is2_windowing_sub, is2_chunking_sub,
% is2_sfc_detect_sub, is2_sig_count_sub


% Call data - needs to be updated once larger data volumes are available
fname_atll03 = 'ATL03_20190102184312_00810210_001_01.h5';
try
    fname_atl06 = ['ATL06', fname_atll03(6:end)];
catch
    error(['Equivalent ATL06 file not found for: ', fname_atll03])
end

% Beam Data (csb = Central strong beam, cwb = Central weak beam)
time_csb = h5read(fname_atll03, '/gt2l/heights/delta_time');
elev_csb = h5read(fname_atll03, '/gt2l/heights/h_ph');
lat_csb = h5read(fname_atll03, '/gt2l/heights/lat_ph');
lon_csb = h5read(fname_atll03, '/gt2l/heights/lon_ph');
dist = h5read(fname_atll03, '/gt2l/heights/dist_ph_along'); % In meters
seglen = h5read(fname_atll03, '/gt2l/geolocation/segment_length'); % Only used to properly convert dist

time06 = h5read(fname_atl06, '/gt2l/land_ice_segments/delta_time');
elev06 = h5read(fname_atl06, '/gt2l/land_ice_segments/h_li');
lat06 = h5read(fname_atl06, '/gt2l/land_ice_segments/latitude');
lon06 = h5read(fname_atl06, '/gt2l/land_ice_segments/longitude');
snr06 = h5read(fname_atl06, '/gt2l/land_ice_segments/fit_statistics/snr');
dist06 = h5read(fname_atl06, '/gt2l/land_ice_segments/ground_track/x_atc');

% Removing surface type classification, filler values
class_csb = h5read(fname_atll03, '/gt2l/heights/signal_conf_ph');
class_consol_csb = is2_class_merge(class_csb);
snr06(snr06 > 1e6) = NaN;
snr06 = 10*log(snr06); % Convert to dBZ
elev06(elev06 > 1e6) = NaN;

% Converting time, distance relative to start of track
for j = 1:length(time_csb)
    delta_time_csb(j) = time_csb(j) - time_csb(1);
end

total_dist = sum(seglen) + dist(end);
dist_corrected = linspace(0,total_dist, length(dist))./1000; % Also convert to km

% Separating the data into multiple windows; will need to be updated
% according to the length of each data granule
window_size_csb = floor(length(elev_csb)/25);
window_count_csb = ceil(length(elev_csb)/window_size_csb);

window_size06 = floor(length(elev06)/25);
window_count06 = ceil(length(elev06)/window_size_csb);

%% Lake Surface-Bed Separation
for j = 1:window_count_csb
    % Windowing Subroutines
    [lat_bin_csb,lon_bin_csb,time_bin_csb,elev_bin_csb,class_bin_csb,dist_bin] = is2_windowing_sub(lat_csb,lon_csb,delta_time_csb,elev_csb,class_consol_csb,dist_corrected,window_size_csb,j);
    [lat_bin06,lon_bin06,time_bin06,elev_bin06,snr_bin06,dist_bin06] = atl06_windowing_sub(lat06,lon06,time06,elev06,snr06,dist06,window_size06,j);
    high_bin_csb = elev_bin_csb;
    high_bin_csb(class_bin_csb~=4) = NaN;
    
    % Surface-/Bed-Finding Subroutines 
    window_lake_sfc = is2_sfc_detect_sub(high_bin_csb,class_bin_csb,elev_bin06,snr_bin06);
    [~,window_lake_btm,lake_btm_mean,lake_btm_error] = is2_ph_dist(dist_bin,high_bin_csb);
    
    % Lake Bounds
    lake_start = find(~isnan(window_lake_sfc), 1, 'first');
    lake_end = find(~isnan(window_lake_sfc), 1, 'last');
    
    % Boundary conditions
    window_lake_btm(lake_start) = window_lake_sfc(lake_start); 
    window_lake_btm(lake_end) = window_lake_sfc(lake_end);
    lat_lake_start03 = lat_bin_csb(lake_start);
    lat_lake_end03 = lat_bin_csb(lake_end);
    
    if ~isempty(lake_start)
        lake_start06 = find(round(lat06,3) == round(lat_lake_start03,3),1,'first');
        lake_end06 = find(round(lat06,3) == round(lat_lake_end03,3),1,'last');
        disp(['Mean SNR over ''lake'': ', num2str(nanmean(snr06(lake_start06:lake_end06)))]);
        lake_range = lake_start:lake_end;
    else
        lake_range = -9999; % Filler value to stop polyfitting
    end
    
    % Filtering Statistics
    lake_sfc_mean = movmean(window_lake_sfc, 1600, 'omitnan');
    high_bin_mean = movmean(high_bin_csb, 3000, 'omitnan');
    high_bin_std = movstd(high_bin_csb, 10000, 'omitnan');
    number_of_btm_nans = movsum(~isnan(lake_btm_mean),3000,'omitnan');
    
    for k = 1:length(window_lake_sfc)
        if k < lake_start % Remove points before start of lake
            lake_sfc_mean(k) = NaN;
            lake_btm_mean(k) = NaN;
        elseif k > lake_end % Remove points after end of lake
            lake_sfc_mean(k) = NaN;
            lake_btm_mean(k) = NaN;
        elseif high_bin_mean(k) > lake_sfc_mean(k) % Remove points above lake
            lake_sfc_mean(k) = NaN;
            lake_btm_mean(k) = NaN;
        elseif lake_sfc_mean(k)-lake_btm_mean(k) < 1 % Remove shallow points
            lake_sfc_mean(k) = NaN;
            lake_btm_mean(k) = NaN;
        elseif number_of_btm_nans(k) < 1800 % Remove small false positives
            lake_btm_mean(k) = NaN;
            lake_sfc_mean(k) = NaN;
        elseif high_bin_std(k) >= 2.1 % Remove points on rough surfaces
            lake_btm_mean(k) = NaN;
            lake_sfc_mean(k) = NaN;
        end
    end
    lake_btm_mean(lake_start) = lake_sfc_mean(lake_start);
    lake_btm_mean(lake_end) = lake_sfc_mean(lake_end);
    
    % Filtering and the start and end points
    if lake_btm_mean(lake_start+1)-lake_btm_mean(lake_start) < -1
        lake_btm_mean(lake_start) = NaN;
    elseif lake_btm_mean(lake_end)-lake_btm_mean(lake_end-1) > 1
        lake_btm_mean(lake_end) = NaN;
    end
    
    % Polynomial Fitting, with lake endpoints as boundary
    % conditions
    lake_btm_fitted = is2_polyfit_sub(lon_bin_csb,window_lake_btm);

    % Filtering the polyfit curves to look nicer
    lake_btm_fitted(isnan(window_lake_sfc)) = NaN;
    
    % Lake depth estimates
    sfc_btm_diff_polyfit = lake_sfc_mean - lake_btm_fitted;
    lake_depth_polyfit = max(sfc_btm_diff_polyfit(lake_start:lake_end));
    sfc_btm_diff_means = lake_sfc_mean - lake_btm_mean;
    lake_depth_means = max(sfc_btm_diff_means(lake_start:lake_end));
    mean_depth = mean([lake_depth_means lake_depth_polyfit]);
    
    % Depth markers for plots
    try
        if lake_depth_means > lake_depth_polyfit
            marked_int = find(sfc_btm_diff_means==lake_depth_means, 1, 'first');
            depth_marker = lake_btm_mean(marked_int);
        elseif lake_depth_polyfit > lake_depth_means
            marked_int = find(sfc_btm_diff_polyfit==lake_depth_polyfit, 1, 'first');
            depth_marker = lake_btm_fitted(marked_int);
        end
    catch
        warning('No lake detected for this data window.')
        lake_depth_polyfit = NaN;
    end
    
    % Plotting
    if any(lake_sfc_mean)
        figure;
        plot(dist_bin, high_bin_csb, '.', 'MarkerSize', 2, 'Color', rgb('sky blue'))
        hold on; plot(dist_bin, lake_sfc_mean, '.', 'MarkerSize', 2, 'Color', rgb('green'))
        plot(dist_bin, lake_btm_fitted, 'LineWidth', 2, 'Color', rgb('rose'))
        plot(dist_bin, lake_btm_mean, 'LineWidth', 2)
        plot(dist_bin(marked_int), depth_marker, 'r*', 'MarkerSize',12)
        xlabel('Along-track distance [km]', 'FontSize',14, 'FontWeight','bold');
        ylabel('Elevation [m]', 'FontSize',14, 'FontWeight','bold')
        legend('Raw', 'Signal Surface', 'Polyfit Bed', 'Signal Bed', 'Max Depth')
        pause; close all;
    end
    
end
   
