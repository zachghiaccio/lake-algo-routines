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
time_csb = h5read(fname_atll03, '/gt1l/heights/delta_time');
elev_csb = h5read(fname_atll03, '/gt1l/heights/h_ph');
lat_csb = h5read(fname_atll03, '/gt1l/heights/lat_ph');
lon_csb = h5read(fname_atll03, '/gt1l/heights/lon_ph');

time06 = h5read(fname_atl06, '/gt1l/land_ice_segments/delta_time');
elev06 = h5read(fname_atl06, '/gt1l/land_ice_segments/h_li');
lat06 = h5read(fname_atl06, '/gt1l/land_ice_segments/latitude');
lon06 = h5read(fname_atl06, '/gt1l/land_ice_segments/longitude');
snr06 = h5read(fname_atl06, '/gt1l/land_ice_segments/fit_statistics/snr');
dist06 = h5read(fname_atl06, '/gt1l/land_ice_segments/ground_track/x_atc');

% Removing surface type classification, filler values
class_csb = h5read(fname_atll03, '/gt1l/heights/signal_conf_ph');
class_consol_csb = is2_class_merge(class_csb);
snr06(snr06 > 1e6) = NaN;
snr06 = 10*log(snr06); % Convert to dBZ
elev06(elev06 > 1e6) = NaN;

dist = ((1:length(elev_csb)).*0.7)./1000; % ATL03 Along track distance in km

% Converting time relative to start of track
for j = 1:length(time_csb)
    delta_time_csb(j) = time_csb(j) - time_csb(1);
end

% Separating the data into multiple windows; will need to be updated
% according to the length of each data granule
window_size_csb = floor(length(elev_csb)/25);
window_count_csb = ceil(length(elev_csb)/window_size_csb);

window_size06 = floor(length(elev06)/25);
window_count06 = ceil(length(elev06)/window_size_csb);

for j = 1:window_count_csb
    [lat_bin_csb,lon_bin_csb,time_bin_csb,elev_bin_csb,class_bin_csb,dist_bin] = is2_windowing_sub(lat_csb,lon_csb,delta_time_csb,elev_csb,class_consol_csb,dist,window_size_csb,j);
    [lat_bin06,lon_bin06,time_bin06,elev_bin06,snr_bin06,dist_bin06] = atl06_windowing_sub(lat06,lon06,time06,elev06,snr06,dist06,window_size06,j);
    

    high_bin_csb = elev_bin_csb;
    high_bin_csb(class_bin_csb~=4) = NaN;
    [window_lake_sfc, window_lake_btm] = is2_sfc_detect_sub(high_bin_csb,class_bin_csb,elev_bin06,snr_bin06);
    lake_start = find(~isnan(window_lake_sfc), 1, 'first');
    lake_end = find(~isnan(window_lake_sfc), 1, 'last');
    
    
    window_lake_btm(lake_start) = window_lake_sfc(lake_start); % Boundary conditions
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
    
    lake_sfc_mean = movmean(window_lake_sfc, 1600, 'omitnan');
    lake_btm_max = movmax(window_lake_btm, 1600, 'omitnan');
    high_bin_mean = movmean(high_bin_csb, 3000, 'omitnan');
    high_bin_std = movstd(high_bin_csb, 10000, 'omitnan');
    number_of_btm_nans = movsum(~isnan(lake_btm_max),3000,'omitnan');
    
    for k = 1:length(window_lake_sfc)
        if k < lake_start
            %lake_btm_fitted(k) = NaN;
            lake_sfc_mean(k) = NaN;
            lake_btm_max(k) = NaN;
        elseif k > lake_end
            %lake_btm_fitted(k) = NaN;
            lake_sfc_mean(k) = NaN;
            lake_btm_max(k) = NaN;
        elseif high_bin_mean(k) > lake_sfc_mean(k)
            %lake_btm_fitted(k) = NaN;
            lake_sfc_mean(k) = NaN;
            lake_btm_max(k) = NaN;
        elseif lake_sfc_mean(k)-lake_btm_max(k) < 1
            %lake_btm_fitted(k) = NaN;
            lake_sfc_mean(k) = NaN;
            lake_btm_max(k) = NaN;
        elseif number_of_btm_nans(k) < 1800
            %lake_btm_fitted(k) = NaN;
            lake_btm_max(k) = NaN;
            lake_sfc_mean(k) = NaN;
        elseif high_bin_std(k) >= 2.1
            %lake_btm_fitted(k) = NaN;
            lake_btm_max(k) = NaN;
            lake_sfc_mean(k) = NaN;
        end
    end
    lake_btm_max(lake_start) = lake_sfc_mean(lake_start);
    lake_btm_max(lake_end) = lake_sfc_mean(lake_end);
    
    % Further filtering, for cleanliness
    if lake_btm_max(lake_start+1)-lake_btm_max(lake_start) < -1
        lake_btm_max(lake_start) = NaN;
    elseif lake_btm_max(lake_end)-lake_btm_max(lake_end-1) > 1
        lake_btm_max(lake_end) = NaN;
    end
    
    % Polynomial Fitting, with lake endpoints as boundary
    % conditions
    lake_btm_fitted = is2_polyfit_sub(lon_bin_csb,window_lake_btm);
    
%     if lake_start == 1 % If lake is at beginning of data file
%         idx = isnan(lake_btm_max);
%         [p,~,mu] = polyfit(lon_bin_csb(~idx), window_lake_btm(~idx), 3);
%         lake_btm_fitted = polyval(p, lon_bin_csb, [], mu);
%     elseif any(length(high_bin_csb)-10:length(high_bin_csb) == lake_end) % If lake is at end of data file
%         idx = isnan(window_lake_btm);
%         [p,~,mu] = polyfit(lon_bin_csb(~idx), window_lake_btm(~idx), 3);
%         lake_btm_fitted = polyval(p, lon_bin_csb, [], mu);
%     else
%         idx = isnan(window_lake_btm);
%         p = polyfix(lon_bin_csb(~idx),window_lake_btm(~idx), 3, ...
%             [lon_bin_csb(lake_start) lon_bin_csb(lake_end)], ...
%             [window_lake_sfc(lake_start) window_lake_sfc(lake_end)]);
%         lake_btm_fitted = polyval(p, lon_bin_csb);
%     end
    
%     for k = 1:length(lake_btm_fitted)
%         if k<lake_start | k>lake_end
%             lake_btm_fitted(k) = NaN;     
%         end
%     end

    % Filtering the polyfit curves to look nicer
    lake_btm_fitted(isnan(window_lake_sfc)) = NaN;
    
    sfc_btm_diff_polyfit = lake_sfc_mean - lake_btm_fitted;
    lake_depth_polyfit = max(sfc_btm_diff_polyfit(lake_start:lake_end));
    sfc_btm_diff_means = lake_sfc_mean - lake_btm_max;
    lake_depth_means = max(sfc_btm_diff_means(lake_start:lake_end));
    mean_depth = mean([lake_depth_means lake_depth_polyfit]);
    
    
    try
        if lake_depth_means > lake_depth_polyfit
            marked_int = find(sfc_btm_diff_means==lake_depth_means, 1, 'first');
            depth_marker = lake_btm_max(marked_int);
        elseif lake_depth_polyfit > lake_depth_means
            marked_int = find(sfc_btm_diff_polyfit==lake_depth_polyfit, 1, 'first');
            depth_marker = lake_btm_fitted(marked_int);
        end
    catch
        warning('No lake detected for this data window.')
        lake_depth_polyfit = NaN;
    end
    
    if any(lake_sfc_mean)
        figure;
        plot(lon_bin_csb, high_bin_csb, '.', 'MarkerSize', 2, 'Color', rgb('sky blue'))
        hold on; plot(lon_bin_csb, lake_sfc_mean, '.', 'MarkerSize', 2, 'Color', rgb('green'))
        plot(lon_bin_csb, lake_btm_fitted, 'LineWidth', 2, 'Color', rgb('rose'))
        plot(lon_bin_csb, lake_btm_max, '.', 'MarkerSize', 2)
        plot(lon_bin_csb(marked_int), depth_marker, 'r*', 'MarkerSize',12)
        xlabel('Longitude', 'FontSize',14, 'FontWeight','bold');
        ylabel('Elevation [m]', 'FontSize',14, 'FontWeight','bold')
        legend('Raw', 'Polyfit Bed', 'Signal Surface', 'Signal Bed', 'Max Depth')
        %keyboard
        pause; close all;
    end
    
    
%     elev_signal_high_csb = elev_bin_csb; elev_noise_csb = elev_bin_csb;
%     elev_signal_high_csb(class_bin_csb~=4) = NaN;
%     elev_noise_csb(class_bin_csb~=0) = NaN;
%     
%     % SNR Computation
%     time_bin_csb_interval = time_bin_csb(end)/420;
%     for k = 1:420
%         signal_high_count_csb(k) = length(elev_signal_high_csb(~isnan(elev_signal_high_csb) & time_bin_csb'<time_bin_csb_interval*k & time_bin_csb'>=time_bin_csb_interval*(k-1)));
%         noise_count_csb(k) = length(elev_noise_csb(~isnan(elev_noise_csb) & time_bin_csb'<time_bin_csb_interval*k & time_bin_csb'>=time_bin_csb_interval*(k-1)));
%     end
%     
%     figure;
%     subplot(1,2,1)
%     plot(signal_high_count_csb, 'LineWidth', 2)
%     hold on; plot(noise_count_csb, 'LineWidth', 2)
%     subplot(1,2,2)
%     histogram(elev_signal_high_csb)
%     
%     figure;
%     plot(elev_signal_high_csb, '.', 'MarkerSize', 2)
%     hold on; plot(elev_noise_csb, 'r.', 'MarkerSize', 2)
%     pause; close all
     clear('signal_high_count_csb', 'noise_count_csb')
end
   
