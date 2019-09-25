clear; clc

% The lake detection algorithm for ICESat-2 (first iteration)
% Currently, only the center strong beam (gt2r) is considered
% Subroutines used: is2_class_merge, is2_windowing_sub, is2_chunking_sub,
% is2_sfc_detect_sub, is2_sig_count_sub

% ----------------------------------------------------------------------- %

% Call data - needs to be updated once larger data volumes are available
fname_atll03 = 'ATL03_20190102184312_00810210_001_01.h5';

% Plotting Booleans
indyplot = 0;
multiplot = 1;

% Number of windows per granule
window_number = 25;

% ----------------------------------------------------------------------- %

% Beam Data
time = h5read(fname_atll03, '/gt2l/heights/delta_time');
elev = h5read(fname_atll03, '/gt2l/heights/h_ph');
lat = h5read(fname_atll03, '/gt2l/heights/lat_ph');
lon = h5read(fname_atll03, '/gt2l/heights/lon_ph');
dist = h5read(fname_atll03, '/gt2l/heights/dist_ph_along'); % In meters
seglen = h5read(fname_atll03, '/gt2l/geolocation/segment_length'); % Only used to properly convert dist
theta_elev = h5read(fname_atll03, '/gt2l/geolocation/ref_elev');
theta_elev(theta_elev>2*pi) = NaN;
theta_elev_mean = nanmean(theta_elev); % Elevation angle change is small
theta_inc = pi/2 - theta_elev_mean; % Incidence angle, for refraction correction
class = h5read(fname_atll03, '/gt2l/heights/signal_conf_ph');
class_consol = is2_class_merge(class); % Removing surface type classification
total_dist = sum(seglen) + dist(end);
dist_corrected = linspace(0,total_dist, length(dist))./1000; % Also convert to km

% Separating the data into multiple windows (Currently 
window_size = floor(length(elev)/window_number);
window_count = ceil(length(elev)/window_size);
lake_count = 0;

% Lake Surface-Bed Separation
for j = 1:window_count
    %% Windowing Subroutines
    [lat_bin,lon_bin,time_bin,elev_bin,class_bin,dist_bin] = is2_windowing_sub(lat,lon,delta_time,elev,class_consol,dist_corrected,window_size,j);
    high_bin = elev_bin; sig_bin = elev_bin;
    high_bin(class_bin~=4) = NaN; % High-confidence photons only
    sig_bin(class_bin<1) = NaN; % High/medium/low/buffer photons
    
    %% Surface-/Bed-Finding Subroutines 
    window_lake_sfc = is2_sfc_detect_sub(high_bin,class_bin);
    [~,window_lake_btm,lake_btm_mean,lake_btm_error] = is2_ph_dist(dist_bin,high_bin);
    
    % Lake Bounds
    lake_start = find(~isnan(window_lake_sfc), 1, 'first');
    lake_end = find(~isnan(window_lake_sfc), 1, 'last');
    if ~isempty(lake_start)
        lake_range = lake_start:lake_end;
    else
        lake_range = -9999; % Filler value to stop polyfitting
    end
    
    % Boundary conditions
    window_lake_btm(lake_start) = window_lake_sfc(lake_start); 
    window_lake_btm(lake_end) = window_lake_sfc(lake_end);
    
    % Filtering Statistics
    lake_sfc_mean = movmean(window_lake_sfc, 1600, 'omitnan');
    high_bin_mean = movmean(high_bin, 3000, 'omitnan');
    high_bin_std = movstd(high_bin, 10000, 'omitnan');
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
    
    % Filtering of the start and end points
    if lake_btm_mean(lake_start+1)-lake_btm_mean(lake_start) < -1
        lake_btm_mean(lake_start) = NaN;
    elseif lake_btm_mean(lake_end)-lake_btm_mean(lake_end-1) > 1
        lake_btm_mean(lake_end) = NaN;
    end
    
    %% Lake Surface/Bed Corrections
    [window_lake_sfc_corr,lake_sfc_mean_corr,~,lake_sfc_slope] = is2_sfc_find(dist_bin,high_bin,window_lake_sfc);
    
    % Lake depth estimates
    sfc_btm_diff_means = lake_sfc_mean - lake_btm_mean;
    lake_depth_means = max(sfc_btm_diff_means(lake_start:lake_end));
    
    % Refraction Correction
    [dist_bin_corr,signal_depth_corr,dz] = depth_refrac_fix(dist_bin,sfc_btm_diff_means,theta_inc);
    window_lake_btm_corr = window_lake_btm + dz;
    idx = (high_bin == window_lake_btm);
    high_bin(idx) = high_bin(idx) + dz(idx);
    
    lake_btm_mean_corr = movmean(window_lake_btm_corr,3000,'omitnan');
    for i = 2:length(lake_btm_mean_corr)
        if abs(lake_btm_mean_corr(i)-lake_btm_mean_corr(i-1))>0.2
            lake_btm_mean_corr(i) = NaN;
        end
    end
    
    
    %% Polynomial Fitting
    lake_btm_fitted = is2_polyfit_sub(lon_bin,window_lake_btm);

    % Polynomial Filtering
    lake_btm_fitted(isnan(window_lake_sfc)) = NaN;
    lake_btm_fitted = lake_btm_fitted + dz;
    lake_btm_fitted(lake_btm_fitted>lake_sfc_mean) = NaN;
    sfc_btm_diff_polyfit = lake_sfc_mean - lake_btm_fitted;
    lake_depth_polyfit = max(sfc_btm_diff_polyfit(lake_start:lake_end)); % Fitted depth
    
    %% Error Statistics
    sigma_btm = movstd(window_lake_btm_corr,100,'omitnan'); % Bottom uncertainty
    idx = isnan(window_lake_btm_corr);
    resid = window_lake_btm_corr - lake_btm_fitted;
    delta = sqrt(nansum(resid(~idx).^2)/length(window_lake_btm(~idx))); % Polynomial standard error
    
    %% Plotting
    if any(lake_sfc_mean)
        lake_count = lake_count + 1;
        disp(lake_count)
        
        if indyplot
            figure;
            plot(dist_bin, high_bin, '.', 'MarkerSize', 6, 'Color', rgb('sky blue'))
            hold on; plot(dist_bin, lake_sfc_mean_corr, 'LineWidth', 2, 'Color', rgb('green'))
            plot(dist_bin, lake_btm_fitted, 'LineWidth', 2, 'Color', rgb('rose'))
            plot(dist_bin_corr, lake_btm_mean_corr, 'b', 'LineWidth', 2)
            plot(dist_bin(marked_int), depth_marker, 'r*', 'MarkerSize',12)
            xlabel('Along-track distance [km]', 'FontSize',14, 'FontWeight','bold');
            ylabel('Elevation [m]', 'FontSize',14, 'FontWeight','bold')
            legend('Raw', 'Signal Surface', 'Polyfit Bed', 'Signal Bed', 'Max Depth')
            pause; close all;
        elseif multiplot
            
            if lake_count == 2
                figure;
                x = 445797:473763;
                marked_int = 461876; depth_marker = lake_btm_mean_corr(marked_int);
                subplot(2,2,1)
                plot(dist_bin(x), high_bin(x), '.', 'MarkerSize', 3, 'color', rgb('sky blue'))
                hold on; plot(dist_bin(x),lake_sfc_mean_corr(x),'LineWidth',2,'color',rgb('green'))
                plot(dist_bin(x), lake_btm_fitted(x), 'LineWidth', 2, 'Color', rgb('rose'))
                plot(dist_bin_corr(x), lake_btm_mean_corr(x), 'b', 'LineWidth', 2)
                plot(real(dist_bin(marked_int)), depth_marker, 'r*', 'MarkerSize',12)
                xlim([min(dist_bin(x)) max(dist_bin(x))])
                ylim([min(high_bin(x))-2 max(high_bin(x))])
                
                marked_int = find(lake_btm_fitted==min(lake_btm_fitted(x),[],'omitnan'));
                keyboard;
            elseif lake_count == 4
                subplot(2,2,2)
                x = 131034:203691;
                marked_int = 188177; depth_marker = lake_btm_mean_corr(marked_int);
                plot(dist_bin(x), high_bin(x), '.', 'MarkerSize', 3, 'color', rgb('sky blue'))
                hold on; plot(dist_bin(x),lake_sfc_mean_corr(x),'LineWidth',2,'color',rgb('green'))
                plot(dist_bin(x), lake_btm_fitted(x), 'LineWidth', 2, 'Color', rgb('rose'))
                plot(dist_bin_corr(x), lake_btm_mean_corr(x), 'b', 'LineWidth', 2)
                plot(real(dist_bin(marked_int)), depth_marker, 'r*', 'MarkerSize',12)
                xlim([min(dist_bin(x)) max(dist_bin(x))])
                ylim([min(high_bin(x))-2 max(high_bin(x))])
                
                marked_int = find(lake_btm_fitted==min(lake_btm_fitted(x),[],'omitnan'));
                keyboard;
                
                subplot(2,2,3)
                x = 251643:291204;
                lake_sfc_mean_corr(251643:271258) = NaN;
                lake_btm_mean(251643:271258) = NaN;
                marked_int = 279070; depth_marker = lake_btm_mean_corr(marked_int);
                plot(dist_bin(x), high_bin(x), '.', 'MarkerSize', 3, 'color', rgb('sky blue'))
                hold on; plot(dist_bin(x),lake_sfc_mean_corr(x),'LineWidth',2,'color',rgb('green'))
                plot(dist_bin(x), lake_btm_fitted(x), 'LineWidth', 2, 'Color', rgb('rose'))
                plot(dist_bin_corr(x), lake_btm_mean_corr(x), 'b', 'LineWidth', 2)
                plot(real(dist_bin(marked_int)), depth_marker, 'r*', 'MarkerSize',12)
                xlim([min(dist_bin(x)) max(dist_bin(x))])
                ylim([min(high_bin(x))-2 max(high_bin(x))])
                
                marked_int = find(lake_btm_fitted==min(lake_btm_fitted(x),[],'omitnan'));
                keyboard;
                
                subplot(2,2,4)
                x = 326860:371236;
                lake_sfc_mean_corr(348647:371236) = NaN;
                lake_btm_mean(348647:371236) = NaN;
                marked_int = 333702; depth_marker = lake_btm_mean_corr(marked_int);
                plot(dist_bin(x), high_bin(x), '.', 'MarkerSize', 3, 'color', rgb('sky blue'))
                hold on; plot(dist_bin(x),lake_sfc_mean_corr(x),'LineWidth',2,'color',rgb('green'))
                plot(dist_bin(x), lake_btm_fitted(x), 'LineWidth', 2, 'Color', rgb('rose'))
                plot(dist_bin_corr(x), lake_btm_mean_corr(x), 'b', 'LineWidth', 2)
                plot(real(dist_bin(marked_int)), depth_marker, 'r*', 'MarkerSize',12)
                xlim([min(dist_bin(x)) max(dist_bin(x))])
                ylim([min(high_bin(x))-2 max(high_bin(x))])
                
                marked_int = find(lake_btm_fitted==min(lake_btm_fitted(x),[],'omitnan'));
                keyboard;
                pause;
            end
        end
    end
    
end
   
