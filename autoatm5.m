clear; clc; tic;
% An update to autoatm4.m, improving the polynomial interpolation scheme
% and making the process more automated.
% The struct "lake" is produced at the end, giving the latitudes,
% longitudes, and depths of all detected supraglacial lakes.

addpath(genpath('/Users/zhfair/Documents/'))


%-------------------------------------------------------------------------%
%% User Input

% List of data files; change the date to view different data sets
fnames = importdata('fnames_20170719.txt');
if ~strcmp(fnames{1}(3:5), 'ATM')
    error('Input file is not from ATM.')
end

% for i = 1:length(fnames)
%     disp([num2str(i), ': ', fnames{i}])
% end
% disp(' ')
% stamp = input(['Which file do you want to look at? Choose 1-', num2str(length(fnames)), ': ']);
% fname = fnames{stamp};

% Number of photons needed in an elevation bin to be classified as a lake
lake_thresh = 12000;

indyplot = 1; % Display plots = 1, no plots = 0
multiplot = 0;

%-------------------------------------------------------------------------%
lat_list = []; lon_list = [];
mean_depth_list = []; fitted_depth_list = [];
lake_vs_list = []; lake_vp_list = [];
lake_radius_list = [];
lake_count = 0;
for k = 1:length(fnames(1:9))
    fname = fnames{k};
    time = h5read(fname, '/instrument_parameters/rel_time');
    for j = 1:length(time)
        delta_time(j) = time(j) - time(1);
    end
    
    lat = h5read(fname, '/latitude');
    lon = h5read(fname, '/longitude');
    lon = lon-360;
    
    time_hhmmss = h5read(fname, '/instrument_parameters/time_hhmmss');
    time_str = [fname(18:19), ':', fname(20:21), ':', fname(22:23)];
    date_str = fname(9:16);
    elev = h5read(fname, '/elevation');
    
    % Separating the data into iterative windows
    bin_window = floor(length(elev)/9);
    bin_count = round(length(elev)/bin_window);
    
    for i = 1:bin_count
        [lat_bin,lon_bin,time_bin_hhmmss,elev_bin] = atm_windowing_sub(lat,lon,time,elev,bin_window,i);
	
	
        % Setting up the histogram bins
        elev_sort = sort(elev_bin);
        [n, edges] = histcounts(elev_bin);
        
        if any(n > lake_thresh)
            
            elev_diff = 0;
            for j = 2:length(elev_bin)
                elev_diff(j) = -abs(elev_bin(j) - elev_bin(j-1)); % Difference transform anomalies for verification purposes
            end
            elev_diff_min = min(elev_diff);
            
            %-------------------------------------------------------------------------%
            % Lake Surface Identification - See atm_chunking_sub.m and
            % atm_sfc_detect_sub.m for debugging
            [window_lake_sfc,window_lake_btm] = atm_sfc_detect_sub(elev_bin);
            lake_start = find(~isnan(window_lake_sfc), 1, 'first');
            lake_end = find(~isnan(window_lake_sfc), 1, 'last');
            lon_bin_mean = movmean(lon_bin, 3000); % Remove scanning redundancy
            lat_bin_mean = movmean(lat_bin, 3000);
            
            % Great Circle Formula, for distance calculation
            lat1 = min([lat(1) lat(end)]);
            lon1 = min([lon(1) lon(end)]);
            for j = 1:length(lon_bin_mean)
                cent(j) = acosd(sind(lat1)*sind(lat_bin_mean(j)) + cosd(lat1)*cosd(lat_bin_mean(j))*cosd(abs(lon1-lon_bin_mean(j))));
            end
            cent_rad = cent*(3.14/180);
            dist = 6378*cent_rad'; % Units of km
            
            window_lake_btm(lake_start) = window_lake_sfc(lake_start); % Boundary conditions
            window_lake_btm(lake_end) = window_lake_sfc(lake_end);
            
            % Polynomial Fitting, with lake endpoints as boundary
            % conditions
            if lake_start == 1 % If lake is at beginning of data file
                idx = isnan(window_lake_btm);
                [p,~,mu] = polyfit(dist(~idx), window_lake_btm(~idx), 3);
                lake_btm_fitted = polyval(p, dist, [], mu);
            elseif any(length(elev_bin)-10:length(elev_bin) == lake_end) % If lake is at end of data file
                idx = isnan(window_lake_btm);
                [p,~,mu] = polyfit(dist(~idx), window_lake_btm(~idx), 3);
                lake_btm_fitted = polyval(p, dist, [], mu);
            else
                idx = isnan(window_lake_btm);
                p = polyfix(dist(~idx),window_lake_btm(~idx), 3, ...
                    [dist(lake_start) dist(lake_end)], ...
                    [window_lake_sfc(lake_start) window_lake_sfc(lake_end)]);
                lake_btm_fitted = polyval(p, dist);
            end
            
            
            lake_sfc_mean = movmean(window_lake_sfc, 1600, 'omitnan');
            lake_btm_mean = movmean(window_lake_btm, 1600, 'omitnan');
            % Remove Points Outside of Lake
            lake_btm_fitted(lake_btm_fitted>max(window_lake_sfc)) = NaN;
            for j = 1:length(window_lake_sfc)
                if j < lake_start
                    lake_btm_fitted(j) = NaN;
                    lake_sfc_mean(j) = NaN;
                    lake_btm_mean(j) = NaN;
                elseif j > lake_end
                    lake_btm_fitted(j) = NaN;
                    lake_sfc_mean(j) = NaN;
                    lake_btm_mean(j) = NaN;
                end
            end
            lake_btm_mean(lake_start) = lake_sfc_mean(lake_start);
            lake_btm_mean(lake_end) = lake_sfc_mean(lake_end);
            % Further filtering, for cleanliness
            if lake_btm_mean(lake_start+1)-lake_btm_mean(lake_start) < -1
                lake_btm_mean(lake_start) = NaN;
            elseif lake_btm_mean(lake_end)-lake_btm_mean(lake_end-1) > 1
                lake_btm_mean(lake_end) = NaN;
            end
            
            sfc_btm_diff_polyfit = lake_sfc_mean - lake_btm_fitted;
            lake_depth_polyfit = max(sfc_btm_diff_polyfit(lake_start:lake_end));
            sfc_btm_diff_means = lake_sfc_mean - lake_btm_mean;
            lake_depth_means = max(sfc_btm_diff_means(lake_start:lake_end));
            mean_depth = mean([lake_depth_means lake_depth_polyfit]);
            
            
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
            
            

            
            
            %-------------------------------------------------------------------------%
            
            if ~isempty(lake_depth_polyfit)
                if ~isnan(lake_depth_polyfit)
                    lake_count = lake_count + 1;
                    disp(['Lake found! The total lake count is now: ', num2str(lake_count)])
                    
                    % Lake volume (half spheroid approximation)
                    lake_radius = 0.5*abs(dist(lake_end)-dist(lake_start))*1000; % Units of m
                    if lake_radius > max([lake_depth_means lake_depth_polyfit])
                        lake_vs = 0.5*(4/3)*lake_radius^2*lake_depth_means;
                        lake_vp = 0.5*(4/3)*lake_radius^2*lake_depth_polyfit;
                    else
                        lake_vs = 0.5*(4/3)*lake_radius*lake_depth_means^2;
                        lake_vp = 0.5*(4/3)*lake_radius*lake_depth_polyfit^2;
                    end

                    lat_list = [lat_list lat_bin(marked_int)];
                    lon_list = [lon_list lon_bin(marked_int)];
                    mean_depth_list = [mean_depth_list lake_depth_means];
                    fitted_depth_list = [fitted_depth_list lake_depth_polyfit];
                    lake_vs_list = [lake_vs_list lake_vs];
                    lake_vp_list = [lake_vp_list lake_vp];
                    lake_radius_list = [lake_radius_list lake_radius];
                    
                    lake.Raw(i).Data(1).Name = 'binned elevation data';
                    lake.Raw(i).Data(2).Name = 'binned longitudes';
                    lake.Raw(i).Data(3).Name = 'binned latitudes';
                    
                    lake.Raw(i).Data(1).Dataset = elev_bin;
                    lake.Raw(i).Data(2).Dataset = lon_bin;
                    lake.Raw(i).Data(3).Dataset = lat_bin;
                    
                    if indyplot
                        % Individual Plots
                        figure;
                        plot(real(dist), elev_bin, '.', 'MarkerSize', 3, 'color', rgb('sky blue'))
                        hold on;
                        plot(real(dist), lake_btm_fitted, 'LineWidth', 2, 'color', rgb('rose'))
                        plot(real(dist), lake_sfc_mean,'LineWidth', 2, 'color', rgb('green'))
                        plot(real(dist), lake_btm_mean, 'b', 'LineWidth', 2)
                        plot(real(dist(marked_int)), depth_marker, 'r*', 'MarkerSize',12)
                        xlabel('Distance from start of data path [m]', 'FontSize',14, 'FontWeight','bold');
                        ylabel('Elevation [m]', 'FontSize',14, 'FontWeight','bold')
                        try
                            xlim([real(dist(1)) real(dist(end))])
                        catch
                            warning('X-axis limits need to be reversed.')
                            xlim([real(dist(end)) real(dist(1))])
                        end
                        ylim([min(elev_bin)-2 max(elev_bin)])
                        title(['ILATM1B - 20170719, ', fname(18:19), ':', fname(20:21), ':', fname(22:23)], 'FontSize', 14, 'FontWeight', 'bold')
                        grid on;
                        legend('Raw', 'Polyfit Bed', 'Signal Surface', 'Signal Bed', 'Max Depth')
                        pause; close all
                        close all
                    end
                    
                    if multiplot    
                        
                        % Removing unwanted profiles; can be removed or
                        % changed depending on the data used
                        if lake_count~=10 && lake_count~=13 && lake_count~=15
                            plot_count = lake_count;
                            if plot_count>10 && plot_count<13
                                plot_count = plot_count-1;
                            elseif plot_count==14
                                plot_count = 12;
                            end
                            
                            % Multi-Panel Plots
                            subplot(3, 4, plot_count)
                            plot(real(dist), elev_bin, '.', 'MarkerSize', 3, 'color', rgb('sky blue'))
                            hold on;
                            plot(real(dist), lake_btm_fitted, 'LineWidth', 2, 'color', rgb('rose'))
                            plot(real(dist), lake_sfc_mean,'LineWidth', 2, 'color', rgb('green'))
                            plot(real(dist), lake_btm_mean, 'b', 'LineWidth', 2)
                            plot(real(dist(marked_int)), depth_marker, 'r*', 'MarkerSize',12)
                            try
                                xlim([real(dist(1)) real(dist(end))])
                            catch
                                warning('X-axis limits need to be reversed.')
                                xlim([real(dist(end)) real(dist(1))])
                            end
                            ylim([min(elev_bin)-2 max(elev_bin)])
                            grid on;
                            title(['ILATM1B - 20170719, ', fname(18:23)], 'FontWeight', 'bold')
                        end

                    end
                end
            end
            
        end
        clear('cent')
    end
end

%% Lake Struct Development
lake.Detections.Data(1).Name = 'lake latitude';
lake.Detections.Data(2).Name = 'lake longitude';
lake.Detections.Data(3).Name = 'mean signal depth';
lake.Detections.Data(4).Name = 'polynomial-fitted depth';
lake.Detections.Data(5).Name = 'lake volume (signal)';
lake.Detections.Data(6).Name = 'lake volume (fitted)';
lake.Detections.Data(7).Name = 'lake radius';

lake.Detections.Data(1).Dataset = lat_list;
lake.Detections.Data(2).Dataset = lon_list;
lake.Detections.Data(3).Dataset = mean_depth_list;
lake.Detections.Data(4).Dataset = fitted_depth_list;
lake.Detections.Data(5).Dataset = lake_vs_list;
lake.Detections.Data(6).Dataset = lake_vp_list;
lake.Detections.Data(7).Dataset = lake_radius_list;

