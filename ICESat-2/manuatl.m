function [ds_mean,ds_max,dp_mean,dp_max,sigma_btm,delta_p] = manuatl(fl_atl03,lat_in,lon_in,extent)
% If a supraglacial lake location is known, approximate the lake depth
% using available photon data.
% May use a lat (recommended) or lon as input, with 5 maximum significant
% digits (e.g. 79.123)

% Lat/lon input must be a single value
if length(lat_in)>1 || length(lon_in)>1
    error('Lat/lon input must be a single number.')
end

% Check for extent input
if ~exist('extent', 'var')
    warning('Input for swath extent missing. Defaulting to 10000.')
    extent = 10000;
end

% Photon height information
lat = h5read(fl_atl03,'/gt2l/heights/lat_ph');
lon = h5read(fl_atl03,'/gt2l/heights/lon_ph');
class = h5read(fl_atl03,'/gt2l/heights/signal_conf_ph');
elev = h5read(fl_atl03, '/gt2l/heights/h_ph');

% Photon geolocation information
theta_elev = h5read(fl_atl03, '/gt2l/geolocation/ref_elev');
theta_elev(theta_elev>2*pi) = NaN;
theta_elev_mean = nanmean(theta_elev); % Elevation angle change is small
theta_inc = pi/2 - theta_elev_mean; % Incidence angle, for refraction correction
dist = h5read(fl_atl03, '/gt2l/heights/dist_ph_along'); % In meters
seglen = h5read(fl_atl03, '/gt2l/geolocation/segment_length'); % Only used to properly convert dist
total_dist = sum(seglen) + dist(end);
dist_corr = linspace(0,total_dist, length(dist))./1000; % Also convert to km

% Remove surface type from photon classification
class_consol = is2_class_merge(class);

% Remove solar background photons
signal = elev;
signal(class_consol<1) = NaN; % Signal and buffer photons
hml = signal; % High-confidence photons only
hml(class_consol<4) = NaN;

% Define a range of photons centered at the input lat/lon
if exist('lat_in','var')
    idx = find(round(lat,3)==lat_in,1,'first');
    latx = lat(idx-extent:idx+extent);
    lonx = lon(idx-extent:idx+extent);
    sigx = signal(idx-extent:idx+extent);
    hx = hml(idx-extent:idx+extent);
    dx = dist_corr(idx-extent:idx+extent);
    
    lon_in = -999; % Filler value to prevent input error
elseif exist('lon_in','var')
    idx = find(round(lon,3)==lon_in,1,'first');
    latx = lat(idx-extent:idx+extent);
    lonx = lon(idx-extent:idx+extent);
    sigx = signal(idx-extent:idx+extent);
    hx = hml(idx_extent:idx+extent);
    dx = dist_corr(idx-extent:idx+extent);
    
    lat_in = -999; % Filler value
else
    error('No lat/lon input(s).')
end

%% Surface, Bed Determination
sfc = is2_sfc_detect_sub(sigx);

% Define the start and end of lake
lake_start = find(~isnan(sfc),1,'first');
lake_end = find(~isnan(sfc),1,'last');
before_lake = 1:lake_start-1;
after_lake = lake_start+1:length(sfc);

% Lake bottom approximation
window_mean = movmean(sigx,1600,'omitnan');
window_std = movstd(sigx,1600,'omitnan');
sig_btm = sigx;
sig_btm(sig_btm>=window_mean-0.75*window_std) = NaN;
sig_btm(sig_btm<=window_mean-1.8*window_std) = NaN;
sigma_btm = movstd(sig_btm,100,'omitnan'); % Lake bottom uncertainty

if sum(~isnan(sig_btm)) > 0
    lake_flag = 1;
else
    lake_flag = 0;
end

% [~,~,subs] = unique(high);
% freq = accumarray(subs,subs,[],@numel);
% [~,sortie] = sort(freq(subs));
% sigsort = high(sortie);
% 
% sig_sort_std = movstd(sigsort,600,'omitnan');
% lake_points = sigsort(sig_sort_std<0.15);
% lake_points_std = std(lake_points);
% lake_points_mean = mean(lake_points);
% 
% std_bound_upper = lake_points_mean + 3*lake_points_std;
% std_bound_lower = lake_points_mean - 3*lake_points_std;
% lake_sfc = high;
% lake_sfc(lake_sfc > std_bound_upper) = NaN;
% lake_sfc(lake_sfc < std_bound_lower) = NaN;


%% Depth Estimation
if lake_flag
    idx = (sigx == sig_btm); % Needed for visualization purposes
    
    lake_sfc_mean = movmean(sfc,600,'omitnan');
    lake_btm_mean = movmean(sig_btm,600,'omitnan');
    depth = lake_sfc_mean-lake_btm_mean;
    
    % Refraction correction
    [~,final_depth,dz] = depth_refrac_fix(dx,depth,theta_inc);
    
    % Correct lake bed arrays
    sig_btm = sig_btm + dz;
    lake_btm_mean = movmean(sig_btm,600,'omitnan');
    sigx(idx) = sigx(idx) + dz(idx);
    
    % Signal depth statistics
    ds_mean = nanmean(final_depth);
    ds_max = nanmax(final_depth);
    
    % Polynomial fitting
    idx = isnan(sig_btm);
    [p,S,mu] = polyfit(dx(~idx)',sig_btm(~idx),3);
    [poly_btm,~] = polyval(p,dx,S,mu);
    poly_btm(poly_btm>max(lake_sfc_mean)) = NaN; 
    
    % Polynomial error
    resid = sig_btm - poly_btm';
    delta_p = sqrt(nansum(resid(~idx).^2)/length(sig_btm(~idx)));
    
    % Polynomial-fitted depths and statistics
    dp = lake_sfc_mean - poly_btm';
    dp_mean = nanmean(dp);
    dp_max = nanmax(dp);
    
    sigx(sigx<=window_mean-1.8*window_std) = NaN; %For visualization
end


figure;
plot(dx,sigx, '.', 'MarkerSize', 6, 'color',rgb('sky blue'))
hold on; plot(dx,lake_sfc_mean, 'LineWidth', 3, 'color', rgb('forest green'))
plot(dx,lake_btm_mean, 'b', 'LineWidth', 2)
plot(dx(lake_start:lake_end),poly_btm(lake_start:lake_end),'LineWidth',2, 'color',rgb('orange'))
legend('Raw', 'Surface', 'Bed', 'Polyfit')

itlooksgood = input('Does the plot provide a good lake profile? 0 if no, 1 if yes: ');
if itlooksgood
    disp(' ')
    disp(['Mean Depth: ', num2str(ds_mean)])
    disp(['Max Depth: ', num2str(ds_max)])
    disp(['Mean Depth (Polyfit): ', num2str(dp_mean)])
    disp(['Max Depth (Polyfit): ', num2str(dp_max)])
    disp(['Mean Depth Uncertainty: ', num2str(nanmean(sigma_btm))])
    disp(['Polyfit Depth Error: ', num2str(delta_p)])
else
    %% Photon Refinement
    close all;
    refine_flag = input('What changes are needed? 1 = Surface only, 2 = Bed only, 3 = Both ');
    
    % Surface Refinement only %
    %-------------------------%
    if refine_flag == 1
        [sfc_ref,~,~,~] = is2_sfc_find(dx,sigx,sfc); % Surface refinement
        lake_sfc_mean = movmean(sfc_ref,600,'omitnan');
        depth = lake_sfc_mean-lake_btm_mean;
        
        % Refraction correction
        [~,final_depth,dz] = depth_refrac_fix(dx,depth,theta_inc);
        
        % Correct lake bed arrays
        idx = (sigx == sig_btm); % Needed for visualization purposes
        sig_btm = sig_btm + dz;
        lake_btm_mean = movmean(sig_btm,600,'omitnan');
        sigx(idx) = sigx(idx) + dz(idx);
        
        % Signal depth statistics
        ds_mean = nanmean(final_depth);
        ds_max = nanmax(final_depth);
        
        % Polynomial fitting
        idx = isnan(sig_btm);
        [p,S,mu] = polyfit(dx(~idx)',sig_btm(~idx),3);
        [poly_btm,~] = polyval(p,dx,S,mu);
        poly_btm(poly_btm>max(lake_sfc_mean)) = NaN;
        
        % Polynomial error
        resid = sig_btm - poly_btm';
        delta_p = sqrt(nansum(resid(~idx).^2)/length(sig_btm(~idx)));
        
        % Polynomial-fitted depths and statistics
        dp = lake_sfc_mean - poly_btm';
        dp_mean = nanmean(dp);
        dp_max = nanmax(dp);
        
        lake_btm_mean(lake_btm_mean>max(lake_sfc_mean)) = NaN;
        sigx(sigx<=window_mean-1.8*window_std) = NaN; %For visualization
        
        % Bed Refinement only %
        %---------------------%
    elseif refine_flag == 2
        [~,sig_btm,~,~] = is2_ph_dist(dx,sigx);
        lake_btm_mean = movmean(sig_btm,300,'omitnan');
        lake_btm_mean(lake_btm_mean>max(lake_sfc_mean)) = NaN;
        sigma_btm = movstd(sig_btm,100,'omitnan');
        depth = lake_sfc_mean-lake_btm_mean;
        
        idx = (sigx == sig_btm); % Needed for visualization purposes
        % Refraction correction
        [~,final_depth,dz] = depth_refrac_fix(dx,depth,theta_inc);
        
        % Correct lake bed arrays
        sig_btm = sig_btm + dz;
        lake_btm_mean = movmean(sig_btm,600,'omitnan');
        lake_btm_mean(lake_btm_mean>nanmean(lake_sfc_mean)) = NaN;
        sigx(idx) = sigx(idx) + dz(idx);
        
        % Signal depth statistics
        ds_mean = nanmean(final_depth);
        ds_max = nanmax(final_depth);
        
        % Polynomial fitting
        idx = isnan(sig_btm);
        [p,S,mu] = polyfit(dx(~idx)',sig_btm(~idx),3);
        [poly_btm,~] = polyval(p,dx,S,mu);
        poly_btm(poly_btm>max(lake_sfc_mean)) = NaN;
        
        % Polynomial error
        resid = sig_btm - poly_btm';
        delta_p = sqrt(nansum(resid(~idx).^2)/length(sig_btm(~idx)));
        
        % Polynomial-fitted depths and statistics
        dp = lake_sfc_mean - poly_btm';
        dp_mean = nanmean(dp);
        dp_max = nanmax(dp);
        
        sigx(sigx<=window_mean-2*window_std) = NaN; %For visualization
        lake_btm_mean(lake_btm_mean==max(lake_btm_mean)) = NaN;
        
        % Surface and Bed Refinement %
        %----------------------------%
    elseif refine_flag>=3 || refine_flag<1
        if refine_flag > 3
            warning('Invalid refinement flag input. Refining both surface and bed by default.')
        end
        
        % Refinement algorithms
        [sfc_ref,~,~,~] = is2_sfc_find(dx,sigx,sfc);
        [~,sig_btm,~,sigma_btm] = is2_ph_dist(dx,sigx);
        lake_sfc_mean = movmean(sfc_ref,600,'omitnan');
        lake_btm_mean = movmean(sig_btm,600,'omitnan');
        lake_btm_mean(lake_btm_mean>max(lake_sfc_mean)) = NaN;
        depth = lake_sfc_mean-lake_btm_mean;
        
        % Refraction correction
        [~,final_depth,dz] = depth_refrac_fix(dx,depth,theta_inc);
        
        % Correct lake bed arrays
        sig_btm = sig_btm + dz;
        lake_btm_mean = movmean(sig_btm,600,'omitnan');
        lake_btm_mean(lake_btm_mean>mean(lake_sfc_mean)) = NaN;
        sigx(idx) = sigx(idx) + dz(idx);
        
        % Signal depth statistics
        ds_mean = nanmean(final_depth);
        ds_max = nanmax(final_depth);
        
        % Polynomial fitting
        idx = isnan(sig_btm);
        [p,S,mu] = polyfit(dx(~idx)',sig_btm(~idx),3);
        [poly_btm,~] = polyval(p,dx,S,mu);
        poly_btm(poly_btm>max(lake_sfc_mean)) = NaN;
        
        % Polynomial error
        resid = sig_btm - poly_btm';
        delta_p = sqrt(nansum(resid(~idx).^2)/length(sig_btm(~idx)));
        
        % Polynomial-fitted depths and statistics
        dp = lake_sfc_mean - poly_btm';
        dp_mean = nanmean(dp);
        dp_max = nanmax(dp);
        
        sigx(sigx<=window_mean-1.8*window_std) = NaN; %For visualization
    end
    
    
    figure;
    plot(dx,sigx, '.', 'MarkerSize', 6, 'color',rgb('sky blue'))
    hold on; plot(dx,lake_sfc_mean, 'LineWidth', 3, 'color',rgb('forest green'))
    plot(dx,lake_btm_mean, 'b', 'LineWidth', 2)
    plot(dx(lake_start:lake_end),poly_btm(lake_start:lake_end),'LineWidth',2, 'color',rgb('orange'))
    legend('Raw', 'Surface', 'Bed', 'Polyfit')
    
    disp(' ')
    disp(['Mean Depth: ', num2str(ds_mean)])
    disp(['Max Depth: ', num2str(ds_max)])
    disp(['Mean Depth (Polyfit): ', num2str(dp_mean)])
    disp(['Max Depth (Polyfit): ', num2str(dp_max)])
    disp(['Mean Depth Uncertainty: ', num2str(nanmean(sigma_btm))])
    disp(['Polyfit Depth Error: ', num2str(delta_p)])
end

pause;
%% H5 File Exportation
lake_label = input('Please input a label for the lake of interest. ');

% Generate data structure
% h5create('is2_depth_data.h5', ['/',lake_label,'/raw/height'], size(sigx));
% h5create('is2_depth_data.h5', ['/',lake_label,'/raw/latitude'], size(sigx));
% h5create('is2_depth_data.h5', ['/',lake_label,'/raw/longitude'], size(sigx));
% h5create('is2_depth_data.h5', ['/',lake_label,'/raw/along_track_distance'], size(dx));

% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/raw_surface'], size(sfc));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/raw_bed'], size(sig_btm));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/mean_surface'], size(lake_sfc_mean));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/mean_bed'], size(lake_btm_mean));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/polyfit_bed'], size(poly_btm));
% 
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/depths/signal_depth'], size(final_depth));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/depths/polyfit_depth'], size(dp));
% 
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/mean_signal_depth'], size(ds_mean));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/mean_polyfit_depth'], size(dp_mean));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/max_signal_depth'], size(ds_max));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/max_polyfit_depth'], size(dp_max));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/sigma_bed'], size(sigma_btm));
% h5create('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/polyfit_rmse'], size(delta_p));


% Write data to file
% h5write('is2_depth_data.h5', ['/',lake_label,'/raw/height'], sigx);
% h5write('is2_depth_data.h5', ['/',lake_label,'/raw/latitude'], latx);
% h5write('is2_depth_data.h5', ['/',lake_label,'/raw/longitude'], lonx);
% h5write('is2_depth_data.h5', ['/',lake_label,'/raw/along_track_distance'], dx);

% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/raw_surface'], sfc);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/raw_bed'], sig_btm);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/mean_surface'], lake_sfc_mean);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/mean_bed'], lake_btm_mean);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/lake_parameters/polyfit_bed'], poly_btm);
% 
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/depths/signal_depth'], final_depth);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/depths/polyfit_depth'], dp);
% 
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/mean_signal_depth'], ds_mean);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/mean_polyfit_depth'], dp_mean);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/max_signal_depth'], ds_max);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/max_polyfit_depth'], dp_max);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/sigma_bed'], sigma_btm);
% h5write('is2_depth_data.h5', ['/',lake_label,'/lsbs/statistics/polyfit_rmse'], delta_p);

end
