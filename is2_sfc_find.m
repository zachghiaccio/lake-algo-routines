function [window_lake_sfc_corr,lake_sfc_mean,lake_sfc_error] = is2_sfc_find(d_bin,elev_bin,window_lake_sfc)
% A script that attempts to fill gaps in the lake surface that are caused
% by filters. This is mostly for cleanliness, but it also allows for more
% of the lake bottom to be identified.

d_bin = d_bin.*1000; % Convert to m

% Initial Conditions
segment_length = 0;
lower_bound = 1;
upper_bound = lower_bound + 1;
h_w_prev = 1;

% Presets
sfc_nans = isnan(window_lake_sfc);
chunk_size = floor(length(window_lake_sfc)/15);
number_of_chunks = round(length(window_lake_sfc)/chunk_size);
window_lake_sfc_corr = window_lake_sfc;
lake_sfc_mean = NaN(length(window_lake_sfc),1);
lake_sfc_error = NaN(length(window_lake_sfc),1);

for chunk = 1:number_of_chunks
    bound = chunk_size*chunk;
    bound_range = (bound-chunk_size+1):bound;
    
    lake_sfc_bounded = window_lake_sfc(bound_range);
    if sum(sfc_nans(bound_range)) == length(bound_range)
        there_is_data = 0;
    else
        there_is_data = 1;
        elev_bin_bounded = elev_bin(bound_range);
        d_bin_bounded = d_bin(bound_range);
    end
    
    if there_is_data
        lake_start = find(~isnan(lake_sfc_bounded),1,'first');
        lake_end = find(~isnan(lake_sfc_bounded),1,'last');
        
        while segment_length<40 && upper_bound<length(bound_range)
            segment_length = d_bin_bounded(upper_bound) - d_bin_bounded(lower_bound);
            
            % Only runs statistics if there is a 40 m segment of photons (N~190)
            if segment_length >= 40
                elev_segment = elev_bin_bounded(lower_bound:upper_bound);
                d_segment = d_bin_bounded(lower_bound:upper_bound);
                first_ph = find(~isnan(elev_segment),1,'first');
                last_ph = find(~isnan(elev_segment),1,'last');
                mid_ph = round( (first_ph+last_ph)/2 );
                if isempty(mid_ph)
                    ph_count = -999;
                    d_ph_max = -999;
                else
                    ph_count = length(elev_segment(~isnan(elev_segment)));
                    d_ph = d_segment - d_segment(mid_ph);
                    d_ph = d_ph.';
                    d_ph_max = d_segment(last_ph) - d_segment(first_ph);
                end
                
                if ph_count >=10 && d_ph_max>=20
                    %% Linear Least-Squares Regression
                    % Outputs: Residuals (r), Mean elevation (elev_lsqf_mean), Slope (elev_lsqf_slope)
                    idx = isnan(elev_segment);
                    [p,~,mu] = polyfit(d_ph(~idx),elev_segment(~idx),1);
                    elev_lsqf = polyval(p,d_ph,[],mu);
                    r = elev_segment - elev_lsqf;
                    elev_lsqf_slope = ( elev_lsqf(last_ph)-elev_lsqf(first_ph) )/d_ph_max;
                    
                    %% Residual Analysis
                    % Outputs: Median residual (r_med), Observed spread (robust_spread), Predicted spread (expected_spread)
                    r_med = median(r, 'omitnan');
                    [robust_spread,~,~] = is2_robust_spread(elev_segment);
                    expected_spread = (3e8/2)*sqrt( (0.68e-9)^2 + ((17/(24e8))*tan(elev_lsqf_slope))^2 ); % From Smith et al., (2019)
                    
                    %% Lake Surface Finding
                    % Outputs: Lake surface (window_lake_sfc_corr), Mean lake surface (lake_sfc_mean), Lake surface error (lake_sfc_error)
                    h_w = max([6*robust_spread,6*expected_spread,0.75*h_w_prev,1]);
                    r_diff = abs(r-r_med);
                    elev_segment(r_diff>0.75*h_w) = NaN; % Note the inequality difference from bottom-finding
                    
                    idx = isnan(elev_segment);
                    [p,S,mu] = polyfit(d_ph(~idx),elev_segment(~idx),1);
                    [elev_lsqf,delta] = polyval(p,d_ph,S,mu);
                    elev_lsqf_mean = nanmean(elev_lsqf);
                    
                    % Function outputs
                    window_lake_sfc_corr(lower_bound:upper_bound) = elev_segment;
                    lake_sfc_mean(lower_bound:upper_bound) = elev_lsqf_mean;
                    lake_sfc_error(lower_bound:upper_bound) = delta;
                    
                    % Reset segment bounds
                    segment_length = 0;
                    lower_bound = upper_bound;
                    upper_bound = lower_bound + 1;
                    
                else
                    disp('Not enough photons found! Moving on to next segment.')
                    lower_bound = upper_bound;
                    upper_bound = lower_bound + 1;
                end
            else
                upper_bound = upper_bound + 1;
            end
        end
        
        % Remove stray data points before/after the lake start/end
        ints = 1:length(lake_sfc_bounded);
        not_in_lake = (ints<lake_start | ints>lake_end);
        window_lake_sfc_corr(not_in_lake) = NaN;
        lake_sfc_mean(not_in_lake) = NaN;
        lake_sfc_error(not_in_lake) = NaN;
        
    end
end

end
