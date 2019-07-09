function [flag,window_lake_btm,lake_btm_mean,lake_btm_error] = is2_ph_dist(d_bin, elev_bin)
% A subroutine of autoatl.m to examine the photon distribution within a 
% 40-meter along-track segment.
% d_bin is the binned along-track distance, in km
% elev_bin is the binned elevation data, in m
% flag is a boolean telling the main script whether there are sufficient 
% photons for refinement.

d_bin = d_bin.*1000; % Convert to m

% Initial conditions for while loop
segment_length = 0;
lower_bound = 1;
upper_bound = lower_bound + 1;
h_w_prev = 1; % Default window height, in meters

% Array presets
window_lake_btm = NaN(length(elev_bin),1);
lake_btm_mean = NaN(length(elev_bin),1);
lake_btm_error = NaN(length(elev_bin),1);

%% Segmentation Sequence
while segment_length<40 && upper_bound<length(d_bin)
    segment_length = d_bin(upper_bound) - d_bin(lower_bound);
    
    % Only runs statistics if there is a 40 m segment of photons (N~190)
    if segment_length >= 40
        elev_segment = elev_bin(lower_bound:upper_bound);
        d_segment = d_bin(lower_bound:upper_bound);
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
        
        
        if ph_count>=10 && d_ph_max>=20
            %% Linear Least-Squares Regression
            % Outputs: Residual (r), Mean elevation (elev_lsqf_mean), Slope (elev_lsqf_slope)
            idx = isnan(elev_segment);
            [p,~,mu] = polyfit(d_ph(~idx),elev_segment(~idx),1);
            elev_lsqf = polyval(p,d_ph,[],mu);
            r = elev_segment - elev_lsqf;
            elev_lsqf_slope = ( elev_lsqf(last_ph)-elev_lsqf(first_ph) )/d_ph_max;
            
            %% Residual Analysis
            % Outputs: Median residual (r_med), Observed spread
            % (robust_spread), Predicted spread (expected_spread)
            r_med = median(r, 'omitnan');
            [robust_spread,~,~] = is2_robust_spread(elev_segment);
            expected_spread = (3e8/2)*sqrt( (0.68e-9)^2 + ((17/(24e8))*...
                tan(elev_lsqf_slope))^2 ); % From Smith et al., (2019)
            
            %% Lake Bed Finding
            % Outputs: Lake bottom (window_lake_btm), Mean lake bottom
            % (lake_btm_mean), Lake bottom error (lake_btm_error)
            h_w = max([6*robust_spread,6*expected_spread,0.75*h_w_prev,1]);
            r_diff = abs(r-r_med);
            elev_segment(r_diff<0.75*h_w) = NaN;
            
            idx = isnan(elev_segment);
            [p,S,mu] = polyfit(d_ph(~idx),elev_segment(~idx),1);
            [elev_lsqf,delta] = polyval(p,d_ph,S,mu);
            elev_lsqf_mean = nanmean(elev_lsqf);
            
            % Function outputs
            window_lake_btm(lower_bound:upper_bound) = elev_segment;
            lake_btm_mean(lower_bound:upper_bound) = elev_lsqf_mean;
            lake_btm_error(lower_bound:upper_bound) = delta;

            % Reset segment bounds
            flag = 1;
            segment_length = 0;
            lower_bound = upper_bound;
            upper_bound = lower_bound + 1;
            
        else
            flag = 0;
            disp('Not enough photons found! Moving on to next segment.')
            
            segment_length = 0;
            lower_bound = upper_bound;
            upper_bound = lower_bound + 1;
        end
        
    else
        upper_bound = upper_bound + 1;
    end
end


end
