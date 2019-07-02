function [lake_btm_fitted] = is2_polyfit_sub(lon_bin_csb, window_lake_btm)
% Subroutine to fit a polynomial curve to each lake
% This accounts for multiple lakes within a single data window
% Currently using 3rd-order polynomial curves
tic;
chunk_size = floor(length(window_lake_btm)/15);
number_of_chunks = round(length(window_lake_btm)/chunk_size);

lake_btm_fitted = NaN(length(window_lake_btm),1);
for chunk = 1:number_of_chunks
    bound = chunk_size*chunk;
    bound_range = (bound-chunk_size+1):bound;
    
    lon_bounded = lon_bin_csb(bound_range);
    lake_btm_bounded = window_lake_btm(bound_range);
    
    if sum(isnan(lake_btm_bounded)) ~= length(lake_btm_bounded)
        lake_start = find(~isnan(lake_btm_bounded), 1, 'first');
        lake_end = find(~isnan(lake_btm_bounded), 1, 'last');
        
        if lake_start == 1 % If lake is at beginning of data file
            idx = isnan(lake_btm_bounded);
            [p,~,mu] = polyfit(lon_bounded(~idx), lake_btm_bounded(~idx), 3);
            lake_btm_fitted(bound_range) = polyval(p, lon_bounded, [], mu);
        elseif any(length(lake_btm_bounded)-10:length(lake_btm_bounded) == lake_end) % If lake is at end of data file
            idx = isnan(lake_btm_bounded);
            [p,~,mu] = polyfit(lon_bounded(~idx), lake_btm_bounded(~idx), 3);
            lake_btm_fitted(bound_range) = polyval(p, lon_bounded, [], mu);
        else
            idx = isnan(lake_btm_bounded);
            p = polyfix(lon_bounded(~idx),lake_btm_bounded(~idx), 3, ...
                [lon_bounded(lake_start) lon_bounded(lake_end)], ...
                [lake_btm_bounded(lake_start) lake_btm_bounded(lake_end)]);
            lake_btm_fitted(bound_range) = polyval(p, lon_bounded);
        end
    else
        continue;
    end
    toc;
end