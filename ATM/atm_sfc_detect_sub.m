function [window_lake_sfc,window_lake_btm] = atm_sfc_detect_sub(elev_bin)
% A subroutine for autoatm4.m, designed to chunk the data and better
% identify lake surfaces/bottoms in ATM data

addpath(genpath('/Users/zhfair/Documents/'))

plotting = 0; % Yes=1, No=0

% Elevation window information
window_mean = movmean(elev_bin,3000);
window_std = movstd(elev_bin, 3000);
window_lake_sfc = NaN(length(elev_bin),1);
window_lake_btm = NaN(length(elev_bin),1);

% Chunking
for k = 1:9
    [n_hist, ~, elev_chunk,bound_range] = atm_chunking_sub(elev_bin, k);
    perc = prctile(n_hist,1:100);
    [~,index] = min(abs(perc' - 1200));
    
    % Sorting chunked data by frequency of occurrence
    [~,~,subs] = unique(elev_chunk);
    freq = accumarray(subs,subs,[],@numel);
    [~,sortie] = sort(freq(subs));
    elev_chunk_sort = elev_chunk(sortie);
    
    % Finding the lake surface (if any) and statistics
    elev_sort_std = movstd(elev_chunk_sort, 1600); 
    lake_points = elev_chunk_sort(elev_sort_std<0.02); % To find lake flatness
    lake_points_std = std(lake_points); 
    lake_points_mean = mean(lake_points);
    
    % Standard deviation check for the sorted data
    perc_bound = round(0.01*index*length(elev_chunk_sort));
    elev_chunk_std = std(elev_chunk_sort(perc_bound:end), 'omitnan');
    
    if elev_chunk_std < 0.009 && ~isnan(lake_points_std) % Chunks without lakes tend to have std > 0.009
        
        % Acceptable data points for lake surface
        std_bound_upper = lake_points_mean + 3*lake_points_std;
        std_bound_lower = lake_points_mean - 3*lake_points_std;
        elev_chunk_sfc = elev_chunk;
        elev_chunk_sfc(elev_chunk_sfc > std_bound_upper) = NaN;
        elev_chunk_sfc(elev_chunk_sfc < std_bound_lower) = NaN;
        
        window_lake_sfc(bound_range) = elev_chunk_sfc;
        first_sfc_point = find(~isnan(window_lake_sfc),1,'first');
        last_sfc_point = find(~isnan(window_lake_sfc),1,'last');
        
        % Acceptable data points for the lake bottom
        window_lake_btm(min(bound_range):max(bound_range)) = elev_chunk;
        window_lake_btm(window_lake_btm>window_mean-1.8*window_std) = NaN;
        window_lake_btm(movstd(window_lake_btm,1600,'omitnan')<0.04) = NaN;
        
        
        
        % Non-surface data points, for debugging
        non_surface = elev_chunk;
        non_surface(non_surface<std_bound_upper & non_surface>std_bound_lower) = NaN;
        
        if plotting %For debigging and testing purposes
            figure;
            histogram(elev_chunk)
            
            figure;
            plot(1:length(elev_chunk),elev_chunk_sfc, '.')
            hold on;
            plot(1:length(elev_chunk), non_surface, '.')
            pause; close all
        end
    end
end

end
