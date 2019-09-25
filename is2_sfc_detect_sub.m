function [window_lake_sfc] = is2_sfc_detect_sub(elev_bin)
% A subroutine for autoatl.m, designed to chunk the data and better
% identify lake surfaces in IS2 data. Bottom finding is now relegated to
% is2_ph_dist.m

addpath(genpath('/Users/zhfair/Documents/'))

% Elevation window information
window_lake_sfc = NaN(length(elev_bin),1);
elev_bin = round(elev_bin, 3);

%% Chunking
for k = 1:15
    [n_hist, ~, elev_chunk,bound_range] = is2_chunking_sub(elev_bin, k);
    perc = prctile(n_hist,1:100);
    [~,index] = min(abs(perc' - 1200));

    % Sorting chunked data by frequency of occurrence
    [~,~,subs] = unique(elev_chunk);
    freq = accumarray(subs,subs,[],@numel);
    [~,sortie] = sort(freq(subs));
    elev_chunk_sort = elev_chunk(sortie);
    
    % Finding the lake surface (if any) and statistics
    elev_sort_std = movstd(elev_chunk_sort, 600); 
    lake_points = elev_chunk_sort(elev_sort_std<0.07); % To find lake flatness
    lake_points_std = std(lake_points); 
    lake_points_mean = mean(lake_points);
    
    % Standard deviation check for the sorted data
    perc_bound = round(0.01*index*length(elev_chunk_sort));
    elev_chunk_std = std(elev_chunk_sort(perc_bound:end), 'omitnan');

    if elev_chunk_std < 0.05 && ~isnan(lake_points_std) % Chunks without lakes tend to have std > 0.05
        
        % Acceptable data points for lake surface
        std_bound_upper = lake_points_mean + 3*lake_points_std;
        std_bound_lower = lake_points_mean - 3*lake_points_std;
        elev_chunk_sfc = elev_chunk;
        elev_chunk_sfc(elev_chunk_sfc > std_bound_upper) = NaN;
        elev_chunk_sfc(elev_chunk_sfc < std_bound_lower) = NaN;
        
        window_lake_sfc(bound_range) = elev_chunk_sfc;
        window_lake_sfc( (movstd(elev_chunk,1600,'omitnan')<0.71)+bound_range(1) ) = NaN;
        
    end
end

end
