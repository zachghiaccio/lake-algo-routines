function [window_lake_sfc,window_lake_btm] = is2_sfc_detect_sub(elev_bin, class_bin)
% A subroutine for autoatm4.m, designed to chunk the data and better
% identify lake surfaces/bottoms in ATM data

addpath(genpath('/Users/zhfair/Documents/'))

plotting = 0; % Yes=1, No=0

% Elevation window information
window_mean = movmean(elev_bin,3000, 'omitnan');
window_std = movstd(elev_bin, 3000, 'omitnan');
window_lake_sfc = NaN(length(elev_bin),1);
window_lake_btm = NaN(length(elev_bin),1);

elev_bin = round(elev_bin, 3);
% Chunking
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
    disp(k)
    if elev_chunk_std < 0.05 && ~isnan(lake_points_std) % Chunks without lakes tend to have std > 0.05
        
        % Acceptable data points for lake surface
        std_bound_upper = lake_points_mean + 3*lake_points_std;
        std_bound_lower = lake_points_mean - 3*lake_points_std;
        elev_chunk_sfc = elev_chunk;
        elev_chunk_sfc(elev_chunk_sfc > std_bound_upper) = NaN;
        elev_chunk_sfc(elev_chunk_sfc < std_bound_lower) = NaN;
        
        window_lake_sfc(bound_range) = elev_chunk_sfc;
        first_sfc_point = find(~isnan(window_lake_sfc),1,'first');
        last_sfc_point = find(~isnan(window_lake_sfc),1,'last');
        window_lake_sfc( (movstd(elev_chunk,1600,'omitnan')<0.71)+bound_range(1) ) = NaN;
        
        % Acceptable data points for the lake bottom
        window_lake_btm(min(bound_range):max(bound_range)) = elev_chunk;
        window_lake_btm(window_lake_btm>window_mean-1.8*window_std) = NaN;
        window_lake_btm(window_lake_btm > nanmean(window_lake_sfc(bound_range))) = NaN;
        window_lake_btm( (movstd(elev_chunk,1600,'omitnan')<0.71)+bound_range(1) ) = NaN;
        window_lake_sfc( (movstd(elev_chunk,1600,'omitnan')<0.71)+bound_range(1) ) = NaN;
        if isnan(window_lake_btm) % To help account for false positives
            window_lake_sfc(bound_range) = NaN;
        end
        
        % Signal/Noise Count Tests
        class_chunk = class_bin(bound_range);
        high_count = sum(class_chunk == 4);
        lower_count = sum(class_chunk==3 | class_chunk==2);
        noise_count = sum(class_chunk == 0);
        disp(['High SNR: ', num2str(high_count./(high_count + noise_count))])
        disp(['Lower SNR: ', num2str(lower_count./(lower_count + noise_count))])
        
        %keyboard;
        
        % Non-surface data points, for debugging
        non_surface = elev_chunk;
        non_surface(non_surface<std_bound_upper & non_surface>std_bound_lower) = NaN;
        
        if plotting && ~isnan(window_lake_sfc(bound_range)) %For debigging and testing purposes
            figure;
            histogram(elev_chunk)
            
            figure;
            hold on;
            plot(1:length(elev_chunk), non_surface, '.')
            plot(1:length(elev_chunk),elev_chunk_sfc, '.')
            pause; close all
        end
    end
end

end
