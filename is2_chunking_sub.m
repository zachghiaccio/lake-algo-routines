function[n_hist,hist_edges,elev_chunk,bound_range] = is2_chunking_sub(elev_bin, chunk)
% Function to divide binned data into smaller chunks and obtaining
% histogram counts
% elev_bin is the elevation profile of the binned data
% chunk refers to which chunk of the binned data you want to analyze. You
% must input a number between 1-chunk_count (chunk_count = 12 right now)
%-------------------------------------------------------------------------%

chunk_window = floor(length(elev_bin)/15);
chunk_count = round(length(elev_bin)/chunk_window);
if chunk > chunk_count
    error('Binned data does not have that many chunks.')
elseif chunk < 1
    error('Invalid input for chunk.')
end

bound = chunk_window*chunk;
bound_range = (bound-chunk_window+1):bound;
elev_chunk = elev_bin(bound_range);

[n_hist,hist_edges] = histcounts(elev_chunk);

end