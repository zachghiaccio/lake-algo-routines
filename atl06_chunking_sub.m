
function [snr_chunk,elev_chunk,bound_range] = atl06_chunking_sub(elev_bin,snr_bin,chunk);
% The ATL06 equivalent to is2_chunking_sub.m - breaks down windowed data into smaller chunks
% to obtain histogram counts for more reliable surface detection.
% Histograms aren't needed (yet?), so are replaced with SNR

chunk_window = floor(length(elev_bin)/15);
chunk_count = round(length(elev_bin)/chunk_window);
if chunk > chunk_count
	error('Binned data does not have that many chunks available.')
elseif chunk < 1
	error('Invalid input for chunk.')
end

bound = chunk_window*chunk;
bound_range = (bound-chunk_window+1):bound;
elev_chunk = elev_bin(bound_range);
snr_chunk = snr_bin(bound_range);

end
