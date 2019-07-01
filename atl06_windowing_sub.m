function [lat_bin, lon_bin, time_bin, elev_bin, snr_bin, dist_bin] = atl06_windowing_sub(lat,lon,time,elev,snr,dist,window_size,window_count)
% The ATL06 equivalent to is2_windowing_sub.m - separates ATL06 arrays into discrete windows for analysis.
% Instead of the "class" variable, uses "snr" in smaller arrays.

bound = window_size*window_count;

if bound <= length(lat)
	lat_bin = lat(bound-window_size+1:bound);
	lon_bin = lat(bound-window_size+1:bound);
	time_bin = lat(bound-window_size+1:bound);
	elev_bin = lat(bound-window_size+1:bound);
	snr_bin = lat(bound-window_size+1:bound);
	dist_bin = lat(bound-window_size+1:bound);
else
	bound = window_size*(window_count-1);
	lat_bin = lat(bound-window_size+1:bound);
	lon_bin = lat(bound-window_size+1:bound);
	time_bin = lat(bound-window_size+1:bound);
	elev_bin = lat(bound-window_size+1:bound);
	snr_bin = lat(bound-window_size+1:bound);
	dist_bin = lat(bound-window_size+1:bound);
end

end

