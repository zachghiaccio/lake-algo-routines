function [lat_bin,lon_bin,time_bin,elev_bin,class_bin,dist_bin] = is2_windowing_sub(lat,lon,time,elev,class,dist,window_size,window_count)
% A subroutine for autoatl.m, designed to separate a data file into
% windows that are viewed iteratively. is2_chunking_sub further separates
% the output of this function into smaller pieces for is2_sfc_detect_sub.m

bound = window_size*window_count;

if bound <= length(lat)
    lat_bin =  lat(bound-window_size+1:bound);
    lon_bin =  lon(bound-window_size+1:bound);
    time_bin =  time(bound-window_size+1:bound);
    elev_bin =  elev(bound-window_size+1:bound);
    class_bin =  class(bound-window_size+1:bound);
    dist_bin =  dist(bound-window_size+1:bound);
else
    bound = window_size*(window_count-1); % Change this - it is only for testing!
    lat_bin =  lat(bound-window_size+1:bound);
    lon_bin =  lon(bound-window_size+1:bound);
    time_bin =  time(bound-window_size+1:bound);
    elev_bin =  elev(bound-window_size+1:bound);
    class_bin =  class(bound-window_size+1:bound);
    dist_bin =  dist(bound-window_size+1:bound);
end


end
