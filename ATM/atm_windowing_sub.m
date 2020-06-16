function [lat_bin,lon_bin,time_bin_hhmmss,elev_bin] = atm_windowing_sub(lat,lon,time_hhmmss,elev,bin_window,bin_count)
% A subroutine for autoatm4.m, designed to separate a data file into
% windows that are viewed iteratively. atm_chunking_sub further separates
% the output of this function into smaller pieces for atm_sfc_detect_sub.m

bound = bin_window*bin_count;
lat_bin = lat(bound-bin_window+1:bound);
lon_bin = lon(bound-bin_window+1:bound);
time_bin_hhmmss = time_hhmmss(bound-bin_window+1:bound);
elev_bin = elev(bound-bin_window+1:bound);

end