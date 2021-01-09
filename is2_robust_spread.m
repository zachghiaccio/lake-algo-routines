function [robust_spread,z_25,z_75] = is2_robust_spread(elev_segment)
% A quick function designed to find the statistical spread of a segment
% of ATL03 high-confidence photons. If the photon distribution is near-Gaussian,
% then the spread is functionally the same as standard deviation.

% Sorting the data by height
N = sum(~isnan(elev_segment));
elev_seg_filt = elev_segment(~isnan(elev_segment));
elev_seg_sort = sort(elev_seg_filt);
idx = 0.5:1:N;

% Spread calculation
z_25 = max(elev_seg_sort(idx<0.25*N));
z_75 = min(elev_seg_sort(idx>0.75*N));
robust_spread = (z_75-z_25)/0.6745;

end
