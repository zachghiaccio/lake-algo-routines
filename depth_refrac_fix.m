function [dist_bin_corrected, depth_corrected, dz] = depth_refrac_fix(dist_bin,signal_depth,theta_inc)
% Subroutine to apply refraction corrections to depths in ICESat-2
% dist_bin is the uncorrected horizontal distance
% signal_depth is the uncorrected water depth; corresponds to sfc_btm_diff_means in autoatl.m
% K is the reference azimuth angle; only needed if not using along-track distance
% theta_inc is the angle of incidence onto the water; defined in autoatl.m

dist_bin_corrected = NaN(length(dist_bin),1);

% Angle of refraction
n_a = 1.00029;
n_w = 1.33469;
theta_ref = asin( (n_a*sin(theta_inc))/n_w );

% Computing triangle geometries
S = signal_depth/cos(theta_inc);
R = (S*n_a)/n_w;
gamma = (pi/2) - theta_inc;

phi = theta_inc - theta_ref;
P = sqrt(R.^2 + S.^2 - 2.*R.*S.*cos(phi));
alfa = asin( (R.*sin(phi))./P );
beta = gamma - alfa;

% Offsets and corrections
dy = P.*cos(beta);
dz = P.*sin(beta);
depth_corrected = signal_depth - dz;

% Apply distance correction over water only
idx = ~isnan(signal_depth);
dist_bin_corrected(idx) = dist_bin(idx)' - dy(idx);

end
