clear; clc;

addpath(genpath('/Users/zhfair/Documents/IS2_data'))

fname = 'ATL03_20181225185948_13460110_001_01.h5';

for i = 1:3
    track_str1 = ['/gt', num2str(i), 'r/heights/'];
    track_str2 = ['/gt', num2str(i), 'l/heights/'];
    
    elev_weak = h5read(fname, [track_str1, '/h_ph']);
    elev_strong = h5read(fname, [track_str2, '/h_ph']);
    
    class_weak = h5read(fname, [track_str1, '/signal_conf_ph']);
    class_strong = h5read(fname, [track_str2, '/signal_conf_ph']);
    class_weak_consol = is2_class_merge(class_weak);
    class_strong_consol = is2_class_merge(class_strong);
    
    high_weak = elev_weak; high_weak(class_weak_consol~=4) = NaN;
    high_strong = elev_strong; high_strong(class_strong_consol~=4) = NaN;
    med_weak = elev_weak; med_weak(class_weak_consol~=3) = NaN;
    med_strong = elev_strong; med_strong(class_strong_consol~=3) = NaN;
    low_weak = elev_weak; low_weak(class_weak_consol~=2) = NaN;
    low_strong = elev_strong; low_strong(class_strong_consol~=2) = NaN;
    noise_weak = elev_weak; noise_weak(class_weak_consol~=0) = NaN;
    noise_strong = elev_strong; noise_strong(class_strong_consol~=0) = NaN;
    
    lon_weak = h5read(fname, [track_str1, '/lon_ph']);
    lon_strong = h5read(fname, [track_str2, '/lon_ph']);
    lat_weak = h5read(fname, [track_str1, '/lat_ph']);
    lat_strong = h5read(fname, [track_str2, '/lat_ph']);
    
    figure; 
    subplot(1,2,1)
    plot(lon_weak,high_weak, '.', 'MarkerSize', 2)
    hold on; plot(lon_weak,med_weak, '.', 'MarkerSize', 2, 'Color', rgb('green'))
    plot(lon_weak,low_weak, '.', 'MarkerSize', 2, 'Color', rgb('rose'))
    plot(lon_weak,noise_weak, '.', 'MarkerSize', 2, 'Color', rgb('sky blue'))
%     xlim([67.27 67.32]); ylim([175 275])
    title('Amery Ice Shelf, Weak Beam')
    xlabel('Longitude', 'FontWeight', 'bold', 'FontSize', 12)
    ylabel('Elevation [m]', 'FontWeight', 'bold', 'FontSize', 12)
    
    subplot(1,2,2)
    plot(lon_strong,high_strong, '.', 'MarkerSize', 2)
    hold on; plot(lon_strong,med_strong, '.', 'MarkerSize', 2, 'Color', rgb('green'))
    plot(lon_strong,low_strong, '.', 'MarkerSize', 2, 'Color', rgb('rose'))
    plot(lon_strong,noise_strong, '.', 'MarkerSize', 2, 'Color', rgb('sky blue'))
%     xlim([67.27 67.32]); ylim([175 275])
    title('Amery Ice Shelf, Strong Beam')
    xlabel('Longitude', 'FontWeight', 'bold', 'FontSize', 12)
    ylabel('Elevation [m]', 'FontWeight', 'bold', 'FontSize', 12)
    
    figure;
    antarctica = shaperead('landareas', 'UseGeoCoords', true,...
            'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'}); 
    worldmap('antarctica')
    patchm(antarctica.Lat, antarctica.Lon, [0.5 1 0.5])
    hold on;
    plotm(lat_weak,lon_weak, 'LineWidth',3)
    plotm(lat_strong,lon_strong, 'LineWidth',3)
    pause; close all
    
end