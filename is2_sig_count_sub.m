function [hicount,medcount,lowcount,bkgcount] = is2_sig_count_sub(elev_bin, class_bin)
% A function to extract photon counts for high/medium/low confidence
% signals and background noise

addpath(genpath('/Users/zhfair/Documents/'))

plotting = 0; % Yes=1, No=0

hicount = []; medcount = [];
lowcount = []; bkgcount = [];
for k = 1:12
    [~, ~, elev_chunk,~] = is2_chunking_sub(elev_bin, k);
    [~, ~, class_chunk, ~] = is2_chunking_sub(class_bin,k);
    
    high_conf = elev_chunk;
    med_conf = elev_chunk;
    low_conf = elev_chunk;
    bckgrd = elev_chunk;
    
    high_conf(class_chunk~=4) =  NaN;
    med_conf(class_chunk~=3) = NaN;
    low_conf(class_chunk~=2) = NaN;
    bckgrd(class_chunk~=0) = NaN;
    
    ints = 1:23456;
    for i = 1:83
        high(i) = length(high_conf(~isnan(high_conf(ints>(i-1)*285 & ints<=i*285))));
        medium(i) = length(med_conf(~isnan(med_conf(ints>(i-1)*285 & ints<=i*285))));
        low(i) = length(low_conf(~isnan(low_conf(ints>(i-1)*285 & ints<=i*285))));
        bkg(i) = length(bckgrd(~isnan(bckgrd(ints>(i-1)*285 & ints<=i*285))));
        
        if high(i)>=81 && high(i)<87
            high(i) = NaN;
        end
    end
        
        
    
    hicount = [hicount high];
    medcount = [medcount medium];
    lowcount = [lowcount low];
    bkgcount = [bkgcount bkg];
    
end