function [avg_freq] = getAvgFreq(interval,h,k)
h_thresh = 0;
k_thresh = 0;
make_plot = true;

time = 0:interval:(length(h)-1)*interval;

[h_pks, h_locs] = findpeaks(h,time,'Threshold',.1);
mask = h_pks > h_thresh;
h_pks = h_pks(mask);
h_locs = h_locs(mask);

if (make_plot)
    tiledlayout(2,1)
    nexttile
    plot(time,h);
    hold on
    plot(h_locs,h_pks,'g*')
end
h_pers = diff(h_locs);
h_freq = mean(1./h_pers);


[k_pks, k_locs] = findpeaks(k,time,'Threshold',.1);
mask = k_pks > k_thresh;
k_pks = k_pks(mask);
k_locs = k_locs(mask);
if (make_plot)
    nexttile
    plot(time,k);
    hold on
    plot(k_locs,k_pks,'g*')
end
k_pers = diff(k_locs);
k_freq = mean(1./k_pers);

avg_freq = (k_freq*length(k_pers) + h_freq*length(h_pers))/(length(k_pers)+length(h_pers));
end