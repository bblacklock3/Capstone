function [hip_avg_max, knee_avg_max] = getAvgPeak(interval,h,k,make_plot)
%h_thresh is the velocity threshold for the hip velocity to be counted
h_thresh = 30;
k_thresh = 30;
h_vel = gradient(h,interval);
k_vel = gradient(k,interval);

%h_vel = h;
%k_vel = k;

time = 0:interval:(length(h_vel)-1)*interval;

[h_pks, h_locs] = findpeaks(h_vel,time,'Threshold',0);
mask = h_pks > h_thresh;
h_pks = h_pks(mask);
h_locs = h_locs(mask);

if (make_plot)
    tiledlayout(2,1)
    nexttile
    plot(time,h_vel);
    hold on
    plot(h_locs,h_pks,'g*')
end
hip_avg_max = mean(h_pks);


[k_pks, k_locs] = findpeaks(k_vel,time,'Threshold',0);
mask = k_pks > k_thresh;
k_pks = k_pks(mask);
k_locs = k_locs(mask);
if (make_plot)
    nexttile
    plot(time,k_vel);
    hold on
    plot(k_locs,k_pks,'g*')
end
knee_avg_max = mean(k_pks);
end