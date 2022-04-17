function [h_rom,k_rom] = getROM(h,k,make_plot)
%h is a vector of hip angles
%k is a vector of knee angles
%make_plot is a boolean which decides whether to show the plot
%of time spent at each angle
h_d = sort(round(h));
k_d = sort(round(k));

h_t = [];
k_t = [];

for i = min(h_d):1:max(h_d)
    cur = sum(h_d == i);
    h_t = [h_t cur];
end
h_t(h_t == 0) = [];

for i = min(k_d):1:max(k_d)
    cur = sum(k_d == i);
    k_t = [k_t cur];
end
k_t(k_t == 0) = [];

if (make_plot)
    tiledlayout(2,2)
    nexttile
    bar(unique(h_d),h_t)
    xlabel('Hip Angle (discretized to 1 degree)')
    ylabel('Number of frames at angle')
    nexttile
    bar(unique(h_d),h_t/sum(h_t)*100)
    xlabel('Hip Angle (discretized to 1 degree)')
    ylabel('% of time spent at angle')
    
    nexttile
    bar(unique(k_d),k_t)
    xlabel('Knee Angle (discretized to 1 degree)')
    ylabel('Number of frames at angle')
    nexttile
    bar(unique(k_d),k_t/sum(k_t)*100)
    xlabel('Knee Angle (discretized to 1 degree)')
    ylabel('% of time spent at angle')
end
threshold_bot = .05;
threshold_top = .05;

size_bot = round(length(h_d)*threshold_bot);
size_top = round(length(h_d)*threshold_top);



bot_s = 0;
i = 1;
while (bot_s < size_bot)
    if (bot_s + h_t(i) > size_bot)
        break
    end
    bot_s = bot_s + h_t(i);
    i = i + 1;
end

top_s = 0;
j = length(h_t);
while (top_s < size_top)
    if (top_s + h_t(j) > size_top)
        break
    end
    top_s = top_s + h_t(j);
    j = j - 1;
end
h_un = unique(h_d);
h_top_avg = (sum(h_t(j:end).*h_un(j:end)) + (size_top-top_s)*h_un(j-1))/size_top;
h_bot_avg = (sum(h_t(1:i).*h_un(1:i)) + (size_bot-bot_s)*h_un(i+1))/size_bot;


bot_s = 0;
i = 1;
while (bot_s < size_bot)
    if (bot_s + k_t(i) > size_bot)
        break
    end
    bot_s = bot_s + k_t(i);
    i = i + 1;
end

top_s = 0;
j = length(k_t);
while (top_s < size_top)
    if (top_s + k_t(j) > size_top)
        break
    end
    top_s = top_s + k_t(j);
    j = j - 1;
end
k_un = unique(k_d);
k_top_avg = (sum(k_t(j:end).*k_un(j:end)) + (size_top-top_s)*k_un(j-1))/size_top;
k_bot_avg = (sum(k_t(1:i).*k_un(1:i)) + (size_bot-bot_s)*k_un(i+1))/size_bot;

h_rom = [h_bot_avg, h_top_avg];
k_rom = [k_bot_avg, k_top_avg];
end