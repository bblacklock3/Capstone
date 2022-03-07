function [p_Hip, p_Knee] = HipKneePlots(Data, T_Hip, Hip, T_Knee, Knee, Bounds, Fig_Num, n)

% Obtaining data
T_Hip = Data(:,T_Hip);
T_Hip = T_Hip(~isnan(T_Hip));
Hip = Data(:,Hip);
Hip = Hip(~isnan(Hip));
T_Knee = Data(:,T_Knee);
T_Knee = T_Knee(~isnan(T_Knee));
Knee = Data(:,Knee);
Knee = Knee(~isnan(Knee));

% Trimming data
mask = T_Hip > Bounds(1) & T_Hip < Bounds(2);
T_Hip = T_Hip(mask);
Hip = Hip(mask);
mask = T_Knee > Bounds(1) & T_Knee < Bounds(2);
T_Knee = T_Knee(mask);
Knee = Knee(mask);

% Min and max time/angles
T_min = min([T_Hip(1), T_Knee(1)]);
T_max = max([T_Hip(end), T_Knee(end)]);
Angle_min = min([min(Hip), min(Knee)]);
Angle_max = max([max(Hip), max(Knee)]);

% Standardizing time
T_Hip = [T_min; T_Hip; T_max];
T_Hip = T_Hip - T_Hip(1);
Hip = [Hip(1); Hip; Hip(end)];
T_Knee = [T_min; T_Knee; T_max];
T_Knee = T_Knee - T_Knee(1);
Knee = [Knee(1); Knee; Knee(end)];

% Min and max time/angles
T_min = min([T_Hip(1), T_Knee(1)]);
T_max = max([T_Hip(end), T_Knee(end)]);
Angle_min = min([min(Hip), min(Knee)]);
Angle_max = max([max(Hip), max(Knee)]);

% Creating poly fits (n = degree)
warning('off');
p_Hip = polyfit(T_Hip, Hip, n);
x1 = linspace(min(T_Hip), max(T_Hip));
y1 = polyval(p_Hip, x1);
p_Knee = polyfit(T_Knee, Knee, n);
x2 = linspace(min(T_Knee), max(T_Knee));
y2 = polyval(p_Knee, x2);
warning('on');

% Plotting
colors1 = {[0 0.4470 0.7410], [0.9290 0.6940 0.1250]};
colors2 = {[0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560]};
figure(Fig_Num);
hold on;
plot(T_Hip, Hip, '-', 'Color', colors1{1});
plot(x1, y1, '--', 'Color', colors2{1});
plot(T_Knee, Knee, '-', 'Color', colors1{2});
plot(x2, y2, '--', 'Color', colors2{2});

% Labeling plot
xlabel('Time (s)');
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Hip Fit', 'Knee Data', 'Knee Fit');
lgd.Location = 'southeast';
xlim([T_min, T_max]);
ylim([Angle_min - 10, Angle_max + 10]);

end