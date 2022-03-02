% Angle Plot

% Resetting MATLAB
%clear, clc, close;

% Obtaining data
Data = xlsread('Angle_Data.xlsx');
T_Hip = Data(:,1);
T_Hip = T_Hip(~isnan(T_Hip));
Hip = Data(:,2);
Hip = Hip(~isnan(Hip));
T_Ankle = Data(:,4);
T_Ankle = T_Ankle(~isnan(T_Ankle));
Ankle = Data(:,5);
Ankle = Ankle(~isnan(Ankle));
T_Knee = Data(:,7);
T_Knee = T_Knee(~isnan(T_Knee));
Knee = Data(:,8);
Knee = Knee(~isnan(Knee));

%Standardizing time
T_min = min([T_Hip(1), T_Ankle(1), T_Knee(1)]);
T_max = max([T_Hip(end), T_Ankle(end), T_Knee(end)]);

T_Hip = [T_min; T_Hip; T_max];
T_Hip = T_Hip - T_Hip(1);
Hip = [Hip(1); Hip; Hip(end)];

T_Ankle = [T_min; T_Ankle; T_max];
T_Ankle = T_Ankle - T_Ankle(1);
Ankle = [Ankle(1); Ankle; Ankle(end)];

T_Knee = [T_min; T_Knee; T_max];
T_Knee = T_Knee - T_Knee(1);
Knee = [Knee(1); Knee; Knee(end)];

% Creating poly fits (n = degree)
n = 17;
warning('off');
p1 = polyfit(T_Hip, Hip, n);
x1 = linspace(min(T_Hip), max(T_Hip));
y1 = polyval(p1, x1);

p2 = polyfit(T_Ankle, Ankle, n);
x2 = linspace(min(T_Ankle), max(T_Ankle));
y2 = polyval(p2, x2);

p3 = polyfit(T_Knee, Knee, n);
x3 = linspace(min(T_Knee), max(T_Knee));
y3 = polyval(p3, x3);
warning('on');

% Plotting
colors1 = {[0 0.4470 0.7410], [0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880]};
colors2 = {[0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560], [0.6350 0.0780 0.1840]};

figure(1);
hold on;
plot(T_Hip, Hip, '-', 'Color', colors1{1});
plot(x1, y1, '--', 'Color', colors2{1});
plot(T_Ankle, Ankle, '-', 'Color', colors1{2});
plot(x2, y2, '--', 'Color', colors2{2});
plot(T_Knee, Knee, '-', 'Color', colors1{3});
plot(x3, y3, '--', 'Color', colors2{3});

% Labeling plot
xlabel('Frames (100/sec)');
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Hip Fit', 'Ankle Data', 'Ankle Fit', 'Knee Data', 'Knee Fit');
lgd.Location = 'southeast';
