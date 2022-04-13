%% Continuous demo

% Resetting MATLAB
clear, clc, close all;

%% Changeable variables
Length_Bounds = [0.9, 1.5]; % Range of kick times
smooth_Bounds = [0.1, 0.1]; % Ratio of smoothing to data
% For following bounds first two are for hip, last two are for knee
Flat_Diff_Bounds = [10, 15, 10, 15]; % How much flat area can deviate
Angle_End_Bounds = [150, 135, 140, 125]; % Starting/ending angle for kicks
Angle_Diff_Bounds = [90, 75, 80, 65]; % Minimum angle for kicks
Vel_Bounds = 400; % Maximum allowable velocity
Accel_Bounds = 2000; % Maximum allowable acceleration
visual = false; % true for animation, false for no animation
Total_Time = 60; % Total time of data collection

%% Setup
% Loading fits
load('Flat_Fits.mat');
P_Flat = P;
load('Kick_Fits.mat');

% Initializing variables
x = linspace(0,1);
X = x;
Hip = mean(Angle_End_Bounds([1,2]));
Knee = mean(Angle_End_Bounds([3,4]));
End_Slope_Hip = 0;
End_Slope_Knee = 0;
X_Total = [];
Hip_Total = [];
Knee_Total = [];
Hip_Diff_Total = [];
Knee_Diff_Total = [];
checker = false;
T = table(X_Total, Hip_Total, Knee_Total, Hip_Diff_Total, Knee_Diff_Total);

% Setting up animated plot
figure(1);
subplot(1,2,1);
h1 = animatedline;
h2 = animatedline;
h1.Color = '#0072BD';
h2.Color = '#D95319';
subplot(1,2,2);
leg = animatedline('MaximumNumPoints', 3);
leg.LineWidth = 3;
xlim([-0.5,1.75]);
ylim([-0.75,1.75]);
xticks([]);
yticks([]);
title('Infant Joint Positions over Time');
subplot(1,2,1);
ylim([50, 180]);
set(gcf, 'Position', [500, 500, 1000, 400]);
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Knee Data');
xticks([]);

%% Data collection
while checker == false
    
    % Data collection main functions
    [X, Hip, Knee, End_Slope_Hip, End_Slope_Knee, T] = Continuous_Demo_Plot(X, Hip, Knee, Length_Bounds, smooth_Bounds, Flat_Diff_Bounds, Angle_End_Bounds, P_Flat, true, End_Slope_Hip, End_Slope_Knee, h1, h2, leg, visual, T);
    [X, Hip, Knee, End_Slope_Hip, End_Slope_Knee, T] = Continuous_Demo_Plot(X, Hip, Knee, Length_Bounds, smooth_Bounds, Flat_Diff_Bounds, Angle_End_Bounds, P_Flat, true, End_Slope_Hip, End_Slope_Knee, h1, h2, leg, visual, T);
    [X, Hip, Knee, End_Slope_Hip, End_Slope_Knee, T] = Continuous_Demo_Plot(X, Hip, Knee, Length_Bounds, smooth_Bounds, Angle_Diff_Bounds, Angle_End_Bounds, P, false, End_Slope_Hip, End_Slope_Knee, h1, h2, leg, visual, T);
    
    % Continues until over max time
    checker = max(T.X_Total) > Total_Time .* 1.05;
    
end

%% Data processing

% Deleting duplicate time values
[n, bin] = histc(T.X_Total, unique(T.X_Total));
multiple = find(n > 1);
index = find(ismember(bin, multiple));
B = [];
for i=1:length(index)
    if (mod(i,2) ~= 0)
        B = vertcat(B,index(i,:));
    end
end
T(B,:) = [];

% Manual lowpass filter
alpha = 0.98;
for i = 1:length(T.X_Total)-1
    T.Hip_Total(i+1) = (alpha .* T.Hip_Total(i)) + ((1 - alpha) .* T.Hip_Total(i+1));
    T.Knee_Total(i+1) = (alpha .* T.Knee_Total(i)) + ((1 - alpha) .* T.Knee_Total(i+1));
end

% Debug plot
figure(4);
hold on;
plot(T.X_Total, T.Hip_Total);
plot(T.X_Total, T.Hip_Diff_Total);
plot(T.X_Total, T.Knee_Total);
plot(T.X_Total, T.Knee_Diff_Total);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Angle Data', 'Hip Velocity Data', 'Knee Angle Data', 'Knee Velocity Data');

% Zeroing time values
T.X_Total = T.X_Total - T.X_Total(1);

% Eliminating extreme velocities
[Time, Hip, Knee] = Continuous_Demo_Trim(T.X_Total, T.Hip_Total, T.Knee_Total, Vel_Bounds, Accel_Bounds);

% Finding velocities
Hip_Vel = diff(Hip) ./ diff(Time);
Hip_Vel = [Hip_Vel; Hip_Vel(end)];
Knee_Vel = diff(Knee) ./ diff(Time);
Knee_Vel = [Knee_Vel; Knee_Vel(end)];

% Debug plot
figure(5);
hold on;
plot(Time, Hip);
plot(Time, Hip_Vel);
plot(Time, Knee);
plot(Time, Knee_Vel);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Angle Data', 'Hip Velocity Data', 'Knee Angle Data', 'Knee Velocity Data');

% Fitting data
xq1 = [0:0.01:Total_Time];
P_Hip = pchip(Time, Hip, xq1);
Vel_P_Hip = diff(P_Hip) ./ diff(xq1);
Vel_P_Hip = [Vel_P_Hip, Vel_P_Hip(end)];
Accel_P_Hip = diff(Vel_P_Hip) ./ diff(xq1);
Accel_P_Hip = [Accel_P_Hip, Accel_P_Hip(end)];
P_Knee = pchip(Time, Knee, xq1);
Vel_P_Knee = diff(P_Knee) ./ diff(xq1);
Vel_P_Knee = [Vel_P_Knee, Vel_P_Knee(end)];
Accel_P_Knee = diff(Vel_P_Knee) ./ diff(xq1);
Accel_P_Knee = [Accel_P_Knee, Accel_P_Knee(end)];
T = table(xq1', P_Hip', P_Knee');

% Plotting fitted data
figure(6);
hold on;
plot(xq1, P_Hip);
plot(xq1, Vel_P_Hip);
plot(xq1, P_Knee);
plot(xq1, Vel_P_Knee);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Angle Data', 'Hip Velocity Data', 'Knee Angle Data', 'Knee Velocity Data');

% Plotting acceleration fitted data
figure(7);
hold on;
plot(xq1, Accel_P_Hip);
plot(xq1, Accel_P_Knee);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees/s^2)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Acceleration Data', 'Knee Acceleration Data');

% Final position plot
figure(8);
hold on;
plot(xq1, P_Hip);
plot(xq1, P_Knee);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Knee Data');




x_new = linspace(xq1(1),xq1(end),xq1(end).*20);
Hip_New = pchip(xq1, P_Hip, x_new);
Knee_New = pchip(xq1, P_Knee, x_new);

x_new2 = linspace(x_new(1),x_new(end),x_new(end).*100);
Hip_New = spline(x_new, Hip_New, x_new2);
Knee_New = spline(x_new, Knee_New, x_new2);
x_new = x_new2;

x_new2 = linspace(x_new(1),x_new(end),x_new(end).*200);
Hip_New = pchip(x_new, Hip_New, x_new2);
Knee_New = pchip(x_new, Knee_New, x_new2);
x_new = x_new2;

Hip_Velocity_New = gradient(Hip_New) ./ gradient(x_new);
Hip_Acceleration_New = gradient(Hip_Velocity_New) ./ gradient(x_new);
Knee_Velocity_New = gradient(Knee_New) ./ gradient(x_new);
Knee_Acceleration_New = gradient(Knee_Velocity_New) ./ gradient(x_new);

% Get x_new, Hip_New, Knee_New as mat file

% Final position plot 2
figure(9);
hold on;
plot(x_new, Hip_New);
plot(x_new, Knee_New);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data new', 'Knee Data new');

% Final velocity plot 2
figure(10);
hold on;
plot(x_new, Hip_Velocity_New);
plot(x_new, Knee_Velocity_New);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees/s)');
title('Infant Joint Angles Velocities vs Time');
lgd = legend('Hip Velocity new', 'Knee Velocity new');

% Final acceleration plot 2
figure(11);
hold on;
plot(x_new, Hip_Acceleration_New);
plot(x_new, Knee_Acceleration_New);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees/s^2)');
title('Infant Joint Angles Accelerations vs Time');
lgd = legend('Hip Acceleration new', 'Knee Acceleration new');



[xhip, yhip, xknee, yknee, xfoot, yfoot] = Plot_Position(xq1, P_Hip, P_Knee);

n=3;
xq1(n:n:end) = [];
P_Hip(n:n:end) = [];
P_Knee(n:n:end) = [];
xhip(n:n:end) = [];
yhip(n:n:end) = [];
xknee(n:n:end) = [];
yknee(n:n:end) = [];
xfoot(n:n:end) = [];
yfoot(n:n:end) = [];

% Plotting flat area
if visual == true
    for i = 1:length(xq1)
        
        % Angle plot
        addpoints(h1, xq1(i), P_Hip(i));
        addpoints(h2, xq1(i), P_Knee(i));
        drawnow limitrate;
        figure(1);
        subplot(1,2,1);
        xlim([xq1(i) - 3,xq1(i) + 1]);
        
        % Position plot
        addpoints(leg, xhip(i), yhip(i));
        addpoints(leg, xknee(i), yknee(i));
        addpoints(leg, xfoot(i), yfoot(i));
        
    end
end

