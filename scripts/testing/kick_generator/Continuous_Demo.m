%% Continuous demo

% Resetting MATLAB
clear, clc, close all;

% Loading fits
load('Flat_Fits.mat');
P_Flat = P;
load('Kick_Fits.mat');

% Initializing variables
x = linspace(0,1);
X = x;
Hip = 140;
Knee = 130;
End_Slope_Hip = 0;
End_Slope_Knee = 0;
X_Total = [];
Hip_Total = [];
Knee_Total = [];
Hip_Diff_Total = [];
Knee_Diff_Total = [];
checker = false;
T = table(X_Total, Hip_Total, Knee_Total, Hip_Diff_Total, Knee_Diff_Total);
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

% Labeling plot
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Knee Data');
xticks([]);

% Parameter bounds
Length_Bounds = [0.9, 1.5];
smooth_Bounds = [0.3, 0.3];
Flat_Diff_Bounds = [10, 15, 10, 15];
Angle_Start_Bounds = [150, 120, 130, 100];
Angle_Diff_Bounds = [100, 70, 90, 60];
Angle_End_Bounds = [140, 140, 120, 120];
visual = true;
Total_Time = 40;

% Main loop
while checker == false
    
    
    [X, Hip, Knee, End_Slope_Hip, End_Slope_Knee, T] = Continuous_Demo_Plot(X, Hip, Knee, Length_Bounds, smooth_Bounds, Flat_Diff_Bounds, Angle_End_Bounds, P_Flat, true, End_Slope_Hip, End_Slope_Knee, h1, h2, leg, visual, T);
    [X, Hip, Knee, End_Slope_Hip, End_Slope_Knee, T] = Continuous_Demo_Plot(X, Hip, Knee, Length_Bounds, smooth_Bounds, Flat_Diff_Bounds, Angle_End_Bounds, P_Flat, true, End_Slope_Hip, End_Slope_Knee, h1, h2, leg, visual, T);
    [X, Hip, Knee, End_Slope_Hip, End_Slope_Knee, T] = Continuous_Demo_Plot(X, Hip, Knee, Length_Bounds, smooth_Bounds, Angle_Diff_Bounds, Angle_End_Bounds, P, false, End_Slope_Hip, End_Slope_Knee, h1, h2, leg, visual, T);
    
    
%     % Initializing variables
%     X_Start = X(end);
%     isflat = true;
%     Length = rand(1).* (Length_Bounds(2) - Length_Bounds(1)) + Length_Bounds(1);
%     smooth = rand(1).* (smooth_Bounds(2) - smooth_Bounds(1)) + smooth_Bounds(1);
%     
%     % Hip joint
%     Joint = 1;
%     Angle_Start = Hip(end);
%     Angle_Diff = rand(1).* (Flat_Diff_Bounds(2) - Flat_Diff_Bounds(1)) + Flat_Diff_Bounds(1);
%     Angle_End = Angle_Start;
%     [~, Hip, End_Slope_Hip] = Continuous_Demo_Math(P_Flat, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Hip);
%     
%     % Knee joint
%     Joint = 2;
%     Angle_Start = Knee(end);
%     Angle_Diff = rand(1).* (Flat_Diff_Bounds(2) - Flat_Diff_Bounds(1)) + Flat_Diff_Bounds(1);
%     Angle_End = Angle_Start;
%     [X, Knee, End_Slope_Knee] = Continuous_Demo_Math(P_Flat, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Knee);
%     
%     % Plotting flat area
%     for i = 1:length(X)
%         addpoints(h1, X(i), Hip(i));
%         addpoints(h2, X(i), Knee(i));
%         drawnow;
%         xlim([X(i) - 3,X(i) + 1]);
%     end
%     
%     % Storing values for large plot
%     X_Total = [X_Total, X];
%     Hip_Total = [Hip_Total, Hip];
%     Diff_Hip = diff(Hip) ./ diff(X);
%     Knee_Total = [Knee_Total, Knee];
%     Diff_Knee = diff(Knee) ./ diff(X);
%     Hip_Diff_Total = [Hip_Diff_Total, [Diff_Hip, Diff_Hip(end)]];
%     Knee_Diff_Total = [Knee_Diff_Total, [Diff_Knee, Diff_Knee(end)]];
%     
%     
%     
%     % Initializing variables
%     X_Start = X(end);
%     isflat = true;
%     Length = rand(1).* (Length_Bounds(2) - Length_Bounds(1)) + Length_Bounds(1);
%     smooth = rand(1).* (smooth_Bounds(2) - smooth_Bounds(1)) + smooth_Bounds(1);
%     
%     % Hip joint
%     Joint = 1;
%     Angle_Start = Hip(end);
%     Angle_Diff = rand(1).* (Flat_Diff_Bounds(2) - Flat_Diff_Bounds(1)) + Flat_Diff_Bounds(1);
%     Angle_End = Angle_Start;
%     [~, Hip, End_Slope_Hip] = Continuous_Demo_Math(P_Flat, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Hip);
%     
%     % Knee joint
%     Joint = 2;
%     Angle_Start = Knee(end);
%     Angle_Diff = rand(1).* (Flat_Diff_Bounds(2) - Flat_Diff_Bounds(1)) + Flat_Diff_Bounds(1);
%     Angle_End = Angle_Start;
%     [X, Knee, End_Slope_Knee] = Continuous_Demo_Math(P_Flat, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Knee);
%     
%     % Plotting flat area
%     for i = 1:length(X)
%         addpoints(h1, X(i), Hip(i));
%         addpoints(h2, X(i), Knee(i));
%         drawnow;
%         xlim([X(i) - 3,X(i) + 1]);
%     end
%     
%     % Storing values for large plot
%     X_Total = [X_Total, X];
%     Hip_Total = [Hip_Total, Hip];
%     Diff_Hip = diff(Hip) ./ diff(X);
%     Knee_Total = [Knee_Total, Knee];
%     Diff_Knee = diff(Knee) ./ diff(X);
%     Hip_Diff_Total = [Hip_Diff_Total, [Diff_Hip, Diff_Hip(end)]];
%     Knee_Diff_Total = [Knee_Diff_Total, [Diff_Knee, Diff_Knee(end)]];
%     
%     
%     
%     % Initializing variables
%     X_Start = X(end);
%     isflat = false;
%     Length = rand(1).* (Length_Bounds(2) - Length_Bounds(1)) + Length_Bounds(1);
%     smooth = rand(1).* (smooth_Bounds(2) - smooth_Bounds(1)) + smooth_Bounds(1);
%     
%     % Hip joint
%     Joint = 1;
%     Angle_Start = Hip(end);
%     Angle_Diff = rand(1).* (Angle_Diff_Bounds(2) - Angle_Diff_Bounds(1)) + Angle_Diff_Bounds(1);
%     Angle_End = rand(1).* (Angle_End_Bounds(2) - Angle_End_Bounds(1)) + Angle_End_Bounds(1);
%     [~, Hip, End_Slope_Hip] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Hip);
%     
%     % Knee joint
%     Joint = 2;
%     Angle_Start = Knee(end);
%     Angle_Diff = rand(1).* (Angle_Diff_Bounds(2) - Angle_Diff_Bounds(1)) + Angle_Diff_Bounds(1);
%     Angle_End = rand(1).* (Angle_End_Bounds(2) - Angle_End_Bounds(1)) + Angle_End_Bounds(1);
%     [X, Knee, End_Slope_Knee] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Knee);
%     
%     % Plotting kicks
%     for i = 1:length(X)
%         addpoints(h1, X(i), Hip(i));
%         addpoints(h2, X(i), Knee(i));
%         drawnow;
%         xlim([X(i) - 3,X(i) + 1]);
%     end
%     
%     % Storing values for large plot
%     X_Total = [X_Total, X];
%     Hip_Total = [Hip_Total, Hip];
%     Diff_Hip = diff(Hip) ./ diff(X);
%     Knee_Total = [Knee_Total, Knee];
%     Diff_Knee = diff(Knee) ./ diff(X);
%     Hip_Diff_Total = [Hip_Diff_Total, [Diff_Hip, Diff_Hip(end)]];
%     Knee_Diff_Total = [Knee_Diff_Total, [Diff_Knee, Diff_Knee(end)]];
    
    checker = max(T.X_Total) > Total_Time;
    %checker = max(X_Total) > 100;
    
end

% Debug plots
figure(2);
hold on;
plot(T.X_Total, T.Hip_Total);
plot(T.X_Total, T.Hip_Diff_Total);
xlabel('Time (s)');
ylabel('Hip Angle (Degrees)');
title('Infant Hip Angles vs Time');
lgd = legend('Angle Data', 'Velocity Data');
figure(3);
hold on;
plot(T.X_Total, T.Knee_Total);
plot(T.X_Total, T.Knee_Diff_Total);
xlabel('Time (s)');
ylabel('Knee Angle (Degrees)');
title('Infant Knee Angles vs Time');
lgd = legend('Angle Data', 'Velocity Data');
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





