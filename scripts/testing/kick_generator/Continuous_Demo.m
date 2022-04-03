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
h1 = animatedline;
h2 = animatedline;
h1.Color = '#0072BD';
h2.Color = '#D95319';

% Plot settings
figure(1);
ylim([50, 170]);

% Labeling plot
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Knee Data');
xticks([]);

% Parameter bounds
X_Start = X(end);
Length_Bounds = [0.6, 2];
smooth_Bounds = [0.03, 0.1];
Angle_Start_Bounds = [100, 150];
Angle_Min_Bounds = [55, 90];
Angle_End_Bounds = [100, 150];

while true
    
    % Plotting flat area
    for i = 1:length(X)
        addpoints(h1, X(i), Hip(end));
        addpoints(h2, X(i), Knee(end));
        drawnow;
        xlim([X(i) - 3,X(i) + 1]);
        pause(0.001);
    end
    
    % Initializing variables
    X_Start = X(end);
    Length = rand(1).* (Length_Bounds(2) - Length_Bounds(1)) + Length_Bounds(1);
    smooth = rand(1).* (smooth_Bounds(2) - smooth_Bounds(1)) + smooth_Bounds(1);
    
    % Hip joint
    Joint = 1;
    Angle_Start = rand(1).* (Angle_Start_Bounds(2) - Angle_Start_Bounds(1)) + Angle_Start_Bounds(1);
    Angle_Start = Hip(end);
    Angle_Min = rand(1).* (Angle_Min_Bounds(2) - Angle_Min_Bounds(1)) + Angle_Min_Bounds(1);
    Angle_End = rand(1).* (Angle_End_Bounds(2) - Angle_End_Bounds(1)) + Angle_End_Bounds(1);
    [~, Hip] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Min, Angle_End, smooth);
    
    % Knee joint
    Joint = 2;
    Angle_Start = rand(1).* (Angle_Start_Bounds(2) - Angle_Start_Bounds(1)) + Angle_Start_Bounds(1);
    Angle_Start = Knee(end);
    Angle_Min = rand(1).* (Angle_Min_Bounds(2) - Angle_Min_Bounds(1)) + Angle_Min_Bounds(1);
    Angle_End = rand(1).* (Angle_End_Bounds(2) - Angle_End_Bounds(1)) + Angle_End_Bounds(1);
    [X, Knee] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Min, Angle_End, smooth);
    
    % Plotting kicks
    for i = 1:length(X)
        addpoints(h1, X(i), Hip(i));
        addpoints(h2, X(i), Knee(i));
        drawnow;
        xlim([X(i) - 3,X(i) + 1]);
        pause(0.001);
    end
    
    % Iterating time
    X = X + Length;
    
end









% for i = 1:length(P)
%
%     figure(1);
%     hold on;
%     x = linspace(0,1);
%     P_Hip = P{i}(:,1);
%     P_Knee = P{i}(:,2);
%
%     %Hip = P_Hip(1)*sin(P_Hip(2)*x+P_Hip(3)) + P_Hip(4)*sin(P_Hip(5)*x+P_Hip(6)) + P_Hip(7)*sin(P_Hip(8)*x+P_Hip(9)) + P_Hip(10)*sin(P_Hip(11)*x+P_Hip(12)) + P_Hip(13)*sin(P_Hip(14)*x+P_Hip(15)) + P_Hip(16)*sin(P_Hip(17)*x+P_Hip(18)) + P_Hip(19)*sin(P_Hip(20)*x+P_Hip(21)) + P_Hip(22)*sin(P_Hip(23)*x+P_Hip(24));
%     %Knee = P_Knee(1)*sin(P_Knee(2)*x+P_Knee(3)) + P_Knee(4)*sin(P_Knee(5)*x+P_Knee(6)) + P_Knee(7)*sin(P_Knee(8)*x+P_Knee(9)) + P_Knee(10)*sin(P_Knee(11)*x+P_Knee(12)) + P_Knee(13)*sin(P_Knee(14)*x+P_Knee(15)) + P_Knee(16)*sin(P_Knee(17)*x+P_Knee(18)) + P_Knee(19)*sin(P_Knee(20)*x+P_Knee(21)) + P_Knee(22)*sin(P_Knee(23)*x+P_Knee(24));
%
%     Hip = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end)) + P_Hip(8)*cos(4*x*P_Hip(end)) + P_Hip(9)*sin(4*x*P_Hip(end)) + P_Hip(10)*cos(5*x*P_Hip(end)) + P_Hip(11)*sin(5*x*P_Hip(end)) + P_Hip(12)*cos(6*x*P_Hip(end)) + P_Hip(13)*sin(6*x*P_Hip(end)) + P_Hip(14)*cos(7*x*P_Hip(end)) + P_Hip(15)*sin(7*x*P_Hip(end)) + P_Hip(16)*cos(8*x*P_Hip(end)) + P_Hip(17)*sin(8*x*P_Hip(end));
%     Knee = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end)) + P_Knee(8)*cos(4*x*P_Knee(end)) + P_Knee(9)*sin(4*x*P_Knee(end)) + P_Knee(10)*cos(5*x*P_Knee(end)) + P_Knee(11)*sin(5*x*P_Knee(end)) + P_Knee(12)*cos(6*x*P_Knee(end)) + P_Knee(13)*sin(6*x*P_Knee(end)) + P_Knee(14)*cos(7*x*P_Knee(end)) + P_Knee(15)*sin(7*x*P_Knee(end)) + P_Knee(16)*cos(8*x*P_Knee(end)) + P_Knee(17)*sin(8*x*P_Knee(end));
%
%     %     dT = x(2) - x(1);
%     %     Fs = 1/dT;
%     %     Hip = lowpass(Hip, 0.01, Fs);
%     %     Knee = lowpass(Knee, 0.01, Fs);
%
%
%
%     plot(x, Hip, 'b', x, Knee, 'r');
%     hold off;
%
% end



