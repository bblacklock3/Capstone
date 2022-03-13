%% Automatically generates plots of n random kicks

% Resetting MATLAB
clear, clc, close all;

Starting = [140 140]; %[120 160]
Ending = [140 140]; %[120 160]
Minimum = [80 80]; %[60 100]
Length = 1.1; %[0.6, 1.5]

for i = 1:1
[Time, Hip, Knee] = Kick_Generator_Math(Starting, Ending, Minimum, Length);

figure(1);
hold on;
plot(Time, Hip, 'b');
plot(Time, Knee, 'r');

% Labeling plot
xlabel('Time (sec)');
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Knee Data');
xlim([min(Time), max(Time)]);

end