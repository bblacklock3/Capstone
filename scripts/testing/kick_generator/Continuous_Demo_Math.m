

% clear, clc, close all;
% Joint = 1;
% load('Kick_Fits.mat');
%
% % % Percentage of normalized kicks to be smoothing on each side
% % P = ;- Fits coefficients
% Joint = 1; % Hip(1) or Knee(2)
% X_Start = 0;
% % Length - Length of kick in seconds
% % Angle_Start/Min/End - Angle parameters
% Angle_Start = 120;
% Angle_Min = 70;
% Angle_End = 110;
% smooth = 0.05; %Percentage normalized of smoothing on each side of data
% Length = 0.7;

function[x, Angle, End_Slope] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, Init_Slope)



% Setting normalized coords and gathering delta t
x = linspace(0,1,(100.*Length));
xdiff = x(2) - x(1);

% Gathering dataset
r = randi([1,length(P)],1,1);
P_Angle = P{r}(:,Joint);
Angle = P_Angle(1) + P_Angle(2)*cos(x*P_Angle(end)) + P_Angle(3)*sin(x*P_Angle(end)) + P_Angle(4)*cos(2*x*P_Angle(end)) + P_Angle(5)*sin(2*x*P_Angle(end)) + P_Angle(6)*cos(3*x*P_Angle(end)) + P_Angle(7)*sin(3*x*P_Angle(end)) + P_Angle(8)*cos(4*x*P_Angle(end)) + P_Angle(9)*sin(4*x*P_Angle(end)) + P_Angle(10)*cos(5*x*P_Angle(end)) + P_Angle(11)*sin(5*x*P_Angle(end)) + P_Angle(12)*cos(6*x*P_Angle(end)) + P_Angle(13)*sin(6*x*P_Angle(end)) + P_Angle(14)*cos(7*x*P_Angle(end)) + P_Angle(15)*sin(7*x*P_Angle(end)) + P_Angle(16)*cos(8*x*P_Angle(end)) + P_Angle(17)*sin(8*x*P_Angle(end));
Angle_Slope = diff(Angle) ./ diff(x);

% Smoothing beginning of angles
x_smooth = [-smooth:xdiff:0];
n = length(x_smooth);
Angle_Slope_smooth = [Init_Slope:(Angle_Slope(1) - Init_Slope)/n:Angle_Slope(1)];
y_smooth = 0;
Y_smooth = [];
for i = 1:length(Angle_Slope_smooth)
    y_smooth = y_smooth + (Angle_Slope_smooth(i) .* (xdiff));
    Y_smooth = [Y_smooth, y_smooth];
end
Y_smooth = Y_smooth(1:end-1);
Y_smooth = Y_smooth + Angle(1) - Y_smooth(end) + (Angle_Slope(1) .* x_smooth(end));
warning off;
x_p = [x_smooth, x(2:end/10)];
Angle_p = [Y_smooth, Angle(2:end/10)];
warning on;
p = pchip(x_p, Angle_p, x_p);
n = abs(n);
x_p = x_p(1:n);
p = p(1:n);
x = [-smooth(2:end), x_p, x];
Angle = [p, Angle];

% Smoothing end of angles
x_smooth = [1+xdiff:xdiff:1+smooth];
n = length(x_smooth);
Angle_Slope_smooth = [Angle_Slope(end):-Angle_Slope(end)/n:0];
y_smooth = Angle(end);
Y_smooth = [];
for i = 1:length(Angle_Slope_smooth)
    y_smooth = y_smooth + (Angle_Slope_smooth(i) .* (xdiff));
    Y_smooth = [Y_smooth, y_smooth];
end
Y_smooth = Y_smooth(1:end-1);
warning off;
x_p = [x(9*end/10:end), x_smooth];
Angle_p = [Angle(9*end/10:end), Y_smooth];
warning on;
p = pchip(x_p,Angle_p,x_p);
x_p = x_p(end-n+1:end);
p = p(end-n+1:end);
x = [x, x_p(1:end-1)];
Angle = [Angle, p(1:end-1)];

% Formatting angles
if isflat == true
    % Conforming flat angles to input parameters
    Angle = Angle - Angle(1);
    Angle_Diff_T = Angle_Diff;
    Angle_Diff_A = max(abs(max(Angle) - Angle(1)), abs(min(Angle) - Angle(1)));
    Offset = Angle(length(x_smooth)+1);
    Angle(length(x_smooth)+1:end) = Angle(length(x_smooth)+1:end) - Offset;
    Angle(length(x_smooth)+1:end) = Angle(length(x_smooth)+1:end) .* (Angle_Diff_T ./ Angle_Diff_A);
    Angle(length(x_smooth)+1:end) = Angle(length(x_smooth)+1:end) + Offset;
    Angle = Angle + Angle_Start;
    x = x - x(1);
    x = x ./ x(end) .* Length;
    x = x + X_Start;
else
    % Conforming kick angles to input parameters
    Angle = Angle - Angle(1);
    Angle_Diff_T = Angle_Start - Angle_Diff;
    Angle_Diff_A = Angle(1) - min(Angle);
    Angle = Angle .* (Angle_Diff_T ./ Angle_Diff_A);
    Angle = Angle - min(Angle);
    Angle_Diff_T = Angle_End - Angle_Diff;
    Angle_Diff_A = Angle(end) - min(Angle);
    [~,idx] = min(Angle);
    Angle = [Angle(1:idx), Angle(idx+1:end) .* (Angle_Diff_T ./ Angle_Diff_A)];
    Angle = Angle - Angle(1) + Angle_Start;
    x = x - x(1);
    x = x ./ x(end) .* Length;
    x = x + X_Start;
end


Angle_Slope = diff(Angle) ./ diff(x);
End_Slope = Angle_Slope(end);
x(1) = [];
Angle(1) = [];

% % Debug plot
% figure(1);
% plot(x, Angle);

end





