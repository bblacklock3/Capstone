%% Generate polynomials for kick data sets

% Resetting MATLAB
clear, clc, close all;

% Initializing poly vector
P = {};
type = 'sin8';
num = 1;

% Data set 1 - Optotrak data
Data = readmatrix('an01#001.abc.csv');
Time = [1:length(Data)]';
Time = Time - Time(1);
Time = Time ./ Time(end);
Hip = Data(:,3);
Knee = Data(:,2);
%c_Hip = fit(Time, Hip, type);
%p_Hip = coeffvalues(c_Hip);
%c_Knee = fit(Time, Knee, type);
%p_Knee = coeffvalues(c_Knee);
%P = [P; [{p_Hip}, {p_Knee}]];

[P_Hip] = Curve_Splitting(Time, Hip, type, num);
[P_Knee] = Curve_Splitting(Time, Knee, type, num);
P = [P; [{P_Hip}, {P_Knee}]];

% Data set 2 - Angle data xlsx 1
Data = readmatrix('Angle_Data.xlsx');
Time_Hip = Data(:,1);
Time_Hip = Time_Hip(~isnan(Time_Hip));
Time_Hip = Time_Hip - Time_Hip(1);
Time_Hip = Time_Hip ./ Time_Hip(end);
Hip = Data(:,2);
Hip = Hip(~isnan(Hip));
Time_Knee = Data(:,7);
Time_Knee = Time_Knee(~isnan(Time_Knee));
Time_Knee = Time_Knee - Time_Knee(1);
Time_Knee = Time_Knee ./ Time_Knee(end);
Knee = Data(:,8);
Knee = Knee(~isnan(Knee));
% c_Hip = fit(Time_Hip, Hip, type);
% p_Hip = coeffvalues(c_Hip);
% c_Knee = fit(Time_Knee, Knee, type);
% p_Knee = coeffvalues(c_Knee);
% P = [P; [{p_Hip}, {p_Knee}]];

[P_Hip] = Curve_Splitting(Time_Hip, Hip, type, num);
[P_Knee] = Curve_Splitting(Time_Knee, Knee, type, num);
P = [P; [{P_Hip}, {P_Knee}]];

% Data set 3 & 4 
Data = readmatrix('S1_fig2_both.xlsx');
Time_Hip = Data(:,1);
Time_Hip = Time_Hip(~isnan(Time_Hip));
Time_Hip = Time_Hip - Time_Hip(1);
Time_Hip = Time_Hip ./ Time_Hip(end);
mask_Hip = Time_Hip < 0.35;
mask2_Hip = Time_Hip > 0.55 & Time_Hip < 0.95;
Hip = Data(:,2);
Hip = Hip(~isnan(Hip));
Time_Knee = Data(:,7);
Time_Knee = Time_Knee(~isnan(Time_Knee));
Time_Knee = Time_Knee - Time_Knee(1);
Time_Knee = Time_Knee ./ Time_Knee(end);
mask_Knee = Time_Knee < 0.35;
mask2_Knee = Time_Knee > 0.55 & Time_Knee < 0.95;
Knee = Data(:,8);
Knee = Knee(~isnan(Knee));
Time_Hip_1 = Time_Hip(mask_Hip);
Time_Hip_1 = Time_Hip_1 ./ Time_Hip_1(end);
Hip_1 = Hip(mask_Hip);
Time_Knee_1 = Time_Knee(mask_Knee);
Time_Knee_1 = Time_Knee_1 ./ Time_Knee_1(end);
Knee_1 = Knee(mask_Knee);
Time_Hip_2 = Time_Hip(mask2_Hip);
Time_Hip_2 = Time_Hip_2 - Time_Hip_2(1);
Time_Hip_2 = Time_Hip_2 ./ Time_Hip_2(end);
Hip_2 = Hip(mask2_Hip);
Time_Knee_2 = Time_Knee(mask2_Knee);
Time_Knee_2 = Time_Knee_2 - Time_Knee_2(1);
Time_Knee_2 = Time_Knee_2 ./ Time_Knee_2(end);
Knee_2 = Knee(mask2_Knee);
% c_Hip = fit(Time_Hip_1, Hip_1, type);
% p_Hip = coeffvalues(c_Hip);
% c_Knee = fit(Time_Knee_1, Knee_1, type);
% p_Knee = coeffvalues(c_Knee);
% P = [P; [{p_Hip}, {p_Knee}]];
% c_Hip = fit(Time_Hip_2, Hip_2, type);
% p_Hip = coeffvalues(c_Hip);
% c_Knee = fit(Time_Knee_2, Knee_2, type);
% p_Knee = coeffvalues(c_Knee);
% P = [P; [{p_Hip}, {p_Knee}]];


[P_Hip] = Curve_Splitting(Time_Hip_1, Hip_1, type, num);
[P_Knee] = Curve_Splitting(Time_Knee_1, Knee_1, type, num);
P = [P; [{P_Hip}, {P_Knee}]];


[P_Hip] = Curve_Splitting(Time_Hip_2, Hip_2, type, num);
[P_Knee] = Curve_Splitting(Time_Knee_2, Knee_2, type, num);
P = [P; [{P_Hip}, {P_Knee}]];

Hip_1_Weights = ones(1,length(Hip_1)) .* 0.5;
Time_Hip_1_Test = [Time_Hip_1; 0.4];
Hip_1_Test = [Hip_1; 130];
Hip_1_Weights = [Hip_1_Weights 1];







