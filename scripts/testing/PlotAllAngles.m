%% Plot all Angles

% Resetting MATLAB
clear, clc, close;

% Data set 1
Data = readmatrix('S1_fig2_both.xlsx');
T_Hip = 1;
Hip = 2;
T_Knee = 7;
Knee = 8;
Bounds = [1.5, 2.8];
Fig_Num = 1;
n = 10;
[p_Hip, p_Knee] = HipKneePlots(Data, T_Hip, Hip, T_Knee, Knee, Bounds, Fig_Num, n);














