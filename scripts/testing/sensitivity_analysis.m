%% Load constansts
prop = load_prop('full_system');

syms g L_1 L_2
syms L__i_g1 L__i_g2 m__i_1 m__i_2 I__i_1 I__i_2
syms psi_1 psi_dot_1 psi_ddot_1 psi_2 psi_dot_2 psi_ddot_2
syms L__d_g1 L__d_g2 m__d_1 m__d_2 I__d_1 I__d_2 I__d_m1 n_1 I__d_m2 n_2

% State varialbes: u
u_psi = [psi_1 psi_dot_1 psi_ddot_1 psi_2 psi_dot_2 psi_ddot_2];
u_const = [g L_1 L_2 L__i_g1 L__i_g2 m__i_1 m__i_2 I__i_1 I__i_2];
u_const_0 = [prop.g prop.L_1 prop.L_2 prop.L__i_g1 prop.L__i_g2 prop.m__i_1 prop.m__i_2 prop.I__i_1 prop.I__i_2];

% Design variables: x
x1 = [L__d_g1 L__d_g2 m__d_1 m__d_2 I__d_1 I__d_m1 n_1];
x1_0 = [prop.L__d_g1 prop.L__d_g2 prop.m__d_1 prop.m__d_2 prop.I__d_1 prop.I__d_m1 prop.n_1];
x1_names = {'L__d_g1' 'L__d_g2' 'm__d_1' 'm__d_2' 'I__d_1' 'I__d_m1' 'n_1'};

x2 = [L__d_g2 m__d_2 I__d_2 I__d_m2 n_2];
x2_0 = [prop.L__d_g2 prop.m__d_2 prop.I__d_2 prop.I__d_m2 prop.n_2];
x2_names = {'L__d_g2' 'm__d_2' 'I__d_2' 'I__d_m2' 'n_2'};

% Full dynamics
full_eqn = full_system_eom;
J1 = full_eqn(1,:);
J2 = full_eqn(2,:);

% Compute gradient
for i = 1:length(x1)
    grad_J1(i,1) = diff(J1,x1(i));
end
for i = 1:length(x2)
    grad_J2(i,1) = diff(J2,x2(i));
end

grad_J1;
grad_J2;

%% Sensitivity Test Cases
% Resisting gravity when horizontal
u_psi_0 = [0 0 0 0 0 0];
J1_0 = double(subs(J1,[u_psi u_const x1],[u_psi_0 u_const_0 x1_0]));
J2_0 = double(subs(J2,[u_psi u_const x2],[u_psi_0 u_const_0 x2_0]));

grad_J1_0 = double(subs(grad_J1,[u_psi u_const x1],[u_psi_0 u_const_0 x1_0]));
grad_J2_0 = double(subs(grad_J2,[u_psi u_const x2],[u_psi_0 u_const_0 x2_0]));

norm_sens_1 = (x1_0./J1_0).'.*(grad_J1_0);
norm_sens_2 = (x2_0./J2_0).'.*(grad_J2_0);

figure(1)
subplot(1,2,1)
barh(norm_sens_1)
title('Joint 1 Torque Sensitivity: Horizontal Static Case')
yticklabels(x1_names)

subplot(1,2,2)
barh(norm_sens_2)
title('Joint 2 Torque Sensitivity: Horizontal Static Case')
yticklabels(x2_names)

% Import Joint Data
% Obtaining data
Data = xlsread('Angle_Data.xlsx');
T_Hip = Data(:,1);
T_Hip = T_Hip(~isnan(T_Hip))/100;
Hip = Data(:,2);
Hip = Hip(~isnan(Hip));
T_Knee = Data(:,7);
T_Knee = T_Knee(~isnan(T_Knee))/100;
Knee = Data(:,8);
Knee = Knee(~isnan(Knee));

Knee = Knee-180;

%Standardizing time
T_min = min([T_Hip(1), T_Knee(1)]);
T_max = max([T_Hip(end), T_Knee(end)]);

T_Hip = [T_min; T_Hip; T_max];
T_Hip = T_Hip - T_Hip(1);
Hip = [Hip(1); Hip; Hip(end)];

T_Knee = [T_min; T_Knee; T_max];
T_Knee = T_Knee - T_Knee(1);
Knee = [Knee(1); Knee; Knee(end)];

% Creating poly fits (n = degree)
n = 17;
warning('off');
t = linspace(min(T_Hip), max(T_Hip),1000);
p1 = polyfit(T_Hip, Hip, n);
p2 = polyfit(T_Knee, Knee, n);

hip_deg = [polyval(p1,t); polyval(polyder(p1),t); polyval(polyder(polyder(p1)),t)];
knee_deg = [polyval(p2,t); polyval(polyder(p2),t); polyval(polyder(polyder(p2)),t)];
knee_deg = hip_deg+knee_deg;
warning('on');

% Typical kick (average and max sensitivity)
trim_len = 50;
t(:,1:trim_len) = [];
t(:,end-trim_len:end) = [];
hip_deg(:,1:trim_len) = [];
hip_deg(:,end-trim_len:end) = [];
knee_deg(:,1:trim_len) = [];
knee_deg(:,end-trim_len:end) = [];

num_sp = 5;
figure(10)
subplot(num_sp,1,1)
cla, hold on
plot(t,hip_deg(1,:),'-r')
plot(t,knee_deg(1,:),'-b')
title('Joint Position in Absolute Coordinates')
xlabel('Time (s)')
ylabel('Position (Deg)')
xlim([min(t) max(t)])

subplot(num_sp,1,2)
cla, hold on
plot(t,hip_deg(2,:)/360,'--r')
plot(t,knee_deg(2,:)/360,'--b')
title('Joint Velocity in Absolute Coordinates')
xlabel('Time (s)')
ylabel('Velocity (rev/s)')
xlim([min(t) max(t)])

subplot(num_sp,1,3)
cla, hold on
plot(t,hip_deg(3,:)/360,'.-r')
plot(t,knee_deg(3,:)/360,'.-b')
title('Joint Acceleration in Absolute Coordinates')
xlabel('Time (s)')
ylabel('Accel (rev/s^2)')
xlim([min(t) max(t)])

hip_rad = hip_deg*pi/180;
knee_rad = knee_deg*pi/180;

% u_psi_0 = [0 0 0 0 0 0];
u_psi = {psi_1 psi_dot_1 psi_ddot_1 psi_2 psi_dot_2 psi_ddot_2};
u_psi_0 = {hip_rad(1,:)' hip_rad(2,:)' hip_rad(3,:)' knee_rad(1,:)' knee_rad(2,:)' knee_rad(3,:)'};

J1_0 = subs(J1,[u_const x1],[u_const_0 x1_0]);
J1_0 = double(subs(J1_0,u_psi,u_psi_0));
J2_0 = subs(J2,[u_const x2],[u_const_0 x2_0]);
J2_0 = double(subs(J2_0,u_psi,u_psi_0));

figure(10)
subplot(num_sp,1,4)
cla, hold on
plot(t,J1_0,'-r')
plot(t,J2_0,'-b')
title('Joint Torque')
xlabel('Time (s)')
ylabel('Torque (Nm)')
xlim([min(t) max(t)])

J1_0 = abs(J1_0);
J2_0 = abs(J2_0);

del_idx_1 = find(J1_0 < 0.15);
del_idx_2 = find(J2_0 < 0.05);

t1 = t;
t1(del_idx_1) = [];
t2 = t;
t2(del_idx_2) = [];

u_psi_1_trim = u_psi_0;
u_psi_2_trim = u_psi_0;
for i = 1:length(u_psi_0)
    psi_trim_1 = u_psi_1_trim{i};
    psi_trim_1(del_idx_1) = [];
    u_psi_1_trim(i) = {psi_trim_1};
    psi_trim_2 = u_psi_2_trim{i};
    psi_trim_2(del_idx_2) = [];
    u_psi_2_trim(i) = {psi_trim_2};
end
J1_0_trim = J1_0;
J2_0_trim = J2_0;

J1_0_trim(del_idx_1) = [];
J2_0_trim(del_idx_2) = [];

subplot(num_sp,1,5)
cla, hold on
plot(t,J1_0,'-r')
plot(t(del_idx_1),J1_0(del_idx_1),'or')
plot(t,J2_0,'-b')
plot(t(del_idx_2),J2_0(del_idx_2),'ob')
title('Joint Torque Zero Trimmed')
xlabel('Time (s)')
ylabel('Torque (Nm)')
xlim([min(t) max(t)])

grad_J1_0 = subs(grad_J1.',[u_const x1],[u_const_0 x1_0]);
grad_J1_0 = double(subs(grad_J1_0,u_psi,u_psi_1_trim))';
grad_J2_0 = subs(grad_J2.',[u_const x2],[u_const_0 x2_0]);
grad_J2_0 = double(subs(grad_J2_0,u_psi,u_psi_2_trim))';

sens_1 = grad_J1_0;
sens_2 = grad_J2_0;
norm_sens_1 = (x1_0./J1_0_trim).'.*sens_1;
norm_sens_2 = (x2_0./J2_0_trim).'.*sens_2;

figure(12)
subplot(2,1,1)
cla, hold on
for i = 1:size(norm_sens_1,1)
    plot(t1,norm_sens_1(i,:))
end
warning('off')
legend(x1_names)
subplot(2,1,2)
cla, hold on
for i = 1:size(norm_sens_2,1)
    plot(t2,norm_sens_2(i,:))
end
legend(x2_names)
warning('on')

%% 
mean_norm_sens_1 = mean(norm_sens_1,2)
mean_norm_sens_2 = mean(norm_sens_2,2)

max_norm_sens_1 = max(norm_sens_1,[],2)
max_norm_sens_2 = max(norm_sens_2,[],2)

figure(2)
subplot(1,2,1)
barh(mean_norm_sens_1)
title('Joint 1 Torque Sensitivity Mean: Typical Kick')
yticklabels(x1_names)

subplot(1,2,2)
barh(mean_norm_sens_2)
title('Joint 2 Torque Sensitivity Mean: Typical Kick')
yticklabels(x2_names)

%% 
% subplot(2,2,3)
% barh(max_norm_sens_1)
% title('Joint 1 Torque Sensitivity Max: Typical Kick')
% yticklabels(x1_names)
% 
% subplot(2,2,4)
% barh(max_norm_sens_2)
% title('Joint 2 Torque Sensitivity Max: Typical Kick')
% yticklabels(x2_names)

%% 
