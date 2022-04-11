plot_on = false;

%%
load('ContinuousKicks1.mat');
T( any(ismissing(T),2), :) = [];
t = T.X_Total.';
t = t - t(1);
theta_1 = 180-T.Hip_Total.';
theta_2 = 180-T.Knee_Total.';
delta_t = diff(t);
theta_dot_1 = [0 diff(theta_1)./delta_t];
theta_dot_2 = [0 diff(theta_2)./delta_t];
theta_ddot_1 = [0 diff(theta_dot_1)./delta_t];
theta_ddot_2 = [0 diff(theta_dot_2)./delta_t];

hip_deg = [theta_1; theta_dot_1; theta_ddot_1];
knee_deg = [theta_2; theta_dot_2; theta_ddot_2];
knee_deg = hip_deg-knee_deg;

% Plot constants
C = linspecer(6);
red = C(1,:);
blue = C(2,:);
green = C(3,:);
orange = C(4,:);
brown = C(5,:);
yellow = C(6,:);

black = [0 0 0];
font_size = 10;
text_spacing = 0.2;
kinematics_box = [10 3];
dynamics_box = [10 6];
% kinematics_box = [4 1];
% dynamics_box = [7 3];

% Plot kinematics
if plot_on
    figure(10)
    clf, hold on
    plot(t,hip_deg(1,:))
    plot(t,knee_deg(1,:))
    colors = {red; blue;};
    plt_pos = Plot();
    plt_pos.Title = 'Typical Kick: Absolute Position';
    plt_pos.XLim = [min(t) max(t)];
    % plt_pos.XLabel = 'Time (s)';
    plt_pos.YLabel = 'Position (Deg)';
    plt_pos.BoxDim = kinematics_box;
    plt_pos.FontSize = font_size;
    plt_pos.Colors = colors;
    plt_pos.LineWidth = [2, 2];
    plt_pos.LineStyle = {'-','-'};
    plt_pos.Legend = {'Hip','Knee'};
    plt_pos.LegendLoc = 'bestoutside';
    plt_pos.Resolution = 300;
end

%%
if plot_on
    figure(11)
    clf, hold on
    plot(t,hip_deg(2,:)*pi/180)
    plot(t,knee_deg(2,:)*pi/180)
    colors = {red; blue;};
    plt_vel = Plot();
    plt_vel.Title = 'Typical Kick: Absolute Velocity';
    plt_vel.XLim = [min(t) max(t)];
    % plt_vel.XLabel = 'Time (s)';
    plt_vel.YLabel = 'Velocity (rev/s)';
    plt_vel.BoxDim = kinematics_box;
    plt_vel.FontSize = font_size;
    plt_vel.Colors = colors;
    plt_vel.LineWidth = [2, 2];
    plt_vel.LineStyle = {'--','--'};
    plt_vel.Legend = {'Hip','Knee'};
    plt_vel.LegendTextColor = [0 0 0];
    plt_vel.LegendLoc = 'bestoutside';
    plt_vel.Resolution = 300;
end

if plot_on
    figure(12)
    clf, hold on
    plot(t,hip_deg(3,:)*pi/180)
    plot(t,knee_deg(3,:)*pi/180)
    colors = {red; blue;};
    plt_accel = Plot();
    plt_accel.Title = 'Typical Kick: Absolute Acceleration';
    plt_accel.XLim = [min(t) max(t)];
    plt_accel.XLabel = 'Time (s)';
    plt_accel.YLabel = 'Acceleration (rev/s^2)';
    plt_accel.BoxDim = kinematics_box;
    plt_accel.FontSize = font_size;
    plt_accel.Colors = colors;
    plt_accel.LineWidth = [2, 2];
    plt_accel.LineStyle = {'-.','-.'};
    plt_accel.Legend = {'Hip','Knee'};
    plt_accel.LegendLoc = 'bestoutside';
    plt_accel.Resolution = 300;
end

hip_rad = hip_deg*pi/180;
knee_rad = knee_deg*pi/180;

psi = [hip_rad(1,:); knee_rad(1,:)];
psi_dot = [hip_rad(2,:); knee_rad(2,:)];
psi_ddot = [hip_rad(3,:); knee_rad(3,:)];

%% Inverse dynamics
prop = load_prop('prototype');
prop.max_hip_torque = 0.4;
prop.max_knee_torque = 0.15;
tau_i_1 = [];
tau_i_2 = [];
for i = 1:length(psi)
    tau = infant_id(psi(:,i),psi_dot(:,i),psi_ddot(:,i),prop);
    tau_i_1 = [tau_i_1; tau(1,:)];
    tau_i_2 = [tau_i_2; tau(2,:)];
end
tau_i = [tau_i_1(:,1) tau_i_2(:,1)];

%% Run simulation
psi_sol_arr = zeros([length(t) 4]);
tau_d = zeros([length(t) 2]);

opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9,'MaxStep',0.001);
psi_sol = [psi(1,1) psi_dot(1,1) psi(2,1) psi_dot(2,1)];
for i = 1:length(t)-1
    tau_G = gravity_comp_main([psi_sol(1); psi_sol(3)],prop);


    tau_d(i,:) = ([0.3; 0.3].*tau_G).';

    %     tau_sim = [tau_d(i,1) tau_d(i,2) tau_i(i,1) tau_i(i,2)];
    tau_sim = [0 0 tau_i(i,1) tau_i(i,2)];
    [t_sol,psi_sol] = ode45(@(t,psi) full_system_ode(t,psi,tau_sim,prop),[t(i) t(i+1)],psi_sol,opts_1);
    psi_sol = psi_sol(end,:);
    % Endstop collision
    [j1_clamped,theta_1_lim] = clamp(psi_sol(1),0,pi/2);
    [j2_clamped, theta_2_lim] = clamp(psi_sol(3) - theta_1_lim,-150*pi/180,0);
    psi_sol(1) = theta_1_lim;
    psi_sol(3) = theta_1_lim+theta_2_lim;
    c_r = 0.8;
    if j1_clamped
        psi_sol(2) = -c_r*psi_sol(2);
    end
    if j2_clamped
        psi_sol(4) = -c_r*psi_sol(4);
    end
    psi_sol_arr(i,:) = psi_sol;

    if mod(i,100) == 0
        figure(2)
        subplot(1,3,1)
        plot(t,psi_sol_arr(:,1)*180/pi,'-r',t,psi_sol_arr(:,3)*180/pi,'-b')
        axis tight, grid on, box on
        title('Position vs Time')
        xlabel('Time (s)')
        ylabel('Angular Position (deg)')
        legend('Hip','Knee')
        subplot(1,3,2)
        plot(t,psi_sol_arr(:,2)*180/pi,'--r',t,psi_sol_arr(:,4)*180/pi,'--b')
        axis tight, grid on, box on
        title('Velocity vs Time')
        xlabel('Time (s)')
        ylabel('Angular Veloctiy (deg/s)')
        legend('Hip','Knee')
        subplot(1,3,3)
        plot(t,tau_d(:,1),'-.r',t,tau_d(:,2),'-.b')
        axis tight, grid on, box on
        title('Torque vs Time')
        xlabel('Time (s)')
        ylabel('Torque (Nm)')
        legend('Hip','Knee')
    end
end

%%
figure(20)
clf, hold on
plot(t,tau_i_1(:,1))
plot(t,tau_i_1(:,2))
plot(t,tau_i_1(:,3))
plot(t,tau_i_1(:,4))
colors = {black; red; green; blue;};
plt_min_hip_infant = Plot();
plt_min_hip_infant.Title = 'Infant Torques';
plt_min_hip_infant.TickLength = [0.015 0.015];
plt_min_hip_infant.YLim = [-0.5 0.65];
plt_min_hip_infant.XLabel = 'Time (s)';
plt_min_hip_infant.YLabel = 'Torque (Nm)';
plt_min_hip_infant.BoxDim = dynamics_box;
plt_min_hip_infant.FontSize = font_size;
plt_min_hip_infant.Colors = colors;
plt_min_hip_infant.LineWidth = [3, 2, 2, 2];
plt_min_hip_infant.LineStyle = {'-','-','-','-'};
plt_min_hip_infant.Legend = {'Net','Mass','Cent.','Grav.'};
plt_min_hip_infant.LegendBox = 'on';
x_start = 0.4; y_start = -0.13; y_space = diff(plt_min_hip_infant.YLim)*text_spacing/plt_min_hip_infant.BoxDim(2);
% text(x_start,y_start,0,sprintf('Avg Net = %0.3f (Nm)',min_avg_i_1(1)),'FontSize',font_size,'Color',colors{1})
% text(x_start,y_start-y_space,0,sprintf('Avg Mass = %0.3f (Nm)',min_avg_i_1(2)),'FontSize',font_size,'Color',colors{2})
% text(x_start,y_start-2*y_space,0,sprintf('Avg Cent. = %0.3f (Nm)',min_avg_i_1(3)),'FontSize',font_size,'Color',colors{3})
% text(x_start,y_start-3*y_space,0,sprintf('Avg Grav. = %0.3f (Nm)',min_avg_i_1(4)),'FontSize',font_size,'Color',colors{4})
plt_min_hip_infant.Resolution = 300;


function [clamped,new_val] = clamp(val,min_val,max_val)
new_val = min(max(val,min_val),max_val);
clamped = val ~= new_val;
end
