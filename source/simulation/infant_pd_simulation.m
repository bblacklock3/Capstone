%% Load Joint Data
Data = xlsread('Angle_Data.xlsx');
T_Hip = Data(:,1);
T_Hip = T_Hip(~isnan(T_Hip))/(100);
Hip = Data(:,2);
Hip = Hip(~isnan(Hip));
T_Knee = Data(:,7);
T_Knee = T_Knee(~isnan(T_Knee))/(100);
Knee = Data(:,8);
Knee = Knee(~isnan(Knee));

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
t = linspace(min(T_Hip), max(T_Hip),100);
p1 = polyfit(T_Hip, Hip, n);
p2 = polyfit(T_Knee, Knee, n);

hip_deg = [polyval(p1,t); polyval(polyder(p1),t); polyval(polyder(polyder(p1)),t)];
knee_deg = [polyval(p2,t); polyval(polyder(p2),t); polyval(polyder(polyder(p2)),t)];
warning('on');

trim_len = round(0.05*length(t));
t(:,1:trim_len) = [];
t(:,end-trim_len:end) = [];
t = t - t(1);
hip_deg(:,1:trim_len) = [];
hip_deg(:,end-trim_len:end) = [];
knee_deg(:,1:trim_len) = [];
knee_deg(:,end-trim_len:end) = [];

% figure(100)
% cla, hold on
% plot(t,theta_1,'-r',t,theta_2,'-b')
% plot(t,psi_1,'-r',t,psi_2,'-b')
% axis tight, grid on, box on
% title('Position vs Time')
% xlabel('Time (s)')
% ylabel('Angular Position (deg)')
% legend('Hip','Knee')

hip_rad = hip_deg*pi/180;
knee_rad = knee_deg*pi/180;

theta = [hip_rad(1,:); knee_rad(1,:)].';
theta_dot = [hip_rad(2,:); knee_rad(2,:)].';
theta_ddot = [hip_rad(3,:); knee_rad(3,:)].';

psi = infant2global(theta);
psi_dot = infant2global(theta_dot);
psi_ddot = infant2global(theta_ddot);

%%
% load('ContinuousKicks1.mat');
% T( any(ismissing(T),2), :) = [];
% t = T.X_Total.';
% t = t - t(1);
% theta_1 = 180-T.Hip_Total.';
% theta_2 = 180-T.Knee_Total.';
% delta_t = diff(t);
% theta_dot_1 = [0 diff(theta_1)./delta_t];
% theta_dot_2 = [0 diff(theta_2)./delta_t];
% theta_ddot_1 = [0 diff(theta_dot_1)./delta_t];
% theta_ddot_2 = [0 diff(theta_dot_2)./delta_t];
% 
% hip_deg = [theta_1; theta_dot_1; theta_ddot_1];
% knee_deg = [theta_2; theta_dot_2; theta_ddot_2];
% knee_deg = hip_deg-knee_deg;
% 
% hip_rad = hip_deg*pi/180;
% knee_rad = knee_deg*pi/180;
% 
% psi = [hip_rad(1,:); knee_rad(1,:)].';
% psi_dot = [hip_rad(2,:); knee_rad(2,:)].';
% psi_ddot = [hip_rad(3,:); knee_rad(3,:)].';

% Inverse dynamics
prop = load_prop('prototype');
prop.max_hip_torque = 2;
prop.max_knee_torque = 2;
tau_i_1 = [];
tau_i_2 = [];
for i = 1:length(psi)
    tau = infant_id(psi(i,:),psi_dot(i,:),psi_ddot(i,:),prop);
    tau_i_1 = [tau_i_1; tau(1,:)];
    tau_i_2 = [tau_i_2; tau(2,:)];
end
tau_i_traj = [tau_i_1(:,1) tau_i_2(:,1)];

figure(30)
subplot(1,3,1)
plot(t,psi(:,1)*180/pi,'-r',t,psi(:,2)*180/pi,'-b')
axis tight, grid on, box on
title('Position vs Time')
xlabel('Time (s)')
ylabel('Angular Position (deg)')
legend('Hip','Knee')
subplot(1,3,2)
plot(t,psi_dot(:,1)*180/pi,'--r',t,psi_dot(:,2)*180/pi,'--b')
axis tight, grid on, box on
title('Velocity vs Time')
xlabel('Time (s)')
ylabel('Angular Veloctiy (deg/s)')
legend('Hip','Knee')
subplot(1,3,3)
cla, hold on
plot(t,tau_i_traj(:,1),'-r',t,tau_i_traj(:,2),'-b')
axis tight, grid on, box on
title('Torque vs Time')
xlabel('Time (s)')
ylabel('Torque (Nm)')
legend('Hip','Knee')

% Run simulation
state = 'off';
flex_start = 50*pi/180;
flex_peak = 75*pi/180;
flex_hip_end = 120*pi/180;
flex_knee_end = 10*pi/180;

flex_assist_hip_start = 0;

hip_scaling_factor = 0;

psi_error_1 = 0;
psi_error_prev_1 = 0;
psi_error_dot_1 = 0;

psi_error_2 = 0;
psi_error_prev_2 = 0;
psi_error_dot_2 = 0;

psi_sol_arr = zeros([length(t) 4]);
tau_i = zeros([length(t) 2]);
tau_d = zeros([length(t) 2]);

dt = mean(diff(t));

finish_time = 0;
interp_dur = 0.15;

psi_sol = [psi(1,1) psi_dot(1,1) psi(1,2) psi_dot(1,2)];
for i = 1:length(t)

    tau_i_G = gravity_comp_infant([psi_sol(1); psi_sol(3)],prop);
    tau_d_G = gravity_comp_device([psi_sol(1); psi_sol(3)],prop);
    k_hip = 0.5;
    k_knee = 0.5;
    
    if (strcmp(state,'off'))
        if (finish_time ~= 0)
            curr_time = toc(finish_time);
            if curr_time < interp_dur
                interp_factor = (curr_time - toc(finish_time))/(interp_dur - toc(finish_time));
            end
        else
            interp_factor = 1;
        end
        interp_factor
        % Fake infant torque for simulation only
        tau_d(i,:) = ([k_hip k_knee].*tau_i_G.') + tau_d_G.';
        tau_i_1 = interp_factor*(0.8*(psi(i,1) - psi_sol(1)) + 0.05*(psi_dot(i,1) - psi_sol(2)));
        tau_i_2 = interp_factor*(0.8*(psi(i,2) - psi_sol(3)) + 0.05*(psi_dot(i,2) - psi_sol(4)));
        tau_i(i,:) = [tau_i_1 tau_i_2];
        if (psi(i,1) < flex_start && psi_sol(2) >= flex_start)
            state = 'flex_assist_1';
            flex_hip_assist_start = psi_sol(1);
        end
    end
    if (strcmp(state,'flex_assist_1'))
        Kp_1 = 0.1;
        Kd_1 = 0.1;
        Kp_2 = 0.1;
        Kd_2 = 0.01;
        % Real infant torque
        tau_i_1 = tau_i_traj(i,1);
        tau_i_2 = tau_i_traj(i,2);
        tau_i(i,:) = [tau_i_1 tau_i_2];
        % Device torque: Only assist knee trying to increase angle
        psi_error_2 = flex_peak - psi_sol(3);
        psi_error_dot_2 = (psi_error_2 - psi_error_prev_2)/dt;
        hip_scaling_factor = (psi_sol(1)-flex_hip_assist_start)/(flex_peak - flex_hip_assist_start);
        
        psi_error_1 = flex_hip_end - psi_sol(1);
        psi_error_dot_1 = (psi_error_1 - psi_error_prev_1)/dt;

        tau_d_flex_1_vel = hip_scaling_factor*(Kd_1*psi_error_dot_1)
        tau_d_flex_2_pos = hip_scaling_factor*(Kp_2*psi_error_2)
        tau_d_flex_2_vel = hip_scaling_factor*(Kd_2*psi_error_dot_2)

        tau_d(i,:) = ([k_hip k_knee].*tau_i_G.') + tau_d_G.' + [tau_d_flex_1_vel tau_d_flex_2_pos+tau_d_flex_2_vel];

        if (psi_sol(1) >= flex_peak)
            state = 'flex_assist_2';
            psi_error_1 = 0;
            psi_error_dot_1 = 0;
            psi_error_2 = 0;
            psi_error_dot_2 = 0;
        end
    end
    if (strcmp(state,'flex_assist_2'))
        Kp_1 = 0.1;
        Kd_1 = 0.1;
        Kp_2 = 0.1;
        Kd_2 = 0.01;
        % Real infant torque
        tau_i_1 = tau_i_traj(i,1);
        tau_i_2 = tau_i_traj(i,2);
        tau_i(i,:) = [tau_i_1 tau_i_2];
        % Device torque: Assist hip and knee to reach final state
        psi_error_1 = flex_hip_end - psi_sol(1);
        psi_error_2 = flex_knee_end - psi_sol(3);
        psi_error_dot_1 = (psi_error_1 - psi_error_prev_1)/dt;
        psi_error_dot_2 = (psi_error_2 - psi_error_prev_2)/dt;
        tau_d_flex_1_pos = (Kp_1*psi_error_1)
        tau_d_flex_1_vel = (Kd_1*psi_error_dot_1)
        tau_d_flex_2_pos = (Kp_2*psi_error_2)
        tau_d_flex_2_vel = (Kd_2*psi_error_dot_2)
        tau_d(i,:) = ([k_hip k_knee].*tau_i_G.') + tau_d_G.' + [tau_d_flex_1_pos+tau_d_flex_1_vel tau_d_flex_2_pos+tau_d_flex_2_vel];

        if (psi_sol(2)/abs(psi_sol(2)) == -1)
            finish_time = tic();
            state = 'off';
        end
    end
    psi_error_prev_1 = psi_error_1;
    psi_error_prev_2 = psi_error_2;

    tau_d_hip_max = 1.2;
    if (abs(tau_d(i,1)) > tau_d_hip_max)
        tau_d(i,1) = tau_d(i,1)/abs(tau_d(i,1))*tau_d_hip_max;
    end
    tau_d_knee_max = 0.5;
    if (abs(tau_d(i,1)) > tau_d_knee_max)
        tau_d(i,2) = tau_d(i,2)/abs(tau_d(i,2))*tau_d_knee_max;
    end
%     tau_d(i,:) = ([k_hip k_knee].*tau_i_G.') + tau_d_G.';
%     tau_d(i,:) = [0 0];

%     tau_i_1 = tau_i_traj(i,1);
%     tau_i_2 = tau_i_traj(i,2);
%     tau_i_1 = tau_i_traj(i,1) + (psi(i,1) - psi_sol(1)) + 0.05*(psi_dot(i,1) - psi_sol(2));
%     tau_i_2 = tau_i_traj(i,2) + (psi(i,2) - psi_sol(3)) + 0.05*(psi_dot(i,2) - psi_sol(4));
%     tau_i_1 = 0.3*(psi(i,1) - psi_sol(1)).^2 + 0.05*(psi_dot(i,1) - psi_sol(2));
%     tau_i_2 = 0.3*(psi(i,2) - psi_sol(3)).^2 + 0.05*(psi_dot(i,2) - psi_sol(4));

    tau_sim = [tau_d(i,1) tau_d(i,2) tau_i(i,1) tau_i(i,2)];
    if i < length(t)
        [t_sol,psi_sol] = ode45(@(t,psi) full_system_ode(t,psi,tau_sim,prop),[t(i) t(i+1)],psi_sol);
    else
        [t_sol,psi_sol] = ode45(@(t,psi) full_system_ode(t,psi,tau_sim,prop),[t(i) t(i)+mean(diff(t))],psi_sol);
    end
    psi_sol = psi_sol(end,:);

    state
    hip_scaling_factor

    figure(1)
    subplot(1,2,1)
    cla, hold off
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    axis equal
    xlim([-0.1 0.3])
    ylim([-0.025 0.2])
    cla, hold on
    psi_1 = psi_sol(1);
    psi_2 = psi_sol(3);
    plot_sys(psi_1,psi_2)
    drawnow
    subplot(1,2,2)
    cla, hold on
    plot(t(1:i),tau_i(1:i,1),'-r',t(1:i),tau_i(1:i,2),'-b')

%     % Endstop collision
%     [j1_clamped,theta_1_lim] = clamp(psi_sol(1),0,150*pi/180);
%     [j2_clamped, theta_2_lim] = clamp(psi_sol(3) - theta_1_lim,-150*pi/180,0);
%     psi_sol(1) = theta_1_lim;
%     psi_sol(3) = theta_1_lim+theta_2_lim;
%     c_r = 0.3;
%     if j1_clamped
%         psi_sol(2) = -c_r*psi_sol(2);
%     end
%     if j2_clamped
%         psi_sol(4) = -c_r*psi_sol(4);
%     end
    psi_sol_arr(i,:) = psi_sol;
end

% Calculate inverse dynamics again
psi = [psi_sol_arr(:,1) psi_sol_arr(:,3)];
psi_dot = [psi_sol_arr(:,2) psi_sol_arr(:,2)];
psi_ddot = [0 0 ; diff(psi_dot)./mean(diff(t))];
tau_f_1 = [];
tau_f_2 = [];
for i = 1:length(psi)
    tau = full_system_id(psi(i,:),psi_dot(i,:),psi_ddot(i,:),prop);
    tau_f_1 = [tau_f_1; tau(1,:)];
    tau_f_2 = [tau_f_2; tau(2,:)];
end
tau_f = [tau_f_1(:,1) tau_f_2(:,1)];
tau_i = tau_f - tau_d;

figure(31)
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
cla, hold on
plot(t,tau_d(:,1),'-.r',t,tau_d(:,2),'-.b')
plot(t,tau_i(:,1),'-r',t,tau_i(:,2),'-b')
axis tight, grid on, box on
title('Torque vs Time')
xlabel('Time (s)')
ylabel('Torque (Nm)')
legend('Hip','Knee')

%%
% v = VideoWriter('TypicalKick.gif');
% v.Quality = 100;
% v.FrameRate = 30;
% open(v);

delta_t = 1/30;
interp_index = 1;
for i = 1:length(t)
    if t(i)-t(interp_index(end)) > delta_t
        interp_index = [interp_index i];
    end
end

figure(1)
clf, hold off
set(gca,'XColor','none')
set(gca,'YColor','none')
axis equal
xlim([-0.1 0.3])
ylim([-0.025 0.2])
for i = 1:length(interp_index)-1
    t_start = tic();
    cla, hold on
    psi_1 = psi_sol_arr(interp_index(i),1);
    psi_2 = psi_sol_arr(interp_index(i),3);
    plot_sys(psi_1,psi_2)
    drawnow
    t_plot = toc(t_start);
    pause(delta_t-t_plot)
    %     frame = getframe(gca);
    %     writeVideo(v,frame);
end
% v.Duration
% close(v);

% figure(1)
% clf, hold off
% axis equal
% xlim([-0.2 0.3])
% ylim([-0.2 0.3])
% for c = 1:3
%     for i = 1:length(t)-1
%         t_start = tic();
%         cla
%         plot_sys_rect(psi(1,i),psi(2,i))
%         drawnow
%         delta_t = t(i+1)-t(i);
%         t_plot = toc(t_start);
%         pause(delta_t-t_plot)
%     end
% end


function plot_sys(psi_1,psi_2)
L_1 = 0.108;
L_2 = 0.1;
W = 0.015;
r_O_A = L_1/2*[cos(psi_1); sin(psi_1);];
r_O_B = L_1*[cos(psi_1); sin(psi_1);];
r_O_C = [L_1*cos(psi_1)+L_1/2*cos(psi_2); L_1*sin(psi_1)+L_1/2*sin(psi_2);];
r_O_D = [L_1*cos(psi_1)+L_2*cos(psi_2); L_1*sin(psi_1)+L_2*sin(psi_2);];
[X_1, Y_1] = create_rect(r_O_A,L_1,W,psi_1);
[X_2, Y_2] = create_rect(r_O_C,L_2,W,psi_2);
patch('Vertices',[X_1; Y_1]','Faces',[1 2 3 4],'Facecolor','w','Linewidth',1.2);
patch('Vertices',[X_2; Y_2]','Faces',[1 2 3 4],'Facecolor','w','Linewidth',1.2);
r = 0.01;
circle(0,0,r);
circle(r_O_B(1),r_O_B(2),r);
circle(r_O_D(1),r_O_D(2),r);
end

function [X,Y] = create_rect(center,L,H,psi)
c1 = center(1);
c2 = center(2);
R = [cos(psi), -sin(psi); sin(psi), cos(psi)];
X = [-L/2, L/2, L/2, -L/2];
Y =[-H/2, -H/2, H/2, H/2];
for i=1:4
    T(:,i)=R*[X(i); Y(i)];
end
X = [c1+T(1,1) c1+T(1,2) c1+T(1,3) c1+T(1,4)];
Y = [c2+T(2,1) c2+T(2,2) c2+T(2,3) c2+T(2,4)];
end

function circle(x,y,r)
p = nsidedpoly(1000, 'Center', [x y], 'Radius', r);
plot(p, 'FaceColor', 'w', 'FaceAlpha', 1)
end

function [clamped,new_val] = clamp(val,min_val,max_val)
new_val = min(max(val,min_val),max_val);
clamped = val ~= new_val;
end
