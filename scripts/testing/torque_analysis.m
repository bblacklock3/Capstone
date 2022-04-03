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

Knee = -Knee;

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
knee_deg = hip_deg+knee_deg;
warning('on');

trim_len = 5;
t(:,1:trim_len) = [];
t(:,end-trim_len:end) = [];
hip_deg(:,1:trim_len) = [];
hip_deg(:,end-trim_len:end) = [];
knee_deg(:,1:trim_len) = [];
knee_deg(:,end-trim_len:end) = [];

%% Plotting constants
black = [0 0 0];
blue = [57 106 177]./255;
red = [204 37 41]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
font_size = 10;
text_spacing = 0.2;
kinematics_box = [4 1];
dynamics_box = [7 3];

%% Plot kinematics
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

hip_rad = hip_deg*pi/180;
knee_rad = knee_deg*pi/180;

psi = [hip_rad(1,:); knee_rad(1,:)];
psi_dot = [hip_rad(2,:); knee_rad(2,:)];
psi_ddot = [hip_rad(3,:); knee_rad(3,:)];

%% Min infant values
min_age_weeks = 8;
[min_tau_i_1, min_tau_i_2] = calc_infant_torque(psi,psi_dot,psi_ddot,min_age_weeks);
min_avg_i_1 = mean(abs(min_tau_i_1)).';
min_avg_i_2 = mean(abs(min_tau_i_2)).';

% Plot min hip values
figure(20)
clf, hold on
plot(t,min_tau_i_1(:,1))
plot(t,min_tau_i_1(:,2))
plot(t,min_tau_i_1(:,3))
plot(t,min_tau_i_1(:,4))
colors = {black; red; green; blue;};
plt_min_hip_infant = Plot();
plt_min_hip_infant.Title = [num2str(min_age_weeks) ' Week Old Infant: Hip Torque for Typical Kick'];
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
text(x_start,y_start,0,sprintf('Avg Net = %0.3f (Nm)',min_avg_i_1(1)),'FontSize',font_size,'Color',colors{1})
text(x_start,y_start-y_space,0,sprintf('Avg Mass = %0.3f (Nm)',min_avg_i_1(2)),'FontSize',font_size,'Color',colors{2})
text(x_start,y_start-2*y_space,0,sprintf('Avg Cent. = %0.3f (Nm)',min_avg_i_1(3)),'FontSize',font_size,'Color',colors{3})
text(x_start,y_start-3*y_space,0,sprintf('Avg Grav. = %0.3f (Nm)',min_avg_i_1(4)),'FontSize',font_size,'Color',colors{4})
plt_min_hip_infant.Resolution = 300;

% Plot min knee values
figure(21)
clf, hold on
plot(t,min_tau_i_2(:,1))
plot(t,min_tau_i_2(:,2))
plot(t,min_tau_i_2(:,3))
plot(t,min_tau_i_2(:,4))
colors = {black; red; green; blue;};
plt_min_knee_infant = Plot();
plt_min_knee_infant.Title = [num2str(min_age_weeks) ' Week Old Infant: Knee Torque for Typical Kick'];
plt_min_knee_infant.YLim = [-0.15 0.15];
plt_min_knee_infant.TickLength = [0.015 0.015];
plt_min_knee_infant.XLabel = 'Time (s)';
plt_min_knee_infant.YLabel = 'Torque (Nm)';
plt_min_knee_infant.BoxDim = dynamics_box;
plt_min_knee_infant.FontSize = font_size;
plt_min_knee_infant.Colors = colors;
plt_min_knee_infant.LineWidth = [3, 2, 2, 2];
plt_min_knee_infant.LineStyle = {'-','-','-','-'};
plt_min_knee_infant.Legend = {'Net','Mass','Cent.','Grav.'};
plt_min_knee_infant.LegendBox = 'on';
x_start = 0.4; y_start = -0.05; y_space = diff(plt_min_knee_infant.YLim)*text_spacing/plt_min_knee_infant.BoxDim(2);
text(x_start,y_start,0,sprintf('Avg Net = %0.3f (Nm)',min_avg_i_2(1)),'FontSize',font_size,'Color',colors{1})
text(x_start,y_start-y_space,0,sprintf('Avg Mass = %0.3f (Nm)',min_avg_i_2(2)),'FontSize',font_size,'Color',colors{2})
text(x_start,y_start-2*y_space,0,sprintf('Avg Cent. = %0.3f (Nm)',min_avg_i_2(3)),'FontSize',font_size,'Color',colors{3})
text(x_start,y_start-3*y_space,0,sprintf('Avg Grav. = %0.3f (Nm)',min_avg_i_2(4)),'FontSize',font_size,'Color',colors{4})
plt_min_knee_infant.Resolution = 300;

%% Max infant values
max_age_weeks = 20;
[max_tau_i_1, max_tau_i_2] = calc_infant_torque(psi,psi_dot,psi_ddot,max_age_weeks);
max_avg_i_1 = mean(abs(max_tau_i_1)).';
max_avg_i_2 = mean(abs(max_tau_i_2)).';

% Plot min hip values
figure(30)
clf, hold on
plot(t,max_tau_i_1(:,1))
plot(t,max_tau_i_1(:,2))
plot(t,max_tau_i_1(:,3))
plot(t,max_tau_i_1(:,4))
colors = {black; red; green; blue;};
plt_max_hip_infant = Plot();
plt_max_hip_infant.Title = [num2str(max_age_weeks) ' Week Old Infant: Hip Torque for Typical Kick'];
plt_max_hip_infant.TickLength = [0.015 0.015];
plt_max_hip_infant.XLabel = 'Time (s)';
plt_max_hip_infant.YLabel = 'Torque (Nm)';
plt_max_hip_infant.BoxDim = dynamics_box;
plt_max_hip_infant.FontSize = font_size;
plt_max_hip_infant.Colors = colors;
plt_max_hip_infant.LineWidth = [3, 2, 2, 2];
plt_max_hip_infant.LineStyle = {'-','-','-','-'};
plt_max_hip_infant.Legend = {'Net','Mass','Cent.','Grav.'};
plt_max_hip_infant.LegendBox = 'on';
x_start = 0.4; y_start = -0.2; y_space = diff(plt_max_hip_infant.YLim)*text_spacing/plt_max_hip_infant.BoxDim(2);
text(x_start,y_start,0,sprintf('Avg Net = %0.3f (Nm)',max_avg_i_1(1)),'FontSize',font_size,'Color',colors{1})
text(x_start,y_start-y_space,0,sprintf('Avg Mass = %0.3f (Nm)',max_avg_i_1(2)),'FontSize',font_size,'Color',colors{2})
text(x_start,y_start-2*y_space,0,sprintf('Avg Cent. = %0.3f (Nm)',max_avg_i_1(3)),'FontSize',font_size,'Color',colors{3})
text(x_start,y_start-3*y_space,0,sprintf('Avg Grav. = %0.3f (Nm)',max_avg_i_1(4)),'FontSize',font_size,'Color',colors{4})
plt_max_hip_infant.Resolution = 300;

% Plot min knee values
figure(31)
clf, hold on
plot(t,max_tau_i_2(:,1))
plot(t,max_tau_i_2(:,2))
plot(t,max_tau_i_2(:,3))
plot(t,max_tau_i_2(:,4))
colors = {black; red; green; blue;};
plt_max_knee_infant = Plot();
plt_max_knee_infant.Title = [num2str(max_age_weeks) ' Week Old Infant: Knee Torque for Typical Kick'];
plt_max_knee_infant.YLim = [-0.18 0.2];
plt_max_knee_infant.TickLength = [0.015 0.015];
plt_max_knee_infant.XLabel = 'Time (s)';
plt_max_knee_infant.YLabel = 'Torque (Nm)';
plt_max_knee_infant.BoxDim = dynamics_box;
plt_max_knee_infant.FontSize = font_size;
plt_max_knee_infant.Colors = colors;
plt_max_knee_infant.LineWidth = [3, 2, 2, 2];
plt_max_knee_infant.LineStyle = {'-','-','-','-'};
plt_max_knee_infant.Legend = {'Net','Mass','Cent.','Grav.'};
plt_max_knee_infant.LegendBox = 'on';
x_start = 0.4; y_start = -0.05; y_space = diff(plt_max_knee_infant.YLim)*text_spacing/plt_max_knee_infant.BoxDim(2);
text(x_start,y_start,0,sprintf('Avg Net = %0.3f (Nm)',max_avg_i_2(1)),'FontSize',font_size,'Color',colors{1})
text(x_start,y_start-y_space,0,sprintf('Avg Mass = %0.3f (Nm)',max_avg_i_2(2)),'FontSize',font_size,'Color',colors{2})
text(x_start,y_start-2*y_space,0,sprintf('Avg Cent. = %0.3f (Nm)',max_avg_i_2(3)),'FontSize',font_size,'Color',colors{3})
text(x_start,y_start-3*y_space,0,sprintf('Avg Grav. = %0.3f (Nm)',max_avg_i_2(4)),'FontSize',font_size,'Color',colors{4})
plt_max_knee_infant.Resolution = 300;

%% Combined values
[tau_d_1, tau_d_2] = calc_device_torque(psi,psi_dot,psi_ddot);
avg_d_1 = mean(abs(tau_d_1)).';
avg_d_2 = mean(abs(tau_d_2)).';

min_hip_torque_ratio = min_avg_i_1(1)/avg_d_1(1);
min_knee_torque_ratio = min_avg_i_2(1)/avg_d_2(1);
max_hip_torque_ratio = max_avg_i_1(1)/avg_d_1(1);
max_knee_torque_ratio = max_avg_i_2(1)/avg_d_2(1);

text_spacing_1 = 0.38;
text_spacing_2 = 0.25;

% Plot combined hip values
figure(40)
clf, hold on
plot(t,tau_d_1(:,1))
plot(t,min_tau_i_1(:,1))
plot(t,max_tau_i_1(:,1))
colors = {red; green; blue;};
plt_hip_combined = Plot();
plt_hip_combined.Title = 'Hip Torque for Typical Kick';
plt_hip_combined.TickLength = [0.015 0.015];
plt_hip_combined.XLabel = 'Time (s)';
plt_hip_combined.YLabel = 'Torque (Nm)';
plt_hip_combined.BoxDim = dynamics_box;
plt_hip_combined.FontSize = font_size;
plt_hip_combined.Colors = colors;
plt_hip_combined.LineWidth = [3, 2, 2, 2];
plt_hip_combined.LineStyle = {'-','-','-','-'};
plt_hip_combined.Legend = {'Device','Infant Min','Infant Max'};
plt_hip_combined.LegendBox = 'on';
x_start = 0.58; y_start = -0.15;
y_space_1 = diff(plt_hip_combined.YLim)*text_spacing_1/plt_hip_combined.BoxDim(2);
y_space_2 = diff(plt_hip_combined.YLim)*text_spacing_2/plt_hip_combined.BoxDim(2);
text(x_start,y_start,0,'Torque Ratio - $$R = \frac{\tau_{i}}{\tau_{d}}$$','FontSize',14,'Color',black,'Interpreter','latex','HorizontalAlignment','center')
text(x_start,y_start-y_space_1,0,sprintf('Min Avg Ratio $$(R_{min})$$ = %0.1f',min_hip_torque_ratio),'FontSize',14,'Color',black,'Interpreter','latex','HorizontalAlignment','center')
text(x_start,y_start-y_space_1-y_space_2,0,sprintf('Max Avg Ratio $$(R_{max})$$ = %0.1f',max_hip_torque_ratio),'FontSize',14,'Color',black,'Interpreter','latex','HorizontalAlignment','center')
plt_hip_combined.Resolution = 300;

% Plot combined knee values
figure(41)
clf, hold on
plot(t,tau_d_2(:,1))
plot(t,min_tau_i_2(:,1))
plot(t,max_tau_i_2(:,1))
colors = {red; green; blue;};
plt_knee_combined = Plot();
plt_knee_combined.Title = 'Knee Torque for Typical Kick';
plt_knee_combined.YLim = [-0.18 0.2];
plt_knee_combined.TickLength = [0.015 0.015];
plt_knee_combined.XLabel = 'Time (s)';
plt_knee_combined.YLabel = 'Torque (Nm)';
plt_knee_combined.BoxDim = dynamics_box;
plt_knee_combined.FontSize = font_size;
plt_knee_combined.Colors = colors;
plt_knee_combined.LineWidth = [3, 2, 2, 2];
plt_knee_combined.LineStyle = {'-','-','-','-'};
plt_knee_combined.Legend = {'Device','Infant Min','Infant Max'};
plt_knee_combined.LegendBox = 'on';
x_start = 0.58; y_start = -0.04;
y_space_1 = diff(plt_knee_combined.YLim)*text_spacing_1/plt_knee_combined.BoxDim(2);
y_space_2 = diff(plt_knee_combined.YLim)*text_spacing_2/plt_knee_combined.BoxDim(2);
text(x_start,y_start,0,'Torque Ratio - $$R = \frac{\tau_{i}}{\tau_{d}}$$','FontSize',14,'Color',black,'Interpreter','latex','HorizontalAlignment','center')
text(x_start,y_start-y_space_1,0,sprintf('Min Avg Ratio $$(R_{min})$$ = %0.1f',min_knee_torque_ratio),'FontSize',14,'Color',black,'Interpreter','latex','HorizontalAlignment','center')
text(x_start,y_start-y_space_1-y_space_2,0,sprintf('Max Avg Ratio $$(R_{max})$$ = %0.1f',max_knee_torque_ratio),'FontSize',14,'Color',black,'Interpreter','latex','HorizontalAlignment','center')
plt_knee_combined.Resolution = 300;

%%
colors = {red; green; blue; black; black; black};

% Plot combined hip values
figure(50)
clf, hold on
plot(t,tau_d_1(:,1))
plot(t,min_tau_i_1(:,1))
plot(t,max_tau_i_1(:,1))
yline(0.0484,'--k')
yline(0.0185,'-.k')
yline(0.007,'-ok')
plt_hip_combined_backdrive = Plot();
plt_hip_combined_backdrive.Title = 'Hip Torque for Typical Kick';
plt_hip_combined_backdrive.XLim = [min(t) max(t)];
plt_hip_combined_backdrive.TickLength = [0.015 0.015];
plt_hip_combined_backdrive.XLabel = 'Time (s)';
plt_hip_combined_backdrive.YLabel = 'Torque (Nm)';
plt_hip_combined_backdrive.BoxDim = dynamics_box;
plt_hip_combined_backdrive.FontSize = font_size;
plt_hip_combined_backdrive.Colors = colors;
plt_hip_combined_backdrive.LineWidth = [2, 2, 2, 2, 1, 1, 1];
plt_hip_combined_backdrive.LineStyle = {'-','-','-','-','--','-.','-o'};
plt_hip_combined_backdrive.Legend = {'Device','Infant Min','Infant Max','Backdrive: I=0','Backdrive: Calib','Backdrive: Calib+Bias'};
plt_hip_combined_backdrive.LegendBox = 'on';
plt_hip_combined_backdrive.Resolution = 300;

% Plot combined knee values
figure(51)
clf, hold on
plot(t,tau_d_2(:,1))
plot(t,min_tau_i_2(:,1))
plot(t,max_tau_i_2(:,1))
yline(0.0484)
yline(0.0185)
yline(0.007)
plt_knee_combined_backdrive = Plot();
plt_knee_combined_backdrive.Title = 'Knee Torque for Typical Kick';
plt_knee_combined_backdrive.XLim = [min(t) max(t)];
plt_knee_combined_backdrive.YLim = [-0.05 0.2];
plt_knee_combined_backdrive.TickLength = [0.015 0.015];
plt_knee_combined_backdrive.XLabel = 'Time (s)';
plt_knee_combined_backdrive.YLabel = 'Torque (Nm)';
plt_knee_combined_backdrive.BoxDim = dynamics_box;
plt_knee_combined_backdrive.FontSize = font_size;
plt_knee_combined_backdrive.Colors = colors;
plt_knee_combined_backdrive.LineWidth = [2, 2, 2, 2, 1, 1, 1];
plt_knee_combined_backdrive.LineStyle = {'-','-','-','-','--','-.','-o'};
plt_knee_combined_backdrive.Legend = {'Device','Infant Min','Infant Max','Backdrive: I=0','Backdrive: Calib','Backdrive: Calib+Bias'};
plt_knee_combined_backdrive.LegendBox = 'on';
plt_knee_combined_backdrive.Resolution = 300;

%% Save plots
% plt_pos.export('TypicalKickPos.png');
% plt_vel.export('TypicalKickVel.png');
% plt_accel.export('TypicalKickAccel.png');
% plt_min_hip_infant.export('MinHipInfant.png'); 
% plt_min_knee_infant.export('MinKneeInfant.png'); 
% plt_max_hip_infant.export('MaxHipInfant.png'); 
% plt_max_knee_infant.export('MaxKneeInfant.png'); 
% plt_hip_combined.export('HipCombined.png'); 
% plt_knee_combined.export('KneeCombined.png');
plt_hip_combined_backdrive.export('HipCombinedBackdrive.png'); 
plt_knee_combined_backdrive.export('KneeCombinedBackdrive.png');

%% Animate plot
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

%% Export animation
% v = VideoWriter('TypicalKick.gif');
% v.Quality = 100;
% v.FrameRate = 30;
% open(v);
% 
% delta_t = 1/30;
% interp_index = 1;
% for i = 1:length(t)
%     if t(i)-t(interp_index(end)) > delta_t
%         interp_index = [interp_index i];
%     end
% end
% 
% figure(1)
% clf, hold off
% set(gca,'XColor','none')
% set(gca,'YColor','none')
% axis equal
% xlim([-0.1 0.2])
% ylim([-0.025 0.2])
% for i = 1:length(interp_index)-1
%     cla, hold on
%     psi_1 = psi(1,interp_index(i));
%     psi_2 = psi(2,interp_index(i));
%     plot_sys(psi_1,psi_2)
%     drawnow
%     frame = getframe(gca);
%     writeVideo(v,frame);
% end
% v.Duration
% close(v);


%% Functions
function [tau_i_1, tau_i_2] = calc_infant_torque(psi,psi_dot,psi_ddot,age_weeks)
prop = load_prop('full');
[prop.m__i_1, prop.m__i_2, prop.I__i_1, prop.I__i_2] = SJ_mass(age_weeks);
tau_i_1 = [];
tau_i_2 = [];
for i = 1:length(psi)
    tau_i = infant_id(psi(:,i),psi_dot(:,i),psi_ddot(:,i),prop);
    tau_i_1 = [tau_i_1; tau_i(1,:)];
    tau_i_2 = [tau_i_2; tau_i(2,:)];
end
end

function [tau_d_1, tau_d_2] = calc_device_torque(psi,psi_dot,psi_ddot)
prop = load_prop('full');
tau_d_1 = [];
tau_d_2 = [];
for i = 1:length(psi)
    tau_d = device_id(psi(:,i),psi_dot(:,i),psi_ddot(:,i),prop);
    tau_d_1 = [tau_d_1; tau_d(1,:)];
    tau_d_2 = [tau_d_2; tau_d(2,:)];
end
end

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
