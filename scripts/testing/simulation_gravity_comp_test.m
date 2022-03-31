%% Constants
% anthropometric data, from John's .xls file
age = 0.147945205479452;            % age in years
age_weeks = age * 365/7;
bodyMass = 3.4;                     % weight in kg

length_thigh = 0.108;               %
length_shank = 0.1;                 % segment lengths in m
length_foot  = 0.091;               %

circ_thigh = 0.187;                 %
circ_shank = 0.135;                 % segment circumference in m
circ_foot = 0.097;                  %

width_foot = 0.029;     

mass_thigh_XLS = 0.247359160273973;     %
mass_shank_XLS = 0.109713428;           % segment mass in kg, from the .xls file
mass_foot_XLS  = 0.079837736;           %

mass_thigh_SZ = 6.9126e-02*age      + 2.9582e-00*length_thigh + 3.1541e-00*circ_thigh - 6.7217e-01; %
mass_shank_SZ = 6.5138e-03*bodyMass + 1.8158e-00*length_shank + 1.8743e-00*circ_shank - 3.5460e-01; % estimates from Schneider and Zernicke (1992)
mass_foot_SZ  = 2.9331e-03*bodyMass + 1.2405e-00*length_foot  + 1.9337e-00*width_foot - 1.0250e-01; %

mass_thigh_SJ = 2.39e-01 + 1.24e-02*age_weeks; %
mass_shank_SJ = 1.69e-01 + 4.16e-03*age_weeks; % estimates from Sun and Jensen (1994)
mass_foot_SJ  = 5.22e-02 + 2.12e-03*age_weeks; %

mass_thigh = mass_thigh_SZ;
mass_shank = mass_shank_SZ;
mass_foot  = mass_foot_SZ;

com_thigh = 0.4859*length_thigh;    % distance of center of mass from proximal joint
com_shank = 0.4377*length_shank;    % estimation according to Schneider and Zernicke (1992)
com_foot = 0.3469*length_foot;      %

% estimates for moments of inertia
% according to Schneider and Zernicke (1992)
tvInertia_thigh_1 = 0.017943*length_thigh + 0.005699*circ_thigh - 0.0027078;
tvInertia_shank_1 = 0.000018660*bodyMass + 0.0085431*length_shank + 0.0016127*circ_shank - 0.0011192;

% estimates for moments of inertia
% according to Sun and Jensen (1994)
% units are kg m^2
tvInertia_thigh_2 = -1.90e-04 + 4.35e-05*age_weeks;
tvInertia_shank_2 =  8.74e-05 + 9.89e-06*age_weeks;

% Real values
L_1 = length_thigh;
L_2 = length_shank;
G_1 = com_thigh;
G_2 = (com_shank*mass_shank+com_foot*mass_foot)/(mass_shank+mass_foot);
m_1 = 2*mass_thigh;
m_2 = 2*(mass_shank+mass_foot);
I_1 = tvInertia_thigh_1;
I_2 = tvInertia_shank_1;
prop = [L_1 L_2 G_1 G_2 m_1 m_2 I_1 I_2];

% Simple values
% L_1 = 0.1;
% L_2 = 0.1;
% G_1 = 0.05;
% G_2 = 0.05; 
% m_1 = 1/9.81;
% m_2 = 1/9.81;
% I_1 = 0.001;
% I_2 = 0.001;
% prop = [L_1 L_2 G_1 G_2 m_1 m_2 I_1 I_2];

%% 
figure(1)
cla, hold off
axis equal
xlim([-0.1 0.4])
ylim([-0.25 0.25])
opts_1 = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',0.025);

t_sol_arr = [];
psi_sol_arr = [];
tau_arr = [];

gc_factor = 0.8;
t_sol = 0;
psi_sol = [pi/4 0 -pi/4 0];
delta_t = 0.01;
while t_sol < 1
    t1_start = tic();
    cla, hold on
    plot_two_link_coupled(psi_sol(1),psi_sol(3),prop)
    drawnow
    [tau_1, tau_2] = gravity_comp_coupled(psi_sol(1),psi_sol(3), prop);
    [t_sol,psi_sol] = ode45(@(t,psi) two_link_coupled_ode(t,psi,gc_factor*tau_1,gc_factor*tau_2,prop),[t_sol t_sol+delta_t],psi_sol,opts_1);
    t_sol = t_sol(end);
    psi_sol = psi_sol(end,:);
    t1 = toc(t1_start);
    t_sol_arr = [t_sol_arr; t_sol];
    psi_sol_arr = [psi_sol_arr; psi_sol;];
    tau_arr = [tau_arr; tau_1 tau_2;];
    pause(delta_t-t1)
end

psi_sol_arr = psi_sol_arr*180/pi;
figure(2)
subplot(1,3,1)
plot(t_sol_arr,psi_sol_arr(:,1),'-r',t_sol_arr,psi_sol_arr(:,3),'-b')
axis tight, grid on, box on
title('Position vs Time')
xlabel('Time (s)')
ylabel('Angular Position (deg)')
legend('Hip','Knee')
subplot(1,3,2)
plot(t_sol_arr,psi_sol_arr(:,2),'--r',t_sol_arr,psi_sol_arr(:,4),'--b')
axis tight, grid on, box on
title('Velocity vs Time')
xlabel('Time (s)')
ylabel('Angular Veloctiy (deg/s)')
legend('Hip','Knee')
subplot(1,3,3)
plot(t_sol_arr,tau_arr(:,1),'-.r',t_sol_arr,tau_arr(:,2),'-.b')
axis tight, grid on, box on
title('Torque vs Time')
xlabel('Time (s)')
ylabel('Torque (Nm)')
legend('Hip','Knee')

function plot_two_link_decoupled(psi_1,psi_2,prop)
L_1 = prop(1);
L_2 = prop(2);
G_1 = prop(3);
G_2 = prop(4);
r_O = [0; 0];
r_O_A = G_1*[cos(psi_1); sin(psi_1);];
r_O_B = L_1*[cos(psi_1); sin(psi_1);];
r_O_C = [L_1*cos(psi_1)+G_2*cos(psi_2); L_1*sin(psi_1)+G_2*sin(psi_2);];
r_O_D = [L_1*cos(psi_1)+L_2*cos(psi_2); L_1*sin(psi_1)+L_2*sin(psi_2);];
link_1 = [r_O r_O_B];
plot(link_1(1,:),link_1(2,:))
link_2 = [r_O_B r_O_D];
plot(link_2(1,:),link_2(2,:))

plot(r_O_A(1,:),r_O_A(2,:),'b*')
plot(r_O_C(1,:),r_O_C(2,:),'r*')
end

function plot_two_link_coupled(theta_1,theta_2,prop)
L_1 = prop(1);
L_2 = prop(2);
G_1 = prop(3);
G_2 = prop(4);
r_O = [0; 0];
r_O_A = G_1*[cos(theta_1); sin(theta_1);];
r_O_B = L_1*[cos(theta_1); sin(theta_1);];
r_O_C = [L_1*cos(theta_1)+G_2*cos(theta_1+theta_2); L_1*sin(theta_1)+G_2*sin(theta_1+theta_2);];
r_O_D = [L_1*cos(theta_1)+L_2*cos(theta_1+theta_2); L_1*sin(theta_1)+L_2*sin(theta_1+theta_2);];
link_1 = [r_O r_O_B];
plot(link_1(1,:),link_1(2,:))
link_2 = [r_O_B r_O_D];
plot(link_2(1,:),link_2(2,:))
plot(r_O_A(1,:),r_O_A(2,:),'b*')
plot(r_O_C(1,:),r_O_C(2,:),'r*')
end