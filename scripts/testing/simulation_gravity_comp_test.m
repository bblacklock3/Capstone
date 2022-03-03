%% Constants
prop = load_prop('full_system');

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
    plot_two_link_decoupled(psi_sol(1),psi_sol(3),prop)
    drawnow
    [tau_1, tau_2] = gravity_comp_coupled(psi_sol(1),psi_sol(3), prop);
    [t_sol,psi_sol] = ode45(@(t,psi) full_system_ode(t,psi,gc_factor*tau,prop),[t_sol t_sol+delta_t],psi_sol,opts_1);
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