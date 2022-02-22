%% Constants
L_1 = 0.15;
L_2 = 0.1;
G_1 = 0.08;
G_2 = 0.06; 
m_1 = 0.1;
m_2 = 0.1;
I_1 = 0.01;
I_2 = 0.01;
prop = [L_1 L_2 G_1 G_2 m_1 m_2 I_1 I_2];

%%
tspan = [0 5];
psi_0 = [0 0 0 0];
% psi_0 = [-pi/2 0 -pi/2 0];
tau_1 = 0;
tau_2 = 0;
opts_1 = odeset('RelTol',1e-3,'AbsTol',1e-4,'MaxStep',0.025);
[t_sol,psi_sol] = ode45(@(t,psi) two_link_ode(t,psi,tau_1,tau_2,prop),tspan,psi_0,opts_1);
%% 
delta_t = 0.1;
interp_index = 1;
for i = 1:length(t_sol)
    if t_sol(i)-t_sol(interp_index(end)) > delta_t
        interp_index = [interp_index i];
    end
end

figure(1)
cla, hold off
xlim([-0.2 0.4])
ylim([-0.25 0.25])
for i = 1:length(interp_index)
    cla, hold on
    t1_start = tic();
    plot_two_link(psi_sol(interp_index(i),1),psi_sol(interp_index(i),3))
    t2_start = tic();
    drawnow
    pause(delta_t)
    t1 = toc(t2_start);
end

function plot_two_link(psi_1,psi_2)
L_1 = 0.15; L_2 = 0.1; G_1 = 0.08; G_2 = 0.06;
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