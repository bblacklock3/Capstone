function [tau_1, tau_2] = gravity_comp_coupled(theta_1, theta_2, prop)
% Gravitational torque calculation for two links with different properties
% psi_1 - absolute angle for hip
% psi_2 - absolute angle for knee
% prop  - vector of constants for the system
% tau_1 - gravitational torque for hip joint
% tau_2 - gravitational torque for knee joint

L_1 = prop(1);
L_2 = prop(2);
G_1 = prop(3);
G_2 = prop(4);
m_1 = prop(5);
m_2 = prop(6);
I_1 = prop(7);
I_2 = prop(8);
g = 9.81;
tau_1 = G_1*cos(theta_1)*m_1*g + (L_1*cos(theta_1)+G_2*cos(theta_1+theta_2))*m_2*g;
tau_2 = G_2*cos(theta_1+theta_2)*m_2*g;
end