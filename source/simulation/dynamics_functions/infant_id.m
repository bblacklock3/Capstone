function tau = infant_id(psi,psi_dot,psi_ddot,prop)
assert(length(psi)==2)
assert(length(psi_dot)==2)
assert(length(psi_ddot)==2)

I__i_1 = prop.I__i_1;
I__i_2 = prop.I__i_2;
L_1 = prop.L_1;
L__i_g_1 = prop.L__i_g_1;
L__i_g_2 = prop.L__i_g_2;
g = prop.g;
m__i_1 = prop.m__i_1;
m__i_2 = prop.m__i_2;

M = [m__i_2*L_1^2 + L__i_g_2*m__i_2*cos(psi(1) - psi(2))*L_1 + m__i_1*L__i_g_1^2 + I__i_1, m__i_2*L__i_g_2^2 + L_1*m__i_2*cos(psi(1) - psi(2))*L__i_g_2 + I__i_2; L_1*L__i_g_2*m__i_2*cos(psi(1) - psi(2)), m__i_2*L__i_g_2^2 + I__i_2];
C = [-L_1*L__i_g_2*m__i_2*sin(psi(1) - psi(2)), L_1*L__i_g_2*m__i_2*sin(psi(1) - psi(2)); -L_1*L__i_g_2*m__i_2*sin(psi(1) - psi(2)), 0];
G = [g*(L_1*m__i_2*cos(psi(1)) + L__i_g_1*m__i_1*cos(psi(1)) + L__i_g_2*m__i_2*cos(psi(2))); L__i_g_2*g*m__i_2*cos(psi(2))];

if size(psi_dot,1) ~= 2
    psi_dot = psi_dot.';
end

if size(psi_ddot,1) ~= 2
    psi_ddot = psi_ddot.';
end

tau_M = M*psi_ddot;
tau_M(tau_M > prop.max_hip_torque) = prop.max_hip_torque;
tau_M(tau_M < -prop.max_hip_torque) = -prop.max_hip_torque;

tau_C = C*psi_dot.^2;
tau_G = G;
tau_total = tau_M + tau_C + tau_G;
tau = [tau_total tau_M tau_C tau_G];
end

