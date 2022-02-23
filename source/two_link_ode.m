function dpsi = two_link_ode(t,psi,tau_1,tau_2,prop)
% Full ODE for two links with different properties
% t     - time variable
% psi   - state variable (psi(1) psi(2) psi(3) psi(4))
% tau_1 - torque input for hip
% tau_2 - torque input for knee
% prop  - vector of constants for the system

L_1 = prop(1);
L_2 = prop(2);
G_1 = prop(3);
G_2 = prop(4);
m_1 = prop(5);
m_2 = prop(6);
I_1 = prop(7);
I_2 = prop(8);

g = 9.81;
dpsi = [psi(2);
    -(I_2*tau_2 - I_2*tau_1 - G_2^2*m_2*tau_1 + 2*G_2^2*m_2*tau_2 - G_2^3*g*m_2^2*cos(psi(3)) + G_2^3*L_1*m_2^2*psi(2)^2*sin(psi(1) - psi(3)) + G_2^3*L_1*m_2^2*psi(4)^2*sin(psi(1) - psi(3)) + G_2^2*L_1*g*m_2^2*cos(psi(1)) + G_1*I_2*g*m_1*cos(psi(1)) + I_2*L_1*g*m_2*cos(psi(1)) + G_2*L_1*m_2*tau_2*cos(psi(1) - psi(3)) + 2*G_2^3*L_1*m_2^2*psi(2)*psi(4)*sin(psi(1) - psi(3)) - G_2^2*L_1*g*m_2^2*cos(psi(1) - psi(3))*cos(psi(3)) + G_1*G_2^2*g*m_1*m_2*cos(psi(1)) + 2*G_2^2*L_1^2*m_2^2*psi(2)^2*cos(psi(1) - psi(3))*sin(psi(1) - psi(3)) - G_2*I_2*L_1*m_2*psi(2)^2*sin(psi(1) - psi(3)) + G_2*I_2*L_1*m_2*psi(4)^2*sin(psi(1) - psi(3)) + 2*G_2*I_2*L_1*m_2*psi(2)*psi(4)*sin(psi(1) - psi(3)))/(I_1*I_2 + G_2^2*L_1^2*m_2^2 + G_1^2*I_2*m_1 + G_2^2*I_1*m_2 + G_2^2*I_2*m_2 + I_2*L_1^2*m_2 + G_1^2*G_2^2*m_1*m_2 - G_2^2*L_1^2*m_2^2*cos(psi(1) - psi(3))^2 + 2*G_2*I_2*L_1*m_2*cos(psi(1) - psi(3)));
    psi(4);
    (I_1*tau_2 + G_1^2*m_1*tau_2 - G_2^2*m_2*tau_1 + 2*G_2^2*m_2*tau_2 + L_1^2*m_2*tau_2 - G_2^3*g*m_2^2*cos(psi(3)) + 2*G_2*L_1^3*m_2^2*psi(2)^2*sin(psi(1) - psi(3)) + G_2^3*L_1*m_2^2*psi(2)^2*sin(psi(1) - psi(3)) + G_2^3*L_1*m_2^2*psi(4)^2*sin(psi(1) - psi(3)) + G_2^2*L_1*g*m_2^2*cos(psi(1)) - G_2*L_1^2*g*m_2^2*cos(psi(3)) - G_2*I_1*g*m_2*cos(psi(3)) - G_2*L_1*m_2*tau_1*cos(psi(1) - psi(3)) + 3*G_2*L_1*m_2*tau_2*cos(psi(1) - psi(3)) + 2*G_2^3*L_1*m_2^2*psi(2)*psi(4)*sin(psi(1) - psi(3)) + G_2*L_1^2*g*m_2^2*cos(psi(1) - psi(3))*cos(psi(1)) - 2*G_2^2*L_1*g*m_2^2*cos(psi(1) - psi(3))*cos(psi(3)) + G_1*G_2^2*g*m_1*m_2*cos(psi(1)) - G_1^2*G_2*g*m_1*m_2*cos(psi(3)) + 3*G_2^2*L_1^2*m_2^2*psi(2)^2*cos(psi(1) - psi(3))*sin(psi(1) - psi(3)) + G_2^2*L_1^2*m_2^2*psi(4)^2*cos(psi(1) - psi(3))*sin(psi(1) - psi(3)) + 2*G_2*I_1*L_1*m_2*psi(2)^2*sin(psi(1) - psi(3)) + 2*G_1^2*G_2*L_1*m_1*m_2*psi(2)^2*sin(psi(1) - psi(3)) + 2*G_2^2*L_1^2*m_2^2*psi(2)*psi(4)*cos(psi(1) - psi(3))*sin(psi(1) - psi(3)) + G_1*G_2*L_1*g*m_1*m_2*cos(psi(1) - psi(3))*cos(psi(1)))/(I_1*I_2 + G_2^2*L_1^2*m_2^2 + G_1^2*I_2*m_1 + G_2^2*I_1*m_2 + G_2^2*I_2*m_2 + I_2*L_1^2*m_2 + G_1^2*G_2^2*m_1*m_2 - G_2^2*L_1^2*m_2^2*cos(psi(1) - psi(3))^2 + 2*G_2*I_2*L_1*m_2*cos(psi(1) - psi(3)));
    ];
end