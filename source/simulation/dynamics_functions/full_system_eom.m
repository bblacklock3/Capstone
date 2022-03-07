function full_eqn = full_system_eom
syms g L_1 L_2

syms L__d_g1 L__d_g2 m__d_1 m__d_2 I__d_1 I__d_2 I__d_m1 n_1 I__d_m2 n_2
syms tau__d_1 tau__d_2

syms L__i_g1 L__i_g2 m__i_1 m__i_2 I__i_1 I__i_2
syms tau__i_1 tau__i_2

syms psi_1 psi_dot_1 psi_ddot_1 psi_2 psi_dot_2 psi_ddot_2

full_eqn = [- I__i_2*(psi_ddot_1 - psi_ddot_2) + psi_ddot_1*(m__i_2*L_1^2 + 2*m__i_2*cos(psi_1 - psi_2)*L_1*L__i_g2 + m__i_1*L__i_g1^2 + m__i_2*L__i_g2^2 + I__i_1 + I__i_2) + (m__d_2*L_1^2 + m__d_1*L__d_g1^2 + I__d_m1*n_1^2 + I__d_1)*psi_ddot_1 - L__i_g2^2*m__i_2*(psi_ddot_1 - psi_ddot_2) + L_1*L__d_g2*m__d_2*sin(psi_1 - psi_2)*psi_dot_2^2 + L_1*g*m__d_2*cos(psi_1) + L__d_g1*g*m__d_1*cos(psi_1) + L_1*L__d_g2*m__d_2*psi_ddot_2*cos(psi_1 - psi_2) + L_1*g*m__i_2*cos(psi_1) + L__i_g1*g*m__i_1*cos(psi_1) + L__i_g2*g*m__i_2*cos(psi_2) - L_1*L__i_g2*m__i_2*cos(psi_1 - psi_2)*(psi_ddot_1 - psi_ddot_2) + L_1*L__i_g2*m__i_2*sin(psi_1 - psi_2)*(psi_dot_1 - psi_dot_2)^2 - 2*L_1*L__i_g2*m__i_2*psi_dot_1*sin(psi_1 - psi_2)*(psi_dot_1 - psi_dot_2);
    I__i_2*psi_ddot_1 - (m__i_2*L__i_g2^2 + I__i_2)*(psi_ddot_1 - psi_ddot_2) + (m__d_2*L__d_g2^2 + I__d_m2*n_2^2 + I__d_2)*psi_ddot_2 - L_1*L__d_g2*m__d_2*sin(psi_1 - psi_2)*psi_dot_1^2 + L__d_g2*g*m__d_2*cos(psi_2) + L_1*L__d_g2*m__d_2*psi_ddot_1*cos(psi_1 - psi_2) + L__i_g2^2*m__i_2*psi_ddot_1 + L__i_g2*g*m__i_2*cos(psi_2) - L_1*L__i_g2*m__i_2*psi_dot_1^2*sin(psi_1 - psi_2) + L_1*L__i_g2*m__i_2*psi_ddot_1*cos(psi_1 - psi_2)];
end

