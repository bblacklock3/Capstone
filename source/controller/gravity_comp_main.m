function tau_G = gravity_comp_main(psi,prop)
assert(length(psi)==2)
L_1 = prop.L_1;
L__d_g_1 = prop.L__d_g_1;
L__d_g_2 = prop.L__d_g_2;
L__i_g_1 = prop.L__i_g_1;
L__i_g_2 = prop.L__i_g_2;
g = prop.g;
m__d_1 = prop.m__d_1;
m__d_2 = prop.m__d_2;
m__i_1 = prop.m__i_1;
m__i_2 = prop.m__i_2;

tau_G = [g*(L_1*m__d_2*cos(psi(1)) + L_1*m__i_2*cos(psi(1)) + L__d_g_1*m__d_1*cos(psi(1)) + L__i_g_1*m__i_1*cos(psi(1)) + L__i_g_2*m__i_2*cos(psi(2))); g*cos(psi(2))*(L__d_g_2*m__d_2 + L__i_g_2*m__i_2)];
end

