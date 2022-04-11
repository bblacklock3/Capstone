function tau_G = gravity_comp_device(psi,prop)
assert(length(psi)==2)

L_1 = prop.L_1;
L__d_g_1 = prop.L__d_g_1;
L__d_g_2 = prop.L__d_g_2;
g = prop.g;
m__d_1 = prop.m__d_1;
m__d_2 = prop.m__d_2;

tau_G = [g*cos(psi(1))*(L_1*m__d_2 + L__d_g_1*m__d_1); L__d_g_2*g*m__d_2*cos(psi(2))];
end

