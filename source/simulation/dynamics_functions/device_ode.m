function dy = device_ode(t,y,tau,prop)
I__d_1 = prop.I__d_1;
I__d_2 = prop.I__d_2;
I__d_m1 = prop.I__d_m1;
I__d_m2 = prop.I__d_m2;
L_1 = prop.L_1;
L__d_g1 = prop.L__d_g1;
L__d_g2 = prop.L__d_g2;
g = prop.g;
m__d_1 = prop.m__d_1;
m__d_2 = prop.m__d_2;
n_1 = prop.n_1;
n_2 = prop.n_2;
dy = ...
    [y(2);
    -(L_1*L__d_g2^3*m__d_2^2*y(4)^2*sin(y(1) - y(3)) - I__d_m2*n_2^2*tau(1) - L__d_g2^2*m__d_2*tau(1) - I__d_2*tau(1) + L_1*L__d_g2^2*g*m__d_2^2*cos(y(1)) + I__d_2*L_1*g*m__d_2*cos(y(1)) + I__d_2*L__d_g1*g*m__d_1*cos(y(1)) + L_1*L__d_g2*m__d_2*tau(2)*cos(y(1) - y(3)) - L_1*L__d_g2^2*g*m__d_2^2*cos(y(1) - y(3))*cos(y(3)) + L_1^2*L__d_g2^2*m__d_2^2*y(2)^2*cos(y(1) - y(3))*sin(y(1) - y(3)) + I__d_m2*L_1*g*m__d_2*n_2^2*cos(y(1)) + I__d_m2*L__d_g1*g*m__d_1*n_2^2*cos(y(1)) + L__d_g1*L__d_g2^2*g*m__d_1*m__d_2*cos(y(1)) + I__d_2*L_1*L__d_g2*m__d_2*y(4)^2*sin(y(1) - y(3)) + I__d_m2*L_1*L__d_g2*m__d_2*n_2^2*y(4)^2*sin(y(1) - y(3)))/(I__d_1*I__d_2 + L_1^2*L__d_g2^2*m__d_2^2 + I__d_2*I__d_m1*n_1^2 + I__d_1*I__d_m2*n_2^2 + I__d_2*L_1^2*m__d_2 + I__d_2*L__d_g1^2*m__d_1 + I__d_1*L__d_g2^2*m__d_2 + I__d_m1*I__d_m2*n_1^2*n_2^2 + I__d_m2*L_1^2*m__d_2*n_2^2 + I__d_m1*L__d_g2^2*m__d_2*n_1^2 + I__d_m2*L__d_g1^2*m__d_1*n_2^2 + L__d_g1^2*L__d_g2^2*m__d_1*m__d_2 - L_1^2*L__d_g2^2*m__d_2^2*cos(y(1) - y(3))^2);
    y(4);
    (I__d_1*tau(2) + I__d_m1*n_1^2*tau(2) + L_1^2*m__d_2*tau(2) + L__d_g1^2*m__d_1*tau(2) + L_1^3*L__d_g2*m__d_2^2*y(2)^2*sin(y(1) - y(3)) - L_1^2*L__d_g2*g*m__d_2^2*cos(y(3)) - I__d_1*L__d_g2*g*m__d_2*cos(y(3)) - L_1*L__d_g2*m__d_2*tau(1)*cos(y(1) - y(3)) + L_1^2*L__d_g2*g*m__d_2^2*cos(y(1) - y(3))*cos(y(1)) + L_1^2*L__d_g2^2*m__d_2^2*y(4)^2*cos(y(1) - y(3))*sin(y(1) - y(3)) - I__d_m1*L__d_g2*g*m__d_2*n_1^2*cos(y(3)) - L__d_g1^2*L__d_g2*g*m__d_1*m__d_2*cos(y(3)) + I__d_1*L_1*L__d_g2*m__d_2*y(2)^2*sin(y(1) - y(3)) + I__d_m1*L_1*L__d_g2*m__d_2*n_1^2*y(2)^2*sin(y(1) - y(3)) + L_1*L__d_g1^2*L__d_g2*m__d_1*m__d_2*y(2)^2*sin(y(1) - y(3)) + L_1*L__d_g1*L__d_g2*g*m__d_1*m__d_2*cos(y(1) - y(3))*cos(y(1)))/(I__d_1*I__d_2 + L_1^2*L__d_g2^2*m__d_2^2 + I__d_2*I__d_m1*n_1^2 + I__d_1*I__d_m2*n_2^2 + I__d_2*L_1^2*m__d_2 + I__d_2*L__d_g1^2*m__d_1 + I__d_1*L__d_g2^2*m__d_2 + I__d_m1*I__d_m2*n_1^2*n_2^2 + I__d_m2*L_1^2*m__d_2*n_2^2 + I__d_m1*L__d_g2^2*m__d_2*n_1^2 + I__d_m2*L__d_g1^2*m__d_1*n_2^2 + L__d_g1^2*L__d_g2^2*m__d_1*m__d_2 - L_1^2*L__d_g2^2*m__d_2^2*cos(y(1) - y(3))^2);
    ];
end