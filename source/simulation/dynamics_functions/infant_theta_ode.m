function dy = infant_theta_ode(t,y,tau,prop)
I__i_1 = prop.I__i_1;
I__i_2 = prop.I__i_2;
L_1 = prop.L_1;
L__i_g1 = prop.L__i_g1;
L__i_g2 = prop.L__i_g2;
g = prop.g;
m__i_1 = prop.m__i_1;
m__i_2 = prop.m__i_2;
dy = ...
    [y(2);
    (I__i_2*tau(1) - I__i_2*tau(2) + L__i_g2^2*m__i_2*tau(1) - L__i_g2^2*m__i_2*tau(2) - L_1*L__i_g2^2*g*m__i_2^2*cos(y(1)) - I__i_2*L_1*g*m__i_2*cos(y(1)) - I__i_2*L__i_g1*g*m__i_1*cos(y(1)) - L_1*L__i_g2*m__i_2*tau(2)*cos(y(3)) + L_1*L__i_g2^3*m__i_2^2*y(2)^2*sin(y(3)) + L_1*L__i_g2^3*m__i_2^2*y(4)^2*sin(y(3)) + L_1*L__i_g2^2*g*m__i_2^2*cos(y(1) + y(3))*cos(y(3)) + L_1^2*L__i_g2^2*m__i_2^2*y(2)^2*cos(y(3))*sin(y(3)) + I__i_2*L_1*L__i_g2*m__i_2*y(2)^2*sin(y(3)) + I__i_2*L_1*L__i_g2*m__i_2*y(4)^2*sin(y(3)) - L__i_g1*L__i_g2^2*g*m__i_1*m__i_2*cos(y(1)) + 2*L_1*L__i_g2^3*m__i_2^2*y(2)*y(4)*sin(y(3)) + 2*I__i_2*L_1*L__i_g2*m__i_2*y(2)*y(4)*sin(y(3)))/(I__i_1*I__i_2 + L_1^2*L__i_g2^2*m__i_2^2 + I__i_2*L_1^2*m__i_2 + I__i_2*L__i_g1^2*m__i_1 + I__i_1*L__i_g2^2*m__i_2 + L__i_g1^2*L__i_g2^2*m__i_1*m__i_2 - L_1^2*L__i_g2^2*m__i_2^2*cos(y(3))^2);
    y(4);
    -(I__i_2*tau(1) - I__i_1*tau(2) - I__i_2*tau(2) - L_1^2*m__i_2*tau(2) - L__i_g1^2*m__i_1*tau(2) + L__i_g2^2*m__i_2*tau(1) - L__i_g2^2*m__i_2*tau(2) + L_1^2*L__i_g2*g*m__i_2^2*cos(y(1) + y(3)) - L_1*L__i_g2^2*g*m__i_2^2*cos(y(1)) + I__i_1*L__i_g2*g*m__i_2*cos(y(1) + y(3)) - I__i_2*L_1*g*m__i_2*cos(y(1)) - I__i_2*L__i_g1*g*m__i_1*cos(y(1)) + L_1*L__i_g2*m__i_2*tau(1)*cos(y(3)) - 2*L_1*L__i_g2*m__i_2*tau(2)*cos(y(3)) + L_1*L__i_g2^3*m__i_2^2*y(2)^2*sin(y(3)) + L_1^3*L__i_g2*m__i_2^2*y(2)^2*sin(y(3)) + L_1*L__i_g2^3*m__i_2^2*y(4)^2*sin(y(3)) + L_1*L__i_g2^2*g*m__i_2^2*cos(y(1) + y(3))*cos(y(3)) + 2*L_1^2*L__i_g2^2*m__i_2^2*y(2)^2*cos(y(3))*sin(y(3)) + L_1^2*L__i_g2^2*m__i_2^2*y(4)^2*cos(y(3))*sin(y(3)) - L_1^2*L__i_g2*g*m__i_2^2*cos(y(1))*cos(y(3)) + L__i_g1^2*L__i_g2*g*m__i_1*m__i_2*cos(y(1) + y(3)) + I__i_1*L_1*L__i_g2*m__i_2*y(2)^2*sin(y(3)) + I__i_2*L_1*L__i_g2*m__i_2*y(2)^2*sin(y(3)) + I__i_2*L_1*L__i_g2*m__i_2*y(4)^2*sin(y(3)) - L__i_g1*L__i_g2^2*g*m__i_1*m__i_2*cos(y(1)) + 2*L_1*L__i_g2^3*m__i_2^2*y(2)*y(4)*sin(y(3)) + 2*L_1^2*L__i_g2^2*m__i_2^2*y(2)*y(4)*cos(y(3))*sin(y(3)) + L_1*L__i_g1^2*L__i_g2*m__i_1*m__i_2*y(2)^2*sin(y(3)) + 2*I__i_2*L_1*L__i_g2*m__i_2*y(2)*y(4)*sin(y(3)) - L_1*L__i_g1*L__i_g2*g*m__i_1*m__i_2*cos(y(1))*cos(y(3)))/(I__i_1*I__i_2 + L_1^2*L__i_g2^2*m__i_2^2 + I__i_2*L_1^2*m__i_2 + I__i_2*L__i_g1^2*m__i_1 + I__i_1*L__i_g2^2*m__i_2 + L__i_g1^2*L__i_g2^2*m__i_1*m__i_2 - L_1^2*L__i_g2^2*m__i_2^2*cos(y(3))^2);
    ];
end