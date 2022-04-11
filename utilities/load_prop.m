function prop = load_prop(name)
prop = [];
prop.g = 9.81;
prop.L_1 = 0.1080;
prop.L_2 = 0.1000;
prop.L__i_g_1 = 0;
prop.L__i_g_2 = 0;
prop.m__i_1 = 0;
prop.m__i_2 = 0;
prop.I__i_1 = 0;
prop.I__i_2 = 0;
prop.L__d_g_1 = 0;
prop.L__d_g_2 = 0;
prop.m__d_1 = 0;
prop.m__d_2 = 0;
prop.I__d_1 = 0;
prop.I__d_2 = 0;
prop.I__d_m_1 = 0;
prop.I__d_m_2 = 0;
prop.n_1 = 1;
prop.n_2 = 1;
switch (name)
    case 'ones'
        prop.L_1 = 1;
        prop.L_2 = 1;
        prop.L__d_g_1 = 1;
        prop.L__d_g_2 = 1;
        prop.m__d_1 = 1;
        prop.m__d_2 = 1;
        prop.I__d_1 = 1;
        prop.I__d_2 = 1;
        prop.I__d_m_1 = 1;
        prop.I__d_m_2 = 1;
        prop.n_1 = 1;
        prop.n_2 = 1;
        prop.L__i_g_1 = 1;
        prop.L__i_g_2 = 1;
        prop.m__i_1 = 1;
        prop.m__i_2 = 1;
        prop.I__i_1 = 1;
        prop.I__i_2 = 1;

    case 'infant'
        prop.L__i_g_1 = 0.0525;
        prop.L__i_g_2 = 0.0385;
        prop.m__i_1 = 0.2474;
        prop.m__i_2 = 0.1786;
        prop.I__i_1 = 1.9471e-04;
        prop.I__i_2 = 7.3091e-05;

    case 'device'
        prop.L__d_g_1 = 0.044;
        prop.L__d_g_2 = 0.025;
        prop.m__d_1 = 0.032;
        prop.m__d_2 = 0.0067;
        prop.I__d_1 = 4.7e-5;
        prop.I__d_2 = 3.0e-6;
        prop.I__d_m_1 = 2e-5;
        prop.I__d_m_2 = 2e-5;
        prop.n_1 = 1;
        prop.n_2 = 1;

    case 'full'
        prop.L__i_g_1 = 0.0525;
        prop.L__i_g_2 = 0.0385;
        prop.m__i_1 = 0.2474;
        prop.m__i_2 = 0.1786;
        prop.I__i_1 = 1.9471e-04;
        prop.I__i_2 = 7.3091e-05;
        prop.L__d_g_1 = 0.044;
        prop.L__d_g_2 = 0.025;
        prop.m__d_1 = 0.032;
        prop.m__d_2 = 0.0067;
        prop.I__d_1 = 4.7e-5;
        prop.I__d_2 = 3.0e-6;
        prop.I__d_m_1 = 2e-5;
        prop.I__d_m_2 = 2e-5;
        prop.n_1 = 1;
        prop.n_2 = 1;
    case 'prototype'
        prop.L__i_g_1 = 0.044;
        prop.L__i_g_2 = 0.038;
        prop.m__i_1 = 0.308;
        prop.m__i_2 = 0.1607;
        prop.I__i_1 = 2.301e-04;
        prop.I__i_2 = 1.77e-04;
        prop.L__d_g_1 = 0.044;
        prop.L__d_g_2 = 0.025;
        prop.m__d_1 = 0.036;
        prop.m__d_2 = 0.0067;
        prop.I__d_1 = 4.7e-5;
        prop.I__d_2 = 3.0e-6;
        prop.I__d_m_1 = 2e-5;
        prop.I__d_m_2 = 2e-5;
        prop.n_1 = 1;
        prop.n_2 = 1;
end
end

