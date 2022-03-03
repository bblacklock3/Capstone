function prop = load_prop(name)
prop = [];
prop.g = 9.81;

switch (name)
    case 'ones'
        prop.L_1 = 1;
        prop.L_2 = 1;
        prop.L__d_g1 = 1;
        prop.L__d_g2 = 1;
        prop.m__d_1 = 1;
        prop.m__d_2 = 1;
        prop.I__d_1 = 1;
        prop.I__d_2 = 1;
        prop.I__d_m1 = 1;
        prop.I__d_m2 = 1;
        prop.n_1 = 1;
        prop.n_2 = 1;
        prop.L__i_g1 = 1;
        prop.L__i_g2 = 1;
        prop.m__i_1 = 1;
        prop.m__i_2 = 1;
        prop.I__i_1 = 1;
        prop.I__i_2 = 1;

    case 'infant_only'
        prop.L_1 = 0.1080;
        prop.L_2 = 0.1000;
        prop.L__i_g1 = 0.0525;
        prop.L__i_g2 = 0.0385;
        prop.m__i_1 = 0.2474;
        prop.m__i_2 = 0.1786;
        prop.I__i_1 = 2.9576e-04;
        prop.I__i_2 = 1.6268e-05;
        prop.L__d_g1 = 0;
        prop.L__d_g2 = 0;
        prop.m__d_1 = 0;
        prop.m__d_2 = 0;
        prop.I__d_1 = 0;
        prop.I__d_2 = 0;
        prop.I__d_m1 = 0;
        prop.I__d_m2 = 0;
        prop.n_1 = 1;
        prop.n_2 = 1;

    case 'device_only'
        prop.L_1 = 0.1080;
        prop.L_2 = 0.1000;
        prop.L__i_g1 = 0;
        prop.L__i_g2 = 0;
        prop.m__i_1 = 0;
        prop.m__i_2 = 0;
        prop.I__i_1 = 0;
        prop.I__i_2 = 0;
        prop.L__d_g1 = 0.05;
        prop.L__d_g2 = 0.05;
        prop.m__d_1 = 0.300;
        prop.m__d_2 = 0.150;
        prop.I__d_1 = 3e-4;
        prop.I__d_2 = 1e-4;
        prop.I__d_m1 = 1e-5;
        prop.I__d_m2 = 1e-5;
        prop.n_1 = 1;
        prop.n_2 = 1;

    case 'full_system'
        prop.L_1 = 0.1080;
        prop.L_2 = 0.1000;
        prop.L__i_g1 = 0.0525;
        prop.L__i_g2 = 0.0385;
        prop.m__i_1 = 0.2474;
        prop.m__i_2 = 0.1786;
        prop.I__i_1 = 2.9576e-04;
        prop.I__i_2 = 1.6268e-05;
        prop.L__d_g1 = 0.05;
        prop.L__d_g2 = 0.05;
        prop.m__d_1 = 0.300;
        prop.m__d_2 = 0.150;
        prop.I__d_1 = 3e-4;
        prop.I__d_2 = 1e-4;
        prop.I__d_m1 = 1e-5;
        prop.I__d_m2 = 1e-5;
        prop.n_1 = 1;
        prop.n_2 = 1;        

    otherwise
        prop.L_1 = 0;
        prop.L_2 = 0;
        prop.L__i_g1 = 0;
        prop.L__i_g2 = 0;
        prop.m__i_1 = 0;
        prop.m__i_2 = 0;
        prop.I__i_1 = 0;
        prop.I__i_2 = 0;
        prop.L__d_g1 = 0;
        prop.L__d_g2 = 0;
        prop.m__d_1 = 0;
        prop.m__d_2 = 0;
        prop.I__d_1 = 0;
        prop.I__d_2 = 0;
        prop.I__d_m1 = 0;
        prop.I__d_m2 = 0;
        prop.n_1 = 0;
        prop.n_2 = 0;

end
end

