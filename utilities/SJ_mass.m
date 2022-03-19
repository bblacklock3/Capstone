function [m_1, m_2, I_1, I_2] = SJ_mass(age_weeks)
L_2 = 0.100; L_g_2 = 0.0385;
mass_thigh_SJ = 2.39e-01 + 1.24e-02*age_weeks;      %
mass_shank_SJ = 1.69e-01 + 4.16e-03*age_weeks;      % estimates from Sun and Jensen (1994)
mass_foot_SJ  = 5.22e-02 + 2.12e-03*age_weeks;      %
lgInertia_thigh =  8.44e-05 + 1.43e-05*age_weeks;   %
lgInertia_shank =  5.45e-05 + 2.41e-06*age_weeks;   % estimates from Sun and Jensen (1994)
lgInertia_foot  =  8.07e-06 + 1.15e-06*age_weeks;   %
m_1 = mass_thigh_SJ;
m_2 = mass_shank_SJ + mass_foot_SJ;
I_1 = lgInertia_thigh;
I_2 = lgInertia_shank + lgInertia_foot + mass_foot_SJ*(L_2-L_g_2)^2;
end

