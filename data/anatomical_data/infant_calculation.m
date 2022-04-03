%% Constants
% anthropometric data, from John's .xls file
age = 0.147945205479452;            % age in years
age_weeks = age * 365/7;
bodyMass = 3.4;                     % weight in kg

length_thigh = 0.108;               %
length_shank = 0.1;                 % segment lengths in m
length_foot  = 0.091;               %

circ_thigh = 0.187;                 %
circ_shank = 0.135;                 % segment circumference in m
circ_foot = 0.097;                  %

width_foot = 0.029;     

mass_thigh_XLS = 0.247359160273973;     %
mass_shank_XLS = 0.109713428;           % segment mass in kg, from the .xls file
mass_foot_XLS  = 0.079837736;           %

mass_thigh_SZ = 6.9126e-02*age      + 2.9582e-00*length_thigh + 3.1541e-00*circ_thigh - 6.7217e-01; %
mass_shank_SZ = 6.5138e-03*bodyMass + 1.8158e-00*length_shank + 1.8743e-00*circ_shank - 3.5460e-01; % estimates from Schneider and Zernicke (1992)
mass_foot_SZ  = 2.9331e-03*bodyMass + 1.2405e-00*length_foot  + 1.9337e-00*width_foot - 1.0250e-01; %

mass_thigh_SJ = 2.39e-01 + 1.24e-02*age_weeks; %
mass_shank_SJ = 1.69e-01 + 4.16e-03*age_weeks; % estimates from Sun and Jensen (1994)
mass_foot_SJ  = 5.22e-02 + 2.12e-03*age_weeks; %

mass_thigh = mass_thigh_SZ;
mass_shank = mass_shank_SZ;
mass_foot  = mass_foot_SZ;

com_thigh = 0.4859*length_thigh;    % distance of center of mass from proximal joint
com_shank = 0.4377*length_shank;    % estimation according to Schneider and Zernicke (1992)
com_foot = 0.3469*length_foot;      %

% estimates for moments of inertia
% according to Sun and Jensen (1994)
% units are kg m^2
lgInertia_thigh =  8.44e-05 + 1.43e-05*age_weeks
lgInertia_shank =  5.45e-05 + 2.41e-06*age_weeks
lgInertia_foot  =  8.07e-06 + 1.15e-06*age_weeks

% Real values
L_1 = length_thigh
L_2 = length_shank
G_1 = com_thigh
G_2 = (com_shank*mass_shank+com_foot*mass_foot)/(mass_shank+mass_foot)
m_1 = mass_thigh
m_2 = (mass_shank+mass_foot)
I_1 = lgInertia_thigh
I_2 = lgInertia_shank