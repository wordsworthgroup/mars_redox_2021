%% model_setup.m
% 
% set up model parameters
%
% checked by RDW 16/2/20

% set up t structure containing all time parameters
time.y2s      = 86400*365.25;        % 1 Earth year [s]
time.My       = 1e6*time.y2s;        % 1 million Earth years [s]
time.Gy       = 1e3;                 % 1 billion Earth years [My]
time.T        = 4.5e3-1;             % total time interval [My]
% -1 is simply to avoid NaN in climate interpolation when F >= Ftoday
time.nt       = 4.5e4-10;            % number of timesteps []
time.dt       = time.T/time.nt;      % time interval [My]
time.t_a      = (1:time.nt)*time.dt; % time array [My]       
time.tNoach   = 400;                 % start of Noachian (4.1 Ga) [My]

% set up p structure containing all fixed parameters
params.rho_ice  = 920;             % density of water ice (approx.) [kg/m3]
params.G        = 6.67408e-11;     % gravitational constant [m3/kg/s2]
params.muH2O    = 18.015/1e3;      % H2O molar mass [kg/mol]
params.muCO2    = 44.01/1e3;       % CO2 molar mass [kg/mol]
params.muH2     = 2.016/1e3;       % H2 molar mass [kg/mol]
params.muO2     = 31.999/1e3;      % O2 molar mass [kg/mol]
params.kB       = 1.3806503e-23;   % Boltzman constant [m2 kg s-2 K-1]
params.Rstar    = 8.314463;        % ideal gas constant [J/K/mol]
params.NA       = params.Rstar/params.kB; % Avogadro's constant [1/mol]
params.m_u      = 1.660539e-27;    % atomic mass unit [kg]
params.exa      = 1e18;            % exa conversion factor []
params.Tmelt    = 273.15;          % melting point of liquid water [K]
params.Tseas    = 263.15;          % melting point of liquid water (seasonal melting) [K]
params.fresh_i  = 1;               % re-calculate stochastic H input at each timestep?

% set up mars structure containing all planet-specific parameters
% https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
mars.M     = 0.64171e24;          % mass [kg]
mars.r     = 3389.5e3;            % radius [m]
mars.g     = params.G*mars.M/mars.r^2; % gravity [m/s/s]
mars.Area  = 4*pi*mars.r^2;       % surface area [m2]
mars.d     = 1.524;               % semi-major axis [AU]
mars.bar   = 1e5;                 % bar-->Pa. Should be in params structure, but I can't resist writing mars.bar.

% total moles reducing power in the H in 1 m GEL H2O (factor 2 from stoichometry)
mars.N_e        = 2*mars.Area*1.0*params.rho_ice/params.muH2O; % [mol]

% moles of CO2 given uCO2 in kg/m2 (bulk gas calculation)
mars.N_CO2      = @(uCO2) uCO2*mars.Area/params.muCO2; % [mol]

% H loss to space today [mol/My]
% Jakosky et al. 2018 Table 2 gives present-day rate as 1.6e26 to 11e26 atoms/s
mars.H_loss_now = time.My*1.0e27/params.NA; 




