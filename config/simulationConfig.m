function params = simulationConfig()
    % SIMULATIONCONFIG Returns a struct with all simulation parameters.
    
    %% Radar Parameters
    params.radarParams.BW_az = 0.169 * 2;
    params.radarParams.BW_az = round(params.radarParams.BW_az,1);
    params.radarParams.BW_az_rad = deg2rad(params.radarParams.BW_az);
    params.radarParams.D_az = 5;            % Azimuth radar length (m)
    params.radarParams.D_el = 0.5;          % Elevation radar height (m)
    params.radarParams.PW = 250e-9;         % Pulse width (s)
    params.radarParams.C = 3e8;             % Speed of light (m/s)
    params.radarParams.P_peak = 25000;      % Peak transmitted power (W)
    params.radarParams.f = 9e9;             % Frequency (Hz)
    params.radarParams.lamb = params.radarParams.C / params.radarParams.f;
    params.radarParams.G_t = 1000;
    params.radarParams.G_r = 1000;
    params.radarParams.epsilon = 1e-30;
    params.radarParams.h_r = 20;             % Radar height (m)
    params.radarParams.max_radius = 100000; % Maximum radar range (m)
    params.radarParams.radar_pos = [0, 0];
    
    %% Clutter and Environmental Parameters
    params.clutterParams.L_t = 1;
    params.clutterParams.L_bs = 1;
    params.clutterParams.L_atm = 1;
    params.clutterParams.L_r = 1;
    params.clutterParams.G_pc = 1;
    params.clutterParams.s_s = 3;
    params.clutterParams.BW_el = 1.5;
    params.clutterParams.azimuth_direction = 0;
    params.clutterParams.q = 1.1 / (params.radarParams.lamb + 0.015)^0.4;
    
    %% Object Parameters (for additional targets, e.g. a boat)
    params.objectParams.object_sigma = 10;
    params.objectParams.h_t = 10; % obeject height
    
    %% Grid Configuration
    params.gridParams.num_rows = 10;           % Number of turbines (downwind) 10
    params.gridParams.num_cols = 5;            % Number of turbines (crosswind) 5
    params.gridParams.blade_length = 80;       % Blade length (m)
    params.gridParams.S_cw = 4 * params.gridParams.blade_length; % Crosswind spacing (m) 4
    params.gridParams.S_dw = 4 * params.gridParams.blade_length; % Downwind spacing (m) 8
    params.gridParams.wind_dir = 55;           % Wind direction (deg clockwise from North)
    params.gridParams.R_grid = 3000;           % Distance from radar to grid midpoint (m)
    params.gridParams.theta_grid = 20;         % Azimuth angle from North to grid midpoint (deg) 17.74
    
    %% Turbine Parameters
    params.turbineParams.pole_height = 100;         % Pole height (m)
    params.turbineParams.N_pole = 100;              % Number of pole segments
    params.turbineParams.pole_radius = 5;
    %% Shadow Paramters
    params.shadowParams.radar_pos = [0, 0];
    params.shadowParams.shadow_length = 60000;
end
