function [radarData, radarData_dB, PD_approx, object_Pr_dB, X, Y, seaClutter, turbineReturns, radarNoise, noise_and_clutter] = radarSimulation(params, turbineSegments)
    % RADARSIMULATION computes the radar returns including noise, sea clutter,
    % and power from turbine segments (and an additional object) using a coherent
    % sum (per turbine) of the segment contributions.
    
    % Unpack radar parameters
    rp = params.radarParams;
    epsilon = rp.epsilon;

    radar_pos = params.shadowParams.radar_pos;
    
    % Compute range resolution and create range vector (grid)
    cell_radius = (rp.PW * rp.C) / 2;
    radii = 0:cell_radius:rp.max_radius;
    angles = 0:rp.BW_az:360;
    angles_shifted = mod(angles, 360);
    [Theta, R_grid_mesh] = meshgrid(deg2rad(angles_shifted), radii);
    [X, Y] = pol2cart(Theta, R_grid_mesh);
    
    % Noise calculation
    % k = 1.380649e-23;
    % T = 295;
    % tau = rp.PW;
    % B = 1/tau;
    % N_noise = k * T * B;
    % radarData = N_noise * ones(size(R_grid_mesh));
    % radarNoise = radarData;

    % Noise calculation using AWGN model
    k = 1.380649e-23; % Boltzmann's constant
    T = 295;          % Temperature in Kelvin
    tau = rp.PW;      % Pulse width
    B = 1/tau;        % Receiver bandwidth
    N0 = k * T * B;   % Thermal noise power
    
    % Noise Figure??
    NF_dB = 5; 
    NF = 10^(NF_dB/10);
    N0 = N0 * NF;
    
    % Generate complex AWGN noise:
    radarData = sqrt(N0/2) * (randn(size(R_grid_mesh)) + 1i * randn(size(R_grid_mesh)));
    radarNoise = radarData;
    disp(max(10*log10(abs(radarNoise.^2)), [],"all"));
    
    % Sea clutter
    num_angles = size(R_grid_mesh, 2);
    [seaClutter, ~, ~] = calculateSeaClutter(rp, params.clutterParams, rp.max_radius, num_angles);
    nu = 2;
    [P_cs_k, P_cs_k_dB, ranges] = computeSeaClutterK(rp, params.clutterParams, rp.max_radius, num_angles, nu);
    radarData = radarData + seaClutter;
    noise_and_clutter = abs(radarNoise).^2 + P_cs_k;

    
    %% Coherent processing per turbine
    % Assume each turbine has N_pole segments
    N_pole = params.turbineParams.N_pole;
    num_segments = size(turbineSegments, 1);
    N_turbines = num_segments / N_pole;
    
    lambda = rp.lamb;
    pole_radius = params.turbineParams.pole_radius;
    pole_height = params.turbineParams.pole_height;
    
    % Preallocate per-turbine arrays
    sigma_total_array = zeros(N_turbines, 1);  % Coherent RCS per turbine
    avgRange_array = zeros(N_turbines, 1);       % Average (3D) range from radar
    prop_factor_array = zeros(N_turbines, 1);    % Propagation factor per turbine
    turbine_base = zeros(N_turbines, 2);         % (X,Y) position for each turbine
    turbine_avg_z = zeros(N_turbines, 1);        % Average z (height) for each turbine
    object_real_dist_turbine = zeros(N_turbines, 1); % Radar-to-turbine range (from base)
    
    % Loop over turbines
    for t = 1:N_turbines
        idx_start = (t-1)*N_pole + 1;
        idx_end   = t*N_pole;
        segs = turbineSegments(idx_start:idx_end, :);  % [N_pole x 3]
        
        % Assume all segments in a turbine share the same (X,Y) (turbine base)
        turbine_base(t, :) = segs(1, 1:2);
        turbine_avg_z(t) = mean(segs(:,3));
        
        % Coherent summation over segments (using full 3D range adjusted to surface)
        sum_Re = 0;
        sum_Im = 0;
        distances = zeros(N_pole, 1);
        for s = 1:N_pole
            x_seg = segs(s, 1);
            y_seg = segs(s, 2);
            z_seg = segs(s, 3);

            % Compute the horizontal distance from radar to segment center
            d_center_xy = sqrt(x_seg^2 + y_seg^2);
            % Subtract the pole radius to get the horizontal distance to the surface
            d_surface_xy = max(d_center_xy - pole_radius, 0);
            % Compute the 3D distance from the radar ([0,0,rp.h_r]) to the surface point
            d_seg = sqrt(d_surface_xy^2 + (z_seg - rp.h_r)^2);
            distances(s) = d_seg;

            % Compute phase for this segment using the surface distance
            gamma_s = (2*pi*d_seg) / lambda;

            % Assuming monostatic: phi1 = phi2
            phi = atan2(z_seg - rp.h_r, d_surface_xy);

            % Compute per-segment RCS from the cylinder model
            sigma_seg = computeBistaticRCS_PO(2*pi/lambda, pole_radius, pole_height, 0, phi, phi);
            % Convert to field amplitude and apply phase factor
            field_amp = sqrt(sigma_seg);
            sigma_complex = field_amp * exp(-1j * gamma_s);
            sum_Re = sum_Re + real(sigma_complex);
            sum_Im = sum_Im + imag(sigma_complex);
        end
        sigma_total_array(t) = sqrt(sum_Re^2 + sum_Im^2);


        avgRange_array(t) = mean(distances);
        
        % For the radar equation, compute a representative radar-to-turbine range.
        % Use the turbine base (X,Y) and average z.
        x_base = turbine_base(t,1);
        y_base = turbine_base(t,2);
        obj_dist = sqrt(x_base^2 + y_base^2);
        delta_h = turbine_avg_z(t) - rp.h_r;
        object_real_dist_turbine(t) = sqrt(obj_dist^2 + delta_h^2);
        
        % Compute propagation factor (using method of images)
        delta_h_image = rp.h_r + turbine_avg_z(t);
        d_i = sqrt(obj_dist^2 + delta_h_image^2);
        phase_d = (2*pi*object_real_dist_turbine(t)) / lambda;
        phase_i = (2*pi*d_i) / lambda;
        V_d = 1; V_i = 0.7;
        phase_d_complex = V_d * exp(-1j * phase_d);
        phase_i_complex = V_i * exp(-1j * phase_i);
        prop_factor_array(t) = abs(phase_d_complex + phase_i_complex);
    end
    
    %% Compute gains and map turbine returns onto the radar grid
    % Preallocate grid for turbine return power
    P_r = zeros(size(R_grid_mesh));
    beam_angles = Theta(1, :); % vector of azimuth angles from the grid
    % Loop over turbines to compute and map their return power
    for t = 1:N_turbines
        % For this turbine, get base (X,Y) and compute its azimuth angle
        turbine_az = atan2(turbine_base(t,2), turbine_base(t,1));
        % Relative azimuth for each beam
        rel_az_angle = wrapToPi(beam_angles - turbine_az);
        % Limit relative angles to [-pi/2, pi/2]
        rel_az_angle(rel_az_angle < -pi/2) = -pi/2;
        rel_az_angle(rel_az_angle >  pi/2) =  pi/2;
        % Azimuth gain
        u_az = (rp.D_az * pi * sin(rel_az_angle)) / lambda;
        g_az = sinc(u_az/pi);
        G_az = abs(g_az).^2; 
        
        % Elevation gain (using the turbine’s average elevation)
        el_angle = atan2(turbine_avg_z(t) - rp.h_r, sqrt(turbine_base(t,1)^2 + turbine_base(t,2)^2));
        u_el = (rp.D_el * pi * sin(el_angle)) / lambda;
        g_el = sinc(u_el/pi);
        G_el = abs(g_el)^2;  % scalar
        
        % Total antenna gain for this turbine (applied over the azimuth beam)
        G_t_prime = rp.G_t * G_el * G_az;  % row vector
        G_r_prime = rp.G_r * G_el * G_az;
        
        % Use the monostatic radar equation for this turbine:
        % Note: We use the computed object_real_dist_turbine(t) as the range.
        P_r_turbine = (rp.P_peak * G_t_prime .* G_r_prime * lambda^2 * sigma_total_array(t) * (prop_factor_array(t)^4)) ...
                      / ((4*pi)^3 * (object_real_dist_turbine(t)^4));
                  
        % Map this turbine’s return power onto the radar grid.
        % Determine the range bin corresponding to this turbine.
        idx_range = floor(object_real_dist_turbine(t) / cell_radius) + 1;
        if idx_range <= size(P_r, 1)
            % Add the turbine’s power (a 1 x n_angles row vector) to the grid row.
            P_r(idx_range, :) = P_r(idx_range, :) + P_r_turbine;
        end
    end
    turbineReturns = P_r;
    assignin("base", "turbineReturns", turbineReturns);
    [i, j, t_val] = find(turbineReturns);
    non_zero_turbines = 10 * log10(t_val);
    assignin("base", "non_zero_turbines", non_zero_turbines);
    fprintf("Max Turbines: %f\n", max(non_zero_turbines));
    fprintf("Mean Turbines: %f\n", mean(non_zero_turbines));
    
    % Add turbine returns to the radar data
    radarData = radarData + P_r;
    noise = abs(radarData.^2);
    radarData_dB = 10 * log10(noise + epsilon);
    
    %% Additional object return (e.g. boat)
    object_Pr = P_r_boat360(radii, params.objectParams.h_t, rp.h_r, ...
                              params.objectParams.object_sigma, cell_radius, 0.7, num_angles);
    radarData = radarData + object_Pr;
    object_Pr_dB = 10 * log10(object_Pr + epsilon);
    
    %% SNR and Detection Probability (Gaussian approximation)
    SNR = object_Pr ./ noise;
    threshold = sqrt(-2 * log(1e-6));
    PD_approx = 0.5 * erfc((threshold - sqrt(2 * SNR)) / sqrt(2));
end
