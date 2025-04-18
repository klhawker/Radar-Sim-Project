function [P_cs_rotated, P_cs_dB_rotated, ranges] = calculateSeaClutter(radarParams, clutterParams, max_radius, num_angles)
    % CALCULATESEACLUTTER Computes sea clutter power for a radar system.
    
    C = radarParams.C;
    f = radarParams.f;
    PW = radarParams.PW;
    P_peak = radarParams.P_peak;
    G_t = radarParams.G_t;
    G_r = radarParams.G_r;
    D_el = radarParams.D_el;
    epsilon = radarParams.epsilon;
    
    lamb = C / f;
    B = 1 / PW;
    cell_radius = (PW * C) / 2;
    radii = 0:cell_radius:max_radius;
    
    % Compute range midpoints
    ranges = (radii(1:end-1) + radii(2:end)) / 2;
    ranges = [ranges, (cell_radius/2) + ranges(end)];
  
    
    % Environmental parameters
    s_s = clutterParams.s_s;
    v_wind = 3.16 * s_s^0.8;
    h_r = radarParams.h_r;
    h_a = 4.52e-3 * v_wind^2.5;
    
    % Grazing angle
    grazing_angle = asin(h_r ./ ranges) + epsilon;
    grazing_angle(1) = atan(h_r/ranges(1));
    assignin('base', 'grazing', (abs(grazing_angle)));

    
    % Roughness and gain factors
    roughness_factor = ((14.4 * lamb + 5.5) * grazing_angle .* h_a) / lamb;
    G_a = (roughness_factor.^4) ./ (1 + roughness_factor.^4);
    G_u = exp(0.2 * cosd(clutterParams.azimuth_direction) * (1 - 2.8 * grazing_angle) * (lamb + 0.015)^-0.4);
    G_w = ((1.94 * v_wind) / (1 + (v_wind / 15.4)))^clutterParams.q;
    
    % Clutter reflectivity and cross-section
    clutter_reflectivity_dB = 10 * log10(abs(3.9e-6 * lamb .* grazing_angle.^0.4 .* G_a .* G_u .* G_w) + epsilon);
    clutter_reflectivity = 10.^(clutter_reflectivity_dB / 10);
    theta_az = deg2rad(radarParams.BW_az);
    clutter_cross_section = clutter_reflectivity .* ranges .* theta_az .* (C ./ (2 * B * cos(grazing_angle)));
    
    % Elevation gain
    u_el = (D_el * pi * sin(grazing_angle)) / lamb;
    g_elevation = sinc(u_el / pi);
    
    % Clutter power (linear scale) and conversion to dB
    P_cs = (P_peak * G_t * G_r * 1 .* g_elevation.^4 .* lamb^2 .* clutter_cross_section) ./ ((4*pi)^3 .* ranges.^4);
    P_cs_dB = 10 * log10(max(P_cs, epsilon));
    assignin("base",'P_cs', P_cs);

    % Replicate across the number of azimuth angles
    P_cs_rotated = repmat(P_cs(:), 1, num_angles);
    P_cs_dB_rotated = repmat(P_cs_dB(:), 1, num_angles);
end