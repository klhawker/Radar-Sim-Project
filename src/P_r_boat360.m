function P_r_rotated = P_r_boat360(radii, h_t, h_r, sigma, cell_radius, reflection_coef, num_angles)
    % P_R_BOAT360 computes the rotated power returned matrix for an additional object.
    
    P_peak = 25000;
    f = 9e9;
    C = 3e8;
    lamb = C / f;
    PW = 250e-9;
    G_t = 1000;
    G_r = 1000;
    D_el = 0.5;
    
    ranges = (radii(1:end-1) + radii(2:end)) / 2;
    ranges = [ranges, (cell_radius/2) + ranges(end)];
    
    d_d = sqrt(h_r^2 + ranges.^2);
    d_i = 2 * sqrt(ranges.^2 + (h_r + h_t)^2);
    
    phase_d = (2*pi*d_d)/lamb;
    phase_i = (2*pi*d_i)/lamb;
    
    V_d = 1;
    V_i = V_d * reflection_coef;
    phase_d_complex = V_d * exp(-1j * phase_d);
    phase_i_complex = V_i * exp(-1j * phase_i);
    
    prop_factor = abs(phase_d_complex + phase_i_complex);
    
    elevation_angle = atand(h_t./ranges);
    u_el = (D_el*pi*sind(elevation_angle))/lamb;
    g_elevation = sinc(u_el / pi);
    G_elevation = abs(g_elevation).^2;
    
    G_prime = G_t * G_r .* G_elevation.^2;
    P_r = (P_peak .* G_prime .* lamb^2 .* sigma .* prop_factor.^4) ./ ((4*pi)^3 * ranges.^4);
    epsilon = 1e-30;
    P_r_rotated = repmat(P_r(:), 1, num_angles);
end
