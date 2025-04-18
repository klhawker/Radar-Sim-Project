function turbine_segments = generateTurbineSegmentsFarm(R, pos_angle, dir_angle, pole_height, N_pole, ~, ~, ~, ~, ~)
    % GENERATETURBINSEGMENTSFARM generates only the pole segments for a wind turbine.
    %
    % Inputs:
    %   R         - Distance of the turbine from the radar (m)
    %   pos_angle - Azimuth angle of the turbine relative to North (degrees)
    %   dir_angle - (Unused here; kept for interface compatibility)
    %   pole_height - Height of the turbine pole (m)
    %   N_pole    - Number of pole segments
    %
    % Outputs:
    %   turbine_segments - Matrix of turbine pole segment coordinates (N_pole x 3)
    
    % Calculate turbine base position
    turbine_x = R * cosd(pos_angle);
    turbine_y = R * sind(pos_angle);
    
    % Generate only the pole segments
    z_pole = ((1:N_pole) - 0.5) * (pole_height / N_pole);  % Divide the pole height evenly
    x_pole = turbine_x * ones(1, N_pole);  % All segments share the same X position
    y_pole = turbine_y * ones(1, N_pole);  % All segments share the same Y position
    
    turbine_segments = [x_pole', y_pole', z_pole'];
    assignin('base','turbine_segments',turbine_segments);
end
