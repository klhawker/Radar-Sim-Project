function [turbinePositions, turbineSegments] = generateTurbineGridPositions(gridParams, turbineParams)
    % GENERATETURBINEGRIDPOSITIONS computes the turbine positions and
    % generates their segments (now only the turbine pole segments).
    
    num_rows = gridParams.num_rows;
    num_cols = gridParams.num_cols;
    S_dw = gridParams.S_dw;
    S_cw = gridParams.S_cw;
    wind_dir = gridParams.wind_dir;
    R_grid = gridParams.R_grid;
    theta_grid = gridParams.theta_grid;
    
    % Create rotation matrix for wind direction
    alpha_rad = deg2rad(wind_dir);
    R_rot = [cos(alpha_rad), -sin(alpha_rad); sin(alpha_rad), cos(alpha_rad)];
    
    % Global grid midpoint (with North aligned along positive Y)
    x_grid_mid = R_grid * sind(theta_grid);
    y_grid_mid = R_grid * cosd(theta_grid);
    
    turbinePositions = [];
    turbineSegments = [];
    
    for row = 1:num_rows
        downwind_offset = (row - 1) * S_dw - ((num_rows - 1) * S_dw)/2;
        for col = 1:num_cols
            crosswind_offset = (col - 1) * S_cw - ((num_cols - 1) * S_cw)/2;

            
            % Local grid coordinates
            local_pos = [downwind_offset; crosswind_offset];
            global_pos = R_rot * local_pos;
            
            % Global turbine position
            X = x_grid_mid + global_pos(1);
            Y = y_grid_mid + global_pos(2);
            turbinePositions = [turbinePositions; X, Y];
            
            % Calculate range and azimuth angle for this turbine
            R_dist = sqrt(X^2 + Y^2);
            pos_angle = atan2d(Y, X);
            
            % Generate only the pole segments for this turbine
            segments = generateTurbineSegmentsFarm(R_dist, pos_angle, wind_dir, ...
                turbineParams.pole_height, turbineParams.N_pole, [], [], [], [], []);
            turbineSegments = [turbineSegments; segments];
        end
    end
end
