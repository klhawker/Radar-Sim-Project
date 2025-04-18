function shadowMask = calculateShadowCones(turbinePositions, X, Y, shadowParams, pole_radius)
% CALCULATESHADOWCONES computes a binary mask for grid points in the shadow 
% of turbine poles using a geometrically calculated half-angle.
%
% Inputs:
%   turbinePositions : Nx2 matrix with [X, Y] coordinates of each turbine pole center.
%   X, Y             : Matrices of grid coordinates (e.g., from pol2cart).
%   shadowParams     : Structure with fields:
%                        - radar_pos: 1x2 vector [x, y] for the radar's position.
%                        - shadow_length: Maximum distance (m) that the shadow extends.
%   pole_radius      : The physical radius of each turbine pole.
%
% Output:
%   shadowMask       : A matrix (same size as X and Y) with 1's indicating shadowed zones,
%                      and 0's elsewhere.

    % Extract required shadow parameter
    radar_pos = shadowParams.radar_pos;
    shadow_length = shadowParams.shadow_length;
    
    % Initialize the mask
    shadowMask = false(size(X));
    
    % Loop over each turbine pole
    for i = 1:size(turbinePositions, 1)
        % Turbine center coordinates
        turbine_center = turbinePositions(i, :);
        
        % Compute the direction vector from the radar to the turbine center.
        v = turbine_center - radar_pos;
        vNorm = norm(v);
        if vNorm == 0
            continue; % Skip if the turbine is at the radar location.
        end
        
        % Unit vector along the cone axis (from radar to turbine)
        d = v / vNorm;
        
        % Calculate the geometrical half-angle for this turbine
        theta_min = deg2rad(0.25); % For example, minimum 0.25°
        theta_geo = max(atan(pole_radius / vNorm), theta_min);
        %disp(rad2deg(theta_geo));
        
        % A perpendicular unit vector in 2D
        d_perp = [-d(2), d(1)];
        
        % Shift grid coordinates so that the turbine center is at the origin.
        X_rel = X - turbine_center(1);
        Y_rel = Y - turbine_center(2);
        
        % Project each grid point onto the cone axis:
        x_proj = X_rel * d(1) + Y_rel * d(2);
        
        % Compute the perpendicular distance from the cone axis:
        y_proj = X_rel * d_perp(1) + Y_rel * d_perp(2);
        y_abs = abs(y_proj);
        
        % Compute the allowed half-width at a given x_proj:
        % Starting width is pole_radius and increases linearly with x_proj.
        allowed_width = pole_radius + x_proj * tan(theta_geo);
        
        % A grid point is in the shadow if:
        %   1. It lies behind the turbine (x_proj > 0),
        %   2. Its distance along the axis is less than the shadow length,
        %   3. Its perpendicular distance is within the allowed width.
        inShadow = (x_proj > 0) & (x_proj <= shadow_length) & (y_abs <= allowed_width);
        
        % Combine with the overall shadow mask.
        shadowMask = shadowMask | inShadow;
    end
    
    % Convert the logical mask to double (0's and 1's)
    shadowMask = double(shadowMask);
end
