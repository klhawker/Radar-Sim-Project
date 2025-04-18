function plotRadarHeatMapWithGhosts(X, Y, radarData, ghostReturnMatrix, ghostLocations, epsilon, turbineSegments, turbineReturns, radarNoise)
% PLOTRADARHEATMAPWITHGHOSTS overlays ghost returns onto the radar heat map by
% adding the multipath ghost heat map (in linear scale) to the radar data.
%
% Inputs:
%   X, Y             : Matrices of grid coordinates.
%   radarData        : Radar data in linear scale.
%   ghostReturnMatrix: Ghost heat map (multipath returns) in linear power.
%   epsilon          : Small constant to avoid log(0) (e.g., 1e-30).
%   turbineSegments  : Coordinates of turbine segments.
%   turbineReturns   : Turbine returns in linear scale.
%   radarNoise       : Radar noise in linear scale.
%
% The function computes the total radar power (radarData + ghost returns + turbine returns + radar noise),
% converts it to dB, and plots the resulting heat map.

    assignin("base","ghostReturnMatrix", ghostReturnMatrix);
    assignin("base","turbineReturns", turbineReturns);
    assignin("base","radarNoise", radarNoise);
    assignin("base","radarData", radarData);
    % Sum the ghost returns, turbine returns, and radar noise with the radar data.
    static_noise = ones(size(radarNoise)) * 1.6293e-14;
    totalRadarData = ghostReturnMatrix + turbineReturns + abs(radarNoise).^2;
    
    % Convert the total power to dB scale.
    totalRadarData_dB = 10 * log10(totalRadarData + epsilon);
    
    % Create the pseudocolor plot.
    figure;
    ph = pcolor(X, Y, totalRadarData_dB);
    shading flat;
    colormap('jet');
    colorbar;
    title('Radar Heat Map with Ghost Returns Added');
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal;
    hold on;
    
    % Plot ghost target locations as red circles.
    % if ~isempty(ghostLocations)
    %     plot(ghostLocations(:,1), ghostLocations(:,2), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    % end
    % plot radar at 0,0
    % Plot radar at 0,0 and capture its handle.
    hRadar = plot(0, 0, 'ko', 'MarkerSize', 5, 'LineWidth', 2, 'MarkerFaceColor','k');
    
    % Plot turbine segment positions and capture their handle.
    if ~isempty(turbineSegments)
        hTurbine = plot(turbineSegments(:,1), turbineSegments(:,2), 'ko', 'MarkerSize', 5, 'LineWidth', 1, 'MarkerFaceColor','white');
    else
        hTurbine = [];
    end
    
    % Create the legend only for these two handles.
    legend([hRadar, hTurbine], {'Radar', 'Turbine'}, 'Location','best');
    hold off;
    
    % Set the ZData of the surface to match its CData.
    % This ensures that the data cursor displays the actual dB values.
    ph.ZData = ph.CData;


end
