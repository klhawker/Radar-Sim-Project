function plotRadarHeatMap(X, Y, radarData_dB, turbineSegments)
    % PLOTRADARHEATMAP visualizes the radar heat map along with turbine locations.
    
    figure;
    pcolor(X, Y, radarData_dB);
    shading flat;
    colormap('jet');
    colorbar;
    title('Combined Returns of Noise, Sea Clutter and Turbines (dB)');
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal;
    hold on;
    % Plot turbine segment positions as black circles
    plot(turbineSegments(:,1), turbineSegments(:,2), 'ko', 'MarkerSize', 2);
    legend('Turbine Segments');
    hold off;
end
