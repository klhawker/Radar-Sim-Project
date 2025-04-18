function plotTurbineReturns(X, Y, turbineReturns, turbineSegments,radarNoise)
    % PLOTOBJECTRETURN visualizes the power returned from turbines.
    radarNoise = abs(radarNoise.^2);
    epsilon = 1e-30;
    turbineReturns_dB = 10 * log10(turbineReturns + radarNoise + epsilon);
    figure;
    ph = pcolor(X, Y, turbineReturns_dB);
    shading flat;
    colormap('jet');
    colorbar;
    title('Power Returned From the Turbines');
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal;
    hold on;
        % Plot turbine segment positions as black circles
    plot(turbineSegments(:,1), turbineSegments(:,2), 'ko', 'MarkerSize', 1);
    hold off;
    ph.ZData = ph.CData;
end