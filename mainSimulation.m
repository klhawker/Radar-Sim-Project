function mainSimulation()
    % MAINSIMULATION orchestrates the simulation, from configuration to plotting.
    
    clear; clc;
    
    % Load configuration parameters
    params = simulationConfig();
    
    % Generate turbine grid positions and segments
    [turbinePositions, turbineSegments] = generateTurbineGridPositions(params.gridParams, params.turbineParams);
    
     
    % Run the radar simulation using turbine segments and clutter
    [radarData, radarData_dB, PD_approx, object_Pr_dB, X, Y, seaClutter, turbineReturns, radarNoise, noise_and_clutter] = radarSimulation(params, turbineSegments);
    
    % Calc and get shadow mask  
    pole_radius = params.turbineParams.pole_radius;
    shadowMask = calculateShadowCones(turbinePositions, X, Y, params.shadowParams, pole_radius);

    % Calculate multipath ghost heat map and ghost locations
    %[multipathHeatMap, ghostLocations, debugInfo] = complete_ghosts(params, turbinePositions, X, Y, turbineSegments, true);


    nu = 2;  % example shape parameter; adjust based on sea state
    plotKSeaClutter(params.radarParams, params.clutterParams, params.radarParams.max_radius, 1201, nu);
    %plotSeaClutterPlots(params.radarParams, params.clutterParams, 100000, 1201)
    %Plot the radar heat map, detection probability, and object return
    % plotRadarHeatMap(X, Y, radarData_dB, turbineSegments);
    % plotDetectionProbability(X, Y, PD_approx);
    % plotPDThreshold(X, Y, PD_approx)
    % plotObjectReturn(X, Y, object_Pr_dB);
    %plotSeaClutter(X, Y, seaClutter);
    %plotTurbineReturns(X, Y, turbineReturns, turbineSegments, radarNoise);
    % plotPDShadows(X, Y, PD_approx, shadowMask, turbineSegments);
    %plotRadarHeatMapWithGhosts(X, Y, radarData, multipathHeatMap, ghostLocations, params.radarParams.epsilon, turbineSegments, turbineReturns, radarNoise);
    %plotBistaticRCS();
    %plotNoise(params);

    %     %% Create the PPI-Style Display
    % figure;
    % pcolor(X, Y, radarData_dB);
    % shading flat;
    % colormap('jet');
    % colorbar;
    % title('Radar Display (PPI Style)');
    % xlabel('X (m)');
    % ylabel('Y (m)');
    % axis equal;
    % hold on;
    % 
    % % Place a filled black marker at (0,0) and add the label "Radar"
    % plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    % text(10, 10, 'Radar', 'Color', 'k', 'FontWeight', 'bold');  % adjust the offset (10,10) as needed
    % hold off;
end
    