function plotKSeaClutter(radarParams, clutterParams, max_radius, num_angles, nu)
    % plotKSeaClutter generates separate plots (each in its own window) for
    % K-distribution-based sea clutter, including:
    %   1. A pcolor plot (2D display) of clutter (in dB).
    %   2. A line plot of average clutter (in dB) versus range.
    %   3. A histogram of the clutter values (in dB).
    %   4. A 3D surface plot of clutter (in dB).
    %
    % Inputs:
    %   radarParams   - structure with radar parameters.
    %   clutterParams - structure with sea clutter parameters.
    %   max_radius    - maximum range to compute clutter (m).
    %   num_angles    - number of azimuth angles.
    %   nu            - shape parameter for the Gamma texture (K-distribution).
    
    % Calculate K-distributed sea clutter using your custom function:
    [P_cs_k, P_cs_k_dB, ranges] = computeSeaClutterK(radarParams, clutterParams, max_radius, num_angles, nu);
    
    %AWGN Noise
    zeross = ones(2667,1201);
    radarNoise = computeAWGN(radarParams.PW, zeross, 5);


    noise = abs(radarNoise).^2;
    noise_dB = 10*log10(noise);

    P_cs = load('/Users/kaihawker/Documents/MATLAB/Radar-Project/Modularised_PD_Map_Shadows/P_cs.mat');
    P_cs_dB = 10 * log10(P_cs.P_cs);
    % Retrieve necessary parameters for grid creation:
    C = radarParams.C;
    PW = radarParams.PW;
    cell_radius = (PW * C) / 2;
    % Define azimuth angles (in degrees) using the radar's BW_az parameter:
    angles = 0:radarParams.BW_az:360;
    angles_shifted = mod(angles, 360);
    
    % Create a meshgrid for the polar coordinate grid (using 'ranges' as row vector):
    [Theta, R_grid] = meshgrid(deg2rad(angles_shifted), ranges);
    [X, Y] = pol2cart(Theta, R_grid);
    
    
    %% Figure 1: Average clutter vs. range
    % Average the clutter power (in dB) over all azimuth angles for each range.
    Kclutter = P_cs_k_dB(:,1);
    
    Kclutter_masked = Kclutter;
    Kclutter_masked(Kclutter < noise_dB(:,1)) = NaN;

    line_mask = P_cs_dB;
    line_mask(P_cs_dB < -137.8803) = NaN;

    figure;
    plot(ranges, Kclutter_masked, 'b:o', 'LineWidth', 2, 'Marker','none'); 
    hold on;
    plot(ranges, line_mask, 'b', 'LineWidth',2);
    hold on;
    plot(ranges, noise_dB(:,1), 'r:', 'LineWidth', 2); % Red dashed line for noise floor
    yline(-137.8803, 'r', 'LineWidth', 2); % Optional horizontal reference
    
    grid on;
    xlabel('Range (m)');
    ylabel('Power (dB)');
    legend('K Sea Clutter', 'Empirical Sea Clutter', 'AWGN Thermal Noise', 'Thermal Noise');
    xlim([0 10000]);
    ylim([-160 -80]);

    max_dB = max(P_cs_k_dB, noise_dB);

figure;

% Define zoom window (adjust as needed)
xZoom = [-5000, 5000]; % X range of zoomed-in area
yZoom = [-5000, 5000]; % Y range of zoomed-in area


%% Figure 1: Full range with zoom box
figure;
pcolor(X, Y, max_dB);
shading flat;
colormap('winter');
colorbar;
title('Full Range - K-Distribution Sea Clutter (dB)');
xlabel('X (m)');
ylabel('Y (m)');
axis equal;
xlim([-100000 100000]);
ylim([-100000 100000]);

% Draw rectangle showing zoom region
hold on;
rectangle('Position', [xZoom(1), yZoom(1), diff(xZoom), diff(yZoom)], ...
          'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
hold off;

%% Figure 2: Zoomed-in view
figure;
pcolor(X, Y, max_dB);
shading flat;
colormap('winter');
colorbar;
title('Zoomed In - K-Distribution Sea Clutter (dB)');
xlabel('X (m)');
ylabel('Y (m)');
axis equal;
xlim(xZoom);
ylim(yZoom);
% Add RADAR marker
hold on;
plot(0, 0, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % white circle with black fill
text(0, 0, '  RADAR', 'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10, 'VerticalAlignment', 'bottom');
hold off;


    
    %% Figure 3: Histogram of K-distributed clutter values
    % Flatten the clutter dB values to build a histogram.
    % clutterValues = P_cs_k_dB(:);
    % figure;
    % histogram(clutterValues, 50);  % Using 50 bins; adjust as needed
    % grid on;
    % xlabel('Sea Clutter (dB)');
    % ylabel('Frequency');
    % title('Histogram of K-Distributed Sea Clutter (dB)');
    % 
    % %% Figure 4: 3D Surface Plot of K-distributed Sea Clutter
    % figure;
    % surf(X, Y, P_cs_k_dB);
    % shading interp;
    % colormap('jet');
    % colorbar;
    % title('3D Surface Plot of K-Distributed Sea Clutter (dB)');
    % xlabel('X (m)');
    % ylabel('Y (m)');
    % zlabel('Clutter (dB)');
    % view(45,30);
    % axis tight;
end
