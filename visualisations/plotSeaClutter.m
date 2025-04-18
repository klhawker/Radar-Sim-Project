function plotSeaClutter(X, Y, seaClutter)
    % plotSeaClutter visualizes the sea clutter power using two plots:
    % 1. A pcolor plot of sea clutter (in dB) over the (X,Y) grid.
    % 2. A plot of average sea clutter (in dB) versus range.
    %
    % Inputs:
    %   X, Y      - Matrices representing the grid coordinates.
    %   seaClutter - Matrix of sea clutter power (linear scale).
    
    epsilon = 1e-30;
    seaClutter_dB = 10 * log10(seaClutter + epsilon);
    
    % Create a new figure with two subplots.
    figure;
    
    % %% Subplot 1: pcolor plot of Sea Clutter in dB
    % pcolor(X, Y, seaClutter_dB);
    % shading flat;
    % colormap('jet');
    % colorbar;
    % title('Power Returned From the Sea Clutter');
    % xlabel('X (m)');
    % ylabel('Y (m)');
    % axis equal;
    
    %% Subplot 2: Average Sea Clutter vs. Range
    % Compute the radial distance for each grid point.
    R = sqrt(X.^2 + Y.^2);
    
    % Assuming the grid was constructed in polar coordinates,
    % each row corresponds to a fixed range. We average clutter over azimuth.
    % For clarity, we use the first column (or unique row values) as the range vector.
    ranges = X(:,1);  % if X was built from meshgrid using polar coordinates, X(:,1) gives the range.
    % Alternatively, you can compute it as:
    % ranges = sqrt(X(:,1).^2 + Y(:,1).^2);
    
    % % Compute the average clutter in dB for each range (row).
    % avgSeaClutter_dB = mean(seaClutter_dB, 2);
    % 
    % subplot(1,2,2);
    % plot(ranges, avgSeaClutter_dB, 'b-o', 'LineWidth',1.5);
    % grid on;
    % xlabel('Range (m)');
    % ylabel('Average Sea Clutter (dB)');
    % title('Sea Clutter vs. Range');

    plot(ranges, seaClutter_dB(:,1), 'b-o', 'LineWidth',1.5);
    grid on;
    xlabel('Range (m)');
    ylabel('Sea Clutter (dB)');
    title('Sea Clutter vs. Range');
end
