function plotPDThreshold(X, Y, PD_approx)
    % PLOTDETECTIONPROBABILITY displays the PD with thershold using a colormap.
    threshold = 0.7;
        % Create a binary mask for highlighting
    PD_highlight = PD_approx > threshold;
    
    % Create a custom colormap: black for PD <= 0.7, green for PD > 0.7
    custom_colormap = [0, 0, 0;  % Black for PD <= 0.7
                       0, 1, 0]; % Green for PD > 0.7
    
    figure;
    pcolor(X, Y, double(PD_highlight)); % Visualize binary mask
    shading flat;
    colormap(custom_colormap);
    colorbar('Ticks', [0, 1], 'TickLabels', {'PD \leq 0.7', 'PD > 0.7'});
    title('Highlighted Probability of Detection (PD)');
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal;
    hold off;
end