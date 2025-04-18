function plotPDShadows(X, Y, PD_approx, shadowMask, turbineSegments)
% PLOTPDWITHSHADOWS plots the probability of detection (PD) thresholded at 0.7,
% and marks regions in shadow (from turbine poles) as undetectable.
%
% Inputs:
%   X, Y         : Matrices of grid coordinates.
%   PD_approx    : Matrix of probability of detection values.
%   shadowMask   : Binary mask (same size as X and Y) where 1 indicates a shadowed area.
%
% The function applies a threshold of 0.7 to PD_approx and then forces any grid
% cell that is shadowed to be marked as undetectable (i.e. below threshold).

    % Define the PD threshold.
    threshold = 0.7;
    
    % Create a binary PD mask: true where PD > threshold.
    PD_highlight = PD_approx > threshold;
    
    % Force shadowed regions to be undetectable.
    % That is, if a grid cell is in shadow, mark it as below threshold.
    PD_highlight(shadowMask == 1) = false;
    
    % Define a custom colormap:
    % - Black (0,0,0) for undetectable (PD <= threshold or shadowed)
    % - Green (0,1,0) for detectable (PD > threshold)
    custom_colormap = [0, 0, 0; 0, 1, 0];
    
    % Plot the combined PD mask.
    figure;
    pcolor(X, Y, double(PD_highlight));
    shading flat;
    colormap(custom_colormap);
    colorbar('Ticks', [0, 1], 'TickLabels', {'Undetectable','Detected'});
    title('Probability of Detection (PD) with Shadows (Threshold = 0.7)');
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal;
    hold on;

    if ~isempty(turbineSegments)
        plot(turbineSegments(:,1), turbineSegments(:,2), 'ro', 'MarkerSize', 5, 'LineWidth', 1);
        legend('PD', 'Turbines', 'Location','best');
    end
    hold off;
end
