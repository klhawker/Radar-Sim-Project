function plotDetectionProbability(X, Y, PD_approx)
    % PLOTDETECTIONPROBABILITY displays the detection probability using a colormap.
    
    figure;
    pcolor(X, Y, PD_approx);
    shading flat;
    colormap('jet');
    colorbar;
    title('Probability of Detection (PD) - Gaussian Approximation');
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal;
end
