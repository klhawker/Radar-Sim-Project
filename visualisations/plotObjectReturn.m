function plotObjectReturn(X, Y, object_Pr_dB)
    % PLOTOBJECTRETURN visualizes the power returned from an additional object.
    
    figure;
    pcolor(X, Y, object_Pr_dB);
    shading flat;
    colormap('jet');
    colorbar;
    title('Power Returned From the Object (dB)');
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal;
end
