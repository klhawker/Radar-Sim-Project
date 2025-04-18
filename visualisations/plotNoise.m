function [X, Y, noise_dB] = plotNoise(params)
% PLOTTHERMALNOISE Generate and plot only the thermal noise field.
%   Inputs:
%     params — full simulation parameter struct (must contain params.radarParams)
%   Outputs:
%     X, Y      — Cartesian coordinate grids
%     noise_dB  — noise power in dB

    rp = params.radarParams;

    % Build polar grid
    cell_radius = (rp.PW * rp.C)/2;
    radii = 0:cell_radius:rp.max_radius;
    angles = 0:rp.BW_az:360;
    [Theta, R_mesh] = meshgrid(deg2rad(mod(angles,360)), radii);
    [X, Y] = pol2cart(Theta, R_mesh);

    % Thermal noise power (linear)
    k = 1.380649e-23; T = 295; B = 1/rp.PW;
    N0 = k * T * B * 10^(5/10);  % includes 5 dB noise figure

    % Generate complex AWGN noise
    noise_complex = sqrt(N0/2)*(randn(size(R_mesh)) + 1i*randn(size(R_mesh)));

    [P_cs_k, P_cs_k_dB, ranges] = computeSeaClutterK(rp, params.clutterParams, rp.max_radius, 1201, 2);

    noise_complex = noise_complex + P_cs_k;

    % Convert to power (linear → dB)
    noise_power = abs(noise_complex).^2;
    noise_dB = 10*log10(noise_power + rp.epsilon);


   % Zoom range
    zoom_x_range = [-1000, 1000];
    zoom_y_range = [-1000, 1000];
    
    % 1. Full View with zoom box
    figure;
    p1 = pcolor(X, Y, noise_dB);
    set(p1, 'EdgeColor', 'none');
    axis equal tight;
    colormap winter;
    colorbar;
    title('Full View');
    xlabel('X (m)');
    ylabel('Y (m)');
    
    % Draw zoom box
    hold on;
    rectangle('Position', ...
        [zoom_x_range(1), zoom_y_range(1), ...
         diff(zoom_x_range), diff(zoom_y_range)], ...
        'EdgeColor', 'r', 'LineWidth', 1);
    hold off;
    
    % 2. Zoomed-In View with RADAR marker
    figure;
    p2 = pcolor(X, Y, noise_dB);
    set(p2, 'EdgeColor', 'none');
    axis equal tight;
    xlim(zoom_x_range);
    ylim(zoom_y_range);
    colormap winter;
    colorbar;
    title('Zoomed-In View');
    xlabel('X (m)');
    ylabel('Y (m)');
    
    % Add RADAR marker
    hold on;
    plot(0, 0, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % white circle with black fill
    text(0, 0, '  RADAR', 'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10, 'VerticalAlignment', 'bottom');
    hold off;



end
