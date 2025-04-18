%% shadowMaskPlot_polar_dynamic.m
% This script generates and plots turbine shadow masks on a PD map
% using a polar PPI display with pcolor.
%
% For each turbine:
% - The shadow is now cast in the direction from the radar (at [0,0]) to 
%   the turbine’s position.
% - The half-angle for shadow overlap is computed as:
%       halfAngle = atan(pole_radius / R)
%   where R is the distance from the radar to the turbine.
%
% - Turbines that fall within the shadow cone of an upwind turbine (in the 
%   dynamic shadow direction) do not cast their own shadow.
%
% The PD map is formatted so that any PD value greater than or equal
% to 0.7 is green and values below 0.7 are black.
% Also, the base of each shadow (where it originates at the turbine)
% is drawn with a width equal to the turbine pole’s diameter.

%% Load simulation parameters and turbine data
params = simulationConfig();
[turbinePositions, turbineSegments] = generateTurbineGridPositions(params.gridParams, params.turbineParams);

%% Run the radar simulation to obtain the PD map and polar grid (X,Y)
[radarData, radarData_dB, PD_approx, object_Pr_dB, X, Y, seaClutter, turbineReturns, radarNoise, noise_and_clutter] = radarSimulation(params, turbineSegments);

noise_and_clutter_dB = 10*log10(noise_and_clutter);
% norm_noise_and_seaclutter = noise_and_clutter_dB/-137.8804;

%% Create a binary PD map based on the threshold
PD_threshold = 0.0 ;
PD_binary = double(PD_approx >= PD_threshold);

%% Create a custom binary colormap:
% 0 (below threshold) -> black, 1 (above threshold) -> green.
binaryColormap = [0 0 0; 0 1 0];

%% (Removed fixed shadow direction)
% The fixed shadow direction has been removed. Now, for each turbine the shadow 
% will extend in the direction of the vector from the radar (0,0) to the turbine.
% (old code:)
% fixedShadowAngle = pi/2 - deg2rad(30);

%% Compute the half-angle for the shadow cone for each turbine (for shadow-checking)
numTurbines = size(turbinePositions, 1);
halfAngles = zeros(numTurbines,1);
for i = 1:numTurbines
    R_i = norm(turbinePositions(i,:));
    halfAngles(i) = atan(params.turbineParams.pole_radius / R_i);
end

%% Determine which turbines cast a shadow
% A turbine will not cast a shadow if it lies within the cone (shadow) 
% of another turbine that is upwind relative to that turbine’s own radial direction.
%
% For each turbine j we now compute its dynamic shadow direction from the radar:
%   shadowAngle_j = atan2(turbinePositions(j,2), turbinePositions(j,1));
%
% Then turbine i is considered to be in the shadow of j if:
%   1. Its projection along the direction from the radar to turbine j is greater 
%      than that of turbine j (i.e. turbine i is further out).
%   2. The angular difference between the relative position vector (i from j) and 
%      j’s shadow direction is less than the half-angle for turbine j.
castShadow = true(numTurbines, 1);
for i = 1:numTurbines
    pos_i = turbinePositions(i,:);
    for j = 1:numTurbines
        if i == j
            continue;
        end
        pos_j = turbinePositions(j,:);
        % Compute dynamic shadow direction for turbine j:
        shadowAngle_j = atan2(pos_j(2), pos_j(1));
        d_j = [cos(shadowAngle_j), sin(shadowAngle_j)];
        % Only consider turbine j if i is farther in turbine j’s shadow direction.
        if dot(pos_i, d_j) <= dot(pos_j, d_j)
            continue;
        end
        % Compute the relative vector from turbine j to i.
        rel_vec = pos_i - pos_j;
        theta_rel = atan2(rel_vec(2), rel_vec(1));
        % Compute the minimal angular difference relative to turbine j's shadow direction.
        angle_diff = abs(mod(theta_rel - shadowAngle_j + pi, 2*pi) - pi);
        % If this difference is less than the half-angle for turbine j, turbine i is in j's shadow.
        if angle_diff <= halfAngles(j)
            castShadow(i) = false;
            break;
        end
    end
end

%% Define a very large extension length for the shadow (indefinite extension)
L = 1e5;

xZoom = [-800, 2875];
yZoom = [1225, 4500];


% Create figure with subplots
figure;

% --- Subplot 1: Full PD Map ---
pcolor(X, Y, PD_binary);
shading flat;
colormap(binaryColormap);
caxis([0 1]);
xlabel('X (m)');
ylabel('Y (m)');
axis equal;
xlim([-35000 35000]);
ylim([-500 60000])
hold on;
axis off;

% Overlay concentric circles and axis lines
ax = axis;
x_min = ax(1);
x_max = ax(2);
y_min = ax(3);
y_max = ax(4);

radii = 1000:1000:70000;
theta_circle = linspace(0, 2*pi, 360);
for r = radii
    x_circle = r * cos(theta_circle);
    y_circle = r * sin(theta_circle);
    plot(x_circle, y_circle, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

plot([x_min, x_max], [0, 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot([0, 0], [y_min, y_max], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
rectangle('Position', [xZoom(1), yZoom(1), diff(xZoom), diff(yZoom)], ...
          'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', ':');

% Plot turbines and radar
plot(turbinePositions(:,1), turbinePositions(:,2), 'o', ...
     'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 0.5);
plot(0, 0, 'mo', 'MarkerSize', 6, 'LineWidth', 2, 'MarkerFaceColor', 'm');

% Draw shadows (with widening and indefinite length) using dynamic shadow directions:
for i = 1:numTurbines
    if ~castShadow(i)
        continue;
    end
    base = turbinePositions(i,:);
    pole_radius = params.turbineParams.pole_radius;
    
    % Compute the dynamic shadow direction for turbine i:
    shadowAngle = atan2(base(2), base(1));
    
    % The cone boundaries are computed relative to this dynamic shadow direction:
    angle_left  = shadowAngle + halfAngles(i);
    angle_right = shadowAngle - halfAngles(i);
    
    % Compute far end of the shadow along the two boundaries:
    far_left  = base + L * [cos(angle_left), sin(angle_left)];
    far_right = base + L * [cos(angle_right), sin(angle_right)];
    
    % Define the base of the shadow (turbine's physical extent) rotated about shadowAngle:
    base_left  = base + pole_radius * [cos(shadowAngle + pi/2), sin(shadowAngle + pi/2)];
    base_right = base + pole_radius * [cos(shadowAngle - pi/2), sin(shadowAngle - pi/2)];
    
    coneX = [base_left(1), far_left(1), far_right(1), base_right(1)];
    coneY = [base_left(2), far_left(2), far_right(2), base_right(2)];
    patch(coneX, coneY, 'k', 'EdgeColor', 'k', 'FaceAlpha', 1);
end

% --- Subplot 2: Zoomed-In View ---
figure;
pcolor(X, Y, PD_binary);
colorbar off;
shading flat;
colormap(binaryColormap);
caxis([0 1]);
xlabel('X (m)');
ylabel('Y (m)');
axis equal
xlim(xZoom);
ylim(yZoom)
axis off;
hold on;

ax = axis;
x_min = ax(1);
x_max = ax(2);
y_min = ax(3);
y_max = ax(4);

radii = 1000:1000:max([abs(x_min), abs(x_max), abs(y_min), abs(y_max)]);
theta_circle = linspace(0, 2*pi, 360);
for r = radii
    x_circle = r * cos(theta_circle);
    y_circle = r * sin(theta_circle);
    plot(x_circle, y_circle, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

plot([x_min, x_max], [0, 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot([0, 0], [y_min, y_max], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot(turbinePositions(:,1), turbinePositions(:,2), 'o', ...
     'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 0.5);
plot(0, 0, 'mo', 'MarkerSize', 2, 'LineWidth', 2, 'MarkerFaceColor', 'm');

for i = 1:numTurbines
    if ~castShadow(i)
        continue;
    end
    base = turbinePositions(i,:);
    pole_radius = params.turbineParams.pole_radius;
    
    shadowAngle = atan2(base(2), base(1));
    
    angle_left  = shadowAngle + halfAngles(i);
    angle_right = shadowAngle - halfAngles(i);
    
    far_left  = base + L * [cos(angle_left), sin(angle_left)];
    far_right = base + L * [cos(angle_right), sin(angle_right)];
    
    base_left  = base + pole_radius * [cos(shadowAngle + pi/2), sin(shadowAngle + pi/2)];
    base_right = base + pole_radius * [cos(shadowAngle - pi/2), sin(shadowAngle - pi/2)];
    
    coneX = [base_left(1), far_left(1), far_right(1), base_right(1)];
    coneY = [base_left(2), far_left(2), far_right(2), base_right(2)];
    patch(coneX, coneY, 'k', 'EdgeColor', 'k');
end

% --- Create Dummy Handles for Legend ---
hDetectable = patch(NaN, NaN, [0 1 0], 'EdgeColor', 'none');
hUndetectable = patch(NaN, NaN, [0 0 0], 'EdgeColor', 'none');
hTurbine = plot(NaN, NaN, 'o', 'MarkerSize', 6, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 0.5);
hRadar = plot(NaN, NaN, 'mo', 'MarkerSize', 6, 'LineWidth', 2, 'MarkerFaceColor', 'm');

legend([hDetectable, hUndetectable, hTurbine, hRadar], ...
       {'Detectable', 'Undetectable', 'Turbine', 'Radar'}, ...
       'Location', 'north', 'Orientation', 'horizontal');
