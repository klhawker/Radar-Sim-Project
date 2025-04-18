function [P_cs_k, P_cs_k_dB, ranges] = computeSeaClutterK(radarParams, clutterParams, max_radius, num_angles, nu)
    % computeSeaClutterK computes sea clutter power incorporating a K-distribution model.
    %
    % Inputs:
    %   radarParams   - structure with radar parameters.
    %   clutterParams - structure with clutter parameters.
    %   max_radius    - maximum range to compute clutter.
    %   num_angles    - number of azimuth angles.
    %   nu            - shape parameter for the Gamma-distributed texture.
    %
    % Outputs:
    %   P_cs_k       - Clutter power (linear scale) after applying the K-distribution.
    %   P_cs_k_dB    - Clutter power in dB.
    %   ranges       - Vector of range midpoints.
    
    % First, compute the deterministic clutter using your existing model:
    [P_cs_rotated, ~, ranges] = calculateSeaClutter(radarParams, clutterParams, max_radius, num_angles);
    % P_cs_rotated is of size (nRanges x num_angles)
    P_cs_rotated = ones(size(P_cs_rotated));

    % Generate the random texture component (Gamma distributed) for each cell.
    texture = gamrnd(nu, 1/nu, size(P_cs_rotated));
    
    % Generate the speckle component (Exponential, mean = 1) for each cell.
    speckle = exprnd(1, size(P_cs_rotated));
    
    % Multiply to get the K-distributed random factor for each cell.
    kFactor = texture .* speckle;
    
    % Apply the multiplicative factor to the deterministic clutter.
    P_cs_k = P_cs_rotated .* kFactor;
    
    % Convert to dB
    epsilon = radarParams.epsilon;
    P_cs_k_dB = 10 * log10(max(P_cs_k, epsilon));
end
