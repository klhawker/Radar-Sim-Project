function [multipathHeatMap, ghostLocations, debugInfo] = computeTwoBounceGhost(params, turbinePositions, X, Y, turbineSegments, debugFlag)
% computeTwoBounceGhost computes a ghost return heat map using a modified
% two‐bounce approach that separates the process into three stages.
%
% For each turbine pair (i,j):
%
% Stage 1 – Radar to turbine i segments:
%   [Compute each segment's surface point and the corresponding received
%    power contribution (P_ref).]
%
% Stage 2 – Turbine i midpoint to turbine j segments:
%   [Compute each segment's contribution (P_ref_prime) from turbine i to j.]
%
% Stage 3 – Turbine j midpoint to radar:
%   [Compute the final on‐boresight received power (P_received_max).]
%
% Then, candidate ghost returns are computed using candidate TX/RX gains.
%
% Inputs:
%   params           : Structure with radar and turbine parameters.
%   turbinePositions : Nx2 array of [x, y] positions for each turbine.
%   X, Y             : Matrices of grid coordinates.
%   turbineSegments  : Matrix of turbine pole segments (each turbine has N_pole rows).
%   debugFlag        : (Optional) If true, debugInfo is populated.
%
% Outputs:
%   multipathHeatMap : 2D matrix with summed ghost return power.
%   ghostLocations   : Mx2 array of ghost target locations.
%   debugInfo        : Structure with intermediate variables (if debugFlag true).

    if nargin < 6
        debugFlag = false;
    end
    debugInfo = struct();
    if debugFlag
        % Create a structure array to store data for each turbine pair.
        debugInfo.turbinePairs = [];
        pairCounter = 0;
    end

    %% Unpack parameters
    rp           = params.radarParams;
    radarPos2D   = rp.radar_pos;   % [x, y]
    radarHeight  = rp.h_r;
    radarPos3D   = [radarPos2D, radarHeight];
    P_t          = rp.P_peak;
    baseG_t      = rp.G_t;  % Full (on-boresight) TX gain.
    baseG_r      = rp.G_r;
    lambda       = rp.lamb;
    
    % Turbine parameters.
    r_pole  = params.turbineParams.pole_radius;
    h_pole  = params.turbineParams.pole_height;
    N_pole  = params.turbineParams.N_pole;
    
    N_turbines = size(turbinePositions, 1);
    
    % Define turbine midpoints (for example, at half the pole height).
    midPoints = zeros(N_turbines, 3);
    for i = 1:N_turbines
        midPoints(i,:) = [turbinePositions(i,1), turbinePositions(i,2), h_pole/2];
    end
    
    multipathHeatMap = zeros(size(X));
    ghostLocations   = [];
    
    % Define candidate boresight offsets (in degrees) for the ±5° window.
    offsets_deg = -5:rp.BW_az:5;
    offsets_rad = deg2rad(offsets_deg);  % convert to radians
    
    %% Loop over all turbine pairs (i,j)
    for i = 1:N_turbines
        % Retrieve segments and midpoint for turbine i.
        idx_i = (i-1)*N_pole + 1;
        idx_i_end = i*N_pole;
        segs_i = turbineSegments(idx_i:idx_i_end, :);
        mid_i = midPoints(i, :);
        
        % Compute the central (TX) direction for turbine i: from radar to turbine i.
        vecRadarToTi = mid_i(1:2) - radarPos2D;
        central_tx = atan2(vecRadarToTi(2), vecRadarToTi(1));
        
        % (Pre-loop stage: if needed, one could accumulate a global total here.)
        
        % Now, for each turbine pair (i,j), recompute stage 1 with the proper scattering angle.
        for j = 1:N_turbines
            if i == j, continue; end

            
            
            % Retrieve segments and midpoint for turbine j.
            idx_j = (j-1)*N_pole + 1;
            idx_j_end = j*N_pole;
            segs_j = turbineSegments(idx_j:idx_j_end, :);
            mid_j = midPoints(j, :);
            
            % Define the turbine-to-turbine vector (from turbine i to turbine j)
            v_ij = mid_j - mid_i;
            
            %% Stage 1 (Revisited): Radar → Turbine i segments using central TX
            total_P_ref = 0;
            sum_R1 = 0;
            for s_i = 1:size(segs_i, 1)
                seg_i = segs_i(s_i, :);
                seg_i_surface = seg_i;
                direction_i = radarPos2D - turbinePositions(i,:);
                direction_i = direction_i / norm(direction_i);
                seg_i_surface(1:2) = seg_i(1:2) + r_pole * direction_i;
                
                R1 = norm(seg_i_surface - radarPos3D);
                sum_R1 = sum_R1 + R1;
                
                v1 = seg_i_surface - radarPos3D;
                phi1 = atan2(v1(3), norm(v1(1:2)));
                theta_bi = angleBetween(v1, v_ij);
                sigma_bi = computeBistaticRCS_PO(2*pi/lambda, r_pole, h_pole, theta_bi, phi1, 0);
                
                % Use full TX gain (on-boresight)
                P_ref = (P_t * baseG_t * sigma_bi) / (4*pi * R1^2);
                total_P_ref = total_P_ref + P_ref;
            end
            avg_R1 = sum_R1 / size(segs_i,1);
            % Turbine i now acts as a source with power total_P_ref.
            
            %% Stage 2: Turbine i midpoint → Turbine j segments
            total_P_ref_prime = 0;
            sum_R2 = 0;
            for s_j = 1:size(segs_j, 1)
                seg_j = segs_j(s_j, :);
                seg_j_surface = seg_j;
                direction_j = (mid_i(1:2) - turbinePositions(j,:));
                direction_j = direction_j / norm(direction_j);
                seg_j_surface(1:2) = seg_j(1:2) + r_pole * direction_j;
                
                R2 = norm(seg_j_surface - mid_i);
                sum_R2 = sum_R2 + R2;
                
                v2 = seg_j_surface - mid_i;
                phi2 = atan2(v2(3), norm(v2(1:2)));
                v_rj = radarPos3D - mid_j;
                theta_bj = angleBetween(v_ij, v_rj);
                sigma_bj = computeBistaticRCS_PO(2*pi/lambda, r_pole, h_pole, theta_bj, 0, phi2);
                
                P_ref_prime = (total_P_ref * sigma_bj) / (4*pi * R2^2);
                total_P_ref_prime = total_P_ref_prime + P_ref_prime;
            end
            avg_R2 = sum_R2 / size(segs_j,1);
            % Turbine j acts as a source with power total_P_ref_prime.
            
            %% Stage 3: Turbine j midpoint → Radar (central RX)
            R3 = norm(mid_j - radarPos3D);
            P_received_max = (total_P_ref_prime * lambda^2 * baseG_r) / ((4*pi * R3)^2);
            
            % Define the average ghost distance.
            ghost_distance_avg = (avg_R1 + avg_R2 + R3) / 2;
            
            % Compute the central RX angle (from radar to turbine j)
            central_rx = atan2(mid_j(2) - radarPos2D(2), mid_j(1) - radarPos2D(1));
            
            % Prepare to store candidate information for this turbine pair.
            candidateDebug = [];
            
            % Iterate candidate boresight offsets relative to the central TX direction.
            for k = 1:length(offsets_rad)
                offset = offsets_rad(k);
                candidate_boresight = central_tx + offset;
                
                % Compute candidate TX and RX gains.
                candidate_tx = computeAntennaPattern(rp, lambda, abs(offset), 'TX');
                candidate_rx = computeAntennaPattern(rp, lambda, abs(offset), 'RX');
                candidateGain = candidate_tx * candidate_rx;
                
                % Compute candidate ghost return power.
                candidate_P_received = P_received_max * candidateGain;
                
                % Compute candidate ghost location.
                ghost_loc_candidate = radarPos2D + ghost_distance_avg * [cos(candidate_boresight), sin(candidate_boresight)];
                
                % Map candidate ghost power onto the heat map.
                dMat = sqrt((X - ghost_loc_candidate(1)).^2 + (Y - ghost_loc_candidate(2)).^2);
                [~, idx] = min(dMat(:));
                multipathHeatMap(idx) = multipathHeatMap(idx) + candidate_P_received;
                
                % Record candidate ghost location.
                ghostLocations = [ghostLocations; ghost_loc_candidate];
                
                % If debugging, store candidate details.
                if debugFlag
                    candidateDebug(end+1).offset_rad = offset;
                    candidateDebug(end).offset_deg = rad2deg(offset);
                    candidateDebug(end).candidate_boresight = candidate_boresight;
                    candidateDebug(end).candidate_tx = candidate_tx;
                    candidateDebug(end).candidate_rx = candidate_rx;
                    candidateDebug(end).candidateGain = candidateGain;
                    candidateDebug(end).candidate_P_received = 10*log10(candidate_P_received);
                    candidateDebug(end).ghost_loc_candidate = ghost_loc_candidate;
                end
            end
            
            % If debugging, store all intermediate data for this turbine pair.
            if debugFlag
                pairCounter = pairCounter + 1;
                debugInfo.turbinePairs(pairCounter).turbine_i = i;
                debugInfo.turbinePairs(pairCounter).turbine_j = j;
                debugInfo.turbinePairs(pairCounter).total_P_ref = 10*log10(total_P_ref);
                debugInfo.turbinePairs(pairCounter).avg_R1 = avg_R1;
                debugInfo.turbinePairs(pairCounter).total_P_ref_prime = 10*log10(total_P_ref_prime);
                debugInfo.turbinePairs(pairCounter).avg_R2 = avg_R2;
                debugInfo.turbinePairs(pairCounter).R3 = R3;
                debugInfo.turbinePairs(pairCounter).P_received_max = 10*log10(P_received_max);
                debugInfo.turbinePairs(pairCounter).ghost_distance_avg = ghost_distance_avg;
                debugInfo.turbinePairs(pairCounter).central_tx = central_tx;
                debugInfo.turbinePairs(pairCounter).central_rx = central_rx;
                debugInfo.turbinePairs(pairCounter).candidates = candidateDebug;
            end
        end
    end
    
    % Assign key results to the base workspace.
    assignin("base", "multipathHeatMap", multipathHeatMap);
    assignin("base", "ghostLocations", ghostLocations);
    [~,~,non_zeros_ghost_returns] = find(multipathHeatMap);
    non_zeros_ghost_returns = 10*log10(non_zeros_ghost_returns);
    assignin("base", "non_zeros_ghost_returns", non_zeros_ghost_returns);
    assignin("base", "debugInfo", debugInfo);
    
    fprintf("Max Ghosts: %f\n", max(non_zeros_ghost_returns));

    V = max(multipathHeatMap, [], 2);
    V(V == 0) = [];
    m = mean(V);
    fprintf("Mean Ghosts: %f\n", 10*log10(m));
end

%% ------------------------------------------------------------------------
% Helper function: angleBetween (3D)
function angleRad = angleBetween(a, b)
    dotVal = dot(a, b);
    norms  = norm(a) * norm(b);
    cosTheta = dotVal/(norms + eps);
    cosTheta = max(min(cosTheta, 1), -1);
    angleRad = acos(cosTheta);
end

% Helper function: angleBetween2D
function angleRad = angleBetween2D(a2d, b2d)
    a2d = a2d/(norm(a2d)+eps);
    b2d = b2d/(norm(b2d)+eps);
    c   = dot(a2d, b2d);
    c   = max(min(c, 1), -1);
    angleRad = acos(c);
end

% Helper function: computeAntennaPattern
% Returns a gain factor (between 0 and 1) for a given angular deviation.
% For TX mode, applies a cutoff for abs(angle) > 5°. For RX mode, no cutoff.
function patGain = computeAntennaPattern(rp, lambda, angleRad, mode)
    if nargin < 4
        mode = 'TX'; % default to TX mode
    end

    if strcmpi(mode, 'TX')
        if abs(angleRad) <= deg2rad(5)
            u_az = (rp.D_az * pi * sin(angleRad)) / lambda;
            g_az_val = sinc(u_az/pi);
            u_el = (rp.D_el * pi * sin(angleRad)) / lambda;
            g_el_val = sinc(u_el/pi);
            patGain = (abs(g_az_val)^2) * (abs(g_el_val)^2);
        else
            patGain = 0;
        end
    else  % RX mode: no cutoff applied.
        u_az = (rp.D_az * pi * sin(angleRad)) / lambda;
        g_az_val = sinc(u_az/pi);
        u_el = (rp.D_el * pi * sin(angleRad)) / lambda;
        g_el_val = sinc(u_el/pi);
        patGain = (abs(g_az_val)^2) * (abs(g_el_val)^2);
    end
end
