function radarNoise = computeAWGN(PW, R_grid_mesh, NF_dB)
%GENERATERADARNOISE Generate complex AWGN noise for radar simulation.
%
%   radarNoise = generateRadarNoise(PW, R_grid_mesh, NF_dB) generates complex
%   additive white Gaussian noise (AWGN) for a given radar simulation.
%
%   Inputs:
%       PW          - Pulse width in seconds.
%       R_grid_mesh - Matrix (or array) whose dimensions determine the noise 
%                     matrix size.
%       NF_dB       - (Optional) Noise figure in dB. Default value is 5 dB.
%
%   Output:
%       radarNoise  - Complex AWGN noise matrix.
%
%   The function calculates the thermal noise power using the AWGN model:
%       N0 = k * T * B * NF
%   where:
%       k   = Boltzmann's constant (1.380649e-23 J/K),
%       T   = Temperature (295 K),
%       B   = Receiver bandwidth (approximated as 1/PW),
%       NF  = Noise figure (linear scale).
%

    % Boltzmann's constant (J/K)
    k = 1.380649e-23;
    
    % Temperature in Kelvin
    T = 295;
    
    % Use pulse width to calculate the receiver bandwidth
    tau = PW;
    B = 1 / tau;
    
    % Thermal noise power (in Watts)
    N0 = k * T * B;
    
    % Set default noise figure if not provided
    if nargin < 3
        NF_dB = 5;
    end
    
    % Convert noise figure from dB to linear scale and update noise power
    NF = 10^(NF_dB/10);
    N0 = N0 * NF;
    
    % Generate complex AWGN noise matrix
    radarNoise = sqrt(N0/2) * (randn(size(R_grid_mesh)) + 1i * randn(size(R_grid_mesh)));
end
