% validateAWGN.m
% This script validates the AWGN noise model based on a constant noise power of -137.8803 dB.

% Clear workspace and command window
clear; clc;

%% Define the constant noise power in dB and convert to linear scale
noise_dB_const = -137.8803;       % Given constant noise power (dB)
N0 = 10^(noise_dB_const/10);        % Linear noise power

%% Generate complex AWGN noise samples
% Since the noise is complex, each dimension (real and imaginary) has variance N0/2.
numSamples = 1e5;                 % Number of noise samples to generate
noiseSamples = sqrt(N0/2) * (randn(numSamples,1) + 1i * randn(numSamples,1));

% Extract the in-phase (real) component
noiseReal = real(noiseSamples);

%% Plot the histogram of the real component
figure;
histogram(noiseReal, 'Normalization', 'pdf');
hold on;

%% Overlay the theoretical Gaussian PDF
% The real component is distributed as N(0, N0/2).
sigma = sqrt(N0/2);  % Standard deviation of the real part
x = linspace(min(noiseReal), max(noiseReal), 1000);
pdf_theoretical = (1/(sigma * sqrt(2*pi))) * exp(-x.^2 / (2*sigma^2));

plot(x, pdf_theoretical, 'r', 'LineWidth', 2);
xlabel('In-phase Noise Component');
ylabel('Probability Density');
title('Histogram of AWGN Noise (Real Component) with Theoretical PDF');
legend('Simulated Histogram', 'Theoretical N(0, N_0/2) PDF');
hold off;

%% Validate the noise level by computing the measured average noise power (in dB)
measuredNoise_dB = 10*log10(mean(abs(noiseSamples).^2));
fprintf('Measured average noise power: %.4f dB\n', measuredNoise_dB);
