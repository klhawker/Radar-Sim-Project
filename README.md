# Radar Simulation Project for Offshore Wind Farms

[Repository → https://github.com/klhawker/Radar-Sim-Project](https://github.com/klhawker/Radar-Sim-Project)

## Introduction

This MATLAB‑based project simulates the radar detection environment around offshore wind farms. It models the key contributors to a radar return—including additive white Gaussian noise (AWGN), sea clutter (K‑distribution or Gamma‑Speckle models), coherent/incoherent summation of turbine pole returns, turbine shadowing, and multipath “ghost” reflections. The goal is to provide researchers and engineers with a modular framework to explore detection probability, false alarms, and radar performance impacts of turbine layouts and environmental conditions.

## Installation

1. **Prerequisites**  
   - MATLAB R2020b or later  
   - Signal Processing Toolbox  

2. **Clone the repository**  
   ```bash
   git clone https://github.com/your-username/radar-simulation-project.git
   cd radar-simulation-project
   ```
3. **Add to MATLAB path**
   ```bash
   addpath(genpath(pwd));
   savepath;
   ```
4. **Configure parameters**
   Open simulationConfig.m and set your radar parameters (pulse width, bandwidth, antenna gains, etc.), turbine layout (rows, columns, spacing, wind direction), clutter parameters, shadow and ghost parameters, and any additional object (e.g., a boat).

## Usage

### 1. Basic run

In the MATLAB Command Window:

```matlab
mainSimulation
```
This will:
- Generate turbine grid positions and pole segments
- Compute AWGN noise and sea clutter
- Calculate coherent turbine returns and map onto a PPI grid
- Compute dynamic shadow cones
- Plot detection probability, clutter, shadows, and (optionally) ghost multipath heat maps

### 2. Simple example
```matlab
params = simulationConfig();                  % load defaults
[pos, segs] = generateTurbineGridPositions(params.gridParams, params.turbineParams);
[rd, rd_dB, PD, obj_dB, X, Y, clutter, tReturns, noise, nc] = radarSimulation(params, segs);
figure;
pcolor(X, Y, PD);
shading flat;
colormap(jet);
colorbar;
title('Detection Probability');
```
### 3. Custom demos
Uncomment any of the plotting lines in `mainSimulation.m` to view:
- Sea clutter curves (`plotKSeaClutter`, `plotSeaClutterPlots`)
- Radar heat maps, PD thresholds, object returns
- Ghost heat maps (`complete_ghosts`, `plotRadarHeatMapWithGhosts`)

## Known Issues & Future Improvements
- UI: develop a simple UI for parameter sweeps and live visualization.
- Real Terrain: integrate real‑world bathymetry and wind‑farm layouts.
- Modular Packaging: split into MATLAB classes or functions for easier unit testing.
