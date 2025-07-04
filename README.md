# Plant-Pollinator Network Coexistence Analysis
*Accompanies the manuscript "Equalizing Effect of Pollinator Adaptive Foraging on Plant Network Coexistence" by Valdovinos, Kaur, and Marsland III*

## Overview
This code analyzes plant-pollinator network dynamics, focusing on how adaptive foraging behavior affects species coexistence in ecological networks. It processes multiple network configurations and evaluates network properties, species persistence, and coexistence criteria.

## Data Requirements

### Input Files
1. `network_properties_1200m.csv`: Contains network metrics including:
   - Basic network properties (P, A, S)
   - Connectance (C)
   - Nestedness metrics (NODFst)
   - Degree distributions
   - Other network-level indices

2. `1200m.mat`: Contains the adjacency matrices for the networks

## Core Components

### 1. Network Classification
- Networks are classified based on:
  - Nestedness (NODFst > 0.2)
  - Size categories (S > 10, S > 60, S > 150)
  - Structural properties

### 2. Simulation Parameters
- `dataset`: Network size (default: 1200)
- `muAP`: Mutualistic interaction strength (default: 2)
- `frG`: Foraging behavior flag (0 = fixed, 1 = adaptive)

### 3. Key Analyses
For each network, the code:
1. Simulates dynamics under fixed and adaptive foraging
2. Calculates coexistence metrics
3. Evaluates species persistence
4. Computes niche overlap (rho) between species

### 4. Output Metrics
For each plant species:
- Growth rates (k_i)
- Resource allocation (rho_i)
- Population dynamics (P, R)
- Persistence predictions
- Coexistence criteria

## Output Files

### Network-Level Results
`Coex_mean_muAP{muAP}_{dataset}m.csv`:
- Network properties
- Mean niche overlap
- Fraction of persistent species
- Criterion check results

### Species-Level Results
`Coex_AllPsp_muAP{muAP}_{dataset}m.csv`:
- Individual species metrics
- Network context
- Persistence outcomes
- Coexistence criteria

## Implementation Details

### Main Process Flow
1. Load network data and properties
2. Classify networks
3. For each network:
   ```matlab
   for i = 1:dataset
       % Process both fixed and adaptive foraging
       for frG = [0,1]
           % Run dynamics
           % Calculate metrics
           % Store results
       end
   end
   ```

### Key Functions
- `create_metadata`: Sets up simulation parameters
- `IntegrateValdovinos2013_dataset`: Runs dynamic simulations
- `calCoexist`: Calculates coexistence metrics

## Dependencies
- MATLAB (with Statistics and Machine Learning Toolbox)
- Required data files in working directory

## Usage
1. Ensure input files are present
2. Set desired parameters (dataset size, muAP)
3. Run the main script
4. Results are saved automatically as CSV files

## Notes
- Processing large networks may be computationally intensive
- Results combine both structural and dynamical analyses
- Supports comparison between fixed and adaptive foraging strategies

## References
Valdovinos, F. S., Kaur, T., & Marsland III, R. "Equalizing Effect of Pollinator Adaptive Foraging on Plant Network Coexistence"
