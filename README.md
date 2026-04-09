# Beyond the Next Step: A Multi-Criteria Generative Validation Framework for Step Selection Functions
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19485572.svg)](https://doi.org/10.5281/zenodo.19485572)

This repository contains the code required to reproduce the analyses, synthetic stress tests, and empirical case studies presented in the manuscript:

> Nicosia, A. *Beyond the Next Step: A Multi-Criteria Generative Validation Framework for Step Selection Functions.* (Methods in Ecology and Evolution)

## Repository Structure
This repository contains the canonical code associated with the paper, isolated into four categories:

- **`Functions/`**: Contains `diagnose_issf.R`. This script defines the four-pillar validation framework (Wasserstein, MSD, sinuosity, barrier crossing) applied to standard iSSA.
- **`Simulations/`**: Contains the R scripts used to run the synthetic stress-test scenarios (Scenarios 1 to 6) evaluating structural failure modes.
- **`Data_Analysis/`**: Contains `red_deer_empirical_analysis_final.R`, the empirical application of the framework to the red deer dataset.
- **`Supplementary/`**: Contains `barrier_crossing_sensitivity.R` and associated sensitivity analyses from Supplementary Information S3.
- **`docs/`**: Includes environment/documentation files like `sessionInfo.txt`.

*(Note: Red deer data utilized in the empirical application are acquired programmatically via the `amt` R package. There is no redundant raw `.csv` or `.rds` data file included in this repository.)*

## Dependencies
The code requires **R** and standard movement ecology and spatial packages, notably:

- `amt` (for iSSA fitting and the built-in red deer dataset)
- `transport` (for Wasserstein distance calculation via the network simplex algorithm)
- `terra`, `sf` (for spatial operations)
- `ggplot2`, `patchwork` (for visualization)
