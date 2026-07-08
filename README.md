# Beyond the Next Step: A Multi-Criteria Generative Validation Framework for Step Selection Functions
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19485572.svg)](https://doi.org/10.5281/zenodo.19485572)

This repository contains the code required to reproduce the analyses, synthetic stress tests, and empirical case studies presented in the manuscript:

> Nicosia, A. *Beyond the Next Step: A Multi-Criteria Generative Validation Framework for Step Selection Functions.* (Methods in Ecology and Evolution)

## Part of the research ecosystem

This repository is part of Aurélien Nicosia's open research and teaching ecosystem in computational statistics, scientific R software, reproducible data science and statistical education.

* Research Lab: [https://aureliennicosiaulaval.github.io/web_site/research-lab.html](https://aureliennicosiaulaval.github.io/web_site/research-lab.html)
* GitHub profile: [https://github.com/AurelienNicosiaULaval](https://github.com/AurelienNicosiaULaval)
* Related projects: [`gmov`](https://github.com/AurelienNicosiaULaval/gmov), [`evalue-HMM`](https://github.com/AurelienNicosiaULaval/evalue-HMM)

## Repository Structure
This repository contains the canonical code associated with the paper, isolated into four categories:

- **`Functions/`**: Contains `diagnose_issf.R`. This script defines the four-pillar validation framework (Wasserstein, MSD, sinuosity, barrier crossing) applied to standard iSSA.
- **`Simulations/`**: Contains the R scripts used to run the synthetic stress-test scenarios (Scenarios 1 to 6) evaluating structural failure modes.
- **`Data_Analysis/`**: Contains `red_deer_empirical_analysis_final.R`, the empirical application of the framework to the red deer dataset.
- **`Supplementary/`**: Contains `barrier_crossing_sensitivity.R` and associated sensitivity analyses from Supplementary Information S3.
- **`docs/`**: Includes environment/documentation files like `sessionInfo.txt`.

*(Note: Red deer data utilized in the empirical application are acquired programmatically via the `amt` R package. There is no redundant raw `.csv` or `.rds` data file included in this repository.)*

## Reproducibility overview

This repository provides reproducible research material for the generative validation framework applied to SSF/iSSA analyses. It is not an R package; the scripts are intended to be run from the repository root.

- Shared functions are stored in `Functions/`.
- Synthetic stress-test simulations are stored in `Simulations/`.
- The empirical red deer analysis is stored in `Data_Analysis/`.
- Supplementary analyses are stored in `Supplementary/`.
- Generated outputs are organized in `output/` or in script-specific result and figure folders created during execution.

The main dependencies are listed below. The file `docs/sessionInfo.txt` records one reference R environment, but it does not replace a lockfile such as `renv.lock`. Some analyses may be time-consuming or depend on the local R and spatial software environment. No single-command workflow for full reproduction is provided.

## Dependencies
The code requires **R** and standard movement ecology and spatial packages, notably:

- `amt` (for iSSA fitting and the built-in red deer dataset)
- `transport` (for Wasserstein distance calculation via the network simplex algorithm)
- `terra`, `sf` (for spatial operations)
- `ggplot2`, `patchwork` (for visualization)
