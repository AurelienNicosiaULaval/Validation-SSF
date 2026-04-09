# Supplementary Information S3

This folder contains the full code for the sensitivity analysis discussed in Supplementary Information S3 (barrier-crossing diagnostic power).

## File

- `barrier_crossing_sensitivity.R`
  - Reproduces the Monte Carlo grid experiment used to quantify power across trajectory length (`N`) and diffusivity (`sigma`) as described in SI3.
  - Generates figures used in the manuscript (`barrier_sensitivity_N.png` and `barrier_sensitivity_sigma.png`).
  - Writes summary CSV files in `results/`:
    - `barrier_sensitivity_summary.csv`
    - `barrier_sensitivity_details.csv`

Run from repository root with:

```r
Rscript MEE_Code/Supplementary/S3/barrier_crossing_sensitivity.R
```
