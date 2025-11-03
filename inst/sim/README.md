# Replication Scripts for "The R package geeVerse for Ultra-high Dimensional Heterogeneous Data Analysis with Generalized Estimating Equations"

This folder contains the R scripts necessary to replicate the simulations and data analyses presented in the accompanying manuscript.

## Manuscript

Zu, T., Green, B., & Yu, Y. (2025+). The R package geeVerse for Ultra-high Dimensional Heterogeneous Data Analysis with Generalized Estimating Equations. *Journal of Data Science, XX*(XX), 1-20.

## Requirements

The scripts require R and the following packages. Please ensure you have the latest versions installed from CRAN or the appropriate source.

* `geeVerse` (the package this paper introduces)
* `quantreg`
* `JM` (for the AIDS dataset)
* `SIS` (for the ultra-high dimensional example)

## File Descriptions

This repository contains the following R scripts:

* `replication_script_for_main.R`: This is the primary, complete script to reproduce all analyses presented in the main manuscript. This includes:
  * The real-data application: CD4 in HIV Patients (Section 5.1).
  * The real-data application: Gene Expression of Yeast Cell (Section 5.2).
  * All simulation examples (Sections 5.3.1 - 5.3.4).
  * The computational time comparisons (Section 5.4).
    This script will generate the raw results and data objects used for the tables and figures.

* `replication_script_for_main_test.R`: This is a lightweight testing version of the main replication script. It likely runs with fewer iterations, a smaller sample size, or a subset of the analyses. It is intended for quickly checking if the code and dependencies are set up correctly before running the full, time-consuming main script.


## Usage / How to Replicate

1. Install all required R packages listed under **Requirements**.
2. (Optional) Run `replication_script_for_main_test.R` to ensure your environment is configured correctly. This should run relatively quickly.
3. Set your working directory to this folder.
4. Run `replication_script_for_main.R` to perform the full replication.
   * **Note: This script is computationally intensive and may take a significant amount of time to complete.**
5. After the main script has finished, run `generate_vis_table.R` to format the results into the final tables as seen in the manuscript.