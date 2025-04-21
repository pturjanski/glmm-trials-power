# Empirical Power Analysis with GLMMs

This repository contains scripts and tools to simulate trial datasets and perform empirical power analysis using Generalized Linear Mixed Models (GLMMs). It includes the generation of synthetic data, model fitting, and visualization of power curves across different experimental conditions.

This code is part of the scientific article:

> **"Why ignoring data structure can lead to flawed statistical analysis: an application in epidemiological research"**

## Overview

The analysis involves two main steps:

1. **Data Simulation**  
   Run the script `build_trials_dataframe_GLMM_for_power_analysis.R` to generate the trial dataset. This script simulates multiple trials under different experimental conditions and stores the results in a structured dataframe.

2. **Power Visualization**  
   Once the data is generated, run `generate_empirical_power_plots_GLMM.R` to fit GLMM models and produce the corresponding empirical power plots.

## Requirements

This project is written in **R**. To run the scripts, you will need the following R packages:
- `methods`
- `boot`
- `glmmTMB`
- `dplyr`
- `ggplot2`
- `gridExtra`
- `patchwork`

You can install missing packages using:

```R
install.packages(c("methods", "boot", "glmmTMB", "dplyr", "ggplot2", "gridExtra", "patchwork"))
