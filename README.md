# Empirical Power Analysis with GLMMs

This repository contains scripts and tools to simulate trial datasets and perform empirical power analysis using Generalized Linear Mixed Models (GLMMs). It includes the generation of synthetic data, model fitting, and visualization of power curves across different experimental conditions.

This code is part of the scientific article:

> **"Why accounting for data clustering matters in parasitology: simulation-based insights from hierarchical infection models"**

## Overview

The analysis involves two main steps:

1. **Data Simulation**
   Run the script `build_trials_dataframe_GLMM_for_power_analysis.R` to generate the trial dataset. This script simulates multiple trials under different experimental conditions and stores the results in a structured dataframe.

2. **Power Visualization**
   Once the data is generated, run `generate_empirical_power_plots_GLMM.R` to fit GLMM models and produce the corresponding empirical power plots.

## Interactive Execution

The script `build_trials_dataframe_GLMM_for_power_analysis.R`
supports interactive execution using `readline()`, allowing users
to modify the main simulation parameters without editing the source code.

To run the script interactively in RStudio:
- Use the **"Source"** button located in the top-right corner
  of the script editor, or
- Use the keyboard combination `Ctrl + Shift + S`.

Another way to run the script interactively is:
1. Start an R session from the terminal.
2. Execute:

```R
source("build_trials_dataframe_GLMM_for_power_analysis.R")
```

In non-interactive execution modes, `readline()` may not pause
execution to wait for user input, and default values may be used
automatically.

Default values:
- `input_bp = 0.12`
- `input_OR1 = 2.117`
- `input_OR2 = 1.492`
- `input_sigmaFarms = 1`
- `input_numberOfFarms = c(10, 20, 50, 80, 120, 160)`
- `input_numberOfPigs = c(10, 20, 50, 80)`

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
install.packages(c(
  "methods",
  "boot",
  "glmmTMB",
  "dplyr",
  "ggplot2",
  "gridExtra",
  "patchwork"
))
```
