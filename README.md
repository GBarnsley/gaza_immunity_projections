# Immunity projections for vaccine preventable diseases in the Gaza Strip post-crisis

This repository contains *R* code to run a model based projection of population immunity to disease and infection from the following vaccine-preventable diseases:
- MMR:
    - Measles
- Pentavalent:
    - Diphtheria
    - Pertussis
    - Hib Disease
- Rotavirus
- PCV:
    - Pneumococcal Disease
- Polio:
    - Polio (Wild-type)
    - Polio (Vaccine-derived)

Based upon data on vaccine coverage from 2000 onwards, disease incidence, and duration of immunity.

## Requirements
This repo requires a modern version of *R*, and the installation of non-*CRAN* *R*-packages, i.e. [IVODE](https://github.com/GBarnsley/IVODE), so may require [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

## Repo Strucure
- `analysis`
  - `x_run.R` *main script for running this repo*
  - `0_setup.R` *installs and loads R-packages and sets up projection parameters*
  - `1_demographics.R` *this and the following script compile parameters for the projection model*
  - `2_diseases.R`
  - `3_projections.R` *runs the immunity projection model*
  - `4_format.R` *formats the output of the model*
  - `5_report.rmd` *Rmarkdown report that summarizes the projections*
  - `6_sensitivity.R` *additional sensitivity analyses*
  - `scenario_funcs.R` *miscellaneous R functions*
- `plots`
- `data`
  - `raw`
  - `derived`
  - `output`
