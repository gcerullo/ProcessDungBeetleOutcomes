# Nature_Revision2

This folder contains full copied scripts for the Nature revision workflow:

0. `00_config.R` (shared path/input-output config for all NR2 scripts)
1. `run_all.R` (single entrypoint to execute scripts 01 to 06 in order)
2. `01_CombineDungBeetleDatasets.R` (copied from `Scripts`)
3. `02_FormatDBForNegBinom.R` (copied from `Scripts`)
4. `03_FitModel.R` (copied from `Scripts`)
5. `04_Nature_PredictAbundanceByHab.R` (copied from `Scripts/Nature_revision`)
6. `05_Nature_ProcessScenarioOutcomes.R` (copied from `Scripts/Nature_revision`)
7. `06_scenario_uncertainty_plots.R` (copied from `Scripts`)

## Replicability Path Rule

Each script now sources shared config:

- `source(file.path("Scripts", "Nature_Revision2", "00_config.R"))`
- In `00_config.R`, set `project_root` if you are not running from the project root.

Update `project_root` in that block when needed.  
All script outputs are routed into:

- `Outputs/NR2/figures`
- `Outputs/NR2/rds`
- `Outputs/NR2/models`

## Script-by-Script Inputs, Outputs, and Summary

### 1) `01_CombineDungBeetleDatasets.R`

- Inputs
  - `RawData/*` (all raw dung beetle datasets used in this script)
  - `Inputs/*` (name backbone, GPS, plantation metadata)
  - `Functions/DB_functions.R`
- Outputs
  - `Outputs/NR2/rds/dungBeetlesForAbundanceAnalysis.csv`
  - `Outputs/NR2/rds/data_for_annna_db_borneo.csv`
- Summary
  - Cleans and harmonizes all raw dung beetle datasets.
  - Joins taxonomy, GPS, and intervention-age metadata into one master analysis table.

### 2) `02_FormatDBForNegBinom.R`

- Inputs
  - `Outputs/NR2/rds/dungBeetlesForAbundanceAnalysis.csv`
- Outputs
  - `Outputs/NR2/rds/full_DB_dataframeFor_BRMS_analysis_withSingletons.csv`
  - `Outputs/NR2/rds/full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv`
- Summary
  - Converts trap/day records to model-ready summaries.
  - Scales time-since-intervention and creates versions with/without singleton/doubleton species.

### 3) `03_FitModel.R`

- Inputs
  - `Outputs/NR2/rds/full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv`
- Outputs
  - `Outputs/NR2/models/DB_zi_full_Nov24.rds` (via `brms::brm(..., file = ...)`)
- Summary
  - Fits the zero-inflated negative binomial model in `brms`.
  - Saves fitted model object for prediction scripts.

### 4) `04_Nature_PredictAbundanceByHab.R`

- Inputs
  - `Outputs/NR2/rds/full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv`
  - `Outputs/NR2/models/DB_zi_full_Nov24.rds`
- Outputs
  - `Outputs/NR2/rds/processedOccBeetles_plateau.rds`
  - `Outputs/NR2/rds/processedOccBeetles_no40yrplateau.rds`
  - `Outputs/NR2/rds/DBs_abundance_by_habAge_iterations.rds`
  - `Outputs/NR2/rds/DBsppCategories.rds`
  - `Outputs/NR2/figures/all_beetle_curves_clipped_40yrs.png`
- Summary
  - Generates species-level posterior predicted abundance draws by habitat and age.
  - Builds species categories and saves both intermediate and plotting outputs.

### 5) `05_Nature_ProcessScenarioOutcomes.R`

- Inputs
  - `Inputs/ScenarioParams.R`
  - `Inputs/MasterAllScenarios.rds`
  - `Inputs/ScenariosWithDelaysCSVs/*.csv`
  - `Outputs/NR2/rds/DBs_abundance_by_habAge_iterations.rds`
  - `Outputs/NR2/rds/DBsppCategories.rds`
  - `Outputs/NR2/rds/SL_occ60yr_perIterationJan25.rds` (created within script and reused later in script)
- Outputs
  - `Outputs/NR2/rds/SL_occ60yr_perIterationJan25.rds`
  - `Outputs/NR2/rds/Ab60perScenarioIterationJan25/*.rds`
  - `Outputs/NR2/rds/RelativeAbundancePerIterationJan25/*.rds`
  - `Outputs/NR2/rds/BestScenarioUncertainty/*.rds`
  - `Outputs/NR2/rds/MasterDBPerformance.rds`
- Summary
  - Combines model-derived species responses with scenario trajectories.
  - Produces scenario-level uncertainty outputs and final biodiversity performance summaries.

### 6) `06_scenario_uncertainty_plots.R`

- Inputs
  - `Outputs/NR2/rds/BestScenarioUncertainty/*.rds`
  - `Outputs/NR2/rds/DBsppCategories.rds`
  - `selected_file_index` inside script (choose which `BestScenario` file to plot)
- Outputs
  - `Outputs/NR2/figures/loser_dungbeetle_uncertainty_plot.pdf`
  - `Outputs/NR2/figures/int1lgrp_dungbeetle_uncertainty_plot.pdf`
  - `Outputs/NR2/figures/all_sp_uncertainty_plot.pdf`
  - `Outputs/NR2/figures/all_sp_chunk1.pdf`
  - `Outputs/NR2/figures/all_sp_chunk2.pdf`
  - Same filenames exported as `.png`
- Summary
  - Summarises posterior uncertainty in whether logging or plantation strategies are best per species.
  - Generates uncertainty figures for loser species, intermediate group, and all species.

## Recommended Run Order

Run scripts in numeric order from this folder:

0. `00_config.R` (auto-sourced in scripts, no need to run manually)
or run everything with:

- `Rscript "Scripts/Nature_Revision2/run_all.R"`

Manual order:

1. `01_CombineDungBeetleDatasets.R`
2. `02_FormatDBForNegBinom.R`
3. `03_FitModel.R`
4. `04_Nature_PredictAbundanceByHab.R`
5. `05_Nature_ProcessScenarioOutcomes.R`
6. `06_scenario_uncertainty_plots.R`
