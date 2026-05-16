# =============================================================================
# Nature_Revision2 / 00_config.R
# -----------------------------------------------------------------------------
# I use this as the single place to set paths for my NR2 workflow. I source it
# from the other NR2 scripts so everything writes under Outputs/NR2.
#
# Inputs:  none (optional: I set `project_root` in the global environment before
#           sourcing if I am not running from the project root).
# Outputs: I only create directory paths; I do not write data files here.
#           - Outputs/NR2/figures
#           - Outputs/NR2/rds
#           - Outputs/NR2/models
# =============================================================================

#---- NR2 path configuration ----
# NR2 INPUT: If running this workflow outside the project root, set `project_root`
# manually before sourcing this file.
if (!exists("project_root")) {
  project_root <- "."
}

raw_data_dir <- file.path(project_root, "RawData")
inputs_dir <- file.path(project_root, "Inputs")
functions_dir <- file.path(project_root, "Functions")

nr2_root <- file.path(project_root, "Outputs", "NR2")
nr2_figures_dir <- file.path(nr2_root, "figures")
nr2_rds_dir <- file.path(nr2_root, "rds")
nr2_models_dir <- file.path(nr2_root, "models")

dir.create(nr2_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nr2_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nr2_models_dir, recursive = TRUE, showWarnings = FALSE)

