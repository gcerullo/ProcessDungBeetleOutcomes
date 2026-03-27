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

