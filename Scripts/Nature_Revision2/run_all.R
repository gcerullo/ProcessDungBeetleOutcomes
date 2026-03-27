# Run the full NR2 workflow from one script.
# Usage (from project root):
#   Rscript "Scripts/Nature_Revision2/run_all.R"

scripts_dir <- file.path("Scripts", "Nature_Revision2")

script_sequence <- c(
  "01_CombineDungBeetleDatasets.R",
  "02_FormatDBForNegBinom.R",
  "03_FitModel.R",
  "04_Nature_PredictAbundanceByHab.R",
  "05_Nature_ProcessScenarioOutcomes.R"
)

cat("Starting NR2 workflow\n")
cat("Scripts directory:", normalizePath(scripts_dir, winslash = "/", mustWork = FALSE), "\n")

for (script_name in script_sequence) {
  script_path <- file.path(scripts_dir, script_name)

  if (!file.exists(script_path)) {
    stop("Missing script: ", script_path, call. = FALSE)
  }

  cat("\n--- Running:", script_name, "---\n")
  source(script_path, echo = TRUE)
}

cat("\nNR2 workflow complete.\n")
