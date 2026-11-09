# Publication simulation entry point.
#
# Runs every simulation needed for the publication figures:
#   1. The three cage-trial scenarios (AC 100, DCD 15, DCD 40), driven by
#      Simulation_Instructions.csv and the matching Parameter Files/<scenario>/
#      folders. Outputs land in Outcomes/<scenario>/.
#   2. The nine confinement runs (Standard / Split / DCD x Pop1 / Pop2 / Pop3),
#      a three-population chain with back-and-forth migration, driven by
#      run_confinement_chain.R. Outputs land in Outcomes/<sim>/.
#
# Designed to be self-contained inside this repository. Run with:
#     setwd("Model/Shared Simulations/Publication")
#     source("Simulation setup.R")
# or with Rscript:
#     Rscript "Model/Shared Simulations/Publication/Simulation setup.R"

rm(list = ls())

library(dplyr)
library(tidyr)
library(Matrix)
library(foreach)
library(doParallel)
library(readr)

# Resolve this script's directory so the simulation runs from the right place
# regardless of how it is invoked.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  current_path <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    current_path <- sub("--file=", "", file_arg)
    current_path <- URLdecode(gsub("~\\+~", "%20", current_path))
    current_path <- normalizePath(current_path)
  } else {
    current_path <- ""
  }
}
if (nzchar(current_path)) setwd(dirname(current_path))
cat("Working directory:", getwd(), "\n")

# project_root points at Model/ (which holds Model Code.R and the top-level
# Default_*.csv fallbacks). From this script's location that is two levels up.
project_root <- normalizePath("../..", mustWork = TRUE)

source(file.path(project_root, "Model Code.R"))

# Cage trial runs ------------------------------------------------------------

# When TRUE, conversion / per-genotype matrices are also written out.
produce_conversion_table_files <- TRUE

Figures_to_make <- as.data.frame(read.csv(
  file = "Simulation_Instructions.csv", header = TRUE,
  stringsAsFactors = FALSE
))

run_simulations(create_figure_instruc(Figures_to_make), project_root = project_root)

# Confinement runs -----------------------------------------------------------

source("run_confinement_chain.R")

cat("\nAll publication simulations complete.\n")
