# Confinement Comparison - Three Linked Populations (chain)
#
# Three populations arranged in a chain, connected by back-and-forth migration
# and run simultaneously (lockstep), one generation at a time:
#
#     Pop 1  <-->  Pop 2  <-->  Pop 3  <-->  (wild-type general population)
#   (release)    (middle)     (far end)
#
# Pop 1 receives the release; Pops 2 and 3 start wild-type. Migration is
# exchanged at the gamete stage each generation (see Main_Model_connected),
# with the far end (Pop 3) supplemented by wild-type gametes. This mirrors the
# two-population model of Verkuijl et al. 2022 (PLoS Genetics; OSF 6a2eg),
# generalised to a three-population chain.
#
# Three drives: Standard_Drive (global), Split_Drive, DCD (Direct Daisy-Chain).
# Each (drive, population) pair is saved to its own folder under Outcomes/,
# which is where Analysis.Rmd reads them from.
#
# This file is sourced by Simulation setup.R, but is also safe to source on
# its own. When sourced standalone it locates Model Code.R two levels up.

library(dplyr)
library(tidyr)
library(readr)
library(Matrix)

if (!exists("project_root") || !nzchar(project_root)) {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    current_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                             error = function(e) NULL)
    if (!is.null(current_path) && nzchar(current_path)) {
      setwd(dirname(current_path))
    }
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      current_path <- sub("--file=", "", file_arg)
      current_path <- URLdecode(gsub("~\\+~", "%20", current_path))
      setwd(dirname(normalizePath(current_path)))
    }
  }
  project_root <- normalizePath("../..", mustWork = TRUE)
}

if (!exists("Main_Model", mode = "function")) {
  source(file.path(project_root, "Model Code.R"))
}

cat("\nConfinement runs - working directory:", getwd(), "\n")
cat("Confinement runs - project_root:", project_root, "\n")

# --- Configuration ----------------------------------------------------------

migration_rate         <- 0.01   # back-and-forth migration each generation
n_populations          <- 3      # chain: Pop 1 (release) <-> Pop 2 <-> Pop 3 <-> WT
pop1_release_fraction  <- 0.50   # release fraction per sex in Pop 1
wt_supplementation     <- TRUE   # far end supplemented with wild-type gametes
simulation_length      <- 70     # generations to simulate (past the far-end
                                  # DCD payload peak at ~gen 54, so its decline
                                  # is visible in Pop 3)

drives <- list(
  list(name = "Standard_Drive", folder = "Standard_Drive"),
  list(name = "Split_Drive",    folder = "Split_Drive"),
  list(name = "DCD",            folder = "DCD")
)

# --- Helpers ----------------------------------------------------------------

# Save the standard set of per-simulation files into Outcomes/<sim_name>/,
# matching what run_simulations() writes for the cage trial scenarios.
save_sim_outputs <- function(sim_name, outcomes_mf, mm, params) {
  out_dir <- file.path("Outcomes", sim_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  outcomes_mf$Factor_F1 <- 1
  outcomes_mf$Factor_F2 <- 1

  Alpha_Allele_Names <- do.call(paste, as.data.frame(mm$Alpha_M[, 1:11], stringsAsFactors = FALSE))
  matrix_named <- function(M) {
    df <- cbind(Alpha_Allele_Names, as.data.frame(as.matrix(M)))
    names(df) <- c("Alleles", Alpha_Allele_Names)
    row.names(df) <- Alpha_Allele_Names
    df
  }

  write.csv(outcomes_mf,
            file.path(out_dir, paste0("Outcomes_", sim_name, "_Alleles.csv")),
            row.names = FALSE)

  write.csv(t(rbind(mm$geno_list$A, mm$geno_list$B, mm$geno_list$C)),
            file.path(out_dir, paste0("Allele_", sim_name, "_Alleles.csv")))

  write.csv(mm$nuc_rates_list,
            file.path(out_dir, paste0("Nuc_", sim_name, "_Alleles.csv")))

  write.csv(mm$Release_File,
            file.path(out_dir, paste0("Release_", sim_name, "_Alleles.csv")))

  write.csv(matrix_named(mm$D),
            file.path(out_dir, paste0("Dep_", sim_name, "_Alleles.csv")))
  write.csv(matrix_named(mm$exp_conv_male),
            file.path(out_dir, paste0("ExpM_", sim_name, "_Alleles.csv")))
  write.csv(matrix_named(mm$exp_conv_female),
            file.path(out_dir, paste0("ExpF_", sim_name, "_Alleles.csv")))
  write.csv(cbind(mm$Alpha_M, as.matrix(mm$H_M)),
            file.path(out_dir, paste0("HaploidM_", sim_name, "_Alleles.csv")))
  write.csv(cbind(mm$Alpha_F, as.matrix(mm$H_F)),
            file.path(out_dir, paste0("HaploidF_", sim_name, "_Alleles.csv")))
  write.csv(cbind(mm$Alpha_M, as.matrix(mm$Fit_Male)),
            file.path(out_dir, paste0("FitnessM_", sim_name, "_Alleles.csv")))
  write.csv(cbind(mm$Alpha_F, as.matrix(mm$Fit_Female)),
            file.path(out_dir, paste0("FitnessF_", sim_name, "_Alleles.csv")))
}

# Build the per-population Release_File list from a base Release_File.
# Pop 1 carries the release; Pops 2..n are wild-type with zero release.
make_release_files <- function(base_release, n_pop, release_fraction, sim_length) {
  lapply(1:n_pop, function(p) {
    rf <- base_release
    rf$Simulation_Length <- sim_length
    if (p == 1) {
      rf$Male_Release_Fraction   <- release_fraction
      rf$Female_Release_Fraction <- release_fraction
    } else {
      rf$Release_Type            <- "WT"
      rf$Male_Release_Fraction   <- 0
      rf$Female_Release_Fraction <- 0
    }
    rf
  })
}

# --- Main loop --------------------------------------------------------------

for (d in drives) {
  cat("\n===", d$name, "===\n")

  params <- all_default_parameters(
    Allele_Fit_File_Name = paste0(d$folder, "/Alleles.csv"),
    Nuc_File_Name        = paste0(d$folder, "/Nuc.csv"),
    Release_File_Name    = paste0(d$folder, "/Release.csv"),
    project_root         = project_root
  )

  params$Release_File$Simulation_Length <- simulation_length
  params$Generations                    <- simulation_length

  # Shared drive matrices (all populations use the same drive)
  mm <- All_Main_Calculations(params)

  # Per-population release files (Pop 1 releases; the rest are wild-type)
  release_files <- make_release_files(params$Release_File, n_populations,
                                      pop1_release_fraction, simulation_length)

  cat("  Chain of", n_populations, "populations, migration rate =", migration_rate,
      if (wt_supplementation) "(far end WT-supplemented)\n" else "\n")

  results <- Main_Model_connected(
    Model_Matrices     = mm,
    release_files      = release_files,
    migration_rate     = migration_rate,
    wt_supplementation = wt_supplementation,
    replicate          = 0
  )

  for (pop_i in 1:n_populations) {
    sim_name <- paste0(d$name, "_Pop", pop_i)

    # Record the per-population release file in the saved outputs
    mm$Release_File <- release_files[[pop_i]]

    alleles <- Outcomes_by_Allele(results[[paste0("pop_", pop_i)]], params, 0)
    save_sim_outputs(sim_name, alleles$Outcomes_MF, mm, params)
    cat("    saved Outcomes/", sim_name, "/\n", sep = "")
  }
}

cat("\nConfinement runs complete.\n")
