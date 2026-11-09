# DCD Cage Trial — Publication Package

Self-contained package for the publication: the population-genetics model
code, the parameter files for every published simulation, the figure-building
R Markdown, and the methods / figure-legend text. Everything is reproducible
from a fresh clone in a single R session.

## Folder layout

```
DCD-Cage-Trial-Publication/
├── README.md
├── LICENSE
├── Methods.Rmd                 ← supplemental mathematical-modelling methods
├── Analysis.Rmd                ← knits the figures + Simulation_Scenarios.xlsx
│                                 (figure legends live inline, below each chunk)
├── Figures/                    ← rendered vector figures, SVG (overwritten on knit)
│   ├── Confinement.*               (Figure 1: 3 drives × 3 chained populations)
│   └── Fitted_r1_Overlay.*         (Figure 2: cage trial overlay with fitted r1)
├── Simulation_Scenarios.xlsx   ← per-replicate trajectories (written on knit)
├── renv.lock + .Rprofile + renv/  ← pinned R + package versions (renv)
└── Model/
    ├── Model Code.R                ← core simulation code
    ├── Default_Alleles.csv         ← library-level allele defaults
    ├── Default_Nuc.csv             ← library-level nuclease defaults
    ├── Default_Release.csv         ← library-level release defaults
    ├── LICENSE
    └── Shared Simulations/Publication/
        ├── Simulation setup.R          ← one-shot entry point for all sims
        ├── Simulation_Instructions.csv ← drives the 3 cage trial runs
        ├── run_confinement_chain.R     ← drives the 9 confinement runs (3-population chain)
        ├── Parameter Files/
        │   ├── Default_Alleles.csv     ← project-level allele overrides
        │   ├── Default_Nuc.csv         ← project-level nuclease overrides
        │   ├── Default_Release.csv     ← project-level release overrides
        │   ├── AC 100 Fitted r1_.1/    (Alleles, Release, Experimental_Data)
        │   ├── DCD 15 Fitted r1_.1/    (Alleles, Release, Experimental_Data)
        │   ├── DCD 40 Fitted r1_.1/    (Alleles, Release, Experimental_Data)
        │   ├── Standard_Drive/         (Alleles only)
        │   ├── Split_Drive/            (Alleles only)
        │   └── DCD/                    (Alleles only)
        └── Outcomes/                   ← simulation outputs (written by setup.R)
            ├── AC 100 Fitted r1_.1/    Outcomes + Experimental + ancillary CSVs
            ├── DCD 15 Fitted r1_.1/    Outcomes + Experimental + ancillary CSVs
            ├── DCD 40 Fitted r1_.1/    Outcomes + Experimental + ancillary CSVs
            ├── Standard_Drive_Pop{1,2,3}/
            ├── Split_Drive_Pop{1,2,3}/
            ├── DCD_Pop{1,2,3}/
            └── Meta_Table*.csv         ← per-run parameter audit trail
```

## Reproducing the figures

The repository ships with an `renv` lockfile pinning the R version (4.5.2) and
every package used to produce the published figures. On a fresh machine:

```r
# 1. Open the project folder in R / RStudio. The .Rprofile auto-activates renv.
# 2. Install the pinned package set:
renv::restore()
```

Then either re-knit straight from the bundled simulation outputs, or rerun the
simulations from scratch.

### Just the figures (uses bundled outputs)

```r
rmarkdown::render("Analysis.Rmd")
```

Writes/overwrites:

- `Figures/Confinement.svg`
- `Figures/Fitted_r1_Overlay.svg`
- `Simulation_Scenarios.xlsx`
- `Analysis.pdf` (the `Analysis.Rmd` YAML sets `output: pdf_document`)

### Rerun every simulation, then knit

```r
setwd("Model/Shared Simulations/Publication")
source("Simulation setup.R")
setwd("../../..")
rmarkdown::render("Analysis.Rmd")
```

Or, equivalently, from a shell:

```sh
Rscript "Model/Shared Simulations/Publication/Simulation setup.R"
```

`Simulation setup.R` writes its outputs directly into `Model/Shared Simulations/Publication/Outcomes/`, which is exactly where `Analysis.Rmd` reads them — no files need to be moved.

## What the simulations cover

| Group | Scenario | Output folder under `Outcomes/` |
|-------|----------|----------------------------------|
| Confinement (Figure 1) | Global drive       | `Standard_Drive_Pop1/2/3` |
| Confinement (Figure 1) | Split drive        | `Split_Drive_Pop1/2/3` |
| Confinement (Figure 1) | Daisy chain drive  | `DCD_Pop1/2/3` |
| Cage trial (Figure 2)  | AC 100% (Cas9-absent control) | `AC 100 Fitted r1_.1` |
| Cage trial (Figure 2)  | DCD 15% release    | `DCD 15 Fitted r1_.1` |
| Cage trial (Figure 2)  | DCD 40% release    | `DCD 40 Fitted r1_.1` |

Each cage trial output folder contains one deterministic + 50 stochastic
simulation replicates (`Outcomes_<scenario>_Alleles.csv`) plus the experimental
cage trial measurements (`Experimental_<scenario>_Alleles.csv`). The
confinement folders each contain a single deterministic time series.

## Parameter lookup order

For each scenario, the model resolves parameters in this order:

1. `Parameter Files/<scenario>/<Alleles|Nuc|Release>.csv` (scenario-specific
   override)
2. `Parameter Files/Default_<Alleles|Nuc|Release>.csv` (project-level defaults
   for the publication — includes the cage-trial-specific fitness / phenotype
   settings)
3. `Model/Default_<Alleles|Nuc|Release>.csv` (library-level fallback)

The publication’s parameter table in `Methods.Rmd` is generated at knit time
from these files, so any update to a parameter file is reflected automatically
in the rendered methods.

## Citation

If you use the model or this package, please cite Verkuijl et al., PLOS
Genetics (original framework) and the present manuscript (cage trial
application and fitted-r1 analysis).
