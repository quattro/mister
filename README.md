# mister
Small simulation study of various Mendelian Randomization (MR) approaches using summary GWAS data

Overview
---------

| File | Description |
|--------|-------------|
| `estimators.R` | Contains various MR estimators (e.g., IVW, Egger, Likelihood, Bayes-Egger) |
| `simulation.R` | Standard simulation of outcome mediated by exposure with shared confounder |
| `simulation.balanced.R` | Standard simulation with IVs having balanced (i.e. 0-centered) direct effects on outcome | 
| `simulation.directional.R` | Standard simulation with IVs having directional (i.e. non-zero on average) direct effects on outcome |

To run a simulation load an R interactive session and enter `source(simulation.R)` (or one of the alternatives). Results are stored in `results` and a figure is generated under `plot`.
