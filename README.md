# mister
Straightforward simulation study of various Mendelian Randomization (MR) approaches using summary GWAS data. Intent is mainly for research and code should not be used in real-world MR studies.

Overview
---------

| File | Description |
|--------|-------------|
| `estimators.R` | Contains various MR estimators (e.g., IVW, Egger, Likelihood, Bayes-Egger) |
| `simulation.R` | Standard simulation of outcome mediated by exposure with shared confounder |
| `simulation.balanced.R` | Standard simulation with IVs having balanced (i.e. 0-centered) direct effects on outcome | 
| `simulation.directional.R` | Standard simulation with IVs having directional (i.e. non-zero on average) direct effects on outcome |

To run a simulation load an R interactive session and enter `source(simulation.R)` (or one of the alternatives). Results are stored in `results` and a figure is generated under `plot`.

Caveats
-------
1. Simulations are currently set up such that no SNPs are in linkage (i.e. LD) with one another.
2. Outcome and exposure GWAS results are from the same population, but no shared environmental factors or modeled (in either data generation or estimation)

Methods References
------------------
1. [Burgess, Stephen, Adam Butterworth, and Simon G. Thompson. "Mendelian randomization analysis with multiple genetic variants using summarized data." Genetic epidemiology 37.7 (2013): 658-665.](https://doi.org/10.1002/gepi.21758)

2. [Bowden, Jack, George Davey Smith, and Stephen Burgess. "Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression." International journal of epidemiology 44.2 (2015): 512-525.](https://doi.org/10.1093/ije/dyv080)
