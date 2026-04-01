# bfboinet

The `bfboinet` package implements the Backfill Bayesian Optimal Interval Design (BF-BOIN-ET), 
a novel clinical trial methodology for dose optimization that simultaneously consider both efficacy 
and toxicity outcome as described in Takeda et al  (2025) < https://doi.org/10.1002/pst.2470>.
The package has been extended to include a seamless two-stage phase I/II trial design with backfill 
and joint efficacy and toxicity monitoring as described in Takeda et al (2016) <DOI>.

## Installation

You can install the development version of bfboinet from GitHub with:
```r
# install.packages("devtools")
devtools::install_github("jingzhuzhuzhu/bfboinet")

---

### **3. Usage Example**
```markdown
## Basic Example

library(bfboinet)

# Simulate a dose-finding trial
result <- get.oc.backboinet(
  target_T = 0.3,
  toxprob,
  target_E = 0.25,
  effprob,
  n.dose,
  startdose,
  ncohort,
  cohortsize,
  pT.saf = 0.6 * target_T,
  pT.tox = 1.4 * target_T,
  pE.saf = 0.6 * target_E,
  alpha.T1 = 0.5,
  alpha.E1 = 0.5,
  tau.T,
  tau.E,
  te.corr = 0.2,
  gen.event.time = "weibull",
  accrual,
  gen.enroll.time = "uniform",
  n.elimination = 6,
  stopping.npts = 12,
  suspend = 0,
  stopping.prob.T = 0.95,
  stopping.prob.E = 0.9,
  ppsi01 = 0,
  ppsi00 = 40,
  ppsi11 = 60,
  ppsi10 = 100,
  n.sim = 1000,
  seed.sim = 100
)

