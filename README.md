
## Overview

This repository contains all data, code, and materials required to reproduce the code example, applied examples, and simulation study reported in:

*Mathur MB, Smith LH, Yoshida K, Ding P, VanderWeele TJ (under review). E-values for effect modification and approximations for causal interaction.*

## How to reproduce the applied example on education and dementia incidence

We obtained estimates and 2 x 2 table cell counts from information reported in [Letenneur et al. (2010)](https://pubmed.ncbi.nlm.nih.gov/10873130/). We used the script [analyze.R](https://osf.io/d6wub/) to enter these estimates and cell counts and to conduct analyses. This script calls helper functions in [helper.R](https://osf.io/nrhx7/). Because those functions are limited working versions of the generalized versions in the R package EValue, users wishing to conduct their own analyses should use the R package. 

