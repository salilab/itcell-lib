This is the source code for ITCell, a method for integrative T-cell epitope
prediction.
See [D. Schneidman-Duhovny et al.](https://doi.org/10.1101/415661) for details.

A [web server](https://salilab.org/itcell/) is also available.

# Prerequisites

ITCell will likely only work on a Linux system.

 - [IMP](https://integrativemodeling.org/) (the scripts expect to find the
   `soap_score` binary in the PATH).
 - [SCRWL 4](http://dunbrack.fccc.edu/scwrl4/) (the scripts expect the SCWRL
   executable to be found as `scripts/scwrl4/Scwrl4`).
 - [MultiProt](http://bioinfo3d.cs.tau.ac.il/MultiProt/) (the scripts expect
   to use `scripts/multiprot.Linux` and were developed with version 1.93)
