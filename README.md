[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3226992.svg)](https://doi.org/10.5281/zenodo.3226992)

This is the source code for ITCell, a method for integrative T-cell epitope
prediction.
See [D. Schneidman-Duhovny et al.](https://doi.org/10.1101/415661) for details.

A [web server](https://salilab.org/itcell/) is also available.

# Installation

ITCell will likely only work on a Linux system.

First, download or clone this repository. Then download the prerequisites:

 - [IMP](https://integrativemodeling.org/) (the scripts expect to find the
   `soap_score` binary in the PATH).
 - [SCWRL 4](http://dunbrack.fccc.edu/scwrl4/) (extract SCWRL into the
   `scripts/scwlr4` subdirectory of this repository; the scripts expect the
   SCWRL executable to be found as `scripts/scwrl4/Scwrl4`).
 - [MultiProt](http://bioinfo3d.cs.tau.ac.il/MultiProt/) (the scripts expect
   to use `scripts/multiprot.Linux` and were developed with version 1.93).

Next, download the support libraries from
<a href="https://salilab.org/itcell-lib/">https://salilab.org/itcell-lib/</a>
and extract into the `libs` directory.

# Usage

Run `scripts/runITCell.pl` giving it the MHC type, the path to a TCR PDB
file, and the path to a file containing the antigen sequence, e.g.

    perl <repo>/scripts/runITCell.pl DRB1*0101 TCR.pdb antigen_seq.txt

# License

ITCell is Copyright 2018 Dina Schneidman-Duhovny, and is available under the
terms of the GNU Lesser GPL; see file `LICENSE`.
