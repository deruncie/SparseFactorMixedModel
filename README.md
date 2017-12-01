## BSFG
by Daniel Runcie, Sayan Mukhergee

Reference:
Runcie, D., & Mukherjee, S. (2013). Dissecting High-Dimensional Phenotypes with Bayesian Sparse Factor Analysis of Genetic Covariance Matrices. Genetics, 194(3), 753–767. http://doi.org/10.1534/genetics.113.151217

This package is free software, you can redistribute it and/or modify it 
under the terms of the GNU General Public License (GPL-3).

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

### Version history

#### V1.0
Published version in MATLAB
on website: [](http://www2.stat.duke.edu/~sayan/bfgr/index.shtml)
- includes Ayroles_et_al_Competitive_fitness, Simulations with half-sib design
- should be able to replicate all analyses from paper (up to Monte-carlo error in Gibbs and in simulations)

#### V1.1
- Fixed calculation of genetic and interaction specific effects. The calculation and corresponding text of the paper missed $A^{-1}$. 
This mistake, as well as other errors in the paper are documented [here](Correction/Correction_items.pdf). 
A re-analysis of the simulations presented in the paper and an updated Appendix are presented [here](Correction/Corrected_Appendix.pdf).

#### V2.0
Nearly complete re-write of the model code, but should maintain identical function \
- variables have been re-named to more closely correspond to the paper
- sampler function has been re-written to only sample.
- A new function initializes the sampler, only run once
- sampler function starts where the previous run left off (including maintaining the random number generator), so should be the same as running one continuous chain

#### R_BSFG V2.0
R clone of V2.0 Matlab code
- Functionality should be identical. Worth checking. Note that the RNG is different.
- embeded in the Gibbs sampler are two versions of each sampler function, a native R version, and a Rcpp version. They should be identical (up to RNG differences). The Rcpp function has the same name and arguments, but with "_c" appended to the function name.
