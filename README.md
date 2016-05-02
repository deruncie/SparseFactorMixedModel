## BSFG
by Daniel Runcie, Sayan Mukhergee

Reference:
Runcie, D., & Mukherjee, S. (2013). Dissecting High-Dimensional Phenotypes with Bayesian Sparse Factor Analysis of Genetic Covariance Matrices. Genetics, 194(3), 753–767. http://doi.org/10.1534/genetics.113.151217

### Version history

#### V1.0
Published version in MATLAB
on website: ![](http://www2.stat.duke.edu/~sayan/bfgr/index.shtml)
- includes Ayroles_et_al_Competitive_fitness, Simulations with half-sib design
- should be able to replicate all analyses from paper (up to Monte-carlo error in Gibbs and in simulations)

#### V1.1
- Fixed calculation of genetic and interaction specific effects. The calculation and corresponding text of the paper missed $A^{-1}$. This should not greatly affect the results of the analyses presented in the paper, but will need to check. It doesn't affect Ayroles analysis.

#### V2
Nearly complete re-write of the model code, but should maintain identical function (I believe)
- variables have been re-named to more closely correspond to the paper
- sampler function has been re-written to only sample.
- A new function initializes the sampler, only run once
- sampler function starts where the previous run left off (including maintaining the random number generator), so should be the same as running one continuous chain