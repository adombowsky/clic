# CLIC: Clustering with Independence Centering

This is the github repository for _Product Centered Dirichlet Processes for Bayesian Multiview Clustering._ The product centered Dirichlet process (PCDP) is a novel prior on the mixing measure that can generate partitions exhibiting CLustering with Independence Centering (CLIC). Here, you will find all of the main functions needed to implement CLIC in practice. 

1. The file `sampling.cpp` contains all RCPP code needed to simulate from the posterior of the partitions; `postprocessing.cpp` contains all code needed to compare the clusterings; and `threeview.cpp` contains special functions for the three view setting.
    * These functions require RCPP and RCPP Armadillo:
```r
install.packages("Rcpp")
install.packages("RcppArmadillo")
```
    * To load our functions, use
```r
  sourceCpp("rcppfuncts/sampling.cpp")
  sourceCpp("rcppfuncts/threeview.cpp")
  sourceCpp("rcppfuncts/postprocessing.cpp")
```
2. To run the simulation studies, see `multiview_simulations.R`, `threeview_simulations.R`, and `samp_size_and_dim_simulations.R`. For all simulations, the random number seed is set to `1` by default.
    * Note that the third simulation study requires parallel computing.
3. The application to the CPP data is implemented in `longnecker.R` (though the data are not included).
