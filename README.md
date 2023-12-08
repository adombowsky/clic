# CLIC: Clustering with Independence Centering

This is the github repository for _Product Centered Dirichlet Processes for Dependent Clustering._ The product centered Dirichlet process (PCDP) is a novel prior on the mixing measure that can generate partitions exhibiting CLustering with Independence Centering (CLIC). Here, you will find all of the main functions needed to implement CLIC in practice. 

1. The file `sampling.cpp` contains all RCPP code needed to simulate from the posterior of the partitions.
2. To run the simulation studies, see `multiview_simulations.R` and `conditional_simulations.R`.
3. The application to the CPP data is implemented in `longnecker.R`. 