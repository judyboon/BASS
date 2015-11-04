# BASS
This repository contains the implementation of _BASS_,  a group factor analysis model with Bayesian structured sparsity priors.
For details of the model, please see http://arxiv.org/abs/1411.2698.

# Introduction
_BASS_ stands for Bayesian group factor Analysis with Structured Sparsity priors.
It extends traditional factor analysis model to _m_ coupled observation matrices (views).

BASS could be used to estimate both covariance specific to each observation matrix 
and shared variation across any subset of observations. 
Applications of BASS to construct condition specific gene co-expression networks and to
find shared topics among different classes of documents have been demonstrated in the paper.


# About the folders
* _simData_ contains simulation data sets used in the original paper.
* _code_ contains all source code of implementation.
* _output_ contains the output files/folders after running BASS.

## Code
### Compilation
BASS is implemented using C++. Use Makefile to compile the source code to an executable file. 
It requires installation of GNU Scientific Library (GSL). 
Please make sure the library has been installed in your OS.

### Running
Shell script _runBASS.sh_ provides an example of running BASS with simulation data set "Yn100sim5.sim".
It generates following commands
> ./main -y ../simData/Yn100sim5.sim -sep space -k 15 
  > -v 10 -p1 50 -p2 50 -p3 50 -p4 50 -p5 50 -p6 50 -p7 50 -p8 50 -p9 50 -p10 50
  > -n 100 -iter 20000 -step 50 -mcmc 0 -px 20 
  > -output ../output/MCMC0PX20ik15n100

Explanation of arguments:
* -y: joint data matrix file (p by n) with variables (features) are represented by row and subjects are in column
* -sep: how the data are separated. arguments: "space" "tab"
* -k: initial number of latent factors
* -v: number of observations (views)
* -p1 up to -pv: the dimension of each observation
* -n: number of subjects
* -iter: maximum iteration allowed
* -step: number of steps to store parameters
* -mcmc: number of MCMC steps
* -px: number of PX-EM steps
* -output: output directory

### Analyzing results
R script _analyze.R_ provides an example of analyzing outputs of BASS.

### Acknowledgement
The source code _RngStream.cpp_, _RngStream.h_ and *gig_par.cpp*, *gig_par.h* are downloaded from
http://jonaswallin.github.io/articles/2013/07/simulation-of-gig-distribution/ to draw samples from generalized inverse Gaussian distribution. 
