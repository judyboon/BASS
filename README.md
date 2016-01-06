# BASS
This repository contains the implementation of _BASS_,  a group factor analysis model with Bayesian structured sparsity priors.
For details of the model, please see http://arxiv.org/abs/1411.2698. 
The draft is published on Journal of Machine Learning Research (JMLR).

## License
This software is distributed under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 
of the License, or (at your option) any later version.

## JMLR Copyright
This source code is supplied “as is” without warranty of any kind, 
and its author and the Journal of Machine Learning Research (JMLR) 
and JMLR’s publishers and distributors, disclaim any and all warranties, 
including but not limited to any implied warranties of merchantability 
and fitness for a particular purpose, and any warranties or non 
infringement. The user assumes all liability and responsibility for use 
of this source code, and neither the author nor JMLR, nor JMLR’s 
publishers and distributors, will be liable for damages of any kind 
resulting from its use. Without limiting the generality of the foregoing, 
neither the author, nor JMLR, nor JMLR’s publishers and distributors, 
warrant that the source code will be error-free, will operate without 
interruption, or will meet the needs of the user.


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
Make sure the library has been installed in your OS and 
the path has been correctly specified in _Makefile_.

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
