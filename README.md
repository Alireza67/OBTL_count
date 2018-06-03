# Optimal Bayesian Transfer Learning for Count Data
This is the Matlab code for the paper "Optimal Bayesian Transfer Learning for Count Data". You can run the file named demo.m to reproduce the resluts for the paper. Before running the demo file, make sure to install MatlabStan and add its path in the first line of the demo file. To install MatlabStan, follow the instructions in https://github.com/brian-lau/MatlabStan/wiki.

The source and target data for the two types of lung cancer (LUAD and LUSC) from The Cancer Genome Atlas (TCGA) have been provided. The data have been partitioned to 50 random training and test sets for two different feature sets, as described in the paper. By setting the  prior correlations of the mean and shape parameters, $\rho_\mu$ and $\rho_r$, in the demo file, the corresponding average classification error for the OBTL classifier and OBC is computed. 
