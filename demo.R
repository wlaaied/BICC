#call for packages, load data, run the function, and see the output

# install packages 
library(tidyverse)
library(MCMCpack)
library(TruncatedNormal)
library(extraDistr)
library(MASS)
library(dplyr)
library(rstan)
library(ordinal)
library(marked)
library(lme4)
library(nlme)
library(truncnorm)
library(chemometrics)
library(rpart)
library(invgamma)

# read BICC function from a file 
source('BICC_function.R')

# read data
data = read.csv('data.csv')

# run function
# Arguments
# data = data frame of data
# niter = number of iterations
# nchain = number of Gibbs chains
# nburn = number of burn-in
BICC(data = data, niter = 50000, nchain = 2, nburn = 1000)
