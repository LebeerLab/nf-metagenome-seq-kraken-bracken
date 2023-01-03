#!/usr/bin/env Rscript

# This scripts filters the reads of a set of shotgun samples based on length and
# quality.

# Author: Tim Van Rillaer, adapted from Stijn Wittouck
# Last modified: 02/01/2023

# dependencies: dada2 v1.20.0, tidyverse v1.3.0

library(dada2)
library(tidyverse)

# set filtering and trimming parameters
truncLen <- c(0, 0)
trimLeft <- c(0, 0)
trimRight <- c(0, 0)
minLen <- c(50, 50)
maxN <- c(2, 2)
maxEE <- c(2, 2)

args = commandArgs(trailingOnly = TRUE)
filterAndTrim(args)