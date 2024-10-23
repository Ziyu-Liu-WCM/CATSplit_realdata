library(tidyverse)
library(ape)
library(GUniFrac)
library(MiRKAT)
library(vegan)
library(glmnet)
library(parallel)
library(doParallel)


## Read files
source("helper_function.R")
source("preprocess.R")
source("Split_fun.R")
source("DS_Mann_fun.R")
source("otutabletoabundance.R")
abundanceTable<-otutabletoabundance(otu_table)
apply(abundanceTable,2,sum) # double check if abundance sum to 1 for each obs
## You can load the result directly if you don't have enough time to run
#if(file.exists("result.RData")) load("result.RData")

## Mann-Whitney
Man_selected <- Mann_WhitU(abundanceTable, taxonomy_table, meta_table, qval_bound = 0.2)
Man_selected

# ## You can load the result directly if you don't have enough time to run
# if(file.exists("result.RData")) load("result.RData")
# 
# ## Mann-Whitney
# Man_selected <- Mann_WhitU(otu_table, taxonomy_table, meta_table, qval_bound = 0.05)
# Man_selected


## LASSO with data splitting
meta_bin <- ifelse(meta_table$BinOutcomes == "R", 1, 0)
DS_selected <- DSBin(otu_table, meta_bin, taxonomy_table, num_split = 50, qval_bound = 0.05)
DS_selected

## CATSplit
metric <- "Weighted UniFrac"
parallel <- TRUE
nCore <- 10
nReps <- 50

CAT_selected <- CATSplit_parallel(otu_table, taxonomy_table, meta_table, tree_data,
                                  metric, parallel, nCore, nReps, qval_bound = 0.05)
CAT_selected


## Save results
result <- list(Man_selected = Man_selected, DS_selected = DS_selected, CAT_selected = CAT_selected)
if(!file.exists("result.RData")) save(result, file = "result.RData")
