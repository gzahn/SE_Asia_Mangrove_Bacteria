##  ###################################################  ##
##           ##
##                                                       ##
##  Author: Geoff Zahn - October 12, 2020                ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  tidyverse v 1.3.0                                    ##
##  phyloseq v 1.32.0                                    ##
##  mikropml v 1.0.0                                     ##
##  microbioome v 1.10.0                                 ##
##  ###################################################  ##

library(mikropml)
library(tidyverse)
library(phyloseq)
library(microbiome); packageVersion("microbiome")

# import genus-level ps object
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")


# build dataframe for mikropml (1st column = feature, all other columns = ASVs, rows = samples)
df <- otu_table(ps_genus) %>% as.data.frame()
meta <- meta(ps_genus)
df$Structure <- meta$Structure

# pre-process data and run mikropml function
df_processed <- preprocess_data(df,outcome_colname = "Structure")
ml_results <- run_ml(df_processed$dat_transformed,
                     "glmnet",
                     outcome_colname = "Structure",
                     seed=2021,
                     find_feature_importance = TRUE)

ml_results$performance
ml_results$feature_importance
ml_results$trained_model
