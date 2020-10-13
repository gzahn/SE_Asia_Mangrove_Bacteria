##  ###################################################  ##
##  Combine phyloseq objects from all runs into one      ##
##  This script processes "run7" - 1912KMI-0006          ##
##                                                       ##
##  Sonneratia alba microbiome                           ##
##                                                       ##
##  Author: Geoff Zahn - October 12, 2020                ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  tidyverse v 1.3.0                                    ##
##  phyloseq v 1.32.0                                    ##
##                                                       ##
##  ###################################################  ##

# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")

# Merge ####

# Grab ps objects
ps06 <- readRDS("./Output/run6/noncontam_ps_object.RDS")
ps07 <- readRDS("./Output/run7/noncontam_ps_object.RDS")
ps12 <- readRDS("./Output/run12/noncontam_ps_object.RDS")
ps13 <- readRDS("./Output/run13/noncontam_ps_object.RDS")

# merge phyloseq objects from all 4 runs
full_ps <- merge_phyloseq(ps06,ps07,ps12,ps13)
saveRDS(full_ps, "./Output/full_ps_object_raw.RDS")

# Clean up taxa and any empty samples ####

# remove non-bacteria
full_ps <- 
  full_ps %>% 
  subset_taxa(Kingdom == "Bacteria")

# Remove "Chloroplast" taxa
full_ps <- 
  full_ps %>% 
  subset_taxa(Order != "Chloroplast")

# remove empty samples
full_ps <- 
  full_ps %>% 
  subset_samples(sample_sums(full_ps) > 1)

# remove empty taxa
full_ps <- 
  full_ps %>% 
  subset_taxa(taxa_sums(full_ps) > 0)

saveRDS(full_ps, "./Output/full_ps_object_cleaned.RDS")