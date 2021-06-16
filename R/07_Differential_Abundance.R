##  ###################################################  ##
##  Differential abundance of taxa                       ##
##                                                       ##
##  Author: Geoff Zahn - June 14, 2021                   ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  tidyverse v 1.3.0                                    ##
##  phyloseq v 1.32.0                                    ##
##  vegan v 2.5.6                                        ##
##  broom v 0.7.1                                        ##
##  patchwork v 1.0.1                                    ##
##  microbiome v 1.10.0                                  ##
##  lme4 v 1.1.23                                        ##
##  lmerTest v 3.1.3                                     ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(lme4); packageVersion("lme4")
library(lmerTest); packageVersion("lmerTest")
library(corncob); packageVersion("corncob")
source("./R/bbdml_helper.R")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585")

# Load ps objects
ps <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

meta(ps_genus) %>% glimpse()

# reassign factor level order
ps_genus@sam_data$Structure <- factor(ps_genus@sam_data$Structure, levels = c( "Sediment","Pneumatophore","Leaf","Fruit"))

# corncob analyses ####

# genus-Level vs "Structure"

set.seed(123)
da_analysis <- differentialTest(formula = ~ Structure, #abundance
                                phi.formula = ~ Structure, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

beepr::beep()
plot(da_analysis)

if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps,
                           mu_predictor = "Acidified",
                           phi_predictor = "Acidified",
                           taxlevels = 2:7)
  length(bbdml_obj)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "Acidified",
                   pointsize = 3)
}

