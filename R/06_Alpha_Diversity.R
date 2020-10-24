##  ###################################################  ##
##  Beta-diversity measures - btwn site and structure    ##
##                                                       ##
##  Author: Geoff Zahn - October 20, 2020                ##
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
##  lmerTest v 3.1.2                                     ##
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

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585")

# Load ps objects
ps <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

names(ps@sam_data)

# Plot alpha diversity ####
plot_richness(ps_genus, x="Structure", measures=c("Shannon")) +
  geom_boxplot(alpha=0) +
  facet_wrap(~Species) +
  theme_minimal() +
  theme(strip.text = element_text(face="italic"),
        axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(size=14,face="bold")) +
  labs(y="Shannon diversity")
ggsave("./Output/Figs/Shannon_Diversity_genus-glom_by_Structure_and_Species.png",
       dpi=300,height = 6,width = 6)

plot_richness(ps_genus, x="Location", measures=c("Shannon")) +
  geom_boxplot(alpha=0) +
  facet_wrap(~Species) +
  theme_minimal() +
  theme(strip.text = element_text(face="italic"),
        axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(size=14,face="bold")) +
  labs(y="Shannon diversity")
ggsave("./Output/Figs/Shannon_Diversity_genus-glom_by_Location.png",
       dpi=300,height = 6,width = 6)

# Model alpha diversity ####
meta <- microbiome::meta(ps_genus)                 
meta$Shannon <- vegan::diversity(otu_table(ps_genus),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps_genus))

shannon_mod <- lmerTest::lmer(data = meta,
                              formula = Shannon ~ Species * Structure + (1|Location))
summary(shannon_mod)
anova(shannon_mod)

sink("./Output/Stats/Shannon_Diversity_Model.txt")
summary(shannon_mod)
anova(shannon_mod)
sink(NULL)





