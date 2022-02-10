##  ###################################################  ##
##  Compare bacteria and fungi diversity patterns

##  Author: Geoff Zahn - September, 2021                 ##
##                                                       ##

# load packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(patchwork)
source("./R/plot_bar2.R")
source("./R/palettes.R")

# load data objects


full <- readRDS("./Output/bact_and_fungi_clean_ps_object.RDS")

full <- full %>% 
  subset_samples(Structure != "Blank")

# remove any lingering empty samples/taxa
full <- full %>% 
  subset_taxa(taxa_sums(full) > 0) %>% 
  subset_samples(sample_sums(full) > 0)

# create two subsets and clean them up
bact <- full %>% 
  subset_taxa(Kingdom == "Bacteria")
bact <- bact %>% subset_samples(sample_sums(bact) > 0)
fung <- full %>% 
  subset_taxa(Kingdom == "Fungi")
fung <- fung %>% subset_samples(sample_sums(fung) > 0)

# relative abundance transformation
full_ra <- transform_sample_counts(full,function(x){x/sum(x)})


# table of kingdom-level assignments
tax_table(full)[,1] %>% table()

# barplot
merge_samples(full,"Structure") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Kingdom") +
  theme_minimal() +
  scale_fill_manual(values=pal.discrete)
ggsave("./Output/Figs/bact_vs_fungi/Kingdom_Relative_Abundance_plant_Structure.png")



# alpha div 
full_ra@sam_data$Shannon <- estimate_richness(full_ra, measures = "Shannon")


full@sam_data$Microbe

# beta div
ord <- ordinate(full_ra,method = "PCoA")
plot_ordination(full_ra,ord,color="Structure") +
  facet_wrap(~Microbe) +
  theme_bw() +
  scale_color_manual(values=pal.discrete)
ggsave("./Output/Figs/bact_vs_fungi/PCoA_Plot_by_Study_and_Structure.png")

bact_ord <- bact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "PCoA")
fung_ord <- fung %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "PCoA")
p1 <- plot_ordination(bact,bact_ord,color="Structure") +
  theme_bw() +
  scale_color_manual(values=pal.discrete) +
  theme(legend.position = "none") + labs(title = "Bacteria")
p2 <- plot_ordination(fung,fung_ord,color="Structure") +
  theme_bw() + labs(title = "Fungi") +
  scale_color_manual(values=pal.discrete)
p1 + p2


# models
asv <- full_ra@otu_table %>% as.matrix() %>% as.data.frame()
structure <- full_ra@sam_data$Structure
location <- full_ra@sam_data$Location
domain <- full_ra@sam_data$Microbe

permanova <- vegan::adonis(asv ~ structure * location, strata = full_ra@sam_data$Study)
permanova
