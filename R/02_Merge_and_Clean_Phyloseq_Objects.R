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
##  vegan v 2.5.6                                        ##
##                                                       ##
##  ###################################################  ##

# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
source("./R/plot_bar2.R")

# Merge ####

# Grab ps objects (and add identifier)
ps06 <- readRDS("./Output/run6/noncontam_ps_object.RDS")
ps06@sam_data$IlluminaRun <- "1912KMI-0006"
ps07 <- readRDS("./Output/run7/noncontam_ps_object.RDS")
ps07@sam_data$IlluminaRun <- "1912KMI-0007"
ps12 <- readRDS("./Output/run12/noncontam_ps_object.RDS")
ps12@sam_data$IlluminaRun <- "1912KMI-0012"
ps13 <- readRDS("./Output/run13/noncontam_ps_object.RDS")
ps13@sam_data$IlluminaRun <- "1912KMI-0013"


# merge phyloseq objects from all 4 runs
full_ps <- merge_phyloseq(ps06,ps07,ps12,ps13)
saveRDS(full_ps, "./Output/full_ps_object_raw.RDS")

# plot of kingdom-level assignments
full_ps %>% merge_samples("Structure") %>% 
  plot_bar2(fill="Kingdom") +
  facet_wrap(~IlluminaRun) +
  theme_bw()
ggsave("./Output/Figs/Overall_Raw_Kingdom-Level_Assignments.png")

# Clean up taxa and any empty samples ####

# remove taxa with reads longer than 300bp
full_ps <- subset_taxa(full_ps, 
                       full_ps %>% otu_table() %>% colnames() %>% nchar() == 300)

# remove non-bacteria
full_ps <- 
  full_ps %>% 
  subset_taxa(Kingdom == "Bacteria")

# Remove "Chloroplast" taxa
full_ps <- 
  full_ps %>% 
  subset_taxa(Order != "Chloroplast")

# Remove "Mitochondria" taxa
full_ps <- 
  full_ps %>% 
  subset_taxa(Family != "Mitochondria")


# remove samples with very few reads (<= 100)
full_ps <- 
  full_ps %>% 
  subset_samples(sample_sums(full_ps) > 100)

# remove empty taxa
full_ps <- 
  full_ps %>% 
  subset_taxa(taxa_sums(full_ps) > 0)

# Save cleaned-up phyloseq object ####
saveRDS(full_ps, "./Output/full_ps_object_cleaned.RDS")



# Summary plots ####

# Plot of taxon-level assignment efficiency 
ps_sp <- full_ps
phy <- !is.na(tax_table(ps_sp)[,2])
cla <- !is.na(tax_table(ps_sp)[,3])
ord <- !is.na(tax_table(ps_sp)[,4])
fam <- !is.na(tax_table(ps_sp)[,5])
gen <- !is.na(tax_table(ps_sp)[,6])
spp <- !is.na(tax_table(ps_sp)[,7])
assignments <- data.frame(Phylum=phy, Class=cla,Order=ord,Family=fam,Genus=gen,Species=spp)

assignments %>% pivot_longer(1:6) %>% mutate(name=factor(name,levels = c("Phylum","Class","Order","Family","Genus","Species"))) %>%
  ggplot(aes(x=name,fill=value)) + geom_bar() + scale_fill_manual(values=c("Gray","Black")) +
  labs(x="Taxonomic level",y="Count",fill="Unambiguous\nassignment")

ggsave("./Output/Figs/SILVA_Taxonomic_Assignment_Efficiency_at_Each_Taxonomic_Rank.png",dpi=300)
rm(phy,cla,ord,fam,gen,spp,assignments,ps_sp)


# Plots of sample sums and taxon sums
png("./Output/Figs/Overall_ESV_Count_Distribution.png")
plot(taxa_sums(full_ps),main = "ESV read count distribution",ylab="Raw ESV Count")
dev.off()

png("./Output/Figs/Overall_Sample_Sum_Distribution.png")
plot(sort(sample_sums(full_ps),decreasing = TRUE),main = "Sample ESV Sums",ylab="Sample Sums")
dev.off()

# Rarefaction curves
png("./Output/Figs/Overall_Rarefaction_Curve.png")
rarecurve(otu_table(full_ps),step = 300,label = FALSE)
dev.off()

# ESV richness distribution
png("./Output/Figs/Overall_ESV_Richness_Distribution.png")
plot(sort(specnumber(otu_table(full_ps)),decreasing = TRUE),
     main="ESV Richness by Sample",ylab="ESV Richness")
dev.off()

# alert that script is finished with beep
beepr::beep(sound=8)


