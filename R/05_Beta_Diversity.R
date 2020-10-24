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
##  purrr v 0.3.4                                        ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585")

# Load ps object glom by genus, and clean up a bit
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")
ps_genus <- ps_genus %>% subset_taxa(taxa_sums(ps_genus)>0)
ps_genus <- ps_genus %>% subset_samples(sample_sums(ps_genus)>0)
ps_genus <- ps_genus %>% subset_samples(Structure != "Blank")
saveRDS(ps_genus,"./Output/full_ps_object_w_tree_genus-glom.RDS")


# MANTEL TEST ####
spatial.dist = vegdist(as.matrix(cbind(ps_genus@sam_data$Lat, ps_genus@sam_data$Lon)),method = "bray")
comm.dist = vegdist(as.matrix(ps_genus@otu_table),method = "bray")
man <- vegan::mantel(spatial.dist,comm.dist)
man
sink("./Output/Stats/Mantel_test_summary.txt")
print("Spatial vs community distance...")
man
sink(NULL)


# Beta-diversity distances and ordinations ####
set.seed(123)
unifrac.dist <- UniFrac(ps_genus,weighted = TRUE,normalized = TRUE,parallel = TRUE)

glimpse(sample_data(ps_genus))
ps_genus@sam_data$Structure %>% unique()
ordu <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate("PCoA","unifrac", weighted=TRUE)

plot_ordination(ps_genus, ordu, color="Structure", shape="Species") +
  geom_point(size=3,alpha=.5) + 
  scale_color_manual(values=pal) +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance") +
  theme_minimal()
  
ggsave("./Output/Figs/W-Unifrac_Ordination_Plot_by_Structure.png",dpi=300)


plot_ordination(ps_genus, ordu, color="Location", shape="Species") +
  geom_point(size=3,alpha=.5) + 
  scale_color_viridis_d() +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance") +
  theme_minimal()

ggsave("./Output/Figs/W-Unifrac_Ordination_Plot_by_Location.png",dpi=300)



# beta-dispersion ####
w <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  otu_table() %>% 
  betadiver("w")
w.disper <- betadisper(w,group = meta(ps_genus)$Structure)
plot(w.disper,main = "Beta-Dispersion")

png("./Output/Figs/Beta-Dispersion_by_Structure.png")
plot(w.disper,main = "Beta-Dispersion")
dev.off()

w.disper <- betadisper(w,group = meta(ps_genus)$Species)
plot(w.disper,main = "Beta-Dispersion")

png("./Output/Figs/Beta-Dispersion_by_Species.png")
plot(w.disper,main = "Beta-Dispersion")
dev.off()




# PermANOVA ####
ps_ra <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)})

set.seed(123)
permanova <- vegan::adonis(otu_table(ps_ra) ~ ps_ra@sam_data$Structure * ps_ra@sam_data$Species * ps_ra@sam_data$Location)

sink("./Output/Stats/PermANOVA_test_summary.txt")
permanova
sink(NULL)

perm_df <- tidy(permanova$aov.tab)
perm_df$term <- str_remove_all(perm_df$term,"ps_ra@sam_data\\$")
perm_df %>% 
  arrange(desc(R2)) %>% 
  filter(!term %in% c("Total","Residuals")) %>% 
  mutate(across(where(is.numeric),function(x){round(x,3)})) %>% 
  write_csv("./Output/Stats/PermANOVA_model_terms_tidy.csv")


# What genera are driving differences? ####

simp <- simper(otu_table(ps_ra),group = ps_ra@sam_data$Structure)
saveRDS(simp,"./Output/simper_structure.RDS")
map(simp,"overall")
summary(simp)

# Indicator species? ####