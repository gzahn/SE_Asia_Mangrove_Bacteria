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
##  corncob v 0.1.0                                      ##
##  indicspecies v 1.7.9                                 ##
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
library(corncob); packageVersion("corncob")
library(indicspecies); packageVersion("indicspecies")

source("./R/bbdml_helper.R")

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
set.seed(123)
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
map_dbl(simp,"overall")
summary(simp)

# Indicator species analysis for different plant structures ####

# set up clusters by Structure
clusters <- rownames(ps_ra@otu_table) %>% str_split("_") %>% map_chr(3)
# convert to matrix
OTU = as(otu_table(ps_ra), "matrix")

# find indicator species for each plant part
flower <- indicators(OTU,cluster = clusters,group = "Fr",verbose = TRUE,func = "IndVal.g",At=.5,Bt=.2)
saveRDS(flower,"./Output/indic_flower.RDS")
leaf <- indicators(OTU,cluster = clusters,group = "Le",verbose = TRUE,func = "IndVal.g",At=.5,Bt=.2)
saveRDS(leaf,"./Output/indic_leaf.RDS")
pneumatophore <- indicators(OTU,cluster = clusters,group = "Pn",verbose = TRUE,func = "IndVal.g",At=.5,Bt=.2)
saveRDS(pneumatophore,"./Output/indic_pneumatophore.RDS")

leaf_indics <- as.data.frame(ps_ra@tax_table[leaf$finalsplist,2:6])
flower_indics <- as.data.frame(ps_ra@tax_table[flower$finalsplist,2:6])
pneumatophore_indics <- as.data.frame(ps_ra@tax_table[pneumatophore$finalsplist,2:6])

# CoRNCOB differential abundance analysis ####

set.seed(123)
da_analysis <- differentialTest(formula = ~ Structure, #abundance
                                phi.formula = ~ Structure, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

length(da_analysis$significant_models)

set.seed(123)
bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps_genus,
                           mu_predictor = "Structure",
                           phi_predictor = "Structure",
                           taxlevels = 2:6)
length(bbdml_obj)
saveRDS(bbdml_obj,"./Output/bbdml_Structure.RDS")

# find corncob diffabund taxa that were also found by indicspecies
indic_species_list <- unique(c(leaf_indics$Genus,flower_indics$Genus,pneumatophore_indics$Genus))
found_by_both <- which(str_split(names(bbdml_obj),"_") %>% map_chr(5) %in% indic_species_list)

# generate bbdml plots
plot_multi_bbdml(bbdml_list = bbdml_obj[found_by_both],
                   color = "Structure",
                   pointsize = 3)

#  Compose bbdml plots ####
scaleFUN <- function(x) sprintf("%.2f", x)

bbdml_plot_1 <- bbdml_plot_1 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "none",
        legend.text = element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal)

bbdml_plot_2 <- bbdml_plot_2 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "none",
        legend.text = element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal)
bbdml_plot_3 <- bbdml_plot_3 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "none",
        legend.text = element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal)
bbdml_plot_4 <- bbdml_plot_4 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "none",
        legend.text = element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal)
bbdml_plot_5 <- bbdml_plot_5 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "none",
        legend.text = element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal)
bbdml_plot_6 <- bbdml_plot_6 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "none",
        legend.text = element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal)
bbdml_plot_7 <- bbdml_plot_7 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "none",
        legend.text = element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal)
bbdml_plot_8 <- bbdml_plot_8 + scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        legend.title = element_text(size=18,face="bold"), legend.position = "bottom", legend.key.size = unit(.5,units = "in"),
        legend.text = element_text(size=14,face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold")) +
  scale_color_manual(values=pal) + labs(color="Sample source: ")

bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3 / bbdml_plot_4 / bbdml_plot_5 / bbdml_plot_6 / bbdml_plot_7 / bbdml_plot_8
ggsave("./Output/Figs/Differentially-Abundant_Taxa_bbdml_plots.png",dpi=300,height = 12,width = 10)

# plot DA analysis
plot(da_analysis)



# Which taxa are found in mangroves but NOT in sediment? ####

# for entire study...
ps <- ps %>% subset_samples(Structure != "Blank")

# add t/f column for whether sample is plant or not
plant <- ps %>% 
  meta %>% 
  mutate(Plant = case_when(Structure == "Sediment" ~ FALSE,
                           Structure != "Sediment" ~ TRUE))
ps@sam_data$Plant <- plant$Plant

# separate ps objects
sediment <- ps %>% 
  subset_samples(Structure == "Sediment")
sediment <- sediment %>% 
  subset_taxa(taxa_sums(sediment) > 0)
  
plant <- ps %>% 
  subset_samples(Structure != "Sediment")
plant <- plant %>% 
  subset_taxa(taxa_sums(plant) > 0)


# N sediment taxa found in plants
taxa_names(sediment) %in% taxa_names(plant) %>% sum()
ntaxa(plant)
taxa_names(plant) %in% taxa_names(sediment) %>% sum()



# for each sampling location, individually...
locations <- meta(plant)$Location %>% unique()

taxa_by_loc_plant <- list()
plant@sam_data
for(i in 1:length(locations)){
  x <- plant %>% 
    subset_samples(Location == locations[i]) 
  taxa_by_loc_plant[[i]] <- x %>% 
    subset_taxa(taxa_sums(x) > 0) %>% 
    taxa_names()
}

taxa_by_loc_sediment <- list()

for(i in 1:length(locations)){
  x <- sediment %>% 
    subset_samples(Location == locations[i]) 
  taxa_by_loc_sediment[[i]] <- x %>% 
    subset_taxa(taxa_sums(x) > 0) %>% 
    taxa_names()
}



# write function to compare taxa names by site between plant and sediment lists

num_taxa_by_loc <- list()
for(i in 1:length(locations)){
num_taxa_by_loc[[i]] <- taxa_by_loc_plant[[i]] %in% taxa_by_loc_sediment[[i]]
}


taxa_numbers <- data.frame(
  plant = lapply(taxa_by_loc_plant, length) %>% unlist(),  
  sediment = lapply(taxa_by_loc_sediment, length) %>% unlist(),
  shared_w_sediment = lapply(num_taxa_by_loc,sum) %>% unlist(),
  row.names = locations
)
taxa_numbers
write_csv(taxa_numbers, "./Output/Stats/taxa_numbers_plant_vs_sediment_by_location.csv")

# every taxon seen in a plant was found in sediment from the same sites!?
ps@sam_data$Location %>% unique

leaf <- ps %>% subset_samples(Structure == "Leaf" & Location == "Chek Jawa") %>% 
  subset_taxa(taxa_sums(ps %>% subset_samples(Structure == "Leaf" & Location == "Chek Jawa"))>0) %>% tax_table() %>% row.names()
  
sediment <- ps %>% subset_samples(Structure == "Sediment" & Location == "Chek Jawa") %>% 
  subset_taxa(taxa_sums(ps %>% subset_samples(Structure == "Sediment" & Location == "Chek Jawa"))>0) %>% tax_table() %>% row.names()

leaf %in% sediment %>% sum
