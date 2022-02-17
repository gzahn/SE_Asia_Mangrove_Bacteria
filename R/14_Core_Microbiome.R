##  ###################################################  ##
##  Core Microbiome Members                              ##
##                                                       ##
##  Author: Geoff Zahn - Feb 16, 2022                    ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  tidyverse v 1.3.0                                    ##
##  phyloseq v 1.32.0                                    ##
##  vegan v 2.5.6                                        ##
##  broom v 0.7.1                                        ##
##  patchwork v 1.0.1                                    ##
##  microbiome v 1.10.0                                  ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
source("./R/plot_bar2.R")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585")

# Load ps objects
ps <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")


# find detection thresholds
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
plot_core(prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")




# find core microbiome members overall
overall_core_taxa <- core_members(ps_genus,
                                  detection = .10,
                                  prevalence = .10,
                                  include.lowest = TRUE)
# subset to just those taxa
ps_core <- ps_genus %>% 
  subset_taxa(taxa_names(ps_genus) %in% overall_core_taxa)

# find core for each plant structure
structures <- ps_genus@sam_data$Structure %>% unique()
  
for(i in structures){
 ps_genus %>% 
    subset_samples(Structure == i) %>% 
    core_members(detection = .10,
                 prevalence = .20,
                 include.lowest = TRUE) %>% 
    assign(paste0(i,"_core_members"),value=.,envir = .GlobalEnv)
}

# combine them
all_core_members <- c(Leaf_core_members,Fruit_core_members,
                      Pneumatophore_core_members,Sediment_core_members) %>% unique()

ps_genus_all_core <- ps_genus %>% 
  subset_taxa(taxa_names(ps_genus) %in% all_core_members)

# rename taxa to family/genus
core_taxa_names <- corncob::otu_to_taxonomy(data=ps_genus_all_core,level = c("Family","Genus"),taxa_names(ps_genus_all_core))
taxa_names(ps_genus_all_core) <- core_taxa_names


psm <- ps_genus_all_core %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt()

structure <- ps_genus_all_core %>% 
  merge_samples("Structure",fun="sum")

structure@sam_data$Structure <- row.names(structure@sam_data)

psm <- psm %>%
  arrange(Structure,sample_Species) %>% 
  mutate(sample_Species = factor(sample_Species,levels = unique(sample_Species)),
         Structure = factor(Structure,levels = unique(Structure)),
         Host_Structure = (paste0(sample_Species,"_",Structure)))

psm$Abundance[is.na(psm$Abundance)] <- 0

psm$SampleID <- factor(psm$SampleID,levels=unique(psm$SampleID))

psm <- psm %>% 
  mutate(Str_color = case_when(Structure == "Leaf" ~ pal[1],
                               Structure == "Fruit" ~ pal[2],
                               Structure == "Pneumatophore" ~ pal[3],
                               Structure == "Sediment" ~ pal[4],))
psm %>%   
ggplot(aes(x=Sample,y=OTU,fill=Abundance)) +
  geom_tile() +
  facet_grid(cols = vars(Structure),
             scales = 'free') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face="bold",size=14),
        axis.text.y = element_text(size=5)) +
  scale_fill_viridis_c() +
  labs(y="Taxon",fill="Relative\nabundance")
ggsave("./Output/Figs/Core_Heatmap_by_Structure.png",dpi=400,height = 8,width = 12)

aa_heatmap <- psm %>% 
  filter(sample_Species == "Avicennia alba") %>% 
  ggplot(aes(x=Sample,y=OTU,fill=Abundance)) +
  geom_tile() +
  facet_grid(cols = vars(Structure),
             scales = 'free') +
  theme_minimal() +
  labs(y="Taxon",fill="Relative\nabundance",title = "Avicennia alba") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(face="bold",size=14),
        axis.text.y = element_text(size=5),
        plot.title = element_text(face="italic",size=14,hjust=.5),
        legend.position = 'none') +
  scale_fill_viridis_c() 

sa_heatmap <- psm %>% 
  filter(sample_Species == "Sonneratia alba") %>% 
  ggplot(aes(x=Sample,y=OTU,fill=Abundance)) +
  geom_tile() +
  facet_grid(cols = vars(Structure),
             scales = 'free') +
  theme_minimal() +
  labs(y="Taxon",fill="Relative\nabundance",title = "Sonneratia alba") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(face="bold",size=14),
        axis.text.y = element_blank(),
        plot.title = element_text(face="italic",size=14,hjust=.5),
        axis.title.y = element_blank()) +
  scale_fill_viridis_c() 

aa_heatmap + sa_heatmap
ggsave("./Output/Figs/Core_Heatmap_by_Structure_and_Host.png",dpi = 400,height = 8,width = 16)
