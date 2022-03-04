##  ###################################################  ##
##  Alpha-diversity measures - btwn site and structure   ##
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
source("./R/plot_bar2.R")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585")

# Load ps objects
ps <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

# Model alpha diversity ####
meta <- microbiome::meta(ps_genus)                 
meta$Shannon <- vegan::diversity(otu_table(ps_genus),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps_genus))
# Set "Sediment" to first factor level of Structure
meta$Structure <- factor(meta$Structure,levels = c("Sediment","Fruit","Leaf","Pneumatophore"))

# add to ps object
ps_genus@sam_data$Richness <- meta$Richness
ps_genus@sam_data$Shannon <- meta$Shannon

# lme models
shannon_mod <- lmerTest::lmer(data = meta,
                              formula = Shannon ~ Species * Structure * (1|Location))
richness_mod <- lmerTest::lmer(data = meta,
                               formula = Richness ~ Species * Structure * (1|Location))
# send to file
sink("./Output/Stats/Shannon_Diversity_Model.txt")
summary(shannon_mod)
anova(shannon_mod)
sink(NULL)

sink("./Output/Stats/Richness_Diversity_Model.txt")
summary(richness_mod)
anova(richness_mod)
sink(NULL)


# Stacked phylum-level barcharts (horizontal) by host and structure ####

# merge samples
ps_genus@sam_data %>% names()
ps_genus@otu_table %>% rowSums()
newmergevar <- paste(ps_genus@sam_data$Species,
                     ps_genus@sam_data$Location,
                     ps_genus@sam_data$Structure,
                     sep = "_")
ps_genus@sam_data$newmergevar <- newmergevar

psm_ra <- ps_genus %>%   
  merge_samples(newmergevar,fun = "sum") %>% 
  transform_sample_counts(function(x){x/sum(x)})

# repair metadata
sp <- psm_ra@sam_data %>% row.names() %>% str_split("_") %>% map_chr(1)
lo <- psm_ra@sam_data %>% row.names() %>% str_split("_") %>% map_chr(2)
st <- psm_ra@sam_data %>% row.names() %>% str_split("_") %>% map_chr(3)
psm_ra@sam_data$Species <- sp
psm_ra@sam_data$Location <- lo
psm_ra@sam_data$Structure <- st

aa_plot <- psm_ra %>% 
  subset_samples(Species == "Avicennia alba") %>% 
  plot_bar2(x="Location",fill="Phylum",facet_grid = ~Structure) +
  coord_flip() +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=8,margin = margin(c(5,10,5,10))),
        axis.text.x = element_text(face="bold",size=8,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face="bold",size=10),
        axis.title = element_text(face="bold",size=12),
        legend.title = element_text(face="bold",size=12),
        legend.text = element_text(size=10),
        plot.title = element_text(face="bold.italic",size=14,hjust=.5),
        legend.position = "none",
        plot.margin = margin(r=20)) +
  scale_fill_viridis_d() +
  scale_y_continuous(breaks=c(0,.5,1)) +
  labs(y="\nRelative abundance",title = "Avicennia alba")

sa_plot <- psm_ra %>% 
  subset_samples(Species == "Sonneratia alba") %>% 
  plot_bar2(x="Location",fill="Phylum",facet_grid = ~Structure) +
  coord_flip() +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=8,margin = margin(c(5,-20,5,-20))),
        axis.text.x = element_text(face="bold",size=8,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_blank(),
        legend.title = element_text(face="bold",size=12),
        legend.text = element_text(size=10),
        plot.title = element_text(face="bold.italic",size=14,hjust=.5)) +
  scale_fill_viridis_d() +
  scale_y_continuous(breaks=c(0,.5,1)) +
  labs(y="\nRelative abundance",title = "Sonneratia alba")

aa_plot + sa_plot
ggsave("./Output/Figs/horizontal_bar_charts.png",dpi=400,width = 16,height = 8)

# Plot alpha diversity ####
plot_richness(ps_genus, x="Structure", measures=c("Shannon")) +
  geom_boxplot(alpha=0.5,aes(fill=Structure)) +
  facet_wrap(~Species) +
  theme_bw() +
  theme(strip.text = element_text(face="bold.italic",size=12),
        axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(size=14,face="bold"),
        strip.background = element_blank(),
        axis.text = element_text(face="bold")) +
  labs(y="Shannon diversity") +
  scale_fill_manual(values=pal)
ggsave("./Output/Figs/Shannon_Diversity_genus-glom_by_Structure_and_Species.png",
       dpi=300,height = 6,width = 6)

plot_richness(ps_genus, x="Location", measures=c("Shannon")) +
  geom_boxplot(alpha=0.5,aes(fill=Structure)) +
  facet_wrap(~Species) +
  theme_bw() +
  theme(strip.text = element_text(face="bold.italic",size=12),
        axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(size=14,face="bold"),
        strip.background = element_blank(),
        axis.text = element_text(face="bold")) +
  labs(y="Shannon diversity") +
  scale_fill_manual(values=pal)
ggsave("./Output/Figs/Shannon_Diversity_genus-glom_by_Location.png",
       dpi=300,height = 6,width = 6)





# Plot phylum diversity (barplots) ####
# merge by species and structure
ps_genus@sam_data$mergevar <- paste0(ps_genus@sam_data$Species,"_",ps_genus@sam_data$Structure)
psm <- ps_genus %>%  merge_samples(group = "mergevar")
# repair metadata
psm@sam_data$Species <- rownames(psm@sam_data) %>% str_split("_") %>% map_chr(1)
psm@sam_data$Structure <- rownames(psm@sam_data) %>% str_split("_") %>% map_chr(2)
psm@sam_data$mergevar <- rownames(psm@sam_data)
psm@sam_data$Tree <- psm@sam_data$Species

# plot relative abundance barcharts
Aa <- psm %>% 
  subset_samples(Species == "Avicennia alba") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Phylum",x="Structure",title = "Avicennia alba") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size=12,face="bold",angle=60,hjust=1),
        axis.title = element_text(size=14,face="bold"),
        legend.title = element_text(size=12,face="bold"),
        title = element_text(size=14,face="bold.italic"),
        legend.position = "none",
        axis.title.x = element_blank()) +
  labs(y="Relative abundance")



Sa <- psm %>% 
  subset_samples(Species == "Sonneratia alba") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Phylum",x="Structure",title = "Sonneratia alba") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size=12,face="bold",angle=60,hjust=1),
        axis.title = element_text(size=14,face="bold"),
        legend.title = element_text(size=12,face="bold"),
        title = element_text(size=14,face="bold.italic"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# combine plots and save
Aa + Sa
ggsave("./Output/Figs/Phylum_Stacked_Barchart_by_Structure.png",dpi=300,width = 12,height = 6)



# Most abundant phyla boxplot ####

# For each phylum, subset to all taxa in that phylum, save otu_table in named list 
phylum_names <- ps_genus %>% tax_table() %>% as.data.frame() %>% select(2)
ps_genus_ra <- transform_sample_counts(ps_genus, function(x){x/sum(x)})
x <- 1
newlist <- list()
for(i in unique(phylum_names$Phylum)){
  tax <- subset_taxa(ps_genus_ra, ps_genus_ra@tax_table[,2]  ==  i)
  OTU = as(otu_table(tax), "matrix")
  if(taxa_are_rows(tax)){OTU <- t(OTU)}
  OTUdf = as.data.frame(OTU)
  newlist[i] <- OTUdf
}
# build data frame
phylum_ra_df <- newlist %>% as.data.frame()
phylum_ra_df$Structure <- ps_genus_ra@sam_data$Structure
phylum_ra_df$Species <- ps_genus_ra@sam_data$Species
phylum_ra_df$Location <- ps_genus_ra@sam_data$Location

# tidy for plotting
tax_abund_order <- phylum_ra_df %>% select(unique(phylum_names$Phylum)) %>% map_dbl(mean) %>% sort(decreasing = TRUE) %>% names() # for all samples
tax_abund_order2 <- phylum_ra_df %>% filter(Structure != "Sediment") %>% select(tax_abund_order) %>% colMeans() %>% sort(decreasing = TRUE) %>% names() # for non-sediment samples
phylum_ra_df_long <- pivot_longer(phylum_ra_df, unique(phylum_names$Phylum),names_to="Phylum",values_to = "Relative_Abundance")
# remove zeroes for plotting  
phylum_ra_df_long$Relative_Abundance[phylum_ra_df_long$Relative_Abundance == 0] <- NA
# plot top 6 phyla
phylum_ra_df_long %>% 
  filter(Phylum %in% tax_abund_order2[1:6]) %>% 
  ggplot(aes(x=Structure,y=Relative_Abundance,fill=Structure)) +
    geom_boxplot() +
    facet_wrap(~Phylum,scales = "free_y",) +
    theme_bw() +
    scale_fill_manual(values = pal) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12,face="bold"),
          legend.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=10,face="bold"),
          axis.title = element_text(size=14,face="bold"),
          axis.text.x = element_text(angle=60,hjust=1,size=10)) +
  labs(y="Relative abundance",caption = "Phyla with highest relative abundance in plants, not sediment")
ggsave("./Output/Figs/Most_abundant_phyla_boxplot_by_Structure.png",dpi=300,height = 8,width = 12)



