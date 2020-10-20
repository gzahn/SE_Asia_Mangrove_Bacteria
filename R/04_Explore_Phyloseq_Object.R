##  ###################################################  ##
##  Investigate and tidy up full phyloseq object         ##
##                                                       ##
##  Author: Geoff Zahn - October 12, 2020                ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  tidyverse v 1.3.0                                    ##
##  phyloseq v 1.32.0                                    ##
##  vegan v 2.5.6                                        ##
##  VennDiagram v 1.6.20                                 ##
##  patchwork v 1.0.1                                    ##
##                                                       ##
##  ###################################################  ##

# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(VennDiagram); packageVersion("VennDiagram")
library(patchwork); packageVersion("patchwork")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585")
# names(pal) <- c("f","l","p","s")

# Load ps object with tree ####
full_ps <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")

# Look at ESV overlap (raw reads) ####
# find ESVs shared between Illumina runs
# Grab ps objects (and add identifier)
ps06 <- readRDS("./Output/run6/noncontam_ps_object.RDS")
ps06@sam_data$IlluminaRun <- "1912KMI-0006"
ps07 <- readRDS("./Output/run7/noncontam_ps_object.RDS")
ps07@sam_data$IlluminaRun <- "1912KMI-0007"
ps12 <- readRDS("./Output/run12/noncontam_ps_object.RDS")
ps12@sam_data$IlluminaRun <- "1912KMI-0012"
ps13 <- readRDS("./Output/run13/noncontam_ps_object.RDS")
ps13@sam_data$IlluminaRun <- "1912KMI-0013"
ps06names <- ps06 %>% otu_table() %>% colnames()
ps07names <- ps07 %>% otu_table() %>% colnames()
ps12names <- ps12 %>% otu_table() %>% colnames()
ps13names <- ps13 %>% otu_table() %>% colnames()
sum(ps06names %in% ps07names)
sum(ps06names %in% ps12names)
sum(ps06names %in% ps13names)
sum(ps07names %in% ps12names)
sum(ps07names %in% ps13names)
sum(ps12names %in% ps13names)
# clean up
rm(list = c("ps06","ps07","ps12","ps13","ps06names","ps07names","ps12names","ps13names"))

# agglomerate taxa at genus level ####
full_ps_genus <- tax_glom(full_ps,"Genus")



# Look again at overlap (at genus level) and build venn diagram of overlap by IlluminaRun ####
A <- full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0006") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0006"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0007") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0007"))>0) %>% tax_table() %>% row.names()
C <- full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0012") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0012"))>0) %>% tax_table() %>% row.names()
D <- full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0013") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(IlluminaRun == "1912KMI-0013"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C,D))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n34 <- sum(full %in% unique(C) & full %in% unique(D))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))
n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

venn.plot <- draw.quad.venn(area1=length((A)),
               area2=length((B)),
               area3=length((C)),
               area4=length((D)),
               n12,n13,n14,n23,n24,n34,n123,n124,n134,n234,n1234,
               category = c("1912KMI-0006","1912KMI-0007","1912KMI-0012","1912KMI-0013"))

dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_IlluminaRun.png")
grid.draw(venn.plot)
dev.off()


# Same VennDiagram, but for Plant Parts ####

A <- full_ps_genus %>% subset_samples(Structure == "Fruit") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Structure == "Fruit"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Structure == "Leaf") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Structure == "Leaf"))>0) %>% tax_table() %>% row.names()
C <- full_ps_genus %>% subset_samples(Structure == "Pneumatophore") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Structure == "Pneumatophore"))>0) %>% tax_table() %>% row.names()
D <- full_ps_genus %>% subset_samples(Structure == "Sediment") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Structure == "Sediment"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C,D))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n34 <- sum(full %in% unique(C) & full %in% unique(D))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))
n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

venn.plot2 <- draw.quad.venn(area1=length((A)),
                            area2=length((B)),
                            area3=length((C)),
                            area4=length((D)),
                            n12,n13,n14,n23,n24,n34,n123,n124,n134,n234,n1234,
                            category = c("Fruit","Leaf","Pneumatophore","Sediment"))
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_Structure.png")
grid.draw(venn.plot2)
dev.off()

# Save genus-level-glom object ####

# output genus-level ps_object for convenience
saveRDS(full_ps_genus, "./Output/full_ps_object_w_tree_genus-glom.RDS")


# Explore phylogenetic tree ####
full_ps_genus@sam_data$Structure %>% unique()
# colored by phylum (blanks removed)
full_ps_genus %>% 
  subset_samples(Structure != "Blank") %>% 
  plot_tree(ladderize="left", color="Structure") +
    scale_color_viridis_d()

p1 <- full_ps_genus %>% 
  subset_samples(Structure == "Fruit") %>% 
  plot_tree(ladderize="left", color="Structure") +
  scale_color_manual(values = pal[1]) +
  theme(legend.position = "none") +
  ggtitle("Fruit")

p2 <- full_ps_genus %>% 
  subset_samples(Structure == "Leaf") %>% 
  plot_tree(ladderize="left", color="Structure") +
  scale_color_manual(values = pal[2]) +
  theme(legend.position = "none") +
  ggtitle("Leaf")

p3 <- full_ps_genus %>% 
  subset_samples(Structure == "Pneumatophore") %>% 
  plot_tree(ladderize="left", color="Structure") +
  scale_color_manual(values = pal[3]) +
  theme(legend.position = "none") +
  ggtitle("Pneumatophore")

p4 <- full_ps_genus %>% 
  subset_samples(Structure == "Sediment") %>% 
  plot_tree(ladderize="left", color="Structure") +
  scale_color_manual(values = pal[4]) +
  theme(legend.position = "none") +
  ggtitle("Sediment")

(p1+p2) / (p3+p4)
ggsave("./Output/Figs/Phylogenetic_Dispersion_by_Plant_Structure.png",dpi=300)

rm(p1);rm(p2);rm(p3);rm(p4)
