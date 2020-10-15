##  ###################################################  ##
##  Build 16S phylogeny and add tree to full ps object   ##
##  This script processes the merged runs                ##
##                                                       ##
##  Author: Geoff Zahn - October 12, 2020                ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  tidyverse v 1.3.0                                    ##
##  phyloseq v 1.32.0                                    ##
##  phangorn                                             ##
##  vegan v 2.5.6                                        ##
##  phangorn v 2.5.5                                     ##
##  DECIPHER v 2.6.1                                     ##
##  ape v 5.4.1                                          ##
##  seqinr 4.2.4                                         ##
##                                                       ##
##  ###################################################  ##


# Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
# library(msa); packageVersion("msa")
library(DECIPHER); packageVersion("DECIPHER")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
theme_set(theme_bw())



# Read in phyloseq object ####
ps <- readRDS("./Output/")

# grab sequences as a DNAStringSet object
seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree
seqs_SS <- DNAStringSet(seqs)

# Multiple sequence alignment  ####
alignment <- AlignSeqs(seqs_SS,refinements = 3,processors = NULL,verbose = TRUE)
saveRDS(alignment,"./Output/Trees/16S_dna_alignment_decipher.RDS")

alignment_staggered <- StaggerAlignment(alignment)
saveRDS(alignment_staggered,"./Output/Trees/16S_dna_alignment_decipher_staggered.RDS")

# Convert to phangorn format
phang.align = as.phyDat(alignment_staggered, type = "DNA")
# write.phyDat(phang.align,"./output/trees/16S_dna_alignment_muscle.nex",format="nexus")

# distance max likelihood ####
dm <- dist.ml(phang.align)

#save
saveRDS(dm,"./Output/Trees/16S_ML_Distance.RDS")

# Initial neighbor-joining tree ####
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeNJ$tip.label <- seqs

#save
saveRDS(treeNJ, "./Output/Trees/16S_treeNJ.RDS")

# Estimate model parameters ####
fit = pml(treeNJ, data=phang.align)

#save
saveRDS(fit,"./Output/Trees/16S_fit_treeNJ.RDS")

# Likelihood of tree ####
fitJC <- optim.pml(fit, TRUE)

# save
saveRDS(fitJC, "./Output/Trees/16S_tree_fitJC.RDS") # This is the new tree using optim.pml
write.tree(fitJC$tree, file = "./Output/Trees/16S_tree_JC.nwk")

names(fitJC$tree$tip.label) <- NULL
identical(fitJC$tree$tip.label, taxa_names(ps))


# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
                phy_tree(fitJC$tree))

# Save updated phyloseq object with tree
saveRDS(ps2, "./output/full_ps_object_w_tree.RDS")
