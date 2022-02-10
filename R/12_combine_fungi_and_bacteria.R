##  ###################################################  ##
##  Combine Bacterial and fungal ps objects

##  Author: Geoff Zahn - September, 2021                 ##
##                                                       ##

# load packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")


# load data objects
fung_2019 <- readRDS("./Output/Lee_2019_noncontam_ps_object.RDS")
fung_2020 <- readRDS("./Output/Lee_2020_noncontam_ps_object.RDS")
bact <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")

# pull "Structure" and "Location" from SampleID in fungal objects
meta_2019 <- fung_2019@sam_data %>% microbiome::meta()
meta_2020 <- fung_2020@sam_data %>% microbiome::meta()
bact_meta <- bact %>% microbiome::meta()
bact_meta$Location %>% unique()
bact_meta$Structure %>% unique()

meta_2020 <- meta_2020 %>% 
  mutate(Location = case_when(str_detect(SampleID, pattern = "_CJ_|CJawa") ~ "Chek Jawa",
                              str_detect(SampleID, pattern = "Ti_") ~ "Tioman",
                              str_detect(SampleID, pattern = "Kran_") ~ "Kranji",
                              str_detect(SampleID, pattern = "LK_") ~ "Langkawi",
                              str_detect(SampleID, pattern = "Mer_") ~ "Merang",
                              str_detect(SampleID, pattern = "Sem_") ~ "Semakau",
                              str_detect(SampleID, pattern = "TBali_") ~ "Tok Bali",
                              str_detect(SampleID, pattern = "Red_") ~ "Redang",
                              str_detect(SampleID, pattern = "PD_") ~ "Port Dickson",
                              str_detect(SampleID, pattern = "BLANK") ~ "NA")) %>% 
  mutate(Structure = case_when(str_detect(SampleID, pattern = "Fr_") ~ "Fruit",
                               str_detect(SampleID, pattern = "Le_") ~ "Leaf",
                               str_detect(SampleID, pattern = "Pn_") ~ "Pneumatophore",
                               str_detect(SampleID, pattern = "So_") ~ "Sediment",
                               TRUE ~ "Blank")) %>% 
  mutate(Study = "Lee_2020")

meta_2019 <- meta_2019 %>% 
  mutate(Location = case_when(str_detect(SampleID, pattern = "_CJ_|CJawa") ~ "Chek Jawa",
                              str_detect(SampleID, pattern = "Ti_") ~ "Tioman",
                              str_detect(SampleID, pattern = "Kran_") ~ "Kranji",
                              str_detect(SampleID, pattern = "LK_") ~ "Langkawi",
                              str_detect(SampleID, pattern = "Mer_") ~ "Merang",
                              str_detect(SampleID, pattern = "Sem_") ~ "Semakau",
                              str_detect(SampleID, pattern = "TBali_") ~ "Tok Bali",
                              str_detect(SampleID, pattern = "Red_") ~ "Redang",
                              str_detect(SampleID, pattern = "PD_") ~ "Port Dickson",
                              str_detect(SampleID, pattern = "BLANK") ~ "NA")) %>% 
  mutate(Structure = case_when(str_detect(SampleID, pattern = "Fr_") ~ "Fruit",
                               str_detect(SampleID, pattern = "Le_") ~ "Leaf",
                               str_detect(SampleID, pattern = "Pn_") ~ "Pneumatophore",
                               str_detect(SampleID, pattern = "So_") ~ "Sediment",
                               TRUE ~ "Blank")) %>% 
  mutate(Study = "Lee_2019")

# rename fields in bact
names(bact@sam_data)[3] <- "Host"
bact@sam_data$Study = "Bacteria"

# reassign metadata
fung_2019@sam_data <- sample_data(meta_2019)
fung_2020@sam_data <- sample_data(meta_2020)

# clean fungal ps objects
fung_2019 <- fung_2019 %>% 
  subset_taxa(Kingdom == "k__Fungi") %>% 
  subset_samples(sample_sums(fung_2019) > 0)

fung_2020 <- fung_2020 %>% 
  subset_taxa(Kingdom == "k__Fungi") %>% 
  subset_samples(sample_sums(fung_2020) > 0)

bact <- bact %>% 
  subset_taxa(taxa_sums(bact) > 0) %>% 
  subset_samples(sample_sums(bact) > 0)

# rename fungal taxa
tax_table(fung_2019)[,1] <- tax_table(fung_2019)[,1] %>% str_remove_all("k__")
tax_table(fung_2019)[,2] <- tax_table(fung_2019)[,2] %>% str_remove_all("p__")
tax_table(fung_2019)[,3] <- tax_table(fung_2019)[,3] %>% str_remove_all("c__")
tax_table(fung_2019)[,4] <- tax_table(fung_2019)[,4] %>% str_remove_all("o__")
tax_table(fung_2019)[,5] <- tax_table(fung_2019)[,5] %>% str_remove_all("f__")
tax_table(fung_2019)[,6] <- tax_table(fung_2019)[,6] %>% str_remove_all("g__")
tax_table(fung_2019)[,7] <- tax_table(fung_2019)[,7] %>% str_remove_all("s__")

tax_table(fung_2020)[,1] <- tax_table(fung_2020)[,1] %>% str_remove_all("k__")
tax_table(fung_2020)[,2] <- tax_table(fung_2020)[,2] %>% str_remove_all("p__")
tax_table(fung_2020)[,3] <- tax_table(fung_2020)[,3] %>% str_remove_all("c__")
tax_table(fung_2020)[,4] <- tax_table(fung_2020)[,4] %>% str_remove_all("o__")
tax_table(fung_2020)[,5] <- tax_table(fung_2020)[,5] %>% str_remove_all("f__")
tax_table(fung_2020)[,6] <- tax_table(fung_2020)[,6] %>% str_remove_all("g__")
tax_table(fung_2020)[,7] <- tax_table(fung_2020)[,7] %>% str_remove_all("s__")


tax_table(bact)[,1] %>% table
tax_table(fung_2019)[,1] %>% table
tax_table(fung_2020)[,1] %>% table

# inspect
fung_2019
fung_2020
bact

# merge_taxa at species level
fung_2019 <- tax_glom(fung_2019,taxrank = "Species")
fung_2020 <- tax_glom(fung_2020,taxrank = "Species")
bact <- tax_glom(bact,taxrank = "Genus")


# rename samples
fung_2020 <- fung_2020 %>% 
  subset_samples(!duplicated(fung_2020@sam_data$SampleID))

sample_names(fung_2019) <- fung_2019@sam_data$SampleID
sample_names(fung_2020) <- fung_2020@sam_data$SampleID

fung_2019 <- fung_2019 %>% subset_samples(sample_names(fung_2019) %in% sample_names(bact))
fung_2020 <- fung_2020 %>% subset_samples(sample_names(fung_2020) %in% sample_names(bact))


# cobmine and clean up environment
full <- merge_phyloseq(fung_2019,fung_2020)
full <- full %>% tax_glom(taxrank = "Species")




full2 <- merge_phyloseq(full,otu_table(bact),tax_table(bact),sample_data(bact))
tax_table(full2)[,1] %>% 
  table()


sample_names(bact)
# rm(fung_2019);rm(fung_2020);rm(bact)

# full <- full %>% 
#   subset_taxa(taxa_sums(full) > 0) %>% 
#   subset_samples(sample_sums(full) > 0)

# add microbe designation (bact vs fung)
full2@sam_data$Microbe <- ifelse(full2@sam_data$Study == "Bacteria",yes = "Bacteria",no="Fungi")

# save full ps object
saveRDS(full2,"./Output/bact_and_fungi_clean_ps_object")


