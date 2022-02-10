##  ###################################################  ##
##  Processing raw 16S seqs into a phyloseq object       ##
##  This script processes "Lee_2020"  PRJNA592423        ##
##                                                       ##
##  Sonneratia alba mycobiome                            ##
##                                                       ##
##  Author: Geoff Zahn - September 2021                  ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  dada2 v 1.16.0                                       ##
##  purrr v 0.3.4                                        ##
##  tidyverse v 1.3.0                                    ##
##  readxl v 1.3.1                                       ##
##  decontam v 1.8.0                                     ##
##  phyloseq v 1.32.0                                    ##
##                                                       ##
##  ###################################################  ##


# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(tidyverse); packageVersion("tidyverse")
library(readxl); packageVersion("readxl")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")

# Load metadata ####
Lee_2019_meta <- readxl::read_xlsx("./Data/Lee_2019/SRA_Metadata_Lee_2019.xlsx")
Lee_2020_meta <- readxl::read_xlsx("./Data/Lee_2020/SRA_Metadata_Lee_2020.xlsx")

bact_ps <- readRDS("Output/full_ps_object_cleaned.RDS") 



# tidy metadata values to match bacterial sample IDs
Lee_2019_meta$SampleID <- Lee_2019_meta$SampleID %>% str_replace_all("Soil","So")
Lee_2020_meta$SampleID <- Lee_2020_meta$SampleID %>% str_replace_all("Soil","So")
Lee_2020_meta$SampleID <- Lee_2020_meta$SampleID %>% str_replace_all("CJ","CJawa")
Lee_2020_meta$SampleID <- Lee_2020_meta$SampleID %>% str_replace_all("TB","TBali")
Lee_2020_meta$SampleID <- Lee_2020_meta$SampleID %>% str_replace_all("TI","Ti")
Lee_2020_meta$SampleID <- Lee_2020_meta$SampleID %>% str_remove_all("HO_")
Lee_2020_meta$SampleID <- Lee_2020_meta$SampleID %>% str_replace_all("LCK","LK")
bact <- bact_ps %>% 
  microbiome::meta()

bact$SampleID <- bact$SampleID %>% str_replace_all("LCK","LK")

bact_ps <- phyloseq(otu_table(bact_ps),
         tax_table(bact_ps),
         sample_data(bact))
bact$SampleID

lee_2020_sampleIDs <- Lee_2020_meta$SampleID %>%
  str_split("_")
lee_2020_sampleIDs <- lapply(lee_2020_sampleIDs, `length<-`, max(lengths(lee_2020_sampleIDs)))
Lee_2020_meta$SampleIDList <- lee_2020_sampleIDs

Lee_2020_meta <- Lee_2020_meta %>% 
  mutate(SampleID=case_when(str_detect(SampleID, pattern = "BLANK") ~ paste(lee_2020_sampleIDs %>% map_chr(1),
                                                                                lee_2020_sampleIDs %>% map_chr(2),
                                                                                lee_2020_sampleIDs %>% map_chr(3),
                                                                                sep = "_"),
                                TRUE ~ paste(lee_2020_sampleIDs %>% map_chr(2),
                                             lee_2020_sampleIDs %>% map_chr(1),
                                             lee_2020_sampleIDs %>% map_chr(3),
                                             lee_2020_sampleIDs %>% map_chr(4),
                                             sep = "_"))) 

# subset lee datasets to only those found in bacteria data set
Lee_2019_meta <- Lee_2019_meta %>% 
  filter(SampleID %in% bact$SampleID)
Lee_2020_meta <- Lee_2020_meta %>% 
  filter(SampleID %in% bact$SampleID)

# Add columns identifying PCR negatives
Lee_2019_meta$PCR_Negative <- FALSE
Lee_2019_meta$PCR_Negative[grep("BLANK",Lee_2019_meta$SampleID)] <- TRUE
Lee_2020_meta$PCR_Negative <- FALSE
Lee_2020_meta$PCR_Negative[grep("BLANK",Lee_2020_meta$SampleID)] <- TRUE


# generate file path columns
Lee_2019_meta$Fwd_Path <- paste0("./Data/Lee_2019/fastq/",Lee_2019_meta$Run,"_pass_1.fastq.gz")
Lee_2019_meta$Rev_Path <- paste0("./Data/Lee_2020/fastq/",Lee_2019_meta$Run,"_pass_2.fastq.gz")

Lee_2020_meta$Fwd_Path <- paste0("./Data/Lee_2020/fastq/",Lee_2020_meta$Run,"_pass_1.fastq.gz")
Lee_2020_meta$Rev_Path <- paste0("./Data/Lee_2020/fastq/",Lee_2020_meta$Run,"_pass_2.fastq.gz")


# Save cleaned metadata files
write_csv(Lee_2019_meta,"./Data/Lee_2019/SRA_Metadata_Lee_2019_clean.csv")
write_csv(Lee_2020_meta,"./Data/Lee_2020/SRA_Metadata_Lee_2020_clean.csv")




