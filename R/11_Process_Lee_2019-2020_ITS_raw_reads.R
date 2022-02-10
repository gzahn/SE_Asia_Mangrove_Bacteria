##  ###################################################  ##
##  Processing raw ITS data from Lee et al (2019, 2020)

##  Author: Geoff Zahn - September, 2021                 ##
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

# load metadata
lee_2019_meta <- read_csv("./Data/Lee_2019/SRA_Metadata_Lee_2019_clean.csv")
lee_2020_meta <- read_csv("./Data/Lee_2020/SRA_Metadata_Lee_2020_clean.csv")


path_2019 <- "./Data/Lee_2019/fastq"
path_2020 <- "./Data/Lee_2020/fastq"

# Find raw fastq files and prepare workspace ####
fqs_2019 <- list.files(path_2019, full.names = TRUE,recursive = TRUE,pattern = "_ITS.fastq.gz$")
fqs_2020 <- list.files(path_2019, full.names = TRUE,recursive = TRUE,pattern = "_ITS.fastq.gz$")

# Get Sample Names
sample.names_2019 <- str_remove(basename(fqs_2019),"_ITS.fastq.gz")
sample.names_2020 <- str_remove(basename(fqs_2020),"_ITS.fastq.gz")


# Peek at quality profiles
plotQualityProfile(fqs_2019[c(1,30)]) # fwd reads
plotQualityProfile(fqs_2020[c(1,30)]) # rev reads

# Make filtered outfile names
filt_2019 <- file.path(path_2019, "filtered", paste0(sample.names_2019, "_filt.fastq.gz"))
filt_2020 <- file.path(path_2020, "filtered", paste0(sample.names_2020, "_filt.fastq.gz"))

# make new directory for filtered files
if(!dir.exists(file.path(path_2019,"filtered"))){
  dir.create(file.path(path_2019,"filtered"))
}
if(!dir.exists(file.path(path_2020,"filtered"))){
  dir.create(file.path(path_2020,"filtered"))
}

# check for duplicated sample names
sum(duplicated(sample.names_2019))
sum(duplicated(sample.names_2020))



# Filter and trim ####

# cut fwd reads at 300 and rev reads at 200
out_2019 <- filterAndTrim(fqs_2019, filt_2019,
                          maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE,verbose = TRUE)
out_2020 <- filterAndTrim(fqs_2020, filt_2020,
                          maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE,verbose = TRUE)

saveRDS(out_2019,"./Output/Lee_2019_out.RDS") # save object
saveRDS(out_2020,"./Output/Lee_2020_out.RDS") # save object


filtpath_2019 <- file.path(path_2019,"filtered")
filtpath_2020 <- file.path(path_2020,"filtered")

# reassign filts for any potentially lost samples
filt_2019 <- list.files(filtpath_2019, pattern = "filt", full.names = TRUE)
filt_2020 <- list.files(filtpath_2020, pattern = "filt", full.names = TRUE)

# Get Sample Names (again, just in case)
sample.names_2019 <- str_remove(basename(filt_2019),"_filt.fastq.gz")
sample.names_2020 <- str_remove(basename(filt_2020),"_filt.fastq.gz")

# Learn error rates ####
err_2019 <- learnErrors(filt_2019, multithread=TRUE)
err_2020 <- learnErrors(filt_2020, multithread=TRUE)
saveRDS(err_2019,"./Output/err_Lee_2019.RDS") # save object
saveRDS(err_2019,"./Output/err_Lee_2020.RDS") # save object


# plot error rates for sanity
plotErrors(err_2019, nominalQ=TRUE)
plotErrors(err_2020, nominalQ=TRUE)


# Infer sequence variants ####

# add names to filts
names(filt_2019) <- sample.names_2019
names(filt_2020) <- sample.names_2020

# Dereplication, sample inferrence, and merging ####
derep_2019 <- derepFastq(filt_2019)
derep_2020 <- derepFastq(filt_2020)

# DADA algorithm
dada_2019 <- dada(derep_2019, err=err_2019, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")
dada_2020 <- dada(derep_2020, err=err_2020, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")

# Construct sequence table ####
seqtab_2019 <- makeSequenceTable(dada_2019)
seqtab_2020 <- makeSequenceTable(dada_2020)

saveRDS(seqtab_2019,"./Output/Lee_2019_seqtab.RDS") # save object
saveRDS(seqtab_2020,"./Output/Lee_2020_seqtab.RDS") # save object

# Remove chimeras ####
seqtab.nochim_2019 <- removeBimeraDenovo(seqtab_2019, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim_2020 <- removeBimeraDenovo(seqtab_2020, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim_2019,"./Output/Lee_2019_seqtab.nochim.RDS")
saveRDS(seqtab.nochim_2020,"./Output/Lee_2020_seqtab.nochim.RDS")

# Assign taxonomy - UNITE+Euk / 80% bootstrap min ####
taxa_2019 <- assignTaxonomy(seqtab.nochim_2019,"./Taxonomy/sh_general_release_dynamic_all_10.05.2021.fasta.gz", minBoot = 80,multithread = TRUE)
taxa_2020 <- assignTaxonomy(seqtab.nochim_2020,"./Taxonomy/sh_general_release_dynamic_all_10.05.2021.fasta.gz", minBoot = 80,multithread = TRUE)

saveRDS(taxa_2019,"./Output/Lee_2019_taxa.RDS")
saveRDS(taxa_2020,"./Output/Lee_2020_taxa.RDS")

# rename seqtab object samples
seqtab.df_2019 <- as.data.frame(seqtab.nochim_2019)
seqtab.df_2020 <- as.data.frame(seqtab.nochim_2020)

# Create phyloseq object ####

# subset metadata to remove missing samples
# in.meta <- which(names(seqtab.nochim_2019[,1]) %in% lee_2019_meta$Run == TRUE)
# seqtab.nochim_2019 <- seqtab.nochim_2019[in.meta,]

in.seqtab_2019 <- which(lee_2019_meta$Run %in% row.names(seqtab.df_2019))
lee_2019_meta <- lee_2019_meta[in.seqtab_2019,]
in.seqtab_2020 <- which(lee_2020_meta$Run %in% row.names(seqtab.df_2020))
lee_2020_meta <- lee_2020_meta[in.seqtab_2020,]

# re-order
lee_2019_meta <- lee_2019_meta[order(lee_2019_meta$Run),]
row.names(lee_2019_meta) <- lee_2019_meta$Run
seqtab.nochim_2019 <- (seqtab.df_2019[row.names(lee_2019_meta),])
identical(row.names(seqtab.nochim_2019), as.character(lee_2019_meta$Run))

lee_2020_meta <- lee_2020_meta[order(lee_2020_meta$Run),]
row.names(lee_2020_meta) <- lee_2020_meta$Run
seqtab.nochim_2020 <- (seqtab.nochim_2020[row.names(lee_2020_meta),])
identical(row.names(seqtab.nochim_2020), as.character(lee_2020_meta$Run))


# make phyloseq object
otu_2019 <- otu_table(seqtab.nochim_2019,taxa_are_rows = FALSE)
met_2019 <- sample_data(lee_2019_meta)
tax_2019 <- tax_table(taxa_2019)

sample_names(met_2019) <- lee_2019_meta$Run
ps_2019 <- phyloseq(otu_2019,met_2019,tax_2019)

otu_2020 <- otu_table(seqtab.nochim_2020,taxa_are_rows = FALSE)
met_2020 <- sample_data(lee_2020_meta)
tax_2020 <- tax_table(taxa_2020)

sample_names(met_2020) <- lee_2020_meta$Run
ps_2020 <- phyloseq(otu_2020,met_2020,tax_2020)

# save it
saveRDS(ps_2019, "./Output/Lee_2019_raw_ps_object.RDS")
saveRDS(ps_2020, "./Output/Lee_2020_raw_ps_object.RDS")


# Identify and remove contaminants
blanks = which(ps_2019@sam_data$PCR_Negative == TRUE)
contamdf.prev <- isContaminant(ps_2019, neg=blanks, threshold = 0.01)
table(contamdf.prev$contaminant)

ps.noncontam_2019 <- prune_taxa(!contamdf.prev$contaminant, ps_2019)

blanks = which(ps_2020@sam_data$PCR_Negative == TRUE)
contamdf.prev <- isContaminant(ps_2020, neg=blanks, threshold = 0.01)
table(contamdf.prev$contaminant)

ps.noncontam_2020 <- prune_taxa(!contamdf.prev$contaminant, ps_2020)

# save contam-free phyloseq objects
saveRDS(ps.noncontam_2019, "./Output/Lee_2019_noncontam_ps_object.RDS")
saveRDS(ps.noncontam_2020, "./Output/Lee_2020_noncontam_ps_object.RDS")



#################################

