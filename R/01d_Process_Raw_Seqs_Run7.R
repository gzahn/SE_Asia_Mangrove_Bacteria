##  ###################################################  ##
##  Processing raw 16S seqs into a phyloseq object       ##
##  This script processes "run7" - 1912KMI-0007          ##
##                                                       ##
##  Avicenia alba microbiome                             ##
##                                                       ##
##  Author: Geoff Zahn - October 12, 2020                ##
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

# Load metadata ####
SA_meta <- readxl::read_xlsx("./Data/Avicenia_alba/AA_16S_MetaData.xlsx")

# note: relevant "BLANK" fastq files were renamed as AA_BLANK_XXX and SA_BLANK_XXX respectively in bash

# edit colnames
names(SA_meta)[1:5] <- c("SampleID","Location","Species","Structure","GPS")

# remove extra cols
full_meta <- SA_meta %>% select(c("SampleID","Location","Species","Structure","GPS"))

# clean up GPS (all sites in N/E decimal-degree coordinates)
full_meta$GPS[full_meta$GPS == "Blank"] <- NA
gps <- str_split(full_meta$GPS,pattern = " ")
gps <- lapply(gps, `length<-`, max(lengths(gps)))

lat <- map_chr(gps,1) %>% 
  str_remove("N") %>% 
  as.numeric()

lon <- map_chr(gps,2) %>% 
  str_remove("E") %>% 
  as.numeric()

full_meta$Lat <- lat
full_meta$Lon <- lon

# Add column identifying PCR negatives
full_meta$PCR_Negative <- FALSE
full_meta$PCR_Negative[grep("BLANK",full_meta$SampleID)] <- TRUE

# Find raw fastq files and prepare workspace ####
path <- "./Data/Avicenia_alba/Seqs/run7" 
fqs <- list.files(path, full.names = TRUE,recursive = TRUE,pattern = "R1.fastq.gz$|R2.fastq.gz$")

# Parse fwd and rev reads
fnFs <- fqs[grep("R1.fastq.gz",fqs)]
fnRs <- fqs[grep("R2.fastq.gz",fqs)]

# Get Sample Names
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# subset metadata to this run's samples
meta <- full_meta %>% filter(SampleID %in% sample.names)
rm(SA_meta); rm(full_meta); rm(gps); rm(lat); rm(lon) # quick cleanup of environment


# Peek at quality profiles
# plotQualityProfile(fnFs[c(1,30)]) # fwd reads
# plotQualityProfile(fnRs[c(1,30)]) # rev reads

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# make new directory for filtered files
if(!dir.exists(file.path(path,"filtered"))){
  dir.create(file.path(path,"filtered"))
}

# check for duplicated sample names
sum(duplicated(sample.names))

# Filter and trim ####

# cut fwd reads at 300 and rev reads at 200
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,200),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,verbose = TRUE)
saveRDS(out,"./Output/run7/out.RDS") # save object
out <- readRDS("./Output/run7/out.RDS") # reload point

filtpath <- file.path(path,"filtered")

# reassign filts for any potentially lost samples
filtFs <- list.files("./Data/Avicenia_alba/Seqs/run7/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("./Data/Avicenia_alba/Seqs/run7/filtered", pattern = "_R_filt", full.names = TRUE)

# Get Sample Names (again, just in case)
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"./Output/run7/errF.RDS") # save object
saveRDS(errR,"./Output/run7/errR.RDS") # save object
errF <- readRDS("./Output/run7/errF.RDS") # reload point
errR <- readRDS("./Output/run7/errR.RDS") # reload point


# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# Infer sequence variants ####

# add names to filts
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Dereplication, sample inferrence, and merging ####

# loop through each pair of fwd and rev reads, one file at a time
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

saveRDS(mergers,"./Output/run7/mergers.RDS")

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab,"./Output/run7/seqtab.RDS") # save object
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./Output/run7/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./Output/run7/seqtab.nochim.RDS")

# Assign taxonomy - Silva v132 exact match / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.nochim,"./Taxonomy/silva_nr99_v138_train_set.fa.gz", minBoot = 80,multithread = TRUE)
saveRDS(taxa,"./Output/run7/taxa.RDS")

taxa <- addSpecies(taxa, "./Taxonomy/silva_species_assignment_v138.fa.gz")
saveRDS(taxa,"./Output/run7/taxa_with_spp.RDS")

# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df)

# Create phyloseq object ####

# subset to remove missing samples
in.meta <- which(names(seqtab.nochim[,1]) %in% meta$SampleID == TRUE)
seqtab.nochim <- seqtab.nochim[in.meta,]
dim(seqtab.nochim)
in.seqtab <- which(meta$SampleID %in% names(seqtab.nochim[,1]))
meta <- meta[in.seqtab,]

# re-order
meta <- meta[order(meta$SampleID),]
row.names(meta) <- meta$SampleID
seqtab.nochim <- (seqtab.nochim[row.names(meta),])
identical(row.names(seqtab.nochim), as.character(meta$SampleID))


# make phyloseq object
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
met <- sample_data(meta)
tax <- tax_table(taxa)

sample_names(met) <- met$SampleID
ps <- phyloseq(otu,met,tax)

# save it
saveRDS(ps, "./Output/run7/raw_ps_object.RDS")


# Identify and remove contaminants
blanks = which(ps@sam_data$PCR_Negative == TRUE)
contamdf.prev <- isContaminant(ps, neg=blanks, threshold = 0.01)
table(contamdf.prev$contaminant)

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./Output/run7/noncontam_ps_object.RDS")


