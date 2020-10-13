###########################################################
##  Processing raw 16S seqs into a phyloseq object       ##
##  This script processes "run1" - 1912KMI-0012          ##
##                                                       ##
##  Sonneratia alba microbiome                           ##
##                                                       ##
##  Author: Geoff Zahn - October 12, 2020                ##
##                                                       ##
##  Dependencies:                                        ##
##  dada2 v 1.14.1                                       ##
##  purrr v 0.3.4                                        ##
##  tidyverse v 1.3.0                                    ##
##  readxl v 1.3.1                                       ##
##  decontam v 1.6.0                                     ##
##  phyloseq v 1.30.0                                    ##
##                                                       ##
###########################################################


# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(tidyverse); packageVersion("tidyverse")
library(readxl); packageVersion("readxl")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")

# Load metadata ####
SA_meta <- readxl::read_xlsx("./Data/Sonneratia_alba/SA_16S_MetaData.xlsx")

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
path <- "./Data/Sonneratia_alba/Seqs/run1" 
fqs <- list.files(path, full.names = TRUE,recursive = TRUE,pattern = ".fastq.gz$")

# Parse fwd and rev reads
fnFs <- fqs[grep("R1.fastq.gz",fqs)]
fnRs <- fqs[grep("R2.fastq.gz",fqs)]

# Get Sample Names
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# subset metadata to this run's samples
meta <- full_meta %>% filter(SampleID %in% sample.names)
rm(SA_meta); rm(full_meta); rm(gps); rm(lat); rm(lon) # quick cleanup of environment


# Peek at quality profiles
plotQualityProfile(fnFs[c(1,30)]) # fwd reads
plotQualityProfile(fnRs[c(1,30)]) # rev reads

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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,160),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,verbose = TRUE)
saveRDS(out,"./Output/run1/out.RDS") # save object
out <- readRDS("./Output/run1/out.RDS") # reload point

filtpath <- file.path(path,"filtered")

# reassign filts for any potentially lost samples
filtFs <- list.files("./Data/Sonneratia_alba/Seqs/run1/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("./Data/Sonneratia_alba/Seqs/run1/filtered", pattern = "_R_filt", full.names = TRUE)

# Get Sample Names (again, just in case)
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"./Output/run1/errF.RDS") # save object
saveRDS(errR,"./Output/run1/errR.RDS") # save object
errF <- readRDS("./Output/run1/errF.RDS") # reload point
errR <- readRDS("./Output/run1/errR.RDS") # reload point


# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs,"./Output/run1/derepFs.RDS") # save object
derepFs <- readRDS("./Output/run1/derepFs.RDS") # reload point

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs,"./Output/run1/derepRs.RDS") # save object
derepRs <- readRDS("./Output/run1/derepRs.RDS") # reload point


# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inferrence ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs,"./Output/run1/dadaFs.RDS") # save object
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs,"./Output/run1/dadaRs.RDS") # save object

# Merge fwd and rev reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab,"./Output/run1/seqtab.RDS") # save object
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./output/run1/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./output/run1/seqtab.nochim.RDS")

# Assign taxonomy - Silva v132 exact match / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.nochim,"./Taxonomy/silva_nr99_v138_train_set.fa.gz", minBoot = 80,multithread = TRUE)
saveRDS(taxa,"./output/run1/taxa.RDS")

taxa <- addSpecies(taxa, "./Taxonomy/silva_species_assignment_v138.fa.gz")
saveRDS(taxa,"./output/run1/taxa_with_spp.RDS")

# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
# row.names(seqtab.df) <- unlist(map(strsplit(row.names(seqtab.df), "_"),1))

# Create phyloseq object ####

# subset to remove missing samples
meta <- full_meta

in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Pooled library ID (from GIS)` == TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]
dim(seqtab.nochim)
in.seqtab = which(meta$`Pooled library ID (from GIS)` %in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

# re-order
meta = meta[order(meta$`Pooled library ID (from GIS)`),]
row.names(meta) <- meta$`Pooled library ID (from GIS)`
seqtab.nochim <- (seqtab.nochim[row.names(meta),])
identical(row.names(seqtab.nochim), as.character(meta$`Pooled library ID (from GIS)`))

# check dimensions
dim(taxa)
dim(seqtab.nochim)
dim(meta)

# make phyloseq object ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))
# save it
saveRDS(ps, "./output/phyloseq_object_16S.RDS")

# Identify and remove contaminants
blanks = which(ps@sam_data$date == "Blank")
contamdf.prev <- isContaminant(ps, neg=blanks, threshold = 0.1)
table(contamdf.prev$contaminant)

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
identical(ps, ps.noncontam)
# no contaminants from blanks found via prevalence method in other samples!

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./output/phyloseq_object_16S_noncontam.RDS")


# Sequence stats:
Fs <- readRDS("./output/dadaFs.RDS")
Rs <- readRDS("./output/dadaRs.RDS")
taxa <- readRDS("./output/taxa.RDS")
seqtab.nochim = readRDS("./output/seqtab.nochim.RDS")
dim(seqtab.nochim)

x <- 1
seqscounts <- c()
for(i in names(Fs)){
  t <- Fs[x]
  r <- t[[1]]
  seqscounts[x] <- sum(unlist(r[1]))
  x=x+1
}

FwdReads <- sum(seqscounts)

x <- 1
seqscounts <- c()
for(i in names(Rs)){
  t <- Fs[x]
  r <- t[[1]]
  seqscounts[x] <- sum(unlist(r[1]))
  x=x+1
}

RevReads <- sum(seqscounts)
RevReads
