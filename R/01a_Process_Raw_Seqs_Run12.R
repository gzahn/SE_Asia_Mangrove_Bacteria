##  ###################################################  ##
##  Processing raw 16S seqs into a phyloseq object       ##
##  This script processes "run12" - 1912KMI-0012         ##
##                                                       ##
##  Sonneratia alba microbiome                           ##
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
library(ShortRead); packageVersion("ShortRead")

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
path <- "./Data/Sonneratia_alba/Seqs/run12" 
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



# CHECK FOR AND REMOVE PRIMER SITES WITH CUTADAPT ####

############# You will need to change these two values to match your data #############

# Here, you should supply the primer sequences used during PCR
# The example data in this workflow was generated using the 515f - 806r primer pair
# to amplify the V4 region of bacterial 16S rDNA
FWD <- "GTGCCAGCMGCCGCGGTAA" # Sequence of FWD primer
REV <- "GGACTACHVGGGTWTCTAAT"  # Sequence of REV primer

######################################################################################################

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients; REV.orients

# Prefilter to remove reads with ambiguous (N) bases ####
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # on Windows, set multithread = FALSE

# Discover primer matches, regardless of orientation ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Run cutadapt ####
# If the following command returns an error, you do not have cutadapt installed correctly
system2("cutadapt", args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, "--minimum-length 100", # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# sanity check
# This should show no occurences in any of the orientations now
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))



# Filter and trim ####

# cut fwd reads at 300 and rev reads at 200
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,200),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,verbose = TRUE)
saveRDS(out,"./Output/run12/out.RDS") # save object
out <- readRDS("./Output/run12/out.RDS") # reload point

filtpath <- file.path(path,"filtered")

# reassign filts for any potentially lost samples
filtFs <- list.files("./Data/Sonneratia_alba/Seqs/run12/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("./Data/Sonneratia_alba/Seqs/run12/filtered", pattern = "_R_filt", full.names = TRUE)

# Get Sample Names (again, just in case)
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"./Output/run12/errF.RDS") # save object
saveRDS(errR,"./Output/run12/errR.RDS") # save object
errF <- readRDS("./Output/run12/errF.RDS") # reload point
errR <- readRDS("./Output/run12/errR.RDS") # reload point


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

saveRDS(mergers,"./Output/run12/mergers.RDS")

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab,"./Output/run12/seqtab.RDS") # save object
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./Output/run12/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./Output/run12/seqtab.nochim.RDS")

# Assign taxonomy - Silva v132 exact match / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.nochim,"./Taxonomy/silva_nr99_v138_train_set.fa.gz", minBoot = 80,multithread = TRUE)
saveRDS(taxa,"./Output/run12/taxa.RDS")

taxa <- addSpecies(taxa, "./Taxonomy/silva_species_assignment_v138.fa.gz")
saveRDS(taxa,"./Output/run12/taxa_with_spp.RDS")

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
saveRDS(ps, "./Output/run12/raw_ps_object.RDS")


# Identify and remove contaminants
blanks = which(ps@sam_data$PCR_Negative == TRUE)
contamdf.prev <- isContaminant(ps, neg=blanks, threshold = 0.01)
table(contamdf.prev$contaminant)

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./Output/run12/noncontam_ps_object.RDS")


#################################

