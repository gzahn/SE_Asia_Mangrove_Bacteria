# -----------------------------------------------------------------------------#
# Reanalysis of Lee 2019 & 2020 Fungal data alongside bacterial reads
# Trimming conserved flanking regions away from ITS region
# Author: Geoffrey Zahn
# Software versions:  ITSxpress v 1.8.0
#                     R v 3.6.3
# -----------------------------------------------------------------------------#

#############################################################
#### This script calls itsxpress to remove flanking      ####
#### conserved regions of the rDNA when using fungal     ####
#### ITS1 or ITS2 metabarcodes.                          ####
#### You must have itsxpress installed on your system    ####
#### and present in your PATH. See itsxpress             ####
#### installation documents for instructions.            ####
#############################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")

#################################
#########  PARAMETERS  ##########

# Set the number of parallel threads to use (depends on your hardware):
nthreads <- 20

# Set the region (must be either "ITS1" or "ITS2")
itsregion <- "ITS1"

#################################
#################################

# EXTRACT ITS REGION FROM ALL "cutadapted" SEQUENCES ####

# load metadata
lee_2019_meta <- read_csv("./Data/Lee_2019/SRA_Metadata_Lee_2019_clean.csv")
lee_2020_meta <- read_csv("./Data/Lee_2020/SRA_Metadata_Lee_2020_clean.csv")


# Run ITSxpress on all forward fastq files (reverse reads often aren't used in fungal analyses)

# find the "cutadapted" files
fwds <- c(lee_2019_meta$Fwd_Path,lee_2020_meta$Fwd_Path) %>% 
  str_remove_all("\\./")
fwds <- paste0(getwd(),"/",fwds)
# build names for outfiles
outs <- paste0(str_remove_all(fwds,"pass_1.fastq.gz"), "ITS.fastq.gz")


# subset to just those that don't have associated ITS1 file already (avoid duplication)
itsfiles <- list.files("./Data/Lee_2019/fastq",pattern = "_ITS.fastq.gz$")
itsfiles2 <- list.files("./Data/Lee_2020/fastq",pattern = "_ITS.fastq.gz$")
itsfiles <- c(itsfiles,itsfiles2) %>% str_split("_") %>% map_chr(1)
non_its_fwds <- fwds %>% str_split("/") %>% map_chr(11) %>% str_split("_") %>% map_chr(1)
fwds <- fwds[!non_its_fwds %in% itsfiles]
non_its_outs <- outs %>% str_split("/") %>% map_chr(11) %>% str_split("_") %>% map_chr(1)
outs <- outs[!non_its_outs %in% itsfiles]


# for Lee 2020, "_1" is not present in file names
lee_2020_files <- paste0(getwd(),list.files("./Data/Lee_2020/fastq",full.names = TRUE) %>% str_remove_all("^\\."))
fwds[fwds %>% str_detect("Lee_2020")] <- fwds[fwds %>% str_detect("Lee_2020")] %>% str_replace_all("_pass_1","_pass")


# build the ITSxpress command and run it on each file in turn
ITSxpress_path <- "/home/gzahn/miniconda3/bin/itsxpress"

for(i in seq_along(fwds)){

  itsxpress_args <- c(
    paste0("--fastq ",fwds[i]),
    paste0("--single_end"),
    paste0("--outfile ",outs[i]),
    paste0("--region ",itsregion),
    paste0("--taxa Fungi"),
    paste0("--threads ",nthreads),
    paste0("--log ",outs[i],".log")
  )  
  
  system2(command = ITSxpress_path,
          args = itsxpress_args)
}


# This will generate a new set of output files, tagged with "ITS" in the filenames that can be used for the subsequent 
# filter-and-trim step in our workflow

# an example of the structure of a call to itsxpress is shown below:
"itsxpress --fastq ./data/filtN/2554_pass_1.fastq.gz --outfile ./data/filtN/2554_ITS.fastq.gz.ITS1.gz --region ITS1 --taxa Fungi --threads4 --log ./data/filtN/2554_ITS.fastq.gz.log --single-end"
