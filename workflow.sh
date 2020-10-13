#!/bin/bash
set -euxo pipefail

#Rscript ./R/01a_Process_Raw_Seqs_Run12.R
Rscript ./R/01b_Process_Raw_Seqs_Run13.R
Rscript ./R/01c_Process_Raw_Seqs_Run6.R
Rscript ./R/01d_Process_Raw_Seqs_Run7.R
Rscript ./R/02_Merge_and_Clean_Phyloseq_Objects.R

