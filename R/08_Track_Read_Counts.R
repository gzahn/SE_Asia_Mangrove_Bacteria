library(phyloseq)
library(tidyverse)


out <- readRDS("./Output/run13/out.RDS")
no.chim <- readRDS("./Output/run13/seqtab.nochim.RDS")
run13 <- out %>% 
  as.data.frame() %>% 
  mutate(no.chim = no.chim %>% rowSums(),sampleID = row.names(.))

out <- readRDS("./Output/run12/out.RDS")
no.chim <- readRDS("./Output/run12/seqtab.nochim.RDS")
run12 <- out %>% 
  as.data.frame() %>% 
  mutate(no.chim = no.chim %>% rowSums(),sampleID = row.names(.))

out <- readRDS("./Output/run6/out.RDS")
no.chim <- readRDS("./Output/run6/seqtab.nochim.RDS")
run6 <- out %>% 
  as.data.frame() %>% 
  mutate(no.chim = no.chim %>% rowSums(),sampleID = row.names(.))

out <- readRDS("./Output/run7/out.RDS")
no.chim <- readRDS("./Output/run7/seqtab.nochim.RDS")
run7 <- out %>% 
  as.data.frame() %>% 
  mutate(no.chim = no.chim %>% rowSums(),sampleID = row.names(.))

list(run12,run13,run6,run7) %>% 
  reduce(full_join) %>% 
  mutate(sampleID = sampleID %>% str_remove_all("_R1.fastq.gz")) %>% 
  mutate(Host = case_when(str_detect(sampleID,"SA") ~ "S. alba",
                          str_detect(sampleID,"AA") ~ "A. alba")) %>% 
  select(sampleID,reads.in,reads.out,no.chim) %>% 
  write_csv("Output/Stats/ReadStats.csv")
  
