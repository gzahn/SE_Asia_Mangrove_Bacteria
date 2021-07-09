
# Simple short function to build df of relative abundance values for each unique phylum in each sample
###########################################
# library(phyloseq)
# library(tidyverse)
# library(purrr)


parse_phylum_relabund <- function(ps){

# For each phylum, subset to all taxa in that phylum, save otu_table in named list 
taxa_calls <- as(tax_table(ps)[,2],"character") # make rank a variable
ps_tsfrm <- transform_sample_counts(ps, function(x){x/sum(x)}) # make transformation formula a variable

# for-loop to build a list for each unique phylum
x <- 1
newlist <- list()
for(i in unique(taxa_calls)){
  tax <- subset_taxa(ps_tsfrm, ps_tsfrm@tax_table[,2]  ==  i) # make taxa rank a variable
  OTU = as(otu_table(tax), "matrix")
  if(taxa_are_rows(tax)){OTU <- t(OTU)}
  OTUdf = as.data.frame(OTU)
  newlist[i] <- OTUdf
}

taxa_df <- newlist %>% as.data.frame()
cat("Add predictor variables from original sample_data to 'phylum_ra_df' for plotting\n")

assign("phylum_ra_df",taxa_df,envir = .GlobalEnv)
return(taxa_calls)
}

#########################################

# parse_phylum_relabund(ps_genus)
# 
# 
# # build data frame
# phylum_ra_df$Structure <- ps_genus_ra@sam_data$Structure
# phylum_ra_df$Species <- ps_genus_ra@sam_data$Species
# phylum_ra_df$Location <- ps_genus_ra@sam_data$Location
# 
# # tidy for plotting
# tax_abund_order <- phylum_ra_df %>% select(unique(taxa_calls)) %>% map_dbl(mean) %>% sort(decreasing = TRUE) %>% names() # for all samples
# phylum_ra_df_long <- pivot_longer(phylum_ra_df, unique(taxa_calls),names_to="Phylum",values_to = "Relative_Abundance")
# 
# # remove zeroes for plotting  
# phylum_ra_df_long$Relative_Abundance[phylum_ra_df_long$Relative_Abundance == 0] <- NA
# 
# # plot top 6 phyla
# phylum_ra_df_long %>% 
#   filter(Phylum %in% tax_abund_order[1:6]) %>% 
#   ggplot(aes(x=Structure,y=Relative_Abundance,fill=Structure)) +
#   geom_boxplot() +
#   facet_wrap(~Phylum,scales = "free_y",) +
#   theme_bw() +
#   scale_fill_manual(values = pal) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size=12,face="bold"),
#         legend.title = element_text(size=14,face="bold"),
#         legend.text = element_text(size=10,face="bold"),
#         axis.title = element_text(size=14,face="bold"),
#         axis.text.x = element_text(angle=60,hjust=1,size=10))