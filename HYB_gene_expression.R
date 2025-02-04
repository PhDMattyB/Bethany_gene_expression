## reanalysing Bethanys hybrid data from her second data chapter of her PhD
## the data from the last data chapter needs to be fixed apparently. 
## using a GLMM to model gene expression differences


# Big package -------------------------------------------------------------

library(tidyverse)
library(lme4)
library(optimx)
library(sjPlot)



# Read the data -----------------------------------------------------------

setwd("/media/shaun/2d49aafa-914e-40e6-a1ea-142b85301ef9/rna_aligned/F1 hybrids")

brain_meta = read_csv('hyb_data.csv')
brain_data = read_tsv('brain_gene_read_counts_table_all_final.tsv')

names(brain_data)

liver_meta = read_csv('hyb_data.csv')
liver_data = read_tsv('liver_gene_read_counts_table_all_final.tsv')

