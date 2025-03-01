##############################
## Genes of interest
##
## Matt Brachmann (PhDMattyB)
##
## 29.01.2025
##
##############################

library(tidyverse)
library(data.table)
library(qvalue)



# Sample metadata ---------------------------------------------------------

setwd('~/Parsons_Postdoc/Bethany_gene_expression/')

brain_data = read_csv('F1_labfish_sampledata_brain_over10M(1).csv')

brain_data %>% 
  group_by(Pop_eco, 
           temp) %>% 
  summarize(num = n())

liver_data = read_csv('F1_labfish_sampledata10M_liver.csv')

liver_data %>% 
  group_by(popeco, 
           temp) %>% 
  summarize(num = n())

# stickleback annotation --------------------------------------------------



setwd('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/')

## extract all of the stickle genome annotation data
gene_annotation = read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
                           col_names = F, 
                           skip = 1) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  na.omit()


gene_metadata = gene_annotation %>% 
  filter(feature == 'gene') %>% 
  select(position,
         chromosome,
         feature)

ensemlbe_annotation_data = gene_annotation %>% 
  filter(feature == 'gene') %>% 
  pull(gene_id) %>% 
  as_tibble() %>% 
  separate(value, 
           into = c('ensemble_id', 
                    'gene_name',
                    'relationship'), 
           sep = ';')

ensemble_annotation_genes = ensemlbe_annotation_data %>% 
  separate(ensemble_id, 
           into = c('trash', 
                    'ensemble_id'), 
           sep = '=') %>%
  separate(gene_name, 
           into = c('trash', 
                    'gene_name'), 
           sep = '=') %>% 
  select(-trash) %>% 
  separate(ensemble_id, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '.CDS') %>% 
  select(-trash) 

annotation_data = bind_cols(gene_metadata, 
                            ensemble_annotation_genes)




# eco*temp qvalue  -----------------------------------------------------------


brain_goi = read_csv('~/Parsons_Postdoc/Bethany_gene_expression/GLMER_gene_expression_ecotemp_pval0.01.csv')

brain_goi_names = brain_goi %>% 
  select(gene_name, 
         mean_expression_relative, 
         ecow_temp18_pval) %>% 
  rename(ensemble_name = gene_name)
# brain eco*temp qvalue FDR -----------------------------------------------------

brain_ecotemp_pvalues <- goi$ecow_temp18_pval

## qvalue FDR
brain_ecotemp_qvalues = qvalue_truncp(p = brain_ecotemp_pvalues)

brain_ecotemp_qvalues = brain_ecotemp_qvalues$qvalues %>% 
  as_tibble() %>% 
  rename(qval = value)

## Benjami hochberg FDR
BH_ecotemp_pvals = p.adjust(brain_ecotemp_pvalues, 
                    method = 'hochberg', 
                    n = length(brain_ecotemp_pvalues))


brain_ecotemp_bh_pvals = BH_ecotemp_pvals %>% 
  as_tibble() %>% 
  rename(BH_pval = value)

brain_qvalue_ecotemp = bind_cols(brain_goi_names, 
                                 brain_ecotemp_qvalues) %>% 
  filter(qval < 0.05) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         ecow_temp18_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash)

brain_qvalue_ecotemp%>% 
  select(gene_name) %>% 
  # write_csv('Expression_gene_names_ecotemp_pval0.01.csv')
  write_tsv('BRAIN_qvalue_FDR_Expression_gene_names_ecotemp_pval0.05.txt')


brain_BH_ecotemp = bind_cols(brain_goi_names, 
                             brain_ecotemp_bh_pvals) %>% 
  filter(BH_pval < 0.05)%>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         ecow_temp18_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash)

brain_BH_ecotemp%>% 
  select(gene_name) %>% 
  # write_csv('Expression_gene_names_ecotemp_pval0.01.csv')
  write_tsv('BRAIN_BH_FDR_Expression_gene_names_ecotemp_pval0.05.txt')


brain_qvalue_ecotemp %>% 
  group_by(chromosome) %>% 
  arrange(desc(chromosome)) %>% 
  View()



brain_BH_ecotemp %>% 
  group_by(chromosome) %>% 
  arrange(desc(chromosome)) %>% 
  View()




# Brain - Graph - ecotype temp interaction  -------------------------------

brain_BH_ecotemp_graph = brain_BH_ecotemp %>% 
  rename(CHR = chromosome, 
         POS = position) %>% 
  stickle_CHR_reorder2() %>% 
  dist_cal()

brain_BH_axis_df = axis_df(brain_BH_ecotemp_graph)

ggplot(brain_BH_ecotemp_graph, 
       aes(x = POS, 
           y = mean_expression_relative))+
  geom_point(aes(colour = ))

ggplot(non_outs, 
       aes(x = {{xval}}, 
           y = {{yval}}))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(chr)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = out_col,
             alpha=0.8, 
             size=1.3)+
  scale_x_continuous(label = axisdf$CHR, 
                     breaks = axisdf$center)+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0,1.0))+
  # geom_hline(yintercept = 0.00043, 
  #            linetype = 2, 
  #            col = 'Black')+
  # ylim(0,1.0)+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = 'Fst', 
       title = plot_letter)+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                   angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))

# ecotype only results ----------------------------------------------------

## read in our data
goi = read_csv('~/Parsons_Postdoc/Bethany_gene_expression/GLMER_gene_expression_ecotemp_pval0.01.csv')

goi_names = goi %>% 
  select(gene_name, 
         mean_expression_relative, 
         ecoW_pval) %>% 
  rename(ensemble_name = gene_name)

## qvalue FDR
ecow_pvalues <- goi$ecoW_pval
qvalues = qvalue_truncp(p = ecow_pvalues)

eco_qvalues = qvalues$qvalues %>% 
  as_tibble() %>% 
  # arrange(value) 
  rename(eco_qvalues = value)

eco_BH_pvals = p.adjust(ecow_pvalues, 
                    method = 'hochberg', 
                    n = length(ecow_pvalues))


eco_BH_pval = eco_BH_pvals %>% 
  as_tibble() %>% 
  rename(eco_BH_pval = value)


goi_names = bind_cols(goi_names, 
                      eco_qvalues, 
                      eco_BH_pval)


eco_BH_data = goi_names %>% 
  filter(eco_BH_pval < 0.05) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         eco_BH_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash)

eco_qvalue_data = goi_names %>% 
  filter(eco_qvalues < 0.05)%>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         eco_qvalues,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash)


eco_BH_data %>% 
  select(gene_name) %>% 
  write_tsv('ECOTYPE_GLMER_genes_BH_FDR_0.05.txt')
eco_qvalue_data %>% 
  select(gene_name) %>% 
  write_tsv('ECOTYPE_GLMER_genes_qvalue_FDR_0.05.txt')

###
# pop*ecotype comparison --------------------------------------------------

goi = read_csv('~/Parsons_Postdoc/Bethany_gene_expression/GLMER_gene_expression_ecotemp_pval0.01.csv')

popeco_data = goi %>% 
  select(gene_name, 
         mean_expression_relative, 
         popMYV_ecoW_pval, 
         popSKR_ecoW_pval) %>% 
  rename(ensemble_name = gene_name)


## qvalue FDR
myv_eco_pval <- goi$popMYV_ecoW_pval
myv_eco_qvalues = qvalue_truncp(p = myv_eco_pval)
myv_eco_qvalues = myv_eco_qvalues$qvalues %>% 
  as_tibble() %>% 
  # arrange(value) 
  rename(myv_eco_qvalues = value)

skr_eco_pval <- goi$popSKR_ecoW_pval
skr_eco_qvalues = qvalue_truncp(p = skr_eco_pval)
skr_eco_qvalues = skr_eco_qvalues$qvalues %>% 
  as_tibble() %>% 
  # arrange(value) 
  rename(skr_eco_qvalues = value)



## bh fdr
myv_eco_BH_pvals = p.adjust(myv_eco_pval, 
                        method = 'hochberg', 
                        n = length(myv_eco_pval))

myv_eco_BH_pval = myv_eco_BH_pvals %>% 
  as_tibble() %>% 
  rename(myv_eco_BH_pval = value)

skr_eco_BH_pvals = p.adjust(skr_eco_pval, 
                            method = 'hochberg', 
                            n = length(skr_eco_pval))

skr_eco_BH_pval = skr_eco_BH_pvals %>% 
  as_tibble() %>% 
  rename(skr_eco_BH_pval = value)




popeco_datas = bind_cols(popeco_data, 
                      myv_eco_qvalues, 
                      skr_eco_qvalues,
                      myv_eco_BH_pval, 
                      skr_eco_BH_pval)


###
# MYV*ECO and SKR*ECO comparison ------------------------------------------
## compare the genes between the two fixed effect interactions

myv_eco_qvalue_genes = popeco_datas %>% 
  filter(myv_eco_qvalues < 0.05)

skr_eco_qvalue_genes = popeco_datas %>% 
  filter(skr_eco_qvalues < 0.05)

## 504 genes
myv_eco_genes = myv_eco_qvalue_genes %>% 
  select(ensemble_name)

## 187 genes
skr_eco_genes = skr_eco_qvalue_genes %>% 
  select(ensemble_name)

## 169 genes overlap
intersect(myv_eco_genes, 
          skr_eco_genes)

qvalue_intersection = inner_join(myv_eco_qvalue_genes, 
                                 skr_eco_qvalue_genes)


myv_eco_BH_genes = popeco_datas %>% 
  filter(myv_eco_BH_pval < 0.05)

skr_eco_BH_genes = popeco_datas %>% 
  filter(skr_eco_BH_pval < 0.05)

## 172 genes
myv_eco_genes = myv_eco_BH_genes %>% 
  select(ensemble_name)

## 47 genes
skr_eco_genes = skr_eco_BH_genes %>% 
  select(ensemble_name)

## 20 genes overlap
intersect(myv_eco_genes, 
          skr_eco_genes)

BH_intersection = inner_join(myv_eco_BH_genes, 
           skr_eco_BH_genes)

###
# MYV*ECO and SKR*ECO annotation results ----------------------------------


myv_eco_qvalue_genes %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         myv_eco_qvalues,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash) %>% 
  select(gene_name) %>% 
  write_tsv('~/Parsons_postdoc/bethany_gene_expression/MYVECO_GLMER_genes_Qvalue_FDR_0.05.txt')

myv_eco_BH_genes %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         myv_eco_BH_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash) %>% 
  select(gene_name) %>% 
  write_tsv('MYVECO_GLMER_genes_BH_FDR_0.05.txt')



skr_eco_qvalue_genes %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         skr_eco_qvalues,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash) %>% 
  select(gene_name) %>% 
  write_tsv('~/Parsons_postdoc/bethany_gene_expression/SKRECO_GLMER_genes_Qvalue_FDR_0.05.txt')

skr_eco_BH_genes %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         skr_eco_BH_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash) %>% 
  select(gene_name) %>% 
  write_tsv('SKRECO_GLMER_genes_BH_FDR_0.05.txt')

qvalue_intersection %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         skr_eco_qvalues,
         myv_eco_qvalues,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash) %>% 
  select(gene_name) %>% 
  write_tsv('SKR_MYV_ECO_Intersect_GLMER_genes_qvalue_FDR_0.05.txt')

BH_intersection%>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         skr_eco_BH_pval,
         myv_eco_BH_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash) %>% 
  select(gene_name) %>% 
  write_tsv('SKR_MYV_ECO_Intersect_GLMER_genes_BH_FDR_0.05.txt')





# LIVER expression annotation ---------------------------------------------

# LIVER eco*temp qvalue FDR -----------------------------------------------------

liver_goi = read_csv('~/Parsons_Postdoc/Bethany_gene_expression/GLMER_LIVER_gene_expression_ecotemp_pval0.01.csv')

liver_goi_names = liver_goi %>% 
  select(gene_name, 
         mean_expression_relative, 
         ecow_temp18_pval) %>% 
  rename(ensemble_name = gene_name)

liver_ecotemp_pval = liver_goi$ecow_temp18_pval
liver_ecotemp_qvalues = qvalue_truncp(p = liver_ecotemp_pval)
liver_ecotemp_qvalues = liver_ecotemp_qvalues$qvalues %>% 
  as_tibble() %>% 
  rename(qval = value)

liver_ecotemp_BH = p.adjust(liver_ecotemp_pval, 
                            method = 'hochberg', 
                            n = length(liver_ecotemp_pval))

liver_ecotemp_BH = liver_ecotemp_BH %>% 
  as_tibble() %>% 
  rename(BH_pval = value)


liver_ecotemp_BH = bind_cols(liver_goi_names, 
                         liver_ecotemp_BH) %>% 
  filter(BH_pval < 0.05) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name')%>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         ecow_temp18_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash)

liver_ecotemp_qvalue = bind_cols(liver_goi_names, 
                            liver_ecotemp_qvalues) %>% 
  filter(qval < 0.05) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name')%>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         ecow_temp18_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash)

liver_ecotemp_BH %>% 
  select(gene_name) %>% 
  write_tsv('~/Parsons_Postdoc/Bethany_gene_expression/LIVER_BH_FDR_GLMER_ecotemp_pval0.05.txt')

liver_ecotemp_qvalue %>% 
  select(gene_name) %>% 
  write_tsv('~/Parsons_Postdoc/Bethany_gene_expression/LIVER_qvalue_FDR_GLMER_ecotemp_pval0.05.txt')



# liver_ecotemp_BH %>% 
#   group_by(chromosome) %>% 
#   arrange(desc(chromosome)) %>% 
#   filter(chromosome == 'chrXXI') %>% 
#   filter(position >= 9963830, 
#          position <= 11574445) %>% 
#   arrange(position) %>% 
#   View()
# 
# liver_ecotemp_qvalue %>% 
#   group_by(chromosome) %>% 
#   arrange(desc(chromosome)) %>% 
#   filter(chromosome == 'chrXXI') %>% 
#   filter(position >= 9963830, 
#          position <= 11574445) %>% 
#   arrange(position) %>%  
#   View()


dim(liver_ecotemp_BH)

dim(liver_ecotemp_qvalue)



# Brain vs liver eco temp comparison --------------------------------------

intersect(brain_qvalue_ecotemp$gene_name, 
          liver_ecotemp_qvalue$gene_name) %>% 
  as_tibble()

Brain_Liver_qval_intersect = inner_join(brain_qvalue_ecotemp, 
           liver_ecotemp_qvalue, 
           by = c('gene_name', 
                  'ensemble_name', 
                  'chromosome'))
Brain_Liver_qval_intersect %>% 
  write_tsv('Qvalue_Brain_Liver_gene_ecotemp_intersection.txt')


 
intersect(brain_BH_ecotemp$gene_name, 
          liver_ecotemp_BH$gene_name) %>% 
  as_tibble()

Brain_liver_BH_intersect = inner_join(brain_BH_ecotemp, 
           liver_ecotemp_BH, 
           by = c('gene_name', 
                  'ensemble_name', 
                  'chromosome'))

Brain_liver_BH_intersect %>% 
  write_tsv('BH_brain_liver_gene_ecotemp_intersection.txt')


Brain_liver_BH_intersect %>% 
  filter(gene_name == 'glula')
