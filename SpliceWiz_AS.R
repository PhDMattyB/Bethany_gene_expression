library(dplyr)
library(readr)
library(SpliceWiz)
library(limma)


# Stuff done on the cluster -----------------------------------------------


## Create file with all your bam names so it loops through them
bam_dir <- "/export/legacy_drives/jessie3/frombuzz/mkb6d/"

bam_files <- sort(list.files(bam_dir,
                             pattern = "\\.bam$",
                             full.names = TRUE,
                             recursive = TRUE))

manifest <- data.frame(sample = sub("\\.bam$", "", basename(bam_files)),
                       bam = bam_files) |> 
  as_tibble()

write_tsv(manifest,
          path.expand("~/SKR_BAM/bam_manifest.tsv"),
          sep = "\t",
          quote = FALSE,
          row.names = FALSE)

## splicewiz run

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

work_dir <- "/export/legacy_drives/jessie3/frombuzz/mkb6d"
ref_path <- "/home/mkb6d/SKR_BAM"

manifest <- read.delim(
  file.path(work_dir, "bam_manifest.tsv"),
  stringsAsFactors = FALSE
)

sample <- manifest$sample[task_id]
bam <- manifest$bam[task_id]

out_dir <- file.path(work_dir, "SpliceWiz_Output_v5")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Processing: ", sample)

processBAM(
  bamfiles = bam,
  sample_names = sample,
  reference_path = ref_path,
  output_path = out_dir,
  n_threads = 4,
  overwrite = FALSE,
  run_featureCounts = TRUE
)

message("Finished: ", sample)


## This next part can be run locally if you want

work_dir = path.extend('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Bim_Bam/Liver_bam/')

nxtse_path <- file.path(work_dir, "NxtSE_full_v5")

expr <- findSpliceWizOutput(out_dir)

collateData(Experiment = expr,
            reference_path = ref_path,
            output_path = nxtse_path,
            novelSplicing = FALSE,
            lowMemoryMode = TRUE)

se <- makeSE(nxtse_path,
             realize = FALSE,
             RemoveOverlapping = FALSE)

# saveRDS(se, file.path(work_dir, "SpliceWiz_full_se.rds"))

# setwd("~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Bim_Bam/Liver_bam/")
## save the nxste file
saveHDF5SummarizedExperiment(se,
                             dir = file.path(work_dir, "SpliceWiz_full_HDF5"),
                             replace = FALSE)

se.filtered <- se[applyFilters(se),]

psi <- SpliceWiz::makeMatrix(se,
                             event_list = rownames(se),
                             sample_list = colnames(se),
                             method = "PSI",
                             depth_threshold = 10,
                             na.percent.max = NA)

saveRDS(psi, file.path(work_dir, "SpliceWiz_PSI_matrix.rds"))

## metadata file
metadata <- read_tsv("SKR_Brain_Metadata.tsv")

metadata <- metadata[match(colnames(se), metadata$sample_id),] 

metadata$ecotemp <- factor(paste0(metadata$ecotype,
                                  "_",
                                  metadata$temp),
                           levels = c( "SKRC_12", "SKRC_18",
                                       "SKRW_12", "SKRW_18",
                                       "SKRHYB_12", "SKRHYB_18"))



## Limma/EdgeR
mm <- model.matrix(~ 0 + ecotemp, data = metadata)
colnames(mm) <- levels(metadata$ecotemp)

rownames(mm) <- metadata$sample_id

## Run limma

contrast <- makeContrasts(eco12 = SKRC_12 - SKRW_12,
                          eco18 = SKRC_18 - SKRW_18,
                          plast_amb = SKRC_12 - SKRC_18,
                          plast_geo = SKRW_12 - SKRW_18,
                          plast_hyb = SKRHYB_12 - SKRHYB_18,
                          am_hyb_12 = SKRC_12 - SKRHYB_12,
                          am_hyb_18 = SKRC_18 - SKRHYB_18,
                          geo_hyb_12 = SKRW_12 - SKRHYB_12,
                          geo_hyb_18 = SKRW_18 - SKRHYB_18,
                          levels = mm)

SummarizedExperiment::colData(se)$ecotemp <-metadata$ecotemp
comparisons <- list(eco12 = c("SKRC_12", "SKRW_12"),
                    eco18 = c("SKRC_18", "SKRW_18"),
                    plast_amb = c("SKRC_12", "SKRC_18"),
                    plast_geo = c("SKRW_12", "SKRW_18"),
                    plast_hyb = c("SKRHYB_12", "SKRHYB_18"),
                    am_hyb_12 = c("SKRC_12", "SKRHYB_12"),
                    am_hyb_18 = c("SKRC_18", "SKRHYB_18"),
                    geo_hyb_12 = c("SKRW_12", "SKRHYB_12"),
                    geo_hyb_18 = c("SKRW_18", "SKRHYB_18"))


results <- lapply(comparisons, function(x) {
  SpliceWiz::ASE_edgeR(se = se.filtered,
                       test_factor = "ecotemp",
                       test_nom = x[1],
                       test_denom = x[2],
                       useQL = TRUE,
                       IRmode = "annotated")
})

names(results) <- names(comparisons)
# 
# results$eco12 |> 
#   as_tibble() |>
#   write_csv('Eco12_div_geo_amb_AS.csv')
# 
# results$eco18 |> 
#   as_tibble() |> 
#   write_csv("Eco18_div_geo_amb_AS.csv")
# 
# results$plast_amb |> 
#   as_tibble() |> 
#   write_csv("plast_amb_AS.csv")
# 
# results$plast_geo |> 
#   as_tibble() |> 
#   write_csv("plast_geo_AS.csv")
# 
# results$plast_hyb |> 
#   as_tibble() |> 
#   write_csv("plast_hyb_AS.csv")
# 
# results$am_hyb_12 |> 
#   as_tibble() |> 
#   write_csv("Div_amb_hyb_12_AS.csv")
# 
# results$am_hyb_18 |> 
#   as_tibble() |> 
#   write_csv("Div_amb_hyb_18_AS.csv")




# Post cluster analysis ---------------------------------------------------
library(tidyverse)

# Brain -------------------------------------------------------------------

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

## SE = Skipped Exon
## MXE = Mutually exclusive exons
## A5SS = Alternative 5' splice site usage
## A3ss = Alternative 3' splice site usage
## AFE/ALE = Alternative first or last exon usage
## IR or RI = Intron retention


volcano_cols = c('#127475',
                 '#f2542d', 
                 '#dedbd2')

theme_set(theme_bw())

# Brain AS divergence -----------------------------------------------------
as_eco12_full = read_csv('Brain_Eco12_div_geo_amb_AS.csv')
as_eco12 = read_csv('Brain_Eco12_div_geo_amb_AS.csv') %>% 
  filter(FDR < 0.05) %>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) 

as_eco18_full = read_csv('Brain_Eco18_div_geo_amb_AS.csv') 
as_eco18 = read_csv('Brain_Eco18_div_geo_amb_AS.csv') %>% 
  filter(FDR < 0.05) %>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) 

as_eco12$EventName == as_eco18$EventName

Brain_AS_eco12_volplot = ggplot(data = as_eco12_full, 
                                 aes(x = logFC, 
                                     y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_eco12, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'A) Ambient - Geothermal 12˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


Brain_AS_eco18_volplot = ggplot(data = as_eco18_full, 
                                aes(x = logFC, 
                                    y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_eco18, 
             aes(col = Regulated), 
             size = 2,)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'B) Ambient - Geothermal 18˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))
## Only 1 AS event is divergence across both temperatures


Brain_AS_Div_plots = Brain_AS_eco12_volplot + Brain_AS_eco18_volplot

ggsave('Alternative_Splicing_EcoDiv_Volcanoe_Plot.svg', 
       plot = Brain_AS_Div_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 10)

# Brain AS plasticity -----------------------------------------------------
as_amb_plast_full = read_csv('Brain_plast_amb_AS.csv')
as_amb_plast = read_csv('Brain_plast_amb_AS.csv') %>% 
  filter(FDR < 0.05)%>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) 
as_geo_plast_full = read_csv('Brain_plast_geo_AS.csv')
as_geo_plast = read_csv('Brain_plast_geo_AS.csv') %>% 
  filter(FDR < 0.05)%>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) 
as_hyb_plast_full = read_csv("Brain_plast_hyb_AS.csv")
as_hyb_plast = read_csv("Brain_plast_hyb_AS.csv") %>% 
  filter(FDR < 0.05)%>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) 

intersect(as_amb_plast$EventName, 
          as_geo_plast$EventName)
## 27 overlapp between amb and geo
intersect(as_amb_plast$EventName, 
          as_hyb_plast$EventName)
## 28 overlapp between amb and hyb

intersect(as_geo_plast$EventName, 
          as_hyb_plast$EventName)
## 50 overlapp between geo and hyb

intersect(as_geo_plast$EventName, 
          as_hyb_plast$EventName) %>% 
  intersect(., 
            as_amb_plast$EventName)
## 25 overlap between all 3 plastic ecotypes




Brain_AS_amb_plast_volplot = ggplot(data = as_amb_plast_full, 
                                aes(x = logFC, 
                                    y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_amb_plast, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'A) Ambient plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  # geom_hline(yintercept = 0.05, 
  #            col = '#f72585', 
  #            size = 1, 
  #            linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


Brain_AS_geo_plast_volplot = ggplot(data = as_geo_plast_full, 
                                    aes(x = logFC, 
                                        y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_geo_plast, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'B) Geothermal plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  # geom_hline(yintercept = 0.05, 
  #            col = '#f72585', 
  #            size = 1, 
  #            linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())


Brain_AS_hyb_plast_volplot = ggplot(data = as_hyb_plast_full, 
                                    aes(x = logFC, 
                                        y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_hyb_plast, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'C) Hybrid plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  # geom_hline(yintercept = 0.05, 
  #            col = '#f72585', 
  #            size = 1, 
  #            linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())

plasticity_plots = Brain_AS_amb_plast_volplot+Brain_AS_geo_plast_volplot+Brain_AS_hyb_plast_volplot


ggsave('Alternative_Splicing_Ecotype_plasticity_Volcanoe_Plot.svg', 
       plot = plasticity_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 25, 
       height = 15)


# AS ecotype hybrid divergence -----------------------------------------
as_div_amb_hyb_12_full = read_csv("Brain_Div_amb_hyb_12_AS.csv")
as_div_amb_hyb_12 = read_csv("Brain_Div_amb_hyb_12_AS.csv") %>% 
  filter(FDR < 0.05)%>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) %>% 
  mutate(Div_pattern = 'Amb_hyb_12')
as_div_geo_hyb_12_full = read_csv('Brain_Div_geo_hyb_12_AS.csv')
as_div_geo_hyb_12 = read_csv('Brain_Div_geo_hyb_12_AS.csv') %>% 
  filter(FDR < 0.05)%>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) %>% 
  mutate(Div_pattern = 'Geo_hyb_12')
as_div_amb_hyb_18_full = read_csv("Brain_Div_amb_hyb_18_AS.csv")
as_div_amb_hyb_18 = read_csv("Brain_Div_amb_hyb_18_AS.csv") %>% 
  filter(FDR < 0.05)%>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) %>% 
  mutate(Div_pattern = 'Amb_hyb_18')
as_div_geo_hyb_18_full = read_csv('Brain_Div_geo_hyb_18_AS.csv')
as_div_geo_hyb_18 = read_csv('Brain_Div_geo_hyb_18_AS.csv') %>% 
  filter(FDR < 0.05)%>% 
  mutate(Regulated = case_when(
    logFC < 0 ~ 'Down-Regulated',
    logFC > 0 ~ 'Up-Regulated', 
    logFC == 0 ~ 'Neutral')) %>% 
  mutate(Div_pattern = 'Geo_hyb_18')


intersect(as_div_amb_hyb_12$EventName, 
          as_div_geo_hyb_12$EventName)
## 260 overlapp between amb and geo
intersect(as_div_amb_hyb_18$EventName, 
          as_div_geo_hyb_18$EventName)
## 433 overlapp between amb and hyb

intersect(as_div_amb_hyb_12$EventName, 
          as_div_amb_hyb_18$EventName)
## 348 overlapp between geo and hyb

intersect(as_div_geo_hyb_12$EventName, 
          as_div_geo_hyb_18$EventName)
## 261 overlapp between geo and hyb

intersect(as_div_amb_hyb_12$EventName, 
          as_div_amb_hyb_18$EventName) %>% 
  intersect(., 
            as_div_geo_hyb_12$EventName) %>% 
  intersect(., 
            as_div_geo_hyb_18$EventName)
## 204 overlap between all 4 conditions





Brain_AS_amb_hyb_12_volplot = ggplot(data = as_div_amb_hyb_12_full, 
                                    aes(x = logFC, 
                                        y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_div_amb_hyb_12, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'A) Ambient - Hybrid divergence at 12°C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  # geom_hline(yintercept = 0.05, 
  #            col = '#f72585', 
  #            size = 1, 
  #            linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))



Brain_AS_geo_hyb_12_volplot = ggplot(data = as_div_geo_hyb_12_full, 
                                     aes(x = logFC, 
                                         y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_div_geo_hyb_12, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'B) Geothermal - Hybrid divergence at 12°C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  # geom_hline(yintercept = 0.05, 
  #            col = '#f72585', 
  #            size = 1, 
  #            linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())




Brain_AS_amb_hyb_18_volplot = ggplot(data = as_div_amb_hyb_18_full, 
                                     aes(x = logFC, 
                                         y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_div_amb_hyb_18, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'D) Ambient - Hybrid divergence at 18°C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  # geom_hline(yintercept = 0.05, 
  #            col = '#f72585', 
  #            size = 1, 
  #            linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))



Brain_AS_geo_hyb_18_volplot = ggplot(data = as_div_geo_hyb_18_full, 
                                     aes(x = logFC, 
                                         y = -log10(FDR)))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             col = '#dedbd2')+
  geom_point(data = as_div_geo_hyb_18, 
             aes(col = Regulated), 
             size = 2)+
  # scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = '-log(FDR)', 
       title = 'E) Geothermal - Hybrid divergence at 18°C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  # geom_hline(yintercept = 0.05, 
  #            col = '#f72585', 
  #            size = 1, 
  #            linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())


eco_hyb_div_combo = (Brain_AS_amb_hyb_12_volplot+Brain_AS_geo_hyb_12_volplot)/(Brain_AS_amb_hyb_18_volplot+Brain_AS_geo_hyb_18_volplot)


ggsave('Alternative_Splicing_Ecotype_Hybrid_div_Volcanoe_Plot.svg', 
       plot = eco_hyb_div_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 45, 
       height = 25)



# hybridization inheritance type ------------------------------------------
## quantify how many AS events are transgressive


# Inheritance pattern 12 degrees ------------------------------------------

div_hyb_12_common_genes = inner_join(as_div_amb_hyb_12, 
           as_div_geo_hyb_12, 
           by = c('EventName', 
                  'EventRegion', 
                  'EventType'))


### 12 degree inheritance pattern
exp_pattern_12 = div_hyb_12_common_genes %>% 
  dplyr::select(EventName,
                EventType, 
                EventRegion, 
                AvgPSI_SKRC_12, 
                AvgPSI_SKRHYB_12.x, 
                AvgPSI_SKRW_12, 
                AvgPSI_SKRHYB_12.y, 
                logFC.x, 
                logFC.y, 
                FDR.x, 
                FDR.y) %>% 
  dplyr::rename(am_hyb_12 = logFC.x, 
                geo_hyb_12 = logFC.y)

Nothing_12 = exp_pattern_12 %>% 
  filter(am_hyb_12 <= 0.32 & am_hyb_12 >= -0.32) %>% 
  filter(geo_hyb_12 < 0.32 & geo_hyb_12 > -0.32) %>% 
  mutate(exp_pattern = 'Neutral')
## 987 genes
dom_geo1_12 = exp_pattern_12 %>%
  filter(am_hyb_12 <= 0.32 &
           am_hyb_12 >= -0.32) %>%
  filter(geo_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')

dom_geo2_12 = exp_pattern_12 %>%
  filter(am_hyb_12 <= 0.32 &
           am_hyb_12 >= -0.32) %>%
  filter(geo_hyb_12 < -0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')


dom_amb1_12 = exp_pattern_12 %>% 
  filter(geo_hyb_12 <= 0.32 & geo_hyb_12 >= -0.32) %>%
  filter(am_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Ambient')
dom_amb2_12 = exp_pattern_12 %>% 
  filter(geo_hyb_12 <= 0.32 & geo_hyb_12 >= 0.32) %>%
  filter(am_hyb_12 < -0.32) %>% 
  mutate(exp_pattern = 'Dominant Ambient')

additive1_12 = exp_pattern_12 %>% 
  filter(am_hyb_12 > 0.32 & geo_hyb_12 < -0.32) %>% 
  mutate(exp_pattern = 'Additive')
additive2_12 = exp_pattern_12 %>% 
  filter(am_hyb_12 < -0.32 & geo_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Additive')

transgressive1_12 = exp_pattern_12 %>% 
  filter(am_hyb_12 > 0.32 & geo_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

transgressive2_12 = exp_pattern_12 %>% 
  filter(am_hyb_12 < -0.32 & geo_hyb_12 < -0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

# bind_rows(transgressive1_12,
#           transgressive2_12) %>%
#   write_csv('Brain_AS_Transgressive_expression_12degrees.csv')

exp_pattern_12_graph = bind_rows(dom_amb1_12, 
                                 dom_amb2_12, 
                                 dom_geo1_12,
                                 dom_geo2_12, 
                                 additive1_12, 
                                 additive2_12, 
                                 transgressive1_12, 
                                 transgressive2_12, 
                                 Nothing_12)

# inheritance_pal = c('#ff9e00',
#                     '#669bbc',
#                     '#c1121f', 
#                     '#adb5bd', 
#                     '#7400b8')

inheritance_pal = c('#fee440', 
                    '#1d58ab',
                    '#00bbf9',
                    '#adb5bd',
                    '#9b5de5')

AS_12_inheritance_pal = c('#9b5de5')

Inheritance_pattern_12 = ggplot(data = exp_pattern_12_graph, 
                                aes(x = geo_hyb_12, 
                                    y = am_hyb_12, 
                                    col = exp_pattern))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  geom_hline(yintercept = 0.32, 
             linetype = 'dashed')+
  geom_hline(yintercept = -0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = -0.32, 
             linetype = 'dashed')+
  labs(x = 'Hybrid 12°C - Geothermal 12°C', 
       y = 'Hybrid 12°C - Ambient 12°C',
       title = 'A)')+
  xlim(-8.5, 8.5)+
  ylim(-8.5, 8.5)+
  scale_color_manual(values = AS_12_inheritance_pal)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')



# Inheritance pattern 18 degrees ------------------------------------------

div_hyb_18_common_genes = inner_join(as_div_amb_hyb_18, 
                                     as_div_geo_hyb_18, 
                                     by = c('EventName', 
                                            'EventRegion', 
                                            'EventType'))


### 12 degree inheritance pattern
exp_pattern_18 = div_hyb_18_common_genes %>% 
  dplyr::select(EventName,
                EventType, 
                EventRegion, 
                AvgPSI_SKRC_18, 
                AvgPSI_SKRHYB_18.x, 
                AvgPSI_SKRW_18, 
                AvgPSI_SKRHYB_18.y, 
                logFC.x, 
                logFC.y, 
                FDR.x, 
                FDR.y) %>% 
  dplyr::rename(am_hyb_18 = logFC.x, 
                geo_hyb_18 = logFC.y)

Nothing_18 = exp_pattern_18 %>% 
  filter(am_hyb_18 <= 0.32 & am_hyb_18 >= -0.32) %>% 
  filter(geo_hyb_18 < 0.32 & geo_hyb_18 > -0.32) %>% 
  mutate(exp_pattern = 'Neutral')
## 987 genes
dom_geo1_18 = exp_pattern_18 %>%
  filter(am_hyb_18 <= 0.32 &
           am_hyb_18 >= -0.32) %>%
  filter(geo_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')

dom_geo2_18 = exp_pattern_18 %>%
  filter(am_hyb_18 <= 0.32 &
           am_hyb_18 >= -0.32) %>%
  filter(geo_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')


dom_amb1_18 = exp_pattern_18 %>% 
  filter(geo_hyb_18 <= 0.32 & geo_hyb_18 >= -0.32) %>%
  filter(am_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Ambient')
dom_amb2_18 = exp_pattern_18 %>% 
  filter(geo_hyb_18 <= 0.32 & geo_hyb_18 >= 0.32) %>%
  filter(am_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Dominant Ambient')

additive1_18 = exp_pattern_18 %>% 
  filter(am_hyb_18 > 0.32 & geo_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Additive')
additive2_18 = exp_pattern_18 %>% 
  filter(am_hyb_18 < -0.32 & geo_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Additive')

transgressive1_18 = exp_pattern_18 %>% 
  filter(am_hyb_18 > 0.32 & geo_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

transgressive2_18 = exp_pattern_18 %>% 
  filter(am_hyb_18 < -0.32 & geo_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

# bind_rows(transgressive1_18,
#           transgressive2_18) %>%
#   write_csv('Brain_AS_Transgressive_expression_18degrees.csv')

exp_pattern_18_graph = bind_rows(dom_amb1_18, 
                                 dom_amb2_18, 
                                 dom_geo1_18,
                                 dom_geo2_18, 
                                 additive1_18, 
                                 additive2_18, 
                                 transgressive1_18, 
                                 transgressive2_18, 
                                 Nothing_18)

# inheritance_pal = c('#ff9e00',
#                     '#669bbc',
#                     '#c1121f', 
#                     '#adb5bd', 
#                     '#7400b8')

# inheritance_pal = c('#fee440', 
#                     '#1d58ab',
#                     '#00bbf9',
#                     '#adb5bd',
#                     '#9b5de5')

AS_18_inheritance_pal = c('#9b5de5')

Inheritance_pattern_18 = ggplot(data = exp_pattern_18_graph, 
                                aes(x = geo_hyb_18, 
                                    y = am_hyb_18, 
                                    col = exp_pattern))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  geom_hline(yintercept = 0.32, 
             linetype = 'dashed')+
  geom_hline(yintercept = -0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = -0.32, 
             linetype = 'dashed')+
  labs(x = 'Hybrid 18°C - Geothermal 18°C', 
       y = 'Hybrid 18°C - Ambient 18c',
       title = 'B)')+
  xlim(-8.5, 8.5)+
  ylim(-8.5, 8.5)+
  scale_color_manual(values = AS_18_inheritance_pal)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')


trans_12 = bind_rows(transgressive1_12, 
                     transgressive2_12)

trans_18 = bind_rows(transgressive1_18, 
                     transgressive2_18)

inner_join(trans_12, 
           trans_18, 
           by = c('EventName', 
                  'EventType', 
                  'EventRegion'))


AS_Inheritance_combo = Inheritance_pattern_12+Inheritance_pattern_18
ggsave('Alternative_Splicing_Inheritance_Volcanoe_Plot.svg', 
       plot = AS_Inheritance_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 25)



# Transgressive gene expression and AS overlap ----------------------------

trans_gene_12 = read_csv('Trans_expression_12.csv')
trans_gene_18 = read_csv('Trans_expression_18.csv')

trans_as_18 = read_csv('Brain_AS_Transgressive_expression_18degrees.csv') %>% 
  separate(col = EventRegion, 
           c('chromosome', 
             'Region'), 
           sep = ':') %>% 
  separate(col = Region, 
           c('start', 
             'end'), 
           sep = '-') %>% 
  separate(col = end, 
           c('end', 
                   'direction'), 
           sep = '/') %>% 
  dplyr::select(-direction) %>% 
  mutate(header = 'chr') %>% 
  unite(col = chromosome, 
        c('header', 
          'chromosome'), 
        sep = '')%>% 
  mutate(start = as.numeric(start), 
         end = as.numeric(end)) %>% 
  mutate(mid = start + (end-start)/2)

trans_as_12 = read_csv('Brain_AS_Transgressive_expression_12degrees.csv')%>% 
  separate(col = EventRegion, 
           c('chromosome', 
             'Region'), 
           sep = ':') %>% 
  separate(col = Region, 
           c('start', 
             'end'), 
           sep = '-') %>% 
  separate(col = end, 
           c('end', 
             'direction'), 
           sep = '/') %>% 
  dplyr::select(-direction) %>% 
  mutate(header = 'chr') %>% 
  unite(col = chromosome, 
        c('header', 
          'chromosome'), 
        sep = '') %>% 
  mutate(start = as.numeric(start), 
         end = as.numeric(end))%>% 
  mutate(mid = start + (end-start)/2)

splice_gene_overlap_12 = trans_as_12 %>%
  inner_join(trans_gene_12,
    by = join_by(chromosome,
      between(mid, start, end)),
    relationship = "many-to-many")%>% 
  dplyr::rename(AS_amb_hyb_12 = am_hyb_12.x, 
                AS_geo_hyb_12 = geo_hyb_12.x, 
                DEG_amb_hyb_12 = am_hyb_12.y, 
                DEG_geo_hyb_12 = geo_hyb_12.y)


splice_gene_overlap_18 = trans_as_18 %>%
  inner_join(trans_gene_18,
             by = join_by(chromosome,
                          between(mid, start, end)),
             relationship = "many-to-many") %>% 
  dplyr::rename(AS_amb_hyb_18 = am_hyb_18.x, 
                AS_geo_hyb_18 = geo_hyb_18.x, 
                DEG_amb_hyb_18 = am_hyb_18.y, 
                DEG_geo_hyb_18 = geo_hyb_18.y)


# splice_gene_overlap_12 %>% 
#   write_csv('AS_DEG_Overlap_12.csv')
# splice_gene_overlap_18 %>% 
#   write_csv('AS_DEG_Overlap_18.csv')

ggplot(data = splice_gene_overlap_12, 
       aes(x = AS_amb_hyb_12, 
           y = DEG_amb_hyb_12))+
  geom_point()+
  geom_hline(yintercept = 0.0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.0, 
             linetype = 'dashed')+
  labs(x = 'Alternative splicing logFC Ambient - Hybrid 12°C', 
       y = 'DEG logFC Ambient - Hybrid 12°C')

ggplot(data = splice_gene_overlap_18, 
       aes(x = AS_amb_hyb_18, 
           y = DEG_amb_hyb_18))+
  geom_point()+
  geom_hline(yintercept = 0.0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.0, 
             linetype = 'dashed')+
  labs(x = 'Alternative splicing logFC Ambient - Hybrid 18°C', 
       y = 'DEG logFC Ambient - Hybrid 18°C')
