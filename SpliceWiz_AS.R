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
             size = 2,)+
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
       title = 'A) Ambient - Geothermal 18˚C')+
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

as_amb_plast = read_csv('Brain_plast_amb_AS.csv') %>% 
  filter(FDR < 0.05)
as_geo_plast = read_csv('Brain_plast_geo_AS.csv') %>% 
  filter(FDR < 0.05)
as_hyb_plast = read_csv("Brain_plast_hyb_AS.csv") %>% 
  filter(FDR < 0.05)



# AS ecotype hybrid divergence -----------------------------------------

as_div_amb_hyb_12 = read_csv("Brain_Div_amb_hyb_12_AS.csv") %>% 
  filter(FDR < 0.05)

as_div_geo_hyb_12 = read_csv('Brain_Div_geo_hyb_12_AS.csv') %>% 
  filter(FDR < 0.05)

as_div_amb_hyb_18 = read_csv("Brain_Div_amb_hyb_18_AS.csv") %>% 
  filter(FDR < 0.05)

as_div_geo_hyb_18 = read_csv('Brain_Div_geo_hyb_18_AS.csv') %>% 
  filter(FDR < 0.05)

