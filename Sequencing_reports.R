

# F1 lab fish sequencing reports ------------------------------------------

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

data = read_csv('F1_labfish_seq_results.csv')


data %>% 
  separate(SAMPLE, 
           into = c('Ecotype', 
                    'Date'), 
           sep = '_') %>% 
  mutate(Poppair = case_when(
    Ecotype == 'ASHNW' ~ 'ASHN', 
    Ecotype == 'ASHNC' ~ 'ASHN', 
    Ecotype == 'MYVW' ~ 'MYV', 
    Ecotype == 'MYVC' ~ 'MYV', 
    Ecotype == 'SKRW' ~ 'SKR', 
    Ecotype == 'SKRC' ~ 'SKR')) %>% 
  mutate(Ecotype = case_when(
    Ecotype == 'ASHNW' ~ 'Warm', 
    Ecotype == 'ASHNC' ~ 'Cold', 
    Ecotype == 'MYVW' ~ 'Warm', 
    Ecotype == 'MYVC' ~ 'Cold', 
    Ecotype == 'SKRW' ~ 'Warm', 
    Ecotype == 'SKRC' ~ 'Cold')) %>% 
  group_by(Ecotype, 
           Poppair, 
           Rearing_temp, 
           Tissue) %>% 
  filter(pass_or_fail == 'FALSE') %>% 
  summarize(min_seq = min(`aligned reads`), 
            mean_seq = mean(`aligned reads`), 
            max_seq = max(`aligned reads`), 
            num_sample = n()) %>% 
  arrange(Poppair, 
          Ecotype, 
          Rearing_temp, 
          Tissue) %>% 
  # write_tsv('F1_labfish_sequencing_fail_report.txt')
  write_csv('F1_labfish_sequencing_fail_report.csv')
  


data %>% 
  separate(SAMPLE, 
           into = c('Ecotype', 
                    'Date'), 
           sep = '_') %>% 
  mutate(Poppair = case_when(
    Ecotype == 'ASHNW' ~ 'ASHN', 
    Ecotype == 'ASHNC' ~ 'ASHN', 
    Ecotype == 'MYVW' ~ 'MYV', 
    Ecotype == 'MYVC' ~ 'MYV', 
    Ecotype == 'SKRW' ~ 'SKR', 
    Ecotype == 'SKRC' ~ 'SKR')) %>% 
  mutate(Ecotype = case_when(
    Ecotype == 'ASHNW' ~ 'Warm', 
    Ecotype == 'ASHNC' ~ 'Cold', 
    Ecotype == 'MYVW' ~ 'Warm', 
    Ecotype == 'MYVC' ~ 'Cold', 
    Ecotype == 'SKRW' ~ 'Warm', 
    Ecotype == 'SKRC' ~ 'Cold')) %>% 
  group_by(Ecotype, 
           Poppair, 
           Rearing_temp, 
           Tissue) %>% 
  filter(pass_or_fail == 'TRUE') %>% 
  summarize(min_seq = min(`aligned reads`), 
            mean_seq = mean(`aligned reads`), 
            max_seq = max(`aligned reads`),
            num_sample = n()) %>% 
  arrange(Poppair, 
          Ecotype, 
          Rearing_temp, 
          Tissue) %>% 
  # write_tsv('F1_labfish_sequencing_pass_report.txt')
  write_csv('F1_labfish_sequencing_pass_report.csv')
