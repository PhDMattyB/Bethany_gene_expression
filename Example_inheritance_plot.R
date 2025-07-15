##############################
## example plot for SKR inhertiance patterns
##
## Matt Brachmann (PhDMattyB)
##
## 14.07.2025
##
##############################


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)

theme_set(theme_bw())


example_pal = c('#fee440', 
                '#1d58ab',
                '#00bbf9',
                '#9b5de5')


example_plot = ggplot()+
  geom_rect(data = data.frame(xmin = c(-0.32,
                                       -2, 
                                       -0.32, 
                                       0.32, 
                                       0.32, 
                                       -0.32), 
                              xmax = c(0.32, 
                                       2, 
                                       -2, 
                                       2, 
                                       2, 
                                       -2), 
                              ymin = c(-2, 
                                       -0.32, 
                                       0.32, 
                                       -0.32, 
                                       0.32, 
                                       -0.32), 
                              ymax = c(2, 
                                       0.32, 
                                       2, 
                                       -2, 
                                       2, 
                                       -2),
                              fill = c('Pure-strain1 dominant', 
                                       'Pure-strain2 dominant', 
                                       'Additive', 
                                       'Additive', 
                                       'Transgressive', 
                                       'Transgressive')), 
            aes(xmin = xmin, 
                xmax = xmax, 
                ymin = ymin, 
                ymax = ymax, 
                fill = fill), 
            alpha = 0.5)+
  scale_fill_manual(values = example_pal)+
  annotate("text", 
           x = c(-1.15, 1.15), 
           y = c(1.15, -1.15), 
           label = "Additive", 
           size = 2, 
           fontface = 'bold')+
  annotate("text", 
           x = c(1.15, -1.15), 
           y = c(1.15, -1.15), 
           label = "Transgressive", 
           size = 2, 
           fontface = 'bold')+
  annotate("text", 
           x = c(-1.15, 1.15), 
           y = c(0,0), 
           label = "Pure strain 2 dominant", 
           size = 2, 
           fontface = 'bold')+
  annotate("text", 
           x = c(0, 0), 
           y = c(1.15, -1.15), 
           label = "Pure strain 1 dominant", 
           angle = 90, 
           size = 2, 
           fontface = 'bold')+
  labs(x = 'log2(hybrid) - log2(pure strain 1)', 
       y = 'log2(hybrid) - log2(pure strain 2)')+
  geom_hline(yintercept = c(0.32, -0.32))+
  geom_vline(xintercept = c(0.32, -0.32))+
  theme(panel.grid = element_blank(), 
        legend.position = 'none')


ggsave('')