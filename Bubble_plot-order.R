library(tidyverse)
library(janitor)
# Read data ---------------------------------------------------------------####
data_to_plot<-read_delim("input_bubblePlot_AbsoluteFinal_2.tsv", 
                         delim = "\t") %>%
  clean_names()
#Example of the data frame 

#Groups  Colors  Glk     PfkABC  PK      PckA    Fbp     PorA    PflD    Ldh     pta     azf     pgl     PGD     edd     eda     CS      acnA    acnB    IDH     korA    korB
#C5 (22) Myxococcota     0       0       0.954545455     0.954545455     1       0       0.681818182     0.045454545     0.363636364     0       0       0.045454545     0
#U2 (32) Myxococcota     0.34375 0       0.90625 0.9375  0.125   0.375   0.375   0       0.09375 0.25    0.25    0.25    0       0.09375 0.90625 0.90625 0.09375 0.78125 0.937
#C1 (1)  Myxococcota     0       1       1       1       0       0       1       0       1       0       0       0       0       0       1       1       1       1       1
#C2 (2)  Myxococcota     0.5     0       1       1       0       0       0       0       0.5     0       0       0       0       0       0.5     1       0       0.5     1
#C3 (4)  Myxococcota     0       0.5     1       1       0.5     0.25    0.75    0       1       0       0       0.5     0       0       0.5     0.75    0.75    0.75    0.75
#C4 (2)  Myxococcota     0       0       1       1       1       0       0       0       1       1       1       1       0       0       1       1       0       1       1

# Transform to long -------------------------------------------------------####
data_long<-data_to_plot %>%
  pivot_longer(!c(groups, colors), names_to = "Genes", values_to = "count") %>%
  filter(count != "0")
# Order ------------------------------------------------------------------####
order_taxa<-unique(data_long$groups)
order_genes<-colnames(data_to_plot[,3:90])
# Color ------------------------------------------------------------------####
unique(data_long$colors)
color_pallet<-c("Myxococcota" = "#89CBED", "B64_69" = "#45AA98", 
                "Desulfuromonadia" = "#506297", "Desulfobacterota" = "#DCCD76",
                "Syntrophia" = "#AA4399")
# create the plot --------------------------------------------------------####
ggplot(data_long,
       aes(x=factor(groups, 
                    levels = order_taxa),
           y= factor(Genes, 
                     levels = rev(order_genes)),
           size= count,
           color=colors)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(1,5)) +
  scale_color_manual(values = color_pallet) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=6, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        axis.text.y = element_text(size=8))
