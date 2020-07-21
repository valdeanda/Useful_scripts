library(ggplot2)
library(scales)
library(magrittr)
library(dplyr)
#setwd("~/your/path")
df<-read.table('barplot.tab',sep='\t',header=TRUE)

df$Phylum<-
  factor(df$Phylum,
         levels= 
           c("Odinarchaeota",
             "Helarchaeota",
             "Thorarchaeota",
             "Unclassified_Asgard",
             "Lokiarchaeota ",
             "Heimdallarchaeota",
             "Parvarchaeota",
             "Diapherotrites",
             "Nanohaloarchaeota",
             "Altiarchaeota",
             "Nanoarchaeota",
             "Aenigmarchaeota",
             "Micrarchaeota",
             "Pacearchaeota",
             "Woesearchaeota",
             "Methanohalarchaeia",
             "Theionarchaea",
             "Hydrotermarchaeota",
             "Pontarchaea (MG-III)",
             "Hadasarchaeota",
             "Methanobacteria",
             "Archaeoglobi",
             "Thermococci",
             "Thalassoarchaea (MG-II)",
             "Unclassified_Euryarchaeota",
             "Halobacteria",
             "Thermoplasmata",
             "Methanomicrobia",
             "Aigarchaeota",
             "Geothermarchaeota",
             "Geoarchaeota",
             "Brockarchaeota",
             "Marsarchaeota",
             "Verstraetarchaoeta",
             "Unclassified_Crenarchaeota",
            "Korarchaeota",
             "Bathyarchaeota ",
             "Thaumarchaeota",
             "Crenarchaeota",
             "Unclassified"))

df%>%

ggplot(data = df, mapping= aes(y=Number.of.genomes, x=Phylum, fill = factor(Superphylum))) + 
  geom_bar(position = "dodge", stat = "identity") + 
  coord_flip()
