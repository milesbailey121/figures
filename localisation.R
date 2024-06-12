library(SingleCellExperiment)
library(imcRtools)
library(ggplot2)
library(ggraph)
library(viridis)
library(tidyverse)
library(pheatmap)
library(scales)
library(BiocParallel)
library(data.table)
library(dplyr)
library(pals)
library(janitor)
library(ggpubr)
library(rstatix)

cells <- read_csv("C:/Users/miles/GitHub/cellpose-quantification/cell_intensity_results_normalised.csv")
markers_columns <- c("Nuclear","CD31","CK5","SMA","Ki67","CK8","ANAX1","ECAD","PCNT","CCASP3")


ggplot(cells,aes(x = Filename, y = ANAX1))+
  geom_boxplot()+
  stat_compare_means(method = "t.test",label.x = 1.5, label = "p.signif")


localisation <- data.frame(cells$ANAX1[cells$Filename == "nuc.tif"], cells$ANAX1[cells$Filename == "expanded.tif"])


localisation <- localisation %>% mutate(ratio = cells.ANAX1.cells.Filename.....expanded.tif.. / cells.ANAX1.cells.Filename.....nuc.tif..)


ggplot(localisation, aes(x = ratio))+
  geom_point()


plot(localisation$ratio)
