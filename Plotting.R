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
library(ggbeeswarm)
library(glue)

set.seed(123)
class_colors <- c(Luminal = "#51abcb", Basal = "red", `Blood Vessel` = "green", Undefined = "purple")


cells <- read_csv("E:/Final Data/quantification/joined_data.csv")
patient_data <- read_csv("E:/[1]Datasets/[5]Quantification/log/patient_metadata.csv")

# Patient case statements to include descriptive columns
patient_data <- patient_data %>%
  mutate(diagnosis = case_when(
    (`ER STATUS` == 1 | `ER STATUS` == 0) & (`PR STATUS` == 1 | `PR STATUS` == 0) & (`HER2 STATUS` == 1 | `HER2 STATUS` == 0)~ "Triple negative",
    (`ER STATUS` == 2 | `PR STATUS` == 2) & (`HER2 STATUS` == 1 | `HER2 STATUS` == 0) ~ "Luminal A",
    (`ER STATUS` == 2 | `PR STATUS` == 2) & (`HER2 STATUS` == 2 |`HER2 STATUS` == 1 | `HER2 STATUS` == 0) ~ "Luminal B",
    (`ER STATUS` == 1 | `ER STATUS` == 0) & (`PR STATUS` == 1 | `PR STATUS` == 0) & `HER2 STATUS` == 2 ~ "HER2 enriched",
    TRUE ~ "Unknown"
  ))

patient_data <- patient_data %>%
  mutate(Histology_description = case_when(
    HISTOLOGY == 0 ~ "DCIS",
    HISTOLOGY == 1 ~ "Invasive ductal",
    HISTOLOGY == 2 ~ "Lobular",
    HISTOLOGY == 3 ~ "Invasive mucinous",
    HISTOLOGY == 4 ~ "Mixed",
    HISTOLOGY == 5 ~ "Papillary"
  ))

# Joining patient and single cell data
cells <- inner_join(cells, patient_data, by = "Full_Filename")


colnames(cells)[8] <- "ANAX1"
markers_columns <- c("Nuclear","CD31","CK5","SMA","Ki67","CK8","ANAX1","ECAD","PCNT","CCASP3")

# Normalisation of columns
cells[,markers_columns] <- cells[,markers_columns] + 1
cells[,markers_columns] <- log(cells[markers_columns],10)


cells <- cells %>% 
        filter(rowSums(pick(markers_columns)) > 6)

#Renaming columns for IMHC tools
names(cells)[names(cells) == "Filename"] <- "imageNb"
names(cells)[names(cells) == "centroid-0"] <- "Pos_X"
names(cells)[names(cells) == "centroid-1"] <- "Pos_Y"


#----------------------------------------------------------------------------------------------------------#
#                                        Classification                                                    #
#----------------------------------------------------------------------------------------------------------#

cells <- cells %>%
  mutate(class = case_when(
    (SMA >= 1.5) & (ECAD >= 1.5)~ "BnL",
    (SMA >= 1.5) & (CK8 >= 1.5)~ "LnB",
    CK8 >= 1.5 ~ "Luminal",
    ECAD >= 1.5 ~ "Luminal",
    SMA >= 1.5 ~ "Basal",
    CD31 >= 1.5 ~ "Blood Vessel",
    TRUE ~ "Undefined"
  ))

cells <- cells %>% 
  mutate(proliferating = case_when(
    Ki67 >= 1 ~ "Positive",
    Ki67 < 1 ~ "Negative"
  ))

cells <- cells %>%
  mutate(Type = case_when(
    grepl("normal", imageNb, ignore.case = TRUE) ~ "Normal",
    grepl("edge", imageNb, ignore.case = TRUE) ~ "Edge",
    grepl("tumour", imageNb, ignore.case = TRUE) ~ "Tumour",
    TRUE ~ "Undefined"
  ))

cells$Type <- factor(cells$Type, levels = c("Normal", "Edge", "Tumour"))
cells$diagnosis <- factor(cells$diagnosis, levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
cells$Histology_description <- factor(cells$Histology_description, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))



#----------------------------------------------------------------------------------------------------------#
#                                        Marker expression Plots                                           #
#----------------------------------------------------------------------------------------------------------#

for (marker in markers_columns){
  label_title <- paste(marker, "expression")
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
          geom_violin(aes(fill = class))+
          geom_boxplot(width=.1, outlier.shape=NA) +
          labs(x = "Patient Classification", y = label_title, title = "Scatter plot of Marker Correlations") +
          stat_compare_means(method = "t.test",label = "p.signif", comparisons = c("Luminal","Basal"))+
          facet_wrap(~cells$Type[(cells$class == "Luminal" | cells$class == "Basal")])+
          theme_pubr()
  )
}

for (marker in markers_columns){
  label_title <- paste(marker, "expression")
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
          geom_violin(aes(fill = class))+
          geom_boxplot(width=.1, outlier.shape=NA) +
          labs(x = "Patient Histology" , y = label_title, title = "Scatter plot of Marker Correlations") +
          facet_wrap(~cells$Histology_description[(cells$class == "Luminal" | cells$class == "Basal")])+
          theme_pubr()
  )
}

for (marker in markers_columns){
  label_title <- paste(marker, "expression")
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = diagnosis)) +
          geom_violin(aes(fill = class),outlier.shape=NA)+
          geom_boxplot(width=.1, outlier.shape=NA) +
          labs(x = "Patient Diagnosis", y = label_title, title = "Violin plot of marker expression in different molecular subtypes") +
          # facet_wrap(~cells$diagnosis[(cells$class == "Luminal" | cells$class == "Basal")])+
          theme_pubr()
  )
}

stat.test <- cells[cells$class == "Basal" | cells$class == "Luminal",] %>%
  group_by(diagnosis) %>%
  t_test(ANAX1 ~ class) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "diagnosis", dodge = 1) 

ggplot(cells[(cells$class == "Basal" | cells$class == "Luminal"),], aes(x = diagnosis, y = ANAX1))+
  geom_violin(aes(fill = class))+
  theme_pubr()+
  scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
  stat_compare_means(method = "anova",label = "p.signif", label.y = 2.6, label.x = 2.45)+
  stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0.008, y.position = 2.5)+
  theme(legend.position = "top")

ggplot(cells[cells$class == "Basal" | cells$class == "Luminal",], aes(x = proliferating, y = ANAX1))+
  geom_boxplot(aes(fill = class))+
  theme_pubr()+
  scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
  theme(legend.position = "top")+
  stat_compare_means(method = "anova",label = "p.signif")+
  facet_wrap(~cells$Histology_description[cells$class == "Basal" | cells$class == "Luminal"])






# Plot violin plot
ggplot(cells[(cells$class == "Basal" | cells$class == "Luminal"),], aes(x = Histology_description, y = ANAX1)) +
  geom_violin(aes(fill = class)) +
  theme_pubr() +
  scale_y_continuous(breaks = seq(0, 2.6, by = 0.5)) +
  stat_compare_means(method = "anova", label = "p.signif", label.y = 2.6, label.x = 2.45) +
  stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0.008, y.position = 2.5) +
  theme(legend.position = "top") 

ggplot(cells[cells$class == "Luminal",], aes(x = Histology_description, y = ANAX1))+
  geom_violin(fill = "#51abcb")+
  geom_boxplot(width = .1,outlier.shape=NA)+
  stat_compare_means(method = "anova", label = "p.signif", label.y = 2.6, label.x = 2.45) +
  theme_pubr()
  


#----------------------------------------------------------------------------------------------------------#
#                                        Cell Count Plots                                                  #
#----------------------------------------------------------------------------------------------------------#


cell_counts <- cells %>%
  group_by(imageNb, Type, class, diagnosis, proliferating) %>%
  summarise(Count = n()) %>%
  ungroup()


ggplot(cell_counts, aes(x = Type, y = Count, fill = class)) +
  geom_bar(stat = "identity") +
  labs(x = "Image", y = "Count", fill = "Cell Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion of Cells in Each Image")

cell_counts <- cell_counts %>%
  group_by(Type) %>%
  mutate(Percentage = Count / sum(Count))

cell_counts$Percentage <- as.numeric(cell_counts$Percentage)


ggplot(cell_counts, aes(x = Type, y = Percentage, fill = class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "Image", y = "Percentage", fill = "Cell Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion of Cells in Each Image")





create_violin_plot <- function(marker, diagnosis_col) {
  # Subset the data for the desired classes
  # subset_cells <- cells[(cells$class == "Basal" | cells$class == "Luminal" | cells$class == "BnL"),]
  subset_cells <- cells
  
  # Convert marker to numeric if it is not
  subset_cells[[marker]] <- as.numeric(as.character(subset_cells[[marker]]))
  
  # Generate the statistical test
  stat.test <- subset_cells %>%
    group_by(.data[[diagnosis_col]]) %>%
    t_test(as.formula(paste(marker, "~ class"))) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = diagnosis_col)
  
  # Create the ggplot
  plot <- ggplot(subset_cells, aes(x = .data[[diagnosis_col]], y = .data[[marker]])) +
    geom_violin(aes(fill = class)) +
    geom_boxplot(aes(fill = class)) +
    labs(x = "Patient Diagnosis", y = paste(marker, "expression"), title = "Violin plot of marker expression in different molecular subtypes") +
    theme_pubr() +
    scale_y_continuous(breaks = seq(0, 2.6, by = 0.5)) +
    stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0.008, y.position = 2.5) +
    theme(legend.position = "top")
  
  print(plot)
}

# Example usage:
create_violin_plot("ANAX1", "Type")


#----------------------------------------------------------------------------------------------------------#
#                                           UMAP Plots                                                     #
#----------------------------------------------------------------------------------------------------------#

library(umap)


umap_df = cells[markers_columns] %>% umap(n_components = 2, n_neighbors = 15) # 2 and 15
marker_kmeans <- kmeans(umap_df$data, centers = 8)
umap_marker_df <- as.data.frame(umap_df$layout) %>% merge(cells, by = "row.names", all = TRUE) 


umap_marker_df <- umap_marker_df[markers_columns]
marker_name <- colnames(umap_marker_df)
marker_name <- markers_columns



umap_marker_df <- data.frame(umap_marker_df,marker_kmeans$cluster)
umap_marker_df$marker_kmeans.cluster <- factor(umap_marker_df$marker_kmeans.cluster)

ggplot(data = umap_marker_df, aes(x = V1, y = V2, color = factor(marker_kmeans.cluster))) +
  geom_point() +
  theme_pubr() +  # Assuming you've defined the theme
  labs(color = "Cluster") +  # Adding legend title
  xlab("UMAP1") +  # Adding x-axis label
  ylab("UMAP2") +  # Adding y-axis label
  guides(color = guide_legend(title = "Cluster"))  # Adding legend title



# Create SCE objects
cells$fname <- NULL
counts <- cells[,markers_columns]
cell_meta <- cells[, !(names(cells) %in% markers_columns)]
sce <- SingleCellExperiment(assays = list(counts = t(counts)))
colData(sce) <- as(cell_meta, "DataFrame")



spe <- buildSpatialGraph(sce, img_id = "imageNb", type = "expansion", threshold = 30)
# 
# 
# 
plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA NORMAL_A08.ome_4032_5376_4544_5888.tiff"],img_id = "imageNb", coords = c("Pos_Y","Pos_X"), colPairName = "expansion_interaction_graph", draw_edges = TRUE,
            node_size_by = "area", node_color_by = "class")


