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
library(RColorBrewer)
library(dplyr)
library(pals)
library(janitor)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(glue)
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)
library(gplots)



set.seed(123)
class_colors <- c(Luminal = "#51abcb", Basal = "red", `Blood Vessel` = "green", Undefined = "grey")
Diag <- brewer.pal(4,'Set2')
Hist <- brewer.pal(5, 'Set2')

cells <- read_csv("E:/Final Data/quantification/joined_data.csv")
patient_data <- read_csv("E:/[1]Datasets/[5]Quantification/log/patient_metadata.csv")


#----------------------------------------------------------------------------------------------------------#
#                                         Patient Data                                                     #
#----------------------------------------------------------------------------------------------------------#


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

patient_data <- patient_data %>%
  mutate(type = case_when(
    grepl("normal", Full_Filename, ignore.case = TRUE) ~ "Normal",
    grepl("edge", Full_Filename, ignore.case = TRUE) ~ "Edge",
    grepl("tumour", Full_Filename, ignore.case = TRUE) ~ "Tumour",
    TRUE ~ "Undefined"
  ))


patient_data$diagnosis <- factor(patient_data$diagnosis, levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
patient_data$Histology_description <- factor(patient_data$Histology_description, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))
patient_data$type <- factor(patient_data$type, levels = c("Normal", "Edge", "Tumour"))

ggplot(data = as.data.frame(table(patient_data$diagnosis[unique(patient_data$`Local Patient ID`)])), aes(x = Var1, y = Freq,fill = Var1))+
  geom_bar(stat = "identity")+
  theme_pubr()+
  labs(x = "Patient Diagnosis", y = "Number of Patients") +
  scale_fill_manual("Diagnosis",values = Diag)+
  scale_y_continuous(breaks = seq(0, 105, by = 5))+
  geom_text(aes(label = paste("n = ", Freq), vjust = 0))

ggplot(data = as.data.frame(table(patient_data$Histology_description)), aes(x = Var1, y = Freq,fill = Var1))+
  geom_bar(stat = "identity")+
  theme_pubr()+
  labs(x = "Patient Histology", y = "Number of Samples") +
  scale_fill_manual("Histology",values = Hist)+
  scale_y_continuous(breaks = seq(0, 700, by = 50))+
  geom_text(aes(label = paste("n = ", Freq), vjust = 0))

ggplot(as.data.frame(table(patient_data$type[!(patient_data$type == "Undefined")])),aes(x = Var1, y = Freq,fill = Var1))+
  geom_bar(stat = "identity")+
  labs(x = "Pathology Classificaiton", y = "Number of Images") +
  theme_pubr()+
  scale_fill_manual("Pathology Classificaiton",values =c("gray","#e8e54f","#e84fe0","lightblue"))+
  scale_y_continuous(breaks = seq(0, 200, by = 10))+
  geom_text(aes(label = paste("n = ", Freq), vjust = 0))

# Joining patient and single cell data
cells <- inner_join(cells, patient_data, by = "Full_Filename")


colnames(cells)[8] <- "ANAX1"
markers_columns <- c("Nuclear","CD31","CK5","SMA","Ki67","CK8","ANAX1","ECAD","PCNT","CCASP3")


cells[,markers_columns] <- cells[,markers_columns] %>% mutate(across(everything(), ~ ( . - min(.)) / (max(.) - min(.))* 10))


# cells[,markers_columns] <- cells[,markers_columns] + 1
# cells[,markers_columns] <- log(cells[markers_columns],10)

unfiltered_cells <- cells 

cells <- cells %>% 
        filter(rowSums(pick(markers_columns)) > 6)

#Renaming columns for IMHC tools
names(cells)[names(cells) == "Filename"] <- "imageNb"
names(cells)[names(cells) == "centroid-0"] <- "Pos_X"
names(cells)[names(cells) == "centroid-1"] <- "Pos_Y"


#----------------------------------------------------------------------------------------------------------#
#                                        Classification                                                    #
#----------------------------------------------------------------------------------------------------------#

# cells <- cells %>%
#   mutate(class = case_when(
#     (SMA >= 1.5) & (ECAD >= 1.5)~ "BnL",
#     (SMA >= 1.5) & (CK8 >= 1.5)~ "LnB",
#     CK8 >= 1.5 ~ "Luminal",
#     ECAD >= 1.5 ~ "Luminal",
#     SMA >= 1.5 ~ "Basal",
#     CD31 >= 1.5 ~ "Blood Vessel",
#     TRUE ~ "Undefined"
#   ))

cells <- cells %>%
  mutate(class = case_when(
    CK8 >= 1 ~ "Luminal",
    ECAD >= 1~ "Luminal",
    SMA >= 1 ~ "Basal",
    CD31 >= 1 ~ "Blood Vessel",
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

ggplot(cells,aes(x = ECAD))+
  geom_density(stat = "bin",fill = "lightblue")+
  scale_y_continuous(name="Cell Counts",labels = comma)+
  theme_pubr()

ggplot(unfiltered_cells,aes(x = ECAD))+
  geom_density(stat = "bin",fill = "lightblue")+
  scale_y_continuous(name="Cell Counts",breaks = seq(0, 400000, by = 50000) ,labels = comma)+
  theme_pubr()

ggplot(cells,aes(x = CK8))+
  geom_histogram(stat = "bin",fill = "grey")


#----------------------------------------------------------------------------------------------------------#
#                                        Marker expression Plots                                           #
#----------------------------------------------------------------------------------------------------------#

for (marker in markers_columns){
  label_title <- paste(marker, "expression")
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
          geom_boxplot(aes(fill = class))+
          # geom_boxplot(width=.1, outlier.shape=NA) +
          labs(x = "Patient Classification", y = label_title, title = "Scatter plot of Marker Correlations") +
          stat_compare_means(method = "t.test",label = "p.signif", comparisons = c("Luminal","Basal"))+
          facet_wrap(~cells$Type[(cells$class == "Luminal" | cells$class == "Basal")])+
          theme_pubr()
  )
}

for (marker in markers_columns){
  label_title <- paste(marker, "expression")
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
          geom_boxplot(aes(fill = class))+
          labs(x = "Patient Histology" , y = label_title, title = "Scatter plot of Marker Correlations") +
          facet_wrap(~cells$Histology_description[(cells$class == "Luminal" | cells$class == "Basal")])+
          theme_pubr()
  )
}

for (marker in markers_columns){
  label_title <- paste(marker, "expression")
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = diagnosis)) +
          geom_boxplot(aes(fill = class),outlier.shape=NA)+
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

create_violin_plot <- function(marker, diagnosis_col) {
  # Subset the data for the desired classes
  subset_cells <- cells[(cells$class == "Basal" | cells$class == "Luminal") & !(is.na(cells$Type)),]
  
  
  # Convert marker to numeric if it is not
  subset_cells[[marker]] <- as.numeric(as.character(subset_cells[[marker]]))
  
  # Generate the statistical test
  stat.test <- subset_cells %>%
    group_by(.data[[diagnosis_col]]) %>%
    t_test(as.formula(paste(marker, "~ Type"))) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = diagnosis_col)
  
  # Create the ggplot
  plot <- ggplot(subset_cells, aes(x = .data[[diagnosis_col]], y = .data[[marker]])) +
    geom_boxplot(aes(fill = class)) +
    labs(x = "Patient Diagnosis", y = paste(marker, "expression"), title = "Violin plot of marker expression in different molecular subtypes") +
    theme_pubr() +
    # scale_y_continuous(breaks = seq(0, 2.6, by = 0.5)) +
    stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0.008) +
    theme(legend.position = "top")
  
  print(plot)
}

# Example usage:
create_violin_plot("ANAX1", "Type")
  

table(cells$class[cells$proliferating == "Positive"])

table(cells$proliferating)
ca#----------------------------------------------------------------------------------------------------------#
#                                        Cell Count Plots                                                  #
#----------------------------------------------------------------------------------------------------------#


cell_counts <- cells[] %>%
  group_by(imageNb, Type, class, diagnosis, Histology_description,proliferating) %>%
  summarise(Count = n()) %>%
  ungroup()


ggplot(cell_counts, aes(x =Type, y = Count, fill = class)) +
  geom_bar(stat = "identity",position = "dodge") +
  labs(x = "Image Type", y = "Cell Count", fill = "Cell type") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Cell Counts",breaks = seq(0, 1000000, by = 100000) ,labels = comma)+
  ggtitle("Proportion of Cells in Each Image")


ggplot(cells, aes(x = Type, fill = class))+
  geom_bar(stat = "count",position = "dodge")


cell_counts <- cell_counts %>%
  group_by(Type) %>%
  mutate(Percentage = Count / sum(Count))

cell_counts$Percentage <- as.numeric(cell_counts$Percentage)


ggplot(cell_counts, aes(x = Type, y = Percentage, fill = class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 100),breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Image", y = "Percentage of Cell Types", fill = "Cell Class") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percentage proportion of cell types")


ggplot(cell_counts, aes(x = Type, y = Percentage, fill = class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Correct percentage scaling
  labs(x = "Diagnosis", y = "Percentage", fill = "Cell Type") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion of Cells in Each Diagnosis") +
  facet_wrap(~ class, scales = "free_y")
  

#-------------------------------------------------------------------------------------------------------#
create_violin_plot <- function(marker, cell_type) {
  plt1 <- ggplot(cells[(cells$class == cell_type & !(is.na(cells$Type))),],aes(x = Type, y = .data[[marker]],fill = Type))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    scale_fill_manual(values = c("gray","#e8e54f","#e84fe0"))+
    labs(x = "Type Classification", y = paste(marker," marker expression(a.u.)")) +
    scale_y_continuous(breaks = seq(0, 15, by = 1)) +
    theme_pubr()+
    stat_compare_means(method = "t.test",label = "p.signif", comparisons = list(c("Normal","Edge"),c("Normal","Tumour"),c("Edge","Tumour")))
  
  
  plt2 <-ggplot(cells[(cells$class == cell_type),],aes(x = diagnosis, y = .data[[marker]],fill = diagnosis))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    labs(x = "Patient Diagnosis", y = paste(marker," marker expression(a.u.)")) +
    scale_fill_manual(values = Diag)+
    scale_y_continuous(breaks = seq(0, max(cells[,marker]), by = (max(cells[,marker]) / 10))) +
    theme_pubr()+
    stat_compare_means(method = "anova",label = "p.signif", label.y = max(cells[,marker]),label.x = 0.5)
  
  plt3 <-ggplot(cells[(cells$class == cell_type),],aes(x = Histology_description, y = .data[[marker]],fill = Histology_description))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    labs(x = "Tumour Histology", y = paste(marker," marker expression(a.u.)")) +
    scale_fill_manual(values = Hist)+
    scale_y_continuous(breaks = seq(0, max(cells[,marker]), by = (max(cells[,marker])) / 10)) +
    theme_pubr()+
    stat_compare_means(method = "anova",label = "p.signif", label.y = max(cells[,marker]),label.x = 0.6)
  
  print(plt1)
  print(plt2)
  print(plt3)
}



for (marker in markers_columns) {
  create_violin_plot(marker,"Luminal")
}

for (marker in markers_columns) {
  create_violin_plot(marker,"Basal")
}


#----------------------------------------------------------------------------------------------------------#
#                                             Heatmaps                                                     #
#----------------------------------------------------------------------------------------------------------#


# Select relevant columns (adjust column names if needed)
markers <- cells %>% select(`Local Patient ID`, CD31, CK5, SMA, Ki67, CK8, ANAX1, ECAD, PCNT, CCASP3)

# Aggregate data: calculate the mean of each marker for each Local Patient ID
marker_means <- markers %>%
  group_by(`Local Patient ID`) %>%
  summarise(across(CD31:CCASP3, mean, na.rm = TRUE)) %>% 
  sort(`Local Patient ID`,decreaing = FALSE)
 

# Convert Local Patient ID to row names
row.names(marker_means) <- marker_means$`Local Patient ID`
marker_means <- marker_means %>% select(-`Local Patient ID`)

# Normalize or scale the data if needed
marker_means_scaled <- marker_means

# Generate the heatmap
pheatmap(marker_means_scaled, 
         main = "Summary Heatmap of Average Marker Values by Local Patient ID", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))









rnames <- cells[1:4000,"imageNb"]
mat_data <- data.matrix(cells[1:4000,markers_columns])
rnames <- make.unique(as.character(rnames))
rownames(mat_data) <- rnames
pheatmap(mat_data)


# Check the unique values in diagnosis and histology descriptions
unique_diagnosis <- unique(cells$diagnosis)
unique_histology <- unique(cells$Histology_description)

# Perform stratified sampling
sampled_data <- cells %>%
  group_by(diagnosis, Histology_description) %>%
  filter(n() > 1) %>% # Ensure groups with more than one sample
  sample_n(min(5, n())) %>% # Adjust the number '5' to ensure adequate sampling
  ungroup()

markers_scaled <- sampled_data[,markers_columns]

for (marker in colnames(markers_scaled)) {
  # Extract the marker data
  marker_data <- markers_scaled[, marker, drop = FALSE]
  
  # Generate the heatmap
  pheatmap(marker_data, 
           main = paste("Heatmap of", marker), 
           clustering_distance_rows = "euclidean", 
           clustering_distance_cols = "euclidean", 
           clustering_method = "complete", 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

#----------------------------------------------------------------------------------------------------------#
#                                           UMAP Plots                                                     #
#----------------------------------------------------------------------------------------------------------#
# 
# library(umap)
# 
# cells <- cells[cells$Type == "Normal",]
# 
# umap_df = cells[markers_columns] %>% umap(n_components = 2, n_neighbors = 15) # 2 and 15
# marker_kmeans <- kmeans(umap_df$data, centers = 8)
# umap_marker_df <- as.data.frame(umap_df$layout) %>% merge(cells, by = "row.names", all = TRUE)
# 
# 
# umap_marker_df <- umap_marker_df[markers_columns]
# marker_name <- colnames(umap_marker_df)
# marker_name <- markers_columns
# 
# 
# 
# umap_marker_df <- data.frame(umap_marker_df,marker_kmeans$cluster)
# umap_marker_df$marker_kmeans.cluster <- factor(umap_marker_df$marker_kmeans.cluster)
# 
# ggplot(data = umap_marker_df, aes(x = V1, y = V2, color = factor(marker_kmeans.cluster))) +
#   geom_point() +
#   theme_pubr() +  # Assuming you've defined the theme
#   labs(color = "Cluster") +  # Adding legend title
#   xlab("UMAP1") +  # Adding x-axis label
#   ylab("UMAP2") +  # Adding y-axis label
#   guides(color = guide_legend(title = "Cluster"))  # Adding legend title


# 
# Create SCE objects
cells <- cells[cells$Type == "Normal",]
df <- cells %>% filter_at(vars(Type), all_vars(!is.na(.)))
cells$fname <- NULL
counts <- cells[,markers_columns]
cell_meta <- cells[, !(names(cells) %in% markers_columns)]
sce <- SingleCellExperiment(assays = list(counts = t(counts)))
colData(sce) <- as(cell_meta, "DataFrame")



spe <- buildSpatialGraph(sce, img_id = "imageNb", type = "expansion", threshold = 30)


# 
# plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA NORMAL_A08.ome_4032_5376_4544_5888.tiff"],img_id = "imageNb", coords = c("Pos_Y","Pos_X"), colPairName = "expansion_interaction_graph", draw_edges = TRUE,
#             node_size_by = "area", node_color_by = "class")
# 
# cur_cells <- head(cells,n = 4000)
# 
# dittoHeatmap(sce[,cur_cells], 
#              genes = rownames(spe)[rowData(spe)$marker_class == "class"],
#              assay = "exprs", 
#              cluster_cols = FALSE, 
#              scale = "none",
#              heatmap.colors = viridis(100), 
#              annot.by = c("class", "diagnosis", "patient_id"),
#              annotation_colors = list(indication = metadata(spe)$color_vectors$diagnosis,
#                                       patient_id = metadata(spe)$color_vectors$patient_id,
#                                       celltype = metadata(spe)$color_vectors$`Local Patient ID`))

# 
# ggplot(cells[(cells$class == "Luminal" & !(is.na(cells$Type))),],aes(x = Type, y = ANAX1,fill = Type))+
#   geom_violin()+
#   geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
#   scale_fill_manual(values = c("gray","#e8e54f","#e84fe0"))+
#   labs(x = "Type Classification", y = "ANAX1 marker expression(a.u.)") +
#   scale_y_continuous(breaks = seq(0, 15, by = 1)) +
#   theme_pubr()+
#   stat_compare_means(method = "t.test",label = "p.signif", comparisons = list(c("Normal","Edge"),c("Normal","Tumour"),c("Edge","Tumour")))
# 
# 
# ggplot(cells[(cells$class == "Luminal"),],aes(x = diagnosis, y = ANAX1,fill = diagnosis))+
#   geom_violin()+
#   geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
#   labs(x = "Patient Diagnosis", y = "ANAX1 marker expression(a.u.)") +
#   scale_fill_manual(values = c("#1f78b4", "#33a02c", "#ff7f00", "#e31a1c"))+
#   scale_y_continuous(breaks = seq(0, 15, by = 1)) +
#   theme_pubr()+
#   stat_compare_means(method = "anova",label = "p.format", label.y = 11,label.x = 0.5)
# 
# ggplot(cells[(cells$class == "Luminal"),],aes(x = Histology_description, y = ANAX1,fill = Histology_description))+
#   geom_violin()+
#   geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
#   labs(x = "Tumour Histology", y = "ANAX1 marker expression(a.u.)") +
#   scale_fill_manual(values = c("#1f78b4", "#33a02c", "#ff7f00", "#e31a1c", "#6a3d9a"))+
#   scale_y_continuous(breaks = seq(0, 15, by = 1)) +
#   theme_pubr()+
#   stat_compare_means(method = "anova",label = "p.format", label.y = 11,label.x = 0.6)




create_violin_plot <- function(marker, cell_type) {
  
  # Filter data for the specific cell type and marker
  subset_data <- cells[cells$class == cell_type & !is.na(cells$Type), ]
  
  # Create a plotting dataframe with "Normal" added
  normal_data <- data.frame(Type = "Normal", marker_value = subset_data[[marker]])
  plotting_df <- rbind(subset_data, normal_data)
  
  plt1 <- ggplot(plotting_df, aes(x = Type, y = marker_value, fill = Type)) +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2) +
    scale_fill_manual(values = c("gray", "#e8e54f", "#e84fe0")) +
    labs(x = "Type Classification", y = paste(marker, "marker expression (a.u.)")) +
    scale_y_continuous(breaks = seq(0, 15, by = 1)) +
    theme_pubr() +
    stat_compare_means(method = "t.test", label = "p.signif", 
                       comparisons = list(c("Normal", "Edge"), c("Normal", "Tumour"), c("Edge", "Tumour")))
  
  plt2 <- ggplot(plotting_df, aes(x = diagnosis, y = marker_value, fill = diagnosis)) +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2) +
    labs(x = "Patient Diagnosis", y = paste(marker, "marker expression (a.u.)")) +
    scale_fill_manual(values = c("#1f78b4", "#33a02c", "#ff7f00", "#e31a1c", "#bdbdbd")) +
    scale_y_continuous(breaks = seq(0, 15, by = 1)) +
    theme_pubr() +
    stat_compare_means(method = "anova", label = "p.format", label.y = 11, label.x = 0.5)
  
  plt3 <- ggplot(plotting_df, aes(x = Histology_description, y = marker_value, fill = Histology_description)) +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2) +
    labs(x = "Tumour Histology", y = paste(marker, "marker expression (a.u.)")) +
    scale_fill_manual(values = c("#1f78b4", "#33a02c", "#ff7f00", "#e31a1c", "#6a3d9a", "#bdbdbd")) +
    scale_y_continuous(breaks = seq(0, 15, by = 1)) +
    theme_pubr() +
    stat_compare_means(method = "anova", label = "p.format", label.y = 11, label.x = 0.6)
  
  # Print or return the plots
  print(plt1)
  print(plt2)
  print(plt3)
}

create_violin_plot("ANAX1","Luminal")





plotting <- data.frame(cells[,markers_columns],cells$Type,cells$diagnosis,cells$Histology_description)
normal_df <- cells[cells$Type == "Normal",]

