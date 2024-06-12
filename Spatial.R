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


# data <- read_csv("C:/Users/miles/GitHub/cellpose-quantification/cell_intensity_results_normalised_boxcox.csv")
# cells <- read_csv("E:/[1]Datasets/[5]Quantification/log/joined.csv")
cells <- read_csv("E:/Final Data/quantification/joined_data.csv")
patient_data <- read_csv("E:/[1]Datasets/[5]Quantification/log/patient_metadata.csv")

# Define case statements
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

cells <- inner_join(cells, patient_data, by = "Full_Filename")


colnames(cells)[8] <- "ANAX1"
markers_columns <- c("Nuclear","CD31","CK5","SMA","Ki67","CK8","ANAX1","ECAD","PCNT","CCASP3")

cells[,markers_columns] <- cells[,markers_columns] + 1
cells[,markers_columns] <- log(cells[markers_columns],10)

names(cells)[names(cells) == "Filename"] <- "imageNb"
names(cells)[names(cells) == "centroid-0"] <- "Pos_X"
names(cells)[names(cells) == "centroid-1"] <- "Pos_Y"

names(cells)[names(cells) == "Filename"] <- "imageNb"
names(cells)[names(cells) == "centroid-0"] <- "Pos_Y"
names(cells)[names(cells) == "centroid-1"] <- "Pos_X"

# cells <- cells[cells$area >= 50,]
# table(rowSums(cells[markers_columns] >= 0.2)==6)
# cells <- cells[!(rowSums(cells[markers_columns] >= 0.2)==6),]


# Create classifications based on marker intensities & other metadata
cells <- cells %>%
  mutate(class = case_when(
    SMA >= 0.75 ~ "Basal",
    CK8 >= 0.75 ~ "Luminal",
    ECAD >= 0.75 ~ "Luminal",
    CD31 >= 0.75 ~ "Blood Vessel",
    TRUE ~ "Undefined"
  ))

y = cells[markers_columns]
summary(y[cells$class == "Undefined",])

ggplot(cells[cells$class == "Undefined",],aes(x = class,y = ECAD))+
  geom_violin()

cells <- cells %>% 
  mutate(proliferating = case_when(
    Ki67 >= 1 ~ "Positive",
    Ki67 < 1 ~ "Negative"
  ))


cells <- cells %>% 
  mutate(ANAX1_percent = (ANAX1 / max(ANAX1)) * 100)

cells <- cells %>%
  mutate(Type = case_when(
    grepl("normal", imageNb, ignore.case = TRUE) ~ "Normal",
    grepl("edge", imageNb, ignore.case = TRUE) ~ "Edge",
    grepl("tumour", imageNb, ignore.case = TRUE) ~ "Tumour",
  ))



cells$Type <- factor(cells$Type, levels = c("Normal", "Edge", "Tumour"))
cells$diagnosis <- factor(cells$diagnosis, levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
cells$Histology_description <- factor(cells$Histology_description, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))

# Subset data for each treatment condition
cells_normal <- cells[grepl("NORMAL", cells$imageNb, ignore.case = TRUE), ]
cells_edge <- cells[grepl("EDGE", cells$imageNb, ignore.case = TRUE), ]
cells_tumour <- cells[grepl("TUMOUR", cells$imageNb, ignore.case = TRUE), ]

# 
# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from=plots.png.paths, to="C:/Users/miles/GitHub/Figures")

#----------------------------- Correlations ----------------------------------------------------------#
# Scatter plots showing coorelation and line of best fit for normal, edge, and tumour as well as classification
library('Hmisc')
subset_cells <- cells[markers_columns]
# subset_cells <- subset_cells[cells$class == "Basal",]
res2 <- rcorr(as.matrix(subset_cells),type = c("pearson"))

# Insignificant correlation are crossed
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

# Insignificant correlations are leaved blank
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

heatmap(res2$r)

################################ PCNT vs ANAX1 ################################ 
ggplot(cells_normal[cells_normal$class %in% c("Luminal", "Basal"),], aes(x = PCNT, y = ANAX1, color = class)) +
  geom_point(size = 0.5) +
  theme_pubr() +
  xlab(label = "PCNT expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")

ggplot(cells_edge[cells_edge$class %in% c("Luminal", "Basal"),], aes(x = PCNT, y = ANAX1, color = class)) +
  geom_point(size = 0.5) +
  theme_pubr() +
  xlab(label = "PCNT expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")

ggplot(cells_tumour[cells_tumour$class %in% c("Luminal", "Basal"),], aes(x = PCNT, y = ANAX1, color = class)) +
  geom_point(size = 0.5) +
  theme_pubr() +
  xlab(label = "PCNT expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")


################################ Ki67 vs ANAX1 ################################ 
ggplot(cells_normal[cells_normal$class %in% c("Luminal", "Basal"),], aes(x = Ki67, y = ANAX1, color = class)) +
  geom_point(size = 1) +
  theme_pubr() +
  xlab(label = "Ki67 expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")

ggplot(cells_edge[cells_edge$class %in% c("Luminal", "Basal"),], aes(x = Ki67, y = ANAX1, color = class)) +
  geom_point(size = 1) +
  theme_pubr() +
  xlab(label = "Ki67 expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")

ggplot(cells_tumour[cells_tumour$class == "Basal" | cells_tumour$class == "Luminal",], aes(x = Ki67, y = ANAX1, color = class)) +
  geom_point(size = 1) +
  theme_pubr() +
  xlab(label = "Ki67 expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")



ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(x = Ki67, y = ANAX1, color = class)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm",se=F)+
  xlab(label = "Ki67 expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom",)+
  facet_wrap(~cells$diagnosis[(cells$class == "Luminal" | cells$class == "Basal")])+
  theme_pubr()

ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(x = PCNT, y = ANAX1, color = class)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm",se=F)+
  xlab(label = "Ki67 expression (a.u.)")+
  ylab(label = "ANAX1 expression (a.u.)")+
  stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")+
  facet_wrap(~cells$Type[(cells$class == "Luminal" | cells$class == "Basal")])+
  theme_pubr()

# for (marker in markers_columns){
#   label_title <- paste(marker, "expression")
#   print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(x = .data[[marker]], y = ANAX1, color = class)) +
#           geom_point() +
#           labs(x = label_title, y = "ANAX1 expresssion", title = "Scatter plot of Marker Correlations") +
#           stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")+
#           facet_wrap(~cells$Type[(cells$class == "Luminal" | cells$class == "Basal")])+
#           geom_smooth(method = "lm",se=F)+
#           theme_pubr()
#           # theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   )
# }
# 
# for (marker in markers_columns){
#   label_title <- paste(marker, "expression")
#   print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(x = .data[[marker]], y = ANAX1, color = class)) +
#           geom_point() +
#           labs(x = label_title, y = "ANAX1 expresssion", title = "Scatter plot of Marker Correlations") +
#           stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")+
#           facet_wrap(~cells$diagnosis[(cells$class == "Luminal" | cells$class == "Basal")])+
#           geom_smooth(method = "lm",se=F)+
#           theme_pubr()
#   )
# }

# for (marker in markers_columns){
#   label_title <- paste(marker, "expression")
#   print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(x = .data[[marker]], y = ANAX1, color = class)) +
#           geom_point() +
#           labs(x = label_title, y = "ANAX1 expresssion", title = "Scatter plot of Marker Correlations") +
#           stat_cor(method = "pearson",label.x.npc = "middle",label.y.npc = "bottom")+
#           facet_wrap(~cells$Histology_description[(cells$class == "Luminal" | cells$class == "Basal")])+
#           geom_smooth(method = "lm",se=F)+
#           theme_pubr()
#   )
# }

for (marker in markers_columns){
  label_title <- paste(marker, "expression")
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
          geom_violin(aes(fill = class))+
          geom_boxplot(width=.1, outlier.shape=NA) +
          labs(x = "Patient Classification", y = label_title, title = "Scatter plot of Marker Correlations") +
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
  print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
          geom_violin(aes(fill = class))+
          geom_boxplot(width=.1, outlier.shape=NA) +
          labs(x = "Patient Diagnosis", y = label_title, title = "Violin plot of marker expression in different molecular subtypes") +
          facet_wrap(~cells$diagnosis[(cells$class == "Luminal" | cells$class == "Basal")])+
          theme_pubr()
  )
}




#----------------------------- Cell Counts ----------------------------------------------------------#

# Count the number of cells for each cell type in each image
cell_counts <- cells %>%
  group_by(imageNb, Type, class, diagnosis, ANAX1,PCNT,Ki67,proliferating) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(Type = factor(Type, levels = c("Normal", "Edge", "Tumour")))


# Plot the cell counts for each cell type
ggplot(cell_counts[cell_counts$class == "Basal" | cell_counts$class == "Luminal",], aes(x = Type, y = Count , fill = class)) +
  geom_bar(stat = "identity")+
  labs(x = "Cell Type", y = "Ki67 expression", fill = "Classification") +
  theme_minimal()


#Plots the cell count as percentage
cell_counts <- cell_counts %>%
  group_by(Type) %>%
  mutate(Percentage = Count / sum(Count))

cell_counts$Percentage <- as.numeric(cell_counts$Percentage)

ggplot(cell_counts, aes(x = Type, y = Count, fill = class)) +
  geom_bar(stat = "identity") +
  labs(x = "Image", y = "Count", fill = "Cell Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion of Cells in Each Image")


ggplot(cell_counts, aes(x = Type, y = Percentage, fill = class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "Image", y = "Percentage", fill = "Cell Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion of Cells in Each Image")


subtypes_cell_counts <- cells %>%
  group_by(imageNb, class,diagnosis, HISTOLOGY,Type,ANAX1,Ki67) %>%
  summarise(Count = n()) %>%
  ungroup()%>%
  mutate(Type = factor(Type, levels = c("Normal", "Edge", "Tumour")))


#Plots the cell count as percentage
subtypes_cell_counts <- subtypes_cell_counts %>%
  group_by(diagnosis) %>%
  mutate(Percentage = Count / sum(Count))

hist_labels <- c("Invasive ductal","lobular","Invasive mucinous","Mixed")

# Plot the stacked bar chart
ggplot(subtypes_cell_counts, aes(x = as.factor(diagnosis), y = Percentage, fill = class)) +
  geom_col() +
  labs(x = "Image Number", y = "Percentage (%)", fill = "Classification") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab(label = "Histology")



stat.test <- cells[cells$class == "Basal" | cells$class == "Luminal",] %>%
  group_by(Type) %>%
  t_test(ANAX1 ~ class) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Type", dodge = 1) 


comparisons <- list(c("Normal","Edge"),c("Normal","Tumour"),c("Edge","Tumour"))

ggplot(cells[(cells$class == "Basal" | cells$class == "Luminal"),], aes(x = Type, y = ANAX1))+
  geom_violin(aes(fill = class))+
  theme_pubr()+
  scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
  stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0.008, y.position = 2.5)+
  theme(legend.position = "top")+
  stat_compare_means(method = "anova",label = "p.signif", label.x = 1.98, label.y = 2.6)

 
ggplot(cells[cells$class == "Basal" | cells$class == "Luminal",], aes(x = proliferating, y = ANAX1))+
  geom_boxplot(aes(fill = class))+
  theme_pubr()+
  scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
  theme(legend.position = "top")+
  stat_compare_means(method = "anova",label = "p.signif")+
  facet_wrap(~cells$diagnosis[cells$class == "Basal" | cells$class == "Luminal"])
  # stat_compare_means(method = "t.test",label = "p.signif", label.x = 1, label.y = 2.6)

ggplot(cells[cells$class == "Basal" | cells$class == "Luminal",], aes(x = proliferating, y = ANAX1))+
  geom_boxplot(aes(fill = class))+
  theme_pubr()+
  scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
  theme(legend.position = "top")+
  stat_compare_means(method = "anova",label = "p.signif")+
  facet_wrap(~cells$Histology_description[cells$class == "Basal" | cells$class == "Luminal"])
  
ggplot(cells[cells$class == "Basal" | cells$class == "Luminal",], aes(x = proliferating, y = ANAX1))+
  geom_boxplot(aes(fill = class))+
  theme_pubr()+
  scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
  theme(legend.position = "top")+
  facet_wrap(~cells$Type[cells$class == "Basal" | cells$class == "Luminal"])+
  stat_compare_means(method = "anova",label = "p.signif")


table(cells$diagnosis)

ggplot(cells,aes(x = diagnosis, y  = class))+
  geom_jitter()

cells[cells$proliferating == "Positive",]

table(cells$proliferating[cells$imageNb == "Mx_BEGIN TMA TUMOUR_A13.ome-1.tif"])

table((cells$proliferating == "Positive"))

cells[cells$class == "Undefined",]

# cell_counts <- cells %>% 
#   group_by(imageNb,class) %>% 
#   summarise(Count = n()) %>% 
#   
#   ggplot(cell_counts, aes(x = imageNb, y = Count, fill = class)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Image", y = "Cell Count", fill = "Classification") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# Create SCE objects
counts <- cells[,markers_columns]
cell_meta <- cells[, !(names(cells) %in% markers_columns)]
sce <- SingleCellExperiment(assays = list(counts = t(counts)))
colData(sce) <- as(cell_meta, "DataFrame")


spe <- buildSpatialGraph(sce, img_id = "imageNb", type = "expansion", threshold = 30)



plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA NORMAL_I16.ome.tif"],img_id = "imageNb", coords = c("Pos_X","Pos_Y"), colPairName = "expansion_interaction_graph", draw_edges = TRUE, 
            node_size_by = "area", node_shape_by = "class", node_color_by = "ANAX1_percent") + scale_color_viridis(option = "C")

for (fname in unique(spe$imageNb)) {
  print(plotSpatial(spe[,spe$imageNb == fname],img_id = "imageNb", coords = c("Pos_X","Pos_Y"), 
                    colPairName = "expansion_interaction_graph", node_size_by = "ANAX1_percent",node_shape_by = "proliferating",node_color_by = "class")
                    +scale_color_manual(values=class_colors)
        )
}


spe <- aggregateNeighbors(spe, colPairName = "expansion_interaction_graph", 
                          aggregate_by = "metadata", count_by = "class")

cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 15)
spe$cn_celltypes <- as.factor(cn_1$cluster)


# Heatmap of cell types in each neighborhood for PDG
for_plot <- colData(spe) %>% as_tibble() %>%
  group_by(cn_celltypes, class) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = cn_celltypes, names_from = class, 
              values_from = freq, values_fill = 0) %>%
  ungroup() %>%
  select(-cn_celltypes)

pheatmap(for_plot,
         color=colorRampPalette(c("darkblue","white","darkred"))(100), 
         scale = "column")


# Get count of cn_clust cells per image (PDG)
cn_count <- as.data.frame(cbind(spe$imageNb, spe$cn_celltypes))

cn_count2 <- cn_count %>%
  group_by(V1,V2) %>%
  add_count(V1,V2) %>%
  distinct(V1, V2, .keep_all = TRUE) %>%
  arrange(desc(V1), V2)

cn_count3 <- cn_count2
cn_count3$Type <- NA

cn_count3$Type[grepl("NORMAL", tolower(cn_count3$V1))] <- "NORMAL"
cn_count3$Type[grepl("EDGE", tolower(cn_count3$V1))] <- "EDGE"
cn_count3$Type[grepl("TUMOUR", tolower(cn_count3$V1))] <- "TUMOUR"

cn_count3 <- cn_count3 %>%
  pivot_wider(names_from = V2, values_from = n, values_fill = 0)


# Subset data for each treatment condition
cells_normal <- cells[grepl("NORMAL", cells$imageNb, ignore.case = TRUE), ]
cells_edge <- cells[grepl("EDGE", cells$imageNb, ignore.case = TRUE), ]
cells_tumour <- cells[grepl("TUMOUR", cells$imageNb, ignore.case = TRUE), ]

# Create SCE objects for eachtreatment condition
counts1 <- cells_normal[,markers_columns]
cell_meta1 <- cells_normal[, !(names(cells_normal) %in% markers_columns)]
sce_normal <- SingleCellExperiment(assays = list(counts1 = t(counts1)))
colData(sce_normal) <- as(cell_meta1, "DataFrame")

counts2 <- cells_edge[,markers_columns]
cell_meta2 <- cells_edge[, !(names(cells_edge) %in% markers_columns)]
sce_edge <- SingleCellExperiment(assays = list(counts2 = t(counts2)))
colData(sce_edge) <- as(cell_meta2, "DataFrame")

counts3 <- cells_tumour[,markers_columns]
cell_meta3 <- cells_tumour[, !(names(cells_tumour) %in% markers_columns)]
sce_tumour <- SingleCellExperiment(assays = list(counts3 = t(counts3)))
colData(sce_tumour) <- as(cell_meta3, "DataFrame")



# Create spatial interaction graph per image for each treatment condition using 'expansion' method set to 30 micron diamter. This step can take hours to a full day

spe_normal <- buildSpatialGraph(sce_normal, img_id = "imageNb", type = "expansion", threshold = 30)
spe_edge <- buildSpatialGraph(sce_edge, img_id = "imageNb", type = "expansion", threshold = 30)
spe_tumour <- buildSpatialGraph(sce_tumour, img_id = "imageNb", type = "expansion", threshold = 30)

# 'Classic' method used to test interactions between cell types. Can take multiple days to complete each of the following steps
out_normal <- testInteractions(spe_normal, 
                           group_by = "imageNb",
                           label = "class", 
                           colPairName = "expansion_interaction_graph",
                           BPPARAM = SerialParam(RNGseed = 123))
head(out_normal)

out_edge <- testInteractions(spe_edge, 
                           group_by = "imageNb",
                           label = "class", 
                           colPairName = "expansion_interaction_graph",
                           iter = 1000,
                           BPPARAM = SerialParam(RNGseed = 123))
head(out_edge)

out_tumour <- testInteractions(spe_tumour, 
                           group_by = "imageNb",
                           label = "class", 
                           colPairName = "expansion_interaction_graph",
                           iter = 1000,
                           BPPARAM = SerialParam(RNGseed = 123))
head(out_tumour)

out_normal %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_edge %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_tumour %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





