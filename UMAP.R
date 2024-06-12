library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggcorrplot)
library("ggpubr")
library("cluster")
library("FactoMineR")
library("factoextra")  
library(reshape)
library(glue)
library(viridis)
library(ggvenn)
library("ggdendroplot")

set.seed(123)

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


# data <- read_csv("C:/Users/miles/GitHub/cellpose-quantification/cell_intensity_results_normalised_boxcox.csv")
data <- read_csv("E:/log/cell_intensity_results_normalised_tumour.csv")


markers_columns <- c("Nuclear","CD31","CK5","SMA","Ki67","CK8","ANAX1","ECAD","PCNT","CCASP3")


filtered_data <- data[data$area >= 200,]
table(rowSums(data[markers_columns] >= 0.2)==6)
filtered_data <- filtered_data[!(rowSums(filtered_data[markers_columns] >= 0.2)==6),]

filtered_data[order(filtered_data$CK8, decreasing = TRUE),]



ggplot(filtered_data, aes(x = Ki67, y = ANAX1))+
  geom_point()

classifcation <- filtered_data %>%
    mutate(Luminal = ifelse(CK8 >= 1.5, 1, 0))

classifcation <- classifcation %>%
  mutate(Basal = ifelse(SMA >= 1.5, 1, 0))

classifcation <- classifcation %>%
  mutate(BloodVessel = ifelse(CD31 >= 1.5, 1, 0))

length(which(classifcation$Luminal==1))

classifcation <- data.frame(as.logical(classifcation$Luminal), as.logical(classifcation$Basal), as.logical(classifcation$BloodVessel))
colnames(classifcation) <- c("Luminal","Basal","Blood Vessel")

count(classifcation,)

ggvenn(classifcation)

ggplot(filtered_data) +
  geom_point(aes(x = CK8, y = SMA, color = ifelse(classifcation$Luminal & classifcation$Basal, "green", ifelse(classifcation$Luminal, "blue", "red")))) +
  theme_pubr()





subset_dataset <- filtered_data[,markers_columns]


cell_size <- data[,"area"]

temp <- subset_dataset
temp$Nuclear <- NULL
temp$CK5 <- NULL

temp["Filename"] <- filtered_data["Filename"]


# Convert data to long format
melted_data <- temp %>%
  pivot_longer(-Filename, names_to = "Marker", values_to = "Expression")


clustered_data <- melted_data %>%
  mutate(Marker = factor(Marker, levels = unique(Marker))) %>%
  group_by(Marker) %>%
  mutate(rank = dense_rank(Expression)) %>%
  ungroup() %>%
  arrange(rank) %>%
  select(-rank)

# Create heatmap
ggplot(clustered_data, aes(x = Filename, y = Marker, fill = Expression)) +
  geom_tile() +
  labs(x = "Image", y = "Marker", title = "Heatmap of Marker Expressions Across Images") +
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(option = "C",name = "Expression")







s




# Create heatmap
ggplot(melted_data, aes(x = Filename, y = Marker, fill = Expression)) +
  geom_tile() +
  labs(x = "Image", y = "Marker", title = "Heatmap of Marker Expressions Across Images") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_viridis(option = "C",name = "Expression")

heatmap(as.matrix(subset_dataset),Colv = TRUE, Rowv = TRUE)


#box plot for each marker and image
for (marker in markers_columns){
  marker <- gsub("-", ".", marker)
  label_title <- paste(marker, "expression")
  print(ggplot(filtered_data, aes(x = filtered_data$Filename, y = .data[[marker]])) +
          geom_boxplot(aes(fill = filtered_data$Filename)) +
          labs(x = "Image", y = glue('{marker}'), title = "Boxplot of Marker Expressions Across Images") +
          theme_pubr()+
          # theme(axis.text.x = element_text(angle = 45, hjust = 1))
          theme(axis.text.x=element_blank())+
          scale_fill_discrete(name = "Core Name")
  )
  # filename <- paste(marker, "boxplot.png")
  # ggsave(path = "plots", filename = filename)
}
# 
# #box plot for each marker and image: LUMINAL
# for (marker in markers_columns){
#   marker <- gsub("-", ".", marker)
#   label_title <- paste(marker, "expression")
#   print(ggplot(filtered_data[classifcation$Luminal,], aes(x = filtered_data$Filename[classifcation$Luminal], y = .data[[marker]])) +
#           geom_boxplot(aes(fill = filtered_data$Filename[classifcation$Luminal])) +
#           labs(x = "Image", y = glue('{marker}'), title = "Boxplot of Marker Expressions Across Images") +
#           theme_pubr()+
#           # theme(axis.text.x = element_text(angle = 45, hjust = 1))
#           theme(axis.text.x=element_blank())+
#           scale_fill_discrete(name = "Core Name")
#   )
#   # filename <- paste(marker, "boxplot.png")
#   # ggsave(path = "plots", filename = filename)
# }
# 
# 
# #box plot for each marker and image: BASAL
# for (marker in markers_columns){
#   marker <- gsub("-", ".", marker)
#   label_title <- paste(marker, "expression")
#   print(ggplot(filtered_data[classifcation$Basal,], aes(x = filtered_data$Filename[classifcation$Basal], y = .data[[marker]])) +
#           geom_boxplot(aes(fill = filtered_data$Filename[classifcation$Basal])) +
#           labs(x = "Image", y = glue('{marker}'), title = "Boxplot of Marker Expressions Across Images") +
#           theme_pubr()+
#           # theme(axis.text.x = element_text(angle = 45, hjust = 1))
#           theme(axis.text.x=element_blank())+
#           scale_fill_discrete(name = "Core Name")
#   )
#   # filename <- paste(marker, "boxplot.png")
#   # ggsave(path = "plots", filename = filename)
# }
# 



# 
# Create a histogram for each marker
# for (marker in markers_columns) {
#   print(ggplot(data, aes(x = .data[[marker]])) +
#     geom_histogram(fill = "white", color = "black") +
#     labs(title = paste("Histogram of", marker),
#          x = marker,
#          y = "Frequency") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5)))
# }


pos <- filtered_data %>%
  mutate(Ki67_Positive = ifelse(Ki67 >= 1, "Positive", "Negative"))


ggplot(pos)+
  geom_point(aes(x=Ki67, y = ANAX1, color = classifcation$Luminal))+
  theme_pubr()

ggplot(pos)+
  geom_point(aes(x=Ki67, y = ANAX1, color = classifcation$Basal))+
  theme_pubr()

ggplot(pos)+
  geom_point(aes(x=Ki67, y  = ANAX1, color = Filename))+
  theme_pubr()

ggplot(pos[classifcation$Luminal,])+
  geom_boxplot(aes(x=Ki67_Positive, y  = ANAX1, fill = Filename))+
  theme_pubr()

table(cells_tumour$proliferating[cells_tumour$imageNb == "Mx_BEGIN TMA TUMOUR_A13.ome-1.tif"])
  

#--------------------------------------------------#
#system("pip install leidenalg python-igraph")
#library("leiden")

library(plotly)
library(umap)


sumap = subset_dataset %>% umap(n_components = 2, n_neighbors = 15) # 2 and 15
marker_kmeans <- kmeans(sumap$data, centers = 8)
umap_marker_df <- as.data.frame(sumap$layout) %>% merge(subset_dataset, by = "row.names", all = TRUE) 


umap_marker_df <- umap_marker_df[2:13]
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


ggplot(data = pos)+
  geom_boxplot(aes(x=Ki67_Positive, y=ANAX1,fill=Filename))+
  theme_pubr()

# 
# for (marker in marker_name) {
#   marker <- gsub("-", ".", marker)
#   label_title <- paste(marker, "expression")
#   print(
#     ggplot(
#       data = umap_marker_df,
#       aes(x = marker_kmeans.cluster, y = umap_marker_df[, marker], group = marker_kmeans.cluster)
#     ) +
#       geom_boxplot(
#         aes(color = factor(marker_kmeans.cluster), fill = factor(marker_kmeans.cluster))
#       ) +
#       theme_pubr() +  # Assuming you've defined the theme
#       ylab(label_title) +
#       xlab("Cluster") +  # Adding x-axis label
#       guides(color = guide_legend(title = "Cluster")) +  # Adding legend title
#       labs(fill = "Cluster")  # Adding fill legend title
#   )
# }


for(marker in marker_name){
  marker <- gsub("-", ".", marker)
  label_title <- paste(marker,"expression")
  print(ggplot(umap_marker_df, aes(x = V1,y = V2)) +
    geom_point(aes(color = .data[[marker]]), size = 3, alpha = 0.8)+
    theme_pubr()+
    theme(legend.position = "right",)+
      #scale_color_gradient(low="gray", high="blue")+
      scale_color_viridis(option = "C")+
      xlim(NA,10)+
      xlab("UMAP_1")+
      ylab("UMAP_2")+
      guides(color = guide_colorbar(title = label_title)))
    # filename <- paste(marker, "UMAP.png")
    # ggsave(path = "plots", filename = filename)
}


# 
# for (marker in marker_name) {
#   marker <- gsub("-", ".", marker)
#   label_title <- paste(marker, "expression")
#   print(
#     ggplot(umap_marker_df, aes(x = V1, y = V2)) +
#       geom_point(aes(color = .data[[marker]], shape = fname$Filename), size = 1, alpha = 0.8) +
#       theme_pubr() +
#       theme(legend.position = "right") +
#       #scale_color_gradient(low = "gray", high = "blue") +
#       scale_color_viridis(option = "C")+
#       xlim(NA, 10) +
#       xlab("UMAP_1") +
#       ylab("UMAP_2") +
#       guides(color = guide_colorbar(title = label_title, barheight = 15))
#   )
#   filename <- paste(marker, "UMAP_file.png")
#   ggsave(path = "plots", filename = filename)
# }



# Now try plotting again
ggplot(umap_marker_df, aes(x = V1, y = V2)) +
  geom_point(aes(color = filtered_data$Filename), size = 1, alpha = 0.8) +
  theme_pubr() +
  theme(legend.position = "right") +
  xlim(NA, 10) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  labs(color="Core Name")

# 
# # Step 1: Aggregate quantified marker expression to obtain image representations
# image_representations <- aggregate(. ~ Filename, data = data, FUN = mean)  # Assuming 'ImageID' is a column identifying each image
# 
# image_representations$Cell_ID <- NULL
# image_representations$`centroid-0` <- NULL
# image_representations$`centroid-1` <- NULL
# image_representations$eccentricity <- NULL
# image_representations$orientation <- NULL
# image_representations$area <- NULL
# image_representations$solidity <- NULL
# image_representations$perimeter <- NULL
# 
# image_representations["Normal"] <- grepl("NORMAL", image_representations$Filename, fixed = TRUE)
# image_representations["Edge"]<- grepl("EDGE", image_representations$Filename, fixed = TRUE)
# image_representations["Tumour"]<- grepl("TUMOUR", image_representations$Filename, fixed = TRUE)
# 
# # Step 2: Dimensionality Reduction using UMAP
# embedding <- umap(image_representations[, -1], n_neighbors = min(15, nrow(image_representations)))
# 
# # Step 3: Visualization
# embedding_df <- as.data.frame(embedding$layout)
# colnames(embedding_df) <- c("UMAP1", "UMAP2")
# ggplot(embedding_df, aes(x = UMAP1, y = UMAP2, color = image_representations$Filename)) +
#   geom_point() +
#   labs(title = "UMAP Visualization of Image Comparisons", x = "UMAP Dimension 1", y = "UMAP Dimension 2",color="Core Name")+
#   theme_pubr()
# 
# ggplot(embedding_df, aes(x = UMAP1, y = UMAP2)) +
#   geom_point(aes(color = factor(image_representations$Edge)), size = 3) +
#   labs(title = "UMAP Visualization of Image Comparisons", 
#        x = "UMAP Dimension 1", 
#        y = "UMAP Dimension 2") +
#   theme_pubr()
# 
# 
# x <- PCA(image_representations[,-1],scale.unit = TRUE, graph = TRUE)
# 
# 
# fviz_pca_var(x)
# 
# fviz_pca_ind(x,col.ind = image_representations$Normal)
# 
# 
# 
