library(ggplot2)
library(ggraph)
library(viridis)
library(tidyverse)
library(scales)
library(dplyr)
library(pals)
library(janitor)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)


raw_data <- read.csv("E:/[1]Datasets/[5]Quantification/gray_cell.csv")

data <- raw_data[raw_data$area >= 40,]

data <- data %>%
  mutate(Type = case_when(
    grepl("normal", patch_name, ignore.case = TRUE) ~ "Normal",
    grepl("edge", patch_name, ignore.case = TRUE) ~ "Edge",
    grepl("tumour", patch_name, ignore.case = TRUE) ~ "Tumour",
    grepl("Block 1", patch_name, ignore.case = TRUE) ~ "Block1",
    grepl("Block 2", patch_name, ignore.case = TRUE) ~ "Block2",
    grepl("Block 3", patch_name, ignore.case = TRUE) ~ "Block3",
    grepl("Block 5", patch_name, ignore.case = TRUE) ~ "Block5",
    grepl("TMA 1", patch_name, ignore.case = TRUE) ~ "TMA1",
    grepl("TMA 2", patch_name, ignore.case = TRUE) ~ "TMA2",
    grepl("TMA 3", patch_name, ignore.case = TRUE) ~ "TMA3",
    grepl("TMA 5", patch_name, ignore.case = TRUE) ~ "TMA5",
  ))


data <- data[data$Type == "Normal" | data$Type == "Edge" | data$Type == "Tumour",]
data$Type <- factor(data$Type,levels = c("Normal","Edge","Tumour"))
  
my_comparisons <- list(c("Normal","Edge"),c("Normal", "Tumour"),c("Edge","Tumour"))

ggplot(data, aes(x = Type, y = area, fill = Type)) +
  geom_violin()+
  geom_boxplot(width=.1, outlier.shape=NA) +
  # stat_summary(fun = "max", geom = "point", size = 3, fill = "black") +
  labs(title = "Distribution of Areas by Type",
       x = "Type",
       y = "Area") +
  scale_fill_manual(values = c("gray","#e8e54f","#e84fe0"))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
  theme_pubr()



# Calculate counts by type
count_by_type <- data %>%
  dplyr::count(Type)

# Calculate the number of unique images for each type
unique_images_by_type <- data %>%
  distinct(patch_name, Type) %>%
  dplyr::count(Type, name = "unique_images")

# Merge counts and unique image counts
count_by_type <- count_by_type %>%
  left_join(unique_images_by_type, by = "Type") %>%
  mutate(normalized_count = n / unique_images)

cell_counts <- data %>%
  group_by(patch_name) %>%
  summarise(cell_count = n())

# Merge the cell counts back to the original data
data <- merge(data, cell_counts, by = "patch_name")


# Plot normalized counts
ggplot(count_by_type, aes(x = Type, y = normalized_count, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Normalised Centrosome Count to number of images",
       x = "Type",
       y = "Average Centrosome Count") +
  coord_cartesian(ylim = c(0, 10))+
  scale_y_continuous(breaks = seq(0, 10, 1))+
  scale_fill_manual(values = c("gray","#e8e54f","#e84fe0"))+
  stat_compare_means(method = "t.test",label = "p.signif", comparisons = my_comparisons)+
  theme_pubr() +
  theme(axis.text.x = element_text())

my_comparisons <- list(c("Normal","Edge"),c("Normal", "Tumour"),c("Edge","Tumour"))

ggplot(data, aes(x = Type, y = cell_count, fill = Type)) +
  geom_violin()+
  geom_boxplot(width=.1, outlier.shape=NA) +
  labs(title = "Centrosome Count for each Sample",
       x = "Type",
       y = "Centrosome Count per Image") +
  scale_fill_manual(values = c("gray","#e8e54f","#e84fe0"))+
  scale_y_continuous(breaks = seq(0, 160, 10))+
  stat_compare_means(method = "t.test",label = "p.signif", comparisons = my_comparisons)+
  theme_pubr()
  # stat_compare_means(method = "anova",label = "p.signif",label.y = 160)


mean(data$area[data$Type == "Normal"])
mean(data$area[data$Type == "Tumour"])
