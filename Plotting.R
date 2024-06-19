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
library(insight)



set.seed(123)
class_colors <- c(Luminal = "#51abcb", BnL = "#f8766d",Basal = "red", `Blood Vessel` = "#a3a500", LnB = "#00b0f6",Undefined = "grey")
Diag <- brewer.pal(4,'Pastel1')
Hist <- brewer.pal(5, 'Pastel2')

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
  scale_fill_manual("Pathology Classificaiton",values =c("gray","#f2edb6","#f2b6f2","lightblue"))+
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



ggplot(cells[cells$Type == "Edge",], aes(x = factor(AGE), y = (ANAX1)))+
  geom_boxplot()



#----------------------------------------------------------------------------------------------------------#
#                                        Classification                                                    #
#----------------------------------------------------------------------------------------------------------#

cells <- cells %>%
  mutate(class = case_when(
    (SMA >= 1) & (ECAD >= 1)~ "BnL",
    (SMA >= 1) & (CK8 >= 1)~ "LnB",
    CK8 >= 1 ~ "Luminal",
    ECAD >= 1 ~ "Luminal",
    SMA >= 1 ~ "Basal",
    CD31 >= 1 ~ "Blood Vessel",
    TRUE ~ "Undefined"
  ))

# cells <- cells %>%
#   mutate(class = case_when(
#     CK8 >= 1 ~ "Luminal",
#     ECAD >= 1~ "Luminal",
#     SMA >= 1 ~ "Basal",
#     CD31 >= 1 ~ "Blood Vessel",
#     TRUE ~ "Undefined"
#   ))

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
cells$class <- factor(cells$class, levels = c("Basal","BnL","Luminal","LnB","Blood Vessel","Undefined"))

ggplot(cells,aes(x = ECAD))+
  geom_density(stat = "bin",fill = "lightblue")+
  scale_y_continuous(name="Cell Counts",labels = comma)+
  scale_y_continuous(name="Cell Counts",breaks = seq(0, 110000, by = 10000) ,labels = comma)+
  theme_pubr()

ggplot(unfiltered_cells,aes(x = ECAD))+
  geom_density(stat = "bin",fill = "lightblue")+
  scale_y_continuous(name="Cell Counts",breaks = seq(0, 400000, by = 50000) ,labels = comma)+
  theme_pubr()

ggplot(unfiltered_cells,aes(x = CK8))+
  geom_density(stat = "bin",fill = "lightblue")+
  scale_y_continuous(name="Cell Counts",breaks = seq(0, 700000, by = 100000) ,labels = comma)+
  theme_pubr()

ggplot(cells,aes(x = CK8))+
  geom_density(stat = "bin",fill = "lightblue")+
  scale_y_continuous(name="Cell Counts",labels = comma)+
  theme_pubr()
#----------------------------------------------------------------------------------------------------------#
#                                           Summary Stats                                                  #
#----------------------------------------------------------------------------------------------------------#

type_summary_stats <- cells[!(is.na(cells$Type)),] %>%
  group_by(Type) %>%
  summarise(
    count = n(),
    mean_value = mean(ANAX1, na.rm = TRUE),
    median_value = median(ANAX1, na.rm = TRUE),
    sd_value = sd(ANAX1, na.rm = TRUE)
  )

histology_summary_stats <- cells %>%
  group_by(Histology_description) %>%
  summarise(
    count = n(),
    mean_value = mean(ANAX1, na.rm = TRUE),
    median_value = median(ANAX1, na.rm = TRUE),
    sd_value = sd(ANAX1, na.rm = TRUE)
  )

diagnosis_summary_stats <- cells %>%
  group_by(diagnosis) %>%
  summarise(
    count = n(),
    mean_value = mean(ANAX1, na.rm = TRUE),
    median_value = median(ANAX1, na.rm = TRUE),
    sd_value = sd(ANAX1, na.rm = TRUE)
  )

summary_stats <- type_summary_stats
differences <- data.frame(
  Comparison = c("Normal - Edge", "Normal - Tumour", "Edge - Tumour"),
  mean_difference = c(
    summary_stats$mean_value[1] - summary_stats$mean_value[2],
    summary_stats$mean_value[1] - summary_stats$mean_value[3],
    summary_stats$mean_value[2] - summary_stats$mean_value[3]
  ),
  median_difference = c(
    summary_stats$median_value[1] - summary_stats$median_value[2],
    summary_stats$median_value[1] - summary_stats$median_value[3],
    summary_stats$median_value[2] - summary_stats$median_value[3]
  ),
  sd_difference = c(
    summary_stats$sd_value[1] - summary_stats$sd_value[2],
    summary_stats$sd_value[1] - summary_stats$sd_value[3],
    summary_stats$sd_value[2] - summary_stats$sd_value[3]
  )
)

# Print the differences
print(differences)

# Function to calculate percentage change
percentage_change <- function(new, old) {
  ((new - old) / old) * 100
}

# Calculate percentage changes between types
percentage_changes <- data.frame(
  Comparison = c("Normal vs Edge", "Normal vs Tumour", "Edge vs Tumour"),
  mean_percentage_change = c(
    percentage_change(summary_stats$mean_value[2], summary_stats$mean_value[1]),
    percentage_change(summary_stats$mean_value[3], summary_stats$mean_value[1]),
    percentage_change(summary_stats$mean_value[3], summary_stats$mean_value[2])
  ),
  median_percentage_change = c(
    percentage_change(summary_stats$median_value[2], summary_stats$median_value[1]),
    percentage_change(summary_stats$median_value[3], summary_stats$median_value[1]),
    percentage_change(summary_stats$median_value[3], summary_stats$median_value[2])
  ),
  sd_percentage_change = c(
    percentage_change(summary_stats$sd_value[2], summary_stats$sd_value[1]),
    percentage_change(summary_stats$sd_value[3], summary_stats$sd_value[1]),
    percentage_change(summary_stats$sd_value[3], summary_stats$sd_value[2])
  )
)

# Print the percentage changes
print(percentage_changes)

#----------------------------------------------------------------------------------------------------------#
#                                        Marker expression Plots                                           #
#----------------------------------------------------------------------------------------------------------#
kruskal.test(ANAX1 ~ class, data = cells)

pairwise.wilcox.test(cells$weight, PlantGrowth$group,
                     p.adjust.method = "BH")

# Pairwise Wilcoxon test with Bonferroni adjustment and interaction
pairwise.wilcox.test(cells$ANAX1, cells$interaction_var, p.adjust.method = "bonferroni")

pairwise.wilcox.test(cells$ANAX1, cells$Histology_description, p.adjust.method = "bonferroni")


ggplot(cells[(cells$class == "Luminal")& !(is.na(cells$Type)),],aes(x = class, y = ANAX1,fill = class))+
  geom_violin()+
  geom_boxplot(width=.1, outlier.shape=NA,alpha = .2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  xlab(label = "Cancer Type")+
  ylab(label = "ANAX1 expression(a.u.)")+
  theme_pubr()+
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("Basal","BnL"),c("Basal","Luminal"),c("BnL","Luminal")))

ggplot(cells[(cells$class == "Luminal") & !(is.na(cells$Type)),],
       aes(x = Type, y = ANAX1, fill = Type)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("gray","#f2edb6","#f2b6f2"))+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  theme_pubr() +
  xlab(label = "Cancer Type")+
  ylab(label = "ANAX1 expression(a.u.)")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("Normal","Edge"),c("Normal","Tumour"),c("Edge","Tumour"))) +
  labs(title = "ANAX1 Expression by Type Faceted by Class")

ggplot(cells[(cells$class == "Luminal"),],aes(x = class, y = ANAX1,fill = class))+
  geom_violin()+
  geom_boxplot(width=.1, outlier.shape=NA,alpha = .2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  ylab(label = "ANAX1 expression(a.u.)")+
  theme_pubr()+
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("Basal","BnL"),c("Basal","Luminal"),c("BnL","Luminal")))

ggplot(cells[(cells$class == "Luminal"),],
       aes(x = diagnosis, y = ANAX1, fill = diagnosis)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = Diag)+
  ylab(label = "ANAX1 expression(a.u.)")+
  labs(fill='Diagnosis')+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  theme_pubr() +
  stat_compare_means(method = "kruskal.test",aes(label = paste0("p",scales::label_pvalue()(..p..),"****")),label.x = 0.7, label.y = 11)+
  labs(title = "ANAX1 Expression by Type Faceted by Class")


ggplot(cells[(cells$class == "Luminal"),],aes(x = class, y = ANAX1,fill = class))+
  geom_violin()+
  geom_boxplot(width=.1, outlier.shape=NA,alpha = .2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  theme_pubr()+
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("Basal","BnL"),c("Basal","Luminal"),c("BnL","Luminal")))

ggplot(cells[(cells$class == "Luminal"),],
       aes(x = Histology_description, y = ANAX1, fill = Histology_description)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = Hist)+
  ylab(label = "ANAX1 expression(a.u.)")+
  labs(fill='Histology')+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  theme_pubr() +
  stat_compare_means(method = "kruskal.test",aes(label = paste0("p",scales::label_pvalue()(..p..),"****")),label.x = 0.7, label.y = 11)+
  labs(title = "ANAX1 Expression by Type Faceted by Class")

ggplot(cells[(cells$class == "BnL" | cells$class == "Luminal" | cells$class == "Basal")& !(is.na(cells$Type)),],aes(x = class, y = ECAD,fill = class))+
  geom_violin()+
  geom_boxplot(width=.1, outlier.shape=NA,alpha = .2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  ylab(label = "ANAX1 expression(a.u.)")+
  theme_pubr()+
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("Basal","BnL"),c("Basal","Luminal"),c("BnL","Luminal")))+
  facet_wrap(~cells$Type[(cells$class == "BnL" | cells$class == "Luminal" | cells$class == "Basal")& !(is.na(cells$Type))])

ggplot(cells[(cells$class == "BnL" | cells$class == "Luminal" | cells$class == "Basal") & !(is.na(cells$Type)),],
       aes(x = Type, y = ANAX1, fill = Type)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("gray","#f2edb6","#f2b6f2"))+
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
ylab(label = "ANAX1 expression(a.u.)")+
  theme_pubr() +
  facet_wrap(~ class) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("Normal","Edge"),c("Normal","Tumour"),c("Edge","Tumour"))) +
  labs(title = "ANAX1 Expression by Type Faceted by Class")




#----------------------------------------------------------------------------------------------------------#  
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
    stat_compare_means(method = "kruskal.test",label = "p.signif", label.y = max(cells[,marker]),label.x = 0.5)
  
  plt3 <-ggplot(cells[(cells$class == cell_type),],aes(x = Histology_description, y = .data[[marker]],fill = Histology_description))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    labs(x = "Tumour Histology", y = paste(marker," marker expression(a.u.)")) +
    scale_fill_manual(values = Hist)+
    scale_y_continuous(breaks = seq(0, max(cells[,marker]), by = (max(cells[,marker])) / 10)) +
    theme_pubr()+
    stat_compare_means(method = "kruskal.test",label = "p.signif", label.y = max(cells[,marker]),label.x = 0.6)
  
  print(plt1)
  print(plt2)
  print(plt3)
}


create_violin_plot("ANAX1","BnL")
create_violin_plot("ANAX1","Luminal")

for (marker in markers_columns) {
  create_violin_plot(marker,"Luminal")
}


for (marker in markers_columns) {
  create_violin_plot(marker,"BnL")
}


for (marker in markers_columns) {
  create_violin_plot(marker,"LnB")
}


for (marker in markers_columns) {
  create_violin_plot(marker,"Basal")
}

for (marker in markers_columns) {
  create_violin_plot(marker,"Undefined")
}


#----------------------------------------------------------------------------------------------------------#
#                                        Cell Count Plots                                                  #
#----------------------------------------------------------------------------------------------------------#


# cell_counts <- cells[] %>%
#   group_by(imageNb, Type, class, diagnosis, Histology_description,proliferating) %>%
#   summarise(Count = n()) %>%
#   ungroup()
# 
# 
# ggplot(cell_counts, aes(x =Type, y = Count, fill = class)) +
#   geom_bar(stat = "identity",position = "dodge") +
#   labs(x = "Image Type", y = "Cell Count", fill = "Cell type") +
#   theme_pubr() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_y_continuous(name="Cell Counts",breaks = seq(0, 1000000, by = 100000) ,labels = comma)+
#   ggtitle("Proportion of Cells in Each Image")
# 
# 
# ggplot(cells, aes(x = Type, fill = class))+
#   geom_bar(stat = "count",position = "dodge")
# 
# 
# cell_counts <- cell_counts %>%
#   group_by(diagnosis) %>%
#   mutate(Percentage = Count / sum(Count))
# 
# cell_counts$Percentage <- as.numeric(cell_counts$Percentage)
# 
# 
# ggplot(cell_counts, aes(x = diagnosis, y = Percentage, fill = class)) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = scales::percent_format(scale = 100),breaks = seq(0, 1, by = 0.1)) +
#   labs(x = "Image", y = "Percentage of Cell Types", fill = "Cell Class") +
#   scale_fill_manual(fill())
#   theme_pubr() +
#   ggtitle("Percentage proportion of cell types")
# 
# 
# ggplot(cell_counts, aes(x = Type, y = Percentage, fill = class)) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Correct percentage scaling
#   labs(x = "Diagnosis", y = "Percentage", fill = "Cell Type") +
#   theme_pubr() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ggtitle("Proportion of Cells in Each Diagnosis") +
#   facet_wrap(~ class, scales = "free_y")
  

diagnosis_summary <- cells %>%
  group_by(diagnosis, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(diagnosis) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(diagnosis, class, percentage)

# Plot for Diagnosis
ggplot(diagnosis_summary, aes(x = diagnosis, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 2) +
  coord_cartesian(ylim = c(0, 100))+
  labs(title = "Percentage of Cell Types by Diagnosis",
       x = "Diagnosis",
       y = "Percentage",
       fill = "Class") +
  scale_y_continuous(name="Percentage of cells",breaks = seq(0, 100, by = 10) ,labels = function(x) paste0(x,"%"))+
  scale_fill_manual(values = class_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_pubr() 

histology_summary <- cells %>%
  group_by(Histology_description, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Histology_description) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(Histology_description, class, percentage)

# Plot for Histology
ggplot(histology_summary, aes(x = Histology_description, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 2) +
  coord_cartesian(ylim = c(0, 100))+
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(name="Percentage of cells",breaks = seq(0, 100, by = 10) ,labels = function(x) paste0(x,"%"))+
  labs(title = "Percentage of Cell Types by Histology",
       x = "Histology",
       y = "Percentage",
       fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_pubr() 

# Generate the summary table for type
type_summary <- cells[!(is.na(cells$Type)),] %>%
  group_by(Type, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Type) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(Type, class, percentage)

# Reshape data for chi-square test
mean_percentages_wide <- type_summary %>%
  pivot_wider(names_from = class, values_from = percentage)

# Perform chi-square test
chi_square_result <- chisq.test(mean_percentages_wide[, -1])  # Exclude Molecular_Subtype column

# Plot for Type
ggplot(type_summary, aes(x = Type, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 3) +
  coord_cartesian(ylim = c(0, 100))+
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(name="Percentage of cells",breaks = seq(0, 100, by = 10) ,labels = function(x) paste0(x,"%"))+
  labs(title = "Percentage of Cell Types by Type",
       x = "Type",
       y = "Percentage",
       fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = 0.9, y = 90,size = 3, label = insight::format_p(chi_square_result$p.value, stars = TRUE,))+
  theme_pubr() 

diagnosis_summary <- cells %>%
  group_by(diagnosis, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(diagnosis) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(diagnosis, class, percentage)

# Reshape data for chi-square test
diagnosis_summary_wide <- diagnosis_summary %>%
  pivot_wider(names_from = class, values_from = percentage)

# Perform chi-square test
chi_square_result_diagnosis <- chisq.test(diagnosis_summary_wide[, -1])  # Exclude diagnosis column

# Plot for Diagnosis
ggplot(diagnosis_summary, aes(x = diagnosis, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 2) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Percentage of Cell Types by Diagnosis",
       x = "Diagnosis",
       y = "Percentage",
       fill = "Class") +
  scale_y_continuous(name = "Percentage of cells", breaks = seq(0, 100, by = 10), labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = class_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = 0.9, y = 90, size = 3, label = insight::format_p(chi_square_result_diagnosis$p.value, stars = TRUE)) +
  theme_pubr()

histology_summary <- cells %>%
  group_by(Histology_description, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Histology_description) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(Histology_description, class, percentage)

# Reshape data for chi-square test
histology_summary_wide <- histology_summary %>%
  pivot_wider(names_from = class, values_from = percentage)

# Perform chi-square test
chi_square_result_histology <- chisq.test(histology_summary_wide[, -1])  # Exclude Histology_description column

# Plot for Histology
ggplot(histology_summary, aes(x = Histology_description, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 2) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Percentage of Cell Types by Histology",
       x = "Histology",
       y = "Percentage",
       fill = "Class") +
  scale_y_continuous(name = "Percentage of cells", breaks = seq(0, 100, by = 10), labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = class_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = 1, y = 90, size = 3, label = insight::format_p(chi_square_result_histology$p.value, stars = TRUE)) +
  theme_pubr()

#----------------------------------------------------------------------------------------------------------#
#                                             Proliferating                                               #
#----------------------------------------------------------------------------------------------------------#


proliferating_summary <- cells[!(is.na(cells$Type)),] %>%
  group_by(Type,proliferating,class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Type) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(Type, proliferating, percentage,class)

# Plot for Diagnosis
ggplot(proliferating_summary[proliferating_summary$proliferating == "Positive",], aes(x = Type, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 3) +
  coord_cartesian(ylim = c(0, 100))+
  labs(title = "Percentage of Proliferating cells by Diagnosis",
       x = "Diagnosis",
       y = "Percentage",
       fill = "Class") +
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(name="Percentage of cells",breaks = seq(0, 100, by = 10) ,labels = function(x) paste0(x,"%"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_pubr() 




#-------------------------------------------------------------------------------------------------------#http://127.0.0.1:21445/graphics/29b151d0-99ae-4286-b8ce-06e54833382d.png
#----------------------------------------------------------------------------------------------------------#
#                                             Heatmaps                                                     #
#----------------------------------------------------------------------------------------------------------#



anova_model <- aov(ANAX1 ~ class * Histology_description, data = cells)

posthoc_class <- glht(anova_model, linfct = mcp(class = "Tukey"))
posthoc_histology <- glht(anova_model, linfct = mcp(Histology_description = "Tukey"))

interaction_term <- interaction(cells$class, cells$Histology_description)
cells$interaction_term <- interaction_term

# Fit a new ANOVA model with the interaction term explicitly
anova_model_interaction <- aov(ANAX1 ~ interaction_term, data = cells)

# Post-hoc analysis for the interaction term
posthoc_interaction <- glht(anova_model_interaction, linfct = mcp(interaction_term = "Tukey"))

# Create a boxplot for visualization
ggplot(cells, aes(x = interaction(class, Histology_description), y = ANAX1, fill = class)) +
  geom_boxplot() +
  labs(title = "ANAX1 Expression by Cell Type and Histology", x = "Cell Type and Histology", y = "ANAX1 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

interaction.plot(cells$class, cells$Histology_description, cells$ANAX1,
                 col = c("red", "blue", "green", "purple", "orange"),
                 legend = TRUE, xlab = "Cell Type", ylab = "ANAX1 Expression",
                 main = "Interaction Plot of ANAX1 Expression by Cell Type and Histology")
ggplot(cells, aes(x = class, y = ANAX1, fill = class)) +
  geom_boxplot() +
  facet_wrap(~ Histology_description) +
  labs(title = "ANAX1 Expression by Cell Type and Histology", x = "Cell Type", y = "ANAX1 Expression") +
  theme_minimal()

# Calculate means for heatmap
mean_expression <- cells[cells$class == "Basal" | cells$class == "BnL" | cells$class == "Luminal",] %>%
  group_by(class, Type) %>%
  summarise(mean_ANAX1 = mean(ANAX1, na.rm = TRUE)) %>%
  spread(Type, mean_ANAX1)

# Convert to matrix for pheatmap
mean_expression_matrix <- as.matrix(mean_expression[, -1])
rownames(mean_expression_matrix) <- mean_expression$class

# Draw heatmap
pheatmap(mean_expression_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Heatmap of Mean ANAX1 Expression",
         color = colorRampPalette(c("blue", "white", "red"))(100))


mean_expressions <- cells %>%
  group_by(class) %>%
  summarise(across(markers_columns, mean, na.rm = TRUE)) %>%
  pivot_longer(cols = markers_columns, names_to = "Marker", values_to = "Mean_Expression")

# Reshape data to wide format for heatmap
mean_expressions_wide <- mean_expressions %>%
  pivot_wider(names_from = class, values_from = Mean_Expression)

# Convert to matrix for pheatmap
mean_expression_matrix <- as.matrix(mean_expressions_wide[,-1])
rownames(mean_expression_matrix) <- mean_expressions_wide$Marker

# Draw heatmap
pheatmap(mean_expression_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Heatmap of Mean Marker Expressions",
         color = inferno(100),
         fontsize_col = 10, angle_col = 45)


mean_expressions <- cells %>%
  group_by(Full_Filename) %>%
  summarise(across(markers_columns, \(x) mean(x, na.rm = TRUE)))



rowna <- mean_expressions[,1]
colna <- colnames(mean_expressions[markers_columns[!(markers_columns == "Nuclear")]])
mean_expressions_matrix <- as.matrix(mean_expressions[markers_columns[!(markers_columns == "Nuclear")]])

rownames(mean_expressions_matrix) <- rowna$Full_Filename
colnames(mean_expressions_matrix) <- colna

subset_patients <- rownames(mean_expressions_matrix)[1:100]
mean_expressions_subset <- mean_expressions_matrix[subset_patients, ]


pheatmap(mean_expressions_matrix,
         clustering_distance_rows = "euclidean",  # Distance metric for row clustering
         clustering_distance_cols = "euclidean",  # Distance metric for column clustering
         clustering_method = "complete",  # Clustering method
         row_dendrogram = TRUE,  # Display row dendrogram
         col_dendrogram = TRUE,  # Display column dendrogram
         scale = "row",  # Scale rows (optional)
         main = "Clustered Heatmap of Mean Marker Expressions by Patient",
         fontsize = 8,
         fontsize_row = 1,
         color = inferno(100)) 

pheatmap(mean_expressions_subset,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         row_dendrogram = TRUE,
         col_dendrogram = TRUE,
         scale = "row",
         main = "Clustered Heatmap of Mean Marker Expressions by Patient (Subset)",
         fontsize = 8,
         fontsize_row = 6,
         color = inferno(100))

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
subset_cells <- cells[cells$Full_Filename == "Mx_BEGIN TMA NORMAL_A08.ome.tif",]
# subset_cells <- cells[cells$Full_Filename == "Mx_BEGIN TMA TUMOUR_A11.ome.tif",]
# df <- cells %>% filter_at(vars(Type), all_vars(!is.na(.)))
cells$fname <- NULL
counts <- subset_cells[,markers_columns]
cell_meta <- subset_cells[, !(names(subset_cells) %in% markers_columns)]
sce <- SingleCellExperiment(assays = list(counts = t(counts)))
colData(sce) <- as(cell_meta, "DataFrame")



spe <- buildSpatialGraph(sce, img_id = "imageNb", type = "expansion", threshold = 30)


# 
plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA NORMAL_A08.ome_4928_4032_5440_4544.tiff"],img_id = "imageNb", coords = c("Pos_Y","Pos_X"), colPairName = "expansion_interaction_graph", draw_edges = FALSE,
            node_size_fix = 5, node_color_by = "class") + scale_color_manual(values = class_colors)



# plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA TUMOUR_A11.ome_4480_1792_4992_2304.tiff"],img_id = "imageNb", coords = c("Pos_Y","Pos_X"), colPairName = "expansion_interaction_graph", draw_edges = FALSE,
#             node_size_fix = 5, node_color_by = class_colors)
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



# for (marker in markers_columns){
#   label_title <- paste(marker, "expression")
#   print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
#           geom_boxplot(aes(fill = class))+
#           # geom_boxplot(width=.1, outlier.shape=NA) +
#           labs(x = "Patient Classification", y = label_title, title = "Scatter plot of Marker Correlations") +
#           stat_compare_means(method = "t.test",label = "p.signif", comparisons = c("Luminal","Basal"))+
#           facet_wrap(~cells$Type[(cells$class == "Luminal" | cells$class == "Basal")])+
#           theme_pubr()
#   )
# }
# 
# for (marker in markers_columns){
#   label_title <- paste(marker, "expression")
#   print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = class)) +
#           geom_boxplot(aes(fill = class))+
#           labs(x = "Patient Histology" , y = label_title, title = "Scatter plot of Marker Correlations") +
#           facet_wrap(~cells$Histology_description[(cells$class == "Luminal" | cells$class == "Basal")])+
#           theme_pubr()
#   )
# }
# 
# for (marker in markers_columns){
#   label_title <- paste(marker, "expression")
#   print(ggplot(cells[(cells$class == "Luminal" | cells$class == "Basal"),], aes(y = .data[[marker]], x = diagnosis)) +
#           geom_boxplot(aes(fill = class),outlier.shape=NA)+
#           labs(x = "Patient Diagnosis", y = label_title, title = "Violin plot of marker expression in different molecular subtypes") +
#           # facet_wrap(~cells$diagnosis[(cells$class == "Luminal" | cells$class == "Basal")])+
#           theme_pubr()
#   )
# }
# 
# stat.test <- cells[cells$class == "BnL" | cells$class == "Luminal",] %>%
#   group_by(diagnosis) %>%
#   t_test(ANAX1 ~ class) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance("p.adj") %>%
#   add_xy_position(x = "diagnosis", dodge = 1) 
# 
# ggplot(cells[(cells$class == "BnL" | cells$class == "Luminal"),], aes(x = diagnosis, y = ANAX1))+
#   geom_violin(aes(fill = class))+
#   theme_pubr()+
#   scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
#   stat_compare_means(method = "anova",label = "p.signif", label.y = 2.6, label.x = 2.45)+
#   stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0.008, y.position = 2.5)+
#   theme(legend.position = "top")
# 
# ggplot(cells[cells$class == "Basal" | cells$class == "Luminal",], aes(x = proliferating, y = ANAX1))+
#   geom_boxplot(aes(fill = class))+
#   theme_pubr()+
#   scale_y_continuous(breaks = seq(0, 2.6, by=0.5))+
#   theme(legend.position = "top")+
#   stat_compare_means(method = "anova",label = "p.signif")+
#   facet_wrap(~cells$Histology_description[cells$class == "Basal" | cells$class == "Luminal"])
# 
# 
# 
# 
# 
# # Plot violin plot
# ggplot(cells[!(cells$class == "Undefined"),], aes(x = class, y = ANAX1)) +
#   geom_violin(aes(fill = class)) +
#   theme_pubr() +
#   scale_y_continuous(breaks = seq(0, 10, by = 2)) +
#   theme(legend.position = "top")+
#   facet_wrap(~cells$Histology_description[!(cells$class == "Undefined")])
# 
# ggplot(cells[cells$class == "Luminal",], aes(x = Histology_description, y = ANAX1))+
#   geom_violin(fill = "#51abcb")+
#   geom_boxplot(width = .1,outlier.shape=NA)+
#   stat_compare_means(method = "anova", label = "p.signif", label.y = 2.6, label.x = 2.45) +
#   theme_pubr()


