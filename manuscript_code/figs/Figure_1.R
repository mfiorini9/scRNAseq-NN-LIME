## This script was used to produce Figure 1.  

salloc -A def-sfarhan --time=0-8 -c 1 --mem=50g
module load StdEnv/2020 
module load r/4.2.2 
R


r_lib = '/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2'
library(Seurat)
library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(SeuratDisk)
library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyr)
library(ggpubr, lib="/home/fiorini9/scratch/mjff/MF/practice_ensemblux/R")
library(stringi)
library(stringr)


##################################################################
# Performance boxplots 
##################################################################
## Load in cell-type specific summary frames
# Astrocytes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Astrocytes_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astrocytes2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
Astrocytes_LR$eval_method[grepl("accuracy", Astrocytes_LR$metric )] <- "Accuracy"
Astrocytes_LR$eval_method[grepl("BA", Astrocytes_LR$metric )] <- "Balanced_Accuracy"
Astrocytes_LR$eval_method[grepl("F1", Astrocytes_LR$metric )] <- "F1_score"
Astrocytes_LR$ML_model <- "LR"
Astrocytes_LR$cell_type <- "Astrocytes"

Astrocytes_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astrocytes2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
Astrocytes_SVM$eval_method[grepl("accuracy", Astrocytes_SVM$metric )] <- "Accuracy"
Astrocytes_SVM$eval_method[grepl("BA", Astrocytes_SVM$metric )] <- "Balanced_Accuracy"
Astrocytes_SVM$eval_method[grepl("F1", Astrocytes_SVM$metric )] <- "F1_score"
Astrocytes_SVM$ML_model <- "SVM"
Astrocytes_SVM$cell_type <- "Astrocytes"

Astrocytes_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astrocytes2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
Astrocytes_RF$eval_method[grepl("accuracy", Astrocytes_RF$metric )] <- "Accuracy"
Astrocytes_RF$eval_method[grepl("BA", Astrocytes_RF$metric )] <- "Balanced_Accuracy"
Astrocytes_RF$eval_method[grepl("F1", Astrocytes_RF$metric )] <- "F1_score"
Astrocytes_RF$ML_model <- "RF"
Astrocytes_RF$cell_type <- "Astrocytes"

Astrocytes_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/astrocytes2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
Astrocytes_DNN$eval_method[grepl("accuracy", Astrocytes_DNN$metric )] <- "Accuracy"
Astrocytes_DNN$eval_method[grepl("BA", Astrocytes_DNN$metric )] <- "Balanced_Accuracy"
Astrocytes_DNN$eval_method[grepl("F1", Astrocytes_DNN$metric )] <- "F1_score"
Astrocytes_DNN$ML_model <- "DNN"
Astrocytes_DNN$cell_type <- "Astrocytes"

Astrocytes <- rbind(Astrocytes_LR, Astrocytes_SVM, Astrocytes_RF, Astrocytes_DNN)

# DaN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DaN_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/DaN2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
DaN_LR$eval_method[grepl("accuracy", DaN_LR$metric )] <- "Accuracy"
DaN_LR$eval_method[grepl("BA", DaN_LR$metric )] <- "Balanced_Accuracy"
DaN_LR$eval_method[grepl("F1", DaN_LR$metric )] <- "F1_score"
DaN_LR$ML_model <- "LR"
DaN_LR$cell_type <- "DaN"

DaN_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/DaN2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
DaN_SVM$eval_method[grepl("accuracy", DaN_SVM$metric )] <- "Accuracy"
DaN_SVM$eval_method[grepl("BA", DaN_SVM$metric )] <- "Balanced_Accuracy"
DaN_SVM$eval_method[grepl("F1", DaN_SVM$metric )] <- "F1_score"
DaN_SVM$ML_model <- "SVM"
DaN_SVM$cell_type <- "DaN"

DaN_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/DaN2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
DaN_RF$eval_method[grepl("accuracy", DaN_RF$metric )] <- "Accuracy"
DaN_RF$eval_method[grepl("BA", DaN_RF$metric )] <- "Balanced_Accuracy"
DaN_RF$eval_method[grepl("F1", DaN_RF$metric )] <- "F1_score"
DaN_RF$ML_model <- "RF"
DaN_RF$cell_type <- "DaN"

DaN_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/DaN2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
DaN_DNN$eval_method[grepl("accuracy", DaN_DNN$metric )] <- "Accuracy"
DaN_DNN$eval_method[grepl("BA", DaN_DNN$metric )] <- "Balanced_Accuracy"
DaN_DNN$eval_method[grepl("F1", DaN_DNN$metric )] <- "F1_score"
DaN_DNN$ML_model <- "DNN"
DaN_DNN$cell_type <- "DaN"

DaN <- rbind(DaN_LR, DaN_SVM, DaN_RF, DaN_DNN)

# Endothelial ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Endothelial_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Endothelial2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
Endothelial_LR$eval_method[grepl("accuracy", Endothelial_LR$metric )] <- "Accuracy"
Endothelial_LR$eval_method[grepl("BA", Endothelial_LR$metric )] <- "Balanced_Accuracy"
Endothelial_LR$eval_method[grepl("F1", Endothelial_LR$metric )] <- "F1_score"
Endothelial_LR$ML_model <- "LR"
Endothelial_LR$cell_type <- "Endothelial"

Endothelial_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Endothelial2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
Endothelial_SVM$eval_method[grepl("accuracy", Endothelial_SVM$metric )] <- "Accuracy"
Endothelial_SVM$eval_method[grepl("BA", Endothelial_SVM$metric )] <- "Balanced_Accuracy"
Endothelial_SVM$eval_method[grepl("F1", Endothelial_SVM$metric )] <- "F1_score"
Endothelial_SVM$ML_model <- "SVM"
Endothelial_SVM$cell_type <- "Endothelial"

Endothelial_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Endothelial2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
Endothelial_RF$eval_method[grepl("accuracy", Endothelial_RF$metric )] <- "Accuracy"
Endothelial_RF$eval_method[grepl("BA", Endothelial_RF$metric )] <- "Balanced_Accuracy"
Endothelial_RF$eval_method[grepl("F1", Endothelial_RF$metric )] <- "F1_score"
Endothelial_RF$ML_model <- "RF"
Endothelial_RF$cell_type <- "Endothelial"

Endothelial_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Endothelial2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
Endothelial_DNN$eval_method[grepl("accuracy", Endothelial_DNN$metric )] <- "Accuracy"
Endothelial_DNN$eval_method[grepl("BA", Endothelial_DNN$metric )] <- "Balanced_Accuracy"
Endothelial_DNN$eval_method[grepl("F1", Endothelial_DNN$metric )] <- "F1_score"
Endothelial_DNN$ML_model <- "DNN"
Endothelial_DNN$cell_type <- "Endothelial"

Endothelial <- rbind(Endothelial_LR, Endothelial_SVM, Endothelial_RF, Endothelial_DNN)

# Microglia ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Microglia_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Microglia2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
Microglia_LR$eval_method[grepl("accuracy", Microglia_LR$metric )] <- "Accuracy"
Microglia_LR$eval_method[grepl("BA", Microglia_LR$metric )] <- "Balanced_Accuracy"
Microglia_LR$eval_method[grepl("F1", Microglia_LR$metric )] <- "F1_score"
Microglia_LR$ML_model <- "LR"
Microglia_LR$cell_type <- "Microglia"

Microglia_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Microglia2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
Microglia_SVM$eval_method[grepl("accuracy", Microglia_SVM$metric )] <- "Accuracy"
Microglia_SVM$eval_method[grepl("BA", Microglia_SVM$metric )] <- "Balanced_Accuracy"
Microglia_SVM$eval_method[grepl("F1", Microglia_SVM$metric )] <- "F1_score"
Microglia_SVM$ML_model <- "SVM"
Microglia_SVM$cell_type <- "Microglia"

Microglia_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Microglia2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
Microglia_RF$eval_method[grepl("accuracy", Microglia_RF$metric )] <- "Accuracy"
Microglia_RF$eval_method[grepl("BA", Microglia_RF$metric )] <- "Balanced_Accuracy"
Microglia_RF$eval_method[grepl("F1", Microglia_RF$metric )] <- "F1_score"
Microglia_RF$ML_model <- "RF"
Microglia_RF$cell_type <- "Microglia"

Microglia_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Microglia2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
Microglia_DNN$eval_method[grepl("accuracy", Microglia_DNN$metric )] <- "Accuracy"
Microglia_DNN$eval_method[grepl("BA", Microglia_DNN$metric )] <- "Balanced_Accuracy"
Microglia_DNN$eval_method[grepl("F1", Microglia_DNN$metric )] <- "F1_score"
Microglia_DNN$ML_model <- "DNN"
Microglia_DNN$cell_type <- "Microglia"

Microglia <- rbind(Microglia_LR, Microglia_SVM, Microglia_RF, Microglia_DNN)

# OPC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OPC_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/OPC2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
OPC_LR$eval_method[grepl("accuracy", OPC_LR$metric )] <- "Accuracy"
OPC_LR$eval_method[grepl("BA", OPC_LR$metric )] <- "Balanced_Accuracy"
OPC_LR$eval_method[grepl("F1", OPC_LR$metric )] <- "F1_score"
OPC_LR$ML_model <- "LR"
OPC_LR$cell_type <- "OPC"

OPC_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/OPC2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
OPC_SVM$eval_method[grepl("accuracy", OPC_SVM$metric )] <- "Accuracy"
OPC_SVM$eval_method[grepl("BA", OPC_SVM$metric )] <- "Balanced_Accuracy"
OPC_SVM$eval_method[grepl("F1", OPC_SVM$metric )] <- "F1_score"
OPC_SVM$ML_model <- "SVM"
OPC_SVM$cell_type <- "OPC"

OPC_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/OPC2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
OPC_RF$eval_method[grepl("accuracy", OPC_RF$metric )] <- "Accuracy"
OPC_RF$eval_method[grepl("BA", OPC_RF$metric )] <- "Balanced_Accuracy"
OPC_RF$eval_method[grepl("F1", OPC_RF$metric )] <- "F1_score"
OPC_RF$ML_model <- "RF"
OPC_RF$cell_type <- "OPC"

OPC_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/OPC2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
OPC_DNN$eval_method[grepl("accuracy", OPC_DNN$metric )] <- "Accuracy"
OPC_DNN$eval_method[grepl("BA", OPC_DNN$metric )] <- "Balanced_Accuracy"
OPC_DNN$eval_method[grepl("F1", OPC_DNN$metric )] <- "F1_score"
OPC_DNN$ML_model <- "DNN"
OPC_DNN$cell_type <- "OPC"

OPC <- rbind(OPC_LR, OPC_SVM, OPC_RF, OPC_DNN)

# nDaN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nDaN_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/nDaN2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
nDaN_LR$eval_method[grepl("accuracy", nDaN_LR$metric )] <- "Accuracy"
nDaN_LR$eval_method[grepl("BA", nDaN_LR$metric )] <- "Balanced_Accuracy"
nDaN_LR$eval_method[grepl("F1", nDaN_LR$metric )] <- "F1_score"
nDaN_LR$ML_model <- "LR"
nDaN_LR$cell_type <- "nDaN"

nDaN_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/nDaN2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
nDaN_SVM$eval_method[grepl("accuracy", nDaN_SVM$metric )] <- "Accuracy"
nDaN_SVM$eval_method[grepl("BA", nDaN_SVM$metric )] <- "Balanced_Accuracy"
nDaN_SVM$eval_method[grepl("F1", nDaN_SVM$metric )] <- "F1_score"
nDaN_SVM$ML_model <- "SVM"
nDaN_SVM$cell_type <- "nDaN"

nDaN_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/nDaN2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
nDaN_RF$eval_method[grepl("accuracy", nDaN_RF$metric )] <- "Accuracy"
nDaN_RF$eval_method[grepl("BA", nDaN_RF$metric )] <- "Balanced_Accuracy"
nDaN_RF$eval_method[grepl("F1", nDaN_RF$metric )] <- "F1_score"
nDaN_RF$ML_model <- "RF"
nDaN_RF$cell_type <- "nDaN"

nDaN_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/nDaN2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
nDaN_DNN$eval_method[grepl("accuracy", nDaN_DNN$metric )] <- "Accuracy"
nDaN_DNN$eval_method[grepl("BA", nDaN_DNN$metric )] <- "Balanced_Accuracy"
nDaN_DNN$eval_method[grepl("F1", nDaN_DNN$metric )] <- "F1_score"
nDaN_DNN$ML_model <- "DNN"
nDaN_DNN$cell_type <- "nDaN"

nDaN <- rbind(nDaN_LR, nDaN_SVM, nDaN_RF, nDaN_DNN)

# all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/all2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
all_LR$eval_method[grepl("accuracy", all_LR$metric )] <- "Accuracy"
all_LR$eval_method[grepl("BA", all_LR$metric )] <- "Balanced_Accuracy"
all_LR$eval_method[grepl("F1", all_LR$metric )] <- "F1_score"
all_LR$ML_model <- "LR"
all_LR$cell_type <- "all"

all_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/all2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
all_SVM$eval_method[grepl("accuracy", all_SVM$metric )] <- "Accuracy"
all_SVM$eval_method[grepl("BA", all_SVM$metric )] <- "Balanced_Accuracy"
all_SVM$eval_method[grepl("F1", all_SVM$metric )] <- "F1_score"
all_SVM$ML_model <- "SVM"
all_SVM$cell_type <- "all"

all_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/all2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
all_RF$eval_method[grepl("accuracy", all_RF$metric )] <- "Accuracy"
all_RF$eval_method[grepl("BA", all_RF$metric )] <- "Balanced_Accuracy"
all_RF$eval_method[grepl("F1", all_RF$metric )] <- "F1_score"
all_RF$ML_model <- "RF"
all_RF$cell_type <- "all"

all_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/all2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
all_DNN$eval_method[grepl("accuracy", all_DNN$metric )] <- "Accuracy"
all_DNN$eval_method[grepl("BA", all_DNN$metric )] <- "Balanced_Accuracy"
all_DNN$eval_method[grepl("F1", all_DNN$metric )] <- "F1_score"
all_DNN$ML_model <- "DNN"
all_DNN$cell_type <- "all"

all <- rbind(all_LR, all_SVM, all_RF, all_DNN)                       
                                                
# Oligo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Oligo_LR <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Oligo2/LR_10_fold_summary_eval.csv", header = T, sep = ",")
Oligo_LR$eval_method[grepl("accuracy", Oligo_LR$metric )] <- "Accuracy"
Oligo_LR$eval_method[grepl("BA", Oligo_LR$metric )] <- "Balanced_Accuracy"
Oligo_LR$eval_method[grepl("F1", Oligo_LR$metric )] <- "F1_score"
Oligo_LR$ML_model <- "LR"
Oligo_LR$cell_type <- "Oligo"

Oligo_SVM <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Oligo2/SVM_10_fold_summary_eval.csv", header = T, sep = ",")
Oligo_SVM$eval_method[grepl("accuracy", Oligo_SVM$metric )] <- "Accuracy"
Oligo_SVM$eval_method[grepl("BA", Oligo_SVM$metric )] <- "Balanced_Accuracy"
Oligo_SVM$eval_method[grepl("F1", Oligo_SVM$metric )] <- "F1_score"
Oligo_SVM$ML_model <- "SVM"
Oligo_SVM$cell_type <- "Oligo"

Oligo_RF <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Oligo2/RF_10_fold_summary_eval.csv", header = T, sep = ",")
Oligo_RF$eval_method[grepl("accuracy", Oligo_RF$metric )] <- "Accuracy"
Oligo_RF$eval_method[grepl("BA", Oligo_RF$metric )] <- "Balanced_Accuracy"
Oligo_RF$eval_method[grepl("F1", Oligo_RF$metric )] <- "F1_score"
Oligo_RF$ML_model <- "RF"
Oligo_RF$cell_type <- "Oligo"

Oligo_DNN <- read.delim("/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/Oligo2/DNN_10_fold_summary_eval.csv", header = T, sep = ",")
Oligo_DNN$eval_method[grepl("accuracy", Oligo_DNN$metric )] <- "Accuracy"
Oligo_DNN$eval_method[grepl("BA", Oligo_DNN$metric )] <- "Balanced_Accuracy"
Oligo_DNN$eval_method[grepl("F1", Oligo_DNN$metric )] <- "F1_score"
Oligo_DNN$ML_model <- "DNN"
Oligo_DNN$cell_type <- "Oligo"

Oligo <- rbind(Oligo_LR, Oligo_SVM, Oligo_RF, Oligo_DNN)

## combine all dataframes
total_bind <- rbind(Astrocytes, DaN, Endothelial, Microglia, OPC, Oligo, all, nDaN )
total_bind_BA <- subset(total_bind, eval_method == "Balanced_Accuracy")

total_bind_BA$process[str_detect(total_bind_BA$method, "m1")] <- "Method 1"
total_bind_BA$process[str_detect(total_bind_BA$method, "m2")] <- "Method 2"
total_bind_BA$process[str_detect(total_bind_BA$method, "m3")] <- "Method 3"
total_bind_BA$process[str_detect(total_bind_BA$method, "m4")] <- "Method 4"

blups <- brewer.pal(11, "BrBG")

my_comparisons <- list( c("Method 1", "Method 2"), c("Method 1", "Method 3"), c("Method 1", "Method 4"),
                        c("Method 2", "Method 3"), c("Method 2", "Method 4"),
                        c("Method 3", "Method 4") )

ggplot(total_bind_BA, aes(x = process, y = X0, fill = process)) + 
geom_boxplot( outlier.size = .5, outliers = F, size = 0.4) +
theme_classic() +
theme(panel.grid = element_blank(),
axis.title = element_text(face = "bold"),
axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
axis.text.y = element_text(colour = "black"),
strip.background =element_rect(fill="white", colour = "white"),
legend.position = "none") +
ylab("Balanced accracy") + xlab("Feature selection\nmethod") +
scale_x_discrete(labels=c("HVGs", "PCs", "Factors",  "Topics")) +
scale_fill_manual(values=c("#F6E8C3", "#DFC27D", "#BF812D",  "#8C510A")) +
stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.001, hide.ns = T) # Add pairwise comparisons p-value
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf', sep=""), height = 3, width = 3) 

my_comparisons <- list( c("DNN", "LR"), c("DNN", "RF"), c("DNN", "SVM"), 
                            c("LR", "RF"), c("LR", "SVM"),
                            c("SVM", "RF") )

ggplot(total_bind_BA, aes(x = ML_model, y = X0)) + 
geom_boxplot( outlier.size = .5) +
theme_classic() +
theme(panel.grid = element_blank(),
axis.title = element_text(face = "bold"),
axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
axis.text.y = element_text(colour = "black"),
strip.background =element_rect(fill="white", colour = "white"),
legend.position = "bottom") +
ylab("Balanced accracy") + xlab("Machine learning\nclassifier") +
stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value
ggsave(paste('/home/fiorini9/scratch/PD_machine_lerning/figures/BA_ML_eval.pdf', sep=""), height = 5, width = 5) 


## compute mean
mean_results <- total_bind_BA %>%
group_by(method) %>%
summarise(mean_X0 = mean(X0, na.rm = TRUE))

## compute mean
mean_results <- total_bind_BA %>%
group_by(ML_model) %>%
summarise(mean_X0 = mean(X0, na.rm = TRUE))

##