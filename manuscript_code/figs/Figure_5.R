## This script was used to produce Figure 5.

salloc -A def-sfarhan --time=0-8 -c 1 --mem=40g

module load StdEnv/2020 
module load r/4.2.2 
R

library(Seurat)
library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(rstatix)
library(inflection, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(SeuratData, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(viridis)
library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(circlize)
library(MAST, lib ='/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2')
library(robustbase)
library(purrr)
library(GenomeInfoDb, lib ='/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2')
library(GenomicRanges, lib ='/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2')
library(SummarizedExperiment, lib ='/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2')
library(SingleCellExperiment, lib ='/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2')
library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2" )

############################################################################ Summary heatmap
############################################################################
############################################################################
############################################################################

######################################## Get genes
########################################
########################################
########################################
######################################## 
celltypelist <- c("astro", "endo", "mg", "olig", "opc", "peri", "nonda", "da")

feature <- "fill"
celltype<- "fill"
mean_z_score_weighted    <- "fill"

fill <- data.frame(feature,celltype, mean_z_score_weighted)
for (i in celltypelist){
meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/table_genes_',i,'.csv'), header = T, sep = ",")
meta_data <- meta_data %>% dplyr::select(feature, celltype, mean_z_score_weighted)
fill <- rbind(fill, meta_data)
}

fill <- subset(fill, celltype != "fill")
length(unique(fill$feature))
genes <- unique(fill$feature)

filtered_genes <- genes[!grepl("^MT-", genes)]
filtered_genes <- filtered_genes[!grepl("^RP", filtered_genes)]

freq_df <- data.frame(table(fill$feature))
freq_df <- subset(freq_df, Var1 %in% filtered_genes)


######################################## Z-score
########################################
########################################
########################################
######################################## 
celltypelist <- c("astro", "endo", "mg", "olig", "opc", "peri", "nonda")

feature <- "fill"
celltype<- "fill"
mean_z_score_weighted    <- "fill"

fill <- data.frame(feature,celltype, mean_z_score_weighted)
for (i in celltypelist){
meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/New_focused_', i, '_lim_FC.csv'), header = T, sep = ",")
meta_data$celltype <- i
meta_data <- meta_data %>% dplyr::select(feature, celltype, mean_z_score_weighted)
fill <- rbind(fill, meta_data)
}

fill <- subset(fill, celltype != "fill")
fill$mean_z_score_weighted <- as.numeric(fill$mean_z_score_weighted)

## add Da
meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/New_focused_da_lim_FC.csv'), header = T, sep = ",")
meta_data$celltype <- "da"
meta_data <- meta_data %>% dplyr::select(feature, celltype, mean_z_score_weighted)
fill <- rbind(fill, meta_data)

unique(fill$celltype)

## keep genes that meat the criteria
fill_lim <- subset(fill, feature %in% filtered_genes)

## Set cell type factor levels
fill_lim$celltype <- factor(fill_lim$celltype, levels = c("astro", "endo", "mg", "olig",  "opc", "peri", "nonda", "da"))

fill_lim2 <- merge(fill_lim, freq_df, by. = "feature", by.y = "Var1" )
nrow(fill_lim) == nrow(fill_lim2)

# Custom order for celltype
celltype_order <- c("astro", "endo", "mg", "olig", "opc", "peri", "nonda", "da")

# Order the celltype factor levels based on the custom order
fill_lim2$celltype <- factor(fill_lim2$celltype, levels = celltype_order)

## add real frequency column
fill_lim2 <- fill_lim2 %>%
    group_by(feature) %>%
    mutate(Frequency = n()) %>%
    ungroup()
## add mean across all cell types
fill_lim2 <- fill_lim2 %>%
    group_by(feature) %>%
    mutate(mean = mean(mean_z_score_weighted)) %>%
    ungroup()

# Create a new feature factor ordered by Freq, then by celltype
fill_lim2 <- fill_lim2 %>%
arrange(desc(Frequency), celltype, desc(mean)) %>%  
mutate(feature = factor(feature, levels = unique(feature)))  

## Set color scheme
blups <- brewer.pal(9, "Oranges")

## Plot
LIME <- ggplot(fill_lim2, aes(y = feature, x = celltype, fill = as.numeric(mean_z_score_weighted))) +
geom_tile() +
theme_bw() +
theme(panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "white", colour = "white"),
    axis.text.y = element_text( face = "italic", colour = c("black")),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = c("black"))) +
labs(fill = "") +
facet_wrap(~"LIME Z-score") +
scale_fill_gradientn(colors = blups, na.value = "grey") +
scale_x_discrete(expand = c(0,0), labels = c("astro" = "Astrocytes", "all" = "All cells", "endo" = "Endothelial", "mg" = "Microglia", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "nonda" = "Neurons", "da" = "DaNeurons")) +
scale_y_discrete(expand = c(0,0))
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),height = 13, width = 3)


########################################
########################################
########################################
########################################
######################################## FC

## Load in expression for all cell types
celltypelist <- c("astro", "endo", "mg", "olig", "opc", "peri", "nonda")
X <- "fill" 
avg_log2FC <- "fill"  
p_val_adj <- "fill" 
celltype <- "fill" 
fill <- data.frame(X, avg_log2FC, p_val_adj, celltype)

for(i in celltypelist){
temp1 <- paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/',i,'_intergrated_MAST_DGE.csv')
temp1 <- read.delim(temp1, header = T, sep = ",")
temp1 <- temp1 %>% dplyr::select(X, avg_log2FC, p_val_adj)
temp1$celltype <- i
fill <- rbind(fill, temp1)
}

## add da
temp1 <- paste0('/home/fiorini9/scratch/machine_learning/DEG/MAST_kamath/Cell_based_celltype_groups/da/da_DEG.csv')
temp1 <- read.delim(temp1, header = T, sep = ",")
temp1 <- temp1 %>% dplyr::select(X, avg_log2FC, p_val_adj)
temp1$celltype <- "da"
fill <- rbind(fill, temp1)

fill <- subset(fill, celltype != "fill")
fill$avg_log2FC <- as.numeric(fill$avg_log2FC)*-1
colnames(fill) <- c("feature", "avg_log2FC", "p_val_adj", "celltype")

complete_df <- merge(fill_lim2, fill, by = c("feature","celltype"), all.x = T, all.y = F)
nrow(fill_lim2) == nrow(complete_df)

unique(complete_df$celltype)

complete_df <- subset(complete_df, feature %in% filtered_genes)
complete_df$feature <- factor(complete_df$feature, levels = unique(fill_lim2$feature))
complete_df$celltype <- factor(complete_df$celltype, levels = c("astro", "endo", "mg", "olig",  "opc", "peri", "nonda", "da"))

blups <- rev(brewer.pal(11, "RdBu"))
head(fill)

min <- min(complete_df$avg_log2FC)
max <- max(complete_df$avg_log2FC)

## add missing entries
# Create a dataframe with all combinations of feature and celltype
all_combinations <- expand.grid(
feature = unique(complete_df$feature),
celltype = unique(complete_df$celltype),
stringsAsFactors = FALSE
)

# Merge with the existing dataframe, filling in missing values with NA
result_df <- all_combinations %>%
left_join(complete_df, by = c("feature", "celltype"))

result_df$feature <- factor(result_df$feature, levels = unique(fill_lim2$feature))
result_df$celltype <- factor(result_df$celltype, levels = c("astro", "endo", "mg", "olig",  "opc", "peri", "nonda", "da"))

Expression <- ggplot(result_df, aes(y = celltype, x = feature, fill = as.numeric(avg_log2FC))) +
geom_tile(colour = "white", size = 1) +
theme_bw() +
theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = c("black"), size = 14)) +
labs(fill = "")  +
scale_fill_gradientn(colors = blups, values = scales::rescale(c(min,  0, max)), na.value = "lightgrey")+
scale_y_discrete(expand = c(0,0), labels=c("astro" = "Astrocytes", "all" = "All cells",  "endo" = "Endothelial", "mg" = "Microglia", "olig" = "Oligodendrocytes", "opc" = "OPC", "peri" = "Pericyte", "nonda" = "Neurons", "da" = "DaNeurons")) +
scale_x_discrete(expand = c(0,0)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),width = 13, height = 3)

########################
## remake the LIME plot with Greyed NA values
########################
blups <- brewer.pal(9, "Oranges")

LIME <- ggplot(result_df, aes(x = feature, y = celltype, fill = as.numeric(mean_z_score_weighted))) +
geom_tile(colour = "white", size = 1) +
theme_bw() +
theme(panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", colour = c("black"), size = 14),
    axis.text.y = element_text(colour = c("black"), size = 14)) +
labs(fill = "") +
scale_fill_gradientn(colors = blups, na.value = "lightgrey") +
scale_y_discrete(expand = c(0,0), labels = c("astro" = "Astrocytes", "all" = "All cells", "endo" = "Endothelial", "mg" = "Microglia", "olig" = "Oligodendrocytes", "opc" = "OPC", "peri" = "Pericyte", "nonda" = "Neurons", "da" = "DaNeurons")) +
scale_x_discrete(expand = c(0,0)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),height = 3, width = 15)

    
######################################## PD GWAS
########################################
########################################
########################################
######################################## 
########################################
## Identify PD genes
########################################
## label PD GWAS genes
## Nalls
all_PD_gwas_loci_genes <- read.delim(paste('/home/fiorini9/scratch/PD_machine_lerning/Nalls_GWAS.csv', sep=""), header = T, sep = ",")
all_PD_gwas_loci_genes <- all_PD_gwas_loci_genes[!is.na(all_PD_gwas_loci_genes$P..all.studies),]
nrow(all_PD_gwas_loci_genes)
unique(all_PD_gwas_loci_genes$Nearest.Gene)
Nalls_genes <- unique(all_PD_gwas_loci_genes$Nearest.Gene)
Nalls <- all_PD_gwas_loci_genes %>% dplyr::select(Nearest.Gene, P..all.studies)
colnames(Nalls) <- c("feature", "P")
Nalls$P <- as.numeric(Nalls$P)
Nalls$study <- "Nalls"
Nalls <- Nalls %>%
group_by(feature) %>%
filter(P == min(P)) %>%
ungroup()
Nalls <- data.frame(Nalls)
Nalls <- subset(Nalls, P <= 5e-8)

## Kim
all_PD_gwas_loci_genes <- read.delim(paste('/home/fiorini9/scratch/machine_learning/gene_lists/KIM_GWAS_2023.csv', sep=""), header = T, sep = ",")
subset(all_PD_gwas_loci_genes, Nearest.Gene.Feature == "RIMS1")
Kim_genes <- unique(all_PD_gwas_loci_genes$Nearest.Gene.Feature)
Kim <- all_PD_gwas_loci_genes %>% dplyr::select(Nearest.Gene.Feature, P.RE., P.MR.MEGA.)
Kim <- Kim %>%
        mutate(Min_P = pmin(P.RE., P.MR.MEGA., na.rm = TRUE))
Kim <- Kim %>%
    group_by(Nearest.Gene.Feature) %>%
    filter(Min_P == min(Min_P)) %>%
    ungroup()
nrow(Kim)
unique(Kim$Nearest.Gene.Feature)
Kim <- data.frame(Kim)
table(Kim$Nearest.Gene.Feature)
Kim <- Kim[!duplicated(Kim$Nearest.Gene.Feature),]
nrow(Kim) == length(unique(Kim$Nearest.Gene.Feature))
Kim <- Kim %>% dplyr::select(Nearest.Gene.Feature, Min_P)
colnames(Kim) <- c("feature", "P")
Kim$study <- "Kim"
Kim <- subset(Kim, P <= 5e-8)
subset(Kim, feature == "RNU4-66P")
Kim$feature[Kim$feature == "RNU4-66P"] <- "RIMS1"
subset(Kim, feature == "RIMS1")

## Combine Nalls and Kim
GWAS_bind <- rbind(Nalls, Kim)
PD_genes <- GWAS_bind$feature
unique_genes <- unique(as.character(fill_lim2$feature))
temp_df <- data.frame(unique_genes)
temp_df$Nalls <- "No"
temp_df$Kim <- "No"
temp_df$Nalls[temp_df$unique_genes %in% Nalls$feature] <- "Yes"
temp_df$Kim[temp_df$unique_genes %in% Kim$feature] <- "Yes"

# Reshape from wide to long format
long_df <- temp_df %>%
pivot_longer(
    cols = c("Nalls", "Kim"),
    names_to = "metric",
    values_to = "value"
) %>%
arrange(unique_genes, metric)
long_df <- data.frame(long_df)

long_df$unique_genes <- factor(long_df$unique_genes, levels = unique(fill_lim2$feature))
long_df$metric <- factor(long_df$metric, levels = c("Nalls", "Kim"))

## Plot
GWAS <- ggplot(long_df, aes(y = metric, x = unique_genes, fill = value)) +
    geom_tile(colour = "white", size = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background =element_rect(fill=c("white"), colour = "white"),
    strip.text = element_text(colour = "black", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 14)) +
    labs(fill = "")  +
    scale_fill_manual(values = c("lightgrey", "#af4b91"))+
    scale_y_discrete(expand = c(0,0), labels=c("Nalls" = "Nalls et al. 2019", "Kim" = "Kim et al. 2023")) +
    scale_x_discrete(expand = c(0,0)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),width = 12, height = 1.5)


######################################## NDKP
########################################
########################################
########################################
######################################## 
########################################
## HuGE
########################################
GWAS_PD_HUGE <- read.delim('/home/fiorini9/scratch/machine_learning/gene_lists/HUGE_PD.csv', header = T, sep = ",")
GWAS_PD_HUGE <- GWAS_PD_HUGE %>% dplyr::select(gene, bf_common, bf_rare, huge)
GWAS_PD_HUGE <- subset(GWAS_PD_HUGE, gene %in% filtered_genes)

# Identify genes that are not in the DataFrame
HUGE_gene <- GWAS_PD_HUGE$gene

########################################
## common
########################################
GWAS_PD_common <- read.delim('/home/fiorini9/scratch/machine_learning/gene_lists/NDKP_meta_PD_common.csv', header = T, sep = ",")
GWAS_PD_common <- GWAS_PD_common %>% dplyr::select(gene, Parkinsons.pValue, Parkinsons.zStat)
GWAS_PD_common <- subset(GWAS_PD_common, gene %in% filtered_genes)

# Identify genes that are not in the DataFrame
common_gene <- GWAS_PD_common$gene

########################################
## rare
########################################
GWAS_PD_rare <- read.delim('/home/fiorini9/scratch/machine_learning/gene_lists/NDKP_meta_PD_rare.csv', header = T, sep = ",")
GWAS_PD_rare <- GWAS_PD_rare %>% dplyr::select(gene, Parkinsons.pValue, Parkinsons.zStat)
GWAS_PD_rare <- subset(GWAS_PD_rare, gene %in% filtered_genes)

# Identify genes that are not in the DataFrame
rare_gene <- GWAS_PD_rare$gene

########################################
## Create dataframe and plot
########################################
temp_df <- data.frame(unique_genes)

temp_df$common <- "No"
temp_df$rare <- "No"
temp_df$HuGE <- "No"

temp_df$common[temp_df$unique_genes %in% common_gene] <- "common"
temp_df$rare[temp_df$unique_genes %in% rare_gene] <- "rare"
temp_df$HuGE[temp_df$unique_genes %in% HUGE_gene] <- "HuGe"

## conver to long
# Reshape from wide to long format
long_df <- temp_df %>%
pivot_longer(
    cols = c("common", "rare", "HuGE"),
    names_to = "metric",
    values_to = "value"
) %>%
arrange(unique_genes, metric)
long_df <- data.frame(long_df)

long_df$unique_genes <- factor(long_df$unique_genes, levels = unique(fill_lim2$feature))
long_df$metric <- factor(long_df$metric, levels = c("common", "rare", "HuGE"))

## modify the dataframe to only have one entry for PD
long_df <- long_df %>%
    group_by(unique_genes) %>%
    mutate(PD = if_else(all(value == "No"), "No", "Yes")) %>%
    ungroup()
long_df <- data.frame(long_df)
long_df <- long_df %>% dplyr::select(unique_genes, PD)
long_df <- aggregate(. ~ unique_genes, data = long_df, FUN = function(x) x[1])

## read in file for the remaining diseases
NDKP_ND <- read.delim('/home/fiorini9/scratch/machine_learning/gene_lists/ND_NDKP_entries.txt', header = T, sep = "\t")
NDKP_ND <- NDKP_ND %>% dplyr::select(Gene, Phenotype)
NDKP_ND$Gene[NDKP_ND$Gene == "LPHN2 (ADGRL2)"] <- "LPHN2"

unique(NDKP_ND$Phenotype)
NDKP_ND$Phenotype[NDKP_ND$Phenotype == "Alzheimer's disease"] <- "AD"
NDKP_ND$Phenotype[NDKP_ND$Phenotype == "Amyotrophic lateral sclerosis"] <- "ALS"
NDKP_ND$Phenotype[NDKP_ND$Phenotype == "Multiple sclerosis"] <- "MS"
NDKP_ND$Phenotype[NDKP_ND$Phenotype == "Lewy body dementia"]<- "LBD"
NDKP_ND$Phenotype[NDKP_ND$Phenotype == "Lewy Body Dementia"]<- "LBD"

# Transform the dataframe
df_transformed <- NDKP_ND %>%
    # Create a column for each unique phenotype
    mutate(value = "Yes") %>%
    pivot_wider(names_from = Phenotype, values_from = value, values_fill = list(value = "No"))
df_transformed <- data.frame(df_transformed)

long_df_merge <- merge(long_df, df_transformed, by.x = "unique_genes", by.y = "Gene", all.x = T)

long_df_merge <- long_df_merge %>%
    mutate(across(everything(), ~replace_na(.x, "No")))

long_df_merge <- long_df_merge %>%
    pivot_longer(
        cols = c("PD", "AD", "ALS", "MS", "LBD"),
        names_to = "metric",
        values_to = "value"
    ) %>%
    arrange(unique_genes, metric)

long_df_merge <- data.frame(long_df_merge)

long_df_merge$metric <- factor(long_df_merge$metric, levels = c("PD", "AD", "ALS", "LBD", "MS"))

long_df_merge$label <- NA
long_df_merge$label[long_df_merge$metric == "PD" & long_df_merge$unique_genes == "ARL17B"] <- "*"
long_df_merge$label[long_df_merge$metric == "PD" & long_df_merge$unique_genes == "CHORDC1"] <- "*"
long_df_merge$label[long_df_merge$metric == "PD" & long_df_merge$unique_genes == "BAG3"] <- "*"
long_df_merge$label[long_df_merge$metric == "PD" & long_df_merge$unique_genes == "SIPA1L2"] <- "*"
long_df_merge$label[long_df_merge$metric == "PD" & long_df_merge$unique_genes == "TMEM163"] <- "*"
long_df_merge$label[long_df_merge$metric == "PD" & long_df_merge$unique_genes == "GPC6"] <- "*"

long_df_merge$label[long_df_merge$metric == "AD" & long_df_merge$unique_genes == "ARL17B"] <- "*"
long_df_merge$label[long_df_merge$metric == "AD" & long_df_merge$unique_genes == "STIP1"] <- "*"

long_df_merge$label[long_df_merge$metric == "ALS" & long_df_merge$unique_genes == "MOBP"] <- "*"

long_df_merge$label[long_df_merge$metric == "MS" & long_df_merge$unique_genes == "HSPA1A"] <- "*"
long_df_merge$label[long_df_merge$metric == "MS" & long_df_merge$unique_genes == "JUND"] <- "*"
long_df_merge$label[long_df_merge$metric == "MS" & long_df_merge$unique_genes == "HSPA1B"] <- "*"
long_df_merge$label[long_df_merge$metric == "MS" & long_df_merge$unique_genes == "HLA-B"] <- "*"

NDKP <- ggplot(long_df_merge, aes(y = metric, x = unique_genes, fill = value, label = label)) +
    geom_tile(colour = "white", size = 1) +
    geom_text(size = 7, hjust = 0.5, vjust =0.75) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background =element_rect(fill=c("white"), colour = "white"),
    strip.text = element_text(colour = "black", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text( colour = "black", size = 14)) +
    labs(fill = "")  +
    scale_fill_manual(values = c("lightgrey", "#e6a532"))+
    scale_y_discrete(expand = c(0,0), labels=c("PD" = "Parkinson's disease", "AD" = "Alzheimer's disease", "ALS" = "Amyotrophic lateral sclerosis", "LBD" = "Lewy body dementia", "MS" = "Multiple sclerosis")) +
    scale_x_discrete(expand = c(0,0)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),width = 15, height = 1.5)

######################################## eQTL
########################################
########################################
########################################
######################################## 
###################################
## Use previous dataframe
###################################
result_df2 <- result_df

###################################
## How many all specific eQTL
###################################
## ADD Byrois eQTL
all_eQTL <- read.delim('/home/fiorini9/scratch/machine_learning/gene_lists/All_cis_eQTL_Bryois.csv', header = T, sep = ",")
PD_eQTL <- read.delim('/home/fiorini9/scratch/machine_learning/gene_lists/PD_cis_eQTL_Bryois.csv', header = T, sep = ",")

unique(all_eQTL$cell_type)
unique(PD_eQTL$tissue)
unique(result_df2$celltype)

all_eQTL <- subset(all_eQTL, adj_p <= 0.05)
all_eQTL_lim <- subset(all_eQTL, Replication == "5% FDR")

unique(all_eQTL_lim$cell_type)
result_df2$all_eQTL <- "No"
result_df2$all_eQTL[result_df2$celltype == "astro" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Astrocytes"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "endo" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Endothelial cells"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "nonda" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Excitatory neurons"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "nonda" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Inhibitory neurons"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "mg" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Microglia"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "olig" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Oligodendrocytes"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "opc" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "OPCs / COPs"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "peri" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Pericytes"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "da" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Excitatory neurons"] ] <- "Yes"
result_df2$all_eQTL[result_df2$celltype == "da" & result_df2$Freq >= 1 & result_df2$feature %in% all_eQTL_lim$symbol[all_eQTL_lim$cell_type == "Inhibitory neurons"] ] <- "Yes"

result_df2$label <- NA
result_df2$label[result_df2$feature == "ARL17B" & result_df2$celltype == "astro"] <- "X"
result_df2$label[result_df2$feature == "ARL17B" & result_df2$celltype == "mg"] <- "X"
result_df2$label[result_df2$feature == "ARL17B" & result_df2$celltype == "olig"] <- "X"
result_df2$label[result_df2$feature == "TMEM163" & result_df2$celltype == "mg"] <- "X"

eQTL <- ggplot(result_df2, aes(y = celltype, x = feature, fill = all_eQTL, label = label)) +
    geom_tile(colour = "white", size = 1) +
    geom_text()+
    theme_bw() +
    theme(panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background =element_rect(fill=c("white"), colour = "white"),
    strip.text = element_text(colour = "black", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 14)) +
    labs(fill = "")  +
    scale_fill_manual(values = c("lightgrey","#41afaa"))+
    scale_y_discrete(expand = c(0,0), labels=c("astro" = "Astrocytes", "all" = "All cells",  "endo" = "Endothelial", "mg" = "Microglia", "olig" = "Oligodendrocytes", "opc" = "OPC", "peri" = "Pericyte", "nonda" = "Neurons", "da" = "DaNeurons")) +
    scale_x_discrete(expand = c(0,0)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),width = 15, height = 1.5)

######################################## arrange
########################################
########################################
########################################
######################################## 
ggarrange(NDKP, GWAS, eQTL, Expression, LIME,   ncol = 1, nrow = 5, align = "v", heights = c(0.32, 0.15, 0.5, 0.5,  0.7))
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temps.pdf',sep=""),width = 19, height = 10)
