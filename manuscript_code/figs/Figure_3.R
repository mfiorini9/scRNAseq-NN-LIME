## This script was used to produce Figure 3.
salloc -A def-sfarhan --time=0-8 -c 6 --mem=40g

module load StdEnv/2020 
module load r/4.2.2 
R

library(Seurat)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(rstatix)
library(inflection, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(SeuratData, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(circlize)
library(ggrepel)
library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(ggridges)

############################################################################ LIME feature importance correlation, Wang and Kamath
############################################################################
############################################################################
############################################################################
## Cell type List
cell_type_list <- c('astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'all')

## Loop through each cell type to produce scatter plot
for(i in unique(cell_type_list)){
    ## first dataset
    temp_1 <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_kamath_',i,'.csv')
    temp_1 <- read.delim(temp_1, header = T, sep = ",")
    temp_1$celltype <- i
    temp_1$dataset <- "Kamath"

    ## second dataset
    temp_2 <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_wang_',i,'.csv')
    temp_2 <- read.delim(temp_2, header = T, sep = ",")
    temp_2$celltype <- i
    temp_2$dataset <- "Wang"

    ## find intersecting genes
    intersect <- intersect(temp_1$feature, temp_2$feature)

    temp_1 <- subset(temp_1, feature %in% intersect)
    length(unique(temp_1$feature))
    temp_2 <- subset(temp_2, feature %in% intersect)
    length(unique(temp_2$feature))

    ## remove duplicated gene entries for each dataset
    temp_1 <- aggregate(. ~ feature, data = temp_1, FUN = function(x) x[1])
    nrow(temp_1)
    temp_2 <- aggregate(. ~ feature, data = temp_2, FUN = function(x) x[1])
    nrow(temp_2)

    ## merge
    temp_1 <- temp_1 %>% dplyr::select(feature, Mean_z_score_weighted)
    colnames(temp_1) <- c("feature", "Kamath_Mean_z_score_weighted")

    temp_2 <- temp_2 %>% dplyr::select(feature, Mean_z_score_weighted)
    colnames(temp_2) <- c("feature", "Wang_Mean_z_score_weighted")

    bind <- merge(temp_1, temp_2, by = "feature")
    nrow(bind)

    ## Top in both Kamath and Wang
    bind <- bind %>%
        mutate(percent_rank_kamath = percent_rank(as.numeric(Kamath_Mean_z_score_weighted)))

    bind <- bind %>%
        mutate(percent_rank_wang = percent_rank(as.numeric(Wang_Mean_z_score_weighted)))

    bind$sum <- as.numeric(bind$percent_rank_wang) + as.numeric(bind$percent_rank_kamath)

    bind <- bind[order(-bind$sum),]
    top <- unique(bind$feature)[1:10]

    bind$label <- NA
    bind$label[bind$feature %in% top] <- bind$feature[bind$feature %in% top]

    ## plot
    ggplot(bind, aes(x = as.numeric(Kamath_Mean_z_score_weighted), y = as.numeric(Wang_Mean_z_score_weighted), label = label)) +
    geom_point(colour = "darkgrey", size = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", size = 0.25, colour = "red") +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.25, colour = "red") +
    theme_bw() +
    theme(axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)) +
    geom_smooth(method='lm', fill = "grey", colour = "black", size = 0.5) +
    xlab("Kamath: Z-score") +
    ylab("Wang: Z-score") +
    ggtitle(i) +
    scale_y_continuous(limits = c(-3,20)) + 
    scale_x_continuous(limits = c(-3,20)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/kamath_wang_correlation_',i,'.pdf',sep=""),width = 1.75, height = 2.5)

}

############################################################################ Ridge plot
############################################################################
############################################################################
############################################################################
##################
## Kamath
##################
cell_type_list <- c('astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'da', 'all')

    X <- "fill"
    feature<- "fill"
    Mean_feature_importance<- "fill"
    Sample_ID<- "fill"
    z_scores<- "fill"
    Mean_z_score<- "fill"
    celltype <- "fill"
    dataset <- "fill"

    fill <- data.frame(X, feature, Mean_feature_importance, Sample_ID, z_scores, Mean_z_score, celltype, dataset)

    for(i in unique(cell_type_list)){
        temp <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_kamath_',i,'.csv')
        temp <- read.delim(temp, header = T, sep = ",")
        temp$celltype <- i
        temp$dataset <- "Kamath"

        fill <- rbind(fill, temp)
    }

    fill <- subset(fill, X != "fill")
    kamath_df <- fill
    
## Factorize cell types
kamath_df$celltype <- factor(kamath_df$celltype, levels = rev(c("all", 'astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'da')))

## Ridge plot
plot_kamath <- ggplot(kamath_df, aes(x=as.numeric(Mean_z_score), y = celltype, fill = celltype)) +
geom_density_ridges() +
theme_ridges() + 
theme(axis.y.title = element_blank(),
        axis.x.title = element_blank()) +
scale_fill_manual(values = c("astro" = "#CAB2D6",
        "endo"="#FDBF6F",
        "olig"="#B2DF8A",
        "mg"="#FB9A99",
        "da"="#A6CEE3",
        "nonda"="#1F78B4",
        "opc"="#33A02C",
        "peri"="#B15928",
        "all" = "black"))+
theme(legend.position = "none") +
scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "all" = "All cells")) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 2.75, height = 6)

##################
## Wang
##################
cell_type_list <- c('astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', "all")

X <- "fill"
feature<- "fill"
Mean_feature_importance<- "fill"
Sample_ID<- "fill"
z_scores<- "fill"
Mean_z_score<- "fill"
celltype <- "fill"
dataset <- "fill"

fill <- data.frame(X, feature, Mean_feature_importance, Sample_ID, z_scores, Mean_z_score, celltype, dataset)

for(i in unique(cell_type_list)){
    temp <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_wang_',i,'.csv')
    temp <- read.delim(temp, header = T, sep = ",")
    temp$celltype <- i
    temp$dataset <- "Wang"

    fill <- rbind(fill, temp)
}

fill <- subset(fill, X != "fill")
wang_df <- fill

## Factorize cell types     
wang_df$celltype <- factor(wang_df$celltype, levels = rev(c("all", 'astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'da')))

## Add DaN filler
X <- "fill"
feature<- "fill"
Mean_feature_importance<- "fill"
Sample_ID<- "fill"
z_scores<- "fill"
Mean_z_score<- NA
celltype <- "da"
dataset <- 'Wang'

fill <- data.frame(X,
feature,
Mean_feature_importance,
Sample_ID,
z_scores,
Mean_z_score,
celltype,
dataset)

wang_df <-rbind(wang_df, fill)

plot_wang <- ggplot(wang_df, aes(x=as.numeric(Mean_z_score), y = celltype, fill = celltype)) +
    geom_density_ridges() +
    theme_ridges() + 
    theme(axis.y.title = element_blank(),
        axis.x.title = element_blank()) +
    scale_fill_manual(values = c("astro" = "#CAB2D6",
            "endo"="#FDBF6F",
            "olig"="#B2DF8A",
            "mg"="#FB9A99",
            "da"="#A6CEE3",
            "nonda"="#1F78B4",
            "opc"="#33A02C",
            "peri"="#B15928",
            "all" = "black"))+
    theme(legend.position = "none") +
    scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "all" = "All cells")) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 2.75, height = 6)

##################
## Kamath and Wang
##################
cell_type_list <- c('astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'all')

X <- "fill"
feature<- "fill"
Mean_z_score_weighted <- "fill"
celltype <- "fill"

fill <- data.frame(X, feature, Mean_z_score_weighted, celltype)

for(i in unique(cell_type_list)){
    temp <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_master_LIME_whole_',i,'.csv')
    temp <- read.delim(temp, header = T, sep = ",")
    temp$celltype <- i

    fill <- rbind(fill, temp)
}

fill <- subset(fill, X != "fill")
combo_df <- fill

## ridge plot
combo_df$celltype <- factor(combo_df$celltype, levels = rev(c("all", 'astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'da')))

## Add DaN filler
X <- "fill"
feature<- "fill"
Mean_z_score_weighted<- "fill"
celltype <- "da"

fill <- data.frame(X,
feature,
Mean_z_score_weighted,
celltype)

combo_df <-rbind(combo_df, fill)

## Combo ridgeplot
plot_combo <- ggplot(combo_df, aes(x=as.numeric(Mean_z_score_weighted), y = celltype, fill = celltype)) +
    geom_density_ridges() +
    theme_ridges() + 
    theme(axis.y.title = element_blank(),
        axis.x.title = element_blank()) +
    scale_fill_manual(values = c("astro" = "#CAB2D6",
            "endo"="#FDBF6F",
            "olig"="#B2DF8A",
            "mg"="#FB9A99",
            "da"="#A6CEE3",
            "nonda"="#1F78B4",
            "opc"="#33A02C",
            "peri"="#B15928",
            "all" = "black"))+
    theme(legend.position = "none") +
    scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "all" = "All cells")) +
    scale_x_continuous(limits = c(-1,15)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 2.75, height = 6)

##################
## Arrange plots
##################
kamath_df <- kamath_df %>% dplyr::select(X, feature, Mean_z_score, celltype)
kamath_df$dataset <- "kamath"
head(subset(kamath_df, celltype == "da"))

wang_df <- wang_df %>% dplyr::select(X, feature, Mean_z_score, celltype)
wang_df$dataset <- "wang"

combo_df$dataset <- "combo"
colnames(combo_df) <- colnames(wang_df)
head(subset(combo_df, celltype == "da"))

total_df <- rbind(kamath_df, wang_df, combo_df)
total_df$dataset <- factor(total_df$dataset, levels = c("kamath", "wang", "combo"))

plot_combo <- ggplot(total_df, aes(x=as.numeric(Mean_z_score), y = celltype, fill = celltype)) +
            geom_density_ridges() +
            theme_ridges() + 
            facet_wrap(~dataset) +
            theme(axis.y.title = element_blank(),
                axis.x.title = element_blank(),
                    strip.background =element_rect(fill=c("white"), colour = "white")) +
            scale_fill_manual(values = c("astro" = "#CAB2D6",
                    "endo"="#FDBF6F",
                    "olig"="#B2DF8A",
                    "mg"="#FB9A99",
                    "da"="#A6CEE3",
                    "nonda"="#1F78B4",
                    "opc"="#33A02C",
                    "peri"="#B15928",
                    "all" = "black"))+
            theme(legend.position = "none") +
            scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "all" = "All cells")) +
            scale_x_continuous(limits = c(-1,15)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 5.5, height = 6)       


############################################################################ Plot permutations and optimal gene set bar plot
############################################################################
############################################################################
############################################################################

optimal_kamath_wang_unweighted <- function(){  
    ############################################################################ all
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "all"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_all_weighted <- max(bind$combo_diff_weighted)
    max_2_all_weighted <- subset(bind, combo_diff_weighted == max_2_all_weighted)
    max_2_all_weighted <- max_2_all_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_all_weighted)
    n_genes_all <- length(unique(rank_list_lim$feature))
    z_thresh_all <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_all_wang = total_wang$Mean[total_wang$method == max_2_all_weighted & total_wang$ident == 'random_gene']
    y4_all_wang = total_wang$Mean[total_wang$method == max_2_all_weighted & total_wang$ident == 'optimal']


    all_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_all_weighted, y = y3_all_wang, xend = max_2_all_weighted, yend = y4_all_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)), labels = scales::number_format(accuracy = 0.01)) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_all_kamath = total_kamath$Mean[total_kamath$method == max_2_all_weighted & total_kamath$ident == 'random_gene']
    y4_all_kamath = total_kamath$Mean[total_kamath$method == max_2_all_weighted & total_kamath$ident == 'optimal']


    all_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_all_weighted, y = y3_all_kamath, xend = max_2_all_weighted, yend = y4_all_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)), labels = scales::number_format(accuracy = 0.01)) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_all_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_all_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_Smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_Smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_Smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_all_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_all <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_all$genes <- factor(bind_all$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_all <- data_summary(bind_all, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_all$genes <- factor(bind_all$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    all_box <- ggplot(bind_all, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
            position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_text(colour = "black"),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_all)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    
    ############################################################################ astro
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "astro"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_astro_weighted <- max(bind$combo_diff_weighted)
    max_2_astro_weighted <- subset(bind, combo_diff_weighted == max_2_astro_weighted)
    max_2_astro_weighted <- max_2_astro_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_astro_weighted)
    n_genes_astro <- length(unique(rank_list_lim$feature))
    z_thresh_astro <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_astro_wang = total_wang$Mean[total_wang$method == max_2_astro_weighted & total_wang$ident == 'random_gene']
    y4_astro_wang = total_wang$Mean[total_wang$method == max_2_astro_weighted & total_wang$ident == 'optimal']


    astro_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_astro_weighted, y = y3_astro_wang, xend = max_2_astro_weighted, yend = y4_astro_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)), labels = scales::number_format(accuracy = 0.01)) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_astro_kamath = total_kamath$Mean[total_kamath$method == max_2_astro_weighted & total_kamath$ident == 'random_gene']
    y4_astro_kamath = total_kamath$Mean[total_kamath$method == max_2_astro_weighted & total_kamath$ident == 'optimal']


    astro_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_astro_weighted, y = y3_astro_kamath, xend = max_2_astro_weighted, yend = y4_astro_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)), labels = scales::number_format(accuracy = 0.01)) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_astro_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_astro_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_astro_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_astro <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_astro$genes <- factor(bind_astro$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_astro <- data_summary(bind_astro, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_astro$genes <- factor(bind_astro$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    astro_box <- ggplot(bind_astro, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_astro)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ############################################################################ endo
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "endo"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_endo_weighted <- max(bind$combo_diff_weighted)
    max_2_endo_weighted <- subset(bind, combo_diff_weighted == max_2_endo_weighted)
    max_2_endo_weighted <- max_2_endo_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_endo_weighted)
    n_genes_endo <- length(unique(rank_list_lim$feature))
    z_thresh_endo <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_endo_wang = total_wang$Mean[total_wang$method == max_2_endo_weighted & total_wang$ident == 'random_gene']
    y4_endo_wang = total_wang$Mean[total_wang$method == max_2_endo_weighted & total_wang$ident == 'optimal']


    endo_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_endo_weighted, y = y3_endo_wang, xend = max_2_endo_weighted, yend = y4_endo_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)), labels = scales::number_format(accuracy = 0.01)  ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_endo_kamath = total_kamath$Mean[total_kamath$method == max_2_endo_weighted & total_kamath$ident == 'random_gene']
    y4_endo_kamath = total_kamath$Mean[total_kamath$method == max_2_endo_weighted & total_kamath$ident == 'optimal']


    endo_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_endo_weighted, y = y3_endo_kamath, xend = max_2_endo_weighted, yend = y4_endo_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)), labels = scales::number_format(accuracy = 0.01)  ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_endo_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_endo_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_endo_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_endo <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_endo$genes <- factor(bind_endo$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_endo <- data_summary(bind_endo, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_endo$genes <- factor(bind_endo$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    endo_box <- ggplot(bind_endo, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_endo)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ############################################################################ mg
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "mg"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_mg_weighted <- max(bind$combo_diff_weighted)
    max_2_mg_weighted <- subset(bind, combo_diff_weighted == max_2_mg_weighted)
    max_2_mg_weighted <- max_2_mg_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_mg_weighted)
    n_genes_mg <- length(unique(rank_list_lim$feature))
    z_thresh_mg <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_mg_wang = total_wang$Mean[total_wang$method == max_2_mg_weighted & total_wang$ident == 'random_gene']
    y4_mg_wang = total_wang$Mean[total_wang$method == max_2_mg_weighted & total_wang$ident == 'optimal']


    mg_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_mg_weighted, y = y3_mg_wang, xend = max_2_mg_weighted, yend = y4_mg_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_mg_kamath = total_kamath$Mean[total_kamath$method == max_2_mg_weighted & total_kamath$ident == 'random_gene']
    y4_mg_kamath = total_kamath$Mean[total_kamath$method == max_2_mg_weighted & total_kamath$ident == 'optimal']


    mg_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_mg_weighted, y = y3_mg_kamath, xend = max_2_mg_weighted, yend = y4_mg_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_mg_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_mg_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_mg_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_mg <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_mg$genes <- factor(bind_mg$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_mg <- data_summary(bind_mg, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_mg$genes <- factor(bind_mg$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    mg_box <- ggplot(bind_mg, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_mg)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ############################################################################ olig
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "olig"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_olig_weighted <- max(bind$combo_diff_weighted)
    max_2_olig_weighted <- subset(bind, combo_diff_weighted == max_2_olig_weighted)
    max_2_olig_weighted <- max_2_olig_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_olig_weighted)
    n_genes_olig <- length(unique(rank_list_lim$feature))
    z_thresh_olig <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_olig_wang = total_wang$Mean[total_wang$method == max_2_olig_weighted & total_wang$ident == 'random_gene']
    y4_olig_wang = total_wang$Mean[total_wang$method == max_2_olig_weighted & total_wang$ident == 'optimal']


    olig_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_olig_weighted, y = y3_olig_wang, xend = max_2_olig_weighted, yend = y4_olig_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_olig_kamath = total_kamath$Mean[total_kamath$method == max_2_olig_weighted & total_kamath$ident == 'random_gene']
    y4_olig_kamath = total_kamath$Mean[total_kamath$method == max_2_olig_weighted & total_kamath$ident == 'optimal']


    olig_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_olig_weighted, y = y3_olig_kamath, xend = max_2_olig_weighted, yend = y4_olig_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_olig_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_olig_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_olig_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_olig <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_olig$genes <- factor(bind_olig$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_olig <- data_summary(bind_olig, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_olig$genes <- factor(bind_olig$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    olig_box <- ggplot(bind_olig, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_olig)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ############################################################################ opc
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "opc"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_opc_weighted <- max(bind$combo_diff_weighted)
    max_2_opc_weighted <- subset(bind, combo_diff_weighted == max_2_opc_weighted)
    max_2_opc_weighted <- max_2_opc_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_opc_weighted)
    n_genes_opc <- length(unique(rank_list_lim$feature))
    z_thresh_opc <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_opc_wang = total_wang$Mean[total_wang$method == max_2_opc_weighted & total_wang$ident == 'random_gene']
    y4_opc_wang = total_wang$Mean[total_wang$method == max_2_opc_weighted & total_wang$ident == 'optimal']


    opc_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_opc_weighted, y = y3_opc_wang, xend = max_2_opc_weighted, yend = y4_opc_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_opc_kamath = total_kamath$Mean[total_kamath$method == max_2_opc_weighted & total_kamath$ident == 'random_gene']
    y4_opc_kamath = total_kamath$Mean[total_kamath$method == max_2_opc_weighted & total_kamath$ident == 'optimal']


    opc_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_opc_weighted, y = y3_opc_kamath, xend = max_2_opc_weighted, yend = y4_opc_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_opc_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_opc_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_opc_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_opc <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_opc$genes <- factor(bind_opc$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_opc <- data_summary(bind_opc, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_opc$genes <- factor(bind_opc$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    opc_box <- ggplot(bind_opc, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_opc)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ############################################################################ peri
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "peri"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_peri_weighted <- max(bind$combo_diff_weighted)
    max_2_peri_weighted <- subset(bind, combo_diff_weighted == max_2_peri_weighted)
    max_2_peri_weighted <- max_2_peri_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_peri_weighted)
    n_genes_peri <- length(unique(rank_list_lim$feature))
    z_thresh_peri <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_peri_wang = total_wang$Mean[total_wang$method == max_2_peri_weighted & total_wang$ident == 'random_gene']
    y4_peri_wang = total_wang$Mean[total_wang$method == max_2_peri_weighted & total_wang$ident == 'optimal']


    peri_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_peri_weighted, y = y3_peri_wang, xend = max_2_peri_weighted, yend = y4_peri_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_peri_kamath = total_kamath$Mean[total_kamath$method == max_2_peri_weighted & total_kamath$ident == 'random_gene']
    y4_peri_kamath = total_kamath$Mean[total_kamath$method == max_2_peri_weighted & total_kamath$ident == 'optimal']


    peri_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_peri_weighted, y = y3_peri_kamath, xend = max_2_peri_weighted, yend = y4_peri_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_peri_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_peri_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_peri_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_peri <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_peri$genes <- factor(bind_peri$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_peri <- data_summary(bind_peri, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_peri$genes <- factor(bind_peri$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    peri_box <- ggplot(bind_peri, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_peri)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ############################################################################ nonda
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "nonda"

    ########################################
    ########################################
    ########################################
    ######################################## Wang
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"


    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_wang <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_wang <- toal_c
    toal_c_wang <- toal_c_wang %>% dplyr::select(method, diff_weighted)

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both
    bind <- merge(toal_c_wang, toal_c_kamath, by = "method")
    bind$combo_diff_weighted <- bind$diff_weighted.x + bind$diff_weighted.y

    max_2_nonda_weighted <- max(bind$combo_diff_weighted)
    max_2_nonda_weighted <- subset(bind, combo_diff_weighted == max_2_nonda_weighted)
    max_2_nonda_weighted <- max_2_nonda_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_master_LIME_Z_1_',i,'.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_nonda_weighted)
    n_genes_nonda <- length(unique(rank_list_lim$feature))
    z_thresh_nonda <- min(rank_list_lim$Mean_z_score_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## dotplot wang

    y3_nonda_wang = total_wang$Mean[total_wang$method == max_2_nonda_weighted & total_wang$ident == 'random_gene']
    y4_nonda_wang = total_wang$Mean[total_wang$method == max_2_nonda_weighted & total_wang$ident == 'optimal']


    nonda_wang <- ggplot(total_wang, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_nonda_weighted, y = y3_nonda_wang, xend = max_2_nonda_weighted, yend = y4_nonda_wang), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_wang$method), mean(total_wang$method), max(total_wang$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_nonda_kamath = total_kamath$Mean[total_kamath$method == max_2_nonda_weighted & total_kamath$ident == 'random_gene']
    y4_nonda_kamath = total_kamath$Mean[total_kamath$method == max_2_nonda_weighted & total_kamath$ident == 'optimal']


    nonda_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_nonda_weighted, y = y3_nonda_kamath, xend = max_2_nonda_weighted, yend = y4_nonda_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## wang
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_wang.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_nonda_weighted)
    smaj$dataset <- "Dataset 2"
    box_df_wang <- smaj

    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_nonda_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj

    ########################################
    ## Smajic
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_optimal_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_HVG_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/Master_LIME_summary_eval_',i,'_random_gene_LIME_perm_lim_10000_smajic.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_nonda_weighted)
    smaj$dataset <- "Dataset 3"
    box_df_smajic <- smaj

    ########################################
    ## bind and plot
    ########################################

    bind_nonda <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_nonda$genes <- factor(bind_nonda$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_nonda <- data_summary(bind_nonda, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_nonda$genes <- factor(bind_nonda$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    nonda_box <- ggplot(bind_nonda, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_nonda)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    bind_nonda2 <- bind_nonda

    ############################################################################ da
    ############################################################################ 
    ############################################################################
    ############################################################################
    i = "da"

    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/LIME_summary_eval_da_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    random <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/LIME_summary_eval_da_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random$test <- "random"

    random_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/LIME_summary_eval_da_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    random_gene$test <- "random_gene"

    ## Optimal
    optimal <- optimal %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    optimal$ident <- "optimal"

    optimal <- optimal %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    optimal <- data.frame(optimal)

    ## Random HVG
    random <- random %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random$ident <- "random"

    random <- random %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random <- data.frame(random)

    ## Random gene
    random_gene <- random_gene %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(Mean = mean(X0, na.rm=TRUE))
    random_gene$ident <- "random_gene"

    random_gene <- random_gene %>%
    mutate(percent_rank = rank(Mean)/length(Mean))

    random_gene <- data.frame(random_gene)

    total <- rbind(optimal, random, random_gene)
    total_kamath <- total

    ## define optimal threshold
    toal_c <- merge(optimal, random_gene, by = "method")
    toal_c$diff <- toal_c$Mean.x - toal_c$Mean.y
    
    #toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    toal_c$diff_weighted <- toal_c$diff*toal_c$percent_rank.x
    #toal_c$diff_weighted <- toal_c$diff
    toal_c_kamath <- toal_c
    toal_c_kamath <- toal_c_kamath %>% dplyr::select(method, diff_weighted)
    
    ########################################
    ########################################
    ########################################
    ######################################## Find threshold that maximizes both

    max_2_da_weighted <- max(toal_c_kamath$diff_weighted)
    max_2_da_weighted <- subset(toal_c_kamath, diff_weighted == max_2_da_weighted)
    max_2_da_weighted <- max_2_da_weighted$method

    ## how many genes at threshold
    rank_list <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_kamath_percent_ranked_da.csv')
    rank_list <- read.delim(rank_list, header = T, sep = "," )

    rank_list_lim <- subset(rank_list, percent_rank >= max_2_da_weighted)
    n_genes_da <- length(unique(rank_list_lim$feature))
    z_thresh_da <- min(rank_list_lim$Mean_z_score_weighted)
    

    ########################################
    ########################################
    ########################################
    ######################################## dotplot kamath
    
    y3_da_kamath = total_kamath$Mean[total_kamath$method == max_2_da_weighted & total_kamath$ident == 'random_gene']
    y4_da_kamath = total_kamath$Mean[total_kamath$method == max_2_da_weighted & total_kamath$ident == 'optimal']


    da_kamath <- ggplot(total_kamath, aes(x = method, y = Mean, colour = ident, group = ident  )) + 
    geom_point(size = 1) +
    geom_line() +
    geom_segment(aes(x = max_2_da_weighted, y = y3_da_kamath, xend = max_2_da_weighted, yend = y4_da_kamath), colour = "black", linetype = "dashed") +
    theme_classic() + 
    theme(legend.position = "none",
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
    scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    scale_x_continuous(breaks = c(min(total_kamath$method), mean(total_kamath$method), max(total_kamath$method)) , labels = scales::number_format(accuracy = 0.01) ) +
    ylim(0.5, 1) +
    ggtitle(i) +
    ylab("Balanced accuracy") +
    xlab("Z-score percentile rank")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)


    ########################################
    ########################################
    ########################################
    ######################################## test combo box plot
    
    ########################################
    ## Kamath
    ########################################
    smaj_optimal <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/LIME_summary_eval_da_optimal_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_hvg <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/LIME_summary_eval_da_random_HVG_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    smaj_gene <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/LIME_summary_eval_da_random_gene_LIME_perm_lim_10000_kamath.csv'), header = T, sep = ",")
    
    ## Process 
    smaj_optimal$genes <- "optimal"
    smaj_hvg$genes <- "random_hvg"
    smaj_gene$genes <- "random_gene"
    
    colnames(smaj_optimal)
    colnames(smaj_hvg)
    colnames(smaj_gene)

    ## bind 
    smaj <- rbind(smaj_optimal, smaj_hvg, smaj_gene)

    ## subset for Z value
    smaj <- subset(smaj, method == max_2_da_weighted)
    smaj$dataset <- "Dataset 1"
    box_df_kamath <- smaj
    
    
    ########################################
    ## wang
    ########################################
    
    box_df_wang <- box_df_kamath
    box_df_wang$dataset <- "Dataset 2"
    box_df_wang$X0 <- 0


    ########################################
    ## Smajic
    ########################################
    box_df_smajic <- box_df_kamath
    box_df_smajic$dataset <- "Dataset 3"
    box_df_smajic$X0 <- 0
    
    ########################################
    ## bind and plot
    ########################################

    bind_nonda <- rbind(box_df_kamath, box_df_wang, box_df_smajic)
    bind_nonda$genes <- factor(bind_nonda$genes, levels = c("optimal", "random_hvg", "random_gene"))
    
    #+++++++++++++++++++++++++
    # Function to calculate the mean and the standard deviation
    # for each group
    #+++++++++++++++++++++++++
    # data : a data frame
    # varname : the name of a column containing the variable
    #to be summariezed
    # groupnames : vector of column names to be used as
    # grouping variables
    data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
    }

    bind_nonda <- data_summary(bind_nonda, varname="X0", 
                groupnames=c("genes", "dataset"))
    # Convert dose to a factor variable
    bind_nonda$genes <- factor(bind_nonda$genes, levels = c("optimal", "random_hvg", "random_gene"))
    

    ## plot
    da_box <- ggplot(bind_nonda, aes(x = dataset , y = X0, fill = genes)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=X0-sd, ymax=X0+sd), width=.2,
                position=position_dodge(.9)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = c('#1b9e77', '#e6ab02', '#7570b3' )),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none")+
    coord_cartesian(ylim=c(0.5,1.03)) +
    ylab("Balanced accuracy") + 
    scale_fill_manual(values = c("#F18F01","#048BA8" , "#2E4057")) +
    ggtitle(paste0(i, '\n',n_genes_da)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 3, height = 3)

    bind_da <- bind_nonda

    #####################
    ## arrange dot for kamath
    #####################
    ggarrange(all_kamath, astro_kamath, endo_kamath, mg_kamath, olig_kamath, opc_kamath, peri_kamath, nonda_kamath, da_kamath,  ncol = 9, nrow = 1, align = "h", widths = c(1.35,1,1,1,1,1,1,1,1))
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/Master_dot_kamath.pdf',sep=""),width = 10.5, height = 1.5)

    #####################
    ## arrange dot for wang
    #####################
    ggarrange(all_wang, astro_wang, endo_wang, mg_wang, olig_wang, opc_wang, peri_wang, nonda_wang, da_kamath, ncol = 9, nrow = 1, align = "h", widths = c(1.35,1,1,1,1,1,1,1, 1))
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/Master_dot_wang.pdf',sep=""),width = 10.5, height = 1.5)

    #####################
    ## arrange box top
    #####################
    ggarrange(astro_box, endo_box, mg_box, olig_box, ncol = 4, nrow = 1, align = "h", widths = c(1.325,1, 1,1,1))
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/Master_box_top.pdf',sep=""),width = 4.5, height = 2.5)

    #####################
    ## arrange box bottom
    #####################
    ggarrange(opc_box, peri_box, nonda_box, da_box, ncol = 4, nrow = 1, align = "h", widths = c(1.325,1, 1,1,1))
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/Master_box_bottom.pdf',sep=""),width = 4.5, height = 2.5)

    #####################
    ## arrange all box
    #####################
    ggarrange(all_box, astro_box, endo_box, mg_box, olig_box, opc_box, peri_box, nonda_box, da_box, ncol = 9, nrow = 1, align = "h", widths = c(1.4, 1,1, 1,1,1,1,1,1))
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/Master_box_all.pdf',sep=""),width = 7.5, height = 2.5)

    #####################
    ## Prin Z-score threshold csv file
    #####################
    z_score_thresh <- c(max_2_all_weighted, max_2_astro_weighted, max_2_endo_weighted, max_2_mg_weighted, max_2_olig_weighted, max_2_opc_weighted, max_2_peri_weighted, max_2_nonda_weighted, max_2_da_weighted)
    cell_type <- c("all", "astro", "endo", "mg", "olig", "opc", "peri", "nonda", "da")
    df <- data.frame(z_score_thresh, cell_type)
    write.csv(df, paste0('/home/fiorini9/scratch/machine_learning/temp_outs/optimal_Z_score_threshold_Master.csv'))  
}

    optimal_kamath_wang_unweighted()

    ## compute means
    bind_total <- rbind(bind_all, bind_astro, bind_mg, bind_endo, bind_olig, bind_opc, bind_peri, bind_nonda2, bind_da)

    bind_total <- subset(bind_total, as.numeric(X0)  != 0)

    bind_total <- bind_total %>%
        dplyr::group_by(genes, dataset) %>%
        dplyr::summarize(mean_X0 = mean(as.numeric(X0)))

    ## across all datasets
    bind_total <- rbind(bind_all, bind_astro, bind_mg, bind_endo, bind_olig, bind_opc, bind_peri, bind_nonda2, bind_da)

    bind_total <- subset(bind_total, as.numeric(X0)  != 0)

    bind_total <- bind_total %>%
        dplyr::group_by(genes) %>%
        dplyr::summarize(mean_X0 = mean(as.numeric(X0)))
