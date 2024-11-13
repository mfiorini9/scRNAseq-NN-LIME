## This script was used to produce Figure 4.

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

############################################################################ Create Pan-dataset, cell type specific objects
############################################################################
############################################################################
############################################################################

########################################
########################################
########################################
########################################
######################################## Load Seurat objects
## load Wang
seu_wang <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Wang.h5Seurat")    
ncol(seu_wang) 
DefaultAssay(seu_wang)
unique(seu_wang@meta.data$Disease_Status)

## load Kamath
seu_kamath <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Kamath.h5Seurat")    
ncol(seu_kamath) 
DefaultAssay(seu_kamath)
unique(seu_kamath@meta.data$Disease_Status) 

## load Smajic
seu_smajic <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Smajic.h5Seurat")    
ncol(seu_smajic)
DefaultAssay(seu_smajic)
unique(seu_smajic@meta.data$Disease_Status) 

## create cell type list
celltype_list <- c("astro", "endo", "mg", "olig", "opc", "peri", "nonda")

for(i in celltype_list) {
    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## Prep Wang
    ## only keep celltype of interest
    # Main meta data
    seu_1 <- seu_wang
    unique(seu_1@meta.data$predicted.id)
    
    ## Subset to only include designated cell type
    unique(seu_1@meta.data$Cell_Type)

    xx <- unique(seu_1@meta.data$Cell_Type)
    xx <- xx[xx %in% i]

    Idents(seu_1) <- "Cell_Type"
    seu_1=subset(seu_1,idents=xx)

    unique(seu_1@meta.data$Cell_Type) 
    nrow(seu_1@meta.data)

    DefaultAssay(seu_1) <- "RNA"
    
    seu_wang_final <- seu_1


    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## Kamath
    ## only keep celltype of interest
    # Main meta data
    seu_1 <- seu_kamath

    ## Subset to only include astro
    unique(seu_1@meta.data$Cell_Type)

    xx <- unique(seu_1@meta.data$Cell_Type)
    xx <- xx[xx %in% i]

    Idents(seu_1) <- "Cell_Type"
    seu_1=subset(seu_1,idents=xx)

    unique(seu_1@meta.data$Cell_Type) 
    nrow(seu_1@meta.data)

    DefaultAssay(seu_1) <- "RNA"

    seu_kamath_final <- seu_1

    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## Smajic
    ## only keep celltype of interest
    # Main meta data
    seu_1 <- seu_smajic

    ## Subset to only include astro
    unique(seu_1@meta.data$Cell_Type)

    xx <- unique(seu_1@meta.data$Cell_Type)
    xx <- xx[xx %in% i]

    Idents(seu_1) <- "Cell_Type"
    seu_1=subset(seu_1,idents=xx)

    unique(seu_1@meta.data$Cell_Type) 
    nrow(seu_1@meta.data)

    DefaultAssay(seu_1) <- "RNA"

    seu_smajic_final <- seu_1

    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## Integrate
    DefaultAssay(seu_kamath_final)
    DefaultAssay(seu_wang_final)
    DefaultAssay(seu_smajic_final)
    ncol(seu_kamath_final)
    ncol(seu_wang_final)
    ncol(seu_smajic_final)

    seu_kamath_final@meta.data$dataset <- "kamath"
    seu_wang_final@meta.data$dataset <- "wang"
    seu_smajic_final@meta.data$dataset <- "smajic"

    ## select inetgration features
    seu_list <- c(seu_kamath_final, seu_wang_final, seu_smajic_final)
    
    seu_int <- Merge_Seurat_List(
        list_seurat = seu_list,
        merge.data = T,
        project = "MergeSeurat"
        )

    ## Need to renormalize and scale. 
    DefaultAssay(seu_int)
    seu_int <- NormalizeData(seu_int)
    seu_int <- FindVariableFeatures(seu_int)
    seu_int <- ScaleData(seu_int)

    ## fix metadata
    unique(seu_int@meta.data$Cell_Type)
    unique(seu_int@meta.data$Disease_Status)

    seu_int@meta.data$Disease_Status[seu_int@meta.data$Disease_Status == "Parkinson disease"] <- "PD"
    seu_int@meta.data$Disease_Status[seu_int@meta.data$Disease_Status == "PD"] <- "PD"
    seu_int@meta.data$Disease_Status[seu_int@meta.data$Disease_Status == "normal"] <- "ctrl"
    seu_int@meta.data$Disease_Status[seu_int@meta.data$Disease_Status == "control"] <- "ctrl"
    seu_int@meta.data$Disease_Status[seu_int@meta.data$Disease_Status == "ctrl"] <- "ctrl"

    unique(seu_int@meta.data$Disease_Status)
    
    ## DEG with MAST test

    par_statistical_method <- "MAST"
        
    Idents(seu_int) <- "Disease_Status"   
    DGE <- FindMarkers(seu_int, ident.1 = "ctrl", ident.2 = "PD",  logfc.threshold = 0, test.use = par_statistical_method)
    
    #write dge
    write.csv(DGE, file = paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/',i,'_intergrated_MAST_DGE.csv'), quote = FALSE, sep = ",")

    ## save RDS
    saveRDS(seu_int, paste0("/home/fiorini9/scratch/machine_learning/top_genes_outs/",i,"_MASTER_LIME_weighted.rds"))
}


############################################################################ Scatter plot
############################################################################
############################################################################
############################################################################
prep <- function(celltype2){
    ## Load LIME objects
    master_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_master_LIME_whole_',celltype2,'.csv'), header = T, sep = ",")
    master_list <- master_list %>% dplyr::select(feature, Mean_z_score_weighted)
    smajic_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_smajic_whole_', celltype2,'.csv'), header = T, sep = ",")
    smajic_list <- smajic_list %>% dplyr::select(feature, Mean_z_score_weighted)
    smajic_list <- aggregate(. ~ feature, data = smajic_list, FUN = function(x) x[1])
    colnames(smajic_list) <- c('feature', 'Smajic_z_score_weighted')

    ## Merge LIME outputs
    merge <- merge(master_list, smajic_list, by = "feature")

    ## Label retained genes
    meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/temp_outs/Master_gene_sets.csv'), header = T, sep = ",")
    meta_data <- subset(meta_data, cell_type == celltype2)

    merge$included <- "No"
    merge$included[merge$feature %in% meta_data$feature] <- "Yes"

    ## Remove genes deemed not biologically relevant
    merge <- subset(merge, feature != "XIST")
    merge <- subset(merge, feature != "TSIX")
    merge <- merge %>% filter(!grepl("^MT-|^LINC", feature)) 
    table(merge$included)
    merge$label <- NA
    merge$label[ merge$included == "Yes"] <- merge$feature[ merge$included == "Yes"]
    merge <- subset(merge, included == "Yes")

    ## Plot
    LIME <- ggplot(merge, aes(y = as.numeric(Smajic_z_score_weighted), x = as.numeric(Mean_z_score_weighted), colour = included, label = label)) +
    geom_point() +
    geom_text_repel(aes(label = label), size = 2.5,
                box.padding = 0.2,
                segment.color = 'grey',
                segment.size = 0.2,
                colour = "black") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background =element_rect(fill=c("white"), colour = "white"),
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text( colour = c("black"), size = 11),
    axis.text.x = element_text( colour = c("black"), size = 11),
    axis.ticks.y = element_blank()) +
    labs(colour = "LIME optimal gene set", x = "Mean LIME feature importance Z-score\n(Kamath et al. and Wang et al.)", y = "LIME feature importance Z-score\n(Smajic et al.)") +
    scale_colour_manual(values = c("orange")) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/scatter_whole_',celltype2,'.pdf',sep=""),width = 7, height = 4)
}
    
## Astro
celltype2 = "astro"

############################################################################ Volcano-type plot
############################################################################
############################################################################
############################################################################
prep <- function(celltype2, top){
    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## prep kamath
    # Test annData object meta data
    meta_data <- paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/HVG_test_meta_',celltype2,'_kamath.csv')
    meta_data <- read.delim(meta_data, header = T, sep = ",")

    test <- meta_data
    
    ## Fold 
    DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_fold',top,'_kamath_',celltype2,'.csv')
    DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_cell_index_fold',top,'_kamath_',celltype2,'.csv')
    LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
    
    LIME$test <- str_count(LIME$X0, ' ')

    LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
    LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])

    ## merge with cell index and remove incorrectly classified cells
    index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
    LIME <- LIME %>% dplyr::select(X0, X1, cell_index, feature )
    nrow(LIME)
    
    index <- index %>% dplyr::select(Disease_Status, predicted_label, cell_index)
    
    LIME <- merge(LIME, index, by = "cell_index", all = T)
    nrow(LIME) 

    ## only keep correctly classified instances
    LIME <- LIME[LIME$Disease_Status == LIME$predicted_label,]
    
    LIME$rep <- top

    test$X %in% LIME$cell_index

    test_2 <- subset(test, X %in% LIME$cell_index)
    test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
    nrow(test_2)
    test_2 <- test_2[!duplicated(test_2$X),]
    nrow(test_2)

    LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
    nrow(LIME_merge) == nrow(LIME)

    table(LIME_merge$Disease_Status, LIME_merge$Sample_ID)
    
    ########################################################################################################################
    ## weigh LIME feature importance scores by the percent of cells that express the gene: Kamath
    ########################################################################################################################

    #### Compute sample specific Z scores
    ## calculate the average importance for each gene  
    LIME_merge_PD <- subset(LIME_merge, Disease_Status == 1)
    
    ## read in pre-computed percent expression
    express_fill <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_kamath_',celltype2,'.csv')
    express_fill <- read.delim(express_fill, header = T, sep = ",")
    express_fill <- express_fill %>% dplyr::select(feature, percent_expressed, normalized_percent_expressed)
    express_fill <- aggregate(. ~ feature, data = express_fill, FUN = function(x) x[1])

    ## Create filler frame
    feature = "fill"
    Mean_feature_importance= "fill"
    Sample_ID= "fill"
    percent_expressed= "fill"
    normalized_percent_expressed= "fill"
    Mean_feature_importance_weighted= "fill"
    z_scores_weighted= "fill"

    filler <- data.frame(feature, Mean_feature_importance, Sample_ID, percent_expressed,normalized_percent_expressed,Mean_feature_importance_weighted, z_scores_weighted)

    for (i in unique(LIME_merge_PD$Sample_ID)) {
        temper <- subset(LIME_merge_PD, Sample_ID == i)
        
        ## Compute average feature importance across cells
        LIME_avg <- temper %>%
        group_by(feature) %>%
        summarize(Mean_feature_importance = mean(X1, na.rm=TRUE))
        LIME_avg$Sample_ID <- i
        nrow(LIME_avg)

        ## weigh by percent expression
        LIME_avg <- merge(LIME_avg,express_fill, by = "feature")
        nrow(LIME_avg)
        LIME_avg$Mean_feature_importance_weighted <- LIME_avg$Mean_feature_importance*LIME_avg$normalized_percent_expressed

        ## compute Z scores
        data <- abs(LIME_avg$Mean_feature_importance_weighted)
        z_scores <- (data-mean(data))/sd(data)
        z_scores

        LIME_avg$z_scores_weighted <- z_scores

        filler <- rbind(filler, LIME_avg) 
    }     
    
    kamath_filler <- subset(filler, feature != "fill")

    kamath_filler <- kamath_filler %>% dplyr::select(feature, z_scores_weighted)
    kamath_list <- kamath_filler %>%
        group_by(feature) %>%
        summarize(Kamath_z_score_weighted = mean(as.numeric(z_scores_weighted), na.rm=TRUE))
    kamath_list <- data.frame(kamath_list)
    subset(kamath_list, feature == "ABCA10")

    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## prep wang
    # Test annData object meta data
    meta_data <- paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/HVG_test_meta_',celltype2,'_wang.csv')
    meta_data <- read.delim(meta_data, header = T, sep = ",")

    test <- meta_data
    
    ## Fold 1
    DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_LIME_fold',top,'_wang_',celltype2,'.csv')
    DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/full_DNN_cell_index_fold',top,'_wang_',celltype2,'.csv')
    LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
    
    LIME$test <- str_count(LIME$X0, ' ')

    LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
    LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])

    ## merge with cell index and remove incorrectly classified cells
    index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
    LIME <- LIME %>% dplyr::select(X0, X1, cell_index, feature )
    nrow(LIME) 
    
    index <- index %>% dplyr::select(Disease_Status, predicted_label, cell_index)
    
    LIME <- merge(LIME, index, by = "cell_index", all = T)
    nrow(LIME) 

    ## only keep correctly classified instances
    LIME <- LIME[LIME$Disease_Status == LIME$predicted_label,]
    
    LIME$rep <- top

    test$X %in% LIME$cell_index

    test_2 <- subset(test, X %in% LIME$cell_index)
    test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
    nrow(test_2)
    test_2 <- test_2[!duplicated(test_2$X),]
    nrow(test_2)

    LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
    nrow(LIME_merge) == nrow(LIME)

    table(LIME_merge$Disease_Status, LIME_merge$Sample_ID)
    
    ########################################################################################################################
    ## weigh LIME feature importance scores by the percent of cells that express the gene: Wang
    ########################################################################################################################

    #### Compute sample specific Z scores
    ## calculate the average importance for each gene  
    LIME_merge_PD <- subset(LIME_merge, Disease_Status == 1)
    
    ## read in pre-computed percent expression
    express_fill <- paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_wang_',celltype2,'.csv')
    express_fill <- read.delim(express_fill, header = T, sep = ",")
    express_fill <- express_fill %>% dplyr::select(feature, percent_expressed, normalized_percent_expressed)
    express_fill <- aggregate(. ~ feature, data = express_fill, FUN = function(x) x[1])

    ## Create filler frame
    feature = "fill"
    Mean_feature_importance= "fill"
    Sample_ID= "fill"
    percent_expressed= "fill"
    normalized_percent_expressed= "fill"
    Mean_feature_importance_weighted= "fill"
    z_scores_weighted= "fill"

    filler <- data.frame(feature, Mean_feature_importance, Sample_ID, percent_expressed,normalized_percent_expressed,Mean_feature_importance_weighted, z_scores_weighted)

    for (i in unique(LIME_merge_PD$Sample_ID)) {
        temper <- subset(LIME_merge_PD, Sample_ID == i)
        
        ## Compute average feature importance across cells
        LIME_avg <- temper %>%
        group_by(feature) %>%
        summarize(Mean_feature_importance = mean(X1, na.rm=TRUE))
        LIME_avg$Sample_ID <- i
        nrow(LIME_avg)

        ## weigh by percent expression
        LIME_avg <- merge(LIME_avg,express_fill, by = "feature")
        nrow(LIME_avg)
        LIME_avg$Mean_feature_importance_weighted <- LIME_avg$Mean_feature_importance*LIME_avg$normalized_percent_expressed

        ## compute Z scores
        data <- abs(LIME_avg$Mean_feature_importance_weighted)
        z_scores <- (data-mean(data))/sd(data)
        z_scores

        LIME_avg$z_scores_weighted <- z_scores

        filler <- rbind(filler, LIME_avg) 
    }     
    
    wang_filler <- subset(filler, feature != "fill")
    wang_filler <- wang_filler %>% dplyr::select(feature, z_scores_weighted)
    wang_list <- wang_filler %>%
        group_by(feature) %>%
        summarize(Wang_z_score_weighted = mean(as.numeric(z_scores_weighted), na.rm=TRUE))
    wang_list <- data.frame(wang_list)
    subset(wang_list, feature == "ABCA10")


    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## prep Smajic
    smajic_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_fold1_smajic_whole_', celltype2,'.csv'), header = T, sep = ",")
    smajic_list <- smajic_list %>% dplyr::select(feature, Mean_z_score_weighted)
    smajic_list <- aggregate(. ~ feature, data = smajic_list, FUN = function(x) x[1])
    colnames(smajic_list) <- c('feature', 'Smajic_z_score_weighted')


    ########################################
    ########################################
    ########################################
    ########################################
    ######################################## combo all 
    
    merge <- merge(kamath_list, wang_list, by = "feature", all.x = T, all.y = T)
    merge <- merge(merge, smajic_list, by = "feature", all.x = T, all.y = T)

    ## label retained genes
    meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/temp_outs/Master_gene_sets.csv'), header = T, sep = ",")
    meta_data <- subset(meta_data, cell_type == celltype2)

    merge$included <- "No"
    merge$included[merge$feature %in% meta_data$feature] <- "Yes"

    
    ######## filter out genes we elected to remove
    remove_genes <- c("LINC00499", "MEG3", "TTTY14", "UTY", "NLGN4Y", "MTRNR2L1", "XIST", "TSIX", "CCDC26"  )
    
    merge <- merge[!(merge$feature %in% remove_genes),]

    filtered_genes <- merge$feature
    filtered_genes <- filtered_genes[!grepl("^MT-", filtered_genes)]
    filtered_genes <- filtered_genes[!grepl("^RP", filtered_genes)]

    merge <- subset(merge, feature %in% filtered_genes)
    ######## ######## ######## ######## ######## ######## ######## ########

    table(merge$included)
    merge$label <- NA
    merge$label[ merge$included == "Yes"] <- merge$feature[ merge$included == "Yes"]

    ## Integrated DEG
    temp1 <- paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/',celltype2,'_intergrated_MAST_DGE.csv')
    temp1 <- read.delim(temp1, header = T, sep = ",")

    bind_FC <- temp1
    bind_FC$mean_FC <- bind_FC$avg_log2FC

    bind_FC <- data.frame(bind_FC)

    bind_FC <- subset(bind_FC, X %in% merge$feature)

    nrow(bind_FC) == nrow(merge)
    
    bind <- merge(merge, bind_FC, by.x = "feature", by.y = "X")

    nrow(bind_FC) == nrow(bind)

    bind$mean_z_score_weighted <- rowMeans(bind[, c("Kamath_z_score_weighted", "Wang_z_score_weighted", "Smajic_z_score_weighted")], na.rm = TRUE)

    #############################################
    ### All genes
    #############################################
    bind$label[bind$included == "Yes"] <- bind$feature[bind$included == "Yes"]

    bind$colour <- "lightgrey"
    bind$colour[as.numeric(-bind$mean_FC) < -0.25] <- "dodgerblue1"
    bind$colour[as.numeric(-bind$mean_FC) > -0.25] <- "indianred2"
    bind$colour[as.numeric(bind$mean_FC) > -0.25 & as.numeric(bind$mean_FC) < 0.25] <- "lightgrey"

    blups <- rev(brewer.pal(11, "RdBu"))

    bind$label <- NA
    bind$label[bind$included == "Yes"] <- bind$feature[bind$included == "Yes"]

    LIME <- ggplot(bind, aes(y = as.numeric(mean_z_score_weighted), x = as.numeric(-mean_FC), colour = colour, label = label, size = abs(mean_FC*mean_z_score_weighted), alpha = included)) +
    geom_point() +
    geom_text_repel(aes(label = label), size = 2,
                box.padding = 0.5,
                segment.color = 'black',
                segment.size = 0.5,
                colour = "black") +
    theme_bw() +
    geom_vline(xintercept = 0.25, linetype = "dashed") +
    geom_vline(xintercept = -0.25, linetype = "dashed") +
    theme(panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background =element_rect(fill=c("white"), colour = "white"),
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text( colour = c("black"), size = 11),
    axis.text.x = element_text( colour = c("black"), size = 11),
    axis.ticks.y = element_blank()) +
    labs(colour = "Mean LIME feature importance\nacross all datasets", x = "Log2FC across all datasets\n", y = "Mean LIME feature importance\nacross all datasets") +
    scale_size(range = c(1, 5)) +
    scale_colour_manual(values = c("#2166AC","#B2182B","lightgrey", "orange" )) +
    scale_alpha_manual(values = c(0.1, 1))
    ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/volcano_whole_', celltype2,'.pdf',sep=""),width = 7, height = 4)

    #############################################
    ### only LIME identified genes
    #############################################
    
    temp1 <- subset(bind, included == "Yes")

    LIME <- ggplot(temp1, aes(y = as.numeric(mean_z_score_weighted), x = as.numeric(-mean_FC), colour = colour, label = label, alpha = included)) +
    geom_point(size = 2.5) +
    geom_text_repel(aes(label = label), size = 2,
                box.padding = 0.5,
                segment.color = 'black',
                segment.size = 0.5,
                colour = "black") +
    theme_bw() +
    geom_vline(xintercept = 0.25, linetype = "dashed") +
    geom_vline(xintercept = -0.25, linetype = "dashed") +
    theme(panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background =element_rect(fill=c("white"), colour = "white"),
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text( colour = c("black"), size = 11),
    axis.text.x = element_text( colour = c("black"), size = 11),
    axis.ticks.y = element_blank()) +
    labs(colour = "Mean LIME feature importance\nacross all datasets", x = "Log2FC across all datasets\n", y = "Mean LIME feature importance Z-score\nacross all datasets") +
    scale_size(range = c(1, 5)) +
    scale_colour_manual(values = c("#2166AC","#B2182B","grey", "orange" )) +
    scale_alpha_manual(values = c(1, 1))
    ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/volcano_lim_', celltype2,'.pdf',sep=""),width = 7, height = 4)

    ## Retain genes with abs(Log2FC >=0.25)
    temp2 <- subset(temp1, abs(mean_FC) >= 0.25)
    nrow(temp2)  

    ## print whole
    write.csv(bind, paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/New_focused_', celltype2, '_whole_FC.csv'))

    ## print lim
    write.csv(temp2, paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/New_focused_', celltype2, '_lim_FC.csv'))
}

## Astro
celltype2 = "astro"
top = 1

############################################################################ Heatmap
############################################################################
############################################################################
############################################################################
prep <- function(keep_cell){
    ## Create cell type list
    celltypelist <- c("astro", "endo", "mg", "olig", "opc", "peri", "nonda")

    ## Load in FC of all cell types for particular optimal gene set
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

    ## keep gene for designated cell type
    j = keep_cell
    fill_cell <- subset(fill, celltype == j)
    fill_lim <- subset(fill, feature %in% fill_cell$feature)
    
    
    ## Factorize cell types and order
    fill_lim$celltype <- factor(fill_lim$celltype, levels = c("astro", "endo", "mg", "olig",  "opc", "peri", "nonda", "da"))
    fill_cell <- fill_cell[order(fill_cell$mean_z_score_weighted),]
    fill_lim$feature <- factor(fill_lim$feature, levels = unique(fill_cell$feature))

    ## add missing entries -- to ensure that tiles aren't grey for the cell types in which the optimal gene set was not identified as HVG
    # Create a dataframe with all combinations of feature and celltype
    all_combinations <- expand.grid(
    feature = unique(fill_lim$feature),
    celltype = unique(fill_lim$celltype),
    stringsAsFactors = FALSE
    )

    # Merge with the existing dataframe, filling in missing values with NA
    result_df <- all_combinations %>%
    left_join(fill_lim, by = c("feature", "celltype"))

    ## Set color scale
    blups <- brewer.pal(9, "Oranges")
    
    # Plot    
    LIME <- ggplot(result_df, aes(x = feature, y = celltype, fill = as.numeric(mean_z_score_weighted))) +
    geom_tile(colour = "white", size = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 8, colour = c("black"),),
        axis.text.y = element_text(colour = c("black"), size = 8)) +
    labs(fill = "") +
    scale_fill_gradientn(colors = blups, na.value = "grey") +
    scale_y_discrete(expand = c(0,0), labels = c("astro" = "Astrocytes", "all" = "All cells", "endo" = "Endothelial", "mg" = "Microglia", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "nonda" = "Neurons", "da" = "DaNeurons")) +
    scale_x_discrete(expand = c(0,0)) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/heatmap_',keep_cell,'.pdf',sep=""),width = 13, height = 2.25)

    ## check how many meet the following criteria: Z score in astro > 5 and appear in > 3 cell types with Z > 2
    criteria <- subset(fill_lim, celltype == keep_cell & mean_z_score_weighted >= 5 )
    Z_keep <- criteria$feature
    Z_keep <- as.character(Z_keep)
    
    criteria <- subset(fill_lim, celltype == keep_cell & mean_z_score_weighted >= 2 )
    criteria_2 <- as.character(criteria$feature)
    criteria <- fill_lim
    criteria <- subset(criteria, feature %in% criteria_2)

    # Count occurrences of each feature
    feature_counts <- table(criteria$feature)

    # Identify features with at least 3 entries
    features_to_keep <- names(feature_counts[feature_counts >= 3])
    total_keep <- c(Z_keep, features_to_keep, PD_genes)
    total_keep <- unique(total_keep)

    # Filter the dataframe
    keep_genes <- subset(fill_lim, feature %in% total_keep)

    # Print the filtered dataframe
    print(length(unique(keep_genes$feature)))
    write.csv(keep_genes, file = paste0('/home/fiorini9/scratch/machine_learning/top_genes_outs/table_genes_',keep_cell,'.csv'), quote = FALSE, sep = ",")
}

## Astro
keep_cell <- "astro"
prep(keep_cell)


