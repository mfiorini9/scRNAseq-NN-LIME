## This script was used to annotate cells and convert the Seurat objects to AnnData

salloc -A def-sfarhan --time=0-8 -c 1 --mem=200g

module load StdEnv/2020 
module load r/4.2.2 
R

library(Seurat)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(SeuratData, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

################################################################################################################
################################################################################################################
################################################################################################################ Kamath annotation and conversion to AD
## code
    
    ## Seurat Kam
    seu_kam <- readRDS('/home/fiorini9/scratch/machine_learning/scRNAbox_kamath2/scRNAbox/step6/objs6/seu_step6.rds')
    DefaultAssay(seu_kam)
    head(colnames(seu_kam))

    ## UMAP integrated res 0.8
    p1 <- DimPlot(seu_kam, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 5, width =5)

    ## dot plot with Smaj marker
    par_select_features_list= c("AQP4", "CADPS2", "TH", "CLDN5", "CD74", "SLC17A6", "GAD2", "GRIK1", "MOBP", "VCAN")

    ## Define default assay   
    DefaultAssay(seu_kam) <- "RNA"

    ## dotplot
    dot_plt <- DotPlot(seu_kam, features = par_select_features_list, group.by = 'integrated_snn_res.0.8') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size =12),
    axis.title = element_text(face="bold", size =12),
    legend.text = element_text( size =10),
    legend.title = element_text( size= 10),
    legend.position = "top") +
    coord_flip() +
    scale_colour_gradient(low = "white", high = "black") +
    scale_size(guide = guide_legend(direction = "horizontal")) +
    ylab("") + xlab("") +
    guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 10, height = 5)

    ## dot plot with Wang marker
    par_select_features_list= c("SKAP1", "COL1A2", "PDGFRB", "FLT1", "VCAN", "MOG", "C3", "AQP4", "RBFOX3", "RIT2", "GALNTL6", "DCC", "SLC18A2", "TH", "SLC6A3", "GAD2", "SLC17A6")

    ## load parameters    
    DefaultAssay(seu_kam) <- "RNA"

    ## dotplot
    dot_plt <- DotPlot(seu_kam, features = par_select_features_list, group.by = 'integrated_snn_res.0.8') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size =12),
    axis.title = element_text(face="bold", size =12),
    legend.text = element_text( size =10),
    legend.title = element_text( size= 10),
    legend.position = "top") +
    coord_flip() +
    scale_colour_gradient(low = "white", high = "black") +
    scale_size(guide = guide_legend(direction = "horizontal")) +
    ylab("") + xlab("") +
    guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 10, height = 5)

    ## UMAP integrated res 0.8
    p1 <- DimPlot(seu_kam, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 5, width =5)

    ## add first round of clustering labels
    df <- data.frame(seu_kam@meta.data)

    df <- df %>% dplyr::select(integrated_snn_res.0.8)
    df$celltype_1 <- as.character(df$integrated_snn_res.0.8)

    df$celltype_1[df$integrated_snn_res.0.8  == 0] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 1] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8 == 2] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 3] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 4] <- "da"
    df$celltype_1[df$integrated_snn_res.0.8  == 5] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 6] <- "mg"
    df$celltype_1[df$integrated_snn_res.0.8  == 7] <- "astro"
    df$celltype_1[df$integrated_snn_res.0.8  == 8] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 9] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 10] <- "astro"
    df$celltype_1[df$integrated_snn_res.0.8  == 11] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 12] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 13] <- "mg"
    df$celltype_1[df$integrated_snn_res.0.8  == 14] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 15] <- "opc"
    df$celltype_1[df$integrated_snn_res.0.8  == 16] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 17] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 18] <- "peri"
    df$celltype_1[df$integrated_snn_res.0.8  == 19] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 20] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 21] <- "endo"
    df$celltype_1[df$integrated_snn_res.0.8  == 22] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 23] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 24] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 25] <- "peri"
    df$celltype_1[df$integrated_snn_res.0.8  == 26] <- "nonda"
    df$celltype_1[df$integrated_snn_res.0.8  == 27] <- "mg"

    ## Add metadata
    seu_kam <- AddMetaData(seu_kam, df)

    ## UMAP with celltype_1
    p1 <- DimPlot(seu_kam, reduction = "umap", group.by = "celltype_1", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 5, width =5)

    ## fix the meta data columns
    df <- data.frame(seu_kam_lim@meta.data)
    colnames(df)
    df$Cell_Type <- df$celltype_1
    df$Sample_ID <- df$Sample_ID
    df$Disease_Status <- df$disease__ontology_label

    seu_kam_lim <- AddMetaData(seu_kam_lim,df)
    ncol(seu_kam_lim) #300513

    saveRDS(seu_kam_lim, "/home/fiorini9/scratch/machine_learning/temp_seu_objects/temp_seu_kam_lim.Rds")

    ## convert to Ann data
    DefaultAssay(seu_kam_lim)
    DefaultAssay(seu_kam_lim) <- "RNA"
    DefaultAssay(seu_kam_lim)

    SaveH5Seurat(seu_kam_lim, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Kamath.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Kamath.h5Seurat", dest = "h5ad")

    max(seu_kam_lim@assays$RNA@data)
##


################################################################################################################
################################################################################################################
################################################################################################################ Wang annotation and conversion to AD
## code
    ## new Seurat
    seu_wang <- readRDS('/home/fiorini9/scratch/machine_learning/scRNAbox_wang/scRNAbox/step6/objs6/seu_step6.rds')
    DefaultAssay(seu_wang)

    ## UMAP plot by clutser
    p1 <- DimPlot(seu_wang, reduction = "umap", group.by = "integrated_snn_res.2", label = T, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

    ## explore marker gene expression
    par_select_features_list= c("SKAP1", "COL1A2", "PDGFRB", "FLT1", "VCAN", "MOG", "C3", "AQP4", "RBFOX3", "RIT2", "GALNTL6", "DCC", "SLC18A2", "TH", "SLC6A3", "GAD2", "SLC17A6")

    ## confirm default assay
    DefaultAssay(seu_wang) <- "RNA"

    ## dot plot
    dot_plt <- DotPlot(seu_wang, features = par_select_features_list, group.by = 'integrated_snn_res.0.8') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size =12),
    axis.title = element_text(face="bold", size =12),
    legend.text = element_text( size =10),
    legend.title = element_text( size= 10),
    legend.position = "top") +
    scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC")) +
    coord_flip() +
    scale_colour_gradient(low = "white", high = "black") +
    scale_size(guide = guide_legend(direction = "horizontal")) +
    ylab("") + xlab("") +
    guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 5, height = 5)
    

    ## add first round of clustering labels
    df <- data.frame(seu_wang@meta.data)

    df <- df %>% dplyr::select(integrated_snn_res.0.8)
    df$celltype_1 <- as.character(df$integrated_snn_res.0.8)

    df$celltype_1[df$integrated_snn_res.0.8  == 0] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 1] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8 == 2] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 3] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 4] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 5] <- "mg"
    df$celltype_1[df$integrated_snn_res.0.8  == 6] <- "opc"
    df$celltype_1[df$integrated_snn_res.0.8  == 7] <- "astro"
    df$celltype_1[df$integrated_snn_res.0.8  == 8] <- "endo"
    df$celltype_1[df$integrated_snn_res.0.8  == 9] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 10] <- 10
    df$celltype_1[df$integrated_snn_res.0.8  == 11] <- "mg"
    df$celltype_1[df$integrated_snn_res.0.8  == 12] <- "astro"
    df$celltype_1[df$integrated_snn_res.0.8  == 13] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 14] <- 14
    df$celltype_1[df$integrated_snn_res.0.8  == 15] <- "peri"
    df$celltype_1[df$integrated_snn_res.0.8  == 16] <- "peri"
    df$celltype_1[df$integrated_snn_res.0.8  == 17] <- "opc"
    df$celltype_1[df$integrated_snn_res.0.8  == 18] <- "t"
    df$celltype_1[df$integrated_snn_res.0.8  == 19] <- "fibro"
    df$celltype_1[df$integrated_snn_res.0.8  == 20] <- 20
    df$celltype_1[df$integrated_snn_res.0.8  == 21] <- "olig"
    df$celltype_1[df$integrated_snn_res.0.8  == 22] <- "astro"
    df$celltype_1[df$integrated_snn_res.0.8  == 23] <- 23
    df$celltype_1[df$integrated_snn_res.0.8  == 24] <- "mg"

    ## Add metadata
    seu_wang <- AddMetaData(seu_wang, df)

    p1 <- DimPlot(seu_wang, reduction = "umap", group.by = "celltype_1", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 10, width =10)

    ## save a temp RDS file
    saveRDS(seu_wang, "/home/fiorini9/scratch/machine_learning/temp_seu_objects/temp_seu_wang.Rds")

    ## reference based annotation with kam
    ## load name of existing Seurat objects
    seu_wang <- readRDS("/home/fiorini9/scratch/machine_learning/temp_seu_objects/temp_seu_wang.Rds")
    DefaultAssay(seu_wang)

    ## set user defined clustering resolution
    seu_wang <- SetIdent(seu_wang, value = 'integrated_snn_res.0.8')

    ## load reference Seurat object
    reference0 <- readRDS('/home/fiorini9/scratch/PD_machine_lerning/figures/Kamath.rds' )
    DefaultAssay(reference0) <- "RNA" 
    unique(reference0@meta.data$disease__ontology_label)
    unique(reference0@meta.data$Cell_type)
    
    # perform standard preprocessing on reference object
    reference0<- NormalizeData(reference0)
    reference0 <- FindVariableFeatures(reference0)
    reference0<- ScaleData(reference0)
    reference0 <- RunPCA(object = reference0, assay = "RNA", npcs = 20)

    ## find transfer anchors between reference and query Seurat objects
    transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu_wang, dims = 1:20, reference.reduction = "pca")

    ## add reference-based annotations to the qeury object
    eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',"Cell_Type" ,',dims = 1:',20,')', sep='')))
    seu_wang <- AddMetaData(object = seu_wang, metadata = predictions)

    # Add metadata column for reference object
    seu_wang$temp_temp_2 <- seu_wang@meta.data$predicted.id
    name_meta <- names(seu_wang@meta.data) 
    length <- length(name_meta)
    name_meta[length] <- paste("Kamath", "_predictions", sep = "")
    names(seu_wang@meta.data) <- name_meta

    ## Print a umap projection showing the predicted cell types on the query object 
    reference0 <- RunUMAP(reference0, dims = 1:20, reduction = "pca", return.model = TRUE)
    seu_wang <- MapQuery(anchorset = transfer.anchors, reference = reference0, query = seu_wang,
        refdata = list(celltype = "Cell_Type"), reference.reduction = "pca", reduction.model = "umap")
    p1 <- DimPlot(reference0, reduction = "umap", group.by = "Cell_Type", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Reference annotations")
    p2 <- DimPlot(seu_wang, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
    p1 + p2
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 10, width =10)

    ## UMAP integrate 2
    p1 <- DimPlot(seu_wang, reduction = "umap", group.by = "integrated_snn_res.2", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 10, width =10)
    
    min(seu_wang@meta.data$prediction.score.max)

    ## explore marker gene expression
    par_select_features_list= c("SKAP1", "COL1A2", "PDGFRB", "FLT1", "VCAN", "MOG", "C3", "AQP4", "RBFOX3", "RIT2", "GALNTL6", "DCC", "SLC18A2", "TH", "SLC6A3", "GAD2", "SLC17A6")

    ## Default assay
    DefaultAssay(seu_wang) <- "RNA"

    ## dotplot    
    dot_plt <- DotPlot(seu_wang, features = par_select_features_list, group.by = 'predicted.celltype') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size =12),
    axis.title = element_text(face="bold", size =12),
    legend.text = element_text( size =10),
    legend.title = element_text( size= 10),
    legend.position = "top") +
    scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC")) +
    coord_flip() +
    scale_colour_gradient(low = "white", high = "black") +
    scale_size(guide = guide_legend(direction = "horizontal")) +
    ylab("") + xlab("") +
    guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 5, height = 5)

    ## Save seu object again
    saveRDS(seu_wang, "/home/fiorini9/scratch/machine_learning/temp_seu_objects/temp_seu_wang.Rds")

    ## remove clusters we dont want
    seu_wang <- readRDS("/home/fiorini9/scratch/machine_learning/temp_seu_objects/temp_seu_wang.Rds")
    DefaultAssay(seu_wang) <- "integrated"

    ## UMAP 
    p1 <- DimPlot(seu_wang, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 10, width =10)

    ## UMAP
    p1 <- DimPlot(seu_wang_lim, reduction = "umap", group.by = "celltype_1", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 5, width =5)

    ## label the remaing Neuron clusters
    df <- data.frame(seu_wang_lim@meta.data)
    df <- df %>% dplyr::select(celltype_1)
    df$celltype_1[df$celltype_1  == "10"] <- "nonda"
    df$celltype_1[df$celltype_1  == "20"] <- "nonda"
    df$celltype_1[df$celltype_1 == "23"] <- "nonda"
    df$celltype_1[df$celltype_1 == "14"] <- "nonda"

    ## Add metadata
    seu_wang_lim <- AddMetaData(seu_wang_lim, df)

    ## UMAP
    p1 <- DimPlot(seu_wang_lim, reduction = "umap", group.by = "celltype_1", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 5, width =5)

    ## remove fibro, and t
    unique(seu_wang_lim@meta.data$celltype_1)
        
        xx <- unique(seu_wang_lim@meta.data$celltype_1)
        xx <- xx[!(xx %in% c("t", "fibro"))]

        Idents(seu_wang_lim) <- "celltype_1"
        seu_wang_lim=subset(seu_wang_lim,idents=xx)

    unique(seu_wang_lim@meta.data$celltype_1) 
    nrow(seu_wang_lim@meta.data)

    ## UMAP
    p1 <- DimPlot(seu_wang_lim, reduction = "umap", group.by = "celltype_1", label = T, pt.size= 0.000001, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""), height = 5, width =5)

    saveRDS(seu_wang_lim, "/home/fiorini9/scratch/machine_learning/temp_seu_objects/temp_seu_wang_lim.Rds")

    ## convert to Ann data
    DefaultAssay(seu_wang_lim)
    DefaultAssay(seu_wang_lim) <- "RNA"
    DefaultAssay(seu_wang_lim)

    ## fix the meta data columns
    df <- data.frame(seu_wang_lim@meta.data)
    colnames(df)
    df$Cell_Type <- df$celltype_1
    df$Sample_ID <- df$Sample_ID
    df$Disease_Status <- df$orig.ident

    seu_wang_lim <- AddMetaData(seu_wang_lim,df)
    ncol(seu_wang_lim)

    SaveH5Seurat(seu_wang_lim, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Wang.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Wang.h5Seurat", dest = "h5ad")

    max(seu_kam_lim@assays$RNA@data)
##

################################################################################################################
################################################################################################################
################################################################################################################ Smajic annotation and conversion to AD
## code
    ## Seurat
    seu_new <- readRDS('/home/fiorini9/scratch/machine_learning/scRNAbox_smajic/step6/objs6/seu_step6.rds')

    df_new <- data.frame(seu_new@meta.data)
    nrow(df)
    df$names <- rownames(df)
    nrow(df_new)
    df_new$names <- rownames(df_new)

    df_all <- merge(df_new, df, by = "names", all = T)
    df_all <- subset(df_all, names %in% df_new$names )
    nrow(df_all) == nrow(df_new)

    colnames(df_all)
    rownames(df_all) <- df_all$names
    df_all <- df_all %>% dplyr::select(Kamath_celltype)
    table(df_all)

    ## Add metadata
    seu_new <- AddMetaData(seu_new, df_all)

    ## print a UMAP
    p1 <- DimPlot(seu_new, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))


    ## UMAP
    p1 <- DimPlot(seu_new, reduction = "umap", group.by = "Kamath_celltype", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    scale_color_manual(values = c("astro" = "#CAB2D6",
    "endo"="#FDBF6F",
    "olig"="#B2DF8A",
    "mg"="#FB9A99",
    "da"="#A6CEE3",
    "nonda"="#1F78B4",
    "opc"="#33A02C",
    "peri" = "#B15928",
    "ependymal" = "#FFFF99"
    ))+
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

    ## explore marker gene expression
    par_select_features_list= c("SKAP1", "COL1A2", "PDGFRB", "FLT1", "VCAN", "MOG", "C3", "AQP4", "RBFOX3", "RIT2", "GALNTL6", "DCC", "SLC18A2", "TH", "SLC6A3", "GAD2", "SLC17A6")

    ## confirm default assay
    DefaultAssay(seu_new) <- "RNA"

    ## dot plot
    dot_plt <- DotPlot(seu_new, features = par_select_features_list, group.by = 'integrated_snn_res.0.8') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size =12),
    axis.title = element_text(face="bold", size =12),
    legend.text = element_text( size =10),
    legend.title = element_text( size= 10),
    legend.position = "top") +
    scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC")) +
    coord_flip() +
    scale_colour_gradient(low = "white", high = "black") +
    scale_size(guide = guide_legend(direction = "horizontal")) +
    ylab("") + xlab("") +
    guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 5, height = 5)

    ## re annotate
    df <- data.frame(seu_new@meta.data)

    df <- df %>% dplyr::select(integrated_snn_res.0.8)
    df$Celltype[df$integrated_snn_res.0.8 == 0 ] <- "olig"
    df$Celltype[df$integrated_snn_res.0.8 == 1 ] <- "olig"
    df$Celltype[df$integrated_snn_res.0.8 == 2] <- "olig"
    df$Celltype[df$integrated_snn_res.0.8 == 3 ] <- "mg"
    df$Celltype[df$integrated_snn_res.0.8 == 4] <- "olig"
    df$Celltype[df$integrated_snn_res.0.8 == 5] <- "astro"
    df$Celltype[df$integrated_snn_res.0.8 == 6] <- "nonda"
    df$Celltype[df$integrated_snn_res.0.8 == 7] <- "olig"
    df$Celltype[df$integrated_snn_res.0.8 == 8] <- "olig"
    df$Celltype[df$integrated_snn_res.0.8 == 9] <- "opc"
    df$Celltype[df$integrated_snn_res.0.8 == 10] <- "endo"
    df$Celltype[df$integrated_snn_res.0.8 == 11] <- "astro"
    df$Celltype[df$integrated_snn_res.0.8 == 12] <- "peri"
    df$Celltype[df$integrated_snn_res.0.8 == 13] <- "epen"
    df$Celltype[df$integrated_snn_res.0.8 == 17] <- "mg"
    df$Celltype[df$integrated_snn_res.0.8 == 19] <- "da"

    seu_new <- AddMetaData(seu_new, df)

    ## UMAP
    p1 <- DimPlot(seu_new, reduction = "umap", group.by = "Celltype", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    scale_color_manual(values = c("astro" = "#CAB2D6",
    "endo"="#FDBF6F",
    "olig"="#B2DF8A",
    "mg"="#FB9A99",
    "da"="#A6CEE3",
    "nonda"="#1F78B4",
    "opc"="#33A02C",
    "peri" = "#B15928",
    "epen" = "#FFFF99"
    ))+
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

    ## save anndata object
    DefaultAssay(seu_new) <- "RNA"
    DefaultAssay(seu_new)

    SaveH5Seurat(seu_new, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Smajic.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Smajic.h5Seurat", dest = "h5ad")

    max(seu_new@assays$RNA@data)

    ## fix the meta data columns
    seu_smaj <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Smajic.h5Seurat")    
    ncol(seu_smaj) ## 40909 cells base
    DefaultAssay(seu_smaj)

    df <- data.frame(seu_smaj@meta.data)
    colnames(df)
    df$Cell_Type <- df$Celltype
    df$Sample_ID <- df$Sample_ID
    df$Disease_Status[df$orig.ident == "Control1" ] <- "ctrl"
    df$Disease_Status[df$orig.ident == "Control2" ] <- "ctrl"
    df$Disease_Status[df$orig.ident == "Control3" ] <- "ctrl"
    df$Disease_Status[df$orig.ident == "Control4" ] <- "ctrl"
    df$Disease_Status[df$orig.ident == "Control5" ] <- "ctrl"
    df$Disease_Status[df$orig.ident == "Control6" ] <- "ctrl"
    df$Disease_Status[df$orig.ident == "Parkinson1" ] <- "PD"
    df$Disease_Status[df$orig.ident == "Parkinson2" ] <- "PD"
    df$Disease_Status[df$orig.ident == "Parkinson3" ] <- "PD"
    df$Disease_Status[df$orig.ident == "Parkinson4" ] <- "PD"
    df$Disease_Status[df$orig.ident == "Parkinson5" ] <- "PD"

    unique(df$Disease_Status)

    seu_smaj <- AddMetaData(seu_smaj,df)
    ncol(seu_smaj) #40909

    ## save object with fixed metadata
    DefaultAssay(seu_smaj) <- "RNA"
    DefaultAssay(seu_smaj)

    SaveH5Seurat(seu_smaj, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Smajic.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Smajic.h5Seurat", dest = "h5ad")

    max(seu_smaj@assays$RNA@counts)
##

