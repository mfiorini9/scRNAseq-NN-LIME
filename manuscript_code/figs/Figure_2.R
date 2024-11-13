## This script was used to produce Figure 2.
salloc -A def-sfarhan --time=0-5 -c 1 --mem=40g

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
library(MetaNeighbor)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(reshape2)
library(ggplot2)

############################################################################ UMAPs
############################################################################
############################################################################
############################################################################
##################
## Kamath
##################
seu_kam <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Kamath.h5Seurat")    

## UMAP by cell type
DefaultAssay(seu_kam) <- "integrated"
    
## print UMAP
DimPlot(seu_kam, reduction = "umap", group.by = "Cell_Type", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
theme_void()+ 
theme(legend.position = "none",
plot.title = element_blank()) +
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

##################
## Wang
##################
seu_wang <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Wang.h5Seurat")    

## UMAP by cell type
DefaultAssay(seu_wang) <- "integrated"
    
## print UMAP
DimPlot(seu_wang, reduction = "umap", group.by = "Cell_Type", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
theme_void()+ 
theme(legend.position = "none",
plot.title = element_blank()) +
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

##################
## smajic
##################
seu_smaj <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Smajic.h5Seurat")    

## UMAP by cell type
DefaultAssay(seu_smaj) <- "integrated"

## print UMAP
DimPlot(seu_smaj, reduction = "umap", group.by = "Celltype", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
theme_void()+ 
theme(legend.position = "none",
plot.title = element_blank()) +
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


############################################################################ n cells and prop in PD and control
############################################################################
############################################################################
############################################################################
##################
## Kamath
##################
df <- data.frame(seu_kam@meta.data)

df2 <- data.frame(table(df$Disease_Status, df$Cell_Type))
df3 <- data.frame(table(df$Cell_Type))
total <- sum(df3$Freq)
df3$Proportion <- df3$Freq / total

df2$Var2 <- factor(df2$Var2, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))
df3$Var1 <- factor(df3$Var1, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))

## Number cells
counts <- ggplot(df3, aes(fill=Var1, y=Proportion, x=Var1)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_classic() + 
    xlab("") + ylab("Number of cells\n(x 1000)") +
    theme(
        axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size =12, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
    scale_x_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte")) +
    scale_y_continuous( limits = c(0,0.7)) +
    coord_flip() +
    scale_fill_manual(values = c("astro" = "#CAB2D6",
        "endo"="#FDBF6F",
        "olig"="#B2DF8A",
        "mg"="#FB9A99",
        "da"="#A6CEE3",
        "nonda"="#1F78B4",
        "opc"="#33A02C",
        "peri" = "#B15928",
        "epen" = "#FFFF99"))
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

## case control proportion
prop <- ggplot(df2, aes( y = Var2, fill = Var1, x = Freq)) +
    geom_bar(position = "fill", stat = 'identity')+
    theme_classic() + 
    theme(
    axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size =12, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
    ylab("") + xlab("Proportion of\ncells") +
    scale_fill_manual(values = c("dodgerblue2", "indianred3")) +
    scale_x_continuous(expand = c(0.01,0.01), breaks = c(0,0.5,1.0), labels = c("0", "0.5", "1.0")) +
    scale_y_discrete(labels=c("astro" = df3$Freq[df3$Var1 == "astro"], 
                            "da" = df3$Freq[df3$Var1 == "da"], 
                            "endo" = df3$Freq[df3$Var1 == "endo"],
                            "mg" = df3$Freq[df3$Var1 == "mg"], 
                            "nonda" = df3$Freq[df3$Var1 == "nonda"], 
                            "olig" = df3$Freq[df3$Var1 == "olig"], 
                            "opc" = df3$Freq[df3$Var1 == "opc"],
                            "peri" = df3$Freq[df3$Var1 == "peri"]))
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

## combine
ggarrange(counts, prop, ncol = 2, widths = c(3,2.7), align = "h")
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 4.5, height = 3.75)

##################
## Wang
##################
df <- data.frame(seu_wang@meta.data)

unique(df$Cell_Type)

df2 <- data.frame(table(df$Disease_Status, df$Cell_Type))
df3 <- data.frame(table(df$Cell_Type))
total <- sum(df3$Freq)
df3$Proportion <- df3$Freq / total

df2$Var2 <- factor(df2$Var2, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))
df3$Var1 <- factor(df3$Var1, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))

## Number cells
counts <- ggplot(df3, aes(fill=Var1, y=Proportion, x=Var1)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_classic() + 
    xlab("") + ylab("Number of cells\n(x 1000)") +
    theme(
        axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size =12, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
    scale_x_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte")) +
    scale_y_continuous( limits = c(0,0.7)) +
    coord_flip() +
    scale_fill_manual(values = c("astro" = "#CAB2D6",
        "endo"="#FDBF6F",
        "olig"="#B2DF8A",
        "mg"="#FB9A99",
        "da"="#A6CEE3",
        "nonda"="#1F78B4",
        "opc"="#33A02C",
        "peri" = "#B15928",
        "epen" = "#FFFF99"))

ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

## case control proportion
prop <- ggplot(df2, aes( y = Var2, fill = Var1, x = Freq)) +
    geom_bar(position = "fill", stat = 'identity')+
    theme_classic() + 
    theme(
    axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size =12, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
    ylab("") + xlab("Proportion of\ncells") +
    scale_fill_manual(values = c("dodgerblue2", "indianred3")) +
    scale_x_continuous(expand = c(0.01,0.01), breaks = c(0,0.5,1.0), labels = c("0", "0.5", "1.0")) +
    scale_y_discrete(labels=c("astro" = df3$Freq[df3$Var1 == "astro"], 
                            "da" = df3$Freq[df3$Var1 == "da"], 
                            "endo" = df3$Freq[df3$Var1 == "endo"],
                            "mg" = df3$Freq[df3$Var1 == "mg"], 
                            "nonda" = df3$Freq[df3$Var1 == "nonda"], 
                            "olig" = "142494", 
                            "opc" = df3$Freq[df3$Var1 == "opc"],
                            "peri" = df3$Freq[df3$Var1 == "peri"]))

ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

## combine
ggarrange(counts, prop, ncol = 2, widths = c(3,2.6), align = "h")
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 4.5, height = 3.75)

##################
## smajic
##################
df <- data.frame(seu_smaj@meta.data)
df2 <- data.frame(table(df$DiseaseStatus, df$Celltype))
df3 <- data.frame(table(df$Celltype))

total <- sum(df3$Freq)
df3$Proportion <- df3$Freq / total

df2$Var2 <- factor(df2$Var2, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))
df3$Var1 <- factor(df3$Var1, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))

## Number cells
counts <- ggplot(df3, aes(fill=Var1, y=Proportion, x=Var1)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_classic() + 
    xlab("") + ylab("Number of cells\n(x 1000)") +
    theme(
        axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size =12, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
    scale_x_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte")) +

    scale_y_continuous( limits = c(0,0.7)) +
    coord_flip() +
    #scale_fill_manual(values = c("#CAB2D6", "#1F78B4",  "#FDBF6F", "#FB9A99","#A6CEE3",  "#B2DF8A", "#33A02C")) 
    scale_fill_manual(values = c("astro" = "#CAB2D6",
        "endo"="#FDBF6F",
        "olig"="#B2DF8A",
        "mg"="#FB9A99",
        "da"="#A6CEE3",
        "nonda"="#1F78B4",
        "opc"="#33A02C",
        "peri" = "#B15928",
        "epen" = "#FFFF99"))

## case control proportion
prop <- ggplot(df2, aes( y = Var2, fill = Var1, x = Freq)) +
    geom_bar(position = "fill", stat = 'identity')+
    theme_classic() + 
    theme(
    axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size =12, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
    ylab("") + xlab("Proportion of\ncells") +
    scale_fill_manual(values = c("dodgerblue2", "indianred3")) +
    scale_x_continuous(expand = c(0.01,0.01), breaks = c(0,0.5,1.0), labels = c("0", "0.5", "1.0")) +
    scale_y_discrete(labels=c("astro" = df3$Freq[df3$Var1 == "astro"], 
                            "da" = df3$Freq[df3$Var1 == "da"], 
                            "endo" = df3$Freq[df3$Var1 == "endo"], 
                            "epen" = df3$Freq[df3$Var1 == "epen"],
                            "mg" = df3$Freq[df3$Var1 == "mg"], 
                            "nonda" = df3$Freq[df3$Var1 == "nonda"], 
                            "olig" = df3$Freq[df3$Var1 == "olig"], 
                            "opc" = df3$Freq[df3$Var1 == "opc"],
                            "peri" = df3$Freq[df3$Var1 == "peri"]))

ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

## combine
ggarrange(counts, prop, ncol = 2, widths = c(3,2.7), align = "h")
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 5, height = 4 )

############################################################################ Marker gene expression
############################################################################
############################################################################
############################################################################
#####################
## Kam and Wang combo
#####################
## features to plot
par_select_features_list= c("AQP4","CLDN5", "FLT1", "CD74", "C3", "MOBP", "MOG", "GRIK1", "VCAN", "DCC",  "COL1A2", "PDGFRB",  "SLC17A6", "GAD2", "RBFOX3", "RIT2", "GALNTL6", "TH")
DefaultAssay(seu_kam) <- "RNA"

## dotplot
kam_dot_plt <- DotPlot(seu_kam, features = par_select_features_list, group.by = 'Cell_Type')

kam_dot_plt2 <- data.frame(kam_dot_plt$data)
kam_dot_plt2$id <- factor(kam_dot_plt2$id, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))
kam_dot_plt2$features.plot <- factor(kam_dot_plt2$features.plot, levels = c("AQP4","CLDN5", "FLT1", "CD74", "C3", "MOBP", "MOG", "GRIK1", "VCAN", "DCC",  "COL1A2", "PDGFRB",  "SLC17A6", "GAD2", "RBFOX3", "RIT2", "GALNTL6", "TH"))

blups <- rev(brewer.pal(11, "RdBu"))

kam_plot <- ggplot(kam_dot_plt2, aes(x = features.plot , y = id, size = pct.exp, colour = avg.exp.scaled)) + 
geom_point() +
theme_classic() +
theme(
axis.text.x = element_text(angle=45, hjust = 1),
axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.y = element_text(size =12, colour = c("#CAB2D6", "#FDBF6F", "#FB9A99", "#B2DF8A", "#33A02C", "#B15928", "#1F78B4", "#A6CEE3")),
axis.title = element_text(face="bold", size =12),
legend.text = element_text( size =10),
legend.title = element_text( size= 10),
legend.position = "right") +
scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericytes")) +
scale_colour_gradientn(colors = blups, values = scales::rescale(c(-1,  0, 3)), na.value = "white")+
scale_size(guide = guide_legend(direction = "vertical")) +
ylab("") + xlab("") +
guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 10, height = 5)

## Wang 
DefaultAssay(seu_wang) <- "RNA"

## dotplot
wang_dot_plt <- DotPlot(seu_wang, features = par_select_features_list, group.by = 'Cell_Type')

wang_dot_plt2 <- data.frame(wang_dot_plt$data)
wang_dot_plt2$id <- factor(wang_dot_plt2$id, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))
wang_dot_plt2$features.plot <- factor(wang_dot_plt2$features.plot, levels = c("AQP4","CLDN5", "FLT1", "CD74", "C3", "MOBP", "MOG", "GRIK1", "VCAN", "DCC",  "COL1A2", "PDGFRB",  "SLC17A6", "GAD2", "RBFOX3", "RIT2", "GALNTL6", "TH"))

blups <- rev(brewer.pal(11, "RdBu"))

wang_plot <- ggplot(wang_dot_plt2, aes(x = features.plot , y = id, size = pct.exp, colour = avg.exp.scaled)) + 
geom_point() +
theme_classic() +
theme(
axis.text.x = element_text(angle=45, hjust = 1),
axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.y = element_text(size =12, colour = c("#CAB2D6", "#FDBF6F", "#FB9A99", "#B2DF8A", "#33A02C", "#B15928", "#1F78B4", "#A6CEE3")),
axis.title = element_text(face="bold", size =12),
legend.text = element_text( size =10),
legend.title = element_text( size= 10),
legend.position = "right") +
scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericytes")) +
scale_colour_gradientn(colors = blups, values = scales::rescale(c(-1,  0, 3)), na.value = "white")+
scale_size(guide = guide_legend(direction = "vertical")) +
ylab("") + xlab("") +
guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 10, height = 5)

## arrange plots
ggarrange(kam_plot, wang_plot, ncol = 1, nrow = 2, heights = c(1,1))
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 6.7, height = 7.5)

##################
## smajic
##################
## features to plot
par_select_features_list= c("AQP4","CLDN5", "FLT1", "CD74", "C3", "MOBP", "MOG", "GRIK1", "VCAN", "DCC",  "COL1A2", "PDGFRB",  "SLC17A6", "GAD2", "RBFOX3", "RIT2", "GALNTL6", "TH")
DefaultAssay(seu_smaj) <- "RNA"

## dotplot
kam_dot_plt <- DotPlot(seu_smaj, features = par_select_features_list, group.by = 'Cell_Type')

kam_dot_plt2 <- data.frame(kam_dot_plt$data)
kam_dot_plt2$id <- factor(kam_dot_plt2$id, levels = c('astro', 'endo',  'mg',  'olig',  'opc',   'peri', 'nonda', 'da'))
kam_dot_plt2$features.plot <- factor(kam_dot_plt2$features.plot, levels = c("AQP4","CLDN5", "FLT1", "CD74", "C3", "MOBP", "MOG", "GRIK1", "VCAN", "DCC",  "COL1A2", "PDGFRB",  "SLC17A6", "GAD2", "RBFOX3", "RIT2", "GALNTL6", "TH"))

blups <- rev(brewer.pal(11, "RdBu"))

kam_plot <- ggplot(kam_dot_plt2, aes(x = features.plot , y = id, size = pct.exp, colour = avg.exp.scaled)) + 
geom_point() +
theme_classic() +
theme(
axis.text.x = element_text(angle=45, hjust = 1),
axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.y = element_text(size =12, colour = c("#CAB2D6", "#FDBF6F", "#FB9A99", "#B2DF8A", "#33A02C", "#B15928", "#1F78B4", "#A6CEE3")),
axis.title = element_text(face="bold", size =12),
legend.text = element_text( size =10),
legend.title = element_text( size= 10),
legend.position = "right") +
scale_y_discrete(labels=c("astro" = "Astrocyte", "da" = "DaN", "endo" = "Endothelial", "mg" = "Microglia", "nonda" = "nDaN", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericytes")) +
scale_colour_gradientn(colors = blups, values = scales::rescale(c(-1,  0, 3)), na.value = "white")+
scale_size(guide = guide_legend(direction = "vertical")) +
ylab("") + xlab("") +
guides(colour = guide_colorbar("Average expression", title.position = "top"), size = guide_legend("% expressing cells", title.position = "top")) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 6.7, height = 3.75)
##


############################################################################ NN performance boxplots
############################################################################
############################################################################
############################################################################
##################
## kamath
##################
cell_list <- c('astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'da','all')
model <- c("DNN")

X0 <- "fill"
metric <- "fill"
Cell_type <- "fill"
dataset <- "fill"
method <- "fill"
fill <- data.frame(X0,metric, Cell_type, dataset,method)

for (j in unique(model)){
for (i in unique(cell_list)){
temp <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/',j,'_5FC_eval_',i,'_kamath.csv'), header = T, sep = ",")
temp <- subset(temp, metric == "BA")
temp$Cell_type <- i
temp$dataset <- "Kamath"
temp$method <- j
temp <- temp %>% dplyr::select(X0,metric, Cell_type,dataset,method)

fill <- rbind(fill, temp)
}
}

fill <- subset(fill,method != "fill")
kamath_df <- fill

unique(kamath_df$Cell_type)

##################
## wang
##################
cell_list <- c('astro', 'endo', 'mg', 'olig', 'opc', 'peri', 'nonda', 'all')
model <- c("DNN")

X0 <- "fill"
metric <- "fill"
Cell_type <- "fill"
dataset <- "fill"
method <- "fill"
fill <- data.frame(X0,metric, Cell_type, dataset,method)

for (j in unique(model)){
for (i in unique(cell_list)){
temp <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/machine_learning_outs/',j,'_5FC_eval_',i,'_wang.csv'), header = T, sep = ",")
temp <- subset(temp, metric == "BA")
temp$Cell_type <- i
temp$dataset <- "Wang"
temp$method <- j
temp <- temp %>% dplyr::select(X0,metric, Cell_type,dataset,method)

fill <- rbind(fill, temp)
}
}

fill <- subset(fill,method != "fill")
wang_df <- fill

unique(wang_df$Cell_type)

##################
## plot
##################
bind <- rbind(kamath_df, wang_df)
bind$dataset[bind$dataset == "Kamath"] <- "Kamath et al."
bind$dataset[bind$dataset == "Wang"] <- "Wang et al."

mean_results <- bind %>%
group_by(dataset) %>%
summarise(mean_X0 = mean(as.numeric(X0), na.rm = TRUE))

## add Da filler for Wang et al. 
X0 <- 9
metric <- "BA"
Cell_type <- "da"
dataset <- "Wang et al."
method <- "DNN"
filler_da <- data.frame(X0, metric, Cell_type, dataset, method)

bind <- rbind(bind, filler_da)

bind$Cell_type <- factor(bind$Cell_type, levels = c("all", "astro","endo","mg","olig","opc","peri","nonda","da"))

ggplot(bind, aes(x = dataset, y = as.numeric(X0), fill = dataset, colour = dataset)) + 
geom_boxplot( alpha = 0.5) +
theme_bw() +
theme(
axis.title.x = element_blank(),
panel.grid = element_blank(),
axis.title = element_text(face = "bold"),
axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
axis.text.y = element_text(colour = "black"),
strip.background =element_rect(fill=c("white"), colour = "white"),
legend.position = "none") +
ylab("Balanced accracy") + xlab("Machine learning\nclassifier") +
facet_wrap(~Cell_type, ncol = 9, scales = "free_x") +
scale_fill_manual(values = c("#1b9e77", "#e6ab02")) +
scale_colour_manual(values = c("#1b9e77", "#e6ab02")) +
scale_y_continuous(limits = c(0.85, 1))
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""),width = 10, height = 3)    

    

