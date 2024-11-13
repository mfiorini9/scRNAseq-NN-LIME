## This script was used to produce Figure 6.

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

############################################################################ DaN volcano-type plot
############################################################################
############################################################################
############################################################################

## Parameters
celltype2 = "da"
top = 1

## Load in LIME outputs
merge <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/temp_outs/LIME_weighted_kamath_percent_ranked_',celltype2,'.csv'), header = T, sep = ",")

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

## FC
temp1 <- paste0('/home/fiorini9/scratch/machine_learning/DEG/MAST_kamath/Cell_based_celltype_groups/',celltype2,'/',celltype2,'_DEG.csv')
temp1 <- read.delim(temp1, header = T, sep = ",")

bind_FC <- temp1
bind_FC$mean_FC <- bind_FC$avg_log2FC

bind_FC <- data.frame(bind_FC)

bind_FC <- subset(bind_FC, X %in% merge$feature)

nrow(bind_FC) == nrow(merge)

bind <- merge(merge, bind_FC, by.x = "feature", by.y = "X")

nrow(bind_FC) == nrow(bind)

bind$mean_z_score_weighted <- bind$Mean_z_score_weighted

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

LIME_identified_genes <- c("GRID2", "HSP90AB1", "SEMA5A", "AK5", 'CACNA2D3', 'CDH8', 'CNTN5', 'GNAS', 'GPC6', 'LPHN2', 'PAM', 'PLCB4', 'PTPRD', 'RSPO2', 'SEZ6L', 'ST6GALNAC5', 'TIAM2', 'UNC5C')
bind$shape <- "circle"
bind$shape[bind$feature %in% LIME_identified_genes] <- "triangle"


LIME <- ggplot(bind, aes(y = as.numeric(mean_z_score_weighted), x = as.numeric(-mean_FC), colour = colour, label = label, alpha = included, shape = shape)) +
geom_point() +
geom_text_repel(aes(label = label), size = 2,
            box.padding = 0.5,
            segment.color = 'black',
            segment.size = 0.5,
            colour = "black") +
theme_classic() +
geom_vline(xintercept = 0.25, linetype = "dashed") +
geom_vline(xintercept = -0.25, linetype = "dashed") +
theme(panel.grid = element_blank(),
legend.position = "none",
axis.title.y = element_text(face = "bold"),
axis.title.x = element_text(face = "bold"),
plot.title = element_text(face = "bold", hjust = 0.5),
strip.background =element_rect(fill=c("white"), colour = "white"),
legend.title = element_text(face = "bold"),
axis.text.y = element_text( colour = c("black"), size = 10),
axis.text.x = element_text( colour = c("black"), size = 10),
axis.ticks.y = element_blank()) +
labs(colour = "LIME feature importance Z-score", x = "Log2FC", y = "LIME feature importance Z-score") +
scale_size(range = c(1, 5)) +
scale_colour_manual(values = c("#2166AC","#B2182B","lightgrey", "orange" )) +
scale_alpha_manual(values = c(0.3, 1)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/volcano_whole_', celltype2,'.pdf',sep=""),width = 2.75, height = 4)


############################################################################ PPMI UMAP 
############################################################################
############################################################################
############################################################################

## Seurat
seu_foudin <- readRDS('/home/fiorini9/nearline/rrg-tdurcan/DATA/FoundIN-PD/download/FOUNDIN/processed/SCRN/SCRN/iPSCsDopaALL_integratedAfterBroadCellType.RDS')
ncol(seu_foudin)
DefaultAssay(seu_foudin) <- "RNA"

## base UMAP plot
p1 <- DimPlot(seu_foudin, reduction = "umap", group.by = "CellType", label = T, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
theme(axis.text.x = element_text(size =10),
axis.text.y = element_text(size =10),
axis.title = element_text(face="bold", size =10),
legend.text = element_text( size =10),
plot.title = element_text( size =10)) + 
xlab("UMAP1") + ylab("UMAP2")
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))

seu_1 <- seu_foudin

## fix cell type label to reflect broad cell types
seu_foudin@meta.data$CellType <- as.character(seu_foudin@meta.data$CellType)
unique(seu_foudin@meta.data$CellType)
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "iDA1"] <- "DA"
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "iDA2"] <- "DA"
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "iDA3"] <- "DA"
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "iDA4"] <- "iDA"
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "eProg1"] <- "eprog"
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "eProg2"] <- "eprog"
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "lProg1"] <- "lProg"
seu_foudin@meta.data$CellType[seu_foudin@meta.data$CellType == "lProg2"] <- "lProg"

## base UMAP plot
p1 <- DimPlot(seu_foudin, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
theme_void()+ 
theme(legend.position = "none",
plot.title = element_blank()) +
scale_color_manual(values = c("eprog" = "#1F78B4",
"Ependymal"="navy",
"PFPP"="goldenrod4",
"NE"="#A6CEE3",
"lProg"="#33A02C",
"DA"="darkgoldenrod1",
"iDA"="firebrick3",
"peri" = "#B15928",
"epen" = "#FFFF99"
))+
xlab("UMAP1") + ylab("UMAP2")
ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/temp.pdf',sep=""))


############################################################################ DaN heatmap 
############################################################################
############################################################################
############################################################################

LIME_identified_genes <- c("GRID2", "HSP90AB1", "SEMA5A", "AK5", 'CACNA2D3', 'CDH8', 'CNTN5', 'GNAS', 'GPC6', 'LPHN2', 'PAM', 'PLCB4', 'PTPRD', 'RSPO2', 'SEZ6L', 'ST6GALNAC5', 'TIAM2', 'UNC5C')

####################################
## iPSC DA
####################################
DA_DEG <- read.delim(paste0('/home/fiorini9/scratch/machine_learning/gene_lists/DA_DEG.csv'), header = T, sep = ",")
DA_DEG <- subset(DA_DEG, X %in% LIME_identified_genes)
blups <- rev(brewer.pal(11, "RdBu"))
min <- min(-DA_DEG$avg_log2FC)
max <- max(-DA_DEG$avg_log2FC)

DA_DEG <- DA_DEG[order(-DA_DEG$avg_log2FC),]

DA_DEG$X <- factor(DA_DEG$X, levels = unique(DA_DEG$X))

DA_DEG$label <- NA
DA_DEG$label[DA_DEG$p_val_adj <= 0.05] <- "*"

Expression <- ggplot(DA_DEG, aes(x = X, y = "iPSC DaNeurons", fill = as.numeric(-avg_log2FC), label = label)) +
geom_tile() +
geom_text(vjust = 0.75)+
theme_bw() +
theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background =element_rect(fill=c("white"), colour = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = c("black"), face = "italic"),
        axis.text.y = element_text(colour = "black")) +
labs(fill = "")  +
scale_fill_gradientn(colors = blups, values = scales::rescale(c(min,  0, max)), na.value = "white")+
scale_y_discrete(expand = c(0,0), labels=c("astro" = "Astrocytes", "all" = "All cells",  "endo" = "Endothelial", "mg" = "Microglia", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "nonda" = "Neurons", "da" = "DaNeurons")) +
scale_x_discrete(expand = c(0,0)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),width = 7, height = 1.5)

####################################
## Kamath DaNs
####################################
i = "da"
DA_kam_DEG <- paste0('/home/fiorini9/scratch/machine_learning/DEG/MAST_kamath/Cell_based_celltype_groups/',i,'/',i,'_DEG.csv')
DA_kam_DEG <- read.delim(DA_kam_DEG, header = T, sep = ",")
DA_kam_DEG <- subset(DA_kam_DEG, X %in% LIME_identified_genes)
blups <- rev(brewer.pal(11, "RdBu"))
min <- min(-DA_kam_DEG$avg_log2FC)
max <- max(-DA_kam_DEG$avg_log2FC)

DA_kam_DEG <- DA_kam_DEG[order(-DA_kam_DEG$avg_log2FC),]

DA_kam_DEG$X <- factor(DA_kam_DEG$X, levels = unique(DA_kam_DEG$X))

DA_kam_DEG$label <- NA
DA_kam_DEG$label[DA_kam_DEG$p_val_adj <= 0.05] <- "*"

####################################
## Combined heatmap
####################################
DA_kam_DEG$celltype <- "DaNeurons"
DA_DEG$celltype <- "iPSC DaNeurons"

merge <- rbind(DA_kam_DEG, DA_DEG)

Expression <- ggplot(merge, aes(y = X, x = celltype, fill = as.numeric(-avg_log2FC), label = label)) +
geom_tile() +
geom_text(vjust = 0.75)+
theme_bw() +
theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        #axis.text.y = element_blank(),
        strip.background =element_rect(fill=c("white"), colour = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = c("black"), size = 11),
        axis.text.y = element_text(colour = "black", face = "italic"), size = 11) +
labs(fill = "")  +
scale_fill_gradientn(colors = blups, values = scales::rescale(c(min,  0, max)), na.value = "white")+
scale_y_discrete(expand = c(0,0), labels=c("astro" = "Astrocytes", "all" = "All cells",  "endo" = "Endothelial", "mg" = "Microglia", "olig" = "Oligo.", "opc" = "OPC", "peri" = "Pericyte", "nonda" = "Neurons", "da" = "DaNeurons")) +
scale_x_discrete(expand = c(0,0)) 
ggsave(paste('/home/fiorini9/scratch/machine_learning/figure_4/temp.pdf',sep=""),width = 3.25, height = 5.5)


############################################################################ ALSKP plot
############################################################################
############################################################################
############################################################################
GPC6 <- fread('/Users/mfiorini/Desktop/machine_learning/from_scratch/Burden_analysis/GPC6_Rare_Variants.csv', sep = ',', header = TRUE)

## Set phenptype factor level
GPC6$Category <- factor(GPC6$Category, levels =c("IMMUNOLOGICAL", "METABOLITE" , "MUSCULOSKELETAL", "NEUROLOGICAL", "PSYCHIATRIC", "SLEEP AND CIRCADIAN", "STROKE"))
GPC6 <- GPC6 %>%
  mutate(logP = -log10(Pvalue)) %>%
  arrange(Category, desc(logP)) %>%
  mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)))

# Reorder by Category
category_order <- c("IMMUNOLOGICAL", "METABOLITE", "MUSCULOSKELETAL", 
                    "NEUROLOGICAL", "PSYCHIATRIC", "SLEEP AND CIRCADIAN", "STROKE")

GPC6 <- GPC6 %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category, desc(logP)) %>%
  mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)))

GPC6 <- GPC6 %>%
  group_by(Category) %>%
  mutate(Xlab = if_else(Pvalue == min(Pvalue), Category, NA_character_)) %>%
  ungroup()

ggplot(GPC6, aes(x = Phenotype, y = -log10(Pvalue), shape = Direction, fill = Category, label = Phenotype)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour="goldenrod1")+
  geom_hline(yintercept = -log10(6.57e-7), linetype = "dashed", colour="goldenrod1")+
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"))+
  scale_shape_manual(values = c(25, 24)) +
  scale_fill_manual(values = c("orchid3", "red", "gold", "blue", "chartreuse3","mediumpurple3", "palevioletred2"))

############################################################################ GPC6 rare variant enrichement
############################################################################
############################################################################
############################################################################

####################################
## Prepare for VEP
####################################
### read in the files downloaded from PD variant browser database
# Define the path to your folder
path <- "/Users/mfiorini/Desktop/machine_learning/from_scratch/Burden_analysis/PD_genetics_data_download"

# List all CSV files in the folder
files <- list.files(path = path, pattern = "*.csv", full.names = TRUE)

# Read all CSV files into a list of dataframes
df_list <- lapply(files, read.csv)

# Combine all dataframes in the list into one dataframe
combined_df <- do.call(rbind, df_list)
combined_df <- combined_df %>% dplyr::select(SNP, ac_case_GENOME, ac_control_GENOME, af_case_GENOME, af_control_GENOME, ac_case_EXOME, ac_control_EXOME, af_case_EXOME, af_control_EXOME)

## PD Genome dataset
PD_genome <- combined_df[!is.na(combined_df$ac_case_GENOME) & 
                  !is.na(combined_df$ac_control_GENOME) & 
                  !is.na(combined_df$af_case_GENOME) & 
                  !is.na(combined_df$af_control_GENOME), ]
variants_PD_genome <- unique(PD_genome$SNP)

## PDGSC dataset
PDGSC <- combined_df[!is.na(combined_df$ac_case_EXOME) | 
                       !is.na(combined_df$ac_control_EXOME) |
                       !is.na(combined_df$af_case_EXOME) |
                       !is.na(combined_df$af_control_EXOME), ]
variants_PDGSC<- unique(PDGSC$SNP)

## prepare for VEP input
combined_df <- combined_df %>%
  separate(SNP, into = c("chromosome", "position", "ref_allele", "alt_allele"), sep = ":", remove = FALSE)
combined_df$variant <- paste(combined_df$chromosome, " ", combined_df$position, " ", combined_df$SNP, " ", combined_df$ref_allele, " ", combined_df$alt_allele, " . . .")
GPC6_forVEP <- combined_df %>% dplyr::select(variant)
nrow(GPC6_forVEP) 

## Write
write.table(GPC6_forVEP$variant, "/Users/mfiorini/Desktop/machine_learning/from_scratch/Burden_analysis/Round2/GPC6_forVEP.csv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


## read in the annotate file
GPC6_annotated <- fread('/Users/mfiorini/Desktop/machine_learning/from_scratch/Burden_analysis/Round2/GPC6_VEPoutput.txt', sep='\t', header=TRUE)

####################################
## Add VEP annotations
####################################
## merge VEP annotate with initial dataframe
GPC6_total <- merge(combined_df, GPC6_annotated, by.x = "SNP", by.y = "#Uploaded_variation")
nrow(combined_df) == nrow(GPC6_total) ## TRUE

## Remove unecessary columns
colnames(GPC6_total)
GPC6_total <- GPC6_total %>% select(SNP, chromosome, position, ref_allele, alt_allele, ac_case_GENOME, ac_control_GENOME, af_case_GENOME,
                                    af_control_GENOME, ac_case_EXOME, ac_control_EXOME, af_case_EXOME, af_control_EXOME, 
                                    Location, Allele, Consequence, SYMBOL, BIOTYPE, REF_ALLELE, UPLOADED_ALLELE, SIFT, PolyPhen, AF, CLIN_SIG,
                                    SOMATIC, MOTIF_NAME, MOTIF_POS, HIGH_INF_POS, MOTIF_SCORE_CHANGE, TRANSCRIPTION_FACTORS, SpliceAI_pred_DP_AG,
                                    SpliceAI_pred_DP_AL, SpliceAI_pred_DP_DG, SpliceAI_pred_DP_DL, SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG,
                                    SpliceAI_pred_DS_DL, SpliceAI_pred_SYMBOL, CADD_PHRED, am_class, am_pathogenicity, gnomADe_AF, gnomADg_AF,Existing_variation)

GPC6_total <- GPC6_total[!duplicated(GPC6_total$SNP), ]

GPC6_total <- GPC6_total %>% rename(gene = SYMBOL)

unique(GPC6_total$Consequence)

## Bin variants by consequence
GPC6_total$Consequence_Bins <- ifelse(GPC6_total$Consequence == "missense_variant" | GPC6_total$Consequence == "missense_variant,splice_region_variant" | 
                                        GPC6_total$Consequence == "inframe_insertion" | GPC6_total$Consequence == "inframe_deletion", "missense", 
                                      ifelse(GPC6_total$Consequence == "frameshift_variant" | GPC6_total$Consequence == "splice_acceptor_variant"
                                             | GPC6_total$Consequence == "splice_donor_variant" | GPC6_total$Consequence == "stop_gained", "lof","synonymous"))
table(GPC6_total$Consequence, GPC6_total$Consequence_Bins)

## We will filter by allele frequency
GPC6_total$gnomADe_AF_numeric <- as.numeric(GPC6_total$gnomADe_AF)

## remove entires > 0.01 in gnomAD all ancestries from VEP
GPC6_total_common <- subset(GPC6_total, gnomADe_AF_numeric > 0.01)
common_name <- unique(GPC6_total_common$SNP)

GPC6_total <- GPC6_total[!(GPC6_total$SNP %in% common_name),]

## Call SNP name to fit the code
colnames(GPC6_total)
colnames(GPC6_total)[1] <- "name"

## From here on out, we will seprate the two databases and perform independent burden test. 

######################################## PD GENOME
########################################
########################################
########################################
######################################## 
## N subjects
n_case <- 2745*2
n_control <- 4071*2
  
GPC6_total_PD_genome <- subset(GPC6_total, name %in% variants_PD_genome)
nrow(GPC6_total_PD_genome) 

## Compute the frequency
temp <- GPC6_total_PD_genome %>% dplyr::select(name, ac_case_GENOME, ac_control_GENOME)

# Separate 'ac_case_GENOME' into three new columns
df_separated <- temp %>%
  separate(ac_case_GENOME, into = c("homo_case", "hetero_case", "non_case"), sep = "/")

# Separate 'ac_control_GENOME' into three new columns
df_separated <- df_separated %>%
  separate(ac_control_GENOME, into = c("homo_ctrl", "hetero_ctrl", "non_ctrl"), sep = "/")

## merge back
merge_temp <- merge(temp, df_separated, by = "name")

## compute allele counts
merge_temp$homo_case <- as.numeric(merge_temp$homo_case)
merge_temp$homo_case <- merge_temp$homo_case*2

merge_temp$homo_ctrl <- as.numeric(merge_temp$homo_ctrl)
merge_temp$homo_ctrl <- merge_temp$homo_ctrl*2

## compute total
# Ensure columns are numeric
merge_temp$homo_case <- as.numeric(merge_temp$homo_case)
merge_temp$hetero_case <- as.numeric(merge_temp$hetero_case)
merge_temp$homo_ctrl <- as.numeric(merge_temp$homo_ctrl)
merge_temp$hetero_ctrl <- as.numeric(merge_temp$hetero_ctrl)

# Add new columns
merge_temp$case_count <- merge_temp$homo_case + merge_temp$hetero_case
merge_temp$control_count <- merge_temp$homo_ctrl + merge_temp$hetero_ctrl


merge_temp <- merge_temp %>% dplyr::select(name, control_count, case_count)
colnames(merge_temp) <- c("name", "Controls", "Cases") 

GPC6_total_PD_genome <- merge(GPC6_total_PD_genome, merge_temp, by = "name")

GPC6_total_PD_genome$am_class[GPC6_total_PD_genome$am_class == "ambiguous"] <- "pathogenic"

## this finds the count per gene for each variant type
geneset_enrichment_counts_PM <- function(newdf){
  # extract rows from als_file containing genes from the gene category
  # create column for synonymous counts
  syn_case <- newdf %>% 
    filter(Consequence_Bins == "synonymous") %>% 
    summarize(syn_case = sum(Cases))
  syn_control <- newdf %>% 
    filter(Consequence_Bins == "synonymous") %>% 
    summarize(syn_control = sum(Controls))
  miss_case <- newdf %>% 
    filter(Consequence_Bins == "missense") %>% 
    summarize(miss_case = sum(Cases))
  miss_control <- newdf %>% 
    filter(Consequence_Bins == "missense") %>% 
    summarize(miss_control = sum(Controls)) 
  miss_dam_case <- newdf %>% 
    filter(am_class == "pathogenic", Consequence_Bins == "missense") %>% 
    summarize(missense_damaging_case = sum(Cases)) 
  miss_dam_control <- newdf %>% 
    filter(am_class == "pathogenic", Consequence_Bins == "missense") %>% 
    summarize(missense_damaging_control = sum(Controls)) 
  miss_benign_case <- newdf %>% 
    filter(am_class == "benign", Consequence_Bins == "missense") %>% 
    summarize(missense_benign_case = sum(Cases)) 
  miss_benign_control <- newdf %>% 
    filter(am_class == "benign",  Consequence_Bins == "missense") %>% 
    summarize(missense_benign_control = sum(Controls)) 
  ptv_case <- newdf %>% 
    filter(Consequence_Bins == "lof") %>% 
    summarize(ptv_case = sum(Cases)) 
  ptv_control <- newdf %>% 
    filter(Consequence_Bins == "lof") %>% 
    summarize(ptv_control = sum(Controls)) 
  merged_df <- merge(syn_case, syn_control, all = TRUE) %>% 
    merge(miss_case, all= TRUE) %>%
    merge(miss_control, all= TRUE) %>%
    merge(miss_dam_case, all = TRUE) %>%
    merge(miss_dam_control, all = TRUE) %>%
    merge(miss_benign_case, all = TRUE) %>%
    merge(miss_benign_control, all = TRUE) %>%
    merge(ptv_case, all= TRUE) %>% 
    merge(ptv_control, all= TRUE) %>% 
    replace(is.na(.), 0)
  return(merged_df)
}

GPC6_total_PD_genome$CADD_PHRED <- as.numeric(GPC6_total_PD_genome$CADD_PHRED)
SUMO_Genes_PM_counts <- geneset_enrichment_counts_PM(GPC6_total_PD_genome)


####################################
## Perform Fisher's exact test
####################################
#syn fisher test 
fisher_syn_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$syn_case[i]
    control_carrier <- df$syn_control[i]
    case_noncarrier <- n_case - df$syn_case[i]
    control_noncarrier <- n_control - df$syn_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$syn_pval <- pval
  df$syn_OR <- OR
  df$syn_lowci <- lowci
  df$syn_highci <- highci
  df[order(df$syn_pval),]}

SUMO_Genes_PM_Fisher_Syn <- fisher_syn_PM(SUMO_Genes_PM_counts)

#missense fisher test 
fisher_miss_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$miss_case[i]
    control_carrier <- df$miss_control[i]
    case_noncarrier <- n_case - df$miss_case[i]
    control_noncarrier <- n_control - df$miss_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$miss_pval <- pval
  df$miss_OR <- OR
  df$miss_lowci <- lowci
  df$miss_highci <- highci
  df[order(df$miss_pval),]}

SUMO_Genes_PM_Fisher_Miss <- fisher_miss_PM(SUMO_Genes_PM_counts)

#missense damaging fisher test 
fisher_miss_dam_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$missense_damaging_case[i]
    control_carrier <- df$missense_damaging_control[i]
    case_noncarrier <- n_case - df$missense_damaging_case[i]
    control_noncarrier <- n_control - df$missense_damaging_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$missense_damaging_pval <- pval
  df$missense_damaging_OR <- OR
  df$missense_damaging_lowci <- lowci
  df$missense_damaging_highci <- highci
  df[order(df$missense_damaging_pval),]}


SUMO_Genes_PM_Fisher_MissDamg <- fisher_miss_dam_PM(SUMO_Genes_PM_counts)

#missense benign fisher test 
fisher_miss_ben_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$missense_benign_case[i]
    control_carrier <- df$missense_benign_control[i]
    case_noncarrier <- n_case - df$missense_benign_case[i]
    control_noncarrier <- n_control - df$missense_benign_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$missense_benign_pval <- pval
  df$missense_benign_OR <- OR
  df$missense_benign_lowci <- lowci
  df$missense_benign_highci <- highci
  df[order(df$missense_benign_pval),]}

SUMO_Genes_PM_Fisher_MissBen <- fisher_miss_ben_PM(SUMO_Genes_PM_counts)

#ptv fisher test 
fisher_ptv_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$ptv_case[i]
    control_carrier <- df$ptv_control[i]
    case_noncarrier <- n_case - df$ptv_case[i]
    control_noncarrier <- n_control - df$ptv_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$ptv_pval <- pval
  df$ptv_OR <- OR
  df$ptv_lowci <- lowci
  df$ptv_highci <- highci
  df[order(df$ptv_pval),]}

SUMO_Genes_PM_Fisher_PTV <- fisher_ptv_PM(SUMO_Genes_PM_counts)

## merging the different fisher results together
SUMO_Genes_PM_Fisher <- merge(SUMO_Genes_PM_Fisher_Syn, SUMO_Genes_PM_Fisher_Miss, all = TRUE) %>%
  merge(SUMO_Genes_PM_Fisher_MissBen, all = TRUE) %>%
  merge(SUMO_Genes_PM_Fisher_MissDamg, all = TRUE) %>%
  merge(SUMO_Genes_PM_Fisher_PTV, all = TRUE) %>%
  mutate(geneset = "Project Mine ALS")

## convert to long format
# Pvalue
pval_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("pval"))

# Convert to long format with column names
pval_long <- pval_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "pval_value")

pval_long <- data.frame(pval_long)
pval_long <- pval_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

pval_long <- ## convert to long format
  ### Pvalue
  pval_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("pval"))

# Convert to long format with column names
pval_long <- pval_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "pval_value")

pval_long <- data.frame(pval_long)
pval_long <- pval_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

pval_long <- pval_long %>% dplyr::select(-column_name)

## OR
OR_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("OR"))

# Convert to long format with column names
OR_long <- OR_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "OR")

OR_long <- data.frame(OR_long)
OR_long <- OR_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

OR_long <- OR_long %>% dplyr::select(-column_name)

## lowci
lowci_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("lowci"))

# Convert to long format with column names
lowci_long <- lowci_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "lowci")

lowci_long <- data.frame(lowci_long)
lowci_long <- lowci_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

lowci_long <- lowci_long %>% dplyr::select(-column_name)

## highci
highci_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("highci"))

# Convert to long format with column names
highci_long <- highci_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "highci")

highci_long <- data.frame(highci_long)
highci_long <- highci_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

highci_long <- highci_long %>% dplyr::select(-column_name)

## merge them all
merge <- merge(pval_long, OR_long, by = "type" )
merge <- merge(merge, lowci_long, by = "type" )
merge <- merge(merge, highci_long, by = "type" )

merge <- subset(merge, type != "ptv")

## set factor level
merge$type <- factor(merge$type, levels = c("synonymous", "missense", "missense_benign", "missense_damaging"))

## plot
ggplot(data = merge, aes(x=type, y=OR, ymin=lowci, ymax=highci)) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=lowci, ymax=highci), width=0.1, cex=1, colour =ifelse(merge$pval_value < 0.05, "black", "darkgrey")) +
  geom_point(size = 3, aes(col="black")) +
  coord_flip(clip = "off")+
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text= element_text(size=12), axis.text = element_text(size=12, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18), axis.title.y = element_blank()) +
  geom_text(aes(y=15, label=ifelse(is.na(pval_value), "", formatC(pval_value, format = "e", digits = 2))), size=4, hjust=0)+
  theme(panel.spacing.x = unit(5, "lines"))+
  scale_color_manual(name="Odds Ratio", breaks = c("OR > 1", "OR < 1"), values = c("#F2846D", "#729EFC")) +
  scale_y_continuous(trans="log10",labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(labels = c("synonymous" = "Synonymous", 
                              "missense" = "Missense",
                              "missense_benign" = "Missense benign",
                              "missense_damaging" = "Missense pathogenic"
  ))

PD_GENOME_counts <- SUMO_Genes_PM_Fisher
PD_genome_total <- GPC6_total_PD_genome


######################################## IPDGC
########################################
########################################
########################################
######################################## 
## N subjects
n_case <- 2110*2
n_control <- 2978*2

GPC6_total_PD_genome <- subset(GPC6_total, name %in% variants_PDGSC)
nrow(GPC6_total_PD_genome)

## Compute the frequency
temp <- GPC6_total_PD_genome %>% dplyr::select(name, ac_case_EXOME, ac_control_EXOME)

# Separate 'ac_case_EXOME' into three new columns
df_separated <- temp %>%
  separate(ac_case_EXOME, into = c("homo_case", "hetero_case", "non_case"), sep = "/")

# Separate 'ac_control_EXOME' into three new columns
df_separated <- df_separated %>%
  separate(ac_control_EXOME, into = c("homo_ctrl", "hetero_ctrl", "non_ctrl"), sep = "/")

## merge back
merge_temp <- merge(temp, df_separated, by = "name")

## compute allele counts
merge_temp$homo_case <- as.numeric(merge_temp$homo_case)
merge_temp$homo_case <- merge_temp$homo_case*2

merge_temp$homo_ctrl <- as.numeric(merge_temp$homo_ctrl)
merge_temp$homo_ctrl <- merge_temp$homo_ctrl*2

## compute total
# Ensure columns are numeric
merge_temp$homo_case <- as.numeric(merge_temp$homo_case)
merge_temp$hetero_case <- as.numeric(merge_temp$hetero_case)
merge_temp$homo_ctrl <- as.numeric(merge_temp$homo_ctrl)
merge_temp$hetero_ctrl <- as.numeric(merge_temp$hetero_ctrl)

# Add new columns
merge_temp$case_count <- merge_temp$homo_case + merge_temp$hetero_case
merge_temp$control_count <- merge_temp$homo_ctrl + merge_temp$hetero_ctrl


merge_temp <- merge_temp %>% dplyr::select(name, control_count, case_count)
colnames(merge_temp) <- c("name", "Controls", "Cases") 

GPC6_total_PD_genome <- merge(GPC6_total_PD_genome, merge_temp, by = "name")

GPC6_total_PD_genome$am_class[GPC6_total_PD_genome$am_class == "ambiguous"] <- "pathogenic"

## this finds the count per gene for each variant type. 
geneset_enrichment_counts_PM <- function(newdf){
  # extract rows from als_file containing genes from the gene category
  # create column for synonymous counts
  syn_case <- newdf %>% 
    filter(Consequence_Bins == "synonymous") %>% 
    summarize(syn_case = sum(Cases))
  syn_control <- newdf %>% 
    filter(Consequence_Bins == "synonymous") %>% 
    summarize(syn_control = sum(Controls))
  miss_case <- newdf %>% 
    filter(Consequence_Bins == "missense") %>% 
    summarize(miss_case = sum(Cases))
  miss_control <- newdf %>% 
    filter(Consequence_Bins == "missense") %>% 
    summarize(miss_control = sum(Controls)) 
  miss_dam_case <- newdf %>% 
    filter(am_class == "pathogenic", Consequence_Bins == "missense") %>% 
    summarize(missense_damaging_case = sum(Cases)) 
  miss_dam_control <- newdf %>% 
    filter(am_class == "pathogenic", Consequence_Bins == "missense") %>% 
    summarize(missense_damaging_control = sum(Controls)) 
  miss_benign_case <- newdf %>% 
    filter(am_class == "benign", Consequence_Bins == "missense") %>% 
    summarize(missense_benign_case = sum(Cases)) 
  miss_benign_control <- newdf %>% 
    filter(am_class == "benign",  Consequence_Bins == "missense") %>% 
    summarize(missense_benign_control = sum(Controls)) 
  ptv_case <- newdf %>% 
    filter(Consequence_Bins == "lof") %>% 
    summarize(ptv_case = sum(Cases)) 
  ptv_control <- newdf %>% 
    filter(Consequence_Bins == "lof") %>% 
    summarize(ptv_control = sum(Controls)) 
  merged_df <- merge(syn_case, syn_control, all = TRUE) %>% 
    merge(miss_case, all= TRUE) %>%
    merge(miss_control, all= TRUE) %>%
    merge(miss_dam_case, all = TRUE) %>%
    merge(miss_dam_control, all = TRUE) %>%
    merge(miss_benign_case, all = TRUE) %>%
    merge(miss_benign_control, all = TRUE) %>%
    merge(ptv_case, all= TRUE) %>% 
    merge(ptv_control, all= TRUE) %>% 
    replace(is.na(.), 0)
  return(merged_df)
}

GPC6_total_PD_genome$CADD_PHRED <- as.numeric(GPC6_total_PD_genome$CADD_PHRED)
SUMO_Genes_PM_counts <- geneset_enrichment_counts_PM(GPC6_total_PD_genome)


####################################
## Perform Fisher's exact test
####################################

#syn fisher test 
fisher_syn_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$syn_case[i]
    control_carrier <- df$syn_control[i]
    case_noncarrier <- n_case - df$syn_case[i]
    control_noncarrier <- n_control - df$syn_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$syn_pval <- pval
  df$syn_OR <- OR
  df$syn_lowci <- lowci
  df$syn_highci <- highci
  df[order(df$syn_pval),]}


SUMO_Genes_PM_Fisher_Syn <- fisher_syn_PM(SUMO_Genes_PM_counts)


#missense fisher test 
fisher_miss_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$miss_case[i]
    control_carrier <- df$miss_control[i]
    case_noncarrier <- n_case - df$miss_case[i]
    control_noncarrier <- n_control - df$miss_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$miss_pval <- pval
  df$miss_OR <- OR
  df$miss_lowci <- lowci
  df$miss_highci <- highci
  df[order(df$miss_pval),]}


SUMO_Genes_PM_Fisher_Miss <- fisher_miss_PM(SUMO_Genes_PM_counts)

#missense damaging fisher test 
fisher_miss_dam_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$missense_damaging_case[i]
    control_carrier <- df$missense_damaging_control[i]
    case_noncarrier <- n_case - df$missense_damaging_case[i]
    control_noncarrier <- n_control - df$missense_damaging_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$missense_damaging_pval <- pval
  df$missense_damaging_OR <- OR
  df$missense_damaging_lowci <- lowci
  df$missense_damaging_highci <- highci
  df[order(df$missense_damaging_pval),]}


SUMO_Genes_PM_Fisher_MissDamg <- fisher_miss_dam_PM(SUMO_Genes_PM_counts)

#missense benign fisher test 
fisher_miss_ben_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$missense_benign_case[i]
    control_carrier <- df$missense_benign_control[i]
    case_noncarrier <- n_case - df$missense_benign_case[i]
    control_noncarrier <- n_control - df$missense_benign_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$missense_benign_pval <- pval
  df$missense_benign_OR <- OR
  df$missense_benign_lowci <- lowci
  df$missense_benign_highci <- highci
  df[order(df$missense_benign_pval),]}

SUMO_Genes_PM_Fisher_MissBen <- fisher_miss_ben_PM(SUMO_Genes_PM_counts)

#ptv fisher test 
fisher_ptv_PM <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$ptv_case[i]
    control_carrier <- df$ptv_control[i]
    case_noncarrier <- n_case - df$ptv_case[i]
    control_noncarrier <- n_control - df$ptv_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$ptv_pval <- pval
  df$ptv_OR <- OR
  df$ptv_lowci <- lowci
  df$ptv_highci <- highci
  df[order(df$ptv_pval),]}

SUMO_Genes_PM_Fisher_PTV <- fisher_ptv_PM(SUMO_Genes_PM_counts)

## merging the different fisher results together
SUMO_Genes_PM_Fisher <- merge(SUMO_Genes_PM_Fisher_Syn, SUMO_Genes_PM_Fisher_Miss, all = TRUE) %>%
  merge(SUMO_Genes_PM_Fisher_MissBen, all = TRUE) %>%
  merge(SUMO_Genes_PM_Fisher_MissDamg, all = TRUE) %>%
  merge(SUMO_Genes_PM_Fisher_PTV, all = TRUE) %>%
  mutate(geneset = "Project Mine ALS")


## convert to long format
### Pvalue
pval_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("pval"))

# Convert to long format with column names
pval_long <- pval_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "pval_value")

pval_long <- data.frame(pval_long)
pval_long <- pval_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

pval_long <- ## convert to long format
  ### Pvalue
  pval_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("pval"))

# Convert to long format with column names
pval_long <- pval_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "pval_value")

pval_long <- data.frame(pval_long)
pval_long <- pval_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

pval_long <- pval_long %>% dplyr::select(-column_name)

## OR
OR_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("OR"))

# Convert to long format with column names
OR_long <- OR_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "OR")

OR_long <- data.frame(OR_long)
OR_long <- OR_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

OR_long <- OR_long %>% dplyr::select(-column_name)

## lowci
lowci_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("lowci"))

# Convert to long format with column names
lowci_long <- lowci_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "lowci")

lowci_long <- data.frame(lowci_long)
lowci_long <- lowci_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

lowci_long <- lowci_long %>% dplyr::select(-column_name)

## highci
highci_columns <- SUMO_Genes_PM_Fisher %>%
  select(contains("highci"))

# Convert to long format with column names
highci_long <- highci_columns %>%
  pivot_longer(everything(), names_to = "column_name", values_to = "highci")

highci_long <- data.frame(highci_long)
highci_long <- highci_long %>%
  mutate(type = case_when(
    grepl("syn", column_name) ~ "synonymous",
    grepl("missense_benign", column_name) ~ "missense_benign",
    grepl("missense_damaging", column_name) ~ "missense_damaging",
    grepl("miss", column_name) ~ "missense",
    grepl("ptv", column_name) ~ "ptv"
  ))

highci_long <- highci_long %>% dplyr::select(-column_name)

## merge them all
merge <- merge(pval_long, OR_long, by = "type" )
merge <- merge(merge, lowci_long, by = "type" )
merge <- merge(merge, highci_long, by = "type" )

merge <- subset(merge, type != "ptv")

## set factor level
merge$type <- factor(merge$type, levels = c("synonymous", "missense", "missense_benign", "missense_damaging"))

## plot
ggplot(data = merge, aes(x=type, y=OR, ymin=lowci, ymax=highci)) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=lowci, ymax=highci), width=0.1, cex=1, colour =ifelse(merge$pval_value < 0.05, "black", "darkgrey")) +
  geom_point(size = 3, aes(col="black")) +
  coord_flip(clip = "off")+
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text= element_text(size=12), axis.text = element_text(size=12, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18), axis.title.y = element_blank()) +
  geom_text(aes(y=45, label=ifelse(is.na(pval_value), "", formatC(pval_value, format = "e", digits = 2))), size=4, hjust=0)+
  theme(panel.spacing.x = unit(5, "lines"))+
  scale_color_manual(name="Odds Ratio", breaks = c("OR > 1", "OR < 1"), values = c("#F2846D", "#729EFC")) +
  scale_y_continuous(trans="log10",labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(labels = c("synonymous" = "Synonymous", 
                              "missense" = "Missense",
                              "missense_benign" = "Missense benign",
                              "missense_damaging" = "Missense pathogenic"
  ))


PDGSC_counts <- SUMO_Genes_PM_Fisher
PDGSC_total <- GPC6_total_PD_genome


######################################## CMH
########################################
########################################
########################################
######################################## 
## create a filler dataframe
type = "fill"
OR = "fill"
lowci = "fill"
highci = "fill"
pval_value = "fill"

fill = data.frame(type, OR, lowci, highci, pval_value)

## missense damaging
PD_genome_case_carrier <- PD_GENOME_counts$missense_damaging_case
PD_genome_control_carrier <- PD_GENOME_counts$missense_damaging_control
PD_genome_case_noncarrier <- (2745*2) - PD_GENOME_counts$missense_damaging_case
PD_genome_control_noncarrier <- (4071*2) - PD_GENOME_counts$missense_damaging_control

PDGSC_case_carrier <- PDGSC_counts$missense_damaging_case
PDGSC_control_carrier <- PDGSC_counts$missense_damaging_control
PDGSC_case_noncarrier <- (2110*2) - PDGSC_counts$missense_damaging_case
PDGSC_control_noncarrier <- (2978*2) - PDGSC_counts$missense_damaging_control
ALStest <- array(c(PD_genome_case_carrier, PD_genome_control_carrier, PD_genome_case_noncarrier, PD_genome_control_noncarrier, 
                   PDGSC_case_carrier, PDGSC_control_carrier, PDGSC_case_noncarrier, PDGSC_control_noncarrier),
                 dim = c(2, 2, 2),
                 dimnames = list(
                   Data_type = c("Cases", "Controls"),
                   Status = c("Carrier", "NC"),
                   Status_2 = c("PD_genome", "PDGSC")))
mod <- mantelhaen.test(ALStest)

type = "missense_damaging"
OR = mod$estimate
lowci = mod$conf.int[1]
highci = mod$conf.int[2]
pval_value = mod$p.value

temp = data.frame(type, OR, lowci, highci, pval_value)
fill <- rbind(fill, temp)

## missense benign
PD_genome_case_carrier <- PD_GENOME_counts$missense_benign_case
PD_genome_control_carrier <- PD_GENOME_counts$missense_benign_control
PD_genome_case_noncarrier <- (2745*2) - PD_GENOME_counts$missense_benign_case
PD_genome_control_noncarrier <- (4071*2) - PD_GENOME_counts$missense_benign_control

PDGSC_case_carrier <- PDGSC_counts$missense_benign_case
PDGSC_control_carrier <- PDGSC_counts$missense_benign_control
PDGSC_case_noncarrier <- (2110*2) - PDGSC_counts$missense_benign_case
PDGSC_control_noncarrier <- (2978*2) - PDGSC_counts$missense_benign_control
ALStest <- array(c(PD_genome_case_carrier, PD_genome_control_carrier, PD_genome_case_noncarrier, PD_genome_control_noncarrier, 
                   PDGSC_case_carrier, PDGSC_control_carrier, PDGSC_case_noncarrier, PDGSC_control_noncarrier),
                 dim = c(2, 2, 2),
                 dimnames = list(
                   Data_type = c("Cases", "Controls"),
                   Status = c("Carrier", "NC"),
                   Status_2 = c("PD_genome", "PDGSC")))
mod <- mantelhaen.test(ALStest)

type = "missense_benign"
OR = mod$estimate
lowci = mod$conf.int[1]
highci = mod$conf.int[2]
pval_value = mod$p.value

temp = data.frame(type, OR, lowci, highci, pval_value)
fill <- rbind(fill, temp)

## missense 
PD_genome_case_carrier <- PD_GENOME_counts$miss_case
PD_genome_control_carrier <- PD_GENOME_counts$miss_control
PD_genome_case_noncarrier <- (2745*2) - PD_GENOME_counts$miss_case
PD_genome_control_noncarrier <- (4071*2) - PD_GENOME_counts$miss_control

PDGSC_case_carrier <- PDGSC_counts$miss_case
PDGSC_control_carrier <- PDGSC_counts$miss_control
PDGSC_case_noncarrier <- (2110*2) - PDGSC_counts$miss_case
PDGSC_control_noncarrier <- (2978*2) - PDGSC_counts$miss_control
ALStest <- array(c(PD_genome_case_carrier, PD_genome_control_carrier, PD_genome_case_noncarrier, PD_genome_control_noncarrier, 
                   PDGSC_case_carrier, PDGSC_control_carrier, PDGSC_case_noncarrier, PDGSC_control_noncarrier),
                 dim = c(2, 2, 2),
                 dimnames = list(
                   Data_type = c("Cases", "Controls"),
                   Status = c("Carrier", "NC"),
                   Status_2 = c("PD_genome", "PDGSC")))
mod <- mantelhaen.test(ALStest)

type = "missense"
OR = mod$estimate
lowci = mod$conf.int[1]
highci = mod$conf.int[2]
pval_value = mod$p.value

temp = data.frame(type, OR, lowci, highci, pval_value)
fill <- rbind(fill, temp)

## syn 
PD_genome_case_carrier <- PD_GENOME_counts$syn_case
PD_genome_control_carrier <- PD_GENOME_counts$syn_control
PD_genome_case_noncarrier <- (2745*2) - PD_GENOME_counts$syn_case
PD_genome_control_noncarrier <- (4071*2) - PD_GENOME_counts$syn_control

PDGSC_case_carrier <- PDGSC_counts$syn_case
PDGSC_control_carrier <- PDGSC_counts$syn_control
PDGSC_case_noncarrier <- (2110*2) - PDGSC_counts$syn_case
PDGSC_control_noncarrier <- (2978*2) - PDGSC_counts$syn_control
ALStest <- array(c(PD_genome_case_carrier, PD_genome_control_carrier, PD_genome_case_noncarrier, PD_genome_control_noncarrier, 
                   PDGSC_case_carrier, PDGSC_control_carrier, PDGSC_case_noncarrier, PDGSC_control_noncarrier),
                 dim = c(2, 2, 2),
                 dimnames = list(
                   Data_type = c("Cases", "Controls"),
                   Status = c("Carrier", "NC"),
                   Status_2 = c("PD_genome", "PDGSC")))
mod <- mantelhaen.test(ALStest)

type = "synonymous"
OR = mod$estimate
lowci = mod$conf.int[1]
highci = mod$conf.int[2]
pval_value = mod$p.value

temp = data.frame(type, OR, lowci, highci, pval_value)
fill <- rbind(fill, temp)

fill <- subset(fill, OR != "fill")
fill$OR <- as.numeric(fill$OR)
fill$lowci <- as.numeric(fill$lowci)
fill$highci <- as.numeric(fill$highci)
fill$pval_value <- as.numeric(fill$pval_value)

## need to apply an BH p-value correction
fill$p_adjust <- p.adjust(fill$pval_value, method = "BH")

## set factor
fill$type <- factor(fill$type, levels = c("synonymous", "missense", "missense_benign", "missense_damaging"))

## plot
ggplot(data = fill, aes(x=type, y=OR, ymin=lowci, ymax=highci)) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=lowci, ymax=highci), width=0.1, cex=1, colour =ifelse(fill$pval_value < 0.05, "black", "darkgrey")) +
  geom_point(size = 3, aes(col="black")) +
  coord_flip(clip = "off")+
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text= element_text(size=12), axis.text = element_text(size=14, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18), axis.title.y = element_blank()) +
  geom_text(aes(y=15, label=ifelse(is.na(pval_value), "", formatC(pval_value, format = "e", digits = 2))), size=4.5, hjust=0)+
  theme(panel.spacing.x = unit(5, "lines"))+
  scale_color_manual(name="Odds Ratio", breaks = c("OR > 1", "OR < 1"), values = c("#F2846D", "#729EFC")) +
  scale_y_continuous(trans="log10",labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(labels = c("synonymous" = "Synonymous", 
                              "missense" = "Missense",
                              "missense_benign" = "Missense benign",
                              "missense_damaging" = "Missense pathogenic"
  ))




