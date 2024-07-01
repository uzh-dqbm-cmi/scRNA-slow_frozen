library(SingleR)
library(dplyr)
library(tidyr)
library(Seurat)
library(reshape)
library(data.table)
library(readr)
library(hash)
library(ggplot2)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(lattice)
library(patchwork)
library(gtable)
library(el)
library(plotly)
library(Scillus) # for plot_heatmap: multibar heatmaps
library(stringr)
library("org.Hs.eg.db", character.only = TRUE)
library(presto)
library(msigdbr)
library(fgsea)
library(tibble)
library(enrichR)
require(graphics)
library(fasano.franceschini.test)
library(cowplot)
library(ggsignif)
library(ggpubr)

# Set variables

R_object_folder <- "/R_objects/"

myColorsheat = c("#E31A1C", "#001CED", "#FB9A99", "#FDBF6F", "#FF7F00", "#1F78B4", "#C24F4F", "#FFFF00", "#B2DF8A",
                 "#838383", "#6A3D9A", "#CAB2D6", "#A6CEE3", "#33A02C", "#F652A0", "#B15928", "#DEF5B0")

myColors = c("#E31A1C", "#001CED", "#FB9A99", "#FDBF6F", "#FF7F00", "#1F78B4", "#C24F4F", "#FFFF00", "#B2DF8A",
             "#46FF3E", "#838383", "#6A3D9A", "#CAB2D6", "#99BA00", "#A6CEE3", "#C9C9C9", "#33A02C", "#F652A0", "#B15928", "#DEF5B0")

myColorsBC = c("#E31A1C", "#A0DEDC", "#FF7F00", "#FDBF6F", "#FFFF00", "#6A3D9A", "#CAB2D6", "#A6CEE3", 
               "#C9C9C9", "#B15928", "#DEF5B0")
memory.limit(9999999999)

# Load Seurat object

samples[["BCC5613_fresh"]] <- "/filtered_matrices/33_STO_actually_BCC5613_fresh/"
samples[["BCC5613_frozen"]] <- "/filtered_matrices/BCC5613_frozen/"
samples[["BCC5679_fresh"]] <- "/filtered_matrices/BCC5679_fresh/"
samples[["BCC5679_frozen"]] <- "/filtered_matrices/BCC5679_frozen/"
samples[["BCC5680_fresh"]] <- "/filtered_matrices/BCC5680_fresh/"
samples[["BCC5680_frozen"]] <- "/filtered_matrices/BCC5680_frozen/"
samples[["Melanoma_FNA_fresh"]] <- "/filtered_matrices/Melanoma_FNA_fresh/"
samples[["Melanoma_FNA_frozen"]] <- "/filtered_matrices/Melanoma_FNA_frozen/"
samples[["p83_CRC_fresh"]] <- "/filtered_matrices/p83_CRC_fresh/"
samples[["p83_CRC_frozen"]] <- "/filtered_matrices/p83_CRC_frozen/"
samples[["p86_CRC_fresh"]] <- "/filtered_matrices/p86_CRC_fresh/"
samples[["p86_CRC_frozen"]] <- "/filtered_matrices/p86_CRC_frozen/"

merged_project <- "all_frozen"
regev_cellcycle_path <- "regev_lab_cell_cycle_genes.txt"

cell_type_merging <- "cell_type_merging.txt"

### The reference dataset used for the celltyping ###

BP_ENCODE <- BlueprintEncodeData()
View(BP_ENCODE)
### Create and merge SingleR objects ###

for (s in ls(samples)) {
  X10 <- Read10X(data.dir = samples[[s]])
  seurat <- CreateSeuratObject(counts = X10, project = s, min.cells = 2, min.features = 100)
  seurat$orig.ident <- s
  singleR <- SingleR(test = X10, ref = BP_ENCODE, labels = BP_ENCODE$label.fine)
  
  celltypes <- data.frame("barcode" = rownames(singleR), "cell_type" = singleR$labels)
  rownames(celltypes) <- celltypes$barcode
  seurat@meta.data <- merge(x = seurat@meta.data, y = celltypes, by = 0, all.x = TRUE)
  rownames(seurat@meta.data) <- seurat@meta.data$barcode
  seurat@meta.data <- subset(seurat@meta.data, select=-c(Row.names))
  assign(s, seurat)
}

ct <- read_delim(cell_type_merging, "\t", escape_double = FALSE, trim_ws = TRUE)

################################# All SAMPLES ####################################
seurat_merged <- merge(eval(parse(text = keys(samples)[1])), 
                       y = c(eval(parse(text = keys(samples)[2])), eval(parse(text = keys(samples)[3])),
                             eval(parse(text = keys(samples)[4])), eval(parse(text = keys(samples)[5])),
                             eval(parse(text = keys(samples)[6])), eval(parse(text = keys(samples)[7])),
                             eval(parse(text = keys(samples)[8])), eval(parse(text = keys(samples)[9])),
                             eval(parse(text = keys(samples)[10])), eval(parse(text = keys(samples)[11])),
                             eval(parse(text = keys(samples)[12]))
                       ),
                       add.cell.ids = c(keys(samples)),
                       project = merged_project)

################################################################################
seurat_merged@meta.data$barcode <- rownames(seurat_merged@meta.data)

seurat_merged@meta.data <- merge(seurat_merged@meta.data, ct, by.x = "cell_type", by.y = "cell_type", all.x=TRUE)
rownames(seurat_merged@meta.data) <- seurat_merged@meta.data$barcode
seurat_merged@meta.data <- seurat_merged@meta.data[order(seurat_merged@meta.data$barcode),]

############################ Labelling samples #################################

seurat_merged@meta.data$sample <- seurat_merged@meta.data$orig.ident
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "BCC5613_fresh"] <- "BCC5613"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "BCC5613_frozen"] <- "BCC5613"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "BCC5679_fresh"] <- "BCC5679"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "BCC5679_frozen"] <- "BCC5679"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "BCC5680_fresh"] <- "BCC5680"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "BCC5680_frozen"] <- "BCC5680"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "Melanoma_FNA_fresh"] <- "Melanoma_FNA"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "Melanoma_FNA_frozen"] <- "Melanoma_FNA"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "p83_CRC_fresh"] <- "p83_CRC"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "p83_CRC_frozen"] <- "p83_CRC"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "p86_CRC_fresh"] <- "p86_CRC"
seurat_merged@meta.data$sample[seurat_merged@meta.data$sample == "p86_CRC_frozen"] <- "p86_CRC"

seurat_merged@meta.data$treatment <- seurat_merged@meta.data$orig.ident
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "BCC5613_fresh"] <- "fresh"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "BCC5613_frozen"] <- "frozen"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "BCC5679_fresh"] <- "fresh"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "BCC5679_frozen"] <- "frozen"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "BCC5680_fresh"] <- "fresh"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "BCC5680_frozen"] <- "frozen"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "Melanoma_FNA_fresh"] <- "fresh"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "Melanoma_FNA_frozen"] <- "frozen"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "p83_CRC_fresh"] <- "fresh"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "p83_CRC_frozen"] <- "frozen"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "p86_CRC_fresh"] <- "fresh"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "p86_CRC_frozen"] <- "frozen"

seurat_merged@meta.data$disease <- seurat_merged@meta.data$orig.ident
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "BCC5613_fresh"] <- "BCC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "BCC5613_frozen"] <- "BCC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "BCC5679_fresh"] <- "BCC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "BCC5679_frozen"] <- "BCC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "BCC5680_fresh"] <- "BCC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "BCC5680_frozen"] <- "BCC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "Melanoma_FNA_fresh"] <- "MEL"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "Melanoma_FNA_frozen"] <- "MEL"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "p83_CRC_fresh"] <- "CRC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "p83_CRC_frozen"] <- "CRC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "p86_CRC_fresh"] <- "CRC"
seurat_merged@meta.data$disease[seurat_merged@meta.data$disease == "p86_CRC_frozen"] <- "CRC"

cc.genes <- readLines(regev_cellcycle_path)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")

seurat_merged <- subset(seurat_merged, subset = nFeature_RNA > 99 & nCount_RNA < 180000,)
seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(seurat_merged)
seurat_merged <- ScaleData(seurat_merged, features = all.genes)
seurat_merged <- CellCycleScoring(seurat_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(object = seurat_merged))
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:20)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:30)

seurat_merged@meta.data$new_type <- seurat_merged@meta.data$merged_type
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 26 | seurat_merged@meta.data$seurat_clusters == 12 |
                                seurat_merged@meta.data$seurat_clusters == 5 | seurat_merged@meta.data$seurat_clusters == 22 | 
                                seurat_merged@meta.data$seurat_clusters == 10] <- "BCC"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 2 | seurat_merged@meta.data$seurat_clusters == 13 |
                                seurat_merged@meta.data$seurat_clusters == 21 | seurat_merged@meta.data$seurat_clusters == 25] <- "Vasc endothel"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 18] <- "Lymph endothel"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 20] <- "Mast cells"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 6] <- "CAFs"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 17] <- "Epithelial cells"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 3 | seurat_merged@meta.data$seurat_clusters == 23] <- "Fibroblasts"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 24] <- "Schwann-cells"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$cell_type == "Neutrophils"] <- "Neutrophils"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 4 | seurat_merged@meta.data$seurat_clusters == 8 |
                                seurat_merged@meta.data$seurat_clusters == 0 | seurat_merged@meta.data$seurat_clusters == 11] <- "MEL"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$seurat_clusters == 7 | seurat_merged@meta.data$seurat_clusters == 15] <- "CRC"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "Keratinocytes"] <- "Epithelial cells"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "Macrophages"] <- "Mono-Macro-DC"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "Monocytes"] <- "Mono-Macro-DC"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "DC"] <- "Mono-Macro-DC"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "Endothelial cells"] <- "Vasc endothel"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "Adipocytes"] <- "Other"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "Melanocytes" & seurat_merged@meta.data$disease == "MEL"] <- "MEL"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "MEL" & seurat_merged@meta.data$disease != "MEL"] <- "Melanocytes"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "CRC" & seurat_merged@meta.data$disease != "CRC"] <- "BCC"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "BCC" & seurat_merged@meta.data$disease == "CRC"] <- "CRC"
seurat_merged@meta.data$new_type[seurat_merged@meta.data$new_type == "BCC" & seurat_merged@meta.data$disease == "MEL"] <- "MEL"
###################################
###### Cell type plot #########

nobc_new <- ggplot() + 
  geom_bar(data = seurat_merged@meta.data, aes(x = treatment, fill = new_type), position = "fill") +  
  theme_classic() + 
  scale_fill_manual(name = "Cell types" ,values = myColors) + 
  scale_y_continuous(labels=percent) +
  xlab("Sample") + 
  ylab("Percentage") +
  labs(title = "Cell-type composition") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(.~sample, scales = "free_x", ncol = 6)
nobc_new
####################
VlnPlot(seurat_merged, "percent.mito", group.by = "orig.ident", pt.size = 0)

markers <- c("MS4A1", "YAP1", "ACTA2", "CD3E", "CD4", "CD8A", "EPCAM", "CA6", "COL1A1",
             "FLT4", "CPA3", "PMEL", "CD163", "CD55", "KLRD1", "CD38", "MPZ", "IL2RA", "CD93")
DotPlot(seurat_merged, features = markers, group.by = "new_type", split.by = "treatment", cols = c("tomato", "cyan")) + RotatedAxis()

################ Heatmap ###################

library(DoMultiBarHeatmap)
cols.use <- list(new_type=myColors)
all_heat <- subset(seurat_merged, subset = new_type != "Mast cells" & new_type != "Neutrophils" & new_type != "Other")
all_heat@meta.data <- subset(all_heat@meta.data, select=-c(barcode))

DoMultiBarHeatmap(object = all_heat, features = merged_df$genes[merged_df$notna>7 & merged_df$max>1],
                  group.by="new_type", additional.group.by=c("treatment"), group.bar = T,
                  hjust = T)

##############################################
################# GSEA ##################

m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
bccmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                         subset.ident = "BCC.BCC", group.by = "treatment")

cafmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                         subset.ident = "BCC.CAFs", group.by = "treatment")
cd4marker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                         subset.ident = "BCC.CD4+ T-cells", group.by = "treatment")
cd8marker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                         subset.ident = "BCC.CD8+ T-cells", group.by = "treatment")
epimarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                         subset.ident = "BCC.Epithelial cells", group.by = "treatment")
fibromarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                           subset.ident = "BCC.Fibroblasts", group.by = "treatment")
melamarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                          subset.ident = "BCC.Melanocytes", group.by = "treatment")
monomarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                          subset.ident = "BCC.Mono-Macro-DC", group.by = "treatment")
vascmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                          subset.ident = "BCC.Vasc endothel", group.by = "treatment")
tregmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                          subset.ident = "BCC.Tregs", group.by = "treatment")

crcbcellmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                              subset.ident = "CRC.B-cells", group.by = "treatment")
crccrcmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                            subset.ident = "CRC.CRC", group.by = "treatment")
crccafmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                            subset.ident = "CRC.CAFs", group.by = "treatment")
crccd4marker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                            subset.ident = "CRC.CD4+ T-cells", group.by = "treatment")
crccd8marker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                            subset.ident = "CRC.CD8+ T-cells", group.by = "treatment")
crcfibromarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                              subset.ident = "CRC.Fibroblasts", group.by = "treatment")
crcmonomarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                             subset.ident = "CRC.Mono-Macro-DC", group.by = "treatment")
crcvascmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                             subset.ident = "CRC.Vasc endothel", group.by = "treatment")
crctregmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                             subset.ident = "CRC.Tregs", group.by = "treatment")

melcd8marker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                            subset.ident = "MEL.CD8+ T-cells", group.by = "treatment")
melmelmarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                            subset.ident = "MEL.MEL", group.by = "treatment")
melmonomarker <- FindMarkers(seurat_obj, ident.1 = "fresh", ident.2 = "frozen",
                             subset.ident = "MEL.Mono-Macro-DC", group.by = "treatment")
"bccmarker"
dflist <- c("cafmarker", "cd4marker", "cd8marker", "epimarker", "fibromarker",
            "melamarker", "monomarker", "vascmarker", "tregmarker", 
            "crcbcellmarker", "crccrcmarker", "crccafmarker", "crccd4marker",
            "crccd8marker", "crcfibromarker", "crcmonomarker","crcvascmarker",
            "crctregmarker", "melcd8marker", "melmelmarker", "melmonomarker")

merged_df <- get("bccmarker")[get("bccmarker")$p_val_adj < 0.001,]
merged_df$genes <- rownames(merged_df)
merged_df[["bccmarker"]] <- merged_df$avg_log2FC
merged_df<- merged_df[, c("genes", "bccmarker")]


bccmarker$gene <- rownames(bccmarker)
cluster.genes <- bccmarker %>%
  #  dplyr::filter(cluster == "3") %>%   # Not necessary because there is only one cluster compared to another one in this dataframe (there is no cluster column)
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC)

ranks<- deframe(cluster.genes)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy[["bccmarker"]] <- fgseaResTidy$NES
fgseaResTidy <- fgseaResTidy[, c("pathway", "bccmarker")]
merged_df <- fgseaResTidy

for (datafra in dflist) {
  tmp_df <- get(datafra)
  tmp_df$gene <- rownames(tmp_df)
  cluster.genes <- tmp_df %>%
    #  dplyr::filter(cluster == "3") %>%   # Not necessary because there is only one cluster compared to another one in this dataframe (there is no cluster column)
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC)
  
  ranks<- deframe(cluster.genes)
  
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy[[datafra]] <- fgseaResTidy$NES
  fgseaResTidy <- fgseaResTidy[, c("pathway", datafra)]
  merged_df <- merge(merged_df, fgseaResTidy, by="pathway", all=TRUE)
}
merged_df$max <- apply(subset(merged_df, select=-c(genes)), 1, max, na.rm=TRUE)
merged_df$min <- apply(subset(merged_df, select=-c(genes)), 1, min, na.rm=TRUE)


merged_df$notna <- apply(merged_df, 1, function(x) ncol(merged_df)-sum(is.na(x)))
table(merged_df$notna)
merged_df$genes[merged_df$notna>9 & (merged_df$min<(-2) | merged_df$max>2)]
merged_df$genes[merged_df$notna>9]
#########################################################
################ Heat-shock response ###############
hsf1 <- as.data.frame(t(merged_df[merged_df$pathway=="REACTOME_HSF1_ACTIVATION",]))
colnames(hsf1) <- "REACTOME_HSF1_ACTIVATION"
hsf1 <- subset(hsf1, REACTOME_HSF1_ACTIVATION!="REACTOME_HSF1_ACTIVATION")

hsf1$disease <- c("BCC", "BCC", "BCC", "BCC", "BCC", "BCC",
                  "BCC", "BCC", "BCC", "BCC", "CRC", "CRC","CRC", "CRC",
                  "CRC", "CRC", "CRC", "CRC", "CRC", "MEL", "MEL", "MEL")
hsf1$celltype <- c("BCC", "CAFs", "CD4+ T-cells", "CD8+ T-cells",
                   "Epithelial cells", "Fibroblasts",
                   "Melanocytes", "Mono-Macro-DC",
                   "Vasc endothel", "Tregs", "B-cells", "CRC", "CAFs",
                   "CD4+ T-cells", "CD8+ T-cells", "Fibroblasts",
                   "Mono-Macro-DC", "Vasc endothel",
                   "Tregs", "CD8+ T-cells", "MEL", "Mono-Macro-DC")

hsf1$disease <- factor(hsf1$disease, levels = c("BCC", "CRC", "MEL"))
myColorsbc = c("#001CED", "#FFFF00", "#FDBF6F", "#FF7F00", "#1F78B4", 
               "#C24F4F", "#FB9A99", "#838383", "#6A3D9A", "#CAB2D6",
               "#B15928", "#DEF5B0")

p <- ggplot(hsf1[hsf1$celltype!="Plasma cells" &
                   hsf1$celltype != "B-cells",], aes(x=celltype, y=as.numeric(REACTOME_HSF1_ACTIVATION), fill = celltype)) + 
  geom_bar(stat='identity') + facet_grid(.~disease, scales = "free_x", space = "free") + theme_classic() +
  scale_fill_manual(name = "Cell types" ,values = myColorsbc) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
#################################################
######## QC plot for supplementary #########
qcdf <- read.csv("QC2.txt", sep = "\t")
cellbar <-ggplot(qcdf, aes(x=sample, y=cells, fill = treatment)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cellbar

readbar <-ggplot(qcdf, aes(x=sample, y=reads, fill = treatment)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
readbar

genebar <-ggplot(qcdf, aes(x=sample, y=genes, fill = treatment)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
genebar
##############################################

############## Cell-type markers supplementary #############
schwannmarkers <- c("SOX10", "MPZ", "NGFR", "S100B", "NCAM1")
DotPlot(seurat_obj, features = schwannmarkers, group.by = "new_type") + RotatedAxis()

vascendomarkers <- c("ACE", "CD93", "PECAM1", "CD34")
DotPlot(seurat_obj, features = vascendomarkers, group.by = "new_type") + RotatedAxis()

lymphendomarkers <- c("FLT4", "PDPN", "PROX1", "LYVE1")
DotPlot(seurat_obj, features = lymphendomarkers, group.by = "new_type") + RotatedAxis()

cafmarkers <- c("SMN1", "S100A4", "FAP", "VIM", "PDGFRA")
DotPlot(seurat_obj, features = cafmarkers, group.by = "new_type") + RotatedAxis()

cafmarkers <- c("ACTA2", "LBH", "TAGLN", "MCAM")
DotPlot(seurat_obj, features = cafmarkers, group.by = "new_type") + RotatedAxis()

cafDEmarkers <- FindMarkers(seurat_obj, ident.1 = "CAFs", group.by = "new_type")
View(cafDEmarkers)

mastmarkers <- c("TPSAB1", "CPA3", "CTSG", "KIT")
DotPlot(seurat_obj, features = mastmarkers, group.by = "new_type") + RotatedAxis()

melanomamarkers <- c("S100A1", "PAX3", "PSAT1", "ABCB5")
melanomamarkers <- c("AXL", "CYR61", "DHRS3", "EGR1", "HMGA2", "MOXD1", "PSAT1", "RAB3B", "PCOLCE", "IER2", "DUSP1")
DotPlot(seurat_obj, features = melanomamarkers, group.by = "new_type") + RotatedAxis()

crcmarkers <- c("CEACAM5", "FABP1", "TFF3", "MKI67")
DotPlot(seurat_obj, features = crcmarkers, group.by = "new_type") + RotatedAxis()

bccmarkers <- c("DMKN", "IVL", "IRF6", "SFN", "KRT5")
DotPlot(seurat_obj, features = bccmarkers, group.by = "new_type") + RotatedAxis()
bccDEmarkers <- FindMarkers(seurat_obj, ident.1 = "BCC", ident.2 = "Epithelial cells", group.by = "new_type")

allmarkers <- c("MPZ", "NGFR", "S100B", "NCAM1", "ACE", "CD93", "PECAM1", "CD34", "FLT4", "PDPN", "PROX1", "LYVE1", 
                "ACTA2", "LBH", "TAGLN", "MCAM", "TPSAB1", "CPA3", "CTSG", "KIT",
                "S100A1", "PAX3", "PSAT1", "ABCB5", "CEACAM5", "FABP1", "TFF3", "MKI67", "DMKN", "IVL", "IRF6", "SFN")
DotPlot(seurat_obj, features = allmarkers, group.by = "new_type") + RotatedAxis()

############################################################
################# Breast cancer samples ####################
rm(samples)
### Create a BC seurat object ###

samples[["UHB129_BC_frozen"]] <- "/filtered_matrices/UHB129_BC_frozen/"
samples[["UHB150_BC_frozen"]] <- "/filtered_matrices/UHB150_BC_frozen/"
samples[["UHB173_BC_frozen"]] <- "/filtered_matrices/UHB173_BC_frozen/"
samples[["UHB182_BC_frozen"]] <- "/filtered_matrices/UHB182_BC_frozen/"
samples[["UHB194_BC_frozen"]] <- "/filtered_matrices/UHB194_BC_frozen/"

### Create and merge SingleR objects ###

for (s in ls(samples)) {
  X10 <- Read10X(data.dir = samples[[s]])
  seurat <- CreateSeuratObject(counts = X10, project = s, min.cells = 2, min.features = 100)
  seurat$orig.ident <- s
  singleR <- SingleR(test = X10, ref = BP_ENCODE, labels = BP_ENCODE$label.fine)
  
  celltypes <- data.frame("barcode" = rownames(singleR), "cell_type" = singleR$labels)
  rownames(celltypes) <- celltypes$barcode
  seurat@meta.data <- merge(x = seurat@meta.data, y = celltypes, by = 0, all.x = TRUE)
  rownames(seurat@meta.data) <- seurat@meta.data$barcode
  seurat@meta.data <- subset(seurat@meta.data, select=-c(Row.names))
  assign(s, seurat)
}

ct <- read_delim(cell_type_merging, "\t", escape_double = FALSE, trim_ws = TRUE)

################################# All SAMPLES ####################################
BC <- merge(eval(parse(text = keys(samples)[1])), 
            y = c(eval(parse(text = keys(samples)[2])), eval(parse(text = keys(samples)[3])),
                  eval(parse(text = keys(samples)[4])), eval(parse(text = keys(samples)[5])),
                  ),
            add.cell.ids = c(keys(samples)),
            project = merged_project)

################################################################################
########### Setting labels #################
BC@meta.data$barcode <- rownames(BC@meta.data)

BC@meta.data <- merge(BC@meta.data, ct, by.x = "cell_type", by.y = "cell_type", all.x=TRUE)
rownames(BC@meta.data) <- BC@meta.data$barcode
BC@meta.data <- BC@meta.data[order(BC@meta.data$barcode),]

BC@meta.data$new_type <- BC@meta.data$merged_type
BC@meta.data$new_type[BC@meta.data$seurat_clusters == 2 | BC@meta.data$seurat_clusters == 3 |
                        BC@meta.data$seurat_clusters == 14 | BC@meta.data$seurat_clusters == 13 |  
                        BC@meta.data$seurat_clusters == 5 | BC@meta.data$seurat_clusters == 8 |
                        BC@meta.data$seurat_clusters == 12] <- "BC"
BC@meta.data$new_type[BC@meta.data$new_type == "Plasma cells"] <- "B-cells"
BC@meta.data$new_type[BC@meta.data$new_type == "Keratinocytes"] <- "Epithelial cells"
BC@meta.data$new_type[BC@meta.data$new_type == "Other"] <- "Epithelial cells"
BC@meta.data$new_type[BC@meta.data$new_type == "Adipocytes"] <- "Vasc endothel"
BC@meta.data$new_type[BC@meta.data$new_type == "Endothelial cells"] <- "Vasc endothel"
BC@meta.data$new_type[BC@meta.data$new_type == "Myocytes"] <- "Fibroblasts"
BC@meta.data$new_type[BC@meta.data$new_type == "Macrophages"] <- "Mono-Macro-DC"
BC@meta.data$new_type[BC@meta.data$new_type == "Monocytes"] <- "Mono-Macro-DC"
BC@meta.data$new_type[BC@meta.data$new_type == "DC"] <- "Mono-Macro-DC"
BC@meta.data$new_type[(BC@meta.data$seurat_clusters == 11 | BC@meta.data$seurat_clusters == 1 |
                         BC@meta.data$seurat_clusters == 10)] <- "Fibroblasts"
BC@meta.data$new_type[BC@meta.data$seurat_clusters == 11] <- "BC"

BC@meta.data$sample <- BC@meta.data$orig.ident
BC@meta.data$sample[BC@meta.data$sample == "UHB129_BC_frozen"] <- "UHB129_BC"
BC@meta.data$sample[BC@meta.data$sample == "UHB150_BC_frozen"] <- "UHB150_BC"
BC@meta.data$sample[BC@meta.data$sample == "UHB173_BC_frozen"] <- "UHB173_BC"
BC@meta.data$sample[BC@meta.data$sample == "UHB182_BC_frozen"] <- "UHB182_BC"
BC@meta.data$sample[BC@meta.data$sample == "UHB194_BC_frozen"] <- "UHB194_BC"
BC@meta.data$treatment <- BC@meta.data$orig.ident
BC@meta.data$treatment[BC@meta.data$treatment == "UHB129_BC_frozen"] <- "frozen"
BC@meta.data$treatment[BC@meta.data$treatment == "UHB150_BC_frozen"] <- "frozen"
BC@meta.data$treatment[BC@meta.data$treatment == "UHB173_BC_frozen"] <- "frozen"
BC@meta.data$treatment[BC@meta.data$treatment == "UHB182_BC_frozen"] <- "frozen"
BC@meta.data$treatment[BC@meta.data$treatment == "UHB194_BC_frozen"] <- "frozen"
BC@meta.data$disease <- BC@meta.data$orig.ident
BC@meta.data$disease[BC@meta.data$disease == "UHB129_BC_frozen"] <- "BC"
BC@meta.data$disease[BC@meta.data$disease == "UHB150_BC_frozen"] <- "BC"
BC@meta.data$disease[BC@meta.data$disease == "UHB173_BC_frozen"] <- "BC"
BC@meta.data$disease[BC@meta.data$disease == "UHB182_BC_frozen"] <- "BC"
BC@meta.data$disease[BC@meta.data$disease == "UHB194_BC_frozen"] <- "BC"
BC[["percent.mt"]] <- PercentageFeatureSet(BC, pattern = "^MT-")
BC <- subset(BC, subset=percent.mt<50)
FeaturePlot(BC, c("percent.mt"))

BC <- NormalizeData(BC, normalization.method = "LogNormalize", scale.factor = 10000)
BC <- FindVariableFeatures(BC, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(BC)
BC <- ScaleData(BC, features = all.genes)
BC <- CellCycleScoring(BC, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
BC <- RunPCA(BC, features = VariableFeatures(object = BC))
BC <- FindNeighbors(BC, dims = 1:20)
BC <- FindClusters(BC, resolution = 0.5)
BC <- RunUMAP(BC, dims = 1:30)

BC@meta.data$tissue[BC@meta.data$new_type == "BC"] <- "tumor"
BC@meta.data$tissue[BC@meta.data$new_type != "BC"] <- "stroma"
BC@meta.data$sample[BC@meta.data$sample == "UHB129_BC"] <- "UHB129"
BC@meta.data$sample[BC@meta.data$sample == "UHB150_BC"] <- "UHB150"
BC@meta.data$sample[BC@meta.data$sample == "UHB173_BC"] <- "UHB173"
BC@meta.data$sample[BC@meta.data$sample == "UHB182_BC"] <- "UHB182"
BC@meta.data$sample[BC@meta.data$sample == "UHB194_BC"] <- "UHB194"
##############################################################
############ Expression of markers ####################

table(BC@meta.data$sample)
dat = data.frame(sample=BC@meta.data$sample,
                 tissue=BC@meta.data$tissue,
                 new_type=BC@meta.data$new_type,
                 readcount=BC@meta.data$nCount_RNA,
                 EPCAM=t(as.array(BC@assays$RNA[c("EPCAM")])), 
                 CD3E=t(as.array(BC@assays$RNA[c("CD3E")])), 
                 CD4=t(as.array(BC@assays$RNA[c("CD4")])), 
                 CD8A=t(as.array(BC@assays$RNA[c("CD8A")])), 
                 ITGAX=t(as.array(BC@assays$RNA[c("ITGAX")])), 
                 CD14=t(as.array(BC@assays$RNA[c("CD14")])), 
                 CD19=t(as.array(BC@assays$RNA[c("CD19")])), 
                 MS4A1=t(as.array(BC@assays$RNA[c("MS4A1")])), 
                 PECAM1=t(as.array(BC@assays$RNA[c("PECAM1")])), 
                 NCAM1=t(as.array(BC@assays$RNA[c("NCAM1")])), 
                 CD68=t(as.array(BC@assays$RNA[c("CD68")])), 
                 KRT5=t(as.array(BC@assays$RNA[c("KRT5")])), 
                 KRT8=t(as.array(BC@assays$RNA[c("KRT8")])), 
                 KRT18=t(as.array(BC@assays$RNA[c("KRT18")])), 
                 KRT14=t(as.array(BC@assays$RNA[c("KRT14")])), 
                 ESR1=t(as.array(BC@assays$RNA[c("ESR1")])),
                 FOXP3=t(as.array(BC@assays$RNA[c("FOXP3")])), 
                 GATA3=t(as.array(BC@assays$RNA[c("GATA3")])), 
                 ERBB2=t(as.array(BC@assays$RNA[c("ERBB2")])), 
                 MKI67=t(as.array(BC@assays$RNA[c("MKI67")])), 
                 PGR=t(as.array(BC@assays$RNA[c("PGR")])), 
                 SMN1=t(as.array(BC@assays$RNA[c("SMN1")])),
                 SOX10=t(as.array(BC@assays$RNA[c("SOX10")])), 
                 VIM=t(as.array(BC@assays$RNA[c("VIM")])))

dat$EPCAMbin <- as.factor(ifelse(dat$EPCAM == 0, 0, 1))
dat$CD3Ebin <- as.factor(ifelse(dat$CD3E == 0, 0, 1))
dat$CD4bin <- as.factor(ifelse(dat$CD4 == 0, 0, 1))
dat$CD8Abin <- as.factor(ifelse(dat$CD8A == 0, 0, 1))
dat$ITGAXbin <- as.factor(ifelse(dat$ITGAX == 0, 0, 1))
dat$CD14bin <- as.factor(ifelse(dat$CD14 == 0, 0, 1))
dat$CD19bin <- as.factor(ifelse(dat$CD19 == 0, 0, 1))
dat$MS4A1bin <- as.factor(ifelse(dat$MS4A1 == 0, 0, 1))
dat$PECAM1bin <- as.factor(ifelse(dat$PECAM1 == 0, 0, 1))
dat$NCAM1bin <- as.factor(ifelse(dat$NCAM1 == 0, 0, 1))
dat$CD68bin <- as.factor(ifelse(dat$CD68 == 0, 0, 1))
dat$KRT5bin <- as.factor(ifelse(dat$KRT5 == 0, 0, 1))
dat$KRT8bin <- as.factor(ifelse(dat$KRT8 == 0, 0, 1))
dat$KRT18bin <- as.factor(ifelse(dat$KRT18 == 0, 0, 1))
dat$KRT14bin <- as.factor(ifelse(dat$KRT14 == 0, 0, 1))
dat$ESR1bin <- as.factor(ifelse(dat$ESR1 == 0, 0, 1))
dat$FOXP3bin <- as.factor(ifelse(dat$FOXP3 == 0, 0, 1))
dat$GATA3bin <- as.factor(ifelse(dat$GATA3 == 0, 0, 1))
dat$ERBB2bin <- as.factor(ifelse(dat$ERBB2 == 0, 0, 1))
dat$MKI67bin <- as.factor(ifelse(dat$MKI67 == 0, 0, 1))
dat$PGRbin <- as.factor(ifelse(dat$PGR == 0, 0, 1))
dat$SMN1bin <- as.factor(ifelse(dat$SMN1 == 0, 0, 1))
dat$SOX10bin <- as.factor(ifelse(dat$SOX10 == 0, 0, 1))
dat$VIMbin <- as.factor(ifelse(dat$VIM == 0, 0, 1))

BCquant <- c("EPCAM", "CD3E", "CD4", "CD8A", "ITGAX", "CD14", "CD19", "MS4A1", "PECAM1", "NCAM1", "CD68", "KRT5",
             "KRT8", "KRT18", "KRT14", "ESR1", "FOXP3", "GATA3", "ERBB2", "MKI67", "PGR", "SMN1", "SOX10", "VIM")
DotPlot(BC, features = BCquant, group.by = "sample", split.by = "tissue") #+ RotatedAxis()
#View(wdata)
wdata <- dat[dat$readcount > 99 & dat$sample %in% c("UHB129", "UHB150"),]
genelist <-c("EPCAMbin", "CD3Ebin", "CD4bin", "CD8Abin", "ITGAXbin", "CD14bin", "MS4A1bin",
             "PECAM1bin", "NCAM1bin", "CD68bin", "KRT5bin", "KRT8bin", "KRT14bin", "ESR1bin", "FOXP3bin",
             "GATA3bin", "ERBB2bin", "MKI67bin", "PGRbin", "SMN1bin", "SOX10bin", "VIMbin")
res <- dcast(wdata, sample + tissue ~ EPCAMbin)
res <- res[, c("sample", "tissue")]
for (gene in genelist) {
  tempdf <- dcast(wdata, sample + tissue ~ get(gene))
  tempdf$frac <- tempdf$`1` / (tempdf$`0` + tempdf$`1`)
  res[gene] <- tempdf$frac
}
res$method <- "scRNA"
res[nrow(res) + 1,] = list("UHB129", "stroma", 0.00082607, 0.129882266, 0.08506586, 0.08197431, 0.05857182, 0.34693878,
                           0.08376521, 0.11940293, 0.01376796, 0.56491168, 0.27185435, 0.11323231, 0.0051061, 0.07706714,
                           0.01237351, 0.21149714, 0.00909134, 0.00286757, 0.00758663, 0.06245158, 0.18648968, 0.1136918, "antibody")
res[nrow(res) + 1,] = list("UHB129", "tumor", 0.11694843, 0.01537269, 0.04272203, 0.06601047, 0.02279095, 0.3458445,
                           0.003194, 0.31250134, 0.01048604, 0.4103249, 0.12379886, 0.74862419, 0.09305929, 0.99348083,
                           0.00294435, 0.86476265, 0.24124531, 0.04737833, 0.36038616, 0.2147687, 0.02768747, 0.60313267, "antibody")
res[nrow(res) + 1,] = list("UHB150", "stroma", 0.01102614, 0.05048063, 0.03507548, 0.03215363, 0.12872626, 0.33554218,
                           0.01039912, 0.20922655, 0.01729599, 0.2550682, 0.00054963, 0.03846013, 0.00695765, 0.00077055,
                           0.01121854, 0.05466813, 0.00684637, 0.00834103, 0.03078876, 0.25245762, 0.06234834, 0.58987199, "antibody")
res[nrow(res) + 1,] = list("UHB150", "tumor", 0.84960335, 0.00178258, 0.00121982, 0.00505735, 0.00075614, 0.17361603,
                           0.02040311, 0.02045189, 0.00033321, 0.01511076, 0.07010333, 0.5118446, 0.04685308, 0.01576427,
                           0.00044443, 0.63383743, 0.00085184, 0.30317532, 0.02268218, 0.06177749, 0.84609512, 0.02146652, "antibody")


long <- melt(setDT(res), id.vars = c("sample","tissue", "method"), variable.name = "gene", value.name = "frac")
View(long)
long$x_axis <- paste(long$sample, long$tissue, long$method, sep="_")
long$x_axis2 <- paste(long$sample, long$tissue, sep="_") 

ggplot(long, aes(x = gene, y = frac * 100, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("% of cells expressing") +
  facet_grid(sample~tissue)

ggplot(long, aes(x = gene, y = frac * 100, fill = x_axis)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("% of cells expressing") + theme_classic() +
  scale_fill_manual(values = myColorsbox)+
  facet_wrap(.~sample) + RotatedAxis()

ggplot(long, aes(x = x_axis, y = frac * 100, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("% of cells expressing") + theme_classic() +
  scale_fill_manual(values = myColorsbox)+
  facet_wrap(.~gene, ncol = 3) + RotatedAxis()

ggplot(long, aes(x = x_axis2, y = frac * 100, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("% of cells expressing") + theme_classic() +
  scale_fill_manual(values = myColorsbox)+
  facet_wrap(.~gene, ncol = 4) + RotatedAxis()

genes <- c("EPCAMbin", "ESR1bin", "PGRbin", "GATA3bin", "MKI67bin", "PECAM1bin", "CD3Ebin", "CD4bin", "CD8Abin", "CD68bin")

cor(long$frac[long$method=="antibody"], long$frac[long$method=="scRNA"])
cor.test(long$frac[long$method=="antibody"], long$frac[long$method=="scRNA"], method = "spearman")
summary(cor(long$frac[long$method=="antibody"], long$frac[long$method=="scRNA"], method = "spearman"))
#########################################
########### Correlation plot ##########
newlong <- data.frame(ID = unique(long$x_axis), antibody = long$frac[long$method=="antibody"], scRNA = long$frac[long$method=="scRNA"])

ggscatter(newlong, "antibody", y = "scRNA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Expression% (antibody)", ylab = "Expression% (scRNA)")
tumlong <- long[long$tissue=="tumor", ]
newtum <- data.frame(ID = unique(tumlong$x_axis), antibody = tumlong$frac[tumlong$method=="antibody"], scRNA = tumlong$frac[tumlong$method=="scRNA"])
ggscatter(newtum, "antibody", y = "scRNA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Expression% (antibody)", ylab = "Expression% (scRNA)")

##### Analyzing only breast cancer cells #####
onlyBC <- subset(BC, subset=new_type == "BC")
onlyBC <- NormalizeData(onlyBC, normalization.method = "LogNormalize", scale.factor = 10000)
onlyBC <- FindVariableFeatures(onlyBC, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(onlyBC)
onlyBC <- ScaleData(onlyBC, features = all.genes)
onlyBC <- CellCycleScoring(onlyBC, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
onlyBC <- RunPCA(onlyBC, features = VariableFeatures(object = onlyBC))
onlyBC <- FindNeighbors(onlyBC, dims = 1:20)
onlyBC <- FindClusters(onlyBC, resolution = 0.5)
onlyBC <- RunUMAP(onlyBC, dims = 1:30)

DotPlot(onlyBC, features = c("ESR1", "PGR", "ERBB2"), group.by = "sample")
VlnPlot(onlyBC, c("ESR1", "PGR", "ERBB2", "EGFR"), group.by = "sample", pt.size = 0)
###########################################################################