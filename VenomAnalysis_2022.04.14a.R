##### venom gland snRNAseq analysis #####

##### updated 4/14/22

### (0) ENVIRONMENT
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(viridis)
library(grid)
library(gridExtra)
library(pheatmap)
library(reshape2)
library(SingleCellExperiment)
library(scater)
library(scran)
library(matrixStats)
library(escape)
library(dittoSeq)
library(GSEABase)
library(scico)
library(ggsci)
library(DESeq2)
library(cowplot)
library(ggcorrplot)

setwd("~/Desktop/snRNAseq_viridis/")

'%notin%' <- function(x,y)!('%in%'(x,y))

venomGenes <- read.csv("toxin_transcripts.csv", header = TRUE)
row.names(venomGenes) <- venomGenes$transcript

### (1) Data processing and clustering #####
##### Data should be gzipped
data <- Read10X(data.dir = "completeMapping/filtered_feature_bc_matrix/")
venom <- CreateSeuratObject(counts = data, project = "venom", min.cells = 3, min.features = 200)

venom[["percent.mt"]] <- PercentageFeatureSet(venom, pattern = "^nbis-")
VlnPlot(venom, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(venom, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

### Normalization
venom <- NormalizeData(venom, normalization.method = "LogNormalize", scale.factor = 10000)

### Identify highly variable features
venom <- FindVariableFeatures(venom, selection.method = "vst", nfeatures = 2000,
                              mean.function = ExpMean,
                              dispersion.function = LogVMR,
                              dispersion.cutoff = c(0.0125, 3))
top10 <- head(VariableFeatures(venom), 10)

plot1 <- VariableFeaturePlot(venom)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

### Scale data
all.genes <- rownames(venom)
venom <- ScaleData(venom, vars.to.regress = c("nCount_RNA", "percent.mt"), features = all.genes)

### Linear dimensional reduction
#### PCA for full set of computationally identified variable features
venom <- RunPCA(venom, features = VariableFeatures(venom))

print(venom[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(venom, dims = 1:2, reduction = "pca")
DimPlot(venom, reduction = "pca")
DimHeatmap(venom, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(venom, dims = 1:15, cells = 500, balanced = TRUE)

### Determine dimensionality
venom <- JackStraw(venom, num.replicate = 100)
venom <- ScoreJackStraw(venom, dims = 1:20)
JackStrawPlot(venom, dims = 1:20)
ElbowPlot(venom) ### elbow seems to occur around 6-7 for venom gland

### Cluster cells
venom <- FindNeighbors(venom, dims = 1:6)
venom <- FindClusters(venom, resolution = 0.5)
head(Idents(venom), 5)

### Non-linear dimensional reduction
venom <- RunTSNE(object = venom, dims.use = 1:6)
DimPlot(venom, reduction = "tsne")

venom <- RunUMAP(venom, dims = 1:6)
DimPlot(venom, reduction = "umap")

saveRDS(venom, file = "venom_Seurat.rds")
venom <- readRDS(file = "venom_Seurat.rds")


#### PCA for only venom genes expressed
venom2 <- RunPCA(venom, features = venomGenes$transcript)

print(venom2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(venom2, dims = 1:2, reduction = "pca")
DimPlot(venom2, reduction = "pca")
DimHeatmap(venom2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(venom2, dims = 1:15, cells = 500, balanced = TRUE)

### Determine dimensionality
venom2 <- JackStraw(venom2, num.replicate = 100)
venom2 <- ScoreJackStraw(venom2)
JackStrawPlot(venom2)
ElbowPlot(venom2) ### elbow seems to occur around 6 for venom gland

### Cluster cells
venom2 <- FindNeighbors(venom2, dims = 1:6)
venom2 <- FindClusters(venom2, resolution = 0.5)
head(Idents(venom2), 5)

### Non-linear dimensional reduction
venom2 <- RunUMAP(venom2, dims = 1:6)
DimPlot(venom2, reduction = "umap")
saveRDS(venom, file = "new.venom_Seurat.rds")
venom <- readRDS(file = "venom_Seurat.rds")

saveRDS(venom2, file = "new.venom_Seurat_venomProteins.rds")
venom2 <- readRDS("venom_Seurat_venomProteins.rds")

venom$seurat_clusters <- as.factor(as.numeric(venom$seurat_clusters))
venom2$seurat_clusters <- as.factor(as.numeric(venom2$seurat_clusters))

### (2) Identify Biomarkers and start plotting ######
layout_matrix <- rbind(c(1, 2, 3),
                       c(4, 4, 5))

venom.markers <- FindAllMarkers(venom, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- venom.markers[which(venom.markers$gene %notin% all.genes[grep("^nbis", all.genes, perl = TRUE)]),] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
top10$cluster <- as.factor(as.numeric(top10$cluster))
write.csv(top10, "newMapping.venom_top10.csv")
top10.names <- read.csv("newMapping.venom_top10.csv", header = TRUE, row.names = 1)
top10$name <- top10.names$name

venom2$orig.clust <- venom$seurat_clusters
venom$ven.clust <- venom2$seurat_clusters

umap.seurat <- UMAPPlot(venom,
         cols = c("#0C0A3E", "#7B1E7A", "#D62246", "#E27D60", "#F3C677"),
         group.by = "seurat_clusters") +
  ggtitle("A. UMAP clusters") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.1, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.margin = margin(c(-10, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

heatmap <- DoHeatmap(venom, features = top10$gene,
                    group.colors = c("#0C0A3E", "#7B1E7A", "#D62246", "#E27D60", "#F3C677"),
                    group.by = "seurat_clusters",
                    group.bar.height = 0.01,
                    size = 2.5, angle = 0, draw.lines = TRUE,
                    hjust = 0.6) + 
  scale_fill_gradient2(low = "blue", mid = "black", high = "red", name = "expression") +
  scale_y_discrete(labels = rev(top10$name)) +
  guides(colour = "none") +
  ggtitle("B. Top cluster markers") +
  theme(text = element_text(size = 4),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        title = element_text(size = 6),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        legend.position = "bottom",
        legend.margin = margin(c(-10, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

umap.seurat_venom <- UMAPPlot(venom,
         #cols = c("#D71B3B", "#E8D71E", "#16ACEA", "#4203C9", "orange"),
         cols = c("#D71B3B", "#E8D71E", "#16ACEA"),
         group.by = "ven.clust") +
  ggtitle("C. UMAP clusters") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.1, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.margin = margin(c(-10, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))


umap.venom <- UMAPPlot(venom2,
                   #cols = c("#D71B3B", "#E8D71E", "#16ACEA", "#4203C9", "orange"),
                   cols = c("#D71B3B", "#E8D71E", "#16ACEA"),
                   group.by = "seurat_clusters") +
  ggtitle("E. Venom-based clustering") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.1, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.margin = margin(c(-10, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

grid.arrange(umap.seurat, umap.seurat_venom, umap.venom, nrow = 1)


UMAPPlot(venom2,
         cols = c("#0C0A3E", "#7B1E7A", "#D62246", "#E27D60", "#F3C677"),
         group.by = "orig.clust")
UMAPPlot(venom2,
         cols = c("#D71B3B", "#E8D71E", "#16ACEA", "#4203C9"),
         group.by = "seurat_clusters") +
  ggtitle("Venom clusters") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.4, 'cm'),
        legend.margin = margin(c(-2, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))



### (2b) stacked violin
expr <- venom@assays$RNA@scale.data %>% as.matrix %>% t %>% as.data.frame
expr <- expr[, venomGenes[which(row.names(venomGenes) %in% colnames(expr)),1]]
expr$Idents <- venom$seurat_clusters
expr$Cell <- row.names(expr)

features <- venomGenes[which(row.names(venomGenes) %in% colnames(expr)),1]
features2 <- factor(venomGenes[features,2], levels = venomGenes[features,2])

expr <- reshape2::melt(expr, id.vars = c("Cell","Idents"), measure.vars = features,
                       variable.name = "Feat", value.name = "Expr")

expr$Feat <- venomGenes[as.character(expr$Feat),2]
expr$Feat <- factor(expr$Feat, levels = features2)

venomClasses <- data.frame(x = features2,
                           group = venomGenes[features, 3])
color <- viridis(12)

ggplot(expr, aes(Feat, Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, lwd = 0.1) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
  theme_bw() +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        text = element_text(size = 6),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("D. Venom gene expression in clusters") + ylab("scaled expression")

g <- ggplot(venomClasses, aes(x = x, y = 1, fill = group, label = group)) + geom_tile() +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = color) + scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        plot.background = element_blank(), 
        plot.margin = margin(0, 7, 7, 7, "pt"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

violin <- plot_grid(f, g, ncol = 1, rel_heights = c(0.8, 0.2), align = "v", axis = "lr")

## 4.7 x 4.02
grid.arrange(umap.seurat, heatmap, umap.seurat_venom,
             violin, umap.venom, layout_matrix = layout_matrix)

### (2c) Venom gene spreads ######
featureTheme <- theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        text = element_text(size = 6),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.4, 'cm'),
        legend.margin = margin(c(-2, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[4], order = TRUE) +
               labs(title = "SVSP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[5], order = TRUE) +
               labs(title = "SVSP2", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[6], order = TRUE) +
               labs(title = "SVSP3", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[7], order = TRUE) + 
               labs(title = "SVSP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[8], order = TRUE) +
               labs(title = "SVSP5", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[9], order = TRUE) +
               labs(title = "SVSP6", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[10], order = TRUE) +
               labs(title = "SVSP7", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[11], order = TRUE) +
               labs(title = "SVSP8", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[12], order = TRUE) +
               labs(title = "SVSP9", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme)


grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[15], order = TRUE) +
               labs(title = "SVMP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[16], order = TRUE) +
               labs(title = "SVMP2", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[17], order = TRUE) +
               labs(title = "SVMP3", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[18], order = TRUE) +
               labs(title = "SVMP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[19], order = TRUE) +
               labs(title = "SVMP5", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[20], order = TRUE) +
               labs(title = "SVMP6", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[21], order = TRUE) + 
               labs(title = "SVMP7", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[22], order = TRUE) + 
               labs(title = "SVMP8", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[23], order = TRUE) + 
               labs(title = "SVMP9", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[24], order = TRUE) + 
               labs(title = "SVMP10", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             nrow = 3)

y.svsp <- grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[7], order = TRUE) + 
                         labs(title = "SVSP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
                       FeaturePlot(venom2, features = venomGenes$transcript[4], order = TRUE) +
                         labs(title = "SVSP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme, 
                       FeaturePlot(venom2, features = venomGenes$transcript[11], order = TRUE) +
                         labs(title = "SVSP8", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme, 
                       FeaturePlot(venom2, features = venomGenes$transcript[10], order = TRUE) +
                         labs(title = "SVSP7", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
                       nrow = 1)

x.svmp <- grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[18], order = TRUE) +
                         labs(title = "SVMP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme, 
                       FeaturePlot(venom2, features = venomGenes$transcript[20], order = TRUE) +
                         labs(title = "SVMP6", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme, 
                       FeaturePlot(venom2, features = venomGenes$transcript[15], order = TRUE) +
                         labs(title = "SVMP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme, 
                       FeaturePlot(venom2, features = venomGenes$transcript[19], order = TRUE) +
                         labs(title = "SVMP5", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
                       nrow = 1)

layout_matrix2 <- rbind(c(1, 9, 9, 9, 9),
                        c(2, 9, 9, 9, 9),
                        c(3, 9, 9, 9, 9),
                        c(4, 5, 6, 7, 8))

svsp4 <- FeaturePlot(venom2, features = venomGenes$transcript[7], order = TRUE) + 
  labs(title = "SVSP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp4$layers[[1]]$aes_params$size <- 0.3

svsp1 <- FeaturePlot(venom2, features = venomGenes$transcript[4], order = TRUE) +
  labs(title = "SVSP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp1$layers[[1]]$aes_params$size <- 0.3

svsp8 <- FeaturePlot(venom2, features = venomGenes$transcript[11], order = TRUE) +
  labs(title = "SVSP8", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp8$layers[[1]]$aes_params$size <- 0.3

svsp7 <- FeaturePlot(venom2, features = venomGenes$transcript[10], order = TRUE) +
  labs(title = "SVSP7", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp7$layers[[1]]$aes_params$size <- 0.3

svsp4$layers[[1]]$aes_params$alpha <- 0.5
svsp1$layers[[1]]$aes_params$alpha <- 0.5
svsp8$layers[[1]]$aes_params$alpha <- 0.5
svsp7$layers[[1]]$aes_params$alpha <- 0.5


svmp4 <- FeaturePlot(venom2, features = venomGenes$transcript[18], order = TRUE) +
  labs(title = "SVMP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svmp6 <- FeaturePlot(venom2, features = venomGenes$transcript[20], order = TRUE) +
  labs(title = "SVMP6", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svmp1 <- FeaturePlot(venom2, features = venomGenes$transcript[15], order = TRUE) +
  labs(title = "SVMP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svmp5 <- FeaturePlot(venom2, features = venomGenes$transcript[19], order = TRUE) +
  labs(title = "SVMP5", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme

svmp4$layers[[1]]$aes_params$size <- 0.3
svmp6$layers[[1]]$aes_params$size <- 0.3
svmp1$layers[[1]]$aes_params$size <- 0.3
svmp5$layers[[1]]$aes_params$size <- 0.3

svmp4$layers[[1]]$aes_params$alpha <- 0.5
svmp6$layers[[1]]$aes_params$alpha <- 0.5
svmp1$layers[[1]]$aes_params$alpha <- 0.5
svmp5$layers[[1]]$aes_params$alpha <- 0.5

pla2a1 <- FeaturePlot(venom2, features = venomGenes$transcript[29], order = TRUE) + labs(title = venomGenes$toxin[29], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
pla2b1 <- FeaturePlot(venom2, features = venomGenes$transcript[26], order = TRUE) + labs(title = venomGenes$toxin[26], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
crisp1 <- FeaturePlot(venom2, features = venomGenes$transcript[34], order = TRUE) + labs(title = venomGenes$toxin[34], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
ohanin <- FeaturePlot(venom2, features = venomGenes$transcript[36], order = TRUE) + labs(title = venomGenes$toxin[36], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
exo3 <- FeaturePlot(venom2, features = venomGenes$transcript[39], order = TRUE) + labs(title = venomGenes$toxin[39], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme

pla2a1$layers[[1]]$aes_params$size <- 0.3
pla2b1$layers[[1]]$aes_params$size <- 0.3
crisp1$layers[[1]]$aes_params$size <- 0.3
ohanin$layers[[1]]$aes_params$size <- 0.3
exo3$layers[[1]]$aes_params$size <- 0.3

pla2a1$layers[[1]]$aes_params$alpha <- 0.5
pla2b1$layers[[1]]$aes_params$alpha <- 0.5
crisp1$layers[[1]]$aes_params$alpha <- 0.5
ohanin$layers[[1]]$aes_params$alpha <- 0.5
exo3$layers[[1]]$aes_params$alpha <- 0.5

umap <- UMAPPlot(venom2,
                       #cols = c("#D71B3B", "#E8D71E", "#16ACEA", "#4203C9", "orange"),
                       cols = c("#D71B3B", "#E8D71E", "#16ACEA"),
                       group.by = "seurat_clusters") +
  ggtitle("Venom-based clustering") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.1, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.margin = margin(c(-10, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))
umap$layers[[1]]$aes_params$alpha <- 0.7

###Look at venoms throughout UMAP
grid.arrange(svsp7, svsp8, svsp1, svsp4,
             svmp4, svmp6, svmp1, svmp5,
             umap.venom,
             layout_matrix = layout_matrix2
)

grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[26], order = TRUE) + labs(title = venomGenes$toxin[26], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[28], order = TRUE) + labs(title = venomGenes$toxin[28], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[29], order = TRUE) + labs(title = venomGenes$toxin[29], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[30], order = TRUE) + labs(title = venomGenes$toxin[30], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[31], order = TRUE) + labs(title = venomGenes$toxin[31], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[33], order = TRUE) + labs(title = venomGenes$toxin[33], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[34], order = TRUE) + labs(title = venomGenes$toxin[34], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme)

grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[36], order = TRUE) + labs(title = venomGenes$toxin[36], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[37], order = TRUE) + labs(title = venomGenes$toxin[37], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[38], order = TRUE) + labs(title = venomGenes$toxin[38], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[39], order = TRUE) + labs(title = venomGenes$toxin[39], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[40], order = TRUE) + labs(title = venomGenes$toxin[40], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[42], order = TRUE) + labs(title = venomGenes$toxin[42], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[43], order = TRUE) + labs(title = venomGenes$toxin[43], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[44], order = TRUE) + labs(title = venomGenes$toxin[44], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme)

### Additional venom genes
grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[30], order = TRUE) + labs(title = venomGenes$toxin[30], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[33], order = TRUE) + labs(title = venomGenes$toxin[33], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[31], order = TRUE) + labs(title = venomGenes$toxin[31], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[42], order = TRUE) + labs(title = venomGenes$toxin[42], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[34], order = TRUE) + labs(title = venomGenes$toxin[34], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[36], order = TRUE) + labs(title = venomGenes$toxin[36], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[37], order = TRUE) + labs(title = venomGenes$toxin[37], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[38], order = TRUE) + labs(title = venomGenes$toxin[38], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[39], order = TRUE) + labs(title = venomGenes$toxin[39], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[40], order = TRUE) + labs(title = venomGenes$toxin[40], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[43], order = TRUE) + labs(title = venomGenes$toxin[43], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom2, features = venomGenes$transcript[44], order = TRUE) + labs(title = venomGenes$toxin[44], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme)


svsp4 <- FeaturePlot(venom, features = venomGenes$transcript[7], order = TRUE) + 
  labs(title = "SVSP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp4$layers[[1]]$aes_params$size <- 0.3

svsp1 <- FeaturePlot(venom, features = venomGenes$transcript[4], order = TRUE) +
  labs(title = "SVSP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp1$layers[[1]]$aes_params$size <- 0.3

svsp8 <- FeaturePlot(venom, features = venomGenes$transcript[11], order = TRUE) +
  labs(title = "SVSP8", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp8$layers[[1]]$aes_params$size <- 0.3

svsp7 <- FeaturePlot(venom, features = venomGenes$transcript[10], order = TRUE) +
  labs(title = "SVSP7", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svsp7$layers[[1]]$aes_params$size <- 0.3

svsp4$layers[[1]]$aes_params$alpha <- 0.5
svsp1$layers[[1]]$aes_params$alpha <- 0.5
svsp8$layers[[1]]$aes_params$alpha <- 0.5
svsp7$layers[[1]]$aes_params$alpha <- 0.5


svmp4 <- FeaturePlot(venom, features = venomGenes$transcript[18], order = TRUE) +
  labs(title = "SVMP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svmp6 <- FeaturePlot(venom, features = venomGenes$transcript[20], order = TRUE) +
  labs(title = "SVMP6", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svmp1 <- FeaturePlot(venom, features = venomGenes$transcript[15], order = TRUE) +
  labs(title = "SVMP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
svmp5 <- FeaturePlot(venom, features = venomGenes$transcript[19], order = TRUE) +
  labs(title = "SVMP5", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme

svmp4$layers[[1]]$aes_params$size <- 0.3
svmp6$layers[[1]]$aes_params$size <- 0.3
svmp1$layers[[1]]$aes_params$size <- 0.3
svmp5$layers[[1]]$aes_params$size <- 0.3

svmp4$layers[[1]]$aes_params$alpha <- 0.5
svmp6$layers[[1]]$aes_params$alpha <- 0.5
svmp1$layers[[1]]$aes_params$alpha <- 0.5
svmp5$layers[[1]]$aes_params$alpha <- 0.5

pla2a1 <- FeaturePlot(venom, features = venomGenes$transcript[29], order = TRUE) + labs(title = venomGenes$toxin[29], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
pla2b1 <- FeaturePlot(venom, features = venomGenes$transcript[26], order = TRUE) + labs(title = venomGenes$toxin[26], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
crisp1 <- FeaturePlot(venom, features = venomGenes$transcript[34], order = TRUE) + labs(title = venomGenes$toxin[34], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
ohanin <- FeaturePlot(venom, features = venomGenes$transcript[36], order = TRUE) + labs(title = venomGenes$toxin[36], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme
exo3 <- FeaturePlot(venom, features = venomGenes$transcript[39], order = TRUE) + labs(title = venomGenes$toxin[39], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme

pla2a1$layers[[1]]$aes_params$size <- 0.3
pla2b1$layers[[1]]$aes_params$size <- 0.3
crisp1$layers[[1]]$aes_params$size <- 0.3
ohanin$layers[[1]]$aes_params$size <- 0.3
exo3$layers[[1]]$aes_params$size <- 0.3

pla2a1$layers[[1]]$aes_params$alpha <- 0.5
pla2b1$layers[[1]]$aes_params$alpha <- 0.5
crisp1$layers[[1]]$aes_params$alpha <- 0.5
ohanin$layers[[1]]$aes_params$alpha <- 0.5
exo3$layers[[1]]$aes_params$alpha <- 0.5

umap.seurat$layers[[1]]$aes_params$alpha <- 0.7

grid.arrange(svsp7, svsp8, svsp1, svsp4,
             svmp4, svmp6, svmp1, svmp5,
             umap.seurat,
             pla2a1, pla2b1, crisp1, ohanin, exo3,
             layout_matrix = layout_matrix2
)

grid.arrange(FeaturePlot(venom, features = venomGenes$transcript[4], order = TRUE) +
               labs(title = "SVSP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[5], order = TRUE) +
               labs(title = "SVSP2", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[6], order = TRUE) +
               labs(title = "SVSP3", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[7], order = TRUE) + 
               labs(title = "SVSP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[8], order = TRUE) +
               labs(title = "SVSP5", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[9], order = TRUE) +
               labs(title = "SVSP6", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[10], order = TRUE) +
               labs(title = "SVSP7", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[11], order = TRUE) +
               labs(title = "SVSP8", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[12], order = TRUE) +
               labs(title = "SVSP9", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme)


grid.arrange(FeaturePlot(venom, features = venomGenes$transcript[15], order = TRUE) +
               labs(title = "SVMP1", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[16], order = TRUE) +
               labs(title = "SVMP2", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[17], order = TRUE) +
               labs(title = "SVMP3", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[18], order = TRUE) +
               labs(title = "SVMP4", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[19], order = TRUE) +
               labs(title = "SVMP5", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[20], order = TRUE) +
               labs(title = "SVMP6", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[21], order = TRUE) + 
               labs(title = "SVMP7", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[22], order = TRUE) + 
               labs(title = "SVMP8", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[23], order = TRUE) + 
               labs(title = "SVMP9", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[24], order = TRUE) + 
               labs(title = "SVMP10", x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             nrow = 3)

grid.arrange(FeaturePlot(venom, features = venomGenes$transcript[26], order = TRUE) + labs(title = venomGenes$toxin[26], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[28], order = TRUE) + labs(title = venomGenes$toxin[28], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[29], order = TRUE) + labs(title = venomGenes$toxin[29], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[30], order = TRUE) + labs(title = venomGenes$toxin[30], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[31], order = TRUE) + labs(title = venomGenes$toxin[31], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[33], order = TRUE) + labs(title = venomGenes$toxin[33], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[34], order = TRUE) + labs(title = venomGenes$toxin[34], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme)

grid.arrange(FeaturePlot(venom, features = venomGenes$transcript[36], order = TRUE) + labs(title = venomGenes$toxin[36], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[37], order = TRUE) + labs(title = venomGenes$toxin[37], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[38], order = TRUE) + labs(title = venomGenes$toxin[38], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[39], order = TRUE) + labs(title = venomGenes$toxin[39], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[40], order = TRUE) + labs(title = venomGenes$toxin[40], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[42], order = TRUE) + labs(title = venomGenes$toxin[42], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[43], order = TRUE) + labs(title = venomGenes$toxin[43], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme,
             FeaturePlot(venom, features = venomGenes$transcript[44], order = TRUE) + labs(title = venomGenes$toxin[44], x = "UMAP 1", y = "UMAP 2") + NoLegend() + featureTheme)


### (3) Gene co-expression correlations ######
venom.sce <- as.SingleCellExperiment(venom2)

sce.hsc <- venom.sce
sce.hsc <- addPerCellQC(sce.hsc)
spike.drop <- quickPerCellQC(colData(sce.hsc))
sce.hsc <- sce.hsc[,!spike.drop$discard]

sce.hsc <- computeSumFactors(sce.hsc)
sce.hsc <- logNormCounts(sce.hsc)

set.seed(100)
var.cor <- correlatePairs(sce.hsc, subset.row = unlist(venomGenes$transcript) %in% row.names(sce.hsc))
gene.cor <- correlateGenes(var.cor)
head(var.cor)

sig.cor <- var.cor[which(var.cor$FDR <= 0.05),]
summary(sig.cor)

var.cor <- readRDS("varcor.RDS")

plotExpression(sce.hsc, features = "fgenesh-scaffold-mi2-venom-gene-3.5", x = "fgenesh-scaffold-mi2-venom-gene-3.7")

vg.cor.all <- as.data.frame(sig.cor[which(sig.cor$gene1 %in% venomGenes$transcript | sig.cor$gene2 %in% venomGenes$transcript),])
all.cor.melt <- melt(vg.cor.all, id.vars = "rho")
all.cor.melt <- all.cor.melt[which(all.cor.melt$variable %in% c("gene1", "gene2")),]
all.cor.melt <- all.cor.melt[which(all.cor.melt$value %in% venomGenes$transcript),]
all.cor.melt$name <- venomGenes[all.cor.melt$value, 2]

vg.cor <- as.data.frame(var.cor[which(var.cor$gene1 %in% venomGenes$transcript & var.cor$gene2 %in% venomGenes$transcript),])
vg.cor <- vg.cor[which(vg.cor$FDR <= 0.05),]
head(vg.cor)
vg.cor$gene1.name <- venomGenes[vg.cor$gene1,2]
vg.cor$gene2.name <- venomGenes[vg.cor$gene2,2]
vg.cor$class1 <- venomClasses[vg.cor$gene1.name, 2]
vg.cor$class2 <- venomClasses[vg.cor$gene2.name, 2]

vg.cor.melt <- melt(vg.cor, id.vars = "rho")
vg.cor.melt <- vg.cor.melt[which(vg.cor.melt$variable %in% c("gene1", "gene2")),]
vg.cor.melt$value
row.names(venomGenes) <- venomGenes$transcript
vg.cor.melt$name <- venomGenes[vg.cor.melt$value,2]

p.vg.cor <- ggplot(vg.cor.melt, aes(x = name, y = rho)) + geom_boxplot() +
  ggtitle("Spearman's rho for venom genes among venom genes") +
  xlab("venom gene") + ylab("rho") +
  theme_bw() +
  theme(text = element_text(size = 8),
        plot.margin = margin(c(2, 2, 2, 2)))
p.all.cor <- ggplot(all.cor.melt, aes(x = name, y = rho)) + geom_boxplot() + #geom_jitter() +
  ggtitle("Spearman's rho for venom genes among all genes") +
  xlab("venom gene") + ylab("rho") +
  theme_bw() +
  theme(text = element_text(size = 8),
        plot.margin = margin(c(2, 2, 2, 2)))

grid.arrange(p.all.cor, p.vg.cor, nrow = 2)

### DOES NOT WORK ON VENOM
library(fcoex)
target <- colData(sce.hsc)
target <- target$seurat_clusters

# Get normalized table from the pre-processing
exprs <- as.data.frame(assay(sce.hsc, 'logcounts'))

# Create fcoex object
fc <- new_fcoex(data.frame(exprs),target)
fc <- discretize(fc, number_of_bins = 8)

fc <- find_cbf_modules(fc,n_genes_selected_in_first_step = 200, verbose = FALSE, is_parallel = FALSE)
fc <- get_nets(fc)

# Taking a look at the first two networks: 
network_plots <- show_net(fc)

### (4) Exploring PLA2 expression patterns ######
pla2.exp <- t(as.data.frame(GetAssayData(venom2, "data")[c("fgenesh-scaffold-mi7-venom-gene-1.1", "fgenesh-scaffold-mi7-venom-gene-1.2", "fgenesh-scaffold-mi7-venom-gene-1.4", "fgenesh-scaffold-mi7-venom-gene-1.5"),]))
pla2.exp <- pla2.exp[which(rowSums2(pla2.exp) > 0),]
pla2.exp <- as.data.frame(pla2.exp)
colnames(pla2.exp) <- c("PLA2G2E", "PLA2B1", "PLA2C1", "PLA2A1")
pla2.exp <- pla2.exp[with(pla2.exp, order(PLA2C1, PLA2B1, PLA2G2E, PLA2A1)),]
pla2.exp$cells <- row.names(pla2.exp)
pla2.exp$cells <- factor(pla2.exp$cells, levels = c(pla2.exp$cells))
pla2.melt <- melt(pla2.exp, id = "cells")

rankLines <- ggplot(pla2.melt, aes(x = `cells`, y = `value`, group = `variable`, color = `variable`)) +
  geom_line() + theme_bw() +
  ggtitle("Cells ranked by expression of PLA2 genes") +
  scale_discrete_manual(values = c("slateblue2", "deepskyblue2", "orchid2", "firebrick2"), aesthetics = c("colour")) +
  ylab("scaled expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(c(-9, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 0, 2)))

pla2g2e <- FeaturePlot(venom2, features = "fgenesh-scaffold-mi7-venom-gene-1.1", order = TRUE) +
  ggtitle("PLA2G2E") +
  scale_color_gradient(low = "#d3d3d3", high = "slateblue3") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.4, 'cm'),
        legend.margin = margin(c(-2, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

pla2b1 <- FeaturePlot(venom2, features = "fgenesh-scaffold-mi7-venom-gene-1.2", order = TRUE) +
  ggtitle("PLA2B1") +
  scale_color_gradient(low = "#d3d3d3", high = "deepskyblue3") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.4, 'cm'),
        legend.margin = margin(c(-2, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

pla2c1 <- FeaturePlot(venom2, features = "fgenesh-scaffold-mi7-venom-gene-1.4", order = TRUE) +
  ggtitle("PLA2C1") +
  scale_color_gradient(low = "#d3d3d3", high = "orchid3") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.4, 'cm'),
        legend.margin = margin(c(-2, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

pla2a1 <- FeaturePlot(venom2, features = "fgenesh-scaffold-mi7-venom-gene-1.5", order = TRUE) +
  ggtitle("PLA2A1") +
  scale_color_gradient(low = "#d3d3d3", high = "firebrick3") + theme_bw() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.4, 'cm'),
        legend.margin = margin(c(-2, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 2, 2)))

p.pla2.genes <-  ggplot(pla2.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = gene)) +
  ggrepel::geom_text_repel(data = pla2.info %>% 
                             filter(str_detect(gene,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=1.75) +
  geom_segment(aes(x=pla2.reg.start,xend=pla2.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  scale_discrete_manual(values = c("firebrick2", "deepskyblue2", "orchid2", "slateblue2", "white"), aesthetics = c("fill")) +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  ggtitle("PLA2 arrangement on microchromosome 7") +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 6),
        plot.margin = margin(c(2, 2, 5, -2)))

layout_matrix <- rbind(c(1, 1, 1, 1),
                       c(2, 2, 2, 2),
                       c(3, 4, 5, 6))

grid.arrange(rankLines,
             p.pla2.genes,
             pla2g2e, pla2b1, pla2c1, pla2a1,
             layout_matrix = layout_matrix,
             heights = c(1, 0.6, 1))

### (5) Various dotplots ######
chromMod <- readxl::read_xlsx("chromMod_transcripts.xlsx", col_names = TRUE)
DotPlot(venom, features = c(chromMod$transcript))

venom.dots <- DotPlot(venom, features = c(venomGenes$transcript)) +
  ggtitle("Venom gene expression by cluster (determined by all expression)") +
  scale_x_discrete(labels = c(venomGenes$toxin[which(venomGenes$transcript %in% row.names(venom))])) +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 6),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        plot.margin = margin(c(2, 2, 5, -2)))

venom2.dots <- DotPlot(venom2, features = c(venomGenes$transcript)) +
  ggtitle("Venom gene expression by cluster (determined by venom expression)") +
  scale_x_discrete(labels = c(venomGenes$toxin[which(venomGenes$transcript %in% row.names(venom))])) +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 6),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        plot.margin = margin(c(2, 2, 5, -2)))

grid.arrange(venom.dots, venom2.dots, ncol = 1)

chrom.dots <- DotPlot(venom, features = c(chromMod$transcript)) +
  ggtitle("Chromatin-modifying genes by seurat cluster") +
  scale_x_discrete(labels = c(chromMod$gene[which(chromMod$transcript %in% row.names(venom))])) +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 6),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        plot.margin = margin(c(2, 2, 5, -2)))

chrom2.dots <- DotPlot(venom2, features = c(chromMod$transcript)) +
  ggtitle("Chromatin-modifying genes by venom cluster") +
  scale_x_discrete(labels = c(chromMod$gene[which(chromMod$transcript %in% row.names(venom2))])) +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 6),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        plot.margin = margin(c(2, 2, 5, -2)))

grid.arrange(chrom.dots, chrom2.dots, nrow = 1)

### (6) SVMP ridges ######
svmp.exp <- t(as.data.frame(GetAssayData(venom2, "data")[c("fgenesh-scaffold-mi1-venom-gene-2.1",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.2",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.3",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.4",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.5",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.6",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.7",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.8",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.9",
                                                           "fgenesh-scaffold-mi1-venom-gene-2.10"),]))
svmp.exp <- svmp.exp[which(rowSums2(svmp.exp) > 0),]
svmp.exp <- as.data.frame(svmp.exp)
colnames(svmp.exp) <- c("SVMP1", "SVMP2", "SVMP3", "SVMP4", "SVMP5",
                        "SVMP6", "SVMP7", "SVMP8", "SVMP9", "SVMP10")
svmp.exp <- svmp.exp[with(svmp.exp, order(SVMP10, SVMP5, SVMP2, SVMP7, SVMP9, SVMP6, SVMP1, SVMP3, SVMP8, SVMP4)),]
svmp.exp$cells <- row.names(svmp.exp)
svmp.exp$cells <- factor(svmp.exp$cells, levels = c(svmp.exp$cells))
svmp.melt <- melt(svmp.exp, id = "cells")

ggplot(svmp.melt, aes(x = `cells`, y = `variable`, height = `value`, group = `variable`, color = `variable`)) +
  geom_ridgeline() + theme_bw() +
  ggtitle("Cells ranked by expression of SVMP genes") +
  #scale_discrete_manual(values = c("slateblue2", "deepskyblue2", "orchid2", "firebrick2"), aesthetics = c("colour")) +
  ylab("scaled expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(c(-9, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 0, 2)))

ggplot(svmp.melt, aes(x = `cells`, y = `value`, group = `variable`, color = `variable`)) +
  geom_line() + theme_bw() +
  ggtitle("Cells ranked by expression of SVMP genes") +
  #scale_discrete_manual(values = c("slateblue2", "deepskyblue2", "orchid2", "firebrick2"), aesthetics = c("colour")) +
  ylab("scaled expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(c(-9, 2, 2, 2)),
        plot.margin = margin(c(2, 2, 0, 2)))

### (7) Exploring transcription factors ######
# These are utilizing lists of genes + transcripts to track down TF transcripts
tf.prot <- read.table("annotations/Cvv_TF.ProtBinding.txt")
tf.prot$V1 <- gsub("_", "-", tf.prot$V1)

tf.dna <- read.table("annotations/Cvv_TF.DNABinding.txt")
tf.dna$V1 <- gsub("_", "-", tf.dna$V1)
tf.dna <- tf.dna[,c(1,3)]
colnames(tf.dna) <- c("V1", "V2")

tf.coreg <- read.table("annotations/Cvv_TF.Coregulators.txt")
tf.coreg$V1 <- gsub("_", "-", tf.coreg$V1)

tf.all <- rbind(tf.prot, tf.dna, tf.coreg)
tf.cand <- read.table("annotations/TF.candidate_list.txt")

cand.dna <- tf.dna[which(tf.dna$V2 %in% tf.cand$V1),]

## adjust FOXO3 annotation bc the same transcript is annotated both FOXO3 and FOXO4
## this was determined based on a grep of the gtf for FOXO3 and FOXO4
cand.dna$V1[which(cand.dna$V2 %in% "FOXO3")] <- "augustus-masked-scaffold-ma1-processed-gene-736.1"

cand.coreg <- tf.coreg[which(tf.coreg$V2 %in% tf.cand$V1),]
cand <- unique(rbind(cand.dna, cand.coreg))

# This is the final candidate list based on TFs from Blair's GR MS``
cand <- cand[with(cand, order(V2)),]

tf.cluster1 <- DotPlot(venom, features = c(cand$V1)) +
  scale_x_discrete(labels = c(cand$V2[which(cand$V1 %in% row.names(venom))])) +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 6),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        plot.margin = margin(c(2, 2, 5, -2)))

tf.cluster2 <- DotPlot(venom2, features = c(cand$V1)) +
  scale_x_discrete(labels = c(cand$V2[which(cand$V1 %in% row.names(venom2))])) +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 6),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        plot.margin = margin(c(2, 2, 5, -2)))

grid.arrange(tf.cluster1, tf.cluster2, nrow = 2)

percentage <- DotPlot(venom2, features= venomGenes$transcript)
View(percentage$data)


length(WhichCells(object = venom2, idents = 1,
                  expression = `fgenesh-scaffold-mi2-venom-gene-3.1` > 0 &
                    `fgenesh-scaffold-mi2-venom-gene-3.3` > 0))/nrow(venom2@meta.data)*100

WhichCells(object = venom2, idents = 3,
           expression = `fgenesh-scaffold-mi2-venom-gene-3.1` > 0 &
             `fgenesh-scaffold-mi2-venom-gene-3.2` > 0)


grid.arrange(UMAPPlot(venom), UMAPPlot(venom, group.by = "venom_clusters") + theme(plot.title = element_blank()), UMAPPlot(venom2), nrow = 1)


UMAPPlot(venom, group.by = "venom_clusters") + theme(plot.title = element_blank())

### (8) Compare TF and venom gene expression ######

sub.var <- var.cor[which(var.cor$gene1 %in% c(venomGenes$transcript, cand$V1) &
                           var.cor$gene2 %in% c(venomGenes$transcript, cand$V1)),]

sub.var2 <- sub.var
colnames(sub.var2) <- c("gene2", "gene1", "rho", "p.value", "FDR")
sub.var3 <- rbind(sub.var, sub.var2)

filtered <- sub.var3[which(sub.var3$gene1 %in% venomGenes$transcript &
                 sub.var3$gene2 %in% cand$V1),]

row.names(cand) <- cand$V1

filtered$gene1 <- venomGenes[filtered$gene1, 2]
filtered$gene2 <- cand[filtered$gene2, 2]

filtered <- as.data.frame(filtered)
ggplot(filtered, aes(x = gene1, y = gene2, fill = rho)) + geom_tile() +
  theme_light() +
  scale_fill_viridis(option = "B") +
  xlab("venom genes") +
  ylab("transcription factors") +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
shaped <- dcast(filtered[,1:3], gene2~gene1, value.var = "rho")
row.names(shaped) <- shaped$gene2
shaped <- shaped[,-1]
cols <- 
pheatmap(shaped,
         color = colorRampPalette(c("white","#FFBFBF", "red"))(40),
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 6,
         cutree_rows = 6
         )

hm.vg_tf <- pheatmap(shaped,
         color = colorRampPalette(c("black", "red"))(40),
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 6,
         cutree_rows = 6,
         treeheight_row = 10,
         treeheight_col = 10,
         cellwidth = 7,
         cellheight = 7,
         gaps_col = c(8,9,10,13,23,32)
)


filtered2 <- sub.var3[which(sub.var3$gene1 %in% venomGenes$transcript &
                             sub.var3$gene2 %in% venomGenes$transcript),]

filtered2$gene1 <- venomGenes[filtered2$gene1, 2]
filtered2$gene2 <- venomGenes[filtered2$gene2, 2]

filtered2 <- as.data.frame(filtered2)
head(filtered2)
pal <- colorRampPalette(c("blue", "black", "red", "yellow"), bias = 1.5)
ggplot(filtered2, aes(x = gene1, y = gene2, fill = rho)) + geom_tile() +
  theme_light() +
  scale_fill_gradientn(colours = pal(20),
                       na.value = "grey50") +
  theme(text = element_text(size = 6),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

shaped2 <- dcast(filtered2[,1:3], gene2~gene1, value.var = "rho")
row.names(shaped2) <- shaped2$gene2
shaped2 <- shaped2[,-1]
vg.order <- read.csv("res.toxin_heatmapOrder.csv", header = FALSE)
shaped2 <- shaped2[vg.order$V1, vg.order$V1]

scale_fill_gradientn(colours = c("blue","black","red"), 
                     values = scales::rescale(c(-.14,0,.51)),
                     guide = "colorbar", limits=c(-.14,.51), breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5))

hm.vg_vg <- pheatmap(shaped2,
         color = colorRampPalette(c("blue", "black", "red", "yellow"), bias = 1.5)(40),
         border_color = NA,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 6,
         treeheight_row = 0,
         treeheight_col = 0,
         cellwidth = 7,
         cellheight = 7
)

grid.arrange(hm.vg_tf[[4]], hm.vg_vg[[4]], nrow = 1)

### (9) Comparing to Blair's MS ######
# in: SVMP6, SVSP5, SVSP7, SVSP10, SVSP8, SVSP6, PLA2A1, PLA2B1, PLA2C1
promoter.tf <- c("ATF4", "CREB3L2", "CREB3L1", "FOXI1", "IRX2", "JUN", "ZBTB26", "ARNT", "SOX10", "TBX3")

enhancer.tf <- c("ARID3A", "ATF4", "TFAP4", "SREBF2", "SREBF1", "DDIT3",
                 "CEBPA", "TFCP2L1", "CREB3L2", "CREB3L1", "CREB3", "SPDEF",
                 "ELF5", "EHF", "FOXO4", "FOXO3", "FOXI1", "FOXC2", "FOS",
                 "GRHL2", "GRHL1", "HES6", "IRX2", "JUN", "ZBTB26", "NFATC1",
                 "NR4A2", "DLX4", "BARX2", "NFIX", "NFIA", "PITX2", "SOX9",
                 "SOX10", "FIGLA", "BHLHA15", "MEIS1", "RORC", "RORA")
