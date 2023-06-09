##### python snRNAseq analysis #####
##### I have no idea what I am doing but we are going to try some things

### (0) ENVIRONMENT
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(viridis)
library(grid)
library(gridExtra)
library(reshape2)
library(pheatmap)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/")

'%notin%' <- function(x,y)!('%in%'(x,y))

venomGenes <- read.csv("/Users/aundreawestfall/Desktop/Venom_snRNA/res.toxin_transcripts.csv", header = TRUE)
row.names(venomGenes) <- venomGenes$transcript


### (1) TUTORIAL/VIGNETTE #####
##### Data should be gzipped
data <- Read10X(data.dir = "outs/filtered_feature_bc_matrix/")
venom <- CreateSeuratObject(counts = data, project = "venom", min.cells = 3, min.features = 200)

venom[["percent.mt"]] <- PercentageFeatureSet(venom, pattern = "^MT-")
VlnPlot(venom, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol = 3))
FeatureScatter(venom, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

### Normalization
venom <- NormalizeData(venom, normalization.method = "LogNormalize", scale.factor = 10000)

### Identify highly variable features
venom <- FindVariableFeatures(venom, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(venom), 10)

plot1 <- VariableFeaturePlot(venom)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

### Scale data
all.genes <- rownames(venom)
venom <- ScaleData(venom, features = all.genes)

### Linear dimensional reduction
venom2 <- RunPCA(venom, features = venomGenes$transcript)
print(venom2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(venom2, dims = 1:2, reduction = "pca")
DimPlot(venom2, reduction = "pca")
DimHeatmap(venom2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(venom2, dims = 1:15, cells = 500, balanced = TRUE)

### Determine dimensionality
venom2 <- JackStraw(venom2, num.replicate = 100)
venom2 <- ScoreJackStraw(venom2, dims = 1:18)
JackStrawPlot(venom2, dims = 1:18)
ElbowPlot(venom2) ### elbow seems to occur around 6 for venom gland

### Cluster cells
venom2 <- FindNeighbors(venom2, dims = 1:6)
venom2 <- FindClusters(venom2, resolution = 0.5)
head(Idents(venom2), 5)

### Non-linear dimensional reduction
venom2 <- RunUMAP(venom2, dims = 1:6)
DimPlot(venom2, reduction = "umap")

saveRDS(venom, file = "seuratClusters_06.23.rds")
venom <- readRDS(file = "~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/seuratClusters_06.23.rds")

venom2 <- ScaleData(venom2, features = row.names(venom2))

saveRDS(venom2, file = "venomClusters_06.23.rds")
venom2 <- readRDS("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/venomClusters_06.23.rds")

### Identify Biomarkers
venom.markers <- FindAllMarkers(venom2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- venom.markers[row.names(venom.markers) %notin% venomGenes$transcript,] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
top25 <- venom.markers[row.names(venom.markers) %notin% venomGenes$transcript,] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 25, wt = avg_log2FC)
top25_sec <- top25[which(top25$cluster %in% c(0,1,2)),]

venom.markers <- FindAllMarkers(venom2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- venom.markers[row.names(venom.markers) %notin% row.names(venomGenes),] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
top10.all <- print(venom.markers) %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "venomClusterMarkers_06.23.csv", quote = FALSE)
write.csv(top10.all, "venomClusterMarkers_wVenom_06.23.csv", quote = FALSE)

top25 <- venom.markers[row.names(venom.markers) %notin% venomGenes$transcript,] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 25, wt = avg_log2FC)
top25_sec <- top25[which(top25$cluster %in% c(0,1,2)),]

DoHeatmap(venom2, features = top10$gene) + NoLegend() +
  scale_fill_gradientn(colors = scico(50, palette = "berlin"))

svsp2 <- FeaturePlot(venom2, features = venomGenes$transcript[3], order = TRUE) +
  scale_fill_gradientn(colors = scico(50, palette = "devon")) +
  ggtitle("SVSP2") + NoLegend()

FeaturePlot(venom2, features = venomGenes$transcript[24:28], order = TRUE)


mt3 <- FeaturePlot(venom2, features = venomGenes$transcript[1]) +
  scale_fill_gradientn(colors = scico(50, palette = "devon")) +
  ggtitle("myotoxin3") + NoLegend()

bpp <- FeaturePlot(venom2, features = "maker-scaffold-un187-augustus-gene-0.1", order = TRUE) +
  scale_fill_gradientn(colors = scico(50, palette = "devon")) +
  ggtitle("BPP") + NoLegend()

grid.arrange(svsp2, mt3, bpp, nrow = 1)

grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[2], order = TRUE) + ggtitle("SVSP1") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[3], order = TRUE) + ggtitle("SVSP2") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[4], order = TRUE) + ggtitle("SVSP3") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[5], order = TRUE) + ggtitle("SVSP4") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[6], order = TRUE) + ggtitle("SVSP5") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[7], order = TRUE) + ggtitle("SVSP6") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[8], order = TRUE) + ggtitle("SVSP7") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[9], order = TRUE) + ggtitle("SVSP8") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[10], order = TRUE) + ggtitle("SVSP9") + NoLegend())

grid.arrange(FeaturePlot(venom2, features = venomGenes$transcript[13], order = TRUE) + ggtitle("SVMP1") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[14], order = TRUE) + ggtitle("SVMP2") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[15], order = TRUE) + ggtitle("SVMP3") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[16], order = TRUE) + ggtitle("SVMP4") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[17], order = TRUE) + ggtitle("SVMP5") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[18], order = TRUE) + ggtitle("SVMP6") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[19], order = TRUE) + ggtitle("SVMP7") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[20], order = TRUE) + ggtitle("SVMP8") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[21], order = TRUE) + ggtitle("SVMP9") + NoLegend(),
             FeaturePlot(venom2, features = venomGenes$transcript[22], order = TRUE) + ggtitle("SVMP10") + NoLegend(),
             nrow = 3)

grid.arrange(FeaturePlot(venom, features = "maker-scaffold-ma1-augustus-gene-945.8", order = TRUE) +
               ggtitle("A. EPCAM in naive clusters") + theme_minimal() + theme(text = element_text(size = 7)),
             FeaturePlot(venom, features = "maker-scaffold-ma2-augustus-gene-502.9", order = TRUE) +
               ggtitle("B. HEMGN in naive clusters") + theme_minimal() + theme(text = element_text(size = 7)),
             FeaturePlot(venom2, features = "maker-scaffold-ma1-augustus-gene-945.8", order = TRUE) +
               ggtitle("C. EPCAM in venom-based clusters") + theme_minimal() + theme(text = element_text(size = 7)),
             FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-502.9", order = TRUE) +
               ggtitle("D. HEMGN in venom-based clusters") + theme_minimal()  + theme(text = element_text(size = 7)),
             nrow = 2)

# correlations
### https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/misc.html
library(SingleCellExperiment)

venom.sce <- as.SingleCellExperiment(venom2)

library(scater)
sce.hsc <- venom.sce
sce.hsc <- addPerCellQC(sce.hsc)
spike.drop <- quickPerCellQC(colData(sce.hsc))
sce.hsc <- sce.hsc[,!spike.drop$discard]

library(scran)
sce.hsc <- computeSumFactors(sce.hsc)
sce.hsc <- logNormCounts(sce.hsc)

set.seed(100)
var.cor <- correlatePairs(sce.hsc)
   
head(var.cor)

nrow(gene.cor[which(gene.cor$FDR <= 0.05),])
sig.cor <- var.cor[which(var.cor$FDR <= 0.05),]
summary(sig.cor)

var.cor <- readRDS("/Users/aundreawestfall/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/Correlations/varcor_06.22.RDS")
gene.cor <- readRDS("/Users/aundreawestfall/Desktop/snRNAseq_viridis/genecor_06.22.RDS")

plotExpression(sce.hsc, features = "fgenesh-scaffold-mi2-venom-gene-3.5", x = "fgenesh-scaffold-mi2-venom-gene-3.7")

vg.cor.all <- as.data.frame(sig.cor[which(sig.cor$gene1 %in% venomGenes$transcript | sig.cor$gene2 %in% venomGenes$transcript),])
all.cor.melt <- melt(vg.cor.all, id.vars = "rho")
all.cor.melt <- all.cor.melt[which(all.cor.melt$variable %in% c("gene1", "gene2")),]
all.cor.melt <- all.cor.melt[which(all.cor.melt$value %in% venomGenes$transcript),]
all.cor.melt$name <- venomGenes[all.cor.melt$value, 2]

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
row.names(venomClasses) <- venomClasses$x

vg.cor <- as.data.frame(var.cor[which(var.cor$gene1 %in% venomGenes$transcript & var.cor$gene2 %in% venomGenes$transcript),])
vg.cor <- vg.cor[which(vg.cor$FDR <= 0.05),]
vg.cor$gene1.name <- venomGenes[vg.cor$gene1,2]
vg.cor$gene2.name <- venomGenes[vg.cor$gene2,2]
vg.cor$class1 <- venomGenes[vg.cor$gene1, 3]
vg.cor$class2 <- venomGenes[vg.cor$gene2, 3]

vg.cor.melt <- melt(vg.cor, id.vars = "rho")
vg.cor.melt <- vg.cor.melt[which(vg.cor.melt$variable %in% c("gene1", "gene2")),]
vg.cor.melt$value
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

library(matrixStats)
library(reshape2)

### make venomClasses object
expr <- venom@assays$RNA@scale.data %>% as.matrix %>% t %>% as.data.frame
expr <- expr[, venomGenes[which(row.names(venomGenes) %in% colnames(expr)),1]]

features <- venomGenes[which(row.names(venomGenes) %in% colnames(expr)),1]
features2 <- factor(venomGenes[features,2], levels = venomGenes[features,2])

venomClasses <- data.frame(x = features2,
                           group = venomGenes[features, 3])

pla2.exp <- t(as.data.frame(GetAssayData(venom2, "data")[c("fgenesh-scaffold-mi7-venom-gene-1.1", "fgenesh-scaffold-mi7-venom-gene-1.2", "fgenesh-scaffold-mi7-venom-gene-1.4", "fgenesh-scaffold-mi7-venom-gene-1.5"),]))
pla2.exp <- pla2.exp[which(rowSums2(pla2.exp) > 0),]
pla2.exp <- as.data.frame(pla2.exp)
colnames(pla2.exp) <- c("PLA2G2E", "PLA2B1", "PLA2C1", "PLA2A1")
pla2.exp <- pla2.exp[with(pla2.exp, order(PLA2A1, PLA2B1, PLA2G2E, PLA2C1)),]
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

### (8) Compare TF and venom gene expression ######
# These are utilizing lists of genes + transcripts to track down TF transcripts
tf.prot <- read.table("/Users/aundreawestfall/Desktop/snRNAseq_viridis/annotations/Cvv_TF.ProtBinding.txt")
tf.prot$V1 <- gsub("_", "-", tf.prot$V1)

tf.dna <- read.table("/Users/aundreawestfall/Desktop/snRNAseq_viridis/annotations/Cvv_TF.DNABinding.txt")
tf.dna$V1 <- gsub("_", "-", tf.dna$V1)
tf.dna <- tf.dna[,c(1,3)]
colnames(tf.dna) <- c("V1", "V2")

tf.coreg <- read.table("/Users/aundreawestfall/Desktop/snRNAseq_viridis/annotations/Cvv_TF.Coregulators.txt")
tf.coreg$V1 <- gsub("_", "-", tf.coreg$V1)

tf.all <- rbind(tf.prot, tf.dna, tf.coreg)
tf.cand <- read.table("/Users/aundreawestfall/Desktop/snRNAseq_viridis/annotations/TF.candidate_list.txt")

cand.dna <- tf.dna[which(tf.dna$V2 %in% tf.cand$V1),]

## adjust FOXO3 annotation bc the same transcript is annotated both FOXO3 and FOXO4
## this was determined based on a grep of the gtf for FOXO3 and FOXO4
cand.dna$V1[which(cand.dna$V2 %in% "FOXO3")] <- "augustus-masked-scaffold-ma1-processed-gene-736.1"

cand.coreg <- tf.coreg[which(tf.coreg$V2 %in% tf.cand$V1),]
cand <- unique(rbind(cand.dna, cand.coreg))

# This is the final candidate list based on TFs from Blair's GR MS``
cand <- cand[with(cand, order(V2)),]
saveRDS(cand, "~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/candidateTFs.rds")
cand <- readRDS("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/candidateTFs.rds")

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

tfs <- c("GRHL1", "TFCP2L1",
         "IRX2", "XBP1", "ATF2", "ATF4", "ATF6",
         "PITX2", "SPDEF", "KLF11", "NCOA2", "SP1", "SREBF2", "NFIX", "SREBF1", "NFIB", "NFIA", "NR4A2")
vg.order <- read.csv("/Users/aundreawestfall/Desktop/snRNAseq_viridis/res.toxin_heatmapOrder.csv", header = FALSE)
vgs <- vg.order$V1

row_breaks <- grep("TFCP2L1|ATF6", tfs)
col_breaks <- grep("myotoxin|SVSP9|SVMP10|PLA2C1", vgs)
expr <- as.matrix(t(shaped[tfs, vgs]))
hm.vg_tf <- pheatmap(expr,
                     color = magma(40),
                     border_color = NA,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     fontsize = 6,
                     cutree_rows = 3,
                     scale = "row",
                     gaps_row = col_breaks,
                     gaps_col = row_breaks,
                     main = NA,
                     legend =FALSE
)
hm.vg_tf

data <- melt(expr)
colnames(data) <- c("toxin", "tf", "rho")

venomGenes2 <- venomGenes
row.names(venomGenes2) <- venomGenes2$toxin
data$toxin_family <- venomGenes2$family[data$toxin]
data$toxin_family <- gsub("BPP|LAAO|CRISP|ohanin|EXO|VEGF|CTL|GC|disintegrin|CTL", "other", data$toxin_family)

tfGroups <- data.frame(tfs = tfs,
                       group = c(rep("pioneer", 2), rep("UPR", 5), rep("ERK", 11)))
row.names(tfGroups) <- tfGroups$tfs

temp <- venomGenes2[as.vector(data$toxin),3]
data$toxin_family <- temp
data$tf_group <- tfGroups$group[data$tf]
data$toxin_family <- gsub("BPP|LAAO|CRISP|ohanin|EXO|VEGF|CTL|GC|disintegrin|CTL", "other", data$toxin_family)

devtools::install_github("tidyverse/ggplot2")
ggplot(data, aes(x = tf, y = toxin, fill = rho)) + geom_tile() +
  theme_minimal() + coord_cartesian(expand = FALSE) +
  ggtitle("TFs covary with venom\ngene expression") +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "black", "red"))(20),
                       na.value = "grey50") +
  facet_grid(rows = vars(toxin_family), cols = vars(tf_group), scales = "free", space = "free", shrink = FALSE) +
  theme(text = element_text(size = 6),
        plot.title = element_text(size = 7, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.position = "right",
        legend.margin = margin(c(0, 0, 0, 0)),
        plot.margin = margin(c(2, 2, 2, 2)),
        strip.text.x=element_text(margin=margin(r=10, l = 10)),
        strip.background = element_rect(fill = "transparent"),
        strip.clip = "off",
        panel.grid = element_blank(),
        panel.border = element_blank())
f3.heatmap

filtered2 <- sub.var3[which(sub.var3$gene1 %in% venomGenes$transcript &
                              sub.var3$gene2 %in% venomGenes$transcript),]

filtered2$gene1 <- venomGenes[filtered2$gene1, 2]
filtered2$gene2 <- venomGenes[filtered2$gene2, 2]

filtered2 <- as.data.frame(filtered2)

library(escape)
vg.sets2 <- list(SVSP = venomGenes$transcript[which(venomGenes$family == "SVSP")],
                  SVMP = venomGenes$transcript[which(venomGenes$family == "SVMP")],
                  PLA2 = c("fgenesh-scaffold-mi7-venom-gene-1.2", "fgenesh-scaffold-mi7-venom-gene-1.4", "fgenesh-scaffold-mi7-venom-gene-1.5"),
                  CTL = venomGenes$transcript[which(venomGenes$family == "CTL")],
                  LAAO = venomGenes$transcript[which(venomGenes$family == "LAAO")],
                  EXO = venomGenes$transcript[which(venomGenes$family == "EXO")])

vg.regs <- list(ERK = pathways$gene_ID[which(pathways$Pathway == "ERK")],
               UPR = pathways$gene_ID[which(pathways$Pathway == "UPR")])

ES.venom <- enrichIt(obj = venom2, 
                      gene.sets = vg.regs, 
                      groups = 1000, cores = 2, 
                      min.size = 5)
venom2 <- Seurat::AddMetaData(venom2, ES.venom)
plot_grid(FeaturePlot(venom2, features = "ERK", order = TRUE),
          FeaturePlot(venom2, features = "UPR", order = TRUE))

VlnPlot(venom2, features = c("ERK", "UPR"), group.by = "seurat_clusters")

dittoScatterPlot(venom2, x.var = "maker-scaffold-ma6-augustus-gene-163.6", 
                 y.var = "myotoxin3", 
                 do.contour = FALSE) 

grid.arrange(
  FeaturePlot(venom2, features = "SVSP") + scale_color_viridis(),
  FeaturePlot(venom, features = "SVSP") + scale_color_viridis(),
  FeaturePlot(venom2, features = "SVMP") + scale_color_viridis(),
  FeaturePlot(venom, features = "SVMP") + scale_color_viridis(),
  nrow = 2)


