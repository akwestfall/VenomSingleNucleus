'%notin%' <- function(x,y)!('%in%'(x,y))

venom_sub <- subset(venom, seurat_clusters %in% c(0:2))

### Normalization
venom_sub <- NormalizeData(venom_sub, normalization.method = "LogNormalize", scale.factor = 10000)

### Identify highly variable features
venom_sub <- FindVariableFeatures(venom_sub, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(venom_sub), 10)

plot1 <- VariableFeaturePlot(venom_sub)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

### Scale data
all.genes <- rownames(venom_sub)
venom_sub <- ScaleData(venom_sub, features = all.genes)

### Linear dimensional reduction
venom_sub <- RunPCA(venom_sub, features = venomGenes$transcript)
print(venom_sub[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(venom_sub, dims = 1:2, reduction = "pca")
DimPlot(venom_sub, reduction = "pca")
DimHeatmap(venom_sub, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(venom_sub, dims = 1:15, cells = 500, balanced = TRUE)

### Determine dimensionality
venom_sub <- JackStraw(venom_sub, num.replicate = 100)
venom_sub <- ScoreJackStraw(venom_sub, dims = 1:18)
JackStrawPlot(venom_sub, dims = 1:18)
ElbowPlot(venom_sub) ### elbow seems to occur around 6 for venom gland

### Cluster cells
venom_sub <- FindNeighbors(venom_sub, dims = 1:6, graph.name = "test")
venom_sub <- FindClusters(venom_sub, resolution = 0.5)
head(Idents(venom_sub), 5)

### Non-linear dimensional reduction
venom_sub <- RunUMAP(venom_sub, dims = 1:6)
DimPlot(venom_sub, reduction = "umap")

venom.markers <- FindAllMarkers(venom2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
top10 <- venom.markers[row.names(venom.markers) %notin% venomGenes$transcript,] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)

venom2$seurat_clusters <- as.numeric(venom2$seurat_clusters)
venom2$seurat_clusters <- factor(venom2$seurat_clusters)
Idents(venom2) <- venom2$seurat_clusters

venom2 <- FindSubCluster(venom2, 1, "RNA_nn",
                         subcluster.name = "clust1",
                         resolution = 0.5,
                         algorithm = 1)
venom2 <- FindSubCluster(venom2, 2, "RNA_nn",
                         subcluster.name = "clust2",
                         resolution = 0.5,
                         algorithm = 1)
venom2 <- FindSubCluster(venom2, 3, "RNA_nn",
                         subcluster.name = "clust3",
                         resolution = 0.5,
                         algorithm = 1)
venom2 <- FindSubCluster(venom2, 4, "RNA_nn",
                         subcluster.name = "clust4",
                         resolution = 0.5,
                         algorithm = 1)


subclust1 <- UMAPPlot(venom2, group.by = "clust1", cols = c("red", "blue", "orange", "grey75", "grey50", "grey25"))
subclust2 <- UMAPPlot(venom2, group.by = "clust2", cols = c("grey75", "red", "blue", "orange", "cyan", "grey50", "grey25"))
subclust3 <- UMAPPlot(venom2, group.by = "clust3", cols = c("grey75", "grey50", "red", "blue", "grey25"))
subclust4 <- UMAPPlot(venom2, group.by = "clust4", cols = c("grey75", "grey50", "grey25", "red", "blue", "orange"))

plot_grid(subclust1 + theme_void() + ggtitle("Cluster 1") +
  theme(plot.title = element_text(hjust = 0.5)),
subclust2 + theme_void() + ggtitle("Cluster 2") +
  theme(plot.title = element_text(hjust = 0.5)),
subclust3 + theme_void() + ggtitle("Cluster 3") +
  theme(plot.title = element_text(hjust = 0.5)),
subclust4 + theme_void() + ggtitle("Cluster 4") +
  theme(plot.title = element_text(hjust = 0.5)))

saveRDS(venom2, "~/Desktop/subclustered_venom2.rds")

test <- venom2@meta.data[, c("clust1", "clust2", "clust3", "clust4")]
test[which(test$clust1 %in% c(2:4)),"clust1"] <- NA
test[which(test$clust2 %in% c(1, 3, 4)),"clust2"] <- NA
test[which(test$clust3 %in% c(1, 2, 4)),"clust3"] <- NA
test[which(test$clust4 %in% c(1:3)),"clust4"] <- NA

idents <- test %>% mutate(mycol = coalesce(clust1, clust2, clust3, clust4)) %>%
  select(mycol)
venom2$subclusts <- idents
Idents(venom2) <- venom2$subclusts
venom2@active.ident <- factor(x = venom2@active.ident, levels = c("1_0", "1_1", "1_2",
                                                                  "2_0", "2_1", "2_2", "2_3",
                                                                  "3_0", "3_1",
                                                                  "4_0", "4_1", "4_2"))

venom2 <- SCTransform(venom2, assay = "RNA", new.assay.name = "SCT")
venom.markers <- FindAllMarkers(venom2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- venom.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
DoHeatmap(venom2, top10$gene)


library(CellChat)
library(patchwork)
data.input <- GetAssayData(venom2, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(venom2)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchat <- addMeta(cellchat, meta = meta, meta.name = "group")
cellchat <- setIdent(cellchat, ident.use = "group")

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) 

umap <- venom2@reductions$umap@cell.embeddings[,1:2]
umap <- as.data.frame(umap)
umap$subs <- venom2@meta.data$subclusts
f3.heatmap <- ggplot(filtered2, aes(x = gene1, y = gene2, fill = rho)) + geom_tile(width = 1.1, height = 1.1) +
  theme_light() + coord_fixed(1, expand = FALSE) +
  ggtitle("Across-cell correlations\nin venom expression") +
  scale_fill_gradientn(colours = pal(20),
                       na.value = "grey50") +
  theme(text = element_text(size = 5),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(0,0,1,0)),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.position = "right",
        legend.margin = margin(c(0, 0, 0, 0)),
        plot.margin = margin(c(2, 2, 2, 2)),
        panel.grid = element_blank(),
        panel.border = element_blank())

subUmap <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = subs)) + geom_point(size = 0.5) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Cellular heterogenenity of\nhigh venom family expression") +
  theme_void() + 
  scale_color_manual( #values = c("red", "pink", "purple", "yellow", "orange", "cyan", "blue", "green"),
                      values = c("#66CD00", "#009ACD", "#FF4500", "#9A32CD", "#F5E400", "#D5006A", "#08585A", "#B1740F"),
                      na.value = "grey75",
                      breaks = c("1_0", "1_1", "1_2", "2_0", "2_1", "2_2",
                                 "3_1", "4_2"),
                      labels = c("SVMP", "SVMP/CTL", "SVMP/CTL/PLA2A1",
                                 "PLA2B1", "SVMP/PLA2B1", "SVSP/LAAO/CRISP", 
                                 "PLA2A1", "myotoxin")) + 
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(angle = 90),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
        plot.margin = margin(0, 0, 2, 0),
        legend.title = element_text(face = "bold"),
        legend.position = "right") + coord_fixed(5/6)
subUmap$layers[[1]]$aes_params$size <- 2
subUmap$layers[[1]]$aes_params$alpha <- 1


plot_grid(subUmap, f3.heatmap, nrow = 2, rel_heights = c(1.1, 1.3))


