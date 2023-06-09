library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(escape) ### http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html#4_Enrichment

setwd('/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/ERK_UPR')

pathways <- read.csv("ERK-UPR_TFs.csv", header = TRUE)
pathways$gene_ID <- gsub("_masked", "-masked-scaffold", gsub("maker-", "maker-scaffold-", pathways$gene_ID))

unique(pathways$gene_ID[grep('UPR', pathways$Pathway)])

vg.regs <- list(ERK = unique(pathways$gene_ID[grep('ERK', pathways$Pathway)]),
                UPR = unique(pathways$gene_ID[grep('UPR', pathways$Pathway)]))

venom2 <- readRDS("venomClusters_06.23.rds")

venomGenes <- read.csv("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/res.toxin_transcripts.csv", header = TRUE)

row.names(venomGenes) <- venomGenes$transcript


unique(pathways$gene_ID[grep('UPR', pathways$Pathway)])

vg.regs <- list(ERK = unique(pathways$gene_ID[grep('ERK', pathways$Pathway)]),
                UPR = unique(pathways$gene_ID[grep('UPR', pathways$Pathway)]))

ES.venom <- enrichIt(obj = venom2, 
                     gene.sets = vg.regs, 
                     groups = 1000, cores = 2, 
                     min.size = 5)
venom2 <- Seurat::AddMetaData(venom2, ES.venom)

# Try Other approach (k-means clustering)
library(factoextra)

d <- data.frame(venom2$ERK, venom2$UPR)
dscale <- na.omit(d)
dscale <- scale(dscale)
fviz_nbclust(dscale, kmeans, method = "wss")


k=3 # k

dm <- as.matrix(d)
hc <- stats::hclust(dist(dm))
clust <- stats::cutree(hc, k=k)

colourCount = length(levels(as.factor(clust)))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
getPalette(colourCount)

d$clust <- as.factor(clust)


venom2$KMeans <- as.factor(clust)
#table(as.factor(clust)) # Check to see if each cluster has a good number of cells
venom2 <- SetIdent(venom2, value = "KMeans")

#levels(as.factor(clust))[table(as.factor(clust)) >= 20]


meanfunc1 <- function(x) {  
  mean(subset(d, clust==x)$venom2.ERK)
}
meanfunc2 <- function(x) {  
  mean(subset(d, clust==x)$venom2.UPR)
}

text_df <- data.frame(cbind(unlist(lapply(levels(as.factor(clust)), FUN=meanfunc1)),unlist(lapply(levels(as.factor(clust)), FUN=meanfunc2))))
text_df$cluster <- row.names(text_df)
colnames(text_df)[1:2] <- c('Xpos','Ypos')

p4 <- DotPlot(
  object = venom2,
  features=venomGenes$transcript[which(venomGenes$present=='yes')],
  scale=T) + 
  scale_x_discrete(limits=venomGenes$transcript[which(venomGenes$present=='yes')], labels=venomGenes$toxin[which(venomGenes$present=='yes')]) +
  #scale_y_discrete(limits=levels(as.factor(clust))[table(as.factor(clust)) >= 20]) + # Plot clusters with at least 20 cells
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  scale_color_gradient2(low = "blue", mid="grey", high="red") +
  xlab(NULL) +
  ylab('Cluster')

#p1$layers[[1]]
data <- p4$data

px <- ggplot(data, aes(x=features.plot, y=id)) +
  geom_tile(aes(fill = avg.exp.scaled)) +
  ylab("Group") + xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size=6, angle=90, hjust=1, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        title = element_text(face = "bold", size = 7),
        axis.title = element_text(size = 6),
        text = element_text(size = 6)) +
  scale_x_discrete(limits=venomGenes$transcript[which(venomGenes$present=='yes')], labels=venomGenes$toxin[which(venomGenes$present=='yes')]) +
  #scale_fill_gradient2(low = 'blue', mid='grey60', high = 'red')
  scale_fill_gradientn(colours = viridisLite::mako(100)) +
  labs(fill='Scaled\nexpression') 


p5 <- ggplot()+
  geom_point(data=d,aes(x=venom2.ERK, y=venom2.UPR, col=clust)) + 
  theme_bw() + 
  xlab('ERK enrichment') +
  ylab('UPR enrichment') +
  ggtitle("Cells clustered by ERK/UPR") +
  theme(title = element_text(face = "bold", size = 7),
        axis.title = element_text(size = 6),
        text = element_text(size = 6)) +
  scale_color_discrete(type=getPalette(colourCount)) + 
  geom_label(data = text_df, aes(label = cluster,x=Xpos, y = Ypos, size=3), position='dodge') +
  theme(legend.position="none")

r1 <- plot_grid(p1, p2, p3, nrow = 1, labels = c("A", "B", "C"))
r2 <- plot_grid(p5,px, rel_widths = c(1,2), labels = c("D", "E"))
plot_grid(r1, r2, nrow = 2)

new.data <- data %>% select(features.plot, id, avg.exp) %>% pivot_wider(names_from = id, values_from = avg.exp)
data.mat <- new.data[,c(2:4)] %>% as.data.frame()
row.names(data.mat) <- new.data$features.plot
