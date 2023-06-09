library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(escape) ### http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html#4_Enrichment

setwd('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/ERK_UPR')

pathways <- read.csv("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/ERK_UPR/ERK-UPR_TFs.csv", header = TRUE)
pathways$gene_ID <- gsub("_masked", "-masked-scaffold", gsub("maker-", "maker-scaffold-", pathways$gene_ID))

unique(pathways$gene_ID[grep('UPR', pathways$Pathway)])

vg.regs <- list(ERK = unique(pathways$gene_ID[grep('ERK', pathways$Pathway)]),
                UPR = unique(pathways$gene_ID[grep('UPR', pathways$Pathway)]))

venom2 <- readRDS("venomClusters_06.23.rds")
naive <- readRDS("seuratClusters_06.23.rds")

venomGenes <- read.csv("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/res.toxin_transcripts.csv", header = TRUE)

row.names(venomGenes) <- venomGenes$transcript

# p1 <- FeaturePlot(venom2, pathways$gene_ID[1], order = TRUE) + ggtitle(pathways$Name[1]) + scale_color_viridis()
# p2 <- FeaturePlot(venom2, pathways$gene_ID[2], order = TRUE) + ggtitle(pathways$Name[2]) + scale_color_viridis()
# p3 <- FeaturePlot(venom2, pathways$gene_ID[3], order = TRUE) + ggtitle(pathways$Name[3]) + scale_color_viridis()
# p4 <- FeaturePlot(venom2, pathways$gene_ID[4], order = TRUE) + ggtitle(pathways$Name[4]) + scale_color_viridis()
# p5 <- FeaturePlot(venom2, pathways$gene_ID[5], order = TRUE) + ggtitle(pathways$Name[5]) + scale_color_viridis()
# p6 <- FeaturePlot(venom2, pathways$gene_ID[6], order = TRUE) + ggtitle(pathways$Name[6]) + scale_color_viridis()
# plot_grid(p1, p2, p3, p4, p5, p6)

unique(pathways$gene_ID[grep('UPR', pathways$Pathway)])

vg.regs <- list(ERK = unique(pathways$gene_ID[grep('ERK', pathways$Pathway)]),
                UPR = unique(pathways$gene_ID[grep('UPR', pathways$Pathway)]))
#both.erk.upr <- list(ERK_UPR = unique(pathways$gene_ID))

# ES.venom <- enrichIt(obj = venom2, 
#                      gene.sets = both.erk.upr, 
#                      groups = 1000, cores = 2, 
#                      min.size = 5)
# venom2 <- Seurat::AddMetaData(venom2, ES.venom)

ES.venom <- enrichIt(obj = venom2, 
                     gene.sets = vg.regs, 
                     groups = 1000, cores = 2, 
                     min.size = 5)
venom2 <- Seurat::AddMetaData(venom2, ES.venom)


vg.sets3 <- list(SVSP = venomGenes$transcript[which(venomGenes$family == "SVSP")],
                 SVMP = venomGenes$transcript[which(venomGenes$family == "SVMP")],
                 PLA2 = c("fgenesh-scaffold-mi7-venom-gene-1.2", "fgenesh-scaffold-mi7-venom-gene-1.4", "fgenesh-scaffold-mi7-venom-gene-1.5"),
                 CTL = venomGenes$transcript[which(venomGenes$family == "CTL")],
                 LAAO = venomGenes$transcript[which(venomGenes$family == "LAAO")],
                 EXO = venomGenes$transcript[which(venomGenes$family == "EXO")])
ES.venom <- enrichIt(obj = venom2, 
                     gene.sets = vg.sets3, 
                     groups = 1000, cores = 2, 
                     min.size = 3)
venom2 <- Seurat::AddMetaData(venom2, ES.venom)

vg.sets2 <- list(VenomGenes=venomGenes$transcript)
ES.venom <- enrichIt(obj = venom2, 
                     gene.sets = vg.sets2, 
                     groups = 1000, cores = 2, 
                     min.size = 3)
venom2 <- Seurat::AddMetaData(venom2, ES.venom)




#### Trying UPR key genes from Aundrea ####
# XBP1, PDIA5, EIF2AK3, NFE2L2, CCND1, ERN1, ATF6, MAPK8, PPP1R15A, ATF4, DDIT3
#upr_aundrea <- list(UPR2=pathways$gene_ID[pathways$Pathway=='UPR_aundrea'])
upr_aundrea <- list(UPR2=unique(pathways$gene_ID[pathways$Name=='EIF2AK3'|pathways$Name=='ERN1'|pathways$Name=='ATF6'|pathways$Name=='XBP1'])) # Only EIF2AK3, ERN1, ATF6 ## AKW: added XBP1
ES.venom <- enrichIt(obj = venom2, 
                     gene.sets = upr_aundrea, 
                     groups = 1000, cores = 2, 
                     min.size = 3)
venom2 <- Seurat::AddMetaData(venom2, ES.venom)


plot_grid(FeaturePlot(venom2, features = "UPR2", order = TRUE) + scale_color_gradientn(colours=c('lightgrey', 'red'), limits = range(venom2@meta.data$UPR2, venom2@meta.data$UPR2)) + ggtitle('UPR'),
          FeaturePlot(venom2, features = "augustus-masked-scaffold-ma3-processed-gene-300.3", order = TRUE) + ggtitle('ATF6'),
          FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-362.6", order = TRUE) + ggtitle('ERN1'),
          FeaturePlot(venom2, features = "maker-scaffold-Z-augustus-gene-255.8", order = TRUE) + ggtitle('EIF2AK3'), ncol=1)

umap.upr <- FeaturePlot(venom2, features = "UPR2", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("D. UPR Enrichment") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "red"), limits = range(venom2@meta.data$UPR2, venom2@meta.data$UPR2)) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")

umap.xbp1 <- FeaturePlot(venom2, features = "maker-scaffold-mi8-augustus-gene-35.54", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("XBP1") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")
umap.atf6 <- FeaturePlot(venom2, features = "augustus-masked-scaffold-ma3-processed-gene-300.3", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ATF6") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")
umap.ern1 <- FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-362.6", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ERN1") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")
umap.eif2ak3 <- FeaturePlot(venom2, features = "maker-scaffold-Z-augustus-gene-255.8", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("EIF2AK3") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")

upr2 <- plot_grid(umap.xbp1, umap.atf6, umap.ern1, umap.eif2ak3)
upr.grid <- plot_grid(umap.upr, upr2, ncol = 1)

# ERK vs UPR
# library(reshape2)
# venom2@meta.data %>% select(ERK,UPR) %>% melt %>% 
#   ggplot(aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))
# venom2@meta.data %>% select(ERK,UPR) %>% melt()

# random_genes <- list(RandomGenes=sample(x=rownames(venom2)[!(rownames(venom2) %in% venomGenes$transcript)], size=5000, replace = FALSE))
# ES.venom <- enrichIt(obj = venom2, 
#                      gene.sets = random_genes, 
#                      groups = 1000, cores = 2, 
#                      min.size = 3)




erk.enrich <- ggplot(venom2@meta.data, aes(x = ERK, 
                                           y = VenomGenes)) + 
  geom_hex(bins = 50) + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  #stat_cor(method = "pearson", label.x = -2000, label.y = 1) +
  xlab('ERK') +
  ylab('venom genes') +
  theme_cowplot(12) +
  scale_fill_viridis_c()

upr.enrich <- ggplot(venom2@meta.data, aes(x = UPR, 
                             y = VenomGenes)) + 
  geom_hex(bins = 50) + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  #stat_cor(method = "pearson", label.x = -2000, label.y = 1) +
  xlab('UPR') +
  ylab('venom genes') +
  theme_cowplot(12) +
  scale_fill_viridis_c()


####
#### Trying EKR key genes #### https://www.rndsystems.com/resources/articles/erk-signal-transduction-pathway
FeaturePlot(venom2, features = "augustus-masked-scaffold-ma1-processed-gene-750.0", order = TRUE) + ggtitle('TBX3')
FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-434.11", order = TRUE) + ggtitle('ELK1') # Marais, R. et al. (1997) J. Biol. Chem. 272:4378.
FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-372.8", order = TRUE) + ggtitle('ATF1') # Gupta, P. & R. Prywes (2002) J. Biol. Chem. 277:50550.
FeaturePlot(venom2, features = "maker-scaffold-ma1-augustus-gene-935.8", order = TRUE) + ggtitle('SRF') # Xing, J. et al. (1996) Science 273:959.

erk_2 <- list(ERK2=unique(pathways$gene_ID[pathways$Name=='ELK1'|pathways$Name=='ATF1'|pathways$Name=='SRF'])) # Only ELK1, ATF1, SRF
ES.venom <- enrichIt(obj = venom2, 
                     gene.sets = erk_2, 
                     groups = 1000, cores = 2, 
                     min.size = 3)
venom2 <- Seurat::AddMetaData(venom2, ES.venom)

plot_grid(FeaturePlot(venom2, features = "ERK2", order = TRUE) + scale_color_gradient(low='lightgrey', high= 'red', limits = range(venom2@meta.data$UPR2, venom2@meta.data$UPR2)) + ggtitle('ERK'),
          FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-434.11", order = TRUE) + ggtitle('ELK1'),
          FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-372.8", order = TRUE) + ggtitle('ATF1'),
          FeaturePlot(venom2, features = "augustus-masked-scaffold-ma1-processed-gene-750.0", order = TRUE) + ggtitle('TBX3'),
          FeaturePlot(venom2, features = "maker-scaffold-ma1-augustus-gene-935.8", order = TRUE) + ggtitle('SRF'), ncol=1)

umap.erk <- FeaturePlot(venom2, features = "ERK2", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("C. ERK Enrichment") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "red"), limits = range(venom2@meta.data$ERK2, venom2@meta.data$ERK2)) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")

umap.elk1 <- FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-434.11", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ELK1") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")
umap.atf1 <- FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-372.8", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ATF1") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")
umap.tbx3 <- FeaturePlot(venom2, features = "augustus-masked-scaffold-ma1-processed-gene-750.0", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("TBX3") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")
umap.srf <- FeaturePlot(venom2, features = "maker-scaffold-ma1-augustus-gene-935.8", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("SRF") +
  coord_cartesian(xlim=c(-5, 5), ylim=c(-4, 4), expand = FALSE, clip = "off") +
  annotate(x=-5, xend=-5, y=-4, yend=0, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-5, xend=-1, y=-4, yend=-4, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(legend.key.width = unit(0.1, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.margin = margin(c(0,0,0,0)),
                          plot.margin = margin(c(0,0,0,0)),
                          legend.position = "right")

p1 <- erk.enrich + theme_bw() + scale_fill_viridis(option = "magma") +
  ggtitle("A. ERK and venom gene enrichment") +
  theme(legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.position = "right",
        legend.margin = margin(c(0,0,0,0)),
        plot.margin = margin(c(0,0,0,0)),
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 6))

p2 <- upr.enrich + theme_bw() + scale_fill_viridis(option = "magma") +
  ggtitle("B. UPR and venom gene enrichment") +
  theme(legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.position = "right",
        legend.margin = margin(c(0,0,0,0)),
        plot.margin = margin(c(0,0,0,0)),
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 6))

p3 <- ggplot(venom2@meta.data, aes(x = ERK, 
                             y = UPR)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  ggpubr::stat_cor(method = "pearson", label.x = 500, label.y = -500, size = 2.14) +
  xlab('ERK enrichment') +
  ylab('UPR gene enrichment') +
  ggtitle("F. ERK and UPR enrichment correlation") +
  theme_bw() +
  theme(legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.position = "right",
        legend.margin = margin(c(0,0,0,0)),
        plot.margin = margin(c(0,5.5,0,0)),
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 6))


erk2 <- plot_grid(umap.elk1, umap.atf1, umap.tbx3, umap.srf)
erk.grid <- plot_grid(p1, umap.erk, erk2, ncol = 1)

upr2 <- plot_grid(umap.xbp1, umap.atf6, umap.ern1, umap.eif2ak3) 
upr.grid <- plot_grid(p2, umap.upr, upr2, ncol = 1) 

enrichment <- plot_grid(erk.grid, upr.grid, nrow = 1)

data <- readRDS("~/Desktop/heatmapData.rds")
data$toxin_family <- gsub("UPR", " UPR", data$toxin_family)
htmp <- ggplot(data, aes(x = tf, y = toxin, fill = rho)) + geom_tile(width = 1.1, height = 1.1) +
  theme_minimal() + coord_cartesian(expand = FALSE) +
  ggtitle("TFs covary with venom\ngene expression") +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "black", "red"))(20), na.value = "grey50") +
  # scale_fill_viridis() +
  facet_grid(rows = vars(toxin_family), cols = vars(tf_group), scales = "free", space = "free", shrink = FALSE) +
  theme(text = element_text(size = 5),
        plot.title = element_text(size = 7, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.position = "right",
        legend.margin = margin(c(0, 0, 0, -5)),
        panel.spacing = unit(0.02, "lines"),
        plot.margin = margin(c(2, 2, 2, 2)),
        strip.text=element_text(margin=margin(t = 2, r = 2, b = 2, l = 0)),
        panel.grid = element_blank(),
        panel.border = element_blank())

p5 <- plot_grid(enrichment, htmp, ncol = 2, rel_widths = c(1.65,1))

ggsave("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/ms/figurePanels/Fig5_TFs.pdf", plot = p5, width = 172, height = 187, units = "mm", dpi = 300)

####
library(dittoSeq)
library(ggpubr)
# ### You can fill in any x.var feature or gene
dittoScatterPlot(venom2, x.var = "ERK",
                 y.var = "UPR",
                 do.contour = FALSE)
# 
# dittoScatterHex(venom2, x.var = "ERK_UPR", 
#                  y.var = "VenomGenes", 
#                  do.contour = FALSE)

options(scipen=0)
p1 <- ggplot(venom2@meta.data, aes(x = ERK, 
                     y = VenomGenes)) + 
  geom_hex(bins = 50) + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, size = 2.14) +
  xlab('ERK enrichment') +
  ylab('venom gene enrichment') +
  ggtitle("A. ERK and venom gene enrichment") +
  theme_cowplot(12) +
  scale_fill_viridis_c()

p2 <- ggplot(venom2@meta.data, aes(x = UPR, 
                                   y = VenomGenes)) + 
  geom_hex(bins = 50) + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  stat_cor(method = "pearson", label.x = 500, label.y = 1000, size = 2.14) +
  xlab('UPR enrichment') +
  ylab('venom gene enrichment') +
  ggtitle("B. UPR and venom gene enrichment") +
  theme_cowplot(12) +
  scale_fill_viridis_c()


# ERK plot
p3 <- plot_grid(FeaturePlot(venom2, features = "ERK2", order = TRUE) + scale_color_gradient(low='lightgrey', high= 'red', limits = range(0, venom2@meta.data$ERK2)) + ggtitle('ERK'),
          FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-434.11", order = TRUE) + ggtitle('ELK1'),
          FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-372.8", order = TRUE) + ggtitle('ATF1'),
          FeaturePlot(venom2, features = "maker-scaffold-ma1-augustus-gene-935.8", order = TRUE) + ggtitle('SRF'), ncol=1)
# UPR plots
p4 <- plot_grid(FeaturePlot(venom2, features = "UPR2", order = TRUE) + scale_color_gradientn(colours=c('lightgrey', 'red'), limits = range(0, venom2@meta.data$UPR2)) + ggtitle('UPR'),
                FeaturePlot(venom2, features = "augustus-masked-scaffold-ma3-processed-gene-300.3", order = TRUE) + ggtitle('ATF6'),
                FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-362.6", order = TRUE) + ggtitle('ERN1'),
                FeaturePlot(venom2, features = "maker-scaffold-Z-augustus-gene-255.8", order = TRUE) + ggtitle('EIF2AK3'), ncol=1)

plot_grid(p1,p2,p3,p4, rel_heights = c(1,4,4,4,4))



#### Naive clusters ####
ES.venom <- enrichIt(obj = naive, 
                     gene.sets = both.erk.upr, 
                     groups = 1000, cores = 2, 
                     min.size = 5)

naive <- Seurat::AddMetaData(naive, ES.venom)

ES.venom <- enrichIt(obj = naive, 
                     gene.sets = vg.regs, 
                     groups = 1000, cores = 2, 
                     min.size = 5)

naive <- Seurat::AddMetaData(naive, ES.venom)

vg.sets2 <- list(VenomGenes=venomGenes$transcript)

ES.venom <- enrichIt(obj = naive, 
                     gene.sets = vg.sets2, 
                     groups = 1000, cores = 2, 
                     min.size = 3)
naive <- Seurat::AddMetaData(naive, ES.venom)

ES.venom <- enrichIt(obj = naive, 
                     gene.sets = vg.sets3, 
                     groups = 1000, cores = 2, 
                     min.size = 3)
naive <- Seurat::AddMetaData(naive, ES.venom)



p1 <- ggplot(naive@meta.data, aes(x = ERK, 
                                   y = VenomGenes)) + 
  geom_hex(bins = 50) + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  stat_cor(method = "pearson", label.x = 500, label.y = 1000) +
  xlab('ERK enrichment') +
  ylab('Venom gene enrichment') +
  theme_cowplot(12) +
  scale_fill_viridis_c()

p2 <- ggplot(naive@meta.data, aes(x = UPR, 
                                   y = VenomGenes)) + 
  geom_hex(bins = 50) + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  stat_cor(method = "pearson", label.x = 500, label.y = 1000) +
  xlab('UPR enrichment') +
  ylab('Venom gene enrichment') +
  theme_cowplot(12) +
  scale_fill_viridis_c()

p3 <- ggplot(naive@meta.data, aes(x = ERK_UPR, 
                                   y = VenomGenes)) + 
  geom_hex(bins = 50) + 
  geom_smooth(method = 'lm', se = F, color = 'red') +
  stat_cor(method = "pearson", label.x = 500, label.y = 1000) +
  xlab('Joint ERK-UPR enrichment') +
  ylab('Venom gene enrichment') +
  theme_cowplot(12) +
  scale_fill_viridis_c()

plot_grid(p1,p2,p3,nrow=1,ncol=3)

p4 <- FeaturePlot(naive, features = c("ERK","UPR","ERK_UPR"), order = TRUE, combine=F)
p5 <- FeaturePlot(naive, features = "UPR", order = TRUE) + scale_fill_gradient(limits = range(naive@meta.data$UPR, naive@meta.data$UPR))
p6 <- FeaturePlot(naive, features = "ERK_UPR", order = TRUE) + scale_fill_gradient(limits = range(naive@meta.data$UPR, naive@meta.data$UPR))

p5 <- lapply(p4, function (x) x + scale_color_gradientn(colors=c('lightgrey','blue'), limits = range(naive@meta.data$UPR, naive@meta.data$UPR)))
CombinePlots(p5)


scale_fill_gradient(limits = range(venom2@meta.data$UPR, venom2@meta.data$UPR))

max(venom2@meta.data$ERK_UPR)

#### Misc. non venom associated plots #### 
library(Seurat)
library(tidyverse)
library(viridis)
library(cowplot)
setwd('/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/ERK_UPR')

misc_genes <- read.csv('sc_nontoxin_genes.csv', header=T)
venom2 <- readRDS("venomClusters_06.23.rds")

p1 <- FeaturePlot(venom2, c(misc_genes$gtf_id[2]), order = TRUE) + ggtitle(misc_genes$gene[2])
p2 <- FeaturePlot(venom2, c(misc_genes$gtf_id[3]), order = TRUE) + ggtitle(misc_genes$gene[3])
p3 <- FeaturePlot(venom2, c(misc_genes$gtf_id[4]), order = TRUE) + ggtitle(misc_genes$gene[4])
p4 <- FeaturePlot(venom2, c(misc_genes$gtf_id[5]), order = TRUE) + ggtitle(misc_genes$gene[5])
p5 <- FeaturePlot(venom2, c(misc_genes$gtf_id[6]), order = TRUE) + ggtitle(misc_genes$gene[6])
p6 <- FeaturePlot(venom2, c(misc_genes$gtf_id[7]), order = TRUE) + ggtitle(misc_genes$gene[7])
p7 <- FeaturePlot(venom2, c(misc_genes$gtf_id[8]), order = TRUE) + ggtitle(misc_genes$gene[8])
p8 <- FeaturePlot(venom2, c(misc_genes$gtf_id[9]), order = TRUE) + ggtitle(misc_genes$gene[9])
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow=3)


#### Plot Stochastic gene choice ####

venom2 <- readRDS("venomClusters_06.23.rds")

subset.matrix <- venom2@assays$RNA[c('fgenesh-scaffold-mi2-venom-gene-3.8',"fgenesh-scaffold-mi2-venom-gene-3.9",
                                     'fgenesh-scaffold-mi7-venom-gene-1.5','fgenesh-scaffold-mi7-venom-gene-1.2',
                                     'fgenesh-scaffold-mi7-venom-gene-1.4', 'fgenesh-scaffold-mi1-venom-gene-2.1',
                                     'fgenesh-scaffold-mi1-venom-gene-2.2', 'fgenesh-scaffold-mi1-venom-gene-2.3',
                                     'fgenesh-scaffold-mi1-venom-gene-2.4', 'fgenesh-scaffold-mi1-venom-gene-2.5',
                                     'fgenesh-scaffold-mi1-venom-gene-2.6','fgenesh-scaffold-mi1-venom-gene-2.7',
                                     'fgenesh-scaffold-mi1-venom-gene-2.8','fgenesh-scaffold-mi1-venom-gene-2.9',
                                     'fgenesh-scaffold-mi1-venom-gene-2.10'),]
subset.matrix <- as.matrix(subset.matrix)
rownames(subset.matrix) <- c('SVSP8','SVSP9','PLA2A1','PLA2B1','PLA2C1','SVMP1','SVMP2',
                              'SVMP3','SVMP4','SVMP5','SVMP6','SVMP7','SVMP8','SVMP9','SVMP10')
subset.matrix <- t(subset.matrix)

SVMP_mat <- subset.matrix[,grep('SVMP', colnames(subset.matrix))]
hist(rowSums(SVMP_mat > 0))
colSums(SVMP_mat>0)/nrow(SVMP_mat) # expression frequencies
