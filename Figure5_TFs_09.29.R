### Plot Grids
col1_row1 <- plot_grid(erk.ven.enrich, upr.ven.enrich, labels = c("A", "B"), nrow = 1)
col1_row2 <- plot_grid(erk.enrich, upr.enrich,
                       labels = c("D", "E"))
col1_row3 <- plot_grid(umap.elk1, NULL, umap.atf1, umap.xbp1, NULL, umap.atf6,
                       umap.tbx3, NULL, umap.srf, umap.ern1, NULL, umap.eif2ak3,
                       rel_widths = c(1, -0.1, 1, 1, -0.1, 1),
                       nrow = 2)

col1 <- plot_grid(col1_row1, col1_row2, col1_row3, nrow = 3, rel_heights = c(1, 1.25, 1))
col2 <- plot_grid(erk.upr.enrich, htmp, nrow = 2, rel_heights = c(2, 3), labels = c("C", "F"))

fig5 <- plot_grid(col1, col2, ncol = 2, rel_widths = c(2, 1))

ggsave("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/ms/figurePanels/Fig5_TFs.pdf",
       plot = fig5, width = 172, height = 134, units = "mm", dpi = 300)

### ERK and UPR Enrichment Correlations ######
erk.ven.enrich <- ggplot(venom2@meta.data, aes(x = ERK, y = VenomGenes)) + 
  geom_hex(bins = 50) +
  geom_smooth(method = 'lm', se = F, color = 'red', size = 0.5) +
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 1, size = 2.14) +
  guides(fill = guide_colorbar(barheight = 2.5, barwidth = 0.5)) +
  xlab('ERK enrichment') +
  ylab('venom gene enrichment') +
  ggtitle("ERK and venom gene enrichment") +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_minimal() +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(5,5,0,5),
        legend.margin = margin(0,5,5,5),
        legend.title = element_text(face = "bold"),
        legend.position = "right")

upr.ven.enrich <- ggplot(venom2@meta.data, aes(x = UPR, y = VenomGenes)) + 
  geom_hex(bins = 50) +
  geom_smooth(method = 'lm', se = F, color = 'red', size = 0.5) +
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 1, size = 2.14) +
  guides(fill = guide_colorbar(barheight = 2.5, barwidth = 0.5)) +
  xlab('UPR enrichment') +
  ylab('venom gene enrichment') +
  ggtitle("UPR and venom gene enrichment") +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_minimal() +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(5,5,0,5),
        legend.margin = margin(0,5,5,5),
        legend.title = element_text(face = "bold"),
        legend.position = "right")

erk.upr.enrich <- ggplot(venom2@meta.data, aes(x = ERK, 
                                               y = UPR)) + 
  ggtitle("UPR and ERK enrichment correlation") +
  geom_point(color = "black", size = 0.4, fill = "gray50", shape = 21, alpha = 0.75) + 
  geom_smooth(method = 'lm', se = F, color = 'red', size = 0.5) +
  ggpubr::stat_cor(method = "pearson", label.x = 400, label.y = -500, size = 2.14) +
  xlab('ERK enrichment') +
  ylab('UPR enrichment') +
  theme_minimal() +
  theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
        plot.margin = margin(5, 10, 5, 5),
        text = element_text(size = 6),
        legend.position = "none")

### Pathway Enrichment UMAP Plots ######
umap <- venom2@reductions$umap@cell.embeddings[,1:2]
umap <- as.data.frame(umap)
umap$UPR <- venom2@meta.data$UPR2
umap$ERK <- venom2@meta.data$ERK2

upr.enrich <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = UPR)) + geom_point(size = 0.5) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("UPR Enrichment") + labs(color = "enrichment") +
  theme_void() + guides(color = guide_colorbar(barwidth = 2.5, barheight = 0.5)) +
  scale_color_distiller(palette = "Reds", direction = 1, labels = c(-1000, 2000, 4000), breaks = c(-1000, 2000, 4000)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(angle = 90),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom") + coord_fixed(5/6)

erk.enrich <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = ERK)) + geom_point(size = 0.5) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ERK Enrichment") + labs(color = "enrichment") +
  theme_void() + guides(color = guide_colorbar(barwidth = 2.5, barheight = 0.5)) +
  scale_color_distiller(palette = "Reds", direction = 1, labels = c(-2000, 0, 2000), breaks = c(-2000, 0, 2000)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(angle = 90),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom") + coord_fixed(5/6)

panel1 <- plot_grid(erk.ven.enrich, upr.ven.enrich, erk.enrich, upr.enrich,
                    labels = c("A", "B", "C", "D"))
panel1

### UPR Gene Expression ######
umap.xbp1 <- FeaturePlot(venom2, features = "maker-scaffold-mi8-augustus-gene-35.54", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("XBP1") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(5, 0, 0, 5),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.xbp1$layers[[1]]$aes_params$size <- 0.5

umap.atf6 <- FeaturePlot(venom2, features = "augustus-masked-scaffold-ma3-processed-gene-300.3", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ATF6") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 0, 0),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.atf6$layers[[1]]$aes_params$size <- 0.5

umap.ern1 <- FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-362.6", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ERN1") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(0, 0, 5, 5),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.ern1$layers[[1]]$aes_params$size <- 0.5

umap.eif2ak3 <- FeaturePlot(venom2, features = "maker-scaffold-Z-augustus-gene-255.8", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("EIF2AK3") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(0, 5, 5, 0),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.eif2ak3$layers[[1]]$aes_params$size <- 0.5


### ERK Gene Enrichment ######
umap.elk1 <- FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-434.11", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ELK1") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(5, 0, 0, 5),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.elk1$layers[[1]]$aes_params$size <- 0.5

umap.atf1 <- FeaturePlot(venom2, features = "maker-scaffold-ma2-augustus-gene-372.8", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("ATF1") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 0, 0),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.atf1$layers[[1]]$aes_params$size <- 0.5

umap.tbx3 <- FeaturePlot(venom2, features = "augustus-masked-scaffold-ma1-processed-gene-750.0", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("TBX3") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(0, 0, 5, 5),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.tbx3$layers[[1]]$aes_params$size <- 0.5

umap.srf <- FeaturePlot(venom2, features = "maker-scaffold-ma1-augustus-gene-935.8", order = TRUE) +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("SRF") + coord_fixed(5/6) +
  theme_void() + guides(color = guide_colorbar(barwidth = 0.25, barheight = 2)) +
  scale_color_gradientn(colours=c('lightgrey', "blue")) +
  theme(text = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(0, 5, 5, 0),
        legend.title = element_text(face = "bold"),
        legend.position = "right")
umap.srf$layers[[1]]$aes_params$size <- 0.5

plot_grid(umap.elk1, umap.atf1, umap.xbp1, umap.atf6,
          umap.tbx3, umap.srf, umap.ern1, umap.eif2ak3,
          nrow = 2)

### Coexpression Heatmap ######
## data from venom grant figures
data$toxin_family <- factor(data$toxin_family, levels = c("myo.", "SVSP", "SVMP", "PLA2", "other"))
htmp <- ggplot(data, aes(x = tf, y = toxin, fill = rho)) + geom_tile(width = 1.1, height = 1.1) +
  theme_minimal() + coord_cartesian(expand = FALSE) +
  ggtitle("TF coexpression with venom genes") +
  scale_fill_gradientn(colours = colorRampPalette(c("#3182BD", "black", "#DE2D26"), bias = 2)(20), na.value = "grey50") +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 2.5)) +
  facet_grid(rows = vars(toxin_family), cols = vars(tf_group), scales = "free", space = "free", shrink = FALSE) +
  theme(text = element_text(size = 6),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95),
        legend.title = element_text(face = "bold"),
        legend.position = "right",
        legend.margin = margin(c(0, 0, 0, -5)),
        panel.spacing = unit(0.05, "lines"),
        strip.text=element_text(margin=margin(t = 2, r = 2, b = 2, l = 0), face = "bold", size = 7),
        panel.grid = element_blank(),
        panel.border = element_blank())
htmp
