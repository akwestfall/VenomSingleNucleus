library(tidyverse)
library(DESeq2)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/Correlations/SSG_correlations/") # change for Aundrea

#### Multitissue RNAseq with STAR####
options(scipen = 999)


countdata <- read.table("/Users/aundreawestfall/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/heatmap_piechart/CrotalusViridis_multitissue_featurecounts.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]

colnames(countdata) <- gsub("\\..multitissue_STAR_outs.", "", colnames(countdata))
colnames(countdata) <- gsub("\\_Aligned.sortedByCoord.out.bam", "", colnames(countdata))

tissues <- c("VG_1DPE","VG_3DPE","VG_U", "AVG","blood","brain","kidney","liver","lung","ovaries","pancreas","rictal_gland",
             "shaker_muscle","skin","spleen","stomach","testes","tongue")

countdata2 <- do.call(cbind, lapply(countdata, function(x) data.frame(x,x)))
rownames(countdata2) <- row.names(countdata)

# gross
t1 <- (rowMeans(countdata2[,grep(tissues[1],colnames(countdata2))]))
t2 <- (rowMeans(countdata2[,grep(tissues[2],colnames(countdata2))]))
t3 <- (rowMeans(countdata2[,grep(tissues[3],colnames(countdata2))]))
t4 <- (rowMeans(countdata2[,grep(tissues[4],colnames(countdata2))]))
t5 <- (rowMeans(countdata2[,grep(tissues[5],colnames(countdata2))]))
t6 <- (rowMeans(countdata2[,grep(tissues[6],colnames(countdata2))]))
t7 <- (rowMeans(countdata2[,grep(tissues[7],colnames(countdata2))]))
t8 <- (rowMeans(countdata2[,grep(tissues[8],colnames(countdata2))]))
t9 <- (rowMeans(countdata2[,grep(tissues[9],colnames(countdata2))]))
t10 <- (rowMeans(countdata2[,grep(tissues[10],colnames(countdata2))]))
t11 <- (rowMeans(countdata2[,grep(tissues[11],colnames(countdata2))]))
t12 <- (rowMeans(countdata2[,grep(tissues[12],colnames(countdata2))]))
t13 <- (rowMeans(countdata2[,grep(tissues[13],colnames(countdata2))]))
t14 <- (rowMeans(countdata2[,grep(tissues[14],colnames(countdata2))]))
t15 <- (rowMeans(countdata2[,grep(tissues[15],colnames(countdata2))]))
t16 <- (rowMeans(countdata2[,grep(tissues[16],colnames(countdata2))]))
t17 <- (rowMeans(countdata2[,grep(tissues[17],colnames(countdata2))]))
t18 <- (rowMeans(countdata2[,grep(tissues[18],colnames(countdata2))]))
collapsed_tissue_counts <- cbind(t1, t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18)
rm(t1, t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18, countdata2)
colnames(collapsed_tissue_counts) <- tissues
collapsed_tissue_counts <- round(collapsed_tissue_counts)

countdata_b <- read.table("/Users/aundreawestfall/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/heatmap_piechart/Bulk_RNAseq_featurecounts.txt", header=TRUE, row.names=1)
countdata_b <- countdata_b[ ,6:ncol(countdata_b)]

colnames(countdata_b) <- gsub("\\..Bulk_STAR_outs.", "", colnames(countdata_b))
colnames(countdata_b) <- gsub("\\_Aligned.sortedByCoord.out.bam", "", colnames(countdata_b))
colnames(countdata_b) <- gsub("\\..*$", "", colnames(countdata_b))

collapsed_tissue_counts <- cbind(countdata_b,collapsed_tissue_counts)

#DESeq2 stuff#
sampleTable <- as.data.frame(colnames(collapsed_tissue_counts))
colnames(sampleTable) <- "sample"
sampleTable$tissue <- ifelse(grepl("VG", sampleTable$sample),"Venom","Body")
dds <- DESeq2::DESeqDataSetFromMatrix(collapsed_tissue_counts,
                                      sampleTable,
                                      ~tissue)

dds <- DESeq2::DESeq(dds)

# heatmap
normalized_counts <- t(counts(dds, normalized=TRUE))
saveRDS(object = dds, "~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/heatmap_piechart/bulkRNA.dds.rds")
dds <- readRDS("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/heatmap_piechart/bulkRNA.dds.rds")
venom_names <- read.csv('/Users/aundreawestfall/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/heatmap_piechart/venomTxNames.csv', header=T)

normalized_counts <- normalized_counts[,grep(paste0(venom_names$altname,collapse = '|'), colnames(normalized_counts))]

venom_names$altname2 <- NA
for (i in 1:nrow(venom_names)){
  venom_names$altname2[i] <- grep(venom_names$altname[i], colnames(normalized_counts), value=T)
}

normalized_counts <- normalized_counts[, venom_names$altname2]
colnames(normalized_counts) <- venom_names$toxin

colnames(normalized_counts) <- gsub("_", "", gsub("ADAM28_", "ADAM28 ", colnames(normalized_counts)))

rownames(normalized_counts) <- c('Accessory (Bulk)','LVG (Bulk)','Liver (Bulk)','Pancreas (Bulk)', 
                                 'RVG (Bulk)', 'Skin (Bulk)', 'VG (1DPE)','VG (3DPE)', 'VG (Unext.)', 'Accessory VG',
                                 'Blood', 'Brain', 'Kidney', 'Liver', 'Lung', 'Ovaries',
                                 'Pancreas', 'Rictal gland', 'Shaker muscle', 'Skin', 
                                 'Spleen', 'Stomach', 'Testes', 'Tongue') # rename tissues

normalized_counts <- normalized_counts[, c(1:28,44:49,29:33,50:51,34:43)]

# Make heatmap
ComplexHeatmap::pheatmap(log(t(normalized_counts)+1), cluster_cols = F, cluster_rows  = F, 
                         color=colorRampPalette(c("black",'#D65A14',"#FDFFCD"))(50),
                         gaps_row = c(11,22, 26),
                         gaps_col = c(3,4),
                         border_color = 'NA', 
                         angle_col = c("315"))

mat <- log(t(normalized_counts[c(5,7:24),])+1)
breaks <- seq(min(mat), max(mat), length = 3)
cols <- colorRamp2(breaks, c("black",'#D65A14',"#FDFFCD"), space = "LAB")

row_gaps <- grep("myotoxin|SVSP11|SVMP11|PLA2K|BPP|CTL6|LAAO3|CRISP4", row.names(mat))
col_gaps <- c(1, 4, 5)

anno <- data.frame(gene = row.names(mat), family = NA)
row.names(anno) <- anno$gene

htmp <- pheatmap::pheatmap(mat, color = colorRampPalette(c("black",'#D65A14',"#FDFFCD"))(50),
                           fontsize = 6,
                           angle_col = 90,
                           cluster_rows = FALSE, gaps_row = row_gaps, 
                           cluster_cols = FALSE, gaps_col = col_gaps,
                           border_color = NA, scale = "none",
                           main = "A. Bulk tissue RNA venom gene expression"
)
grid.arrange(htmp[[4]])



# Make piechart
library(tidyverse)
library(plotrix)

setwd('/Users/sidgopalan/Documents/Projects/Venom_scRNA/Figures')

fractions <- read.csv('VenomFractions.csv')

family_summed <- fractions %>% 
  group_by(Family) %>% 
  summarise(Freq = sum(Fraction)) %>% 
  arrange(desc(Freq)) %>% 
  arrange(Family=="Other") %>% 
  filter(Freq>1) %>% 
  mutate(Family=factor(Family, levels=Family))

fractions[fractions$Family=='Other',]$Fraction <- fractions[fractions$Family=='Other',]$Fraction 

pie3D(family_summed$Freq,
      col = hcl.colors(length(family_summed$Freq), "Spectral"),
      labels=family_summed$Family,
      explode=0.05,
      border='white')

#### Number of cells vs bulk expression ####
library(cowplot)

setwd('/Volumes/SeagatePortableDrive/venom_scRNAseq') # Change for Aundrea

VG_seurat <- readRDS('Figures/ERK_UPR/seuratClusters_06.23.rds')
VG_seurat <- venom2

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

venom_genes <- read.csv('/Users/aundreawestfall/Desktop/snRNAseq_viridis/res.toxin_transcripts.csv', header=T)
venom_genes <- read.csv('/Users/aundreawestfall/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/Other_files-figs/heatmap_piechart/res.toxin_transcripts.csv', header=T)

expression_set <- cbind(PrctCellExpringGene(VG_seurat,genes=venom_genes$transcript, group.by = "all"),venom_genes)
expression_set <- expression_set[!is.na(expression_set$Cell_proportion),]
rownames(expression_set) <- NULL

VG_counts <- as.data.frame(normalized_counts[2,])
VG_counts <- tibble::rownames_to_column(VG_counts, 'gene')
VG_counts[1,1] <- "myotoxin3"

expression_set$bulk_exp <- NA
expression_set$bulk_exp[1] <- VG_counts[1,2] # Myotoxin
expression_set$bulk_exp[2:10] <- VG_counts[c(2:10),2] # SVSP 1 - 9
expression_set$bulk_exp[11:20] <- VG_counts[c(13:22),2] # SVMP 1 - 10
expression_set$bulk_exp[21] <- VG_counts[c(25),2] # PLA2B1
expression_set$bulk_exp[22] <- VG_counts[c(26),2] # PLA2C1
expression_set$bulk_exp[23] <- VG_counts[c(24),2] # PLA2A1
expression_set$bulk_exp[24] <- VG_counts[c(28),2] # BPP
expression_set$bulk_exp[25:26] <- VG_counts[c(37,35),2] # LAAO1,3
expression_set$bulk_exp[27] <- VG_counts[c(38),2] # CRISP1
expression_set$bulk_exp[28] <- VG_counts[c(42),2] # ohanin
expression_set$bulk_exp[29:31] <- VG_counts[c(43:45),2] # EXOs
expression_set$bulk_exp[32] <- VG_counts[c(46),2] # VEGF1
expression_set$bulk_exp[33] <- VG_counts[c(30),2] # CTL2
expression_set$bulk_exp[34:35] <- VG_counts[c(48:49),2] # vQC1,2
expression_set$bulk_exp <- as.numeric(expression_set$bulk_exp)

library(ggpubr)
library(ggrepel)
# Reclassify all others (not PLA2, not SVSP, not SVMP)
expression_set$family2 <- ifelse(expression_set$family!='SVSP'&expression_set$family!='SVMP'&expression_set$family!='PLA2', 'Other' , expression_set$family)
p1 <- expression_set %>% ggplot(aes(x=log(bulk_exp+1),y=Cell_proportion)) + 
  geom_point(aes(color=family2)) + 
  theme_bw() +
  scale_color_brewer(palette="Dark2") +
  stat_smooth(method="glm", se=F, method.args = list(family=binomial), color="red") +
  stat_cor(method = "pearson", label.x = 2.5, label.y = 0.625, size = 2.14) +
  xlab('log(normalized counts + 1)') +
  ylab('proportion of expressing cells') +
  ggtitle("E. Bulk RNAseq expression versus\nproportion of cells with gene") +
  #geom_text_repel(aes(label=ifelse(Cell_proportion>0.90,as.character(toxin),'')),size = 2.14, segment.size = 0.25) +
  geom_text_repel(aes(label = toxin),size = 2.14, segment.size = 0.25) +
  theme(legend.position="top",
        text = element_text(size = 6),
        plot.title = element_text(size = 6, face = "bold"),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.35, 'cm'),
        legend.margin = margin(c(0,0,0,0)),
        legend.title = element_blank())

#### Backcalculated single cell expression vs bulk expression ####
RNA_mat <- as.matrix(VG_seurat@assays$RNA[,])

expression_set$backcalculated_exp <- NA
for (i in 1:nrow(expression_set)){
  expression_set$backcalculated_exp[i] <- mean(RNA_mat[expression_set$transcript[i],]/colSums(RNA_mat),na.rm = T)
}


p2 <- expression_set %>% ggplot(aes(x=log(bulk_exp+1),y=backcalculated_exp)) + 
  geom_point(aes(color=family2)) + 
  theme_bw() +
  scale_color_brewer(palette="Dark2") +
  stat_smooth(method="glm", se=F, method.args = list(family=binomial), color="red") +
  stat_cor(method = "pearson", label.x = 2.5, label.y = 0.002, size = 2.14) +
  xlab('log(normalized counts + 1)') +
  ylab('averaged per cell expression magnitude') +
  ggtitle("F. Bulk RNAseq expression versus\nmagnitude of expression in cell") +
  geom_text_repel(aes(label=ifelse(backcalculated_exp>0.0025,as.character(toxin),'')), size = 2.14, segment.size = 0.25) +
  theme(legend.position="top",
        text = element_text(size = 6),
        plot.title = element_text(size = 6, face = "bold"),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(.35, 'cm'),
        legend.margin = margin(c(0,0,0,0)),
        legend.title = element_blank())

p3 <- plot_grid(p1,p2, nrow=1)
ggsave("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/ms/figurePanels/Fig2_bulk.v.sn.pdf", plot = p3, width = 172, height = 85, units = "mm", dpi = 300)






