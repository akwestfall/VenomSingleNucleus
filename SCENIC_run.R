library(tidyverse)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCENIC)
library(Seurat)

setwd('/Users/sidgopalan/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/SCENIC')
RNA_mat <- read.csv('/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/GRN_stuff/SIGNET/RNA_mat_unimputed.csv',row.names = 1)
colnames(RNA_mat) <- gsub('\\.','-',colnames(RNA_mat))
RNA_mat <- as.matrix(RNA_mat)


#dir.create("SCENIC_Cvv_VG")
setwd("SCENIC_Cvv_VG/")

# Get per cell information metadata
seurat_Obj <- readRDS('/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/ERK_UPR/venomClusters_06.23.rds')
cellInfo <- seurat_Obj@meta.data
colnames(cellInfo)[6] <- 'CellType' # Call Seurat clusters celltype
cellInfo$CellType <- as.character(cellInfo$CellType)
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("0"="forestgreen", 
                           "1"="darkorange", 
                           "2"="magenta4", 
                           "3"="hotpink"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]

# Initialize SCENIC
org <- "hgnc" # or hgnc, or dmel
myDatasetTitle <- "SCENIC example Cvv_VG" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir="/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/GRN_stuff/SCENIC", dbs=dbs, datasetTitle = myDatasetTitle, nCores=4)

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars
saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(colVars, file="int/colVars.Rds")

### Gene filter/selection ###
# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(RNA_mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(RNA_mat),
                           minSamples=ncol(RNA_mat)*.005) # default 0.01

exprMat_filtered <- RNA_mat[genesKept, ]
dim(exprMat_filtered)

### Run Correlation ###
runCorrelation(exprMat_filtered, scenicOptions)

### Run GENIE3 ###
scenicOptions@settings$seed <- 1234 # Random forest, reproduce

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions) # Takes long time

### Run SCENIC ###
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # Takes long time (hrs)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status

###Binarize the network activity (regulon on/off) ###
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered)
savedSelections <- shiny::runApp(aucellApp)
# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Run binarization #
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)



#### Clustering / dimensionality reduction on the regulon activity ####
# Run t-SNE with different settings:
nPcs <- c(5,15,50)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

#cellInfo$CellType <- as.character(cellInfo$CellType)
#scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
#(cellInfo, file="int/cellInfo.Rds")
#par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName = "CellType", cex=.5)



VG_seuratObj <- readRDS('/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/ERK_UPR/venomClusters_06.23.rds')
scenicOptions <- readRDS('int/scenicOptions.Rds')
regulons <- loadInt(scenicOptions, "regulons")
dr_coords <- Embeddings(VG_seuratObj, reduction="umap")
tfs <- c("CTCF","ATF4","RORB", "XBP1")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, exprMat_filtered, cellsAUC=aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("CREB3L2","EP300","HCFC1","NFKB1")],], plots = "AUC")


tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("CREB3L2","EP300","HCFC1","NFKB1")],], plots="Expression", exprCols = c("lightgrey","blue"))
AUCell::AUCell_plotTSNE(dr_coords, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))["FOXP1"],], plots="Expression", exprCols = c("lightgrey","blue"))

FeaturePlot(VG_seuratObj,features = 'maker-scaffold-ma6-augustus-gene-195.2', order = T)


aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("CTCF", "ATF4", "RORB","XBP1")],]
aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("CTCF","ATF4")],]

# Write 96 regulon AUCs (activity scores) to file
regulonAUC_mat <- as.matrix(getAUC(aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))]))
write.csv(regulonAUC_mat, file = '/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/GRN_stuff/SCENIC/regulon_analyses/regulonAUC_mat.csv')



#### Regulators for known cell types or clusters ####
library(Seurat)
library(pheatmap)
scenicOptions <- readRDS('int/scenicOptions.Rds')

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

#ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
pheatmap(regulonActivity_byCellType_Scaled, border_color = "NA", legend_labels = c('Regulon activity'), fontsize = 9, cluster_cols = F)
pheatmap(regulonActivity_byCellType_Scaled, border_color = "NA", legend_labels = c('Regulon activity'), fontsize = 9, cluster_cols = F, cluster_rows = F)

#### Combine heatmap with annotation table ####
regulonActivity_byCellType_Scaled <- as.data.frame(regulonActivity_byCellType_Scaled)
regulonActivity_byCellType_Scaled$shortname <- rownames(regulonActivity_byCellType_Scaled)

regulonActivity_byCellType_Scaled$shortname <- gsub('_extended','',regulonActivity_byCellType_Scaled$shortname)
regulonActivity_byCellType_Scaled$shortname <- gsub(' \\([0-9]*g\\)','',regulonActivity_byCellType_Scaled$shortname)


library(tidyverse)
library(readxl)
library(ggpubr)

regulon_classification <- read_excel('/Volumes/SeagatePortableDrive/venom_scRNAseq/Figures/GRN_stuff/SCENIC/regulon_analyses/96_regulon_classification_reordered_test.xlsx')
ngenes=96 # Global variable, doesn't change

# Modify metadata slightly
regulon_classification <- regulon_classification[,c(1:5)]
colnames(regulon_classification) <- c('regulon','chrm_mod','prev_imp_VG','UPR','ERK')

merge(regulonActivity_byCellType_Scaled,regulon_classification,by=c("shortname","regulon"))
merged_data <- merge(regulonActivity_byCellType_Scaled, regulon_classification, by.x=c("shortname"), by.y=c("regulon"))

merged_data <- merged_data[match(regulon_classification$regulon, merged_data$shortname),]

rownames(merged_data) <- NULL
merged_data <- merged_data %>% column_to_rownames(var = 'shortname')
merged_data <- merged_data %>% select(-(1:4)) %>% type_convert(na='NA')

ggballoonplot(merged_data, size.range = c(0,2.5), shape=23, ggtheme = theme_bw(), 
              fill = c(rep("#E69F00", nrow(merged_data)),
                       rep("#56B4E9", nrow(merged_data)),
                       rep("#D55E00", nrow(merged_data)),
                       rep("#CC79A7", nrow(merged_data)))) + 
  theme(legend.position = "none")

regulonActivity_byCellType_Scaled_2 <- regulonActivity_byCellType_Scaled[match(regulon_classification$regulon, regulonActivity_byCellType_Scaled$shortname),]

pheatmap(regulonActivity_byCellType_Scaled_2[,c(1:4)], border_color = "NA", legend_labels = c('Regulon activity'), fontsize = 9, cluster_cols = F, cluster_rows = F)



# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "0")
plotRSS_oneSet(rss, setName = "1")
plotRSS_oneSet(rss, setName = "2")
plotRSS_oneSet(rss, setName = "3")


