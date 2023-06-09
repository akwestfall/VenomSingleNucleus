library(tidyverse)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/Correlations/SSG_correlations/") # change for Aundrea
options(scipen=0)

VG_cor <- read.csv('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/Correlations/venomGenes_VG_correlation.csv')
VG_cor <- subset(VG_cor, select=-c(X))

# remove myotoxins
VG_cor <- VG_cor[!grepl("myo", VG_cor$gene1),]
VG_cor <- VG_cor[!grepl("myo", VG_cor$gene2),]

# Remove duplicates
VG_cor <- distinct(VG_cor,rho,p.value,FDR,.keep_all = T)

# Find missing SVSP8-9 correlation
obj <- readRDS('varcor_06.22.RDS')
test <- as.data.frame(obj[which(obj$gene2=='fgenesh-scaffold-mi2-venom-gene-3.9'),])
test[6,]$gene1 <- "SVSP8"
test[6,]$gene2 <- "SVSP9"
VG_cor <- rbind(VG_cor,test[6,])
VG_cor <- VG_cor %>% arrange(desc(rho))
rm(test)


## Add chromatin loop info
# All PLA2s in loop, SVSP8-9 loop
VG_cor$loop <- NA
VG_cor[grepl('PLA2', VG_cor$gene1) & grepl('PLA2', VG_cor$gene2), 'loop'] <- 1
VG_cor[grepl('SVSP8', VG_cor$gene1) & grepl('SVSP9', VG_cor$gene2), 'loop'] <- 1
VG_cor[which(is.na(VG_cor$loop)),'loop'] <- 0

# rename genes
VG_cor$gene1 <- sub("(.*)(\\d+)$", "\\1_\\2", VG_cor$gene1)
VG_cor$gene2 <- sub("(.*)(\\d+)$", "\\1_\\2", VG_cor$gene2)
rep_str <- c('PLA2 A_1'='PLA2_A1',
             'PLA2 B_1'='PLA2_B1',
             'PLA2 C_1'='PLA2_C1',
             'PLA2 gIIE'='PLA2_gIIE')
VG_cor$gene1 <- str_replace_all(VG_cor$gene1, rep_str)
VG_cor$gene2 <- str_replace_all(VG_cor$gene2, rep_str)
rm(rep_str)

# add same family (same_fam) column. Are the two genes in the same array?
rownames(VG_cor) <- NULL # reset row indices

VG_cor$same_fam <- NA

VG_cor[as.numeric(rownames(with(VG_cor, VG_cor[ grepl( 'SVMP', gene1) & grepl( 'SVMP', gene2), ]))), 'same_fam'] <- 1
VG_cor[as.numeric(rownames(with(VG_cor, VG_cor[ grepl( 'PLA2', gene1) & grepl( 'PLA2', gene2), ]))), 'same_fam'] <- 1
VG_cor[as.numeric(rownames(with(VG_cor, VG_cor[ grepl( 'SVSP', gene1) & grepl( 'SVSP', gene2), ]))), 'same_fam'] <- 1
VG_cor[which(is.na(VG_cor$same_fam)),'same_fam'] <- 0


#### Incorporate linear distance between genes ####
gtf <- rtracklayer::import('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/data/Correlations/SSG_correlations/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf') # can use regular gtf too
gtf=as.data.frame(gtf)
gtf <- gtf[gtf$type=='gene',]

# Select only important families, choose nbis IDs
gtf_VG <- gtf[grep('SVMP|SVSP|PLA2', gtf$gene_id),]
gtf_VG <- gtf_VG[grep('nbis', gtf_VG$gene_id),]

gtf_VG <- gtf_VG %>% select(seqnames,start,end,Name)
gtf_VG$Name <- str_replace(gtf_VG$Name,'PLA2gIIE','PLA2_gIIE')


# calculate start-end and end-start, take whatever is smaller (intergenic space, which excludes gene lengths)
# crappy readability..
VG_cor$dist <- NA
for (i in 1:nrow(VG_cor)){
  if(VG_cor$same_fam[i]==1){
    VG_cor$dist[i] <- min(c(abs(gtf_VG[(grep(VG_cor$gene1[i], gtf_VG$Name)),]$start-gtf_VG[(grep(VG_cor$gene2[i], gtf_VG$Name)),]$end),
                            abs(gtf_VG[(grep(VG_cor$gene1[i], gtf_VG$Name)),]$end-gtf_VG[(grep(VG_cor$gene2[i], gtf_VG$Name)),]$start)
    ))
  }
}

# Multiple linear regression test
#model1 <- lm(rho~loop+dist,data=VG_cor_samefam)
#summary(model1)

#### Random sampling from microchromosomes ####
set.seed(1234)

# Remove duplicates
obj <- as.data.frame(obj)
obj <- distinct(obj,rho,p.value,FDR,.keep_all = T)

# Subset all rows with genes on the same chrs, takes a minute
a = obj[grepl('mi1-', obj$gene1) & grepl('mi1-', obj$gene2),]
b = obj[grepl('mi2-', obj$gene1) & grepl('mi2-', obj$gene2),]
c = obj[grepl('mi7-', obj$gene1) & grepl('mi7-', obj$gene2),]

n = 1000 # sample size
sample_cor <- rbind(a[sample(nrow(a), n),],
                b[sample(nrow(b), n),],
                c[sample(nrow(c), n),])

rownames(sample_cor) <- NULL # reset row indices
rm(list=c('a','b','c'))

rownames(gtf) <- NULL # reset row indices

gtf$gene_id <- gsub("_", "-", gtf$gene_id) # replace underscores in gene_id

# calculate start-end and end-start, take whatever is smaller (intergenic space, which excludes gene lengths)
# bad readability..Takes about 6 mins. In the future, use gene midpoints
sample_cor$dist <- NA
system.time(for (i in 1:nrow(sample_cor)){
    sample_cor$dist[i] <- min(abs(gtf[grep(sample_cor$gene1[i], gtf$gene_id),]$start-gtf[grep(sample_cor$gene2[i], gtf$gene_id),]$end),
                             abs(gtf[grep(sample_cor$gene1[i], gtf$gene_id),]$end-gtf[grep(sample_cor$gene2[i], gtf$gene_id),]$start))
})



# Add gene family names to dfs
sample_cor$family <- NA
sample_cor[grep('mi1', sample_cor$gene1),'family'] <- 'mi1'
sample_cor[grep('mi2', sample_cor$gene1),'family'] <- 'mi2'
sample_cor[grep('mi7', sample_cor$gene1),'family'] <- 'mi7'

VG_cor_samefam <- VG_cor[which(VG_cor$same_fam==1),]
VG_cor_samefam$family <- NA
VG_cor_samefam[grep('SVMP', VG_cor_samefam$gene1),'family'] <- 'mi1'
VG_cor_samefam[grep('SVSP', VG_cor_samefam$gene1),'family'] <- 'mi2'
VG_cor_samefam[grep('PLA2', VG_cor_samefam$gene1),'family'] <- 'mi7'

# combine VG and background data
combined_data <- bind_rows(VG_cor_samefam,sample_cor,.id = 'id')
combined_data$id <- ifelse(combined_data$id=='1','venom','background')

# correct distance by chr lengths
combined_data <- combined_data %>% mutate(cor_dist=case_when(family=='mi1'~dist/22521304,
                                                            family=='mi2'~dist/19978503,
                                                            family=='mi7'~dist/12380205))

#### Begin plotting ####
### Plot background vs venom genes scatterplot ###
f3.bg_v_vg <- combined_data %>% ggplot(aes(x=-log(cor_dist),y=rho,color=id)) +
  geom_point(size=2) +
  ylab(expression(rho)) +
  xlab('-log(Scaled intergenic distance)') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

### Plot background vs venom genes boxplots ###
f3.bg_boxes <- combined_data %>% ggplot(aes(x=id,y=rho,color=id)) + geom_boxplot() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

f3.family_rho <- combined_data %>% 
  filter(id=='venom') %>%  
  mutate(family=case_when(family=='mi1'~'SVMP',
                            family=='mi2'~'SVSP',
                            family=='mi7'~'PLA2')) %>% 
  ggplot(aes(x=(cor_dist),y=rho,color=family)) +
  geom_point(size=2) +
  ylab(expression(rho)) +
  xlab('Scaled intergenic distance') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Non parametric test
wilcox.test(rho~id, data = combined_data)
