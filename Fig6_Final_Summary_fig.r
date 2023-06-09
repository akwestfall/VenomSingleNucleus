# (0) ENVIRONMENT
library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(scales)
library(data.table)
library(ggnewscale)
library(tidyverse)
library(ggh4x)
library(matrixStats)

'%notin%' <- function(x,y)!('%in%'(x,y))

### LM correlation model ######
### DEFUNCT

all <- read.table("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_Venom_scRNA/FromRich_Modeling/fromSid_Modeling_LargeTFs_8Dec22/ByParalog/PARALOG_Alldata_RESULTS_ModelFit_08Dec22.txt",
                  header = T, row.names = NULL)

sig.tfs <- all %>% pivot_longer(!row.names, names_sep = "_", names_to = c("Paralog", "var"), values_to = "value") %>%
  pivot_wider(names_from = var, values_from = value) %>%
  filter(PValue <= 0.05) %>%
  mutate(Family = substr(Paralog, start = 1, stop = 4)) %>%
  rename("TFs" = "row.names", "value" = "Coefficient") %>%
  mutate(var = "Beta") %>%
  select(c("TFs", "Family", "Paralog", "var", "value", "Importance")) 

### GENIE3 results ######
gen.in <- read.table("/Users/aundreawestfall/Desktop/Venom_snRNA/summaryFigData/GENIE3_adj_mat.txt")
tf.list <- row.names(gen.in)

ranks <- as.data.frame(rowRanks(as.matrix(gen.in)))
row.names(ranks) <- row.names(gen.in)
colnames(ranks) <- colnames(gen.in)
gen <- gen.in %>% rownames_to_column() %>%
  as_tibble() %>%
  pivot_longer(!rowname, names_to = "Paralog", values_to = "value") %>%
  mutate(var = "GENIE3",
         Paralog = str_remove(Paralog, "_.*")) %>%
  rename("TFs" = "rowname") %>%
  mutate(Family = case_when(
    grepl('SVSP', Paralog) ~ 'SVSP',
    grepl('SVMP', Paralog) ~ 'SVMP',
    grepl('PLA2', Paralog) ~ 'PLA2',
    TRUE ~ "Other")) %>%
  select(c("TFs", "Family", "Paralog", "var", "value"))


### ATACseq Footprints #####
files <- list.files(path = "/Users/aundreawestfall/Desktop/Venom_snRNA/summaryFigData/tfbs_bound_tsvs", pattern = "*wBoundInfo*", full.names = T)
fp_res <- rbindlist(sapply(files, read.table, simplify = FALSE, header=T, sep='\t'), idcol = 'feature')

fp_res$CRE <- ifelse(grepl('enhancer', fp_res$feature), 'Enhancer', 'Promoter')

sig.fps <- fp_res %>% mutate(tx_id = str_remove_all(tx_id, ".+-")) %>% 
  mutate(AvgFootprint = rowMeans(select(fp_res, contains("_footprint")), na.rm = TRUE),
         gene_id = str_remove_all(gene_id, "_"),
         tx_id = str_remove_all(tx_id, "PER[0-9]+_|_"),
         tfbs = str_replace(toupper(tfbs), "DDIT3::CEBPA", "DDIT3")) %>% 
  filter(!str_detect(tx_id, "^\\d")) %>%
  mutate(Family = case_when(grepl('SVSP', tx_id) | grepl('SVSP', gene_id) ~ 'SVSP',
                            grepl('SVMP', tx_id) | grepl('SVMP', gene_id) ~ 'SVMP',
                            grepl('PLA2', tx_id) | grepl('PLA2', gene_id) ~ 'PLA2',
                            TRUE ~ "Other")) %>%
  mutate(tx_id = strsplit(tx_id, split = "\\.")) %>%
  unnest(tx_id) %>%
  filter(tx_id != "PLA2K") %>%
  group_by(Family, tfbs, tx_id, CRE) %>%
  summarise(AvgFootprint = mean(AvgFootprint)) %>%
  rename("TFs" = "tfbs", "Paralog" = "tx_id", "var" = "CRE", "value" = "AvgFootprint") %>%
  mutate(Paralog = str_replace(Paralog, "Vespryn", "ohanin")) %>%
  select(c("TFs", "Family", "Paralog", "var", "value"))


### Combine all data ######
gen # GENIE3 results
sig.fps # ATACseq results

dat <- rbind(sig.fps, gen)
dat$var <- factor(dat$var)
dat$TFs <- factor(dat$TFs)
dat$Paralog <- factor(dat$Paralog)



# Relevel Paralog factors to not be alphabetical
dat$Paralog <- factor(dat$Paralog,
                      levels = c("SVSP1","SVSP2","SVSP3","SVSP4","SVSP5","SVSP6","SVSP7","SVSP8","SVSP9","SVSP10","SVSP11",
                                 "SVMP1","SVMP2","SVMP3","SVMP4","SVMP5","SVMP6","SVMP7","SVMP8","SVMP9","SVMP10",
                                 "PLA2A1","PLA2B1","PLA2C1",
                                 "BPP", "CRISP1", "CTL", "EXO1", "EXO2", "EXO3", "LAAO1", "LAAO3", "myotoxin", "ohanin", "VEGF1", "vQC1", "vQC2"))

erk <- c("ATF2", "BLZF1", "CREG1", "EHF", "ERF", "ELK4", "ERF", "FOS", "FOXO3", "FOXO4", "GRHL1", "JUN", "KLF11", "NR4A1", "RARA", "SP1", "TBX3")
upr <- c("ATF6", "BHLHA15", "CCND1", "DDIT3", "XBP1")
both <- c("ATF4", "ATF6B", "CREB3", "CREB3L1", "CREB3L2")
none <- c("CREM", "CTCF", "ELF5", "FOXC2", "GRHL1", "GRHL2", "NCOA2", "NFIA", "NFIB", "NFIX", "NR4A2", "SPDEF", "SREBF2",
          "NFE2L1", "TDG", "NFATC1", "CREBRF", "BMPR1A", "HCFC1", "CLOCK", "RERE", "MEIS1", "CARM1", "NCOA3", "FOSB",
          "SUPT4H1", "CDK7", "TAF13", "NR1H3", "ARNT", "GATAD2B", "GATAD2A", "HDAC1", "HDAC3", "SIRT1")
other <- read.table("~/Desktop/Venom_snRNA/summaryFigData/otherTFs.txt")
other <- other$V1

pars <- levels(dat$Paralog)

dat$Pathways <- NA
dat[which(dat$TFs %in% erk), "Pathways"] <- "ERK"
dat[which(dat$TFs %in% upr), "Pathways"] <- "UPR"
dat[which(dat$TFs %in% both), "Pathways"] <- "both"
dat[which(dat$TFs %in% other), "Pathways"] <- "other"
dat[which(dat$TFs %notin% c(erk, upr, both, other)), "Pathways"] <- "none"

dat$Pathways <- factor(dat$Pathways, levels = c("ERK", "both", "UPR", "other", "none"))

path.dat <- dat %>% filter(TFs %in% c(erk, upr, both))
other.dat <- dat %>% filter(TFs %in% c(other))
none.dat <- dat %>% filter(TFs %notin% c(erk, upr, both, other))

path.dat$Family <- factor(path.dat$Family, levels = c("SVMP", "SVSP", "PLA2", "Other"))
path.dat$var <- factor(path.dat$var, levels = c("GENIE3", "Promoter", "Enhancer"))


### Rescale ATAC peaks to plot with Beta
p.scaleFactor <- max(path.dat[which(path.dat$var %in% c("Promoter", "Enhancer") & path.dat$Paralog != "LAAO3"), "value"]) / max(path.dat[which(path.dat$var == "GENIE3"), "value"])
path.dat[which(path.dat$var %in% c("Promoter", "Enhancer")), "value"] <- path.dat[which(path.dat$var %in% c("Promoter", "Enhancer")), "value"] / p.scaleFactor

## Rescale LAAO3 enhancers to fit scale; DO NOT FORGET TO ADD ASTERISKS IN ILLUSTRATOR
path.dat[which(path.dat$Paralog == "LAAO3" & path.dat$var == "Enhancer"), "value"] <- 0.09

o.scaleFactor <- max(other.dat[which(other.dat$var %in% c("Promoter", "Enhancer")), "value"]) / max(other.dat[which(other.dat$var == "GENIE3"), "value"])
other.dat[which(other.dat$var %in% c("Promoter", "Enhancer")), "value"] <- other.dat[which(other.dat$var %in% c("Promoter", "Enhancer")), "value"] / o.scaleFactor

n.scaleFactor <- max(none.dat[which(none.dat$var %in% c("Promoter", "Enhancer")), "value"]) / max(none.dat[which(none.dat$var == "GENIE3"), "value"])
none.dat[which(none.dat$var %in% c("Promoter", "Enhancer")), "value"] <- none.dat[which(none.dat$var %in% c("Promoter", "Enhancer")), "value"] / n.scaleFactor

###  Plot
path.dat %>%
  ggplot(aes(x = Paralog, y = value)) +
  theme_classic() + 
  theme(strip.placement = "outside",
        strip.text = element_text(size = 7, face = "bold", colour = "black"),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.switch.pad.grid = unit(0.3, "cm"),
        strip.background.y = element_rect(fill = "transparent", colour = "transparent"),
        strip.clip = "off",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = -70),
        panel.grid.major.x = element_line(size = 0.25, color = "gray75"),
        legend.position = "bottom", 
        legend.justification = "left",
        legend.title = element_text(face = "bold")) +
  ggtitle("Summarized adjacency weights and footprints of TFs and venom paralogs") +
  
  ## Correlation bars
  geom_col(aes(fill = var, group = var),
           color = "transparent", width = 1,
           position = position_dodge(preserve = "single")) +
  scale_y_continuous(name = "weight",
                     breaks = scales::pretty_breaks(n = 2),
                     limits = c(0, 0.09),
                     sec.axis = sec_axis(~.*p.scaleFactor,
                                         breaks = scales::pretty_breaks(n = 2),
                                         name = "footprint score")) +

  scale_fill_manual(values = c("salmon", "lightskyblue1", "skyblue3"), na.value = "transparent",
                    name = "CRE",
                    guide = guide_legend(keyheight = 0.75, keywidth = 0.75)) +
  
  geom_hline(yintercept = 0, size = 0.25, color = "black") +
  facet_nested(rows = vars(Pathways, TFs), cols = vars(Family), switch = "both", scales = "free_x", space = "free_x")

other.dat %>%
  ggplot(aes(x = Paralog, y = value)) +
  theme_classic() + 
  theme(strip.placement = "outside",
        strip.text = element_text(size = 7, face = "bold", colour = "black"),
        strip.text.y.left = element_text(angle = 0),
        strip.switch.pad.grid = unit(0.3, "cm"),
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = -42),
        panel.grid.major.x = element_line(size = 0.25, color = "gray75"),
        legend.position = "right",
        legend.title = element_text(face = "bold")) +
  ggtitle("Summarized adjacency weights and footprints of TFs and venom paralogs") +
  
  ## Correlation bars
  geom_col(aes(fill = var, group = var),
           color = "black", size = 0.25, width = 0.8,
           position = position_dodge(preserve = "single")) +
  scale_y_continuous(name = "GENIE3 adjacency weight", 
                     sec.axis = sec_axis(~.*o.scaleFactor,
                                         name = "ATACseq footprint score")) +
  
  scale_fill_manual(values = c("skyblue3", "lightskyblue1", "salmon"), na.value = "transparent",
                    name = "",
                    breaks = c("Promoter", "Enhancer", "GENIE3"),
                    guide = guide_legend(keyheight = 0.75, keywidth = 0.75)) +
  
  geom_hline(yintercept = 0, size = 0.25, color = "black") +
  facet_nested(rows = vars(TFs), cols = vars(Family), switch = "both", scales = "free", space = "free_x")

none.dat %>%
  ggplot(aes(x = Paralog, y = value)) +
  theme_classic() + 
  theme(strip.placement = "outside",
        strip.text = element_text(size = 7, face = "bold", colour = "black"),
        strip.text.y.left = element_text(angle = 0),
        strip.switch.pad.grid = unit(0.3, "cm"),
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = -30),
        panel.grid.major.x = element_line(size = 0.25, color = "gray75"),
        legend.position = "right",
        legend.title = element_text(face = "bold")) +
  ggtitle("Summarized GENIE3 adjacency weights and ATACseq footprints of TFs and venom paralogs") +
  
  ## Correlation bars
  geom_col(aes(fill = var, group = var),
           color = "black", size = 0.25, width = 0.8,
           position = position_dodge(preserve = "single")) +
  scale_y_continuous(name = "weight",
                     sec.axis = sec_axis(~.*n.scaleFactor,
                                         name = "footprint score")) +
  
  scale_fill_manual(values = c("skyblue3", "lightskyblue1", "salmon"), na.value = "transparent",
                    name = "CRE",
                    breaks = c("Promoter", "Enhancer", "GENIE3"),
                    guide = guide_legend(keyheight = 0.75, keywidth = 0.75)) +
  
  geom_hline(yintercept = 0, size = 0.25, color = "black") +
  facet_nested(rows = vars(TFs), cols = vars(Family), switch = "both", scales = "free_x", space = "free_x")

