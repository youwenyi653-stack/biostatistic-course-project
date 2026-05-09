install.packages(c("readxl", "ggplot2", "dplyr", "pheatmap"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano", "org.Hs.eg.db"))
BiocManager::install("RNASeqPower")
BiocManager::install("edgeR")
library(RNASeqPower)
suppressPackageStartupMessages({
  library(DESeq2)
  library(readxl)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(org.Hs.eg.db)
  library(dplyr)
})
setwd('/Users/youw/Documents/course/biostatistics')
save.image('crispr2.RData')
load('crispr2.RData')
## ===================== input data =====================
## counts: rows = genes, cols = samples

## ===================== load counts =====================
counts <- read.table(
  "crispr_gene_counts.rename.txt",
  header = TRUE,
  sep = "",        
  skip = 1,        
  comment.char = "",  
  stringsAsFactors = FALSE,
  check.names = FALSE
)
gene_ids <- counts$Geneid
count_matrix <- counts[, 7:ncol(counts)]
rownames(count_matrix) <- gene_ids
colnames(count_matrix)
summary(count_matrix)
#====================== load meta =====================
meta_condition <- read_excel('Metadata_RNA-Seq_CRISPRi_326E_KD_iPSCs.xlsx', col_names = TRUE, sheet = 3)
meta_guide <- read_excel('Metadata_RNA-Seq_CRISPRi_326E_KD_iPSCs.xlsx', col_names = TRUE, sheet = 1)
#==========power analysis=====================
meta_condition$batch
cata <- meta_condition$Description[7:63]
cata <- factor(cata)
real_compare_count <- count_matrix[,7:63]
dge <- DGEList(counts = real_compare_count, group = cata)
dge <- estimateDisp(dge, design=NULL) 
dge$samples$short_group <- ifelse(grepl("Control", dge$samples$group),
                                  "Control",
                                  sub(".*\\+ (.*) knock-down.*", "\\1", dge$samples$group))
kd_genes <- unique(dge$samples$short_group[dge$samples$short_group != "Control"])

power_results <- data.frame(KD = character(), power = numeric(), stringsAsFactors = FALSE)

for (kd in kd_genes) {
  n_control <- sum(dge$samples$short_group == "Control")
  n_kd      <- sum(dge$samples$short_group == kd)
  n_use     <- min(n_control, n_kd)
  depth <- mean(dge$samples$lib.size[dge$samples$short_group %in% c("Control", kd)])/1e6
  cv <- sqrt(dge$common.dispersion)
  pwr <- rnapower(n = n_use, cv = cv, effect = 1.2, alpha = 0.05, depth=depth)
  
  power_results <- rbind(power_results, data.frame(KD = kd, power = pwr))
}

print(power_results)

by(dge$samples$lib.size, dge$samples$short_group, summary)
#===========power figure with original mean
kd <- kd_genes[1]

n_seq <- 6:20

# baseline depth
depth_base <- mean(dge$samples$lib.size[
  dge$samples$short_group %in% c("Control", kd)
]) / 1e6

# CV
cv <- sqrt(dge$common.dispersion)
# baseline power figure
power_vec <- sapply(n_seq, function(n) {
  rnapower(n = n, cv = cv, effect = 1.2, alpha = 0.05, depth = depth_base)
})

plot(n_seq, power_vec, type = "l", lwd = 2,
     ylim = c(0,1),
     xlab = "Sample size",
     ylab = "Power",
     main = paste0("Power curve of ", kd))

# depth scaling
depth_factors <- c(0.8, 0.85, 0.9, 0.95, 1.05, 1.1, 1.15, 1.2)


par(mfrow = c(2, 4))

for (f in depth_factors) {
  
  depth_i <- depth_base * f
  
  power_vec <- sapply(n_seq, function(n) {
    rnapower(n = n, cv = cv, effect = 1.2, alpha = 0.05, depth = depth_i)
  })
  
  plot(n_seq, power_vec, type = "l", lwd = 2,
       ylim = c(0,1),
       xlab = "Sample size",
       ylab = "Power",
       main = paste0("Depth × ", f))
  

  abline(h = 0.8, col = "red", lty = 2)
}
