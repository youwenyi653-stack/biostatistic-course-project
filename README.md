This script performs RNA-seq power analysis for CRISPRi knockdown experiments using edgeR and RNASeqPower. First, the gene count matrix and sample metadata are imported and processed to generate a count matrix containing gene expression counts for each sample. Experimental conditions are extracted from the metadata and used to define control and knockdown groups.

An edgeR DGEList object is then created to estimate the common dispersion of the dataset, which reflects the biological variability among samples. The coefficient of variation (CV) is calculated from the estimated dispersion and used as an input for downstream power analysis.

For each knockdown condition, the script calculates the statistical power to detect differential expression between knockdown and control samples. Power is estimated based on sample size, sequencing depth, effect size, and dispersion using the RNASeqPower package.

The workflow also generates power curves across a range of sample sizes and evaluates how changes in sequencing depth influence statistical power. Multiple sequencing depth scaling factors are tested to simulate different experimental conditions. The resulting plots help determine the sample size and sequencing depth required to achieve sufficient power, commonly defined as 80%.
