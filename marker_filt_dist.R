########################################################################################
######################## Modified from James McNellie 5/15/2024 ########################
# This function is designed to thin out markers in genomic data by ensuring that no ####
# two markers within 125,000 base pairs of each other are retained unless they are part#
# of a trio where distance dictates a random selection. This approach helps in reducing#
# marker density, which can be beneficial for analyses where high marker density might #
# obscure genetic signals or increase computational burdens. ###########################
########################################################################################

#Iterate through the function several times until the number of SNPs is within 5000-10000 for QTL Mapping 
library(vcfR)

# Clear existing variables from the workspace
rm(list=ls(all=TRUE))

# Set working directory
setwd("/mmfs1/projects/brent.hulke/HaplotypeStructureForGenomicSelection/Nectar_Remapping")

# Load VCF file
vcf_data <- read.vcfR("Nectar_18LabelsAndSeqIdentifiers_withParents_Imputed_biallelicSNPs_ChiSquareFilterP0.1.vcf.gz")

# Extract necessary parts for filtering, linking directly with the original data
chrom_pos_df <- data.frame(chr = vcf_data@fix[, "CHROM"],
                           pos = vcf_data@fix[, "POS"],
                           row.names = row.names(vcf_data@fix))
chrom_pos_df$chr <- gsub("Ha412HOChr", "", chrom_pos_df$chr)  # Adjust chromosome names if necessary
chrom_pos_df$pos <- as.numeric(as.character(chrom_pos_df$pos))  # Ensure pos is numeric

# Function to filter markers based on distance criteria within each chromosome
marker_filt_dist <- function(data_in) {
  mtk <- c()
  for(i in unique(data_in$chr)) {
    temp_df <- data_in[data_in$chr == i,]
    for(j in seq(1, nrow(temp_df) - 2, 3)) {
      td2 <- temp_df[j:(j+2),]
      temp_dist <- td2[nrow(td2), 'pos'] - td2[1, 'pos']
      if(temp_dist < 125000) {
        mtk <- c(mtk, sample(row.names(td2), 1))
      } else {
        mtk <- c(mtk, row.names(td2))
      }
    }
  }
  mtk <- unique(as.numeric(mtk))
  return(mtk)
}

# Apply the custom filtering function and extract the indices
filtered_indices <- marker_filt_dist(chrom_pos_df)

# Subset the original vcfR object using the filtered indices
vcf_filtered <- vcf_data[filtered_indices, ]

# Save the filtered VCF to a file
write.vcf(vcf_filtered, file = "Nectar_18LabelsAndSeqIdentifiers_withParents_Imputed_biallelicSNPs_ChiSquareFilterP0.1_distfiltered2.vcf.gz")
vcf_data <- read.vcfR("Nectar_18LabelsAndSeqIdentifiers_withParents_Imputed_biallelicSNPs_ChiSquareFilterP0.1_distfiltered14.vcf.gz")
