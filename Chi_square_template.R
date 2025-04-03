# Save output

sink("name.Rout")

# Install and load necessary package
#install.packages("vcfR")

library(vcfR)

# Clear existing variables from the workspace
rm(list=ls(all=TRUE))

# Set working directory
setwd("") 

# Load VCF file
vcf <- read.vcfR("")

# Extract the genotype data into a matrix
gt_matrix <- extract.gt(vcf, element="GT")

# Drop the parent HA434 and 456, the 193rd and 194th index (Ignores parents in chi-square test but keeps them for r/qtl) 
gt_matrix <- gt_matrix[, -c(193:194)]

# Total number of markers before filtering
total_markers_before <- nrow(gt_matrix)
print(paste("Total markers before filtering:", total_markers_before))

# Define a modified chi-square test filter to include monomorphic sites
chi_square_filter <- function(row) {

  # Count the number of homozygous and heterozygous genotypes
  homozygous_count0 <- sum(row %in% c("0|0", "0/0"))
  homozygous_count1 <- sum(row %in% c("1|1", "1/1"))
  heterozygous_count <- sum(row %in% c("1|0", "0|1", "1/0", "0/1"))
  total_count <- sum(homozygous_count0, homozygous_count1, heterozygous_count)

  # Calculate the expected count based on F6 inbred expectations -replace with your expected ratio
  expected_homozygous0 <- 0.484375 * total_count
  expected_homozygous1 <- 0.484375 * total_count
  expected_heterozygous <- 0.03125 * total_count

  # Perform the chi-square test
  observed <- c(homozygous_count0, homozygous_count1, heterozygous_count)
  expected <- c(expected_homozygous0, expected_homozygous1, expected_heterozygous)
  chisq_test <- chisq.test(observed, p=expected/total_count)

  # Return TRUE if p-value > 0.1, indicating the marker should be retained
  return(chisq_test$p.value > 0.1)
}

# Apply the chi-square test to the genotype data
filter_results <- apply(gt_matrix, 1, chi_square_filter)
vcf_filtered_by_chisq <- vcf[filter_results, ]

# Total number of markers after filtering
total_markers_after <- sum(filter_results)
polymorphic_count <- total_markers_after

# Print summary statistics
cat("Total markers before filtering:", total_markers_before, "\n")
cat("Total polymorphic markers retained:", polymorphic_count, "\n")
cat("Markers retained: ", total_markers_after, "/", total_markers_before, "(", round(total_markers_after / total_markers_before * 100, 2), "%)\n")

# Save filtered VCF
write.vcf(vcf_filtered_by_chisq, file = "")

# Stop capturing output
sink()
