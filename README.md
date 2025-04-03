# Variant Filters Using Segregation Information

This repository contains the code used for the study:  
*"Variant Filters Using Segregation Information Improve Mapping of Nectar-Production Genes in Sunflower (Helianthus annuus L.)"*.
The study evaluates the impact of biologically informed variant filtering strategies on QTL mapping, demonstrating improved identification of candidate genes related to nectar production.

---

## **Contents**

### **marker_filt_dist.R**  
This R script filters genomic markers from a VCF file by removing markers within 125,000 bp of each other. It optimizes marker density while maintaining genome-wide coverage, ensuring the filtered set is suitable for QTL mapping and identifying genomic regions linked to nectar-production traits in sunflower.

### **thinning_loop.R**  
This R script thins genomic markers based on inter-marker distance thresholds, identifying and removing redundant or closely spaced markers. It helps refine marker sets to balance genome coverage and computational efficiency, improving QTL mapping precision in the study of sunflower nectar-production traits.

### **Chi_square_template.R**
This R script filters genomic markers using a chi-square test based on expected segregation ratios. The script is designed as a template that can be adjusted for different population types by modifying the expected ratios. The default values (48.4375% homozygous for each allele and 3.125% heterozygous) are set for F6 inbred lines, but can be modified to match the segregation expectations of any population being filtered. It retains markers whose observed genotype frequencies do not significantly deviate from expectations (p > 0.1), removing markers with segregation distortion that could interfere with accurate QTL identification.

---

## **Citation**

Barstow, A.C., McNellie, J.P., Smart, B.C., Keepers, K.G., Prasifka, J.R., Kane, N.C., & Hulke, B.S. (2025).  
*Variant filters using segregation information improve mapping of nectar-production genes in sunflower (Helianthus annuus L.).*  
**The Plant Genome**.
