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

---

## **Citation**

Barstow, A.C., McNellie, J.P., Smart, B.C., Keepers, K.G., Prasifka, J.R., Kane, N.C., & Hulke, B.S. (2024).  
*Variant filters using segregation information improve mapping of nectar-production genes in sunflower (Helianthus annuus L.).*  
**The Plant Genome**.
