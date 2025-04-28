# Variant Filters Using Segregation Information

This repository contains the code used for the study:
*"Variant Filters Using Segregation Information Improve Mapping of Nectar-Production Genes in Sunflower (Helianthus annuus L.)"*.
The study evaluates the impact of biologically informed variant filtering strategies on QTL mapping, demonstrating improved identification of candidate genes related to nectar production.

---

## **Contents**

### **CandidateGeneGetter.sh**
This shell script extracts candidate genes from a GFF annotation file (`HAN412_Eugene_curated_v1_1.gff3`) based on genomic regions specified in the `Windows` file. For each region (defined by chromosome, start position, and end position), it identifies all genes falling entirely within that window, counts them, and outputs the region information along with a comma-separated list of gene IDs to `AshleyCandidateGenes.txt`.

### **Chi_square_template.R**
This R script filters genomic markers using a chi-square test based on expected segregation ratios. The script is designed as a template that can be adjusted for different population types by modifying the expected ratios. The default values (48.4375% homozygous for each allele and 3.125% heterozygous) are set for F6 inbred lines, but can be modified to match the segregation expectations of any population being filtered. It retains markers whose observed genotype frequencies do not significantly deviate from expectations (p > 0.1), removing markers with segregation distortion that could interfere with accurate QTL identification.

### **mapping.R**
This R script performs QTL (Quantitative Trait Locus) mapping using the `qtl` package. It includes code for three distinct "Approaches," likely representing analyses performed on different datasets or using varied marker filtering strategies (`Approach1.csv`, `Approach2.csv`, `Approach3.csv`). The script covers data loading, genetic map estimation and refinement (including custom marker thinning functions and visualization of recombination frequencies), calculation of genotype probabilities, performing 1D (`scanone`), Composite Interval (`cim`), and 2D (`scantwo`) QTL scans, significance testing via permutations, and refining QTL models (`fitqtl`, `refineqtl`).

### **marker_filt_dist.R**
This R script filters genomic markers from a VCF file by removing markers within 125,000 bp of each other. It optimizes marker density while maintaining genome-wide coverage, ensuring the filtered set is suitable for QTL mapping and identifying genomic regions linked to nectar-production traits in sunflower.

### **proc freq marker data.sas**
This SAS script filters genetic markers based on segregation patterns. It utilizes `PROC FREQ` to calculate genotype frequencies for biallelic markers (assuming three genotype classes) and performs chi-square tests against expected segregation ratios (e.g., specified test probabilities like 0.484375, 0.03125, 0.484375, corresponding to F6 expectations). Markers significantly deviating from these expectations (p < 0.10 in this script) are identified and potentially excluded from downstream analyses, similar in principle to `Chi_square_template.R` but implemented within the SAS environment for specific datasets (`markers.bialw`).

### **thinning_loop.R**
This R script thins genomic markers based on inter-marker distance thresholds, identifying and removing redundant or closely spaced markers. It helps refine marker sets to balance genome coverage and computational efficiency, improving QTL mapping precision in the study of sunflower nectar-production traits. *(Note: Similar custom functions are also included within `mapping.R`)*.

### **Windows**
This plain text file serves as input for the `CandidateGeneGetter.sh` script. Each line defines a genomic window with three columns: Chromosome, Start Position, and End Position. These windows likely represent regions of interest identified through QTL mapping or other analyses.

---

## **Citation**

Barstow, A.C., McNellie, J.P., Smart, B.C., Keepers, K.G., Prasifka, J.R., Kane, N.C., & Hulke, B.S. (2025).
*Variant filters using segregation information improve mapping of nectar-production genes in sunflower (Helianthus annuus L.).*
**The Plant Genome**.
