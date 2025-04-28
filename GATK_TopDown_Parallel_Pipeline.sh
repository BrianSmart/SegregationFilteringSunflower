#Pipeline for calling variants on many samples with GATK, using parallelization. Pipeline starts with fastqs and finishes with a multi-sample VCF

#Author: Kyle Keepers
#Date: 5/10/2023

#Required: Fastqs must be uniquely named after the sample prefix (NOT the well position) and end with _R1.fastq.gz and _R2.fastq.gz. All samples to be processed in this pipeline must be in a single folder whose path is specified in the variable assignment below.

#Required programs:
##fastp (https://github.com/OpenGene/fastp, or try ("conda install -c bioconda fastp")
##bwamem2 (https://github.com/bwa-mem2/bwa-mem2)
##samtools (https://www.htslib.org/download/)
##gatk (https://github.com/broadinstitute/gatk/releases)
##GNU parallel (https://www.gnu.org/software/parallel/)


#Set these variables before running:
threads=53

#Locations of programs
#pathToGATK=/home/ziv/gatk-4.1.8.1/gatk ##INCLUDE EXECUTABLE IN PATH
module load gatk
module load parallel
module load samtools
module load bcftools
pathToFastp=/mmfs1/projects/brent.hulke/Software/fastp #INCLUDE EXECUTABLE IN PATH
pathTobwamem2=/mmfs1/projects/brent.hulke/Software/bwa-mem2/bwa-mem2 #INCLUDE EXECUTABLE IN PATH
#pathToParallel=/usr/bin/parallel #INCLUDE EXECUTABLE IN PATH

#Locations of files
pathToBED=/mmfs1/projects/brent.hulke/HA412HOv2_w_CPMT/Ha412HOv2_w_CPMT_split53.fa.bed #INCLUDE FILE IN PATH
pathToReference=/mmfs1/projects/brent.hulke/F6_and_F4_Crosses_for_Ashley/Ha412HOv2_w_CPMT.fa #INCLUDE FILE IN PATH
pathToFastqs=/mmfs1/projects/brent.hulke/SMGT_paper/SMGT_box21_22_raw_GBS_fastq/SMGT_box21_22_raw_fastq #NO SLASH AT END OF PATH


#Index reference in proprietary bwa-mem2 format, if not done already
#if [[ -z ${pathToReference}.bwt.2bit.64 ]]; then
#	$pathTobwamem2 index $pathToReference
#else
#	echo "Reference is already indexed in bwa-mem2 format. Skipping indexing step."
#fi


#Get sample prefixes
#ls $pathToFastqs/*q.gz | awk '{print $NF}' | sed 's@'"$pathToFastqs"'@@g' | sed 's/\///g' | sed 's/_[12]\.fq\.gz//g' | sed 's/\.fq//g' | sed 's/\.fastq//g' | sed 's/\.gz//g' | sort | uniq > SamplePrefixes


##########################Step 1: Generate indexed alignment maps##########################
cat ToRedo | while read sample; do

	#Create a readgroup string to be added in the read-alignment step
	rg_string="@RG\tID:"$sample"\tLB:"$sample"\tPL:Illumina\tSM:$sample"

	#Uncomment if reads are paired-end. Perform automated read quality control filtering
#	$pathToFastp -i $pathToFastqs/${sample}*1.f*q.gz -I $pathToFastqs/${sample}*2.f*q.gz --stdout | $pathTobwamem2 mem -R "$rg_string" -t $threads $pathToReference - | samtools view -bhu - | samtools collate -@ $threads -O - | samtools fixmate -@ $threads -m - - | samtools sort -@ $threads - | samtools markdup -@ $threads - ${sample}.sorted.bam
#	rm fastp.html; rm fastp.json	

	#Uncomment if reads are single-end. Perform automated read quality control filtering
	$pathToFastp -i $pathToFastqs/${sample}.fq.gz --stdout | $pathTobwamem2 mem -R "$rg_string" -t $threads $pathToReference - | samtools view -bhu - | samtools collate -@ $threads -O - | samtools fixmate -@ $threads -m - - | samtools sort -@ $threads - | samtools markdup -@ $threads - ${sample}.sorted.bam
	rm fastp.html; rm fastp.json

	#Index the completed alignment map
	samtools index ${sample}.sorted.bam

done


########################################################################################
##########################Step 2: Call variants in individuals##########################
########################################################################################

cat ToRedo | while read sample; do

	YOUR_BAM=$sample
	CONTIGS=$pathToBED

## A function for bcftools where samtools view extracts only the indicated regions:
#Modified from the thread at https://www.biostars.org/p/418748/

	function CALL_HAPLOTYPES {
		BAM=$1
		REGION=$2
		GATK=$3
		REF=$4
		gatk HaplotypeCaller --java-options "-Xmx4G" --native-pair-hmm-threads 1 --intervals $REGION --input ${BAM}.sorted.bam --output output_fromBAM_${BAM}_region_${REGION}.vcf.gz --reference $REF -ERC GVCF

	}; export -f CALL_HAPLOTYPES


	#Take each region in the .bed input file and send it through the MPILEUP function with its own thread
	cat $CONTIGS | awk '{print $1":"$2"-"$3}' | parallel --will-cite -j $threads "CALL_HAPLOTYPES ${YOUR_BAM} {} gatk $pathToReference"

	#Create a list of split .vcf.gz in the correct order to concatenate into the final whole-genome vcf.
	sed "s/QQQQQ/${YOUR_BAM}/g" FileListCopy > FileList

	## Concatenate everything using the FileList
	bcftools concat --file-list FileList -O z > concat_fromBAM_${YOUR_BAM}.g.vcf.gz

	## remove temporary files:
	rm output_fromBAM_${YOUR_BAM}_region_*.vcf.gz*; rm FileList


	#Split the genome into many pieces to run variant calling in parallel	
	gatk IndexFeatureFile --input concat_fromBAM_${sample}.g.vcf.gz
done


########################################################################################
############################Step 3: GenomicsDBImport####################################
########################################################################################
#This script creates a database from the reference genome, I guess. The GenomicsDBImport/GenotypeGVCFs steps are required for a large number of samples, instead of the CombineGVCFs command that works for a smaller number of samples.
#rm -r OutputFromGenomicsDBImport

#Create some required files for GenomicsDBImport step:
#awk '{print $1}' $pathToBED | uniq > chr.list

#>col2
#cat SamplePrefixes | while read sample; do
#	ls concat_fromBAM_${sample}.g.vcf.gz | awk '{print $NF}' >> col2
#done
#paste SamplePrefixes col2 > CohortSample.map
#rm col2

	#This step takes a while
#	gatk --java-options "-Xmx150g -Xms150g" \
#		GenomicsDBImport \
#		--genomicsdb-workspace-path OutputFromGenomicsDBImport \
#		--batch-size 50 \
#		--intervals chr.list \
#		--sample-name-map CohortSample.map \
#		--reader-threads $threads

########################################################################################
#############################Step 4: Combine GVCFs######################################
########################################################################################
## A function for bcftools where samtools view extracts only the indicated regions:

function COMBINEVCFS_AND_FILTER {

	REGION=$1
	GATK=$2
	REF=$3
	$GATK --java-options "-Xms3g -Xmx3g" GenotypeGVCFs --intervals $REGION --reference $REF -V gendb://OutputFromGenomicsDBImport -O Combined_${REGION}.vcf.gz
	#Keep only biallelic SNPs
	#bcftools view --max-alleles 2 --exclude-types indels Combined_${REGION}.vcf.gz --threads 1 -o Combined_${REGION}_snps.vcf.gz
	#rm Combined_${REGION}.vcf.gz

	#Filter out by missingness and minor allele frequency. These two filters need different values for the organellar genomes versus the nuclear genome
	#vcftools --gzvcf Combined_${REGION}_snps.vcf.gz --max-missing 0.4 --maf 0.05 --recode --recode-INFO-all --out Combined_${REGION}_MM_MAF_Q.vcf.gz

	#Reformat missing data
	#sed 's/\t\.:\./\t\.|\.:\./g' Combined_${REGION}_MM_MAF_Q.vcf.gz.recode.vcf |  sed 's/\t\.:0,0,0/\t\.|\.:0,0,0/g' | sed 's/\t\.:0,0,1/\t\.|\.:0,0,0/g' | sed 's/\t\.:1,0,0/\t\.|\.:0,0,0/g' | sed 's/\t\.:0,1,0/\t\.|\.:0,0,0/g' > Combined_${REGION}_reformatted.vcf

	#Housekeeping
	#rm Combined_${REGION}_MM_MAF_Q.vcf.gz.recode.vcf; rm Combined_${REGION}_snps.vcf.gz
	#mv Combined_${REGION}_reformatted.vcf Combined_${REGION}.vcf
	#gzip Combined_${REGION}.vcf

}; export -f COMBINEVCFS_AND_FILTER

#Take each region in the .bed input file and send it through the COMBINEVCFS_AND_FILTER function with its own thread
cat $pathToBED | awk '{print $1":"$2"-"$3}' | parallel --will-cite -j $threads "COMBINEVCFS_AND_FILTER {} gatk $pathToReference"

#Make "Regions" file:
awk '{print "Combined_"$1":"$2"-"$3".vcf.gz"}' $pathToBED > Regions

## Concatenate everything using the FileList
bcftools concat --file-list Regions -O z > Combined.vcf.gz

## remove temporary files:
rm Combined_*

echo "Finished pipeline! Outfile is named Combined.vcf.gz. Rename if desired."
