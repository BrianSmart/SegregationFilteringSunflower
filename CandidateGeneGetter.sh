echo "Chromosome	FlankingMk1Pos	FlankingMk2Pos	Count	Genes" > AshleyCandidateGenes.txt


cat Windows | while read chrom startpos endpos; do

	awk -v chr=$chrom -v strt=$startpos -v nd=$endpos '$1 == chr && $4 > strt && $5 < nd && $3 == "gene" {print $9}' HAN412_Eugene_curated_v1_1.gff3 | cut -d ';' -f2 > GenesTmp

	count=`wc -l GenesTmp | awk '{print $1}'`

	if [[ `ls -tor GenesTmp | awk '{print $4}'` -eq 0 ]]; then
		echo "NoGenesTmp" > UnwrappedGeneListTmp
	else
                #Unwrap gene list so they all fit on one line
        	tr '\n' ',' < GenesTmp | sed 's/ID=gene://g' | sed 's/,$//g' > UnwrappedGeneListTmp
	fi

	echo "$chrom	$startpos	$endpos	$count" > chrompos
	paste chrompos UnwrappedGeneListTmp > newline
	cat AshleyCandidateGenes.txt newline > tmp; mv tmp AshleyCandidateGenes.txt

done
rm GenesTmp; rm UnwrappedGeneListTmp; rm chrompos; rm newline
