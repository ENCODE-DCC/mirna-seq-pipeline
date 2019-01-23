#$ -N log.star
##$ -o
##$ -e /dev/null
#$ -pe openmp 8
#$ -q sam,class,abio128,bio
#$ -ckpt restart

slid=$1

export PATH=$PATH:~/software/STAR/2.5.1b/STAR-2.5.1b/source/:/share/samdata/sorenar/miRNA/pipeline/src/

results_path="/share/samdata/sorenar/miRNA/pipeline/mirna-seq-pipeline/test_data/results/"
miRNA_annotation="/share/samdata/sorenar/miRNA/pipeline/mirna-seq-pipeline/test_data/refs/test.chr19.miRNA.gtf"
genomedir="/share/samdata/sorenar/miRNA/pipeline/mirna-seq-pipeline/test_data/refs/test"
chromref="/share/samdata/sorenar/miRNA/pipeline/mirna-seq-pipeline/test_data/refs/test.chrom.sizes.tsv"

echo $results_path
echo $miRNA_annotation
echo $genomedir
echo $chromref

params='--runThreadN 8
	--alignEndsType EndToEnd
	--outFilterMismatchNmax 1
	--outFilterMultimapScoreRange 0
	--quantMode TranscriptomeSAM GeneCounts
	--outReadsUnmapped Fastx
	--outSAMtype BAM SortedByCoordinate
	--outFilterMultimapNmax 10
	--outSAMunmapped Within
	--outFilterScoreMinOverLread 0
	--outFilterMatchNminOverLread 0
	--outFilterMatchNmin 16
	--alignSJDBoverhangMin 1000
	--alignIntronMax 1 
	--outWigType wiggle
	--outWigStrand Stranded
	--outWigNorm RPM
'
STAR --genomeDir $genomedir --readFilesIn $results_path$slid"_trim.fastq" --outFileNamePrefix $results_path --sjdbGTFfile $miRNA_annotation $params

mv $results_path"Aligned.sortedByCoord.out.bam" $results_path$slid".bam"
mv $results_path"ReadsPerGene.out.tab" $results_path$slid".tsv"

wigToBigWig $results_path"Signal.UniqueMultiple.str1.out.wig" $chromref $results_path$slid".signal.all.plus.bigWig"
wigToBigWig $results_path"Signal.UniqueMultiple.str2.out.wig" $chromref $results_path$slid".signal.all.minus.bigWig"
wigToBigWig $results_path"Signal.Unique.str1.out.wig" $chromref $results_path$slid".signal.unique.plus.bigWig"
wigToBigWig $results_path"Signal.Unique.str2.out.wig" $chromref $results_path$slid".signal.unique.minus.bigWig"

md5sum $results_path$slid".bam" > $results_path"md5"
md5sum $results_path$slid".tsv" >> $results_path"md5"
md5sum $results_path*".bigWig" >> $results_path"md5"

md5sum $results_path* > $results_path"all.md5"
