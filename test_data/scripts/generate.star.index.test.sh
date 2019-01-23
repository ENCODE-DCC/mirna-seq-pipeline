#! /bin/bash
#$ -N star.log
#$ -pe openmp 8
#$ -q sam128,sam,bio
#$ -ckpt blcr
#$ -l kernel=blcr
#$ -r y
#$ -l virtual_free=48G

export PATH=$PATH:~/software/STAR/2.5.1b/STAR-2.5.1b/source/
genomedir="/share/samdata/sorenar/miRNA/pipeline/mirna-seq-pipeline/test_data/refs/test/"
annotation="/share/samdata/sorenar/miRNA/pipeline/mirna-seq-pipeline/test_data/refs/test.chr19.gtf"
refgenome="/share/samdata/sorenar/miRNA/pipeline/mirna-seq-pipeline/test_data/refs/test.ref.fasta"

STAR \
	--runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir $genomedir  \
	--sjdbGTFfile $annotation \
	--sjdbOverhang 1 \
	--genomeFastaFiles $refgenome
