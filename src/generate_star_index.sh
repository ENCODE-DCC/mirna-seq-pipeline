#! /bin/bash
#$ -N star.log
#$ -pe openmp 16
#$ -q sam,bio
#$ -ckpt blcr
#$ -l kernel=blcr
#$ -r y
#$ -l virtual_free=48G

export PATH=$PATH:/share/samdata/sorenar/miRNA/pipeline/src/
genomedir="/share/samdata/sorenar/miRNA/pipeline/refs/test/"
annotation="/share/samdata/sorenar/miRNA/pipeline/refs/test.chr19.gtf"
refgenome="/share/samdata/sorenar/miRNA/pipeline/refs/test.ref.fasta"

STAR \
	--runThreadN 16 \
	--runMode genomeGenerate \
	--genomeDir $genomedir  \
	--sjdbGTFfile $annotation \
	--sjdbOverhang 1 \
	--genomeFastaFiles $refgenome
