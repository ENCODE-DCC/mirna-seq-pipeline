# ENCODE micro rna seq pipeline: main pipeline
# Maintainer: Otto Jolanki

#CAPER docker quay.io/encode-dcc/mirna-seq-pipeline:v1.0
#CAPER singularity docker://quay.io/encode-dcc/mirna-seq-pipeline:v1.0
#CROO out_def https://storage.googleapis.com/encode-pipeline-output-definition/mirna.output_definition.json

workflow mirna_seq_pipeline {
    #File inputs

    #cutadapt
    #Array containing the input fastq files
    Array[File] fastqs
    #Array containing Fasta files with 5' adapter sequence(s), in the same order as the fastqs
    Array[File] five_prime_adapters
    #Fasta file with 3' adapter sequence(s)
    File three_prime_adapters

    #star
    #tar.gz archive that contains star index
    File star_index
    #micro-rna annotations
    File mirna_annotation

    #wigtobigwig
    #tsv with chromosome sizes
    File chrom_sizes

    #common
    #Prefix for outputs (additionally replicate number will be added)
    String experiment_prefix

    #Resources

    #cutadapt
    Int cutadapt_ncpus
    Int cutadapt_ramGB
    String cutadapt_disk

    #star
    Int star_ncpus
    Int star_ramGB
    String star_disk

    #wigtobigwig
    Int wigtobigwig_ncpus
    Int wigtobigwig_ramGB
    String wigtobigwig_disk

    #Pipeline starts here

    scatter (i in range(length(fastqs))) {
        call cutadapt { input:
            fastq = fastqs[i],
            five_prime_adapters = five_prime_adapters[i],
            three_prime_adapters = three_prime_adapters,
            output_prefix = "rep"+(i+1)+experiment_prefix,
            ncpus = cutadapt_ncpus,
            ramGB = cutadapt_ramGB,
            disk = cutadapt_disk,
        }

        call star { input:
            fastq = cutadapt.trimmed_fastq,
            index = star_index,
            annotation = mirna_annotation,
            output_prefix = "rep"+(i+1)+experiment_prefix,
            ncpus = star_ncpus,
            ramGB = star_ramGB,
            disk = star_disk,
            }

        call wigtobigwig { input:
            plus_strand_all_wig = star.plus_strand_all_wig,
            minus_strand_all_wig = star.minus_strand_all_wig,
            plus_strand_unique_wig = star.plus_strand_unique_wig,
            minus_strand_unique_wig = star.minus_strand_unique_wig,
            chrom_sizes = chrom_sizes,
            output_prefix = "rep"+(i+1)+experiment_prefix,
            ncpus = wigtobigwig_ncpus,
            ramGB = wigtobigwig_ramGB,
            disk = wigtobigwig_disk,
        }
    }

    #If there are exactly two replicates, calculate spearman correlation between quants in replicates
    if (length(fastqs) == 2) {
        call spearman_correlation { input:
            quants = star.tsv,
            output_filename = experiment_prefix+"_spearman.json",
            }
    }
}

#Task definitions

task cutadapt {
    File fastq
    File five_prime_adapters
    File three_prime_adapters
    String output_prefix
    Int ncpus
    Int ramGB
    String disk

    command {
        cutadapt \
            -a file:${three_prime_adapters} \
            -e 0.25 \
            --match-read-wildcards \
            --untrimmed-output=${output_prefix + "_NO3AD.fastq"} \
            ${fastq} \
            | cutadapt \
            -e 0.34 \
            --match-read-wildcards \
            --no-indels \
            -m 15 \
            -O 6 \
            -n 1 \
            -g file:${five_prime_adapters} \
            --untrimmed-output=${output_prefix + "_NO5AD.fastq"} \
            --too-short-output=${output_prefix + "_SHORT_FAIL.fastq"} \
            - > ${output_prefix + "_trim.fastq"}
    }

    output {
        File no3ad_untrimmed_fastq = glob("*_NO3AD.fastq")[0]
        File no5ad_untrimmed_fastq = glob("*_NO5AD.fastq")[0]
        File too_short_fastq = glob("*_SHORT_FAIL.fastq")[0]
        File trimmed_fastq = glob("*_trim.fastq")[0]
    }

    runtime {
        cpu: ncpus
        memory: "${ramGB} GB"
        disks: disk
    }
}

task star {
    File fastq
    File index
    File annotation
    String output_prefix
    Int ncpus
    Int ramGB
    String disk

    command {
        tar -xzvf ${index}
        gzip -cd ${annotation} > anno.gtf
        STAR \
            --genomeDir out \
            --readFilesIn ${fastq} \
            --sjdbGTFfile anno.gtf \
            --runThreadN ${ncpus} \
            --alignEndsType EndToEnd \
            --outFilterMismatchNmax 1 \
            --outFilterMultimapScoreRange 0 \
            --quantMode TranscriptomeSAM GeneCounts \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 10 \
            --outSAMunmapped Within \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMatchNmin 16 \
            --alignSJDBoverhangMin 1000 \
            --alignIntronMax 1 \
            --outWigType wiggle \
            --outWigStrand Stranded \
            --outWigNorm RPM
        mv Aligned.sortedByCoord.out.bam ${output_prefix}.bam
        mv ReadsPerGene.out.tab ${output_prefix}.tsv
        mv Log.final.out ${output_prefix}.Log.final.out
        python3  $(which make_star_qc.py) --quants ${output_prefix}.tsv \
                                          --star_log ${output_prefix}.Log.final.out \
                                          --output_filename ${output_prefix}_star_qc.json
    }

    output {
        File bam = glob("${output_prefix}.bam")[0]
        File tsv = glob("*.tsv")[0]
        File plus_strand_all_wig = glob("Signal.UniqueMultiple.str1.out.wig")[0]
        File minus_strand_all_wig = glob("Signal.UniqueMultiple.str2.out.wig")[0]
        File plus_strand_unique_wig = glob("Signal.Unique.str1.out.wig")[0]
        File minus_strand_unique_wig = glob("Signal.Unique.str2.out.wig")[0]
        File star_qc_json = glob("*_star_qc.json")[0]
        File star_qc_log = glob("star_qc.log")[0]
    }

    runtime {
        cpu: ncpus
        memory: "${ramGB} GB"
        disks: disk
    }
}

task wigtobigwig {
    File plus_strand_all_wig
    File minus_strand_all_wig
    File plus_strand_unique_wig
    File minus_strand_unique_wig
    File chrom_sizes
    String output_prefix
    Int ncpus
    Int ramGB
    String disk

    command {
        wigToBigWig ${plus_strand_all_wig} ${chrom_sizes} ${output_prefix}.signal.all.plus.bigWig
        wigToBigWig ${minus_strand_all_wig} ${chrom_sizes} ${output_prefix}.signal.all.minus.bigWig
        wigToBigWig ${plus_strand_unique_wig} ${chrom_sizes} ${output_prefix}.signal.unique.plus.bigWig
        wigToBigWig ${minus_strand_unique_wig} ${chrom_sizes} ${output_prefix}.signal.unique.minus.bigWig
    }

    output {
        File plus_strand_all_bigwig = glob("*.signal.all.plus.bigWig")[0]
        File minus_strand_all_bigwig = glob("*.signal.all.minus.bigWig")[0]
        File plus_strand_unique_bigwig = glob("*.signal.unique.plus.bigWig")[0]
        File minus_strand_unique_bigwig = glob("*.signal.unique.minus.bigWig")[0]
    }

    runtime {
        cpu: ncpus
        memory: "${ramGB} GB"
        disks: disk
    }
}

task spearman_correlation {
    Array[File] quants
    String output_filename

    command {
        python3 $(which calculate_correlation.py) --quants ${sep=' ' quants} --output_filename ${output_filename}
    }

    output {
        File spearman_json = glob("*spearman.json")[0]
    }

    runtime {
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 20 SSD"
    }
}

task bamtosam {
    File bamfile
    String output_sam
    Int ncpus
    Int ramGB
    String disk

    command {
        samtools view ${bamfile} > ${output_sam}
    }

    output {
        File samfile = glob("${output_sam}")[0]
    }

    runtime {
        cpu: ncpus
        memory: "${ramGB} GB"
        disks: disk
    }
}
