version 1.0

# ENCODE micro rna seq pipeline: main pipeline
# Maintainer: Otto Jolanki

import "cutadapt_subworkflow.wdl" as cutadapt_sub

workflow mirna_seq_pipeline {
    meta {
        author: "Otto Jolanki"
        version: "1.2.0"
        caper_docker: "encodedcc/mirna-seq-pipeline:1.2.0"
        caper_singularity: "docker://encodedcc/mirna-seq-pipeline:1.2.0"
        croo_out_def: "https://storage.googleapis.com/encode-pipeline-output-definition/mirna.output_definition.json"
    }

    input {
        #Array containing the input fastq files
        Array[Array[File]] fastqs
        #Array containing Fasta files with 5' adapter sequence(s), in the same order as the fastqs
        Array[Array[File]] five_prime_adapters
        #Fasta file with 3' adapter sequence(s)
        File three_prime_adapters
        #tar.gz archive that contains star index
        File star_index
        #micro-rna annotations
        File mirna_annotation
        #tsv with chromosome sizes
        File chrom_sizes
        #Prefix for outputs (additionally replicate number will be added)
        String experiment_prefix
        Int cutadapt_ncpus
        Int cutadapt_ramGB
        String cutadapt_disk
        Int star_ncpus
        Int star_ramGB
        String star_disk
        Int wigtobigwig_ncpus
        Int wigtobigwig_ramGB
        String wigtobigwig_disk
    }

    scatter (i in range(length(fastqs))) {
        call cutadapt_sub.cutadapt_wf as cutadapt { input:
            fastqs_to_trim=fastqs[i],
            five_prime_adapters=five_prime_adapters[i],
            three_prime_adapters=three_prime_adapters,
            output_prefix="rep"+(i+1)+experiment_prefix,
            ncpus=cutadapt_ncpus,
            ramGB=cutadapt_ramGB,
            disk=cutadapt_disk,
        }
    }

    scatter (i in range(length(fastqs))) {
        call star { input:
            fastq=cutadapt.trimmed_fastq[i],
            index=star_index,
            annotation=mirna_annotation,
            output_prefix="rep"+(i+1)+experiment_prefix,
            ncpus=star_ncpus,
            ramGB=star_ramGB,
            disk=star_disk,
        }

        call wigtobigwig { input:
            plus_strand_all_wig=star.plus_strand_all_wig,
            minus_strand_all_wig=star.minus_strand_all_wig,
            plus_strand_unique_wig=star.plus_strand_unique_wig,
            minus_strand_unique_wig=star.minus_strand_unique_wig,
            chrom_sizes=chrom_sizes,
            output_prefix="rep"+(i+1)+experiment_prefix,
            ncpus=wigtobigwig_ncpus,
            ramGB=wigtobigwig_ramGB,
            disk=wigtobigwig_disk,
        }
    }

    #If there are exactly two replicates, calculate spearman correlation between quants in replicates
    if (length(fastqs) == 2) {
        call spearman_correlation { input:
            quants=star.tsv,
            output_filename=experiment_prefix+"_spearman.json",
            }
    }
}

#Task definitions

task star {
    input {
        File fastq
        File index
        File annotation
        String output_prefix
        Int ncpus
        Int ramGB
        String disk
    }

    command {
        tar -xzvf ~{index}
        gzip -cd ~{annotation} > anno.gtf
        STAR \
            --genomeDir out \
            --readFilesIn ~{fastq} \
            --sjdbGTFfile anno.gtf \
            --runThreadN ~{ncpus} \
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
        mv Aligned.sortedByCoord.out.bam ~{output_prefix}.bam
        mv ReadsPerGene.out.tab ~{output_prefix}.tsv
        mv Log.final.out ~{output_prefix}.Log.final.out
        python3  $(which make_star_qc.py) --quants ~{output_prefix}.tsv \
                                          --star_log ~{output_prefix}.Log.final.out \
                                          --output_filename ~{output_prefix}_star_qc.json
    }

    output {
        File bam = "~{output_prefix}.bam"
        File tsv = "~{output_prefix}.tsv"
        File plus_strand_all_wig = "Signal.UniqueMultiple.str1.out.wig"
        File minus_strand_all_wig = "Signal.UniqueMultiple.str2.out.wig"
        File plus_strand_unique_wig = "Signal.Unique.str1.out.wig"
        File minus_strand_unique_wig = "Signal.Unique.str2.out.wig"
        File star_qc_json = "~{output_prefix}_star_qc.json"
        File star_qc_log = "star_qc.log"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disk
    }
}

task wigtobigwig {
    input {
        File plus_strand_all_wig
        File minus_strand_all_wig
        File plus_strand_unique_wig
        File minus_strand_unique_wig
        File chrom_sizes
        String output_prefix
        Int ncpus
        Int ramGB
        String disk
    }

    command {
        wigToBigWig ~{plus_strand_all_wig} ~{chrom_sizes} ~{output_prefix}.signal.all.plus.bigWig
        wigToBigWig ~{minus_strand_all_wig} ~{chrom_sizes} ~{output_prefix}.signal.all.minus.bigWig
        wigToBigWig ~{plus_strand_unique_wig} ~{chrom_sizes} ~{output_prefix}.signal.unique.plus.bigWig
        wigToBigWig ~{minus_strand_unique_wig} ~{chrom_sizes} ~{output_prefix}.signal.unique.minus.bigWig
    }

    output {
        File plus_strand_all_bigwig = "~{output_prefix}.signal.all.plus.bigWig"
        File minus_strand_all_bigwig = "~{output_prefix}.signal.all.minus.bigWig"
        File plus_strand_unique_bigwig = "~{output_prefix}.signal.unique.plus.bigWig"
        File minus_strand_unique_bigwig = "~{output_prefix}.signal.unique.minus.bigWig"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disk
    }
}

task spearman_correlation {
    input {
        Array[File] quants
        String output_filename
    }

    command {
        python3 $(which calculate_correlation.py) --quants ~{sep=' ' quants} --output_filename ~{output_filename}
    }

    output {
        File spearman_json = output_filename 
    }

    runtime {
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 20 SSD"
    }
}

task bamtosam {
    input {
        File bamfile
        String output_sam
        Int ncpus
        Int ramGB
        String disk
    }

    command {
        samtools view ~{bamfile} > ~{output_sam}
    }

    output {
        File samfile = output_sam
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disk
    }
}
