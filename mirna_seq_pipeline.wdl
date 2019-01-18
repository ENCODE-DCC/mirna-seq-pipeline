# ENCODE micro rna seq pipeline: main pipeline
# Maintainer: Otto Jolanki

workflow mirna_seq_pipeline {
    #File inputs
    #Array containing the input fastq files
    Array[File] fastqs

    #Fasta file with 5' adapter sequence(s)
    File five_prime_adapters 
    
    #Fasta file with 3' adapter sequence(s)
    File three_prime_adapters

    #Prefix for outputs (additionally replicate number will be added)
    String experiment_prefix

    #Resources

    #Cutadapt
    Int cutadapt_ncpus
    Int cutadapt_ramGB
    String cutadapt_disk

    scatter (i in range(length(fastqs))) {
        call cutadapt { input:
            fastq = fastqs[i],
            five_prime_adapters = five_prime_adapters,
            three_prime_adapters = three_prime_adapters,
            output_prefix = "rep"+(i+1)+experiment_prefix,
            ncpus = cutadapt_ncpus,
            ramGB = cutadapt_ramGB,
            disk = cutadapt_disk,
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