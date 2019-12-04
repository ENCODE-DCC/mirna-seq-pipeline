# ENCODE micro rna seq pipeline: cutadapt subworkflow
# Maintainer: Otto Jolanki

workflow cutadapt_wf {
    Array[File] fastqs_to_trim
    Array[File] five_prime_adapters
    File three_prime_adapters
    String output_prefix
    Int ncpus
    Int ramGB
    String disk

    scatter (i in range(length(fastqs_to_trim))) {
        call cutadapt { input:
            fastq = fastqs_to_trim[i],
            five_prime_adapters = five_prime_adapters[i],
            three_prime_adapters = three_prime_adapters,
            output_prefix = output_prefix,
            ncpus = ncpus,
            ramGB = ramGB,
            disk = disk,
        }
    }

    output {
        Array[File] no3ad_untrimmed_fastq = cutadapt.no3ad_untrimmed_fastq
        Array[File] no5ad_untrimmed_fastq = cutadapt.no5ad_untrimmed_fastq
        Array[File] too_short_fastq = cutadapt.too_short_fastq
        Array[File] trimmed_fastq = cutadapt.trimmed_fastq
    }
}

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