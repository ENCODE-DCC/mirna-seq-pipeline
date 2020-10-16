version 1.0

# ENCODE micro rna seq pipeline: cutadapt subworkflow
# Maintainer: Otto Jolanki

workflow cutadapt_wf {
    input {
        Array[File] fastqs_to_trim
        Array[File] five_prime_adapters
        File three_prime_adapters
        String output_prefix
        Int ncpus
        Int ramGB
        String disk
    }

    scatter (i in range(length(fastqs_to_trim))) {
        call cutadapt { input:
            fastq=fastqs_to_trim[i],
            five_prime_adapters=five_prime_adapters[i],
            three_prime_adapters=three_prime_adapters,
            output_prefix=output_prefix,
            ncpus=ncpus,
            ramGB=ramGB,
            disk=disk,
        }
    }

        call merge_fastqs { input:
            no3ad_untrimmed_fastqs_=cutadapt.no3ad_untrimmed_fastq,
            no5ad_untrimmed_fastqs_=cutadapt.no5ad_untrimmed_fastq,
            too_short_fastqs_=cutadapt.too_short_fastq,
            trimmed_fastqs_=cutadapt.trimmed_fastq,
            output_prefix=output_prefix,
            ncpus=ncpus,
            ramGB=ramGB,
            disk=disk,
        }

    output {
        File no3ad_untrimmed_fastq = merge_fastqs.no3ad_untrimmed_fastq
        File no5ad_untrimmed_fastq = merge_fastqs.no5ad_untrimmed_fastq
        File too_short_fastq = merge_fastqs.too_short_fastq
        File trimmed_fastq = merge_fastqs.trimmed_fastq
    }
}

task cutadapt {
    input {
        File fastq
        File five_prime_adapters
        File three_prime_adapters
        String output_prefix
        Int ncpus
        Int ramGB
        String disk
    }

    command {
        cutadapt \
            -a file:~{three_prime_adapters} \
            -e 0.25 \
            --match-read-wildcards \
            --untrimmed-output=~{output_prefix + "_NO3AD.fastq"} \
            ~{fastq} \
            | cutadapt \
            -e 0.34 \
            --match-read-wildcards \
            --no-indels \
            -m 15 \
            -O 6 \
            -n 1 \
            -g file:~{five_prime_adapters} \
            --untrimmed-output=~{output_prefix + "_NO5AD.fastq"} \
            --too-short-output=~{output_prefix + "_SHORT_FAIL.fastq"} \
            - > ~{output_prefix + "_trim.fastq"}
    }

    output {
        File no3ad_untrimmed_fastq = "~{output_prefix}_NO3AD.fastq"
        File no5ad_untrimmed_fastq = "~{output_prefix}_NO5AD.fastq"
        File too_short_fastq = "~{output_prefix}_SHORT_FAIL.fastq"
        File trimmed_fastq = "~{output_prefix}_trim.fastq"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disk
    }
}

task merge_fastqs {
    input {
        Array[File] no3ad_untrimmed_fastqs_
        Array[File] no5ad_untrimmed_fastqs_
        Array[File] too_short_fastqs_
        Array[File] trimmed_fastqs_
        String output_prefix
        Int ncpus
        Int ramGB
        String disk
    }

    command {
        cat ~{sep=' ' no3ad_untrimmed_fastqs_} > ~{output_prefix}_merged_NO3AD.fastq
        cat ~{sep=' ' no5ad_untrimmed_fastqs_} > ~{output_prefix}_merged_NO5AD.fastq
        cat ~{sep=' ' too_short_fastqs_} > ~{output_prefix}_merged_SHORT_FAIL.fastq
        cat ~{sep=' ' trimmed_fastqs_} > ~{output_prefix}_merged_trim.fastq
    }

    output {
        File no3ad_untrimmed_fastq = "~{output_prefix}_merged_NO3AD.fastq"
        File no5ad_untrimmed_fastq = "~{output_prefix}_merged_NO5AD.fastq"
        File too_short_fastq = "~{output_prefix}_merged_SHORT_FAIL.fastq"
        File trimmed_fastq = "~{output_prefix}_merged_trim.fastq"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disk
    }
}
