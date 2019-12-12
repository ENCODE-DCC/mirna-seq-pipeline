# Test workflow for cutadapt task in ENCODE micro rna seq pipeline

import "cutadapt_subworkflow.wdl" as cutadapt

workflow test_cutadapt {
    File fastq
    File five_prime_adapters
    File three_prime_adapters
    String output_prefix
    Int ncpus
    Int ramGB
    String disk

    call cutadapt.cutadapt { input:
        fastq = fastq,
        five_prime_adapters = five_prime_adapters,
        three_prime_adapters = three_prime_adapters,
        output_prefix = output_prefix,
        ncpus = ncpus,
        ramGB = ramGB,
        disk = disk,
    }
}
