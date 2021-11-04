version 1.0

# Test workflow for cutadapt task in ENCODE micro rna seq pipeline

import "../../cutadapt_subworkflow.wdl" as cutadapt

workflow test_cutadapt {
    input {
        File fastq
        File five_prime_adapters
        File three_prime_adapters
        String output_prefix
        Int ncpus
        Int ramGB
        String disk
        String docker
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": ""
    }

    call cutadapt.cutadapt { input:
        fastq=fastq,
        five_prime_adapters=five_prime_adapters,
        three_prime_adapters=three_prime_adapters,
        output_prefix=output_prefix,
        ncpus=ncpus,
        ramGB=ramGB,
        disk=disk,
        runtime_environment=runtime_environment,
    }
}
