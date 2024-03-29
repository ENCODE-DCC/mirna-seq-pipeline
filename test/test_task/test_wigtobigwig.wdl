version 1.0

# Test workflow for wigtobigwig task in ENCODE micro rna seq pipeline

import "../../mirna_seq_pipeline.wdl" as mirna

workflow test_wigtobigwig {
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
        String docker
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": ""
    }

    call mirna.wigtobigwig { input:
        plus_strand_all_wig=plus_strand_all_wig,
        minus_strand_all_wig=minus_strand_all_wig,
        plus_strand_unique_wig=plus_strand_unique_wig,
        minus_strand_unique_wig=minus_strand_unique_wig,
        chrom_sizes=chrom_sizes,
        output_prefix =output_prefix,
        ncpus=ncpus,
        ramGB=ramGB,
        disk=disk,
        runtime_environment=runtime_environment,
    }
}
