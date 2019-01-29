# Test workflow for star task in ENCODE micro rna seq pipeline

import "../../mirna_seq_pipeline.wdl" as mirna

workflow test_star {
    File fastq
    File index
    File mirna_annotation 
    String output_prefix
    Int ncpus
    Int ramGB
    String disk

    call mirna.star {input:
        fastq = fastq,
        index = index,
        annotation = mirna_annotation,
        output_prefix = output_prefix,
        ncpus = ncpus,
        ramGB = ramGB,
        disk = disk
    }
}