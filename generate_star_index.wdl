version 1.0

# ENCODE micro rna seq pipeline: generate star index

workflow generate_STAR_index {
    meta {
        author: "Otto Jolanki"
        version: "1.2.0"
        caper_docker: "encodedcc/mirna-seq-pipeline:1.2.0"
        caper_singularity: "docker://encodedcc/mirna-seq-pipeline:1.2.0"
    }

    input {
        # Reference .fasta
        File reference_sequence
        # Annotation .gtf
        File annotation
        # Output filename
        String output_filename
        Int ncpus
        Int ramGB
        String disks = "local-disk 50 SSD"
    }

    call generate_index { input:
        genome_file=reference_sequence,
        annotation_file=annotation,
        output_file=output_filename,
        ncpus=ncpus,
        ramGB=ramGB,
        disks=disks,
    }
}

task generate_index {
    input {
        File genome_file
        File annotation_file
        String output_file
        Int ncpus
        Int ramGB
        String disks
    }

    command {
        python3 $(which generate_star_index.py) \
            --ncpus ~{ncpus} \
            --annotation_file ~{annotation_file} \
            --genome_file ~{genome_file} \
            --output_file ~{output_file}
    }

    output {
        File star_index = output_file
        File star_log = "generate_star_index.log"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
    }
}
