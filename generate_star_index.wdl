version 1.0

# ENCODE micro rna seq pipeline: generate star index

import "wdl/structs/runtime.wdl"

workflow generate_STAR_index {
    meta {
        author: "Otto Jolanki"
        version: "1.2.2"
        caper_docker: "encodedcc/mirna-seq-pipeline:1.2.2"
        caper_singularity: "docker://encodedcc/mirna-seq-pipeline:1.2.2"
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
        String docker = "encodedcc/mirna-seq-pipeline:1.2.2"
        String singularity = "docker://encodedcc/mirna-seq-pipeline:1.2.2"
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": singularity
    }
                                

    call generate_index { input:
        genome_file=reference_sequence,
        annotation_file=annotation,
        output_file=output_filename,
        ncpus=ncpus,
        ramGB=ramGB,
        disks=disks,
        runtime_environment=runtime_environment,
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
        RuntimeEnvironment runtime_environment
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
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
