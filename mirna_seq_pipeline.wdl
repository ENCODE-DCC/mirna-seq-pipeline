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
    String experiment_prexix

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
            output_prefix = "rep"+(i+1)+experiment_prexix,
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

    }

    output {

    }

    runtime {
        cpu: ncpus
        memory: "${ramGB} GB"
        disks: disk
    }
}