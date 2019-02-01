# REFERENCE

This document contains more detailed information on the inputs, outputs and the software.

# CONTENTS

[Software](reference.md#software)  
[Inputs](reference.md#inputs)  
[Outputs](reference.md#outputs)

## Software

### Ubuntu 16.04

The pipeline docker image is based on [Ubuntu base image](https://hub.docker.com/_/ubuntu/) version `16.04`.

### Python 3.5.2

Python parts of the pipeline are run using [Python 3.5.2](https://www.python.org/download/releases/3.5.2/) that ships with Ubuntu 16.04.

### STAR 2.5.1b

Alignment, quantitation and generation of the .wig files is done using [STAR 2.5.1b](https://github.com/alexdobin/STAR/releases/tag/2.5.1b). For detailed description of the software see [Article by Dobin et al](https://www.ncbi.nlm.nih.gov/pubmed/23104886). Multiple versions of the software have been released since writing the article.

## Inputs

A typical input json file looks like this:

```
{
    "mirna_seq_pipeline.fastqs" : ["gs://mirna-seq-pipeline/inputs/ENCSR569QVM/ENCFF119JCH_rep1.fastq.gz","gs://mirna-seq-pipeline/inputs/ENCSR569QVM/ENCFF756CSN_rep2.fastq.gz"],
    "mirna_seq_pipeline.five_prime_adapters" : ["gs://mirna-seq-pipeline/full_sized_reference_files/adapters/five_prime_adapter_set3.fasta","gs://mirna-seq-pipeline/full_sized_reference_files/adapters/five_prime_adapter_set4.fasta"],
    "mirna_seq_pipeline.three_prime_adapters" : "gs://mirna-seq-pipeline/full_sized_reference_files/adapters/three_prime_adapter.fasta",
    "mirna_seq_pipeline.star_index" : "gs://mirna-seq-pipeline/full_sized_reference_files/Star_index_GRCh38.tar.gz",
    "mirna_seq_pipeline.mirna_annotation" : "gs://mirna-seq-pipeline/full_sized_reference_files/ENCFF628BVT_mirna_anno.gtf",
    "mirna_seq_pipeline.chrom_sizes" : "gs://mirna-seq-pipeline/full_sized_reference_files/GRCh38_EBV.chrom.sizes.tsv",
    "mirna_seq_pipeline.experiment_prefix" : "ENCSR569QVM",
    "mirna_seq_pipeline.cutadapt_ncpus" : "2",
    "mirna_seq_pipeline.cutadapt_ramGB" : "7",
    "mirna_seq_pipeline.cutadapt_disk" : "local-disk 200 SSD",
    "mirna_seq_pipeline.star_ncpus" : "16",
    "mirna_seq_pipeline.star_ramGB" : "60",
    "mirna_seq_pipeline.star_disk" : "local-disk 200 SSD",
    "mirna_seq_pipeline.wigtobigwig_ncpus" : "4",
    "mirna_seq_pipeline.wigtobigwig_ramGB" : "26",
    "mirna_seq_pipeline.wigtobigwig_disk" : "local-disk 200 SSD"
}
```

Following elaborates the meaning of each line in the input file.

* `mirna_seq_pipeline.fastqs` Is a list of input fastq files, one for each replicate.
* `mirna_seq_pipeline.five_prime_adapters` Is a list of 5' adapter fasta files, one for each replicate. Note: order of this list should correspond to the order of `mirna_seq_pipeline.fastqs`.
* `mirna_seq_pipeline.three_prime_adapters` Is the fasta file containing 3' adapters. Same set is used for all the replicates.
* `mirna_seq_pipeline.star_index` Is the gzipped tar archive that contains STAR index. 
* `mirna_seq_pipeline.mirna_annotation` Is the .gtf file containing miRNA annotations.
* `mirna_seq_pipeline.chrom_sizes` Is .tsv that contains chromosome sizes.
* `mirna_seq_pipeline.experiment_prefix` Is prefix that will be added to important output filenames.

#### Example: 
    
    Assume the `mmirna_seq_pipeline.experiment_prefix` is "FOO_BAR_BAZ". Outputs from replicate 1 will get a prefix rep1FOO_BAR_BAZ, outputs from replicate 2 get a prefix rep2FOO_BAR_BAZ etc.

Following inputs define the computational resources given to the pipeline tasks.

* `mirna_seq_pipeline.cutadapt_ncpus` Is the number of cores given to the cutadapt task.
* `mirna_seq_pipeline.cutadapt_ramGB` Is the amount of ram in gigabytes given to the cutadapt task.
* `mirna_seq_pipeline.cutadapt_disk` Is the amount of disk space in gigabytes to the cutadapt task.
* `mirna_seq_pipeline.star_ncpus` Is the number of cores given to the star task.
* `mirna_seq_pipeline.star_ramGB` Is the amount of ram in gigabytes given to the star task.
* `mirna_seq_pipeline.star_disk` Is the amount of disk space in gigabytes to the star task.
* `mirna_seq_pipeline.wigtobigwig_ncpus` Is the number of cores given to the wigtobigwig task.
* `mirna_seq_pipeline.wigtobigwig_ramGB` Is the amount of ram in gigabytes given to the wigtobigwig task.
* `mirna_seq_pipeline.wigtobigwig_disk` Is the amount of disk space in gigabytes to the wigtobigwig task.

## Outputs