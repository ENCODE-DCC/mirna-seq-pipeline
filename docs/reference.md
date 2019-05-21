# REFERENCE

This document contains more detailed information on the inputs, outputs and the software.

# CONTENTS

[Software](reference.md#software)  
[Recommended Software](reference.md#recommended-software)  
[Inputs](reference.md#inputs)  
[Reference files](reference.md#getting-reference-files)  
[Outputs](reference.md#outputs)  
[Output organizer](reference.md#cromweller-output-organizer)

## Software

### Ubuntu 16.04

The pipeline docker image is based on [Ubuntu base image](https://hub.docker.com/_/ubuntu/) version `16.04`.

### Python 3.5.2

Python parts of the pipeline are run using [Python 3.5.2](https://www.python.org/download/releases/3.5.2/) that ships with Ubuntu 16.04.

### Cutadapt 1.7.1

For trimming the adapters we use [cutadapt version 1.7.1](https://cutadapt.readthedocs.io/en/stable/). For publication detailing the software, see [Article by Marcel Martin](http://journal.embnet.org/index.php/embnetjournal/article/view/200).

### WigtoBigWig

Conversion from .wig format to .bigWig format is done using WigToBigWig, downloaded from [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig).

### STAR 2.5.1b

Alignment, quantitation and generation of the .wig files is done using [STAR 2.5.1b](https://github.com/alexdobin/STAR/releases/tag/2.5.1b). For detailed description of the software see [Article by Dobin et al](https://www.ncbi.nlm.nih.gov/pubmed/23104886). Multiple versions of the software have been released since writing the article.

## Recommended Software

To make Cromwell usage more pleasant, we have developed [Caper](https://github.com/ENCODE-DCC/caper) that provides highly enhanced experience when running pipelines, and [croo](https://github.com/ENCODE-DCC/croo) that helps organizing Cromwell outputs. Croo configuration file for this pipeline is `output_definition.json` and it is located in the root of this repository.

## Inputs

### Getting reference files

MiRNA annotations are available in [The ENCODE Portal](https://www.encodeproject.org/files/ENCFF628BVT/). Corresponding STAR Index is also available via [The ENCODE Portal](https://www.encodeproject.org/files/ENCFF033AVX/). It is also possible to build STAR index from reference by using the `generate_star_index.wdl`.

### Pipeline input file

A typical input json file looks like this:

```
{
    "mirna_seq_pipeline.fastqs" : ["gs://mirna-seq-pipeline/inputs/ENCSR569QVM/ENCFF119JCH_rep1.fastq.gz","gs://mirna-seq-pipeline/inputs/ENCSR569QVM/ENCFF756CSN_rep2.fastq.gz"],
    "mirna_seq_pipeline.five_prime_adapters" : ["gs://mirna-seq-pipeline/full_sized_reference_files/adapters/five_prime_adapter_set3.fasta","gs://mirna-seq-pipeline/full_sized_reference_files/adapters/five_prime_adapter_set4.fasta"],
    "mirna_seq_pipeline.three_prime_adapters" : "gs://mirna-seq-pipeline/full_sized_reference_files/adapters/three_prime_adapter.fasta",
    "mirna_seq_pipeline.star_index" : "gs://mirna-seq-pipeline/full_sized_reference_files/Star_index_GRCh38.tar.gz",
    "mirna_seq_pipeline.mirna_annotation" : "gs://mirna-seq-pipeline/full_sized_reference_files/ENCFF628BVT_mirna_anno.gtf.gz",
    "mirna_seq_pipeline.chrom_sizes" : "gs://mirna-seq-pipeline/full_sized_reference_files/GRCh38_EBV.chrom.sizes.tsv",
    "mirna_seq_pipeline.experiment_prefix" : "ENCSR569QVM",
    "mirna_seq_pipeline.cutadapt_ncpus" : "2",
    "mirna_seq_pipeline.cutadapt_ramGB" : "7",
    "mirna_seq_pipeline.cutadapt_disk" : "local-disk 200 SSD",
    "mirna_seq_pipeline.star_ncpus" : "16",
    "mirna_seq_pipeline.star_ramGB" : "60",
    "mirna_seq_pipeline.star_disk" : "local-disk 200 SSD",
    "mirna_seq_pipeline.wigtobigwig_ncpus" : "2",
    "mirna_seq_pipeline.wigtobigwig_ramGB" : "7",
    "mirna_seq_pipeline.wigtobigwig_disk" : "local-disk 200 SSD"
}
```

Following elaborates the meaning of each line in the input file.

* `mirna_seq_pipeline.fastqs` Is a list of input fastq files, one for each replicate.
* `mirna_seq_pipeline.five_prime_adapters` Is a list of 5' adapter fasta files, one for each replicate. Note: order of this list should correspond to the order of `mirna_seq_pipeline.fastqs`.
* `mirna_seq_pipeline.three_prime_adapters` Is the fasta file containing 3' adapters. Same set is used for all the replicates.
* `mirna_seq_pipeline.star_index` Is the gzipped tar archive that contains STAR index. GRCh38 based on Gencode V24 version is available for download on [The ENCODE Portal](https://www.encodeproject.org/files/ENCFF033AVX/).
* `mirna_seq_pipeline.mirna_annotation` Is the gzipped .gtf(.gz) file containing miRNA annotations. Gencode V24 version available for download on [The ENCODE Portal](https://www.encodeproject.org/files/ENCFF628BVT/).
* `mirna_seq_pipeline.chrom_sizes` Is .tsv that contains chromosome sizes. This can be downloaded from [The ENCODE Portal](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/)
* `mirna_seq_pipeline.experiment_prefix` Is prefix that will be added to important output filenames.

#### Example: 
    
    Assume the `mirna_seq_pipeline.experiment_prefix` is "FOO_BAR_BAZ". Outputs from replicate 1 will get a prefix rep1FOO_BAR_BAZ, outputs from replicate 2 get a prefix rep2FOO_BAR_BAZ etc.

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

#### Cromwell output directory structure

Cromwell will store outputs for each task under directory cromwell-executions/[WORKFLOW_ID]/call-[TASK_NAME]/shard-[IDX]. For all tasks [IDX] means a zero-based index for each replicate. Most humans find the Cromwell output structure unwieldy and confusing. If you are like most humans, see [the next section](reference.md#cromwell-output-organizer)

#### Cromwell output organizer

To make the experience of looking at the outputs we recommend using [the cromwell output organizer](https://github.com/ENCODE-DCC/cromwell_output_organizer). The output definition file (`output_definition.json`) that the output organizer requires is provided as a part of this github repo.

#### Task Cutadapt

* `no3ad_untrimmed_fastq` .fastq with the reads in which the 3' adapter was not found.
* `no5ad_untrimmed_fastq` .fastq with the reads in which the 5' adapter was not found.
* `too_short_fastq` .fastq with the reads that are too short.
* `trimmed_fastq` .fastq with successfully trimmed reads.

#### Task Star

* `bam` .bam with the alignments
* `tsv` .tsv with quantitations
* `plus_strand_all_wig` Wiggle from all plus strand reads. Intermediate file used for signal track generation.
* `minus_strand_all_wig` Wiggle from all minus strand reads. Intermediate file used for signal track generation.
* `plus_strand_unique_wig` Wiggle from uniquely mapping plus strand reads. Intermediate file used for signal track generation.
* `minus_strand_unique_wig` Wiggle from uniquely mapping minus strand reads. Intermediate file used for signal track generation.
* `star_qc_json` .json with alignment QC metrics.

#### Task WigToBigWig

* `plus_strand_all_bigwig` Signal track (bigwig) from all plus strand reads.
* `minus_strand_all_bigwig` Signal track (bigwig) from all minus strand reads.
* `plus_strand_unique_bigwig` Signal track (bigwig) from uniquely mapping plus strand reads.
* `minus_strand_unique_bigwig` Signal track (bigwig) from uniquely mapping minus strand reads.

#### Task Spearman_correlation (when there are exactly 2 replicates)

* `spearman_json` .json file containing spearman correlation metric between the replicates
