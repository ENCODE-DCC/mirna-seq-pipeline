# HOWTO

Here are concrete instructions for running analyses on different platforms.
Before following these instructions, make sure you have completed installation and possible account setup detailed in [installation instructions](installation.md). These instructions show how to use Cromwell directly. Consider running the pipeline using [Caper](https://github.com/ENCODE-DCC/caper) which is more user friendly way.

# CONTENTS

## Running Analyses

[Google Cloud](howto.md#google-cloud)
[Local with Docker](howto.md#local-with-docker)
[SLURM](howto.md#slurm-singularity)
[Building STAR index](howto.md#building-star-index)

# RUNNING THE PIPELINE

## Google Cloud

Make sure you have completed the steps for installation and Google Cloud setup described in the [installation instructions](installation.md#google-cloud). The following assumes your Google Cloud project is `[YOUR_PROJECT]`, you have created a bucket `gs://[YOUR_BUCKET_NAME]`, and also directories `inputs`, `output` and `reference` in the bucket.
The goal is to run the pipeline with test data using Google Cloud Platform.

1. Launch a VM into your Google Cloud project, and connect to the instance.

2. Get the code and move into the code directory:

```bash
  git clone https://github.com/ENCODE-DCC/mirna-seq-pipeline.git
  cd mirna-seq-pipeline
```

3. Get the STAR index:

```bash
  curl https://storage.googleapis.com/circle_ci_test_data/mirna-seq-pipeline/star_index_mirna_chr19.tar.gz -o test_data/refs/star_index_mirna_chr19.tar.gz
```

4. Move the reference files and required adapters into your Google Cloud bucket (we will assume you have created the directories on your Google Cloud bucket):

```bash
  gsutil cp test_data/refs/* gs://[YOUR_BUCKET_NAME]/references/
  gsutil cp adapters/five_prime_adapter_set3.fasta gs://[YOUR_BUCKET_NAME]/references/
  gsutil cp adapters/five_prime_adapter_set4.fasta gs://[YOUR_BUCKET_NAME]/references/
  gsutil cp adapters/three_prime_adapter.fasta gs://[YOUR_BUCKET_NAME]/references/
```

5. Move the input files into your Google Cloud bucket (we will assume you have created the directories on your Google Cloud bucket):

```bash
  gsutil cp test_data/data/*.fastq.gz gs://[YOUR_BUCKET_NAME]/input_data/
```

6. Set up the `input.json`. Copy-paste the following onto your favorite text editor:

```
{
    "mirna_seq_pipeline.fastqs" : ["gs://[YOUR_BUCKET_NAME]/input_data/rep1ENCSR569QVM_chr19.fastq.gz", "gs://[YOUR_BUCKET_NAME]/input_data/rep2ENCSR569QVM_chr19.fastq.gz"],
    "mirna_seq_pipeline.five_prime_adapters" : ["gs://[YOUR_BUCKET_NAME]/references/five_prime_adapter_set3.fasta", "gs://[YOUR_BUCKET_NAME]/references/five_prime_adapter_set4.fasta"],
    "mirna_seq_pipeline.three_prime_adapters" : "gs://[YOUR_BUCKET_NAME]/references/three_prime_adapter.fasta",
    "mirna_seq_pipeline.star_index" : "gs://[YOUR_BUCKET_NAME]/references/star_index_mirna_chr19.tar.gz",
    "mirna_seq_pipeline.mirna_annotation" : "gs://[YOUR_BUCKET_NAME]/references/test.chr19.miRNA.gtf.gz",
    "mirna_seq_pipeline.chrom_sizes" : "gs://[YOUR_BUCKET_NAME]/references/test.chrom.sizes.tsv",
    "mirna_seq_pipeline.experiment_prefix" : "TEST_WORKFLOW_2REPS",
    "mirna_seq_pipeline.cutadapt_ncpus" : 1,
    "mirna_seq_pipeline.cutadapt_ramGB" : 2,
    "mirna_seq_pipeline.cutadapt_disk" : "local-disk 20 SSD",
    "mirna_seq_pipeline.star_ncpus" : 2,
    "mirna_seq_pipeline.star_ramGB" : 4,
    "mirna_seq_pipeline.star_disk" : "local-disk 20 SSD",
    "mirna_seq_pipeline.wigtobigwig_ncpus" : 1,
    "mirna_seq_pipeline.wigtobigwig_ramGB" : 2,
    "mirna_seq_pipeline.wigtobigwig_disk" : "local-disk 20 SSD"
}
```

Replace `[YOUR_BUCKET_NAME]` with the actual name of your Google Cloud bucket and save the file as `input.json`.

7. Run the pipeline:

```bash
    $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=google -Dbackend.providers.google.config.project=[YOUR_PROJECT] -Dbackend.providers.google.config.root=gs://[YOUR_BUCKET_NAME]/output cromwell-40.jar run mirna_seq_pipeline.wdl -i input.json -o workflow_opts/docker.json
```

8. You can observe the virtual machines spinning up on your Google Cloud Compute Engine dashboard, and after the pipeline finishes in few minutes you can the output in `gs://[YOUR_BUCKET_NAME]/output`.

## Local with Docker

Make sure you have completed the installation of docker, Java and Cromwell as described in the [installation instructions](installation.md). The goal is to run the pipeline with test data locally using Docker.

1. Get the code and move into the code directory:

```bash
  git clone https://github.com/ENCODE-DCC/mirna-seq-pipeline.git
  cd mirna-seq-pipeline
```

2. Get the STAR index:

```bash
  curl https://storage.googleapis.com/circle_ci_test_data/mirna-seq-pipeline/star_index_mirna_chr19.tar.gz -o test_data/refs/star_index_mirna_chr19.tar.gz
```

3. Run the pipeline:
```
  $ java -jar -Dconfig.file=backends/backend.conf cromwell-40.jar run mirna_seq_pipeline.wdl -i input.json -o workflow_opts/docker.json
```

## SLURM Singularity

For this example you need to have Singularity installed. For details see [installation instructions](installation.md). The goal is to run the pipeline with testdata using Singularity on a SLURM cluster. Login into your cluster first and then follow the instructions.
When running workflows on SLURM (or other) HPC clusters, use [Caper](https://github.com/ENCODE-DCC/caper), it takes care of backend configuration for you.

1. Get the code and move into the code directory:

```bash
  git clone https://github.com/ENCODE-DCC/mirna-seq-pipeline.git
  cd mirna-seq-pipeline
```

2. Get the STAR index:

```bash
  curl https://storage.googleapis.com/circle_ci_test_data/mirna-seq-pipeline/star_index_mirna_chr19.tar.gz -o test_data/refs/star_index_mirna_chr19.tar.gz
```

3. Build the singularity image for the pipeline. The following pulls the pipeline docker image, and uses that to construct the singularity image. The image will be stored in `~/.singularity`. It is bad practice to build images (or do any other intensive work) on login nodes. For this reason we will first invoke an interactive session on a different node by running `sdev` command, and building there (It will take few seconds to get back into the shell after running `sdev`).

```bash
  sdev
  mkdir -p ~/.singularity && cd ~/.singularity && SINGULARITY_CACHEDIR=~/.singularity SINGULARITY_PULLFOLDER=~/.singularity singularity pull --name mirna-seq-pipeline-v1.0.simg -F docker://quay.io/encode-dcc/mirna-seq-pipeline:v1.0
  exit #this takes you back to the login node
```

Note: If you want to store your inputs `/in/some/data/directory1`and `/in/some/data/directory2`you must edit `workflow_opts/singularity.json` in the following way:
```
{
    "default_runtime_attributes" : {
        "singularity_container" : "~/.singularity/mirna-seq-pipeline-v1.0.simg",
        "singularity_bindpath" : "~/, /in/some/data/directory1/, /in/some/data/directory2/"
    }
}
```

4. Install caper. Python 3.4.1 or newer is required.

```bash
  pip install caper
```

5. Follow [Caper configuration instructions](https://github.com/ENCODE-DCC/caper#configuration-file).

Note: In Caper configuration file, you will need to give a value to `--time` parameter by editing `slurm-extra-param` line. For example:
```
  slurm-extra-param=--time=01:00:00
```
to give one hour of runtime.

6. Edit the input file `test/test_workflow/test_workflow_2reps_input.json` so that all the input file paths are absolute.
For example replace `test_data/data/rep1ENCSR569QVM_chr19.fastq.gz` in fastq inputs with `[PATH-TO-REPO]/test_data/data/rep1ENCSR569QVM_chr19.fastq.gz`. You can find out the `[PATH-TO-REPO]` by running `pwd` command in the `mirna-seq-pipeline` directory.

7. Run the pipeline using Caper:

```bash
  caper run -i test/test_workflow/test_workflow_2reps_input.json -o workflow_opts/singularity.json -m metadata.json
```

## Building STAR index

If you want to build your own STAR index (maybe you are running analyses on another species, or want to use custom spikein sequences) you can use `generate_star_index.wdl` workflow. Format your `input.json` using the following example,
```
{
  "generate_STAR_index.reference_sequence" : "gs://mirna-seq-pipeline/full_sized_reference_files/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
  "generate_STAR_index.annotation" : "gs://mirna-seq-pipeline/full_sized_reference_files/gencode.v24.primary_assembly.annotation_fromPortal.gtf",
  "generate_STAR_index.output_filename" : "Star_index_GRCh38_from_primary_anno_to_portal.tar.gz",
  "generate_STAR_index.ncpus" : 16,
  "generate_STAR_index.ramGB" : 60,
  "generate_STAR_index.disks" : "local-disk 250 SSD"
}
```
and run on the platform of your choice, adapting from the previous examples, or using [Caper](https://github.com/ENCODE-DCC/caper).
