# INSTALLATION

To run the pipeline you need to install following software. Running the pipeline on Google Cloud requires additional setup detailed below.

## Java 11 or newer

Java is required to run execution engine [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution).
To check which Java version you already have, run:
```bash
  $ java -version
```

## Docker

Pipeline code is packaged and distributed in Docker containers, and thus Docker installation is needed.
Follow instructions for [mac](https://docs.docker.com/docker-for-mac/install/) or [linux](https://docs.docker.com/install/linux/docker-ce/ubuntu/#upgrade-docker-after-using-the-convenience-script).

## Caper

For running the pipeline we recommend using [Caper](https://github.com/ENCODE-DCC/caper). Caper github page has extensive documentation on how to use caper on various platforms.

## croo

For organizing pipeline outputs we recommend using [croo](https://github.com/ENCODE-DCC/croo) that makes a nicely organized directory from the complicated output tree Cromwell defaults to. The configuration file for `croo` is named `output_definition.json` and can be found in the root of this repository.

## Singularity

It is possible to use Singularity instead of Docker. For details on how to run pipeline using singularity see [instructions](https://github.com/ENCODE-DCC/caper/blob/master/DETAILS.md#singularity).

## Google Cloud

To see how to run the pipeline on Google Cloud platform, see [instructions](https://github.com/ENCODE-DCC/caper/blob/master/scripts/gcp_caper_server/README.md).

## AWS

To see how to run the pipeline on AWS, see [instructions](https://github.com/ENCODE-DCC/caper/blob/master/scripts/aws_caper_server/README.md).

## HPC 

To see how to run the pipeline on a HPC cluster, see [instructions](https://github.com/ENCODE-DCC/caper#running-pipelines-on-hpcs)

## Truwl

You can run this pipeline on [truwl.com](https://truwl.com/workflows/library/ENCODE%20Micro%20RNA-seq%20pipeline/v1.2.0). This provides a web interface that allows you to define inputs and parameters, run the job on GCP, and monitor progress in a ready-to-go environment. To run it you will need to create an account on the platform then request early access by emailing [info@truwl.com](mailto:info@truwl.com) to get the right permissions. You can see an example case [here](https://truwl.com/workflows/library/ENCODE%20Micro%20RNA-seq%20pipeline/v1.2.0/instances/WF_844a7e.de.5f29). The example job (or other jobs) can be forked to pre-populate the inputs for your own job.

If you do not run the pipeline on Truwl, you can still share your use-case/job on the platform by getting in touch at [info@truwl.com](mailto:info@truwl.com) and providing your inputs.json file.
