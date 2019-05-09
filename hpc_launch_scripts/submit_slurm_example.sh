#!/bin/bash

# do not touch these settings
# number of tasks and nodes are fixed at 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1

# job name for pipeline
# this name will appear when you monitor jobs with "squeue -u $USER"
#SBATCH --job-name=MIRNA-TEST

# walltime for your job
# give long time enough to finish your pipeline
#SBATCH --time=8:00:00

# total amount of memory
#SBATCH --mem=20G

# max number of cpus for each pipeline
# should be >= NUM_CONCURRENT_TASK x "mirna.star" in input JSON file
# since star is a bottlenecking task in the pipeline
# "mirna_seq_pipeline.star_ncpus" is a number of cpus per replicate
#SBATCH --cpus-per-task=2

# email notification for job status
#SBATCH --mail-type=END,FAIL

# load java module if it exists
module load java || true

# use input JSON for a small test sample
INPUT=test/test_workflow/test_workflow_2reps_input.json

# pipeline metadata will be written in this file
PIPELINE_METADATA=metadata.json

# limit number of concurrent tasks
# we recommend to use a number of replicates here
# so that all replicates are processed in parellel at the same time.
# make sure that resource settings in your input JSON file
# are consistent with SBATCH resource settings (--mem, --cpus-per-task) 
# in this script
NUM_CONCURRENT_TASK=2
# run pipeline
# you can monitor your jobs with "squeue -u $USER"
java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=singularity \
-Dbackend.providers.singularity.config.concurrent-job-limit=${NUM_CONCURRENT_TASK} \
$HOME/cromwell-35.jar run mirna_seq_pipeline.wdl -i ${INPUT} -o workflow_opts/singularity.json -m ${PIPELINE_METADATA}
