#!/bin/bash

set -Eeuo pipefail

# Nextflow Version
NXF_VER=20.10.0

# Nextflow Configuration File
NXF_CONFIG=nextflow.config

# Workflow to Run (e.g. GitHub repository)
WORKFLOW_REPO=main.nf #shashidhar22/airr_seq_pipelines

# Queue to use for exeution
QUEUE=campus-new

# Define scratch directory
SCRATCH=/fh/scratch/delete90/warren_h/cfar/

# Define params file
PARAMS=params/cfar.yaml

# Load the Nextflow module (if running on rhino/gizmo)
ml Nextflow

# Load the Singularity module (if running on rhino/gizmo with Singularity)
ml Singularity

# Load Conda module (if running on rhino/gizmo with Conda)
ml Anaconda3
# Make sure that the singularity executables are in the PATH
export PATH=$SINGULARITYROOT/bin/:$PATH

# Run the workflow
NXF_VER=$NXF_VER \
nextflow \
    run \
    -c ${NXF_CONFIG} \
    ${WORKFLOW_REPO} \
    -entry gliphRunner \
    -w  ${SCRATCH} \
    -with-tower \
    -params-file ${PARAMS} \
    -resume