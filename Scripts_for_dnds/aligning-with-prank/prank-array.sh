#!/bin/bash -l

# Batch script to run a serial array job under SGE

# Request ten minutes of wallclock time (format hours:minutes:seconds).
# one job ran on frontend was 60 seconds
#$ -l h_rt=0:10:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
# /usr/bin/time --verbose prank gave Maximum resident set size (kbytes): 4656
# 5000 kb = 5 MB
# memory per node, so I am assuming per array task/orthogroup
#$ -l mem=10M

# Request 15 gigabyte of TMPDIR space (default is 10 GB)
# I have no idea about this
#$ -l tmpfs=15G

# Set up the job array.
# wc -l tmp/1to1-orthogroups
#$ -t 1-2516
#$ -tc 50

# Set the name of the job.
#$ -N PrankArrayJob

# Set the working directory to somewhere in your scratch space.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /lustre/scratch/scratch/ucfaeef/multispecies-dnds

orthogroup_name=$(sed -n "${SGE_TASK_ID}p" /lustre/scratch/scratch/ucfaeef/multispecies-dnds/tmp/1to1-orthogroups)

# prevent script to go on if part of pipeline fails
# set -eu -o pipeline

# Run the application.
echo "$JOB_NAME $SGE_TASK_ID"

/home/ucfaeef/programs/wasabi/binaries/prank/prank -d=/lustre/scratch/scratch/ucfaeef/multispecies-dnds/tmp/${orthogroup_name}-sequences.fa -o=/lustre/scratch/scratch/ucfaeef/multispecies-dnds/result/${orthogroup_name}-aligned_prank_output -f=paml -F -codon
