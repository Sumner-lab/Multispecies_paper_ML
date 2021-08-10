#!/bin/bash -l

# Batch script to run a serial array job under SGE

# Request 20 minutes of wallclock time (format hours:minutes:seconds).
# one job ran on frontend was 2 min
#$ -l h_rt=1:0:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
# memory per node, so I am assuming per array task/orthogroup
#$ -l tmem=10M

# Set up the job array.
# wc -l input/aligned-orthogroups
#$ -t 1-1971

# Set the working directory
#$ -wd /SAN/ugi/VespaCrabro/multispecies 

export TMPDIR=/SAN/ugi/VespaCrabro/tmp

# set the orthogroup at start of array
orthogroup=$(sed -n "${SGE_TASK_ID}p" input/aligned-orthogroups)

export PATH=/share/apps/perl/bin:$PATH

num_bp=$( sed -n '1,1p' resultOnScratch/phylip-files/${orthogroup}-aligned_prank_output.best.phy | sed 's/..* //' )

perl subset_species.pl resultOnScratch/phylip-files/${orthogroup}-aligned_prank_output.best.phy subsp.txt $num_bp

mv ${orthogroup}-aligned_prank_output.out.aln resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

