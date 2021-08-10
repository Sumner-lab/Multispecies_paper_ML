#!/bin/bash -l

# Batch script to run a serial array job under SGE

# Request 20 minutes of wallclock time (format hours:minutes:seconds).
# one job ran on frontend was 2 min
#$ -l h_rt=0:20:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
# /usr/bin/time --verbose codeml gave Maximum resident set size (kbytes): 2276
# 5000 kb = 5 MB
# memory per node, so I am assuming per array task/orthogroup
#$ -l mem=10M

# Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=15G

# Set up the job array.
# wc -l tmp/aligned-orthogroups
#$ -t 1-1971

# Set the working directory
#$ -wd /lustre/home/ucfaeef/projects/multispecies-dnds

# set the orthogroup at start of array
orthogroup=$(sed -n "${SGE_TASK_ID}p" /lustre/home/ucfaeef/projects/multispecies-dnds/tmp/aligned-orthogroups)

# vespula is foreground in this test
species="vespula"

# create config file for null model (fix_omega = 1 (fixed))
cp input/branch-site-models/config-files/template_omegaFixed.ctl input/branch-site-models/config-files/${species}-${orthogroup}_omegaFixed.ctl

# update the orthogroup name
sed --in-place "s/OG/${orthogroup}/g" input/branch-site-models/config-files/${species}-${orthogroup}_omegaFixed.ctl

# update the species name in the tree file
sed --in-place "s/agelaia/${species}/g" input/branch-site-models/config-files/${species}-${orthogroup}_omegaFixed.ctl

# copy config file to input directory as codeml has a character limit
cp input/branch-site-models/config-files/${species}-${orthogroup}_omegaFixed.ctl input/.

# run algorithm for null model (need full path of config file)
/lustre/home/ucfaeef/programs/paml4.8/bin/codeml /lustre/home/ucfaeef/projects/multispecies-dnds/input/${species}-${orthogroup}_omegaFixed.ctl

# check result: resultOnScratch/codeml/branch-site-models/${species}-${orthogroup}-omegaFixed
# remove working copy (keeping things tidy)
rm -rf input/${species}-${orthogroup}_omegaFixed.ctl



# create config file for alternative model (fix_omega = 0 (not fixed))
cp input/branch-site-models/config-files/template_omegaNotFixed.ctl input/branch-site-models/config-files/${species}-${orthogroup}_omegaNotFixed.ctl

# update the species name in the tree file
sed --in-place "s/agelaia/${species}/g" input/branch-site-models/config-files/${species}-${orthogroup}_omegaNotFixed.ctl

# update the orthogroup name
sed --in-place "s/OG/${orthogroup}/g" input/branch-site-models/config-files/${species}-${orthogroup}_omegaNotFixed.ctl

# copy config file to input directory as codeml has a character limit
cp input/branch-site-models/config-files/${species}-${orthogroup}_omegaNotFixed.ctl input/.

# run algorithm for alternative model (fix_omega = 0 (not fixed))
/lustre/home/ucfaeef/programs/paml4.8/bin/codeml /lustre/home/ucfaeef/projects/multispecies-dnds/input/${species}-${orthogroup}_omegaNotFixed.ctl

# check result: resultOnScratch/codeml/branch-site-models/${species}-${orthogroup}-omegaNotFixed
# remove working copy (keeping things tidy)
rm -rf input/${species}-${orthogroup}_omegaNotFixed.ctl



