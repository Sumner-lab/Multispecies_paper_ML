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
# I have no idea about this
#$ -l tmpfs=15G

# Set up the job array.
# wc -l tmp/aligned-orthogroups
#$ -t 1-1971

# Set the working directory
#$ -wd /lustre/home/ucfaeef/projects/multispecies-dnds

# set the orthogroup at start of array
orthogroup=$(sed -n "${SGE_TASK_ID}p" /lustre/home/ucfaeef/projects/multispecies-dnds/tmp/aligned-orthogroups)

# subset input .phy for only the species needed

seqlength=`head -n 1 resultOnScratch/${orthogroup}-aligned_prank_output.best.phy | sed "s/9 //g"`

theseLines=`echo $(((1 + seqlength  / 60 )))`
# change the header so that the first integer is 7
echo "7 ${seqlength}" > resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

# remove the sequences not needed (ie vespa and vespula)
grep -A ${theseLines} "Agelaia_cajannensis" resultOnScratch/${orthogroup}-aligned_prank_output.best.phy >> resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

grep -A ${theseLines} "Brachygastra_mellifica" resultOnScratch/${orthogroup}-aligned_prank_output.best.phy >> resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

grep -A ${theseLines} "Mischocyttarus_basimacula" resultOnScratch/${orthogroup}-aligned_prank_output.best.phy >> resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

grep -A ${theseLines} "Polistes_canadensis" resultOnScratch/${orthogroup}-aligned_prank_output.best.phy >> resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

grep -A ${theseLines} "Metapolybia_cingulata" resultOnScratch/${orthogroup}-aligned_prank_output.best.phy >> resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

grep -A ${theseLines} "Polybia_quadracincta" resultOnScratch/${orthogroup}-aligned_prank_output.best.phy >> resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

grep -A ${theseLines} "Angiopolybia_pallens" resultOnScratch/${orthogroup}-aligned_prank_output.best.phy >> resultOnScratch/polistine_${orthogroup}-aligned_prank_output.best.phy

