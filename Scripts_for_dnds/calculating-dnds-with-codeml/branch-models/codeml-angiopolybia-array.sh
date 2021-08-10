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

# angiopolybia is foreground in model 2
species="angiopolybia"
spe="ang"

# create config file for null model (m0, dn/ds ratio without variable between lineages)
cp input/config-files/template_M0.ctl input/config-files/${spe}-${orthogroup}_M0.ctl

# update the orthogroup name
sed --in-place "s/OG/${orthogroup}/g" input/config-files/${spe}-${orthogroup}_M0.ctl

# update the species name
sed --in-place "s/agelaia/${species}/g" input/config-files/${spe}-${orthogroup}_M0.ctl

# run algorithm for model 0 (need full path of confid file)
/lustre/home/ucfaeef/programs/paml4.8/bin/codeml /lustre/home/ucfaeef/projects/multispecies-dnds/input/config-files/${spe}-${orthogroup}_M0.ctl

# create config file for model 2 (m2, dn/ds ratio with the species in foreground branch)
cp input/config-files/template_M2.ctl input/config-files/${spe}-${orthogroup}_M2.ctl

# update the orthogroup name
sed --in-place "s/OG/${orthogroup}/g" input/config-files/${spe}-${orthogroup}_M2.ctl

# update the species name
sed --in-place "s/agelaia/${species}/g" input/config-files/${spe}-${orthogroup}_M2.ctl

# run algorithm for model 2 (tree species if different)
/lustre/home/ucfaeef/programs/paml4.8/bin/codeml /lustre/home/ucfaeef/projects/multispecies-dnds/input/config-files/${spe}-${orthogroup}_M2.ctl

# collect variables
# obtain omega, number of parameters, log-likelihood

w_M0=`grep "omega" resultOnScratch/codeml/${species}-M0-mlc-${orthogroup} | sed "s/ //g" | cut -d "=" -f 2`

np_M0=`grep "lnL" resultOnScratch/codeml/${species}-M0-mlc-${orthogroup} | sed "s/ //g" | cut -d ")" -f 1 | cut -d ":" -f 3`

lnL_M0=`grep "lnL" resultOnScratch/codeml/${species}-M0-mlc-${orthogroup} | sed "s/ //g"  | cut -d "+" -f 1 | cut -d ":" -f 4`

np_M2=`grep "lnL" resultOnScratch/codeml/${species}-M2-mlc-${orthogroup} | sed "s/ //g" | cut -d ")" -f 1 | cut -d ":" -f 3`

lnL_M2=`grep "lnL" resultOnScratch/codeml/${species}-M2-mlc-${orthogroup} | sed "s/ //g"  | cut -d "+" -f 1 | cut -d ":" -f 4`

# - the omega and kappa values:
# the first omega is for branch #0, the background
background_w=`grep "w (dN/dS) for branches" resultOnScratch/codeml/${species}-M2-mlc-${orthogroup} | cut -d " " -f 6`

# the second omega is for branch #1, the foreground
foreground_w=`grep "w (dN/dS) for branches" resultOnScratch/codeml/${species}-M2-mlc-${orthogroup} | cut -d " " -f 7 `

kappa_M2=`grep "kappa" resultOnScratch/codeml/${species}-M2-mlc-${orthogroup} | cut -d " " -f 5`

# D = 2(nLL_0 - nLL)
D=`echo "2*($lnL_M2 - $lnL_M0)" | bc`

DF=`echo "$np_M2 - $np_M0" | bc`

chiTest=`/lustre/home/ucfaeef/programs/paml4.8/bin/chi2 $DF $D | cut -d " " -f 8 -s`


# collect all info for this round and save it in result/omega_results
# w_M0 np_M0 lnL_M0 np_M2 lnL_M2 background_w foreground_w kappaM1 TestStatistic DF chiTest
echo -e "${orthogroup}\t$w_M0\t$np_M0\t$lnL_M0\t$np_M2\t$lnL_M2\t$background_w\t$foreground_w\t$kappa_M2\t$D\t$DF\t$chiTest" >> result/omega_results/${species}
