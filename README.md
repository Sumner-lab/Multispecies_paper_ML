# Multispecies_paper_ML

Scripts used for Multispecies paper (Wyatt et al. Unpublished). These scripts allow a bespoke way to do svm machine learning on non-genome guided RNAseq data (Trinity), across multiple species. This process is not automated, and requires a lengthy process of acquiring Trinity assemblies for your species of interest, then compiling the data in a specific data structure (folder organisation- see later) to make the script perform without error. Please contact the authors if you wish to embark on a similar analysis (before trying out the various scripts), as likely the scripts would have to be changed for another analysis.

The 3 folders contain the scripts used to:
1. Run the SVM analysis (with Trinity assemblies)
2. Refine the orthogroup lists, to reduce isoform representatives and to allow NAs across the species. 
3. Run the DnDs analysis


# (1) SVM

Master.ML.pl is a perl wrapper that creates all the R scripts necessary to run the svm analysis. 

# Set Up

The Master script requires the data to be in a specific folder organisation:

**1.** Template R scripts in the following folder:
```unix
MAIN/scripts/template_scripts/
```
R scripts found in : [template_scripts](https://github.com/Sumner-lab/Multispecies_paper_ML/tree/master/template_scripts)

**2.** **Species** gene expression data and Trinity assemblies in :
```unix
MAIN/DATA/Experimental_data_merged/Polybia_quadracincta/

MAIN/DATA/Experimental_data_merged/Polistes_canadensis/
```
In each of these folders you need:

a Trinity assembly in /Pool/

a Queen and Worker folder with RSEM gene expression counts based on the Trinity assembly (/RSEM.isoforms.results)

Each species should have, for example: 
```unix
MAIN/DATA/Experimental_data_merged/Polistes_canadensis/Pool/Pool_trinity.fnn

MAIN/DATA/Experimental_data_merged/Polistes_canadensis/Queen/RSEM.isoforms.results

MAIN/DATA/Experimental_data_merged/Polistes_canadensis/Worker/RSEM.isoforms.results
```
**3.** Gene orthology information from Orthofinder should be found in:

```unix
MAIN/Orthofinder/Orthogroups.copy.noMac.tsv
```

Orthology must be determined using Orthofinder, explained in methods of paper, using all the species in your analysis (based on Trinity gene names: TRINITY_DN16291_c0_g1). This involved calculation of protein fasta sequences of your species and running Orthofinder as default.

**4.** Folder to be intialised:

```unix
mkdir MAIN/FIGURES
mkdir MAIN/DEGS
```


# **Dependencies on UCL myriad cluster**:

It is recommended to run this script on a High performance cluster, as the script uses alot of RAM, but can be run locally if have a powerful machine. These settings will be different depending on your cluster.

For UCL myriad we need to load R to run the script.
```unix
module unload compilers/intel/2018/update3
module unload mpi/intel/2018/update3/intel
module unload java/1.8.0_45
module add r/recommended
```

# **R must have the following packages installed**:
```R
library(tximport)
library(tximportData)
library(edgeR)
library(seqinr)
library(pheatmap)
library(stringr)
library(e1071)
library(probsvm)
```


# **Running the script**
```unix
Usage: Master.ML.pl -j <Version folder name> -f <Foreground species list (comma sep)>  -b <Background species list (comma sep)> 

options:    -r     Run the scripts through R (default = OFF).
            -e     Choose expression data folder (default= \"Experimental_data_merged\").
            -orth  Choose Orthofinder file to use in this analysis (default= \"DATA/Orthofinder/Orthofinder_24.5.2020/Orthogroups.tsv\"). 
            -CPM   Choose the cut off of minimal expression to be considered. Default =1.  
            -Gamma Choose the gamma parameters to search. Default = 10^(-7:-5)
            -Cost  Choose the cost parameters to search. Default = 2^(3:5)

Must be the top level folder (e.g. Multispecies_FILES), which should have folders DATA and scripts and FIGURES.

*** Questions: 
    Chris Wyatt (former post-doctorate) or Seirian Sumner
```

Must be the top level folder (e.g. MAIN), which should have folders DATA and scripts and FIGURES.

A typical example of the script would be :
```unix
perl Master.ML.pl -j TOP5_Test_Angiopolybia_pallens_3MER -CPM 10 -f Angiopolybia_pallens -b Vespula_vulgaris,Vespa_crabro,Metapolybia_cingulata,Polybia_quadracincta -e Experimental_data_merged -orth Orthofinder/Orthogroups.copy.noMac.tsv.filt2.csv -r
```



# (2) Refining the orthogroups

i) Refine_orthogroups.pl

This step is optional to increase the number of orthogroups investigated in the SVM. To allow the SVM matrix to contain a single isoform of each orthogroup, for each species. 


This script requires the Orthogroups.csv file from Orthofinder. Plus the names in the header should match the file system order mentioned previously. So that on lines 31/32, you can find the expression files for each species (as shown below, lines 31 and 32):
```
	my $input_Q_rsem="../DATA\/Experimental_data_merged\/$species\/Queen\/RSEM.isoforms.results";
	my $input_W_rsem="../DATA\/Experimental_data_merged\/$species\/Worker\/RSEM.isoforms.results";
```

ii) Zeros_orthogroups.v3.pl

This script imputes a fake entry for NAs in the orthogroup.csv file. This inputs a fake gene name of "TRINITY_UK00000_c0_g1". Then the gene expression counts are set to 1 for queen and 1 for worker. Again, this allows more genes to be considered in the SVM analysis.

iii) get_fasta_largest_isoform.TrinityMS.pl 

This script find the largest isoform representative of a gene to be used in the analysis. Removing splice variants of the same gene is necessary in order to compare just a single gene representative for each orthogroup, for each species.





# (3) DnDs

# Positive selection in nine species (wasp social spectrum)

Using the nucleotide fasta of 1:1 orthogroups for each species (see WHATIDID.txt for details), we aligned sequences for the same orthogroup (see scripts in aligning-with-prank). We then calculated dnds ratios (see scripts in calculating-dnds-with-codeml) for all species and for the Polistines only (see scripts in subsetting-phylip-files-for-polistines). We analysed the results with R (see scripts in analysing-results-with-R).

