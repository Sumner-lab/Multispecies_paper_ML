# Multispecies_paper_ML

Scripts used for Multispecies paper (Wyatt et al. Unpublished)

Master.ML.pl is a perl wrapper that creates all the R scripts necessary to run the svm analysis. 

# SET UP

The Master script requires the data to be in a specific folder ordering system with:

**1.** Scripts in the following folder:

MAIN/scripts/template_scripts/

**2.** **Species** gene expression data and Trinity assemblies in :

MAIN/DATA/Experimental_data_merged/**Polybia_quadracincta**/

MAIN/DATA/Experimental_data_merged/**Polistes_canadensis**/

In each of these folders you need:

a Trinity assembly in /Pool/

a Queen and Worker folder with RSEM gene expression counts based on the Trinity assembly (/RSEM.isoforms.results)

Each species should have, for example: 

MAIN/DATA/Experimental_data_merged/Polistes_canadensis/**Pool**/Pool_trinity.fnn

MAIN/DATA/Experimental_data_merged/Polistes_canadensis/**Queen**/RSEM.isoforms.results

MAIN/DATA/Experimental_data_merged/Polistes_canadensis/**Worker**/RSEM.isoforms.results

**3.** Gene orthology information from Orthofinder should be found in:

MAIN/Orthofinder/Orthogroups.copy.noMac.tsv

**4.** Folder to be intialised:

mkdir MAIN/FIGURES

mkdir MAIN/DEGS


# **Dependencies on UCL myriad cluster**:

module unload compilers/intel/2018/update3

module unload mpi/intel/2018/update3/intel

module unload java/1.8.0_45

module add r/recommended


**R must have the following packages installed**:

library(tximport)

library(tximportData)

library(edgeR)

library(seqinr)

library(pheatmap)

library(stringr)

library(e1071)

library(probsvm)



# **Running the script**

Usage: Master.ML.pl -j <Version folder name> -f <Foreground species list (comma sep)>  -b <Background species list (comma sep)> -filter <matrix_run7.filter.scale>

options:    **-r**     Run the scripts through R (default = OFF).
            **-e**     Choose expression data folder (default= \"Experimental_data_merged\").
            **-orth**  Choose Orthofinder file to use in this analysis (default= \"DATA/Orthofinder/Orthofinder_24.5.2020/Orthogroups.tsv\"). 
            **-CPM**   Choose the cut off of minimal expression to be considered. Default =1.  

Must be the top level folder (e.g. MAIN), which should have folders DATA and scripts and FIGURES.
