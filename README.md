# Multispecies_paper_ML
Scripts used for Multispecies paper (Wyatt et al. Unpublished)

Master.ML.pl is a perl wrapper that creates all the R scripts necessary to run the svm analysis. 

The Master script requires the data to be in a specific folder ordering system with:

1. Scripts in the following folder:
MAIN/scripts/template_scripts/

2. Species gene expression data and Trinity assemblies in :
MAIN/DATA/Experimental_data_merged/\bPolybia_quadracincta\b/
MAIN/DATA/Experimental_data_merged/Polistes_canadensis/

In each of these folders you need:
a Trinity assembly in /Pool/
a Queen and Worker folder with RSEM gene expression counts based on the Trinity assembly (/RSEM.isoforms.results)

Each species should have, for example: 
MAIN/DATA/Experimental_data_merged/Polistes_canadensis/Pool/Pool_trinity.fnn
MAIN/DATA/Experimental_data_merged/Polistes_canadensis/Queen/RSEM.isoforms.results
MAIN/DATA/Experimental_data_merged/Polistes_canadensis/Worker/RSEM.isoforms.results

Gene orthology informatio from Orthofinder should be found in:
MAIN/Orthofinder/Orthogroups.copy.noMac.tsv
