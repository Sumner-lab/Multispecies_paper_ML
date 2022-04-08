# Orthogroup refining and DnDs scripts

# (1) Refining the orthogroups

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

