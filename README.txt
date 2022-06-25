# Permutation of the DEGs (differentially expressed genes) in the nine wasp species.

run with `Wasp_permutation.pl Hypergeometric_Wasps 1000`

Need to give an input file name: DEGs file (1; taken from Supp table S1) and a number of permutations (2).


# Steps

The script looks up the UP or DW (down)regulation form Supp Table S1, for each fo the nine species to determine the total number fo DEGs. 

Then calculates if for each orthogroup there is a Trinity gene representative for each species.

Based on this information, the script creates a random list of "DEGS" of the same number from the list of potential orthologs.

Then it calculates the number of genes that should overlap between all 9, or al 8,,, etc, for the 9 species.
