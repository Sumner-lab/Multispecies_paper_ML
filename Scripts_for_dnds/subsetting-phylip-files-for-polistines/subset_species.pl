#!/usr/bin/perl
​
use strict;
use warnings;
 
## Script that reads a PHYLIP alignment and a text file 
## with the species names to be subset. It outputs a new 
## alignment with the species listed in this text file.
##
## Usage:
## <path_to_script>/Get_partitions.pl <path_to_aln1> <path_to_text_file> <num_bp>
##
## Contact information: <sandra.ac93@gmail.com>
​
#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!!
#IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
open (INFILE2, "<$ARGV[1]") or die "Cannot open $ARGV[0] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/\..*//;
$outname =~ s/..*\///;
print "\nParsing alignment ".$outname." ... ... \n";
my $outname2 = $outname."_out.aln";
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";
​
## Get lines of input files
my @aln  = <INFILE1>;
my @text = <INFILE2>;
my $num_sp = scalar(@text);
​
## Print PHYLIP header 
print OUT $num_sp."   ".$ARGV[2]."\n\n";
​
## Create variables 
my $count = 0;
my $species = "";
my $lenseq = 0;
my $count_sp = 0;
​
## Loop over all the lines in text file and 
## then just extract them from aln
foreach my $line (@text){
	
	chomp($line);
	# Just in case, make sure there are no extra 
	# chars that can cause issues 
	$species = $line;
	$species =~ s/\r//;
	$species =~ s/\n//;
	
	foreach my $aln_line (@aln) {
		# Check if line has a sp name and if it 
		# is in the subset
		## NOTE: This regex is customised to the 
		## taxa IDs you have in your file. Change 
		## if needed
		if ( $aln_line =~ /^[A-Z][a-z]..*\_/ ){
			chomp( $aln_line );
			#print $aln_line."\n";
			if ( $aln_line =~ /^$species$/ ){
				$count_sp = 1;
				print "Species to subset: ".$species."\n";
				$count += 1;
			}else{
				$count_sp = 0;
			}
		}
		# If species to be printed detected, print lines 
		# until condition FALSE
		if( $count_sp == 1 ){
			chomp( $aln_line );
			print OUT $aln_line."\n";
		}
	}
	
}
			
print "\nTotal no. species parsed: ".$count."\n";
​
## Close files
close(OUT);
close(INFILE1);
close(INFILE2);