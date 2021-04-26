#!/usr/bin/perl
#merge_bams_in_tophat_folder.pl
use warnings;
use strict;

die "Please give 1. Orthogroups.csv  2. Number of zeros you will allow\n" unless(@ARGV==2);

my $orthogroupscsv=$ARGV[0];
open(my $IN_b, "<", $orthogroupscsv)   or die "Could not open $orthogroupscsv \n";
my $number_to_merge=$ARGV[1];

my $outfile="$orthogroupscsv\.Z$number_to_merge\.csv";
open(my $out, ">", $outfile)   or die "Could not open $outfile \n";

print "Starting\n";
my $n=0;


my $header=<$IN_b>;
my @split_head=split("\t", $header);
my $waste=shift @split_head;

foreach my $names (@split_head){
	print $out "\t$names";
}
#print $out "\n";


print "Starting step 2\n";

while (my $line1= <$IN_b>){
	chomp $line1;
	my $n=0;
	my @a=split("\t", $line1);
	my $OG=shift(@a);
	my $len_ele=scalar (@a);
	print "length of line (columns) $len_ele\n";

	my @newline=();
	my $num_zero=0;
	for (my $i=0; $i<=8; $i++){

		my $Trin =$a[$i];
		my $species_here=$split_head[$n];
		$n++;
		if ($Trin){
			#It already has an entry.
			push (@newline, $Trin);
			#print $out "$Trin\t";
		}
		else{
			#OK, so this entry is empty, so we can fill in the gap
			push (@newline, "TRINITY_UK00000_c0_g1");
			#print $out "\t";
			$num_zero++;
		}
		print "$n\n";
	}

	print "Num zero : $num_zero\n$line1\n";

	if ($num_zero <= $number_to_merge){
		my $join=join("\t",@newline);
		print $out "$OG\t";
		print $out "$join\n";
	}
	else{
		print $out "$line1\n";
	}

	
}
