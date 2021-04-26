#!/usr/bin/perl
#merge_bams_in_tophat_folder.pl
use warnings;
use strict;

die "Please give 1. Orthogroups.csv  2. Number of groups you wish to merge\n" unless(@ARGV==2);

my $orthogroupscsv=$ARGV[0];
open(my $IN_b, "<", $orthogroupscsv)   or die "Could not open $orthogroupscsv \n";
my $number_to_merge=$ARGV[1];

my $outfile="$orthogroupscsv\.filt$number_to_merge\.csv";
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

my %RSEM;
foreach my $species (@split_head){
	chomp $species;
	my $input_Q_rsem="../DATA\/Experimental_data_merged\/$species\/Queen\/RSEM.isoforms.results";
	my $input_W_rsem="../DATA\/Experimental_data_merged\/$species\/Worker\/RSEM.isoforms.results";
	open(my $IN_rsem, "<", $input_Q_rsem)   or die "Could not open $input_Q_rsem \n";
	open(my $IN_rsem2, "<", $input_W_rsem)   or die "Could not open $input_W_rsem \n";

	my $header_rsem1=<$IN_rsem>;
	while (my $line1= <$IN_rsem>){
		chomp $line1;
		my @sp_rsem=split("\t", $line1);
		if ($RSEM{$species}{$sp_rsem[1]}){
			my $old=$RSEM{$species}{$sp_rsem[1]};
			if ($sp_rsem[6] >= $old){
				$RSEM{$species}{$sp_rsem[1]}=$sp_rsem[6];
			}
		}
		else{
			$RSEM{$species}{$sp_rsem[1]}=$sp_rsem[6];
		}
	}
	close $IN_rsem;

	

	my $header_rsem2=<$IN_rsem2>;

	while (my $line1= <$IN_rsem2>){
		chomp $line1;
		my @sp_rsem=split("\t", $line1);
		if ($RSEM{$species}{$sp_rsem[1]}){
			my $old=$RSEM{$species}{$sp_rsem[1]};
			if ($sp_rsem[6] >= $old){
				$RSEM{$species}{$sp_rsem[1]}=$sp_rsem[6];
			}
		}
		else{
			$RSEM{$species}{$sp_rsem[1]}=$sp_rsem[6];
		}
	}
	close $IN_rsem2;

	print "Species processed : $species \n";
}



print "Starting step 2\n";

while (my $line1= <$IN_b>){
	chomp $line1;
	my $n=0;
	my @a=split("\t", $line1);
	my $OG=shift(@a);
	print $out "$OG";
	#print "$line1\n";
	foreach my $Trin (@a){
		my $species_here=$split_head[$n];
		$n++;
		if ($Trin){
			my @split=split("\, ", $Trin);
			my $count=scalar(@split);
			if ($count <= $number_to_merge){
				my $HIGH;
				my $HIGH_value;
				#print "$Trin $n\n";
				foreach my $isoforms (@split){
					my $highest_expr=$RSEM{$species_here}{$isoforms};
					if ($HIGH){
						if ($highest_expr >= $HIGH_value){
							$HIGH=$isoforms;
							$HIGH_value=$highest_expr;
							#print "HEREEEE2 $species_here $HIGH $HIGH_value $highest_expr\n";
						}
					}
					else{
						$HIGH=$isoforms;
						$HIGH_value=$highest_expr;
						#print "HEREEEE $species_here $HIGH $HIGH_value $highest_expr\n";
					}
				}
				print $out "\t$HIGH";
				#print "this is the : $HIGH\n";
			}
			else {
				if ($Trin){
					print $out "\t$Trin";
				}
				else{
					print $out "\t";
				}
				
			}
		}
		else{
			print $out "\t";
		}
		
	}
	print $out "\n";
}
