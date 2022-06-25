#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

die "Need to give an input file name: DEGs file (1) and a number of permutations (2)\n" unless (@ARGV==2);

my $info= $ARGV[0];
my $user_option = $ARGV[1];
open(my $input, "<", $info) or die "Could not open $info\n";

my $outfile="Permutation_data.txt";
open(my $out, '>>', $outfile)   or die "Could not open $outfile \n";

my $head=<$input>;
my @head_a= split("\t", $head);
#shift @head_a;

#print "$head_a[0]\n";

my %species_genes;
my %species_DEGs;

while (my $line = <$input>){
	chomp $line;
	my @l = split("\t", $line);
	my $OG= $l[0];
	for (my $i=1; $i<=9; $i++){
		my $t = $i+9;

		#Test if a gene is present
		if ($l[$t]=~ m/TRINITY/){
			$species_genes{"$head_a[$i]"}{"$OG"}="$OG";
		}

		#Test if the gene is a DEG
		if ( $l[$i]=~ m/UP|DW/ ){
			$species_DEGs{$head_a[$i]}++;
		}

	}
}


print $out "\t1\t2\t3\t4\t5\t6\t7\t8\t9\n";

for (my $p=1; $p<=$user_option; $p++){
	print "In loop $p\n";
	print $out "$p";

	my @all;

	for (my $i=1; $i<=9; $i++){
		my $species= $head_a[$i];
		my $sp_deg= $species_DEGs{$species};
		print "$species\t$sp_deg\n";

		my @values;

		#print Dumper %species_genes;

		foreach my $key ( keys %species_genes ) { 
			foreach my $OGS ( keys %{$species_genes{$key}} ){
				push(@values,$species_genes{$key}{$OGS});
			}
		   
		}

		my @a = @values;

		my @items;
		for ( 1 .. $sp_deg )
		{
		  push @items, splice @a, rand @a, 1;
		}

		my $joint=join("\n", @items);
		#print $out "$species\n$joint\n";

		push (@all,@items);

	}

	my %counts = ();
	for (@all) {
	   $counts{$_}++;
	}

	my @array_h;
	foreach my $keys (keys %counts) {
	   #print "$keys = $counts{$keys}\n";
	   push (@array_h, $counts{$keys});
	}

	my %counts_2 = ();
	for (@array_h) {
		#print "$_\n";
		$counts_2{$_}++;
	}


	for (my $i=1; $i<=9; $i++){
		if ($counts_2{$i}){
			print $out "\t$counts_2{$i}";
		}
		else{
			print $out "\t0";
		}
		
	}

	print $out "\n";
}
