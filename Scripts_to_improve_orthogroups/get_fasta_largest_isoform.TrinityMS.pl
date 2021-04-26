#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;


die "Please specify (1)fasta file\n" unless(@ARGV==1);

my $fastafile = $ARGV[0];
my $outfile="$fastafile\.largestIsoform";

#open(my $fastahandle, "<", $fastafile)   or die "Could not open $fastafile \n"; no need to open here for bioperl
#open(my $fastafile, "<", $fastafile)   or die "Could not open $fastafile \n";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";


use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "$fastafile", '-format' => 'Fasta');
my %fastadictionary=();
my @headersplit=();
while (my $seq = $seqio->next_seq){ ## selects one sequence at a time
   	## set variables for THIS sequence
    my $id = $seq->display_id;
	my $string = $seq->seq;
	my @split=split(/\ /, $id);
	my @spl2=split("\_", $split[0]);
	my $waste=pop @spl2;
	my $new_id=join ("\_", @spl2);
	my $len=length($string);
	#print "length = $len\n";
	if ($fastadictionary{$new_id}){
		my $len_old=length($fastadictionary{$new_id});
		if ($len >= $len_old){
			$fastadictionary{$new_id}=$string;
		}
	}
	else{
		$fastadictionary{$new_id}=$string;
	}
	
}

print "Now print new fasta with one main protein per gene.\n";

foreach my $key ( sort keys %fastadictionary){
	if ($fastadictionary{$key} eq "Sequenceunavailable"){

	}
	else{
		print $outhandle ">$key\n$fastadictionary{$key}\n";
	}
	
}

print "Finished:  input lines, output lines\n";

=cut
my $end=$d/8;
#print "$split[0]\n";
   	if ($Dict{$split[0]}){
   		my @chars = split //, $string;
   		my $number=scalar(@chars);
   		my @new_string_array;
   		for (my $i=0; $i<300; $i++){
			if ($chars[$i]){
				push (@new_string_array, "$chars[$i]");
			}
			else{
				push (@new_string_array, "N");
			}
		}
		my $new_string=join ("", @new_string_array);
		print $outhandle ">$id\n$new_string\n";
		$d++;
   	}

