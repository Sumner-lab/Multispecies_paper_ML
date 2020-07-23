#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long qw(GetOptions);
# Super script Mulispecies machine learning

my $outfile="All_summarised_GOLD\.tsv";
open(my $out, ">", $outfile)   or die "Could not open $outfile \n";

my $version;
my $background;
my $foreground;
my $helpFlag;
my $expr_data="Experimental_data_merged";
my $run=0;
my $cpm=1;
my $INPUT_orthfinder_file="DATA/Orthofinder/Orthofinder_24.5.2020/Orthogroups.tsv";
my $FILT_here;

GetOptions(	  "help=s" => \$helpFlag,
			  'j=s'    => \$version,
			  "f=s"    => \$foreground,
			  "b=s"    => \$background,
			  "e=s"    => \$expr_data,
			  "CPM=s"    => \$cpm,
			  "r"    => \$run,
			  "orth=s"     => \$INPUT_orthfinder_file,
			  "filter=s"	=> \$FILT_here,
);

if($helpFlag){
  die "
Usage: Master.ML.pl -j <Version folder name> -f <Foreground species list (comma sep)>  -b <Background species list (comma sep)> -filter <matrix_run7.filter.scale>

options:    -r     Run the scripts through R (default = OFF).
            -e     Choose expression data folder (default= \"Experimental_data_merged\").
            -orth  Choose Orthofinder file to use in this analysis (default= \"DATA/Orthofinder/Orthofinder_24.5.2020/Orthogroups.tsv\"). 
            -CPM   Choose the cut off of minimal expression to be considered. Default =1.  

Must be the top level folder (e.g. Multispecies_FILES), which should have folders DATA and scripts and FIGURES.

*** Questions: 
    Chris Wyatt (chris.wyatt\@crg.eu)
"
}

if (defined $version && defined $foreground && defined $background){
	print "Data good, starting\n";
}
else{
	die "You need to specify version, foreground and background to run correctly.

Usage: Master.ML.pl -j <Version folder name> -f <Foreground species list (comma sep)>  -b <Background species list (comma sep)> -filter <matrix_run7.filter.scale>

options:    -r   Run the scripts through R (default = OFF)
            -e   Choose expression data folder (default= \"DATA\/Experimental_data\")   
            -orth  Choose Orthofinder file to use in this analysis (default= \"DATA/Orthofinder/Orthofinder_24.5.2020/Orthogroups.tsv\").  
            -CPM   Choose the cut off of minimal expression to be considered. Default =1. 

Must be the top level folder (e.g. Multispecies_FILES), which should have folders DATA and scripts and FIGURES.

*** Questions: 
    Chris Wyatt (chris.wyatt\@crg.eu)
";
}


my $script_whole=" Master.ML.pl -j $version -f $foreground -b $background -filter $FILT_here -e $expr_data -CPM $cpm -orth $INPUT_orthfinder_file\n";
`echo $script_whole > Script_used.txt`;

#Parse the names of the species ebing run, to paste them correctly into the R scripts.
chomp $foreground;
$foreground =~ s/[\r\n]+//g;
chomp $background;
$background =~ s/[\r\n]+//g;
my @FG=split(",", $foreground);
my @BG=split(",", $background);
my $FG_joinR=join("\"\,\"", @FG);
my $BG_joinR=join("\"\,\"", @BG);
my $ALL_joinR="\"$FG_joinR\",\"$BG_joinR\"";
my @combined=(@FG,@BG);
my $tbd_p=join(" ", @combined);
print "$tbd_p\n";
my $len_test=scalar(@FG);

`mkdir -p DATA/Orthofinder/$version`;

#Make Orthofinder files, based on the input given.

my $orthogroupscsv="$INPUT_orthfinder_file";
open(my $IN_b, "<", $orthogroupscsv)   or die "Could not open $orthogroupscsv \n";
my $OUTorthogroupscsv="DATA/Orthofinder/$version\/Orthogroups.csv";
open(my $OUT_b, ">", $OUTorthogroupscsv)   or die "Could not open $OUTorthogroupscsv \n";
my $OUTorthogroupsCOUNT="DATA/Orthofinder/$version\/Orthogroups.GeneCount.csv";
open(my $OUT_COUNT, ">", $OUTorthogroupsCOUNT)   or die "Could not open $OUTorthogroupsCOUNT \n";
my $OUTorthogroupsSINGLE="DATA/Orthofinder/$version\/SingleCopyOrthogroups.txt";
open(my $OUT_SINGLE, ">", $OUTorthogroupsSINGLE)   or die "Could not open $OUTorthogroupsSINGLE \n";

#`tr '\\r' '\\n' < DATA/Orthofinder/Orthofinder_GOLD/Orthogroups.csv > DATA/Orthofinder/Orthofinder_GOLD/Orthogroups.csv`

my $head=<$IN_b>;
$head =~ s/[\r\n]+//g;

chomp $head;
chomp $head;
my @headsp=split("\t", $head);
#print "HERE @headsp DONE \n";
#$headsp[-1] =~ s/[\r\n]+//g;

#print "->$headsp[-1]<-\n";

my @Order;
my $n1=1;
foreach my $res (@combined){
	#print "$n1 = \"$res\"\n";
	my $n=0;
	foreach my $head_file (@headsp){
		#print "M = \"$head_file\"\n";
		if ($res eq $head_file){
			push (@Order, $n);
			print "IN : $head_file\n";
		}
		$n++;
	}
	$n1++;
}

my @shorts;
my %done;
my $rep=1;
foreach my $cols (@Order){
	print $OUT_b "\t$headsp[$cols]";
	print $OUT_COUNT "\t$headsp[$cols]";

	my @name_s=split("\_", $headsp[$cols]);
	my $genus=substr ($name_s[0],0,1);
	my $speci=substr ($name_s[1],0,3);
	my $short="$genus\_$speci";
	if ($done{$short}){
		my $new="$short\_$rep";
		$rep++;
		push (@shorts, $new);
		$done{$new}="DONE";
	}
	else{
		$done{$short}="DONE";
		push (@shorts, $short);
	}
}
print $OUT_b "\n";
print $OUT_COUNT "\n";


my $shortie=join("\"\,\"", @shorts);
my @short_test;
for (my $i=0; $i <= $len_test-1; $i++){
	push (@short_test, $shorts[$i]);
}
my $joint_short_test=join("|", @short_test);
print "All  species: @shorts\n";
print "Test species: @short_test\n";



while (my $line=<$IN_b>){
	chomp $line;
	$line =~ s/[\r\n]+//g;
	my @l_sp=split("\t", $line);
	print $OUT_b "$l_sp[0]";
	print $OUT_COUNT "$l_sp[0]";
	my $switch=1;
	foreach my $cols (@Order){
		if ($l_sp[$cols]){
			print $OUT_b "\t$l_sp[$cols]";
		}
		else{
			print $OUT_b "\t";
		}
		if ($l_sp[$cols]){
			my @count=split("\, ", $l_sp[$cols]);
			my $len=scalar(@count);
			print $OUT_COUNT "\t$len";
			if ($len != 1){
				$switch=0;
			}
		}
		else{
			print $OUT_COUNT "\t0";
			$switch=0;
		}
		
		
	}
	#This is the kind of line we could add the single isoform in cases of duplicates.
	if ($switch){
		print $OUT_SINGLE "$l_sp[0]\n";
	}
	print $OUT_b "\n";
	print $OUT_COUNT "\n";
}

my $Singles= `wc -l $OUTorthogroupsSINGLE`;
print "Number of single orthogroups = $Singles\n";

`mkdir -p scripts/$version`;
`cp scripts/template_scripts/* scripts/$version`;

#Fix Script 1: getGeneSets
print "Making getGene script\n";
`perl -pi.back -e 's/VERSION_RUN/$version\/g;' scripts/$version\/getGeneSetsOrthofinder.noNameReplace.R`;
`perl -pi.back -e 's/HEREspecies.all/species.all <- c\($ALL_joinR\)/g;' scripts/$version\/getGeneSetsOrthofinder.noNameReplace.R`;
`perl -pi.back -e 's/HEREspecies.train/species.train <- c\(\"$BG_joinR\"\)/g;' scripts/$version\/getGeneSetsOrthofinder.noNameReplace.R`;
`perl -pi.back -e 's/HEREspecies.test/species.test <- c\(\"$FG_joinR\"\)/g;' scripts/$version\/getGeneSetsOrthofinder.noNameReplace.R`;


#Fix script 2: 
print "Making getExpression script\n";
`perl -pi.back -e 's/VERSION_RUN/$version\/g;' scripts/$version\/getExpressionAllOrthofinder.noNameReplace.R`;
`perl -pi.back -e 's/EATMYSHORTS/<- c(\"$shortie\")/g;' scripts/$version\/getExpressionAllOrthofinder.noNameReplace.R`;
`perl -pi.back -e 's/FIX_FOLDER/$expr_data\/g;' scripts/$version\/getExpressionAllOrthofinder.noNameReplace.R`;
`perl -pi.back -e 's/CPM_filter/$cpm\/g;' scripts/$version\/getExpressionAllOrthofinder.noNameReplace.R`;

#Fix script 3:
print "Making getResults script\n";
`perl -pi.back -e 's/VERSION_RUN/$version\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;

my @list_counts;
my @list_counts2;
my @list_counts3;
my @list_counts4;
my @list_counts5;
my @list_counts6;
my @list_counts7;
my @list_counts8;

my $nr=1;
foreach my $species (@shorts){
	push (@list_counts, "species.data.counts.log2.quantile.species.scaled.orth[[$nr\]]");
	push (@list_counts2, "species.data.tpm.orth[[$nr\]]");
	push (@list_counts3, "species.data.tpm.log2.quantile.species.scaled.orth[[$nr\]]");
	push (@list_counts4, "species.data.counts.orth.tpm.log2.quantile.species.scaled[[$nr\]]");
	push (@list_counts5, "species.data.counts.orth.tpm[[$nr\]]");
	push (@list_counts6, "species.data.counts.orth.tpm.log2.species.scaled[[$nr\]]");
	push (@list_counts7, "species.data.counts.orth.tpm.log2[[$nr\]]");
	push (@list_counts8, "species.data.counts.orth.log2.quantile.species.scaled[[$nr\]]");
	$nr++;
}
my $join_count=join(",", @list_counts);
my $join_count2=join(",", @list_counts2);
my $join_count3=join(",", @list_counts3);
my $join_count4=join(",", @list_counts4);
my $join_count5=join(",", @list_counts5);
my $join_count6=join(",", @list_counts6);
my $join_count7=join(",", @list_counts7);
my $join_count8=join(",", @list_counts8);


#fix all classi figure:
#my $num_samp_train=scalar(@BG);
my $n_classi=3;
my $n_colur=2;
my @new_lines;
my @new_numbers;
foreach my $num_samples (@combined){
	my $lineheres="lines(rev(predictionAll[1:100\,$n_classi\]), col=colours[$n_colur\])";
	push (@new_lines, "$lineheres");
	push (@new_numbers, "$n_classi");
	$n_classi++;
	$n_classi++;
	$n_colur++;
}
my $joint_pred_lines=join("\n", @new_lines);
my $joint_number_lines=join("\,", @new_numbers);
`perl -pi.back -e 's/CLASSI_LEAVE_ONE_OUT/$joint_pred_lines\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/LEGEND_LINE_FIX/$joint_number_lines\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;


`perl -pi.back -e 's/CREATECOUNTHERE/matrix.counts.log2.quantile.species.scaled.orth <- cbind($join_count\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/CREATECOUNTAGAIN2/matrix.data.tpm.orth <- cbind($join_count2\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/CREATETPMHERE/matrix.tpm.log2.quantile.species.scaled.orth <- cbind($join_count3\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/CREATEorthTPMHERE/matrix.data.counts.orth.tpm.log2.quantile.species.scaled <- cbind($join_count4\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/CREATE_REAL_TPM_HERE/matrix_real.tpm <- cbind($join_count5\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/CREATE_UNQUANTILED/matrix_real.UNQ <- cbind($join_count6\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/CREATE_LOG2_TPM_UNQUANT/matrix_real.TPM_UNQUANT <- cbind($join_count7\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;
`perl -pi.back -e 's/CREA_once_more_HERE/matrix_run7 <- cbind($join_count8\)\/g;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;

`perl -pi.back -e 's/GREPFITHERE/grep(\"$joint_short_test\"/gi;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;

`perl -pi.back -e 's/CHR_HERE_FILT/$FILT_here\/gi;' scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 `;

if ($run){
	print "Running R script 1\n";
	`R --vanilla < scripts/$version\/getGeneSetsOrthofinder.noNameReplace.R > output.ofthis.test_1`;
	print "Running R script 2\n";
	`R --vanilla < scripts/$version\/getExpressionAllOrthofinder.noNameReplace.R > output.ofthis.test_2`;
	print "Running R script 3\n";
	`R --vanilla < scripts/$version\/getResultsOrthofinder_scaled.R_RUN7 > output.ofthis.test_3`;
	print "\nResults are in FIGURES/Figure_of_Classifiers\n";
}
else{
	print "Printed scripts, but did not run them, use -r flag to do this\n";
}
print "Results are in FIGURES/Figure_of_Classifiers/$version\n";
