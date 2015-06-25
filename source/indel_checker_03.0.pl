#!/usr/bin/perl -w
#indel_checker version 03.0 Copyright 2015 Andreas Hapke
#This file is part of GIbPSs Copyright 2015 Andreas Hapke
#
#GIbPSs is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#GIbPSs is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#Contact:
#Andreas Hapke 
#Johannes Gutenberg University Mainz
#Institute of Anthropology
#Anselm-Franz-von-Benzel-Weg 7
#D-55099 Mainz, Germany
#E-Mail: ahapke2@gmail.com

use strict;
use warnings;

#Copyright message to screen
print "indel_checker version 03.0 Copyright 2015 Andreas Hapke\n",
"This program is part of GIbPSs Copyright 2015 Andreas Hapke\n",
"GIbPSs comes with ABSOLUTELY NO WARRANTY. GIbPSs is free software,\n",
"and you are welcome to redistribute it under certain conditions.\n",
"For details see the GNU General Public License version 3.\n",
"You should have obtained a copy of the GNU General Public License\n",
"along with GIbPSs. If not, see <http://www.gnu.org/licenses/>.\n\n";
#Message about usearch to screen
print "This program uses USEARCH. Edgar, RC (2010) Search and clustering\n",
"orders of magnitudes faster than BLAST, Bioinformatics 26(19),2460-2461.\n\n";

######################
#variables for usearch
######################

#{
my %user_settings = (
UP => '',#full path to usearch
DP => '',#full path to main database directory
id => 'n',#identity threshold
minsl => 'n',#min of shorter seqlength / longer seqlength
query_cov => 'n',#min fraction of query sequence aligned
target_cov => 'n'#min fraction of target (centroid) sequence aligned
);
my %defaults = (
UP => '',#full path to usearch
DP => '',#full path to main database directory
id => '0.9',#identity threshold
minsl => '0.9',#min of shorter seqlength / longer seqlength
query_cov => '0.9',#min fraction of query sequence aligned
target_cov => '0.9'#min fraction of target (centroid) sequence aligned
);
my $n_arg = 0;#number of arguments provided by user
my $flag = '';#command flag
my $val = '';#value
my $inpath = '';#full path to infile for usearch
my $outdirpath = '';#full path to outfile directory
my $outpath = '';#full path to outfile from usearch
my $ucall = '';#full text for call of usearch
#}

################
#other variables
################

#{
my %sel_loc = ();# {poplocID}=1 selected loci (by data_selector)
my %poploc_sl = ();# {poplocID} = sl (sequence length)
my %popall = ();# {sl}{poplocID}{popall_ID}=popall_seq
my $poplocID = 0;
my $popall_ID = 0;
my $sl = 0;#sequence length
my $min_sl = 32;#minimum sequence length accepted by usearch
my $seq = '';#allele sequence
my %clust_poploc = ();# {clusterNo}{poplocID} = 1 cluster number out of usearch
my %clust_size = ();# {clusterNo} = size (number of poplocs in a cluster)
my $cluster = 0;#cluster number
my $n_clusters = 0;#number of clusters identified by usearch
my $n_loc = 0;#number of loci
my $n_loc_no_indel = 0;#number of loci, which are not indel variants of another locus
my $n_loc_indel = 0;#number of loci, which are indel variants of another locus
my $n_clust_indel_loc = 0;#number of clusters that contain several loci
my $perc_clusters = 0;#percent of clusters that contain indel loci
my $tempstring1 = '';
my @temparr1 = ();
my $i = 0;#counter
#}

########################
#Take over user settings
########################

#{
if (@ARGV) {#if user provided arguments
	$n_arg = @ARGV;#determine number of arguments
	for ($i = 0; $i < ($n_arg - 1); $i += 2) {
		$flag = $ARGV[$i];
		$flag =~ s/^-//;
		$val = $ARGV[$i+1];
		#if flag is defined, take over
		if (defined $user_settings{$flag}) {
			$user_settings{$flag} = $val;
		}
	}
} else {
	print "Please provide at least two arguments:\n",
	"full pathes to usearch and to main database directory, example:\n",
	"-UP C:\\myprograms\\usearch.exe -DP C:\\data\\db1\n",
	"Exiting..\n";
	exit;	
}
#}

####################################################################
#Test user settings for plausibility
#(limited test) and replace missing/implausible values with defaults
####################################################################

#{
#full path to usearch -UP
unless (length $user_settings{'UP'} > 0 and -e $user_settings{'UP'}) {
	print "I don't find usearch (-UP). Exiting..\n";
	exit;
}
#full path to main database directory
unless (length $user_settings{'DP'} > 0 and -d $user_settings{'DP'}) {
	print "I don't find main database directory (-DP). Exiting..\n";
	exit;
}
#identity threshold -id decimal between 0 and 1
unless((($user_settings{'id'} =~ /^\d$/) or ($user_settings{'id'} =~ /^\d\.\d+$/))
and ($user_settings{'id'} >= 0) and ($user_settings{'id'} <= 1)) {
	$user_settings{'id'} = $defaults{'id'};
}
# -minsl decimal between 0 and 1
unless((($user_settings{'minsl'} =~ /^\d$/) or ($user_settings{'minsl'} =~ /^\d\.\d+$/))
and ($user_settings{'minsl'} >= 0) and ($user_settings{'minsl'} <= 1)) {
	$user_settings{'minsl'} = $defaults{'minsl'};
}
# -query_cov decimal between 0 and 1
unless((($user_settings{'query_cov'} =~ /^\d$/) or ($user_settings{'query_cov'} =~ /^\d\.\d+$/))
and ($user_settings{'query_cov'} >= 0) and ($user_settings{'query_cov'} <= 1)) {
	$user_settings{'query_cov'} = $defaults{'query_cov'};
}
# -target_cov decimal between 0 and 1
unless((($user_settings{'target_cov'} =~ /^\d$/) or ($user_settings{'target_cov'} =~ /^\d\.\d+$/))
and ($user_settings{'target_cov'} >= 0) and ($user_settings{'target_cov'} <= 1)) {
	$user_settings{'target_cov'} = $defaults{'target_cov'};
}
#}

########################################
#Read in selected loci and sequence data
########################################

#{
if (-f "export/sel_loc.txt") {#if a locus selection file exists
	print "Found a locus selection file..\n";
	#read in selected loci
	unless(open(INFILE, "export/sel_loc.txt")) {
		print "Cannot open export/sel_loc.txt. Exiting..\n";
		exit;
	}
	print "Reading in your locus selection..\n";
	while ($poplocID = <INFILE>) {#read in poplocIDs in the file
		chomp $poplocID;
		$sel_loc{$poplocID} = 1;#store in %sel_loc
	}
	close INFILE;
	if (keys %sel_loc == 0) {#if the file was empty
		print "Found 0 loci in locus selection file, exiting..\n";
		exit;
	}
} else {#no locus selection file exists, create %sel_loc from file poploc.txt
	print "I will select all loci in the database..\n";
	#read in file poploc.txt
	unless(open(INFILE, "poploc.txt")) {
		print "Cannot open poploc.txt. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip header line
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		$sel_loc{$temparr1[0]} = 1;#store poplocID in %sel_loc
	}
	close INFILE;
}
#read in file poploc.txt
unless(open(INFILE, "poploc.txt")) {
	print "Cannot open poploc.txt. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip header line
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#if current poploc has sl >= min and if it is in selection
	if(($temparr1[1] >= $min_sl) and (defined $sel_loc{$temparr1[0]})) {
		$poploc_sl{$temparr1[0]} = $temparr1[1];#store poplocID and sl
	}
	#if current poploc has sl < min and is in selection
	elsif(($temparr1[1] < $min_sl) and (defined $sel_loc{$temparr1[0]})) {
		print "Your selection contains loci that are shorter than $min_sl.\n",
		"Please refine your selection and exclude such loci. Exiting..\n";
		exit;
	}
}
close INFILE;

#read in file popall.txt
unless(open(INFILE, "popall.txt")) {
	print "Cannot open popall.txt. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip header line
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#if current poploc is in %poploc_sl (in selection)
	if (defined  $poploc_sl{$temparr1[0]}) {
		$poplocID = $temparr1[0];
		$sl = $poploc_sl{$poplocID};
		$popall_ID = $temparr1[1];
		$seq = $temparr1[4];
		$popall{$sl}{$poplocID}{$popall_ID} = $seq;
	}
}
#}

######################################################################
#Determine pathes, create outdirectory, build text for call of usearch
######################################################################

#{
if ($user_settings{'DP'} =~ m/\//) {#if path contains slash
	$user_settings{'DP'} =~ s/\/$//;#remove trailing slash if any
	$outdirpath = $user_settings{'DP'} . '/indelcheck';
	$inpath = $outdirpath . '/popall.fasta';
	$outpath = $outdirpath . '/u_out.txt';
}
elsif ($user_settings{'DP'} =~ m/\\/) {#if path contains backslash
	$user_settings{'DP'} =~ s/\\$//;#remove trailing backslash if any
	$outdirpath = $user_settings{'DP'} . '\indelcheck';
	$inpath = $outdirpath . '\popall.fasta';
	$outpath = $outdirpath . '\u_out.txt';
}
if (-d $outdirpath) {
	print "$outdirpath already exists, please delete or rename.\nExiting..\n";
	exit;
}
unless (mkdir "$outdirpath") {
	print "Can't create directory $outdirpath. Exiting..\n";
	exit;
}
$ucall = "$user_settings{'UP'} -cluster_smallmem $inpath -id $user_settings{'id'} " .
"-fulldp -minsl $user_settings{'minsl'} -query_cov $user_settings{'query_cov'} " .
"-target_cov $user_settings{'target_cov'} -uc $outpath";
#}

#####################################################################
#Produce fasta infile for usearch with one allele sequence per poploc
#####################################################################

#{
#open outfile
unless(open(OUTFILE, ">$inpath")) {
	print "Cannot open $inpath Exiting..\n";
	exit;
}
#print data to outfile
for $sl (reverse sort {$a <=> $b} keys %popall) {
	for $poplocID (sort {$a <=> $b} keys %{$popall{$sl}}) {
		$popall_ID = (sort {$a <=> $b} keys %{$popall{$sl}{$poplocID}})[0];
		$seq = $popall{$sl}{$poplocID}{$popall_ID};
		print OUTFILE ">","$poplocID","_","$popall_ID\n$seq\n";
	}
}
close OUTFILE;
#}

#############
#Call usearch
#############

system ("$ucall");

#############################
#Read in outfile from usearch
#############################

#{
unless (open(INFILE,"$outpath")) {
	print "Cannot open $outpath. Exiting..\n";
	exit;
}
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	if ($temparr1[0] eq 'H') {#if this is a hit
		$poplocID = $temparr1[8];
		$poplocID =~ s/_\d+$//;#get poplocID out of fasta header
		$clust_poploc{$temparr1[1]}{$poplocID}=1;#store in %clust_poploc
	}
	elsif ($temparr1[0] eq 'C') {#if this is a centroid
		$poplocID = $temparr1[8];
		$poplocID =~ s/_\d+$//;#get poplocID out of fasta header
		$clust_poploc{$temparr1[1]}{$poplocID}=1;#store in %clust_poploc
		$clust_size{$temparr1[1]} = $temparr1[2];#store cluster size in %clust_size
	}
}
close INFILE;
#}

###################################################
#Produce outfiles clusters.txt and no_indel_loc.txt
###################################################

#{
#open outfiles
unless (open(OUTCLUST,">indelcheck/clusters.txt")) {
	print "Cannot open indelcheck/clusters.txt. Exiting..\n";
	exit;
}
unless (open(OUTNOINDEL,">indelcheck/no_indel_loc.txt")) {
	print "Cannot open indelcheck/no_indel_loc.txt. Exiting..\n";
	exit;
}
#print headerline to file OUTCLUST
print OUTCLUST "Cluster_No\tpoplocID\n";

#print data
for $cluster (sort {$a <=> $b} keys %clust_poploc) {
	++$n_clusters;#count cluster
	for $poplocID (sort {$a <=> $b} keys %{$clust_poploc{$cluster}}) {
		++$n_loc;#count poploc
		print OUTCLUST "$cluster\t$poplocID\n";
		if ($clust_size{$cluster} == 1) {#if this cluster contains only one locus
			print OUTNOINDEL "$poplocID\n";
			++$n_loc_no_indel;#count locus that is not indel variant of another
		}
	}
}
$n_clust_indel_loc = $n_clusters - $n_loc_no_indel;#determine number of clusters that contain indel loci
$n_loc_indel = $n_loc - $n_loc_no_indel;#determine number of loci that are indel variants of another locus
$perc_clusters = $n_clust_indel_loc / $n_clusters * 100;#percentage of clusters that contain several loci
close OUTCLUST;
close OUTNOINDEL;
#}

######################
#Produce a report file
######################

#{
#open report file
unless (open(REPORT, ">indelcheck/indelcheck_report.txt")) {
	print "Cannot open indelcheck/indelcheck_report.txt. Exiting..\n";
	exit;
}
print REPORT "INDELCHECKER:\n",
"Identification of loci that could be indel variants of each other.\n",
"This program used USEARCH. Edgar, RC (2010) Search and clustering\n",
"orders of magnitudes faster than BLAST, Bioinformatics 26(19),2460-2461.\n\n";

#If you can open it, read in report file out of data_selector and print contents to report
if (open(INFILE,"export/sel_rep.txt")) {
	print REPORT "Your selection of data before this analysis:\n\n";
	while ($tempstring1 = <INFILE>) {
		print REPORT "$tempstring1";
	}
	print REPORT "indel_checker used your selection of loci.\n";
close INFILE;
} else {#could not open report file out of data_selector
	print REPORT "indel_checker analyzed all loci in the database.\n";
}
print REPORT "From each locus, it exported the sequence of one allele into file\n",
"indelcheck/popall.fasta.\n\n",
"indel_checker called usearch as follows:\n\n";
print REPORT "$ucall\n";
print REPORT "USEARCH produced this outfile:\n",
"indelcheck/u_out.txt\n";
print REPORT "indel_checker analyzed the outfile and produced these files:\n",
"indel_check/clusters.txt, indelcheck/no_indel_loc.txt\n";
print REPORT "depth_analyzer can use file clusters.txt.\n";
print REPORT "You can later use file no_indel_loc.txt to select\n",
"only those loci, that are not indel variants of another locus.\n";
print REPORT "Use data_selector, select loci, list of loci in a file\n",
"and provide indelcheck/no_indel_loc.txt as file name.\n";
print REPORT "\nSummary statistics:\n",
"$n_loc loci\n",
"$n_clusters clusters\n",
"$n_loc_no_indel loci are not indel variants of other loci.\n",
"$n_loc_indel loci could be indel variants of other loci and are clustered into\n",
"$n_clust_indel_loc clusters.\n",
"$perc_clusters\% of clusters contain more than one locus.\n";
#}
