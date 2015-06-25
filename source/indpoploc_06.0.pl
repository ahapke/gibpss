#!/usr/bin/perl -w
#indpoploc version 06.0 Copyright 2015 Andreas Hapke
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

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "indpoploc version 06.0 Copyright 2015 Andreas Hapke\n",
"This program is part of GIbPSs Copyright 2015 Andreas Hapke\n",
"GIbPSs comes with ABSOLUTELY NO WARRANTY. GIbPSs is free software,\n",
"and you are welcome to redistribute it under certain conditions.\n",
"For details see the GNU General Public License version 3.\n",
"You should have obtained a copy of the GNU General Public License\n",
"along with GIbPSs. If not, see <http://www.gnu.org/licenses/>.\n\n";

#######################
#declare and initialize
#######################

#{
#settings:
my $individuals_fname = 'individuals.txt';#name of file with list of all individual IDs in database
my $allelefilesuffix = '_alleles.txt';#suffix of allele-files out of indloc
my $popall_infile = 'popall.txt';
my $poplocID_col = 0;
my $popall_ID_col = 1;
my $popall_seq_col = 4;
my $popallvar_col = 3;
my $indlocID_col = 0;
my $Loc_catID_col = 1;
my $valid = 'valid';
my $indall_ID_col = 3;
my $indall_seq_col = 4;
#other variables
my $infile = '';
my $ind = '';
my %ind_infile = ();# {ind}=infilename
my $outfile = '';
my $all_outsuffix = '_indpopall.txt';#suffix for allele-outfile for each ind
my $loc_outsuffix = '_indpoploc.txt';#suffix for locus-outfile for each ind
my %ind_all_outfile = ();# {ind}=allele-outfilename
my %ind_loc_outfile = ();# {ind}=locus-outfilename
my %popallcat = ();# {popall_seq} = @ [0]: poplocID [1]: popall_ID [2]: popallvar
my %indpoploc = ();# {indlocID}{poplocID}=n_indall (number of alleles of indloc just as a counter)
my $popall_seq = '';
my $poplocID = 0;
my $popall_ID = 0;
my $popallvar = '';
my $indlocID = 0;
my @popalldata = ();#contains [0]: poplocID [1]: popall_ID [2]: popallvar for an allele
my %split_loci = ();# {poplocID}{ind} = indlocID
my $n_split_loci = 0;#total number of split loci
my @temparr1 = ();
my $tempstring1 = '';
#}

########################################
#Determine names of infiles and outfiles
########################################

#{
#get individual IDs and determine names of *_alleles.txt infilenames
unless(open(INFILE,$individuals_fname)) {
	print "Cannot open $individuals_fname. Exiting..\n";
	exit;
}
@temparr1 = <INFILE>;
close INFILE;
for $ind (@temparr1) {
	chomp $ind;
	$infile = $ind.$allelefilesuffix;#build infilename
	$ind_infile{$ind} = $infile;#store in %ind_infile
}
#create outfilenames
for $ind (keys %ind_infile) {
	$ind_all_outfile{$ind} = $ind.$all_outsuffix;
	$ind_loc_outfile{$ind} = $ind.$loc_outsuffix;
}
#}

###################
#read in popall.txt
###################

#{
unless(open(INFILE,$popall_infile)){
	print "Can't open $popall_infile. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip header line
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	$popall_seq = $temparr1[$popall_seq_col];
	$poplocID = $temparr1[$poplocID_col];
	$popall_ID = $temparr1[$popall_ID_col];
	$popallvar = $temparr1[$popallvar_col];
	#add data to %popallcat
	@{$popallcat{$popall_seq}} = ($poplocID,$popall_ID,$popallvar);
}
close INFILE;
#}

#####################################################
#For each individual: Read in allele-infile,
#Produce outfiles *_indpopall.txt and *_indpoploc.txt
#Identify split loci
#####################################################

#{
for $ind (keys %ind_infile) {#loop through individuals
	#open allele-infile for this ind
	$infile = $ind_infile{$ind};
	unless(open(INFILE,$infile)){
		print "Can't open $infile. Exiting..\n";
		exit;
	}
	#open allele-outfile for this ind
	$outfile = $ind_all_outfile{$ind};
	unless(open(OUTFILE, ">$outfile")) {
		print "Can't open $outfile. Exiting..\n";
		exit;
	}
	#print headerline to outfile
	print OUTFILE "indlocID\tindall_ID\tpoplocID\tpopall_ID\tpopallvar\n";
	
	#read in infile
	$tempstring1 = <INFILE>;#skip header line
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		if ($temparr1[$Loc_catID_col] eq $valid) {#if this is a valid locus
			#if individual allele is in catalog
			if(defined $popallcat{$temparr1[$indall_seq_col]}) {
				#lookup current allele in %popallcat
				@popalldata = @{$popallcat{$temparr1[$indall_seq_col]}};
				#print data to allele-outfile
				$tempstring1 = join("\t",@popalldata);
				print OUTFILE "$temparr1[$indlocID_col]\t$temparr1[$indall_ID_col]\t$tempstring1\n";
				#add locus to %indpoploc (even if you already did so...)
				++$indpoploc{$temparr1[$indlocID_col]}{$popalldata[0]};
			} else {#if not in catalog: warn				
				print "$ind locus $temparr1[$indlocID_col] allele $temparr1[$indall_ID_col] is unknown.\n";
			}			
		}
	}
	close INFILE;
	close OUTFILE;
	#open locus-outfile
	$outfile = $ind_loc_outfile{$ind};
	unless(open(OUTFILE, ">$outfile")) {
		print "Can't open $outfile. Exiting..\n";
		exit;
	}
	#print headerline
	print OUTFILE "indlocID\tpoplocID\n";
	#loop through %indpoploc and print data to locus-outfile
	for $indlocID (sort {$a <=> $b} keys %indpoploc) {
		for $poplocID (sort {$a <=> $b} keys %{$indpoploc{$indlocID}}) {
			print OUTFILE "$indlocID\t$poplocID\n";
		}
	}
	close OUTFILE;
	#search for split loci in this individual
	for $indlocID (keys %indpoploc) {
		#if this indlocID is linked to more than one poplocID
		if (keys %{$indpoploc{$indlocID}} > 1) {
			#loop through these poplocIDs
			for $poplocID (keys %{$indpoploc{$indlocID}}) {
				$split_loci{$poplocID}{$ind} = $indlocID;#store in %split_loci
			}
		}
	}	
	%indpoploc = ();#set back and next ind
}
#}

###############################
#Produce outfile split_loci.txt
###############################

#{
$n_split_loci = keys %split_loci;#number of split loci
print "Found $n_split_loci split loci. Printing to file split_loci.txt..\n";
#open split_loci.txt and print data
unless(open(OUTFILE, ">split_loci.txt")) {
	print "Can't open split_loci.txt. Exiting..\n";
	exit;
}
#print headerline
print OUTFILE "poplocID\tind\tindlocID\n";
#print data
for $poplocID (sort {$a <=> $b} keys %split_loci) {
	for $ind (sort keys %{$split_loci{$poplocID}}) {
		$indlocID = $split_loci{$poplocID}{$ind};
		print OUTFILE "$poplocID\t$ind\t$indlocID\n";
	}
}
close OUTFILE;
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "run took $run_s seconds.\n";

exit;
