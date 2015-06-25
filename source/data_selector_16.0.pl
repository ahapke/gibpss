#!/usr/bin/perl -w
#data_selector version 16.0 Copyright 2015 Andreas Hapke
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
use v5.12;

#uses subroutines
#select individuals
#select loci
#select genotypes
#varpos

#Keyword Infilecolumns! Order of columns in an infile is important.

#Copyright message to screen
print "data_selector version 16.0 Copyright 2015 Andreas Hapke\n",
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
#settings
my $individuals_fname = 'individuals.txt';#out of indloc, list of all individual IDs in database
my $poplocfilename = 'poploc.txt';#poploc-outfile from poploc
my $popallfilename = 'popall.txt';#popall-outfile from poploc
my $indpoplocsuffix = '_indpoploc.txt';#suffix of *_indpoploc.txt outfile of indpoploc
my $indpopallsuffix = '_indpopall.txt';#uffix of *_indpopall.txt outfile of indpoploc
my $alleles_suffix = '_alleles.txt';#suffix of *_alleles outfiles of indloc
my $loci_suffix = '_loci.txt';#suffix of *_loci outfiles out of indloc
my $splitlocfilename = 'split_loci.txt';#split_loci outfile of indpoploc
my $depthfilename = 'depth_analysis/da_all_loc.txt';#out of depth_analyzer: depth data
my $pairlocgr_fname = 'pairs/loc_groups.txt';#out of pair_finder: loc_group n_loc poplocID
my $pairselfm_fname = 'pairs/self_match.txt';#out of pair_finder:
											#self-matching loci: poplocID n_seqpairs n_ind
my $pair_rep_fname = 'pairs/pairs_rep.txt';#out of pair_finder: report-file
#other variables
my %poploc = ();# {poploc_ID}=(sl,cons,nSNP,varpos,n_all,nInd,nIndloc,ties,lon_over,lon_com_over,merge_conflict)
my %popall = ();# {poploc_ID}{popall_ID}=(popall_seq,popallvar,popall_seq_notmerged)
my %sel_ind = ();# {ind}=1
my %sel_loc = ();# {poplocID}=1
my %sel_gt = ();# {poplocID}{ind}=1 contains only existing genotypes
my $previous_sel = 0;#1 if previous selections have been found
my $ind = '';#an individual ID
my $poplocID = 0;#population-locus ID
my $popall_ID = 0;#population-allele ID
my $startselection = '0';
my %startmenu = (
1 => 'select individuals',
2 => 'select loci',
3 => 'select genotypes',
4 => 'quit'
);
my $infilename = '';#an infilename
my $key1 = '';
my $val1 = '';
my $entry = '';
my @temparr1 = ();
my $tempstring1 = '';
#}

print "\nDATA SELECTOR\n\n";

#####################################################
#Read in poploc.txt and popall.txt outfiles of poploc
#####################################################

#{
print "Loading data, please wait.\n";
#read in poploc.txt outfile of poploc and populate %poploc
unless(open(INFILE,$poplocfilename)) {
	print "Cannot open $poplocfilename. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	#Infilecolumns!
	@temparr1 = split(/\t/,$tempstring1);
	$poplocID = shift @temparr1;#get poplocID and remove it from array
	@{$poploc{$poplocID}} = @temparr1;#load data into %poploc
}
close INFILE;

#read in popall.txt outfile of poploc and populate %popall
unless(open(INFILE,$popallfilename)) {
	print "Cannot open $popallfilename. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	#Infilecolumns!
	@temparr1 = split(/\t/,$tempstring1);
	$poplocID = shift @temparr1;#get poplocID and remove it from array
	$popall_ID = shift @temparr1;#get popallID and remove it from array
	@{$popall{$poplocID}{$popall_ID}} = @temparr1;#load data into %popall

}
close INFILE;

#}

############################################################################
#Check for existing preselections, read in or initialize from whole database
############################################################################

#{
print "Searching for previous selections, please wait.\n";
if (-d "export") {
	if (-f "export/sel_ind.txt") {#if there is an individual selection-file
		unless(open(INFILE,"export/sel_ind.txt")) {
			print "Cannot open export/sel_ind.txt. Exiting..\n";
			exit;
		}
		print "Found previous selection of individuals, loading, please wait.\n";
		$previous_sel = 1;
		@temparr1 = <INFILE>;
		close INFILE;
		#Infilecolumns!
		for $ind (@temparr1) {
			chomp $ind;
			$sel_ind{$ind} = 1;#store selected inds in %sel_ind
		}
	}
	if (-f "export/sel_loc.txt") {#if there is a locus selection-file
		unless(open(INFILE,"export/sel_loc.txt")) {
			print "Cannot open export/sel_loc.txt. Exiting..\n";
			exit;
		}
		print "Found previous selection of loci, loading, please wait.\n";
		$previous_sel = 1;
		#Infilecolumns!
		while ($poplocID = <INFILE>) {#read in poplocIDs
			chomp $poplocID;
			$sel_loc{$poplocID} = 1;#store in %sel_loc
		}
		close INFILE;
	}
	if (-f "export/sel_gt.txt") {#if there is a genotype selection-file
		unless(open(INFILE,"export/sel_gt.txt")) {
			print "Cannot open export/sel_gt.txt. Exiting..\n";
			exit;
		}
		print "Found previous selection of genotypes, loading, please wait.\n";
		$previous_sel = 1;
		#Infilecolumns!
		while ($tempstring1 = <INFILE>) {#read in
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			$sel_gt{$temparr1[0]}{$temparr1[1]} = 1;#store in %sel_gt
		}
		close INFILE;
	}	
} else {#if directory "export" doesnt exist, create it
	unless (mkdir 'export') {
		print "Cannot create folder export. Exiting.\n";
		exit;
	}
}
#if %sel_ind is still empty, create it from file individuals.txt (using all inds)
if (keys %sel_ind == 0) {
	unless(open(INFILE,"$individuals_fname")) {
		print "Cannot open $individuals_fname. Exiting..\n";
		exit;
	}
	@temparr1 = <INFILE>;
	close INFILE;
	#Infilecolumns!
	for $ind (@temparr1) {
		chomp $ind;
		$sel_ind{$ind} = 1;#store in %sel_ind
	}
}
#if %sel_loc is still empty, create from %poploc (using all loci)
if (keys %sel_loc == 0) {
	for $poplocID (keys %poploc) {
		$sel_loc{$poplocID} = 1;
	}
}
#if %sel_gt is still empty, create from %sel_ind
#and individual files *_indpoploc.txt (out of indpoploc)
#loop through selected individuals
if (keys %sel_gt == 0) {
	for $ind (keys %sel_ind) {
		#open file *_indpoploc.txt for this ind
		$infilename = $ind.$indpoplocsuffix;#create infilename
		unless(open(INFILE,"$infilename")) {
			print "Cannot open $infilename. Exiting..\n";
			exit;
		}
		$tempstring1 = <INFILE>;#skip headerline
		while ($tempstring1 = <INFILE>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#store poplocID for this ind in %sel_gt
			#no problem if same poplocID appears several times in *_indpoploc.txt
			#Infilecolumns!
			$sel_gt{$temparr1[1]}{$ind} = 1;
		}
		close INFILE;	
	}
}
#}

####################################################
#display start-menu until you get a valid user-entry
####################################################

#{
while ($startselection eq '0') {
	print "\nDATA SELECTION START MENU\n\n";
	if ($previous_sel == 1) {
		print "Previous selections have been loaded.\n",
		"Your selection will create a subset of previous selections.\n\n";
	}
	for $key1 (sort {$a <=> $b} keys %startmenu) {
		print "$key1 $startmenu{$key1}\n";
	}
	print "\nPlease select an option (number): ";
	$entry = <STDIN>;
	chomp $entry;
	if (defined $startmenu{$entry}) {
		$startselection = $entry;
	}
}
#}

#####################################
#call the selected subroutine or exit
#####################################

#{
if ($startselection eq '1') {
	select_individuals(\%sel_ind,\%sel_loc,\%sel_gt);
}
elsif ($startselection eq '2') {
	select_loci(\%poploc,\%popall,\%sel_ind,\%sel_loc,\%sel_gt,
	$indpopallsuffix,$splitlocfilename,$depthfilename,$pairlocgr_fname,
	$pairselfm_fname,$pair_rep_fname,$alleles_suffix,$loci_suffix);
}
elsif ($startselection eq '3') {
	select_genotypes(\%sel_ind,\%sel_loc,\%sel_gt,$indpopallsuffix,$alleles_suffix,$loci_suffix);
}
elsif ($startselection eq '4') {
	print "No data selected. Exiting..\n"; exit;
}
#}

##########################
#print new selection files
##########################

#{
#sel_ind.txt
unless(open(OUTFILE, ">export/sel_ind.txt")) {
	print "Cannot open file export/sel_ind.txt, exiting ...\n\n";
	exit;
}
for $ind (sort keys %sel_ind) {
	print OUTFILE "$ind\n";
}
close OUTFILE;
#sel_loc.txt and sel_gt
unless(open(OUTSLOC, ">export/sel_loc.txt")) {
	print "Cannot open file export/sel_loc.txt, exiting ...\n\n";
	exit;
}
unless(open(OUTSGT, ">export/sel_gt.txt")) {
	print "Cannot open file export/sel_gt.txt, exiting ...\n\n";
	exit;
}
for $poplocID (sort {$a <=> $b} keys %sel_loc) {
	print OUTSLOC "$poplocID\n";
	for $ind (sort keys %{$sel_gt{$poplocID}}) {
		print OUTSGT "$poplocID\t$ind\n";
	}
}
close OUTSLOC;
close OUTSGT;
#}

exit;

############
#Subroutines
############

#definition of subroutine select_individuals
#selects individuals and updates %sel_ind, %sel_loc, %sel_gt
#expects references to
#%sel_ind: {ind}=1
#%sel_loc: {poplocID}=1
#%sel_gt: {poplocID}{ind}=1

sub select_individuals {
	#declare and initialize: _r means a reference
	my ($sel_ind_r,$sel_loc_r,$sel_gt_r) = @_;
	my $infilename = '';#an infilename
	my %user_sel_ind = ();#individuals selected by user: {ind}=1
	my $ind = '';#individual ID
	my $poplocID = 0;#a locus ID
	my $n_ind = 0;#number of individuals
	my $n_ind_dataset = 0;#number of individuals in dataset
	my $n_ind_user = 0;#number of individuals selected by user
	my $n_loc_dataset = 0;#number of loci in dataset
	my $n_gt_dataset = 0;#number of genotypes in dataset
	my $tempstring1 = '';
	
	#Welcome:
	print "\nSELECT INDIVIDUALS:\n\n";
	#ask user for selection file
	print "I need a list of selected individuals in a text file.\n",
	"Format: one individual ID per line, no spaces, no blank lines.\n\n",
	"Please enter filename: ";
	$infilename = <STDIN>;
	chomp $infilename;
	#open file or die
	unless(open(INFILE,$infilename)) {
		print "Cannot open $infilename. Exiting..\n";
		exit;
	}
	#read in file and populate %user_sel_ind
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		$user_sel_ind{$tempstring1} = 1;
	}
	close INFILE;
	#determine starting dimensions of database or previous selection
	#n_ind
	$n_ind_dataset = keys %{$sel_ind_r};
	#n_loc
	$n_loc_dataset = keys %{$sel_loc_r};
	#n_gt
	for $poplocID (keys %{$sel_gt_r}) {
		$n_gt_dataset += keys %{$$sel_gt_r{$poplocID}};
	}
	#determine number of individuals selected by user
	$n_ind_user = keys %user_sel_ind;
	#open export/sel_rep.txt for append
	unless(open(OUTREP, ">>export/sel_rep.txt")) {
		print "Cannot open file export/sel_rep.txt, exiting ...\n\n";
		exit;
	}
	print OUTREP
	"Selection of individuals: starting from\n\n",
	"$n_ind_dataset individuals\n$n_loc_dataset loci\n$n_gt_dataset genotypes\n\n",
	"Selected data from $n_ind_user individuals in file $infilename:\n\n";	
	
	#loop through individuals in %sel_ind
	for $ind (keys %{$sel_ind_r}) {
		#If this ind is not defined in user selection
		unless (defined $user_sel_ind{$ind}) {
			delete $$sel_ind_r{$ind};#delete it from %sel_ind
		}
	}
	#loop through loci in %sel_gt
	for $poplocID (keys %{$sel_gt_r}) {
		#loop through individuals in this locus
		for $ind (keys %{$$sel_gt_r{$poplocID}}) {
			#If this ind is not defined in user selection
			unless (defined $user_sel_ind{$ind}) {
				delete $$sel_gt_r{$poplocID}{$ind};#delete it from %sel_gt
			}
		}
		#Determine number of individuals left in the locus
		$n_ind = keys %{$$sel_gt_r{$poplocID}};
		#If none are left
		if ($n_ind == 0) {
			delete $$sel_gt_r{$poplocID};#delete locus from %sel_gt
			delete $$sel_loc_r{$poplocID};#delete locus from %sel_loc
		}
	}
	#determine dimensions of selected dataset
	#n_ind
	$n_ind_dataset = keys %{$sel_ind_r};
	#n_loc
	$n_loc_dataset = keys %{$sel_loc_r};
	#n_gt
	$n_gt_dataset = 0;
	for $poplocID (keys %{$sel_gt_r}) {
		$n_gt_dataset += keys %{$$sel_gt_r{$poplocID}};
	}
	print OUTREP "$n_ind_dataset individuals\n$n_loc_dataset loci\n$n_gt_dataset genotypes\n\n\n";
	close OUTREP;
}

#definition of subroutine select_loci
#selects loci and updates %sel_ind, %sel_loc, %sel_gt
#expects
#ref to %poploc: {poploc_ID}=(sl,cons,nSNP,varpos,n_all,nInd,nIndloc,ties,lon_over,lon_com_over,merge_conflict)
#ref to %popall: {poploc_ID}{popall_ID}=(popall_seq,popallvar,popall_seq_notmerged)
#ref to %sel_ind: {ind}=1
#ref to %sel_loc: {poplocID}=1
#ref to %sel_gt: {poplocID}{ind}=1
#$indpopallsuffix: suffix of indpopall-outfiles of indpoploc
#$splitlocfilename: split_loci outfile of indpoploc
#$depthfilename: out of depth_analyzer: depth data
#$pairlocgr_fname: out of pair_finder: loc_group n_loc poplocID
#$pairselfm_fname: out of pair_finder: self-matching loci: poplocID n_seqpairs n_ind
#$pair_rep_fname: out of pair_finder: report-file
#$alleles_suffix: suffix of alleles-outfile of indloc
#$loci_suffix = '_loci.txt';#suffix of *_loci outfiles out of indloc

sub select_loci {
	#######################
	#declare and initialize
	#######################
	
	#{
	#_r means a reference
	my ($poploc_r,$popall_r,$sel_ind_r,$sel_loc_r,$sel_gt_r,
	$indpopallsuffix,$splitlocfilename,$depthfilename,
	$pairlocgr_fname,$pairselfm_fname,$pair_rep_fname,$alleles_suffix,$loci_suffix) = @_;
	my $entry = '';#user entry
	my $usersel = 'invalid';#valid selection made by user
	my %number_option = (
	1 => 'pairs',
	2 => 'loclist',
	3 => 'sl',
	4 => 'nInd',
	5 => 'nAll',
	6 => 'nSNP',
	7 => 'nAllind',
	8 => 'ties',
	9 => 'splitloc',
	10 => 'depth',
	11 => 'merge',
	12 => 'fracud',
	13 => 'run',
	14 => 'quit'
	);
	my %onoff = (
	pairs => 'off',#selection based on pair_finder analysis
	loclist => 'off',#user provides a list of loci
	sl => 'off',#sequence length with min max
	nInd => 'off',#number of individuals with min max
	nAll => 'off',#number of alleles with min max
	nSNP => 'off',#number of SNPs with min max
	nAllind => 'off',#number of alleles in allele-richest ind.
	ties => 'off',#loci with/ without ties
	splitloc => 'off',#split loci only / no split loci
	depth => 'off',#selection based on results from depth_analyzer
	merge => 'off',#loci without/with merging conflict
	fracud => 'off',##fraction of used read depth: depth of genotype / total read depth
	run => 'off',#run
	quit => 'off',#quit the program
	invalid => 'off'#dummy for invalid user selection
	);
	my %paired_loc = ();# {loc_group} = (poplocID1,poplocID2);
						#loc_group: group-number out of pair_finder in file loc_groups.txt
	my $loc_group = 0;
	my $loclistfilename = '';#name of a file with a list of loci to select
	my %userloclist = ();# {poplocID}=1
	my $sl = 0;#sequence length
	my $sl_min = 0;#sequence length: min
	my $sl_max = 0;#sequence length: max
	my $nInd = 0;#number of individuals
	my $nInd_min = 0;#number of individuals: min
	my $nInd_max = 0;#number of individuals: max
	my $nAll = 0;#number of alleles
	my $nAll_min = 0;#number of alleles: min
	my $nAll_max = 0;#number of alleles: max
	my $nSNP = 0;#number of SNPs
	my $nSNP_min = 0;#number of SNPs: min
	my $nSNP_max = 0;#number of SNPs: max
	my %indpopall = ();# {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
						#Attention: selection procedures do not update this hash!
						#Always access it via loci in %sel_loc (always up to date)
	my $nAllind = 0;#number of alleles within an individual
	my $nAllindmax = 0;#number of alleles within the allele-richest individual
	my $nAllindmax_min = 0;#minimum number of alleles within the allele-richest individual
	my $nAllindmax_max = 0;#maximum number of alleles within the allele-richest individual
	my @nAllinds = ();#numbers of alleles of all selected inds in one locus
	my $poplocID = 0;#population-level locus ID
	my $ties_sel = 2;#0: select loci without ties, 1: select loci with ties
	my $ties_thisloc = 2;#0: no ties, 1: ties
	my $splitloc_sel = 2;#0: no split loci, 1: split loci only
	my %splitloc = ();# {poplocID}{ind} = indlocID
	my $nIndloc = 0;#number of individual loci in a population-level locus
	my $ind = '';#individual ID
	my $thisind = '';#current individual ID
	my $indlocID = 0;#individual locus ID
	my $popall_ID = 0;#population-level allele ID
	my %popall_IDs_thisloc = ();#all population-level alleles of a locus
								# {popall_ID} = counter
	my $indall_ID = 0;#individual allele ID
	my @locali = ();#2d-array,alignment of a locus' alleles, d1: rows, d2: cols
	my @varpos = ();#variable positions of a locus, count starting with 0
	my $infilename = '';#an infilename
	my $fracud = 0;#fraction of used read depth: depth of genotype / total read depth
	my $fracud_min = 0;#minimum for fraction of used read depth
	my $fracud_del = 0;#1 if a locus must be excluded because of fracud
	my $Alldep = 0;#depth of an allele
	my %Alldep_thisind = ();#for one ind: {indall_ID} = Alldep
	my %indlocdep_this_ind = ();# {indlocID} = locdep
	my %indpoplocdep = ();# {poplocID}{ind}=locdep
	my $depth = 0;#depth of a genotype (sum of estimated depth of alleles)
	my $totrd = 0;#total read depth of a poplocus in an individual
	#variables for selection based on depth analysis
	my $daslrNo = 0;#seqlength-range-number
	my $daslmin = 0;#min seqlength
	my $daslmax = 0;#max seqlength 
	my $dacrit = '';#selection criterion, e.g. med_dep_perc min
	my $dacol = 0;#column in file da_all_loc.txt out of depth_analyzer
	my $davalmin = 0;#minimum of value-range
	my $davalmax = 0;#maximum of value-range
	my %da_slr_var = ();# {$daslrNo}: criteria as string for display
	my %da_sl_var = ();# {$sl} = ($dacol,$davalmin,$davalmax) criteria for analysis
	my %number_dacrit = (#user selects a number (key) for a criterion (value):
	1 => 'min_scaldep max',
	2 => 'med_scaldep max',
	3 => 'max_scaldep max',
	4 => 'min_dep_perc min',
	5 => 'med_dep_perc min',
	6 => 'max_dep_perc min'
	);
	my %number_dacol = (#still same user selection of a number(key), val: column in file *_da_allind.txt out of depth_analyzer
	1 => '1',
	2 => '2',
	3 => '3',
	4 => '4',
	5 => '5',
	6 => '6'
	);
	my $da_addrange = 'Y';#Y: add a range with criteria, N: done
	my %depth_data = ();# {poplocID} = (sl,min_scaldep,med_scaldep,max_scaldep,
									#min_dep_perc,med_dep_perc,max_dep_perc)
	my $merge_sel = 2;#0: loci without merging conflict, 1: loci with merging conflict
	my $merge_thisloc = 2;#0: no merging conflict, 1: merging conflict
	#variables about dataset-dimensions:
	my $n_ind_dataset = 0;#number of individuals
	my $n_loc_dataset = 0;#number of loci
	my $n_gt_dataset = 0;#number of genotypes
	my $tempstring1 = '';
	my @temparr1 = ();
	my $i = 0;
	#}

	#####################
	#User selects options
	#####################
	
	#{
	while (($onoff{'run'} eq 'off') and ($onoff{'quit'} eq 'off')) {
		print "\nSELECT LOCI:\n\n";
		print " 1 Paired loci: select one of each pair    $onoff{'pairs'}\n";
		print " 2 List of loci in a file                  $onoff{'loclist'}";
		if ($onoff{'loclist'} eq 'on') {
			print " $loclistfilename\n";
		} else {print "\n"}
		print " 3 Sequence length                         $onoff{'sl'}";
		if ($onoff{'sl'} eq 'on') {
			print " min: $sl_min max: $sl_max\n";
		} else {print "\n"}
		print " 4 Number of individuals                   $onoff{'nInd'}";
		if ($onoff{'nInd'} eq 'on') {
			print " min: $nInd_min max: $nInd_max\n";
		} else {print "\n"}
		print " 5 Number of alleles                       $onoff{'nAll'}";
		if ($onoff{'nAll'} eq 'on') {
			print " min: $nAll_min max: $nAll_max\n";
		} else {print "\n"}
		print " 6 Number of SNPs                          $onoff{'nSNP'}";
		if ($onoff{'nSNP'} eq 'on') {
			print " min: $nSNP_min max: $nSNP_max\n";
		} else {print "\n"}
		print " 7 Numb. of alleles in allele-richest ind. $onoff{'nAllind'}";
		if ($onoff{'nAllind'} eq 'on') {
			print " min: $nAllindmax_min max: $nAllindmax_max\n";
		} else {print "\n"}
		print " 8 Ties                                    $onoff{'ties'}";
		if ($onoff{'ties'} eq 'on') {
			if ($ties_sel == 0) {
				print " loci without ties\n";
			} else {print " loci with ties\n"}
		} else {print "\n"}
		print " 9 Split loci                              $onoff{'splitloc'}";
		if ($onoff{'splitloc'} eq 'on') {
			if ($splitloc_sel == 0) {
				print " no split loci\n";
			} else {print " split loci only\n"}
		} else {print "\n"}
		print "10 Selection based on depth analysis       $onoff{'depth'}";
		if ($onoff{'depth'} eq 'on') {
			print " Select 10 to view criteria.\n";
		} else {print "\n"}
		print "11 Merging conflicts                       $onoff{'merge'}";
		if ($onoff{'merge'} eq 'on') {
			if ($merge_sel == 0) {
				print " loci without merging conflict\n";
			} else {print " loci with merging conflict\n"}
		} else {print "\n"}
		print "12 Used fraction of read depth             $onoff{'fracud'}";
		if ($onoff{'fracud'} eq 'on') {
			print " min: $fracud_min\n";
		} else {print "\n"}
		print "\n13 Run\n";
		print "14 Quit\n\n";
		print "\nPlease select a number to activate/inactivate an option: ";
		$entry = <STDIN>;#get user-entry
		chomp $entry;
		if (defined $number_option{$entry}) {#if it is valid
			$usersel = $number_option{$entry};#take over as selection
		} else {$usersel = 'invalid'};#else set selection to 'invalid'
		
		#switch corresponding switch
		if ($onoff{$usersel} eq 'off') {
			$onoff{$usersel} = 'on';
		} else {$onoff{$usersel} = 'off'};
		
		#if user activated pairs, give info, open files
		if (($usersel eq 'pairs') and ($onoff{'pairs'} eq 'on')) {
			print
			"\nPaired loci: select one of each pair:\n\n",
			"Selection of loci based on analysis with pair_finder:\n",
			"-selects one locus of each group with 2 loci,\n",
			"-unselects all other loci.\n\n";
			#open pairfinder report-file
			unless(open(PAIRSREP, $pair_rep_fname)) {
				print "Can't open $pair_rep_fname.\n";
				$onoff{'pairs'} = 'off';
			}
			#open datafile with locus groups out of pair_finder
			unless(open(PAIRSDATA, $pairlocgr_fname)) {
				print "Can't open $pairlocgr_fname.\n";
				$onoff{'pairs'} = 'off';
			}
			#open datafile with self-matching loci out of pair_finder
			unless(open(PAIRSELFMATCH, $pairselfm_fname)) {
				print "Can't open $pairselfm_fname.\n";
				$onoff{'pairs'} = 'off';
			}			
			#open file with split loci out of indpoploc
			unless(open(SPLITLOCLIST, $splitlocfilename)) {
				print "Can't open $splitlocfilename.\n";
				$onoff{'pairs'} = 'off';
			}			
			#if option still active, continue
			if ($onoff{'pairs'} eq 'on') {
				#Read first line of report, check if pair_finder was executed
				#based on a preselection of data
				$tempstring1 = <PAIRSREP>;
				#Infilecolumns!
				if ($tempstring1 =~ /preselection/) {
					print "Pair_finder was executed based on a preselection of data.\n",
					"It could be important to verify that the preselection has not changed\n",
					"since you executed pair_finder.\n",
					"You could look into these files:\n",
					"export/sel_rep.txt    $pair_rep_fname\n",
					"Please close them again before you continue with data_selector.\n";
					$entry = 'invalid';
					while ($entry eq 'invalid') {
						print "Inactivate selection of paired loci? (Y/N): ";
						$entry = <STDIN>;
						chomp $entry;
						unless (($entry eq 'Y') or ($entry eq 'N')) {
							$entry = 'invalid';
						}
					}
					if ($entry eq 'Y') {
						$onoff{'pairs'} = 'off';
					}
				}
				close PAIRSREP;
			}
			#if option still active, continue
			if ($onoff{'pairs'} eq 'on') {
				print "Self-matching loci identified by pair_finder will be unselected.\n",
				"They do not appear in groups identified by pair_finder.\n";
				#automatically activate splitloc
				print "Automatically activating Split loci: no split loci.\n",
				"Split loci will be unselected.\n",
				"Split loci were ignored by pair_finder\n",
				"and do not appear in groups identified by pair_finder.\n";
				$onoff{'splitloc'} = 'on';
				$splitloc_sel = 0;
			}
		}
		#if user inactivated pairs, close data-files out of pair_finder
		#and split_loci.txt out of indpoploc
		elsif (($usersel eq 'pairs') and ($onoff{'pairs'} eq 'off')) {
			close PAIRSDATA;
			close PAIRSELFMATCH;
			close SPLITLOCLIST;
			#automatically inactivate splitloc
			$onoff{'splitloc'} = 'off';
			print "Inactivated Paired loci and Split loci.\n";
		}
		#if user activated loclist, get filename
		elsif (($usersel eq 'loclist') and ($onoff{'loclist'} eq 'on')) {
			print
			"\nList of loci: I need a textfile:\n",
			"One locus ID per line, no spaces, no blank lines.\n\n",
			"Please enter filename: ";
			$loclistfilename = <STDIN>;
			chomp $loclistfilename;
			#open the file or inactivate the option again
			unless(open(LOCLIST, $loclistfilename)) {
				print "Can't open $loclistfilename.\n";
				$onoff{'loclist'} = 'off';
			}
		}
		#if user inactivated loclist, close file
		elsif (($usersel eq 'loclist') and ($onoff{'loclist'} eq 'off')) {
			close LOCLIST;
		}
		#if user activated sl, get min and max
		elsif (($usersel eq 'sl') and ($onoff{'sl'} eq 'on')){
			print "\nSequence length: please enter min: ";
			$sl_min = <STDIN>;
			chomp $sl_min;
			print "Sequence length: please enter max: ";
			$sl_max = <STDIN>;
			chomp $sl_max;
			#Warn and turn off if values are invalid
			unless(($sl_min =~ /^\d+$/) and ($sl_max =~ /^\d+$/) and ($sl_max >= $sl_min) and ($sl_max > 0)) {
				print "Min must be 0 or positive integer.\n",
				"Max must be positive integer >= min.\n";
				$onoff{'sl'} = 'off';
			}
		}
		#if user activated nInd, get min and max
		elsif (($usersel eq 'nInd') and ($onoff{'nInd'} eq 'on')){
			print "\nNumber of individuals: please enter min: ";
			$nInd_min = <STDIN>;
			chomp $nInd_min;
			print "Number of individuals: please enter max: ";
			$nInd_max = <STDIN>;
			chomp $nInd_max;
			#Warn and turn off if values are invalid
			unless(($nInd_min =~ /^\d+$/) and ($nInd_max =~ /^\d+$/) and ($nInd_max >= $nInd_min) and ($nInd_max > 0)) {
				print "Min must be 0 or positive integer.\n",
				"Max must be positive integer >= min.\n";
				$onoff{'nInd'} = 'off';
			}
		}
		#if user activated nAll, get min and max
		elsif (($usersel eq 'nAll') and ($onoff{'nAll'} eq 'on')){
			print "\nNumber of alleles: please enter min: ";
			$nAll_min = <STDIN>;
			chomp $nAll_min;
			print "Number of alleles: please enter max: ";
			$nAll_max = <STDIN>;
			chomp $nAll_max;
			#Warn and turn off if values are invalid
			unless(($nAll_min =~ /^\d+$/) and ($nAll_max =~ /^\d+$/) and ($nAll_max >= $nAll_min) and ($nAll_max > 0)) {
				print "Min must be 0 or positive integer.\n",
				"Max must be positive integer >= min.\n";
				$onoff{'nAll'} = 'off';
			}
		}
		#if user activated nSNP, get min and max
		elsif (($usersel eq 'nSNP') and ($onoff{'nSNP'} eq 'on')){
			print "\nNumber of SNPs: please enter min: ";
			$nSNP_min = <STDIN>;
			chomp $nSNP_min;
			print "Number of SNPs: please enter max: ";
			$nSNP_max = <STDIN>;
			chomp $nSNP_max;
			#Warn and turn off if values are invalid
			unless(($nSNP_min =~ /^\d+$/) and ($nSNP_max =~ /^\d+$/) and ($nSNP_max >= $nSNP_min) and ($nSNP_max > 0)) {
				print "Min must be 0 or positive integer.\n",
				"Max must be positive integer >= min.\n";
				$onoff{'nSNP'} = 'off';
			}
		}
		#if user activated nAllind, get min and max
		elsif (($usersel eq 'nAllind') and ($onoff{'nAllind'} eq 'on')){
			print "\nNumb. of alleles in allele-richest ind. Please enter min: ";
			$nAllindmax_min = <STDIN>;
			chomp $nAllindmax_min;
			print "Numb. of alleles in allele-richest ind. Please enter max: ";
			$nAllindmax_max = <STDIN>;
			chomp $nAllindmax_max;
			#Warn and turn off if values are invalid
			unless(($nAllindmax_min =~ /^\d+$/) and ($nAllindmax_max =~ /^\d+$/) and ($nAllindmax_max >= $nSNP_min) and ($nAllindmax_max > 0)) {
				print "Min must be 0 or positive integer.\n",
				"Max must be positive integer >= min.\n";
				$onoff{'nAllind'} = 'off';
			}
		}
		#if user activated ties, get selection criterion
		elsif (($usersel eq 'ties') and ($onoff{'ties'} eq 'on')){
			print "\nTies:\n",
			"0 Loci without ties\n", "1 Loci with ties\n\n",
			"Please select an option (number): ";
			$ties_sel = <STDIN>;
			chomp $ties_sel;
			#Warn and turn off if value is invalid
			unless(($ties_sel =~ /^\d+$/) and (($ties_sel == 0) or ($ties_sel == 1))) {
				print "Invalid selection.\n";
				$onoff{'ties'} = 'off';
			}
		}
		#if user activated splitloc, open file or inactivate again
		elsif (($usersel eq 'splitloc') and ($onoff{'splitloc'} eq 'on')){
			print "\nSplit loci:\n",
			"0 No split loci\n", "1 Split loci only\n\n",
			"Please select an option (number): ";
			$splitloc_sel = <STDIN>;
			chomp $splitloc_sel;
			#Warn and turn off if value is invalid
			unless(($splitloc_sel =~ /^\d+$/) and (($splitloc_sel == 0) or ($splitloc_sel == 1))) {
				print "Invalid selection.\n";
				$onoff{'splitloc'} = 'off';
			}
			#if still active, try to open file
			if ($onoff{'splitloc'} eq 'on') {
				unless(open(SPLITLOCLIST, $splitlocfilename)) {
					print "Can't open $splitlocfilename.\n";
					$onoff{'splitloc'} = 'off';
				}
			}			
		}
		#if user wants to inactivate splitloc...
		elsif (($usersel eq 'splitloc') and ($onoff{'splitloc'} eq 'off')){
			#if pairs is active, don't inactivate splitloc
			if ($onoff{'pairs'} eq 'on') {
				print "Paired loci is active. I can't inactivate Split loci.\n";
				$onoff{'splitloc'} = 'on';
			} else {#if pairs is not active, close file split_loci.txt
				close SPLITLOCLIST;
			}
		}		
		#if user activated depth, get criteria
		elsif (($usersel eq 'depth') and ($onoff{'depth'} eq 'on')){
			print "\nSelection based on depth analysis:\n\n",
			"You can define criteria for different sequence length ranges.\n",
			"Ranges may not overlap.\n";
			$daslrNo = 1;#set range number to 1
			#display criteria menu until user selected not to add criteria
			while ($da_addrange eq 'Y') {
				#show already existing criteria
				if (keys %da_slr_var > 0) {
					print "\nAlready defined:\n\n";
				}
				for $i (sort {$a <=> $b} keys %da_slr_var) {
					print "$da_slr_var{$i}\n";
				}
				#ask if user wants to add a range
				print "\nAdd a new range with a criterion? (Y/N): ";
				$entry = <STDIN>;
				chomp $entry;
				if ($entry eq 'N') {#user doesn't want to add
					$da_addrange = 'N';#stopcriterion
					next;#stop
				}
				elsif ($entry ne 'Y') {#invalid entry
					next;
				}				
				print "Please define a sequence-length range:\nmin: ";
				$daslmin = <STDIN>;
				chomp $daslmin;
				print "max: ";
				$daslmax = <STDIN>;
				chomp $daslmax;
				#check if range limits are valid (min non-negative integer, max pos. integer >= min)
				unless (($daslmin =~ /^\d+$/) and ($daslmax =~ /^\d+$/) and ($daslmax >= $daslmin) and ($daslmax > 0)) {
					print "invalid range\n";
					next;
				}
				#check if range does not overlap with already defined range
				if ((defined $da_sl_var{$daslmin}) or (defined $da_sl_var{$daslmax})) {
					print "Range overlaps with already defined range.\n";
					next;
				}				
				#display criteria until user made a valid selection
				$dacrit = 'none';#set criterion to none
				while ($dacrit eq 'none') {
					print "1 min_scaldep max\n2 med_scaldep max\n3 max_scaldep max\n",
					"4 min_dep_perc min\n5 med_dep_perc min\n6 max_dep_perc min\n\n",
					"Please select a criterion (number): ";
					$entry = <STDIN>;
					chomp $entry;
					if (defined $number_dacrit{$entry}) {#if entry is valid
						$dacrit = $number_dacrit{$entry};#take over criterion
						$dacol = $number_dacol{$entry};#determine column to look up in file da_all_loc.txt
					}					
				}
				#get a minimum or maximum, whatever is appropriate
				$davalmin = -1;#set min value to -1
				while ($davalmin == -1) {
					print "Please enter value: $dacrit: ";
					$entry = <STDIN>;
					chomp $entry;
					#check if entry is plausible
					if ($dacrit =~ /scaldep/) {#if the criterion is based on scaled depth (0-1)
						unless ((($entry =~ /^\d$/) or ($entry =~ /^\d\.\d+$/)) and
								($entry >= 0) and ($entry <= 1)) {
							print "invalid entry\n";
							next;
						}
					}
					elsif ($dacrit =~ /dep_perc/) {#if the criterion is based on depth percentile (0-100)
						unless ((($entry =~ /^\d+$/) or ($entry =~ /^\d+\.\d+$/)) and
								($entry >= 0) and ($entry <= 100)) {
							print "invalid entry\n";
							next;
						}
					}
					#entry is valid: take over
					if ($dacrit =~ /scaldep/) {#if the criterion is based on scaled depth (0-1)
						$davalmax = $entry;#value is a maximum
						$davalmin = 0;#minimum is 0
						$da_slr_var{$daslrNo} = 'length '.$daslmin.'-'.$daslmax.' '.$dacrit.' '.$davalmax;
						++$daslrNo;
					}
					elsif ($dacrit =~ /dep_perc/) {#if the criterion is based on scaled depth (0-1)
						$davalmin = $entry;#value is a minimum
						$davalmax = 100;#maximum is 100
						$da_slr_var{$daslrNo} = 'length '.$daslmin.'-'.$daslmax.' '.$dacrit.' '.$davalmin;
						++$daslrNo;
					}
				}
				#populate %da_sl_var
				for ($sl = $daslmin; $sl <= $daslmax; ++$sl) {
					@{$da_sl_var{$sl}} = ($dacol,$davalmin,$davalmax);
				}
			}
			#if user did not define any valid criterion, set to inactive again
			if (keys %da_slr_var == 0) {
				$onoff{'depth'} = 'off';
			} else {#User made valid selections, try to open file with depth data
				unless(open(DEPTHDATA, $depthfilename)) {#open file or die
					print "Can't open $depthfilename.\n",
					"This file is needed for selection of loci based on depth analysis.\n",
					"Please run depth_analyzer to produce it.\n";
					$onoff{'depth'} = 'off';
					$da_addrange = 'Y';#set back to Y to enable new selection
					%da_slr_var = ();#set back criteria for display
					%da_sl_var = ();#set back criteria for analysis					
				}
			}
		}
		#if user inactivated depth, show criteria and ask if he wants to inactivate
		elsif (($usersel eq 'depth') and ($onoff{'depth'} eq 'off')){
			print "\nSelection based on depth analysis:\n";
			#show already existing criteria
			if (keys %da_slr_var > 0) {
				print "\nAlready defined:\n\n";
			}
			for $i (sort {$a <=> $b} keys %da_slr_var) {
				print "$da_slr_var{$i}\n";
			}
			#ask if user wants to inactivate until you get a valid entry
			$entry = 'invalid';
			while ($entry eq 'invalid') {
				print "\nInactivate selection based on depth analysis (Y/N): ";
				$tempstring1 = <STDIN>;
				chomp $tempstring1;
				if (($tempstring1 eq 'Y') or ($tempstring1 eq 'N')) {
					$entry = $tempstring1;
				}
			}
			if ($entry eq 'N') {#user doesn't want to inactivate
				$onoff{'depth'} = 'on';
			} else {#user wants to inactivate the option
				%da_slr_var = ();#set back criteria for display
				%da_sl_var = ();#set back criteria for analysis
				close DEPTHDATA;
				$da_addrange = 'Y';#set back to Y to enable new selection
			}			
		}
		#if user activated merge, get selection criterion
		elsif (($usersel eq 'merge') and ($onoff{'merge'} eq 'on')) {
			print "\nMerging conflicts:\n",
			"0 Loci without merging conflicts\n", "1 Loci with merging conflicts\n\n",
			"Please select an option (number): ";
			$merge_sel = <STDIN>;
			chomp $merge_sel;
			#Warn and turn off if value is invalid
			unless(($merge_sel =~ /^\d+$/) and (($merge_sel == 0) or ($merge_sel == 1))) {
				print "Invalid selection.\n";
				$onoff{'merge'} = 'off';
			}			
		}
		#if user activated fraction of used read depth, get min
		elsif (($usersel eq 'fracud') and ($onoff{'fracud'} eq 'on')){
			print "\nUsed fraction of read depth:\n",
			"I will exclude a locus when any genotype has a value below min\n",
			"please enter min (decimal between 0 and 1): ";
			$fracud_min = <STDIN>;
			chomp $fracud_min;
			#Warn and turn off if value is invalid (must be between 0 and 1)
			unless((($fracud_min =~ /^\d$/) or ($fracud_min =~ /^\d\.\d+$/)) and
			($fracud_min >= 0) and ($fracud_min <= 1)) {
				print "Min must be >= 0 and <= 1.\n";
				$onoff{'fracud'} = 'off';
			}
		}
		#if user selected quit, ask for confirmation and quit if appropriate
		elsif ($usersel eq 'quit') {
			print "Do you really want to quit? (Y/N): ";
			$entry = <STDIN>;
			chomp $entry;
			if ($entry eq 'Y') {
				print "No data selected. Exiting..\n";
				exit;
			} else {$onoff{'quit'} = 'off'};
		}		
	}
	#}
	
	##########################################################
	#If needed, read in *_indpopall.txt files out of indpoploc
	#into %indpopall, remove unselected data
	##########################################################
	
	#{
	if (($onoff{'nAll'} eq 'on') or
	($onoff{'nSNP'} eq 'on') or
	($onoff{'nAllind'} eq 'on') or
	($onoff{'ties'} eq 'on') or
	($onoff{'fracud'} eq 'on')) {
		#loop through selected individuals
		for $ind (keys %{$sel_ind_r}) {
			$infilename = $ind.$indpopallsuffix;#build infilename
			unless(open(INFILE, $infilename)) {#open file or die
				print "Can't open $infilename. Exiting..\n";
				exit;
			}
			$tempstring1 = <INFILE>;#skip headerline
			while ($tempstring1 = <INFILE>) {
				chomp $tempstring1;
				@temparr1 = split(/\t/,$tempstring1);
				#load data into %indpopall {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
				#Infilecolumns!
				push @{$indpopall{$temparr1[2]}{$ind}{$temparr1[0]}{$temparr1[3]}}, $temparr1[1];
			}
			close INFILE;
		}
		#remove unselected data from %indpopall
		#loop through poplocIDs in %indpopall
		for $poplocID (keys %indpopall) {
			if (defined $$sel_loc_r{$poplocID}) {#if locus is selected
				#loop through individuals in %indpopall{$poplocID}
				for $ind (keys %{$indpopall{$poplocID}}) {
					#if this individual genotype (poplocID x ind) is not selected
					unless ((defined $$sel_gt_r{$poplocID}) and (defined $$sel_gt_r{$poplocID}{$ind})) {
						delete $indpopall{$poplocID}{$ind};
					}
				}
			} else {#if locus is not selected
				delete $indpopall{$poplocID};
			}
		}
	}
	#}
	
	#####################################################
	#If needed, read in *_alleles.txt files out of indloc
	#Load Alldep into %indpopall
	#####################################################
	
	#{
	if ($onoff{'fracud'} eq 'on') {
		for $ind (keys %{$sel_ind_r}) {
			$infilename = $ind.$alleles_suffix;#build infilename
			unless(open(INFILE, $infilename)) {#open file or die
				print "Can't open $infilename. Exiting..\n";
				exit;
			}
			$tempstring1 = <INFILE>;#skip headerline
			while ($tempstring1 = <INFILE>) {
				chomp $tempstring1;
				@temparr1 = split(/\t/,$tempstring1);
				#Load data into %Alldep_thisind: {indall_ID} = Alldep
				#Infilecolumns!
				if ($temparr1[1] eq 'valid') {
					$Alldep_thisind{$temparr1[3]} = $temparr1[5];
				}
			}
			close INFILE;
			#loop through poplocIDs in %indpopall
			for $poplocID (keys %indpopall) {
				#loop through indlocIDs
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						#look up indall_ID
						$indall_ID = $indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[0];
						#look up Alldep
						$Alldep = $Alldep_thisind{$indall_ID};
						#add to %indpopall
						$indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[1] = $Alldep;
					}
				}
			}
			#set back and next ind
			%Alldep_thisind = ();
		}
	}
	#}
	
	##################################################
	#If needed, read in *_loci.txt files out of indloc
	#populate %indpoplocdep
	##################################################
	
	#{
	if ($onoff{'fracud'} eq 'on') {
		#loop through selected individuals
		for $ind (keys %{$sel_ind_r}) {
			$infilename = $ind.$loci_suffix;
			unless(open(INFILE, $infilename)) {
				print "Can't open $infilename. Exiting..\n";
				exit;
			}
			$tempstring1 = <INFILE>;#skip headerline
			while ($tempstring1 = <INFILE>) {
				chomp $tempstring1;
				@temparr1 = split(/\t/,$tempstring1);
				#load data into %indlocdep_this_ind
				#Infilecolumns!
				if ($temparr1[1] eq 'valid') {
					$indlocdep_this_ind{$temparr1[0]} = $temparr1[3];
				}
			}
			close INFILE;
			#loop through selected poplocs and corresponding indlocs for this ind
			for $poplocID (keys %indpopall) {
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#load depth into %indpoplocdep
					$indpoplocdep{$poplocID}{$ind} += $indlocdep_this_ind{$indlocID};
				}
			}
			%indlocdep_this_ind = ();#set back and next ind
		}	
	}
	#}
	
	################################################################
	#Determine starting dimensions of database or previous selection
	#Append to file export/sel_rep.txt
	################################################################
	
	#{
	#n_ind
	$n_ind_dataset = keys %{$sel_ind_r};
	#n_loc
	$n_loc_dataset = keys %{$sel_loc_r};
	#n_gt
	for $poplocID (keys %{$sel_gt_r}) {
		$n_gt_dataset += keys %{$$sel_gt_r{$poplocID}};
	}
	#open export/sel_rep.txt for append
	unless(open(OUTREP, ">>export/sel_rep.txt")) {
		print "Cannot open file export/sel_rep.txt, exiting ...\n\n";
		exit;
	}
	print OUTREP
	"Selection of Loci: starting from\n\n",
	"$n_ind_dataset individuals\n$n_loc_dataset loci\n$n_gt_dataset genotypes\n\n",
	"Selected loci according to these criteria:\n\n";	
	#}
	
	########################
	#Pairs: if activated, do
	########################
	
	#{
	if ($onoff{'pairs'} eq 'on') {
		#read in file loc_groups: loc_group n_loc poplocID
		$tempstring1 = <PAIRSDATA>;#skip headerline
		while ($tempstring1 = <PAIRSDATA>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#If current locus ([2]) is in a group with 1 or more than 2 members ([1])
			#Infilecolums!
			if ($temparr1[1] != 2) {
				delete $$sel_loc_r{$temparr1[2]};#delete from selected loci
				delete $$sel_gt_r{$temparr1[2]};#delete from selected genotypes			
			} else {#if current locus is part of a pair
				#store in %paired_loc: {loc_group}=(poplocID1,poplocID2)
				push @{$paired_loc{$temparr1[0]}}, $temparr1[2];
			}
		}
		close PAIRSDATA;
		#loop through %paired_loc
		for $loc_group (keys %paired_loc) {
			#if the first locus in this group is defined in %sel_loc
			if (defined $$sel_loc_r{$paired_loc{$loc_group}[0]}) {
				#delete the second
				delete $$sel_loc_r{$paired_loc{$loc_group}[1]};#delete from selected loci
				delete $$sel_gt_r{$paired_loc{$loc_group}[1]};#delete from selected genotypes
			}
		}
		#unselect self-matching loci identified by pair_finder
		#read in file about self-matching loci
		$tempstring1 = <PAIRSELFMATCH>;#read to skip headerline
		while ($tempstring1 = <PAIRSELFMATCH>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#Infilecolumns! poplocID in first col
			delete $$sel_loc_r{$temparr1[0]};#delete from selected loci
			delete $$sel_gt_r{$temparr1[0]};#delete from selected genotypes		
		}		
		#print to report file
		print OUTREP "Paired loci: select one of each pair\n";
	}
	#}
	
	##########################
	#loclist: if activated, do
	##########################
	
	#{
	if ($onoff{'loclist'} eq 'on') {
		#read in user loclist
		while ($poplocID = <LOCLIST>) {
			chomp $poplocID;
			$userloclist{$poplocID} = 1;
		}
		#loop through selected loci
		for $poplocID (keys %{$sel_loc_r}) {
			#if locus is not on user-list
			unless (defined $userloclist{$poplocID}) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}
		}
		#print to report file
		print OUTREP "List of loci in file $loclistfilename\n";
	}
	#}
	
	#####################
	#sl: if activated, do
	#####################
	
	#{
	if ($onoff{'sl'} eq 'on') {
		for $poplocID (keys %{$sel_loc_r}) {#loop through selected loci
			$sl = $$poploc_r{$poplocID}[0];#look up sequence length in %poploc
			#delete locus if sequence length is not in selected range
			unless (($sl >= $sl_min) and ($sl <= $sl_max)) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}
		}
		#print to report file
		print OUTREP "Sequence length: min $sl_min max $sl_max\n";
	}
	#}
	
	#######################
	#nInd: if activated, do
	#######################
	
	#{
	if ($onoff{'nInd'} eq 'on') {
		for $poplocID (keys %{$sel_loc_r}) {#loop through loci in %sel_loc
			$nInd = keys %{$$sel_gt_r{$poplocID}};#determine nInd in %sel_gt
			#delete locus if nInd is not in selected range
			unless (($nInd >= $nInd_min) and ($nInd <= $nInd_max)) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}			
		}
		#print to report file
		print OUTREP "Number of individuals: min $nInd_min max $nInd_max\n";
	}	
	#}

	#######################
	#nAll: if activated, do
	#######################
	
	#{
	if ($onoff{'nAll'} eq 'on') {
		#loop through selected poplocIDs in %sel_loc: {poplocID}=1
		for $poplocID (keys %{$sel_loc_r}) {
			#loop through individuals in %indpopall: {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
			for $ind (keys %{$indpopall{$poplocID}}) {
				#loop through indlocIDs in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs in %indpopall
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						#collect popall_IDs in %popall_IDs_thisloc: {popall_ID} = counter
						++ $popall_IDs_thisloc{$popall_ID};
					}
				}
			}
			$nAll = keys %popall_IDs_thisloc;#determine number of alleles
			#delete locus if nAll is not in selected range
			unless (($nAll >= $nAll_min) and ($nAll <= $nAll_max)) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}
			if ($nAll == 0) {#If there are no alleles
				print "Found no alleles in locus $poplocID.\n",
				"Something must be wrong with previous selection files.\n";
			}
			%popall_IDs_thisloc = ();#set back and next locus
		}
		#print to report file
		print OUTREP "Number of alleles: min $nAll_min max $nAll_max\n";
	}
	#}
	
	#######################
	#nSNP: if activated, do
	#######################
	
	#{
	if ($onoff{'nSNP'} eq 'on') {
		#loop through selected poplocIDs in %sel_loc: {poplocID}=1
		for $poplocID (keys %{$sel_loc_r}) {
			#loop through individuals in %indpopall: {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
			for $ind (keys %{$indpopall{$poplocID}}) {
				#loop through indlocIDs in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs in %indpopall
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						#collect popall_IDs in %popall_IDs_thisloc: {popall_ID} = counter
						++ $popall_IDs_thisloc{$popall_ID};
					}
				}
			}
			$nAll = keys %popall_IDs_thisloc;#determine number of alleles
			if ($nAll == 1) {#if there is one allele
				$nSNP = 0;#There can't be SNPs
			}
			elsif ($nAll > 1) {#if there is more than one allele
				#build an alignment of the alleles as 2d-array @locali
				#loop through popall_IDs
				for $popall_ID (keys %popall_IDs_thisloc) {
					#look up sequence and store in string
					$tempstring1 = $$popall_r{$poplocID}{$popall_ID}[0];
					@temparr1 = split('',$tempstring1);#split into array
					push @locali, [@temparr1];#add to @locali
				}
				#call sub varpos to determine variable positions (count starting with 0)
				@varpos = varpos(\@locali);
				$nSNP = @varpos;#number of SNPs				
			} else {#Something must be wrong
				print "Found no alleles in locus $poplocID.\n",
				"Something must be wrong with previous selection files.\n";
			}			
			#delete locus if nSNP is not in selected range
			unless (($nSNP >= $nSNP_min) and ($nSNP <= $nSNP_max)) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}	
			#set variables back
			%popall_IDs_thisloc = ();
			@locali = ();
			#next locus
		}		
		#print to report file
		print OUTREP "Number of SNPs: min $nSNP_min max $nSNP_max\n";
	}
	#}
	
	##########################
	#nAllind: if activated, do
	##########################
	
	#{
	if ($onoff{'nAllind'} eq 'on') {
		#loop through selected poplocIDs in %sel_loc: {poplocID}=1
		for $poplocID (keys %{$sel_loc_r}) {
			#loop through individuals in %indpopall: {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
			for $ind (keys %{$indpopall{$poplocID}}) {
				#loop through indlocIDs in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs in %indpopall
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
					#count alleles of this ind
					++$nAllind;
					}
				}
				#add number of alleles for this ind to @nAllinds
				push @nAllinds, $nAllind;
				$nAllind = 0;#set back and next ind
			}
			#get number of alleles of allele-richest ind
			@nAllinds = sort {$a <=> $b} @nAllinds;
			$nAllindmax = pop @nAllinds;
			#delete locus if nAllindmax is not in selected range
			unless (($nAllindmax >= $nAllindmax_min) and ($nAllindmax <= $nAllindmax_max)) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}			
			@nAllinds = ();#set back and next locus
		}	
		#print to report file
		print OUTREP "Numb. of alleles in allele-richest ind.: min $nAllindmax_min max $nAllindmax_max\n";
	}
	#}
	
	#######################
	#ties: if activated, do
	#######################
	
	#{
	if ($onoff{'ties'} eq 'on') {
		$nInd = 0;
		$nIndloc = 0;
		#loop through selected poplocIDs in %sel_loc: {poplocID}=1
		for $poplocID (keys %{$sel_loc_r}) {
			#loop through individuals in %indpopall: {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
			for $ind (keys %{$indpopall{$poplocID}}) {
				++$nInd;#count the individual
				#loop through indlocIDs in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					++$nIndloc;#count the individual locus
				}
			}
			if ($nIndloc > $nInd) {#if there are more indloci than individuals
				$ties_thisloc = 1;#locus ties
			} else {#if not
				$ties_thisloc = 0;#locus doesn't tie
			}
			#delete locus if it doesn't match the selected criterion
			unless($ties_thisloc == $ties_sel) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}
			#set variables back
			$nInd = 0;
			$nIndloc = 0;
			#next locus
		}
		#print to report file
		if ($ties_sel == 0) {
			print OUTREP "Loci without ties.\n";
		} else {
			print OUTREP "Loci with ties.\n";
		}
	}	
	#}
	
	###########################
	#splitloc: if activated, do
	###########################
	
	#{
	if ($onoff{'splitloc'} eq 'on') {
		#read in split_loci file into %splitloc: {poplocID}{ind} = indlocID
		$tempstring1 = <SPLITLOCLIST>;#read to skip headerline
		while ($tempstring1 = <SPLITLOCLIST>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#Infilecolumns!
			$splitloc{$temparr1[0]}{$temparr1[1]} = $temparr1[2];
		}
		#if selected, delete all split loci
		if ($splitloc_sel == 0) {
			for $poplocID (keys %splitloc) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}
		}
		#else if selected, delete all non split loci
		elsif ($splitloc_sel == 1) {
			for $poplocID (keys %{$sel_loc_r}) {
				unless (defined $splitloc{$poplocID}) {
					delete $$sel_loc_r{$poplocID};#delete from selected loci
					delete $$sel_gt_r{$poplocID};#delete from selected genotypes
				}
			}
		}
		else {#Else: error
			print "Error with splitloc.\n";
		}	
		#print to report file
		if ($splitloc_sel == 0) {
			print OUTREP "No split loci.\n";
		}else {
			print OUTREP "Split loci only.\n";
		}
	}
	#}
	
	########################
	#depth: if activated, do
	########################
	
	#{
	if ($onoff{'depth'} eq 'on') {
		#read in file with depth data, populate %depth_data
		$tempstring1 = <DEPTHDATA>;#skip headerline
		while ($tempstring1 = <DEPTHDATA>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			$poplocID = shift @temparr1;#take out poplocID
			@{$depth_data{$poplocID}} = @temparr1;#store data in %depth_data
		}
		#verify loci for depth criteria
		$i = 0;#counts selected loci without available depth data
		for $poplocID (keys %{$sel_loc_r}) {
			if (defined $depth_data{$poplocID}) {#if depth data are available
				$sl = $depth_data{$poplocID}[0];#look up sequence length of the locus
				@temparr1 = @{$da_sl_var{$sl}};#look up criteria for this seq-length
				#look up value in dacol for this locus and check if acceptable
				unless (($depth_data{$poplocID}[$temparr1[0]] >= $temparr1[1]) and
						($depth_data{$poplocID}[$temparr1[0]] <= $temparr1[2])) {
					delete $$sel_loc_r{$poplocID};#delete from selected loci
					delete $$sel_gt_r{$poplocID};#delete from selected genotypes
				}				
			} else {++$i} #if not count locus without available depth data
		}
		#print to report file
		print OUTREP "Selection based on depth analysis:\n";
		for $i (sort {$a <=> $b} keys %da_slr_var) {
			print OUTREP "$da_slr_var{$i}\n";
		}
		#Inform user about selected loci with no available depth data
		if ($i > 0) {
			print "Current selection included $i loci without available depth data.\n",
			"I treated these loci as having acceptable depths.\n";
			print OUTREP "No depth data were available for $i selected loci.\n",
			"These loci were treated as having acceptable depth.\n";
		}		
	}
	#}
	
	########################
	#merge: if activated, do
	########################
	
	#{
	if ($onoff{'merge'} eq 'on') {
		for $poplocID (keys %{$sel_loc_r}) {#loop through selected loci
			$merge_thisloc = $$poploc_r{$poplocID}[10];#look up merge_conflict in %poploc
			#delete locus if it doesn't match the selected criterion
			unless($merge_thisloc == $merge_sel) {
				delete $$sel_loc_r{$poplocID};#delete from selected loci
				delete $$sel_gt_r{$poplocID};#delete from selected genotypes
			}
		}
		#print to report file
		if ($merge_sel == 0) {
			print OUTREP "Loci without merging conflict.\n";
		} else {
			print OUTREP "Loci with merging conflict.\n";
		}
	}
	#}

	##############################################
	#Used fraction of read depth: if activated, do
	##############################################
	
	#{
	if ($onoff{'fracud'} eq 'on') {
		$depth = 0;
		#loop through loci in %sel_gt
		for $poplocID (keys %{$sel_gt_r}) {
			$fracud_del = 0;
			#loop through individuals (selected genotypes)
			for $ind (keys %{$$sel_gt_r{$poplocID}}) {
				#loop through indlocIDs for this genotype in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						#lookup Alldep and add to depth
						$depth += $indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[1];
					}
				}
				$totrd = $indpoplocdep{$poplocID}{$ind};#look up total read depth
				$fracud = $depth / $totrd;#calculated used fraction of read depth
				#if fraction too small, note
				if ($fracud < $fracud_min) {
					$fracud_del = 1;
				}
				$depth = 0;#set back and next ind
			}
			if ($fracud_del == 1) {#if fracud was too small in any ind
				delete $$sel_loc_r{$poplocID};#delete locus from selected loci
				delete $$sel_gt_r{$poplocID};#delete locus from selected genotypes
			}
		}
		#print to report file
		print OUTREP "Used fraction of read depth: min $fracud_min\n";
	}	
	#}	

	##############################################################
	#Delete empty individuals (individuals with no locus selected)
	##############################################################
	
	#{
	#set all individuals in %sel_ind to 0
	for $ind (keys %{$sel_ind_r}) {
		$$sel_ind_r{$ind} = 0;
	}
	#loop through selected individuals
	IND: for $ind (keys %{$sel_ind_r}) {
		#loop through loci in %sel_gt
		POPLOC: for $poplocID (keys %{$sel_gt_r}) {
			#loop through individuals in this locus in %sel_gt
			THISIND: for $thisind (keys %{$$sel_gt_r{$poplocID}}) {
				$$sel_ind_r{$thisind} = 1;#set this ind to 1 in %sel_ind
				last POPLOC if ($$sel_ind_r{$ind} == 1);#go to next ind in outer loop once you found ind
			}			
		}
		#delete the individual if it was not found
		if ($$sel_ind_r{$ind} == 0) {
			delete $$sel_ind_r{$ind};
		}
	}		
	#}
	
	#########################################
	#Determine dimensions of selected dataset
	#Append to file export/sel_rep.txt
	#########################################
	
	#{
	#n_ind
	$n_ind_dataset = keys %{$sel_ind_r};
	#n_loc
	$n_loc_dataset = keys %{$sel_loc_r};
	#n_gt
	$n_gt_dataset = 0;
	for $poplocID (keys %{$sel_gt_r}) {
		$n_gt_dataset += keys %{$$sel_gt_r{$poplocID}};
	}
	print OUTREP "\nSelected data include\n\n";
	print OUTREP "$n_ind_dataset individuals\n$n_loc_dataset loci\n$n_gt_dataset genotypes\n\n\n";
	close OUTREP;
	#}	
}

#definition of subroutine select_genotypes
#selects genotypes, updates %sel_ind, %sel_loc, %sel_gt, appends to report-file
#expects:
#ref to %sel_ind: {ind}=1
#ref to %sel_loc: {poplocID}=1
#ref to %sel_gt: {poplocID}{ind}=1
#$indpopallsuffix: suffix of indpopall-outfiles of indpoploc
#$alleles_suffix: suffix of alleles-outfile of indloc
#$loci_suffix = '_loci.txt';#suffix of *_loci outfiles out of indloc
sub select_genotypes {
	#######################
	#declare and initialize
	#######################
	
	#{
	#_r means a reference
	my ($sel_ind_r,$sel_loc_r,$sel_gt_r,$indpopallsuffix,$alleles_suffix,$loci_suffix) = @_;
	my $entry = '';#user entry
	my $usersel = 'invalid';#valid selection made by user
	my %number_option = (
	1 => 'depth',
	2 => 'totrd',
	3 => 'fracud',
	4 => 'depth_minor',
	5 => 'nAll',
	6 => 'gtlist',
	7 => 'run',
	8 => 'quit'
	);
	my %onoff = (
	depth => 'off',#depth of the genotype with min (sum of estimated depth of alleles)
	totrd => 'off',#total read depth of poploc in ind
	fracud => 'off',#fraction of used read depth: depth of genotype / total read depth
	depth_minor => 'off',#depth of the minor allele (homozygous: depth/2) with min
	nAll => 'off',#number of alleles with min max
	gtlist => 'off',#user provides a list of genotypes
	run => 'off',#run
	quit => 'off',#quit the program
	invalid => 'off'#dummy for invalid user selection
	);
	my $depth = 0;#depth of a genotype (sum of estimated depth of alleles)
	my $depth_min = 0;#min depth of genotype  (sum of estimated depth of alleles)
	my $totrd = 0;#total read depth of a poplocus in an individual
	my $totrd_min = 0;#min total read depth of a locus
	my $fracud = 0;#fraction of used read depth: depth of genotype / total read depth
	my $fracud_min = 0;#min fraction of used read depth
	my $depth_minor = 0;#depth of minor allele
	my $depth_minor_min = 0;#min depth of minor allele
	my @depths_thisgt = ();#depths of all alleles of a genotype
	my $nAll = 0;#number of alleles
	my $nAll_min = 0;#min number of alleles
	my $nAll_max = 0;#max number of alleles
	my $gtlistname = '';#filename: list of genotypes
	my %indpopall = ();# {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
						#Attention: selection procedures do not update this hash!
						#Always access it via loci in %sel_loc or %sel_gt (always up to date)
	my $ind = '';#individual ID
	my $thisind = '';#an individual ID
	my $infilename = '';#an infilename
	my $poplocID = 0;#population-level locus ID
	my $indlocID = 0;#individual locus ID
	my $popall_ID = 0;#population-level allele ID
	my $indall_ID = 0;#individual allele ID
	my $Alldep = 0;#depth of an allele
	my %Alldep_thisind = ();#for one ind: {indall_ID} = Alldep
	my %indlocdep_this_ind = ();# {indlocID} = locdep
	my %indpoplocdep = ();# {poplocID}{ind}=locdep
	my $n_gt = 0;#number of genotypes in a locus
	my %user_sel_gt = ();#according to gtlist provided by user: {poplocID}{ind}=1

	#variables about dataset-dimensions:
	my $n_ind_dataset = 0;#number of individuals
	my $n_loc_dataset = 0;#number of loci
	my $n_gt_dataset = 0;#number of genotypes
	
	my $tempstring1 = '';
	my @temparr1 = ();
	#}
	
	#####################
	#User selects options
	#####################
	
	#{
	while (($onoff{'run'} eq 'off') and ($onoff{'quit'} eq 'off')) {
		print "\nSELECT GENOTYPES:\n\n";
		print "1 Depth of genotype            $onoff{'depth'}";
		if ($onoff{'depth'} eq 'on') {
			print " min: $depth_min\n";
		} else {print "\n"}
		print "2 Total read depth of genotype $onoff{'totrd'}";
		if ($onoff{'totrd'} eq 'on') {
			print " min: $totrd_min\n";
		} else {print "\n"}
		print "3 Used fraction of read depth  $onoff{'fracud'}";
		if ($onoff{'fracud'} eq 'on') {
			print " min: $fracud_min\n";
		} else {print "\n"}
		print "4 Depth of minor allele        $onoff{'depth_minor'}";
		if ($onoff{'depth_minor'} eq 'on') {
			print " min: $depth_minor_min\n";
		} else {print "\n"}
		print "5 Number of alleles            $onoff{'nAll'}";
		if ($onoff{'nAll'} eq 'on') {
			print " min: $nAll_min max: $nAll_max\n";
		} else {print "\n"}
		print "6 List of genotypes in a file  $onoff{'gtlist'}";
		if ($onoff{'gtlist'} eq 'on') {
			print " $gtlistname\n";
		} else {print "\n"}
		print "\n7 Run\n";
		print "8 Quit\n\n";
		print "\nPlease select a number to activate/inactivate an option: ";
		$entry = <STDIN>;#get user-entry
		chomp $entry;
		if (defined $number_option{$entry}) {#if it is valid
			$usersel = $number_option{$entry};#take over as selection
		} else {$usersel = 'invalid'};#else set selection to 'invalid'
		
		#switch corresponding switch
		if ($onoff{$usersel} eq 'off') {
			$onoff{$usersel} = 'on';
		} else {$onoff{$usersel} = 'off'};
		
		#if user activated depth, get min
		if (($usersel eq 'depth') and ($onoff{'depth'} eq 'on')){
			print "\nDepth of genotype: please enter min: ";
			$depth_min = <STDIN>;
			chomp $depth_min;
			#Warn and turn off if value is invalid
			unless(($depth_min =~ /^\d+$/) and ($depth_min >= 0)) {
				print "Min must be 0 or positive integer.\n";
				$onoff{'depth'} = 'off';
			}
		}
		#if user activated total read depth, get min
		elsif (($usersel eq 'totrd') and ($onoff{'totrd'} eq 'on')){
			print "\nTotal read depth of genotype: please enter min: ";
			$totrd_min = <STDIN>;
			chomp $totrd_min;
			#Warn and turn off if value is invalid
			unless(($totrd_min =~ /^\d+$/) and ($totrd_min >= 0)) {
				print "Min must be 0 or positive integer.\n";
				$onoff{'totrd'} = 'off';
			}
		}
		#if user activated fraction of used read depth, get min
		elsif (($usersel eq 'fracud') and ($onoff{'fracud'} eq 'on')){
			print "\nUsed fraction of read depth:\n",
			"please enter min (decimal between 0 and 1): ";
			$fracud_min = <STDIN>;
			chomp $fracud_min;
			#Warn and turn off if value is invalid (must be between 0 and 1)
			unless((($fracud_min =~ /^\d$/) or ($fracud_min =~ /^\d\.\d+$/)) and
			($fracud_min >= 0) and ($fracud_min <= 1)) {
				print "Min must be >= 0 and <= 1.\n";
				$onoff{'fracud'} = 'off';
			}
		}
		#if user activated depth_minor, get min
		elsif (($usersel eq 'depth_minor') and ($onoff{'depth_minor'} eq 'on')){
			print "\nDepth of minor allele: please enter min: ";
			$depth_minor_min = <STDIN>;
			chomp $depth_minor_min;
			#Warn and turn off if value is invalid
			unless(($depth_minor_min =~ /^\d+$/) and ($depth_minor_min >= 0)) {
				print "Min must be 0 or positive integer.\n";
				$onoff{'depth_minor'} = 'off';
			}
		}
		#if user activated nAll, get min and max
		elsif (($usersel eq 'nAll') and ($onoff{'nAll'} eq 'on')){
			print "\nNumber of alleles: please enter min: ";
			$nAll_min = <STDIN>;
			chomp $nAll_min;
			print "Number of alleles: please enter max: ";
			$nAll_max = <STDIN>;
			chomp $nAll_max;
			
			#Warn and turn off if values are invalid
			unless(($nAll_min =~ /^\d+$/) and ($nAll_max =~ /^\d+$/) and ($nAll_max >= $nAll_min) and ($nAll_max > 0)) {
				print "Min must be 0 or positive integer.\n",
				"Max must be positive integer >= min.\n";
				$onoff{'nAll'} = 'off';
			}
		}
		#if user activated gtlist, get filename
		elsif (($usersel eq 'gtlist') and ($onoff{'gtlist'} eq 'on')){
			print
			"\nList of genotypes: I need a textfile:\n",
			"Each line: poplocID TAB ind_ID\n",
			"no headers, spaces, blank lines\n",
			"Please enter filename: ";
			$gtlistname = <STDIN>;
			chomp $gtlistname;
			#open the file or inactivate the option again
			unless(open(GTLIST, $gtlistname)) {
				print"Can't open $gtlistname.\n";
				$onoff{'gtlist'} = 'off';
			}
		}
		#if user inactivated gtlist, close file
		elsif (($usersel eq 'gtlist') and ($onoff{'gtlist'} eq 'off')){
			close GTLIST;
		}
		#if user selected quit, ask for confirmation and quit if appropriate
		elsif ($usersel eq 'quit') {
			print "Do you really want to quit? (Y/N): ";
			$entry = <STDIN>;
			chomp $entry;
			if ($entry eq 'Y') {
				print "No data selected. Exiting..\n";
				exit;
			} else {$onoff{'quit'} = 'off'};
		}		
	}
	#}

	###############################################
	#Read in *_indpopall.txt files out of indpoploc
	#into %indpopall, remove unselected data
	###############################################
	
	#{
	#loop through selected individuals
	for $ind (keys %{$sel_ind_r}) {
		$infilename = $ind.$indpopallsuffix;#build infilename
		unless(open(INFILE, $infilename)) {#open file or die
			print "Can't open $infilename. Exiting..\n";
			exit;
		}
		$tempstring1 = <INFILE>;#skip headerline
		while ($tempstring1 = <INFILE>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#load data into %indpopall {poplocID}{ind}{indlocID}{popall_ID}=(indall_ID,Alldep)
			#Infilecolumns!
			push @{$indpopall{$temparr1[2]}{$ind}{$temparr1[0]}{$temparr1[3]}}, $temparr1[1];
		}
		close INFILE;
	}
	#remove unselected data from %indpopall
	#loop through poplocIDs in %indpopall
	for $poplocID (keys %indpopall) {
		if (defined $$sel_loc_r{$poplocID}) {#if locus is selected
			#loop through individuals in %indpopall{$poplocID}
			for $ind (keys %{$indpopall{$poplocID}}) {
				#if this individual genotype (poplocID x ind) is not selected
				unless ((defined $$sel_gt_r{$poplocID}) and (defined $$sel_gt_r{$poplocID}{$ind})) {
					delete $indpopall{$poplocID}{$ind};
				}
			}
		} else {#if locus is not selected
			delete $indpopall{$poplocID};
		}
	}	
	#}

	##########################################
	#Read in *_alleles.txt files out of indloc
	#Load Alldep into %indpopall
	##########################################
	
	#{
	#loop through selected individuals
	for $ind (keys %{$sel_ind_r}) {
		$infilename = $ind.$alleles_suffix;#build infilename
		unless(open(INFILE, $infilename)) {#open file or die
			print "Can't open $infilename. Exiting..\n";
			exit;
		}
		$tempstring1 = <INFILE>;#skip headerline
		while ($tempstring1 = <INFILE>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#Load data into %Alldep_thisind: {indall_ID} = Alldep
			#Infilecolumns!
			if ($temparr1[1] eq 'valid') {
				$Alldep_thisind{$temparr1[3]} = $temparr1[5];
			}
		}
		close INFILE;
		#loop through poplocIDs in %indpoploc
		for $poplocID (keys %indpopall) {
			#loop through indlocIDs
			for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
				#loop through popall_IDs
				for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
					#look up indall_ID
					$indall_ID = $indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[0];
					#look up Alldep
					$Alldep = $Alldep_thisind{$indall_ID};
					#add to %indpopall
					$indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[1] = $Alldep;
				}
			}
		}
		#set back and next ind
		%Alldep_thisind = ();
	}
	#}
	
	#######################################
	#Read in *_loci.txt files out of indloc
	#populate %indpoplocdep
	#######################################
	
	#{	
	#loop through selected individuals
	for $ind (keys %{$sel_ind_r}) {
		$infilename = $ind.$loci_suffix;
		unless(open(INFILE, $infilename)) {#open file or die
			print "Can't open $infilename. Exiting..\n";
			exit;
		}
		$tempstring1 = <INFILE>;#skip headerline
		while ($tempstring1 = <INFILE>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#load data into %indlocdep_this_ind
			#Infilecolumns!
			if ($temparr1[1] eq 'valid') {
				$indlocdep_this_ind{$temparr1[0]} = $temparr1[3];
			}
		}
		close INFILE;
		#loop through selected poplocs and corresponding indlocs for this ind
		for $poplocID (keys %indpopall) {
			for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
				#load depth into %indpoplocdep
				$indpoplocdep{$poplocID}{$ind} += $indlocdep_this_ind{$indlocID};
			}
		}
		%indlocdep_this_ind = ();#set back and next ind
	}	
	#}	
	
	################################################################
	#Determine starting dimensions of database or previous selection
	#Append to file export/sel_rep.txt
	################################################################
	
	#{
	#n_ind
	$n_ind_dataset = keys %{$sel_ind_r};
	#n_loc
	$n_loc_dataset = keys %{$sel_loc_r};
	#n_gt
	for $poplocID (keys %{$sel_gt_r}) {
		$n_gt_dataset += keys %{$$sel_gt_r{$poplocID}};
	}
	#open export/sel_rep.txt for append
	unless(open(OUTREP, ">>export/sel_rep.txt")) {
		print "Cannot open file export/sel_rep.txt, exiting ...\n\n";
		exit;
	}
	print OUTREP
	"Selection of Genotypes: starting from\n\n",
	"$n_ind_dataset individuals\n$n_loc_dataset loci\n$n_gt_dataset genotypes\n\n",
	"Selected genotypes according to these criteria:\n\n";	
	#}
	
	########################
	#depth: if activated, do
	########################
	
	#{
	if ($onoff{'depth'} eq 'on') {
		$depth = 0;
		#loop through loci in %sel_gt
		for $poplocID (keys %{$sel_gt_r}) {
			#loop through individuals (selected genotypes)
			for $ind (keys %{$$sel_gt_r{$poplocID}}) {
				#loop through indlocIDs for this genotype in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						#lookup Alldep and add to depth
						$depth += $indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[1];
					}
				}
				#if the genotype is not deep enough, delete it
				unless ($depth >= $depth_min) {
					delete $$sel_gt_r{$poplocID}{$ind};#delete from selected genotypes
				}
				$depth = 0;#set back and next ind
			}
			#if the locus is empty now (no genotypes left)
			$n_gt = keys %{$$sel_gt_r{$poplocID}};
			if ($n_gt == 0) {
				delete $$sel_gt_r{$poplocID};#delete locus from selected genotypes
				delete $$sel_loc_r{$poplocID};#delete locus from selected loci
			}
		}
		#print to report file
		print OUTREP "Depth of genotype: min $depth_min\n";
	}
	#}
	
	###################################
	#Total read depth: if activated, do
	###################################
	
	#{
	if ($onoff{'totrd'} eq 'on') {
		#loop through loci in %sel_gt
		for $poplocID (keys %{$sel_gt_r}) {
			#loop through individuals (selected genotypes)
			for $ind (keys %{$$sel_gt_r{$poplocID}}) {
				$totrd = $indpoplocdep{$poplocID}{$ind};#get total read depth
				#if total read depth too small: delete genotype
				unless ($totrd >= $totrd_min) {
					delete $$sel_gt_r{$poplocID}{$ind};#delete from selected genotypes
				}
			}
			#if the locus is empty now (no genotypes left)
			$n_gt = keys %{$$sel_gt_r{$poplocID}};
			if ($n_gt == 0) {
				delete $$sel_gt_r{$poplocID};#delete locus from selected genotypes
				delete $$sel_loc_r{$poplocID};#delete locus from selected loci
			}
		}
		#print to report file
		print OUTREP "Total read depth of genotype: min $totrd_min\n";
	}	
	#}
	
	##############################################
	#Used fraction of read depth: if activated, do
	##############################################
	
	#{
	if ($onoff{'fracud'} eq 'on') {
		$depth = 0;
		#loop through loci in %sel_gt
		for $poplocID (keys %{$sel_gt_r}) {
			#loop through individuals (selected genotypes)
			for $ind (keys %{$$sel_gt_r{$poplocID}}) {
				#loop through indlocIDs for this genotype in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						#lookup Alldep and add to depth
						$depth += $indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[1];
					}
				}
				$totrd = $indpoplocdep{$poplocID}{$ind};#look up total read depth
				$fracud = $depth / $totrd;#calculated used fraction of read depth
				#if fraction too small, delete genotype
				unless ($fracud >= $fracud_min) {
					delete $$sel_gt_r{$poplocID}{$ind};#delete from selected genotypes
				}
				$depth = 0;#set back and next ind
			}
			#if the locus is empty now (no genotypes left)
			$n_gt = keys %{$$sel_gt_r{$poplocID}};
			if ($n_gt == 0) {
				delete $$sel_gt_r{$poplocID};#delete locus from selected genotypes
				delete $$sel_loc_r{$poplocID};#delete locus from selected loci
			}
		}
		#print to report file
		print OUTREP "Used fraction of read depth: min $fracud_min\n";
	}	
	#}	
	
	##############################
	#depth_minor: if activated, do
	##############################
	
	#{
	if ($onoff{'depth_minor'} eq 'on') {
		$depth_minor = 0;
		#loop through loci in %sel_gt
		for $poplocID (keys %{$sel_gt_r}) {
			#loop through individuals (selected genotypes)
			for $ind (keys %{$$sel_gt_r{$poplocID}}) {
				#loop through indlocIDs for this genotype in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						#lookup Alldep and add to @depths_thisgt
						push @depths_thisgt, $indpopall{$poplocID}{$ind}{$indlocID}{$popall_ID}[1];
					}
				}
				#determine depth of minor allele
				if (@depths_thisgt > 1) {#if there is more than one allele
					@depths_thisgt = sort {$a <=> $b} @depths_thisgt;
					$depth_minor = $depths_thisgt[0];
				}
				elsif (@depths_thisgt == 1) {#if there is one allele
					$depth_minor = ($depths_thisgt[0] / 2);
				}
				else {
					print "Error: no allele depths found in locus $poplocID ind $ind\n";
				}
				
				#if the minor allele is not deep enough, delete the genotype
				unless ($depth_minor >= $depth_minor_min) {
					delete $$sel_gt_r{$poplocID}{$ind};#delete from selected genotypes
				}
				#set variables back and next ind
				$depth_minor = 0;
				@depths_thisgt = ();
			}
			#if the locus is empty now (no genotypes left)
			$n_gt = keys %{$$sel_gt_r{$poplocID}};
			if ($n_gt == 0) {
				delete $$sel_gt_r{$poplocID};#delete locus from selected genotypes
				delete $$sel_loc_r{$poplocID};#delete locus from selected loci
			}
		}
		#print to report file
		print OUTREP "Depth of minor allele: min $depth_minor_min\n";
	}
	#}
	
	#######################
	#nAll: if activated, do
	#######################
	
	#{
	if ($onoff{'nAll'} eq 'on') {
		$nAll = 0;
		#loop through loci in %sel_gt
		for $poplocID (keys %{$sel_gt_r}) {
			#loop through individuals (selected genotypes)
			for $ind (keys %{$$sel_gt_r{$poplocID}}) {
				#loop through indlocIDs for this genotype in %indpopall
				for $indlocID (keys %{$indpopall{$poplocID}{$ind}}) {
					#loop through popall_IDs
					for $popall_ID (keys %{$indpopall{$poplocID}{$ind}{$indlocID}}) {
						++$nAll;#count the allele
					}
				}
				#if the number of alleles does not match the criteria, delete genotype
				unless (($nAll >= $nAll_min) and ($nAll <= $nAll_max)) {
					delete $$sel_gt_r{$poplocID}{$ind};#delete from selected genotypes
				}
				$nAll = 0;#set back and next ind
			}
			#if the locus is empty now (no genotypes left)
			$n_gt = keys %{$$sel_gt_r{$poplocID}};
			if ($n_gt == 0) {
				delete $$sel_gt_r{$poplocID};#delete locus from selected genotypes
				delete $$sel_loc_r{$poplocID};#delete locus from selected loci
			}
		}
		#print to report file
		print OUTREP "Number of alleles: min $nAll_min max: $nAll_max\n";
	}
	#}
	
	#########################
	#gtlist: if activated: do
	#########################
	
	#{
	if ($onoff{'gtlist'} eq 'on') {
		#Read in gtlist file: poplocID TAB ind_ID
		#and populate %user_sel_gt: {poplocID}{ind}=1
		while ($tempstring1 = <GTLIST>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			$user_sel_gt{$temparr1[0]}{$temparr1[1]} = 1;
		}
		close GTLIST;
		#Loop through selected loci in %sel_gt: {poplocID}{ind}=1
		for $poplocID (keys %{$sel_gt_r}) {
			#loop through selected individuals (genotypes) for this locus
			for $ind (keys %{$$sel_gt_r{$poplocID}}) {
				#if the genotype is not in user-list (%user_sel_gt): delete
				unless((defined $user_sel_gt{$poplocID}) and (defined $user_sel_gt{$poplocID}{$ind})) {
					delete $$sel_gt_r{$poplocID}{$ind};#delete from selected genotypes
				}
			}
			#if the locus is empty now (no genotypes left)
			$n_gt = keys %{$$sel_gt_r{$poplocID}};
			if ($n_gt == 0) {
				delete $$sel_gt_r{$poplocID};#delete locus from selected genotypes
				delete $$sel_loc_r{$poplocID};#delete locus from selected loci
			}
		}
		#print to report file
		print OUTREP "List of genotypes in file $gtlistname\n";
	}
	#}
	
	##############################################################
	#Delete empty individuals (individuals with no locus selected)
	##############################################################
	
	#{
	#set all individuals in %sel_ind to 0
	for $ind (keys %{$sel_ind_r}) {
		$$sel_ind_r{$ind} = 0;
	}
	#loop through selected individuals
	IND: for $ind (keys %{$sel_ind_r}) {
		#loop through loci in %sel_gt
		POPLOC: for $poplocID (keys %{$sel_gt_r}) {
			#loop through individuals in this locus in %sel_gt
			THISIND: for $thisind (keys %{$$sel_gt_r{$poplocID}}) {
				$$sel_ind_r{$thisind} = 1;#set this ind to 1 in %sel_ind
				last POPLOC if ($$sel_ind_r{$ind} == 1);#go to next ind in outer loop once you found ind
			}			
		}
		#delete the individual if it was not found
		if ($$sel_ind_r{$ind} == 0) {
			delete $$sel_ind_r{$ind};
		}
	}		
	#}

	#########################################
	#Determine dimensions of selected dataset
	#Append to file export/sel_rep.txt
	#########################################
	
	#{
	#n_ind
	$n_ind_dataset = keys %{$sel_ind_r};
	#n_loc
	$n_loc_dataset = keys %{$sel_loc_r};
	#n_gt
	$n_gt_dataset = 0;
	for $poplocID (keys %{$sel_gt_r}) {
		$n_gt_dataset += keys %{$$sel_gt_r{$poplocID}};
	}
	print OUTREP "\nSelected data include\n\n";
	print OUTREP "$n_ind_dataset individuals\n$n_loc_dataset loci\n$n_gt_dataset genotypes\n\n\n";
	close OUTREP;
	#}
}

#definition of subroutine varpos
#expects: ref to sequence alignment as 2d array:
#d1: rows, d2: cols
#returns: an array of variable positions
#if no positions are variable: empty
#allowed characters: ACGTN, only capital, N is a 5th char
#all sequences must have same length
sub varpos {
	#declare and initialize
	my ($seqmat_r) = @_;
	my @seqmat = @{$seqmat_r};
	my @varpos = ();#the variable positions
	my $row = 0;
	my $col = 0;
	my $char1 = '';#a character in the first sequence
	my %count = ();#a hash for counting characters
	
	#loop through cols
	for ($col = 0; $col < @{$seqmat[0]}; ++$col) {
		$char1 = ${$seqmat[0]}[$col];#character of 1st sequence in this col
		#count how many sequences have this character:
		#loop through rows
		for ($row = 0; $row < @seqmat; ++$row) {
			++$count{${$seqmat[$row]}[$col]};
		}
		#if not all have this character:
		if ($count{$char1} < @seqmat) {
			#add this col to @varpos
			push @varpos, $col;
		}
		%count = ();#set counthash back
	}
	return @varpos;
}
