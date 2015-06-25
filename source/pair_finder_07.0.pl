#!/usr/bin/perl -w
#pair_finder version 07.0 Copyright 2015 Andreas Hapke
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
use IO::File; #requires Perl 5.004 or higher
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

#uses subroutine:
#clust_loc

#Keyword Infilecolumns! Order of columns in an infile is important.

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "pair_finder version 07.0 Copyright 2015 Andreas Hapke\n",
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
my $individuals_fname = 'individuals.txt';#list of all individual IDs in database out of indloc
my $poploc_name = 'poploc.txt';#out of poploc: poplocID sl cons nSNP varpos n_all nInd nIndloc ties
my $sel_ind_name = 'export/sel_ind.txt';#out of data_selector: ind
my $sel_loc_name = 'export/sel_loc.txt';#out of data_selector: poplocID
my $sel_gt_name = 'export/sel_gt.txt';#out of data_selector: poplocID ind
my $sel_rep_name = 'export/sel_rep.txt';#out of data_selector: report file
my $splitloc_name = 'split_loci.txt';#out of indpoploc: poplocID ind indlocID
my $indpoploc_suffix = '_indpoploc.txt';#suffix of indpoploc outfile of indpoploc
my $reads_suffix = '_reads.txt';#suffix of reads outfile of indloc
my $fq_header_del = '_';#field-separator in fastq-headers
my $fq_header_fr_field = 4;#field in fastq-header that distinguishes forward and reverse reads
							#count starting with 0
my $outdirname = 'pairs';#directory for outfiles
#other variables
my $ind = '';#an individual ID
my $thisind = '';#current individual ID
my %sel_ind = ();# {ind}=1
my %sel_loc = ();# {poplocID}=1
my %sel_gt = ();# {poplocID}{ind}=1
my $poplocID = 0;#population level locus ID
my $poplocID1 = 0;#poplocID 1 in a pair
my $poplocID2 = 0;#poplocID 2 in a pair
my %splitloc = ();# {poplocID}{ind} = indlocID
my $n_splitloc = 0;#number of split loci
my %indpoploc = ();# {indlocID}=poplocID
my %head_poploc = ();# {fqheader} = poplocIDs
					#holds the first fqheader of a pair of headers and its poplocID
					#each new header is looked up in here, if defined, the pair of headers is found
					#without storing both
my $head = '';#fqheader without field that distinguishes forward and reverse read
my %pair_counts = ();# {poplocID1}{poplocID2} = n_seqpairs
					#poplocIDs are sorted
					#n_seqpairs: number of forward-reverse seqpairs connecting
					#these poplocs in an individual
my %pair_counts_allind = ();# {poplocID1}{poplocID2}=(n_seqpairs,n_ind)
my %self_match_loc = ();# {poplocID} = n_seqpairs
					#self-matching loci (forward matches revcomp of reverse) in an individual
					#These loci have invalid genotypes and must be excluded
my %self_match_loc_allind = ();# {poplocID}=(n_seqpairs,n_ind)
my %clust_loc = ();# {groupNo}{poplocID}=1
my $groupNo = 0;#Number for a group of loc
my $outgroupNo = 0;#new consecutive number for a group of loci
my $n_loc = 0;#number of loci in a group
my $n_self_loc = 0;#total of self-matching loci
my %group_sumcount = ();# {$n_loc} = number of groups

my $infilename = '';
my $outfilename = '';
my $reads_h = '';#handle for reads infile (outfile of indloc)

my @temparr1 = ();
my @temparr2 = ();
my $tempstring1 = '';
my $tempstring2 = '';
my $i = 0;
my $message_counter = 0;
my $message_total = 0;
#}

###########################################################
#If directory outdirname already exists: exit, else: create
###########################################################

#{
if (-d "$outdirname") {
	print "Directory $outdirname already exists.\n",
	"Please delete or rename. Exiting..\n";
	exit;
} else {#if not, try to create
	unless (mkdir "$outdirname") {
		print "Can't create directory $outdirname. Exiting..\n";
		exit;
	}
}
#}

##########################
#Create file pairs_rep.txt
##########################

#{
unless(open(OUTREP, ">$outdirname/pairs_rep.txt")) {#open or die
	print "Cannot open $outdirname/pairs_rep.txt, exiting ...\n\n";
	exit;
}
#}

###############################
#Read in existing preselections
###############################

#{
#if directory "export" exists, check for selection files and read them in
if (-d "export") {
	if (-f "$sel_ind_name") {#if there is an individual selection-file
		print "Loading your selection of individuals.\n";
		unless(open(INFILE,"$sel_ind_name")) {
			print "Can't open $sel_ind_name. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;
		}
		@temparr1 = <INFILE>;
		close INFILE;
		#Infilecolumns!
		for $ind (@temparr1) {
			chomp $ind;
			$sel_ind{$ind} = 1;#store selected inds in %sel_ind
		}
		#if the selection file was empty (%sel_ind still empty now)
		if (keys %sel_ind == 0) {
			print "Your selection evaluates to an empty dataset.\n",
			"I can't continue. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;			
		}
	}
	if (-f "$sel_loc_name") {#if there is a locus selection-file
		print "Loading your selection of loci.\n";
		unless(open(INFILE,"$sel_loc_name")) {
			print "Cannot open $sel_loc_name. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;
		}
		#Infilecolumns!
		while ($poplocID = <INFILE>) {#read in poplocIDs
			chomp $poplocID;
			$sel_loc{$poplocID} = 1;#store in %sel_loc
		}
		close INFILE;
		#if the selection file was empty (%sel_loc still empty now)
		if (keys %sel_loc == 0) {
			print "Your selection evaluates to an empty dataset.\n",
			"I can't continue. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;			
		}
	}
	if (-f "$sel_gt_name") {#if there is a genotype selection-file
		print "Loading your selection of genotypes.\n";
		unless(open(INFILE,"$sel_gt_name")) {
			print "Cannot open $sel_gt_name. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;
		}
		#Infilecolumns!
		while ($tempstring1 = <INFILE>) {#read in
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			$sel_gt{$temparr1[0]}{$temparr1[1]} = 1;#store in %sel_gt
		}
		close INFILE;
		#if the selection file was empty (%sel_gt still empty now)
		if (keys %sel_gt == 0) {
			print "Your selection evaluates to an empty dataset.\n",
			"I can't continue. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;			
		}
	}
	#If no preselection has been loaded now: print information to report file
	if ((keys %sel_ind == 0) and (keys %sel_loc == 0) and (keys %sel_gt == 0)) {
		print OUTREP "PAIR_FINDER: analysis based on whole database.\n\n";
	}
	#If a preselection has been loaded now: print information about preselection to report file
	else {
		print OUTREP "PAIR_FINDER: analysis based on preselection of data.\n\n";
		#open the selection report and copy into pair-finder report
		unless(open(INFILE,"$sel_rep_name")) {
			print "Cannot open $sel_rep_name. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;
		}
		print OUTREP "Your preselection of data according to file $sel_rep_name:\n\n";
		while ($tempstring1 = <INFILE>) {
			chomp $tempstring1;
			print OUTREP "$tempstring1\n";
		}
		close INFILE;
		print OUTREP "\nPAIR_FINDER:\n\n";
	}
}
else {#if no directory "export" exists
	print OUTREP "PAIR_FINDER: analysis based on whole database.\n\n";
}
#}

##############################################################
#If no preselections exist, create selection of whole database
##############################################################

#{
#if %sel_ind is still empty, create it from file individuals.txt (using all inds)
if (keys %sel_ind == 0) {
	print "Found no selection of individuals, select all in the database.\n";
	unless(open(INFILE,"$individuals_fname")) {
		print "Cannot open $individuals_fname. Exiting..\n";
		print OUTREP "Analysis aborted.\n";
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
#if %sel_loc is still empty, read in poploc.txt outfile of poploc
#populate %sel_loc
if (keys %sel_loc == 0) {
	print "Found no selection of loci, select all in the database.\n";
	unless(open(INFILE,$poploc_name)) {
		print "Cannot open $poploc_name. Exiting..\n";
		print OUTREP "Analysis aborted.\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip headerline
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		#Infilecolumns!
		@temparr1 = split(/\t/,$tempstring1);
		$poplocID = shift @temparr1;#get poplocID
		$sel_loc{$poplocID} = 1;#store in %sel_loc
	}
	close INFILE;
}
#if %sel_gt is still empty, create from %sel_loc and %sel_ind
#(using all genotypes from selected inds and locs)
if (keys %sel_gt == 0) {
	print "Found no selection of genotypes, select all in the database.\n";
	for $poplocID (keys %sel_loc) {
		for $ind (keys %sel_ind) {
			$sel_gt{$poplocID}{$ind} = 1;
		}
	}
}
#}

#####################################################
#Read in split_loci.txt and check if any are selected
#If yes: inform user that they will be ignored
#####################################################

#{
unless(open(INFILE,$splitloc_name)) {
	print "Cannot open $splitloc_name. Exiting..\n";
	print OUTREP "Analysis aborted.\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	$splitloc{$temparr1[0]}{$temparr1[1]} = $temparr1[2];
}
#Check how many split loci are in the current selection
for $poplocID (keys %splitloc) {
	if (defined $sel_loc{$poplocID}) {
		++$n_splitloc;
	}	
}
#If you found any, exclude them from analysis, inform user
if ($n_splitloc > 0) {
	print "The current selection includes $n_splitloc split loci.\n",
	"I will exclude these loci from locus-pair analysis.\n";
	print OUTREP "Excluded $n_splitloc split loci.\n";
	for $poplocID (keys %splitloc) {
		delete $sel_loc{$poplocID};#delete from selected loci
		delete $sel_gt{$poplocID};#delete from selected genotypes
	}
	#check for empty individuals and exclude
	#set all individuals in %sel_ind to 0
	for $ind (keys %sel_ind) {
		$sel_ind{$ind} = 0;
	}
	#loop through selected individuals
	IND: for $ind (keys %sel_ind) {
		#loop through loci in %sel_gt
		POPLOC: for $poplocID (keys %sel_gt) {
			#loop through individuals in this locus in %sel_gt
			THISIND: for $thisind (keys %{$sel_gt{$poplocID}}) {
				$sel_ind{$thisind} = 1;#set this ind to 1 in %sel_ind
				last POPLOC if ($sel_ind{$ind} == 1);#go to next ind in outer loop once you found ind
			}			
		}
		#delete the individual if it was not found
		if ($sel_ind{$ind} == 0) {
			delete $sel_ind{$ind};
			print "No data left for individual $ind. I exclude it.\n";
			print OUTREP "No data left for individual $ind: excluded.\n";
		}
	}	
}
#}

#########################################################
#Loop through selected individuals, load data
#populate for this ind:
#%indpoploc %head_poploc %pair_counts %self_match_loc
#for all inds: %pair_counts_allind %self_match_loc_allind
#The latter 2 will later be modified
#########################################################

#{
#determine number of individuals for a message
$message_total = keys %sel_ind;
print "Analyzing data from $message_total individuals:\n";
for $ind (sort keys %sel_ind) {
	#read in *_indpoploc.txt: out of indpoploc: indlocID poplocID
	#populate %indpoploc: {indlocID}=poplocID
	$infilename = $ind.$indpoploc_suffix;
	unless(open(INFILE,"$infilename")) {
		print "Cannot open $infilename. Exiting..\n";
		print OUTREP "Analysis aborted.\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip headerline
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		$indpoploc{$temparr1[0]} = $temparr1[1];#store in %indpoploc
	}
	close INFILE;
	#read in *_reads.txt: out of indloc:
	#indlocID Loc_cat(valid/lost) seql All_ID hapID hapcat(used/discarded) readID
	#populate %head_poploc: {fqheader} = (poplocIDs)
	$infilename = $ind.$reads_suffix;
	if (-f $infilename) {#found file with extension .txt: open as plain text
		$reads_h = IO::File->new("< $infilename");
	} else {
		$infilename .= '.gz';
		if (-f $infilename) {#found file with extension .txt.gz: open as gzip compressed file
			$reads_h = new IO::Uncompress::Gunzip $infilename;
		} else {
			print "Individual $ind: I do not find file *_reads.txt or *_reads.txt.gz. Exiting..\n";
			print OUTREP "Analysis aborted.\n";
			exit;
		}
	}
	unless(defined $reads_h) {
		print "Individual $ind: I cannot open file *_reads.txt or *_reads.txt.gz. Exiting..\n";
		print OUTREP "Analysis aborted.\n";
		exit;
	}
	$tempstring1 = <$reads_h>;#skip headerline
	while ($tempstring1 = <$reads_h>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		if ($temparr1[1] eq 'valid') {#if locus is valid (not lost), i.e. is linked to a poploc at all
			#Infilecolumns!
			$tempstring2 = $temparr1[6];#get fastq-header
			@temparr2 = split($fq_header_del,$tempstring2);#split into fields in array
			#remove the field that distinguishes forward- and reverse reads
			splice (@temparr2,$fq_header_fr_field,1);
			$head = join('_',@temparr2);#join header again
			#take indlocID to look up poplocID
			#Infilecolumns!
			$poplocID1 = $indpoploc{$temparr1[0]};
			#if this poplocID for this ind is in current selection
			if ((exists $sel_gt{$poplocID1}) and (exists $sel_gt{$poplocID1}{$ind})) {
				#if this header (2nd of pair) is already stored in %head_poploc
				if (defined $head_poploc{$head}) {
					$poplocID2 = $head_poploc{$head};#get corresponding poplocID
					if ($poplocID1 ne $poplocID2) {#if the two poplocIDs are different
						@temparr1 = ($poplocID1,$poplocID2);#store them in array
						@temparr1 = sort {$a <=> $b} @temparr1;#sort them
						++$pair_counts{$temparr1[0]}{$temparr1[1]};#store pair in %pair_counts/increment n_seqpairs
					} else {#if the two poplocIDs are equal (locus is self-matching)
						++$self_match_loc{$poplocID1};#store in %self_match_loc
					}
					delete $head_poploc{$head};#this header is done, we don't need it any more
				} else {#if the header is new
					$head_poploc{$head} = $poplocID1;#store header and poplocID in %head_poploc
				}
			}
		}
	}
	close $reads_h;
	undef $reads_h;
	#Loop through %pair counts: {poplocID1}{poplocID2}=n_seqpairs
	#populate %pair_counts_allind: {poplocID1}{poplocID2}=(n_seqpairs,n_ind)
	for $poplocID1 (keys %pair_counts) {
		for $poplocID2 (keys %{$pair_counts{$poplocID1}}) {
			$pair_counts_allind{$poplocID1}{$poplocID2}[0] += $pair_counts{$poplocID1}{$poplocID2};
			$pair_counts_allind{$poplocID1}{$poplocID2}[1] += 1;
		}
	}
	#Loop through %self_match_loc: {poplocID} = n_seqpairs
	#populate %self_match_loc_allind: {poplocID}=(n_seqpairs,n_ind)
	for $poplocID (keys %self_match_loc) {
		$self_match_loc_allind{$poplocID}[0] += $self_match_loc{$poplocID};
		$self_match_loc_allind{$poplocID}[1] += 1;
	}
	#set variables back
	%indpoploc = ();
	%head_poploc = ();
	%pair_counts = ();
	%self_match_loc = ();
	#progress message
	++$message_counter;
	print "completed ind $ind: $message_counter of $message_total individuals\n";
	#next ind
}
print "Analyzing collected data..\n";
#}

###################################################################
#Loop through %self_match_loc_allind: {poplocID}=(n_seqpairs,n_ind)
#remove self-matching loci from %sel_loc and %sel_gt
###################################################################

#{
for $poplocID (keys %self_match_loc_allind) {
	delete $sel_loc{$poplocID};
	delete $sel_gt{$poplocID};
}
#}

###################################################################
#Call sub clust_loc to populate %clust_loc: {loc_group}{poplocID}=1
###################################################################

#{
clust_loc(\%pair_counts_allind,\%sel_loc,\%clust_loc);
#}

#############################################################
#Check if any locus group contains a self-matching locus
#add all loci that are connected with any self-matching locus
#to self-matching loci and remove their group
#############################################################

#{
#loop through locus groups
for $groupNo (keys %clust_loc) {
	$i = 0;#counter for self-matching loci in this group
	#loop through loci in group
	for $poplocID (keys %{$clust_loc{$groupNo}}) {
		if (defined $self_match_loc_allind{$poplocID}) {
			++$i;#count self-matching locus
		}
	}
	#if you found a self-matching locus in this group
	if ($i > 0) {
		#loop through loci in group
		for $poplocID (keys %{$clust_loc{$groupNo}}) {
			#add each locus to self-matching loci without counting n_seqpairs or n_ind
			$self_match_loc_allind{$poplocID}[0] += 0;
			$self_match_loc_allind{$poplocID}[1] += 0;
		}
		#delete group from %clust_loc
		delete $clust_loc{$groupNo};
	}
}
#determine number of self-matching loci and print to report
$n_self_loc = keys %self_match_loc_allind;
print OUTREP "Excluded $n_self_loc self-matching loci.\n";
print "Excluded $n_self_loc self-matching loci.\n";
#}

############################################################################
#Remove loci that are now in the self-matching loci from %pair_counts_allind
#sufficient to check locus1 in each pair (we got them by group,see above)
############################################################################

#{
for $poplocID1 (keys %pair_counts_allind) {
	if (defined $self_match_loc_allind{$poplocID1}) {
		delete $pair_counts_allind{$poplocID1};
	}
}
#}

#################################################
#Produce outfiles about data from all individuals
#Determine group summary counts
#################################################

#{
#self_match.txt: poplocID n_seqpairs n_ind
unless(open(OUTFILE, ">$outdirname/self_match.txt")) {#open or die
	print "Cannot open $outdirname/self_match.txt, exiting ...\n\n";
	print OUTREP "Analysis aborted.\n";
	exit;
}
#print headerline
print OUTFILE "poplocID\tn_seqpairs\tn_ind\n";
#loop through %self_match_loc_allind: {poplocID}=(n_seqpairs,n_ind)
for $poplocID (sort {$a <=> $b} keys %self_match_loc_allind) {
	$tempstring1 = join("\t",@{$self_match_loc_allind{$poplocID}});
	print OUTFILE "$poplocID\t$tempstring1\n";
}
close OUTFILE;

#pairs.txt: poplocID1 poplocID2 n_seqpairs n_ind
unless(open(OUTFILE, ">$outdirname/pairs.txt")) {#open or die
	print "Cannot open $outdirname/pairs.txt, exiting ...\n\n";
	print OUTREP "Analysis aborted.\n";
	exit;
}
#print headerline
print OUTFILE "poplocID1\tpoplocID2\tn_seqpairs\tn_ind\n";
#loop through %pair_counts_allind: {poplocID1}{poplocID2}=(n_seqpairs,n_ind)
for $poplocID1 (sort {$a <=> $b} keys %pair_counts_allind) {
	for $poplocID2 (sort {$a <=> $b} keys %{$pair_counts_allind{$poplocID1}}) {
		$tempstring1 = join("\t",@{$pair_counts_allind{$poplocID1}{$poplocID2}});
		print OUTFILE "$poplocID1\t$poplocID2\t$tempstring1\n";
	}
}
close OUTFILE;

#loc_groups.txt: loc_group n_loc poplocID
unless(open(OUTFILE, ">$outdirname/loc_groups.txt")) {#open or die
	print "Cannot open $outdirname/loc_groups.txt, exiting ...\n\n";
	print OUTREP "Analysis aborted.\n";
	exit;
}
#print headerline
print OUTFILE "loc_group\tn_loc\tpoplocID\n";
#loop through %clust_loc: {groupNo}{poplocID}=1
for $groupNo (sort {$a <=> $b} keys %clust_loc) {
	++$outgroupNo;#new consecutive group number for output
	$n_loc = keys %{$clust_loc{$groupNo}};#determine number of loci in group
	#count in %group_sumcount
	++$group_sumcount{$n_loc};
	for $poplocID (sort {$a <=> $b} keys %{$clust_loc{$groupNo}}) {
		print OUTFILE "$outgroupNo\t$n_loc\t$poplocID\n";
	}
}
close OUTFILE;
#}

#####################################
#Print group summary counts to report
#####################################

#{
print OUTREP "\nClustered loci into\n\n";
print "\nClustered loci into\n\n";
for $n_loc (sort {$a <=> $b} keys %group_sumcount) {
	print OUTREP "$group_sumcount{$n_loc}\tgroups with $n_loc loci\n";
	print "$group_sumcount{$n_loc}\tgroups with $n_loc loci\n";
}
print OUTREP "\n";
print "\n";
#}

print OUTREP "Analysis succesfully completed.\n";
close OUTREP;

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";

exit;

############
#Subroutines
############

#definition of subroutine clust_loc
#clusters loci into groups based on connecting read-pairs
#expects
#ref to %pair_counts_allind: {poplocID1}{poplocID2}=(n_seqpairs,n_ind)
		#needs only the 2 keys: pairs of loci
#ref to %sel_loc: {poplocID}=1 all good loci, including loci that are not in a pair
#ref to %clust_loc: {groupNo}{poplocID}=1
#populates %clust_loc
sub clust_loc {
	#declare and initialize: _r means a reference
	my ($pairs_r,$sel_loc_r,$clust_loc_r) = @_;
	my %loc_group = ();# {poplocID} = groupNo
	my %group_locs = ();# {groupNo} = (poplocIDs)
	my $poplocID = 0;#a locus ID
	my $poplocID1 = 0;# first locus ID in a pair
	my $poplocID2 = 0;# second locus ID in a pair
	my $groupNo = 0;#
	my $loc1group = 0;#groupNo of locus1 in a pair
	my $loc2group = 0;#groupNo of locus2 in a pair
	my @temparr1 = ();
	
	#loop through %sel_loc put all poplocIDs into %loc_group
	#initialize groupNo with 0
	for $poplocID (keys %{$sel_loc_r}) {
		$loc_group{$poplocID} = 0;
	}
	
	#####################
	#build groups of loci
	#####################
	
	#{
	#start with loci in pairs: loop through locus pairs in %pairs
	for $poplocID1 (keys %{$pairs_r}) {
		for $poplocID2 (keys %{$$pairs_r{$poplocID1}}) {
			#look up current groupNos
			$loc1group = $loc_group{$poplocID1};
			$loc2group = $loc_group{$poplocID2};
			#if both loci are not yet in a group
			if (($loc1group == 0) and ($loc2group == 0)) {
				#make up a new group for them
				++$groupNo;
				$loc_group{$poplocID1} = $groupNo;
				$loc_group{$poplocID2} = $groupNo;
				@{$group_locs{$groupNo}} = ($poplocID1,$poplocID2);
			}
			#if locus 1 is in a group but locus 2 is not
			elsif (($loc1group != 0) and ($loc2group == 0)) {
				#add locus 2 to group of locus 1
				$loc_group{$poplocID2} = $loc_group{$poplocID1};
				push @{$group_locs{$loc1group}}, $poplocID2;
			}
			#if locus 1 is not in a group but locus 2 is
			elsif (($loc1group == 0) and ($loc2group != 0)) {
				#add locus 1 to group of locus 2
				$loc_group{$poplocID1} = $loc_group{$poplocID2};
				push @{$group_locs{$loc2group}}, $poplocID1;
			}
			#if both are in different groups
			elsif ($loc1group != $loc2group) {
				#join the 2 groups
				#get all members of group loc2 is in
				@temparr1 = @{$group_locs{$loc2group}};
				#put them into group of loc1
				for $poplocID (@temparr1) {
					$loc_group{$poplocID} = $loc1group;
				}
				push @{$group_locs{$loc1group}}, @temparr1;
				#delete former group loc2 was in
				delete $group_locs{$loc2group};
			}			
		}
	}
	#continue with loci that are not yet in a group because they were not in a pair
	for $poplocID (keys %loc_group) {
		if ($loc_group{$poplocID} == 0) {#if the locus is not yet in a group
			#make up a single member group for it
			++$groupNo;
			$loc_group{$poplocID} = $groupNo;
			$group_locs{$groupNo}[0] = $poplocID;
		}
	}	
	#}
	
	############################
	#Populate %clust_loc in main
	############################
	
	#{
	for $groupNo (keys %group_locs) {
		for $poplocID (@{$group_locs{$groupNo}}) {
			$$clust_loc_r{$groupNo}{$poplocID} = 1;
		}
	}	
	#}
}
