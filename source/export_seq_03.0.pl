#!/usr/bin/perl -w
#export_seq version 03.0 Copyright 2015 Andreas Hapke
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
use List::Util qw(shuffle);

#Uses subroutines
#fasta_out
#nexus_out
#phylip_out

#Keyword Infilecolumns! Order of columns in an infile is important.

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "export_seq version 03.0 Copyright 2015 Andreas Hapke\n",
"This program is part of GIbPSs Copyright 2015 Andreas Hapke\n",
"GIbPSs comes with ABSOLUTELY NO WARRANTY. GIbPSs is free software,\n",
"and you are welcome to redistribute it under certain conditions.\n",
"For details see the GNU General Public License version 3.\n",
"You should have obtained a copy of the GNU General Public License\n",
"along with GIbPSs. If not, see <http://www.gnu.org/licenses/>.\n\n";

#######################
#Declare and initialize
#######################

#{
#user settings
my %user_settings = (
f => '',#format: f: FASTA, n: Nexus, p: PHYLIP
mi => '',#1: include missing individuals, 0: don't
m => '',#missing data character
indf => '',#name of user file to add information to indIDs
rs => ''#random number seed
);
my %defaults = (
f => 'f',#format: f: FASTA, n: Nexus, p: PHYLIP
mi => '0',#1: include missing individuals, 0: don't
m => 'N',#missing data character
indf => '',#name of user file to add information to indIDs
rs => ''#random number seed
);
my $n_arg = 0;#number of arguments provided by user
my $flag = '';#command flag
my $val = '';#value
my $ncols1 =0;#number of columns in first (not empty) line in a user file
my $ncols = 0;#number of columns in a line in user files
my $indfname = '';#infilename specified by flag -indf
#program settings
my $outdirname = 'export/';#name of directory for outfiles
my $settings_fname = 'export_seq_settings.txt';#outfile with settings and random number seed
my $popallfilename = 'popall.txt';#poploc-outfile from poploc
my $individuals_fname = 'individuals.txt';#out of indloc: list of all individual IDs in database
my $indpoplocsuffix = '_indpoploc.txt';#suffix of indpoploc-outfile from indpoploc
my $indpopallsuffix = '_indpopall.txt';#suffix of indpopall-outfile from indpoploc
my $group = 'group1';#used for table-output when user does not define a group for an individual (optional user file)
my $indall_outname = '_indall.txt';#name stem for table outfile
my $nid_outname = '_nid';#name stem for nonidentical sequences outfile
my $indseq_outname = '_indseq';#name stem for individual sequences outfile
#other variables
my $seed = '';
my %popall = ();# {poplocID}{popall_ID}=popallseq (merged)
my $poplocID = 0;
my $popall_ID = 0;
my $popallseq = '';
my %sel_ind = ();# {ind}=1
my %sel_loc = ();# {poplocID}=1
my %sel_gt = ();# {poplocID}{ind}=1
my %sel_all = ();# {poplocID}{popall_ID} = 1
my $ind = '';#individual ID
my $infilename = '';
my $outfilename = '';
my %indIDs = ();# {ind}=indID
my %indgroup = ();# {ind}=group group: from user file
my $indID = '';#indID
my $indID1 = '';#indID for seq 1 of an ind
my $indID2 = '';#indID for seq 2 of an ind
my $seq = '';
my $seq1 = '';
my $seq2 = '';
my %indpopall = ();# {ind}{poplocID}{popall_ID} = 1
my %popallind = ();# {poplocID}{popall_ID}{ind} = number of copies of this allele in this ind
my %nidseq = ();# {popall_ID}=seq; nonidentical allele sequences of one locus
my %indseq = ();# {$indID} = seq; individual sequences of one locus
my $loclen = 0;#length of a locus
my $mis_seq = '';#sequence of missing data symbols
my $success = 0;
my $tempstring1 = '';
my @temparr1 = ();
my $i = 0;
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
		#if flag defined, take over
		if (defined $user_settings{$flag}) {
			$user_settings{$flag} = $val;
		}
	}
}
#}

###################################################
#Test user settings for plausibility (limited test)
###################################################

#{
#-f must be f,n or p
if (length $user_settings{'f'} > 0) {
	unless (($user_settings{'f'} eq 'f') or ($user_settings{'f'} eq 'n') or ($user_settings{'f'} eq 'p')) {
		print "Implausible setting -f $user_settings{'f'}, using default: $defaults{'f'}\n";
		$user_settings{'f'} = $defaults{'f'};
	}
} else {
	$user_settings{'f'} = $defaults{'f'};
}
#-mi must be 0 or 1
if (length $user_settings{'mi'} > 0) {
	unless (($user_settings{'mi'} eq '0') or ($user_settings{'mi'} eq '1')) {
		print "Implausible setting -mi $user_settings{'mi'}, using default: $defaults{'mi'}\n";
		$user_settings{'mi'} = $defaults{'mi'};
	}
} else {
	$user_settings{'mi'} = $defaults{'mi'};
}
#-m must have length 1, test only when missing data will be included
if (($user_settings{'mi'} == 1) and (length $user_settings{'m'} > 0)) {
	unless (length $user_settings{'m'} == 1) {
		print "Implausible setting -m $user_settings{'m'}, using default: $defaults{'m'}\n";
		$user_settings{'m'} = $defaults{'m'};
	}
} else {
	$user_settings{'m'} = $defaults{'m'};
}
#indf: can you open the file?
if (length $user_settings{'indf'} > 0) {
	$indfname = $user_settings{'indf'};
	if(open(INFILE, "<", $indfname)) {
		close INFILE;
	} else {
		print "Setting -indf: Cannot open file $indfname. I will ignore it.\n";
		$user_settings{'indf'} = $defaults{'indf'};
		$indfname = $defaults{'indf'};
	}
}
#rs must be positive integer
if (length $user_settings{'rs'} > 0) {
	unless (($user_settings{'rs'} =~ /^\d+$/) and ($user_settings{'rs'} > 0)) {
		print "Setting -rs: must be a positive integer, exiting..\n";
		exit;
	}
}
#}

#######################
#Get random number seed
#######################

#{
$seed = $user_settings{'rs'};
if (length $seed > 0) {
	srand($seed);
} else {
	$seed = srand();
}
#}

####################################################
#Print used settings to file export_seq_settings.txt
####################################################

#{
$settings_fname = $outdirname . $settings_fname;
unless(open(SETOUT, ">>", $settings_fname)) {
	print "Cannot open file $settings_fname, exiting..\n";
	exit;
}
print SETOUT
"export_seq used these settings:\n",
"-f $user_settings{'f'}\n",
"-mi $user_settings{'mi'}\n",
"-m $user_settings{'m'}\n",
"-indf $user_settings{'indf'}\n",
"-rs $seed\n\n",
"You can use these settings to reproduce the same outfiles\n",
"or modify them to export the same data in another format.\n",
"Some flags above may not have values.\n",
"These flags are undefined per default and you did not define them.\n",
"Do not use these flags when you call the program to reproduce your outfiles.\n",
"The random number seed -rs is necessary to reproduce\n",
"the random phasing of alleles in outfiles with individuals sequences.\n",
"Please make sure to use the same selection of data with data_selector\n",
"to reproduce the same outfiles.\n";
#}

####################################################################
#Read in popall.txt, populate %popall and %sel_all (with 0 as value)
####################################################################

#{
print "Loading data...\n";
unless(open(INFILE,$popallfilename)) {
	print "Cannot open $popallfilename. Exiting..\n";
	print SETOUT "Could not open $popallfilename, did not export data.\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	#Infilecolumns!
	@temparr1 = split(/\t/,$tempstring1);#split into array
	$popall{$temparr1[0]}{$temparr1[1]} = $temparr1[2];
	$sel_all{$temparr1[0]}{$temparr1[1]} = 0;
}
close INFILE;
#}

######################################
#Read in existing preselections if any
#populate %sel_ind, %sel_loc, %sel_gt
######################################

#{
if (-d "export") {#if directory "export" exists
	#and if all three selection files are there
	if ((-f "export/sel_ind.txt") and (-f "export/sel_loc.txt") and (-f "export/sel_gt.txt")) {
		print "Loading your selection of individuals...\n";
		unless(open(INFILE,"export/sel_ind.txt")) {
			print "Cannot open export/sel_ind.txt. Exiting..\n";
			print SETOUT "Could not open export/sel_ind.txt, did not export data\n";
			exit;
		}
		@temparr1 = <INFILE>;
		close INFILE;
		for $ind (@temparr1) {
			chomp $ind;
			$sel_ind{$ind} = 1;#store selected inds in %sel_ind
		}
		#if the selection file was empty (%sel_ind still empty now)
		if (keys %sel_ind == 0) {
			print "Your selection evaluates to an empty dataset.\n",
			"I didn't export any data. Exiting..\n";
			print SETOUT "Your selection evaluates to an empty dataset.\n",
			"I didn't export any data.\n";
			exit;			
		}
		print "Loading your selection of loci...\n";
		unless(open(INFILE,"export/sel_loc.txt")) {
			print "Cannot open export/sel_loc.txt. Exiting..\n";
			print SETOUT "Could not open export/sel_loc.txt, did not export data.\n";
			exit;
		}
		while ($poplocID = <INFILE>) {#read in poplocIDs
			chomp $poplocID;
			$sel_loc{$poplocID} = 1;#store in %sel_loc
		}
		close INFILE;
		#if the selection file was empty (%sel_loc still empty now)
		if (keys %sel_loc == 0) {
			print "Your selection evaluates to an empty dataset.\n",
			"I didn't export data. Exiting..\n";
			print SETOUT "Your selection evaluates to an empty dataset.\n",
			"I didn't export data.\n";
			exit;			
		}
		print "Loading your selection of genotypes...\n";
		unless(open(INFILE,"export/sel_gt.txt")) {
			print "Cannot open export/sel_gt.txt. Exiting..\n";
			print SETOUT "Could not open export/sel_gt.txt, did not export data.\n";
			exit;
		}
		while ($tempstring1 = <INFILE>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			$sel_gt{$temparr1[0]}{$temparr1[1]} = 1;#store in %sel_gt
		}
		close INFILE;
		#if the selection file was empty (%sel_gt still empty now)
		if (keys %sel_gt == 0) {
			print "Your selection evaluates to an empty dataset.\n",
			"I didn't export data. Exiting..\n";
			print SETOUT "Your selection evaluates to an empty dataset.\n",
			"I didn't export data.\n";
			exit;			
		}
	}
	else {#if not all selection files are there
		if	(-f "export/sel_ind.txt") {
			print "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_ind.txt later. I will ignore it.\n";
			print SETOUT "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_ind.txt later. I have ignored it.\n";
		}
		if	(-f "export/sel_loc.txt") {
			print "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_loc.txt later. I will ignore it.\n";
			print SETOUT "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_loc.txt later. I have ignored it.\n";
		}
		if	(-f "export/sel_gt.txt") {
			print "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_gt.txt later. I will ignore it.\n";
			print SETOUT "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_gt.txt later. I have ignored it.\n";
		}
		print "I will export all genotypes in the database.\n";
	}
} else {#if directory "export" doesn't exist, create it
	unless (mkdir 'export') {
		print "Cannot create folder export. Exiting.\n";
		print SETOUT "Could not create folder export, did not export data.\n";
		exit;
	}
	print "I will export all genotypes in the database.\n";
}
#if %sel_ind is still empty, create it from file individuals.txt (using all inds)
if (keys %sel_ind == 0) {
	unless(open(INFILE,$individuals_fname)) {#open it
		print "Cannot open $individuals_fname. Exiting..\n";
		print SETOUT "Could not open $individuals_fname, did not export data.\n";
		exit;
	}
	@temparr1 = <INFILE>;#copy into array
	close INFILE;
	for $ind (@temparr1) {
		chomp $ind;
		$sel_ind{$ind} = 1;#store in %sel_ind
	}
}
#if %sel_loc is still empty, create from %popall (using all loci)
if (keys %sel_loc == 0) {
	for $poplocID (keys %popall) {
		$sel_loc{$poplocID} = 1;
	}
}
#if %sel_gt is still empty, create from individual files *_indpoploc.txt
#(using all genotypes from selected inds and locs)
if (keys %sel_gt == 0) {
	for $ind (keys %sel_ind) {
		#read in file *_indpoploc.txt for this ind (outfile from indpoploc)
		$infilename = $ind.$indpoplocsuffix;
		unless(open(INFILE,$infilename)) {
			print "Cannot open $infilename. Exiting..\n";
			print SETOUT "Could not open $infilename, did not export data.\n";
			exit;
		}
		$tempstring1 = <INFILE>;#skip headerline
		while ($tempstring1 = <INFILE>) {
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#Infilecolumns!
			$sel_gt{$temparr1[1]}{$ind} = 1;
		}
		close INFILE;
	}
}
#}

##########################################################
#If defined, read in file with additional strings
#to append to indIDs in outfiles with individual sequences
#if or if not, initialize %indIDs
##########################################################

#{
#initialize %indIDs and %indgroup
for $ind (keys %sel_ind) {
	$indIDs{$ind} = $ind;
	$indgroup{$ind} = $group;
}

if (length $user_settings{'indf'} > 0) {
	$indfname = $user_settings{'indf'};
	unless(open(INFILE, "<", $indfname)) {
		print "Cannot open $indfname, exiting..\n";
		print SETOUT "Could not open $indfname, did not export data.\n";
		exit;
	}
	#read in the file
	while($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		$tempstring1 =~ s/^\s*//;
		$tempstring1 =~ s/\s*$//;
		@temparr1 = split(/\t/,$tempstring1);
		if (@temparr1 != 2) {#ignore a line that has not 2 entries
			print "File $indfname: I ignore this line:\n$tempstring1\n";
			print SETOUT "File $indfname: I have ignored this line:\n$tempstring1\n";
		} else {
			if (defined $indIDs{$temparr1[0]}) {
				$tempstring1 = join("_",@temparr1);
				$indIDs{$temparr1[0]} = $tempstring1;
				$indgroup{$temparr1[0]} = $temparr1[1];
			} else {
				print "File $indfname: ind $temparr1[0] is not in current selection.\n",
				"I ignore it.\n";
				print SETOUT "File $indfname: ind $temparr1[0] is not in current selection.\n",
				"I have ignored it.\n";
			}
		}
	}
	close INFILE;
}
#}

#############################################################################################
#Read in data from individuals,
#populate %indpopall {ind}{poplocID}{popall_ID} = 1
#populate %popallind {poplocID}{popall_ID}{ind} = number of copies of this allele in this ind
#############################################################################################

#{
print "Loading genotype data from individuals...\n";
for $ind (keys %sel_ind) {#loop through selected individuals
	#read in file *_indpopall.txt for this ind (outfile from indpoploc)
	$infilename = $ind.$indpopallsuffix;#create filename
	unless(open(INFILE,$infilename)) {#open or die
		print "Cannot open $infilename. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip headerline
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		if ((exists $sel_gt{$temparr1[2]}) and (exists $sel_gt{$temparr1[2]}{$ind})) {#if poploc selected for this ind
			$indpopall{$ind}{$temparr1[2]}{$temparr1[3]} = 1;#store in %indpopall
			$popallind{$temparr1[2]}{$temparr1[3]}{$ind} = 1;#store in %popallind
		}
	}
	close INFILE;
}
#Change number of allele copies in %popallind to 2 in homozygotes
for $ind (keys %indpopall) {
	for $poplocID (keys %{$indpopall{$ind}}) {
		if (keys %{$indpopall{$ind}{$poplocID}} == 1) {
			$popall_ID = (sort keys %{$indpopall{$ind}{$poplocID}})[0];
			$popallind{$poplocID}{$popall_ID}{$ind} = 2;
		}
	}
}
#}

################################################
#Remove unselected loci and alleles from %popall
################################################

#{
for $ind (keys %indpopall) {
	for $poplocID (keys %{$indpopall{$ind}}) {
		for $popall_ID (keys %{$indpopall{$ind}{$poplocID}}) {
			$sel_all{$poplocID}{$popall_ID} = 1;
		}
	}
}
for $poplocID (keys %sel_all) {
	for $popall_ID (keys %{$sel_all{$poplocID}}) {
		if ($sel_all{$poplocID}{$popall_ID} == 0) {#if this allele did not occur in any ind
			delete $sel_all{$poplocID}{$popall_ID};#delete it
			delete $popall{$poplocID}{$popall_ID};#delete it
		}
	}
}
for $poplocID (keys %sel_all) {
	if (keys %{$sel_all{$poplocID}} == 0) {#if no alleles left in locus
		delete $sel_all{$poplocID};#delete locus
		delete $popall{$poplocID};#delete locus
	}
}
#}

################################
#For each selected locus:
#Assemble data, produce outfiles
################################

#{
for $poplocID (sort {$a <=> $b} keys %sel_loc) {
	#produce outfile *_indall.txt
	$outfilename = $outdirname . $poplocID . $indall_outname;
	unless(open(OUTTAB, ">", $outfilename)) {
		print "Cannot open $outfilename, exiting..\n";
		print SETOUT "Could not open $outfilename, stopped export.\n";
		exit;
	}
	#print header line
	print OUTTAB "poplocID\tpopall_ID\tind\tgroup\tn_copies\n";
	for $popall_ID (sort {$a <=> $b} keys %{$popallind{$poplocID}}) {
		for $ind (sort keys %{$popallind{$poplocID}{$popall_ID}}) {
			print OUTTAB "$poplocID\t$popall_ID\t$ind\t$indgroup{$ind}\t$popallind{$poplocID}{$popall_ID}{$ind}\n";
		}
	}
	close OUTTAB;
	#assemble nonidentical sequence data
	for $popall_ID (keys %{$popall{$poplocID}}) {
		$nidseq{$popall_ID} = $popall{$poplocID}{$popall_ID};
	}
	#assemble individual sequence data
	@temparr1 = ();
	$popall_ID = (sort {$a <=> $b} keys %{$popall{$poplocID}})[0];#select an allele
	$loclen = length $popall{$poplocID}{$popall_ID};#determine locus length
	$mis_seq = $user_settings{'m'} x $loclen;#build missing data seq
	for $ind (sort keys %sel_ind) {
		$indID1 = $indIDs{$ind} . '_1';
		$indID2 = $indIDs{$ind} . '_2';
		if (defined $sel_gt{$poplocID}{$ind}) {#if locus selected for this ind
			if (keys %{$indpopall{$ind}{$poplocID}} == 1) {#if ind is homozygous
				$popall_ID = (sort keys %{$indpopall{$ind}{$poplocID}})[0];
				$indseq{$indID1} = $popall{$poplocID}{$popall_ID};#store ID and seq in %indseq
				$indseq{$indID2} = $popall{$poplocID}{$popall_ID};#store ID and seq in %indseq
			} elsif (keys %{$indpopall{$ind}{$poplocID}} == 2) {#if ind is heterozygous
				for $popall_ID (sort {$a <=> $b} keys %{$indpopall{$ind}{$poplocID}}) {#collect the two popall_IDs
					push @temparr1, $popall_ID;
				}
				@temparr1 = shuffle @temparr1;#randomize their order
				$indseq{$indID1} = $popall{$poplocID}{$temparr1[0]};#store ID and seq in %indseq
				$indseq{$indID2} = $popall{$poplocID}{$temparr1[1]};#store ID and seq in %indseq
				@temparr1 = ();
			} else {#ind has more than 2 alleles
				print "Your selection contains genotypes with more than 2 alleles,\n",
				"e.g. locus $poplocID ind $ind, I stop export of data here.\n";
				print SETOUT "Selection contained genotypes with more than 2 alleles,\n",
				"e.g. locus $poplocID ind $ind, I stopped export of data.\n";
				exit;
			}
		} elsif ($user_settings{'mi'} == 1) {#locus not selected for this ind, missing data shall be included
			$indseq{$indID1} = $mis_seq;
			$indseq{$indID2} = $mis_seq;
		}		
	}
	#produce sequence outfiles
	if ($user_settings{'f'} eq 'f') {#FASTA format
		$success = fasta_out($poplocID,$outdirname,$nid_outname,\%nidseq);
		unless($success eq '1') {
			print "$success, exiting\n";
			print SETOUT "$success, stopped data export\n";
			exit;
		}
		$success = fasta_out($poplocID,$outdirname,$indseq_outname,\%indseq);
		unless($success eq '1') {
			print "$success, exiting\n";
			print SETOUT "$success, stopped data export\n";
			exit;
		}
	}
	elsif ($user_settings{'f'} eq 'n') {#Nexus format
		$success = nexus_out($poplocID,$outdirname,$nid_outname,$loclen,$user_settings{'m'},\%nidseq);
		unless($success eq '1') {
			print "$success, exiting\n";
			print SETOUT "$success, stopped data export\n";
			exit;
		}
		$success = nexus_out($poplocID,$outdirname,$indseq_outname,$loclen,$user_settings{'m'},\%indseq);
		unless($success eq '1') {
			print "$success, exiting\n";
			print SETOUT "$success, stopped data export\n";
			exit;
		}
	}
	elsif ($user_settings{'f'} eq 'p') {#PHYLIP format
		$success = phylip_out($poplocID,$outdirname,$nid_outname,$loclen,\%nidseq);
		unless($success eq '1') {
			print "$success, exiting\n";
			print SETOUT "$success, stopped data export\n";
			exit;
		}
		$success = phylip_out($poplocID,$outdirname,$indseq_outname,$loclen,\%indseq);
		unless($success eq '1') {
			print "$success, exiting\n";
			print SETOUT "$success, stopped data export\n";
			exit;
		}
	}
	%nidseq = ();
	%indseq = ();
}
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";
print SETOUT "Run took $run_s seconds.\n\n";
close SETOUT;

exit;

############
#Subroutines
############

#definition of sub fasta_out
#Expects:
#$poplocID: a locus ID
#$outdirname: name of outfiledir, ends with slash
#$fname: stem of outfilename
#ref to %ali: {ID}=seq
#
#prints sequences in %ali in FASTA format
#to an outfile named "outdir/poplocID_fname.fas"
#returns $success: everything oK: 1, else: error message
sub fasta_out {
	#declare and initialize: _r means a reference
	my ($poplocID,$outdirname,$fname,$ali_r) = @_;
	my $outfname = '';
	my $outsuff = '.fas';
	my $ID = '';
	my $seq = '';
	my $success = '1';
	
	#build outfilename and open outfile
	$outfname = $outdirname . $poplocID . $fname . $outsuff;
	unless(open(OUTFILE, ">", $outfname)) {
		$success = "Could not open $outfname.\n";
		return $success;
	}
	for $ID (sort keys %{$ali_r}) {
		$seq = $$ali_r{$ID};
		print OUTFILE ">$ID\n$seq\n";
	}
	close OUTFILE;
	return $success;
}

#definition of sub nexus_out
#Expects:
#$poplocID: a locus ID
#$outdirname: name of outfiledir, ends with slash
#$fname: stem of outfilename
#$loclen: sequence length
#$mischar: character for missing data
#ref to %ali: {ID}=seq
#
#prints sequences in %ali in Nexus format
#to an outfile named "outdir/poplocID_fname.nex"
#returns $success: everything oK: 1, else: error message
sub nexus_out {
	#declare and initialize: _r means a reference
	my ($poplocID,$outdirname,$fname,$loclen,$mischar,$ali_r) = @_;
	my $outfname = '';
	my $outsuff = '.nex';
	my $ntax = 0;
	my $gap = '-';
	my $ID = '';
	my $seq = '';
	my $success = '1';
	
	#build outfilename and open outfile
	$outfname = $outdirname . $poplocID . $fname . $outsuff;
	unless(open(OUTFILE, ">", $outfname)) {
		$success = "Could not open $outfname.\n";
		return $success;
	}
	#print section before data
	$ntax = keys %{$ali_r};
	if ($mischar eq '-') {
		$gap = 'N';
	}
	print OUTFILE "#NEXUS\n\nBegin data;\n",
	"\tDimensions ntax=$ntax nchar=$loclen;\n",
	"\tFormat datatype=dna gap=$gap missing=$mischar;\n",
	"\tMatrix\n";
	#print data
	for $ID (sort keys %{$ali_r}) {
		$seq = $$ali_r{$ID};
		print OUTFILE "$ID $seq\n";
	}
	#print section after data
	print OUTFILE "\t;\nEnd;\n";
	close OUTFILE;
	return $success;
}

#definition of sub phylip_out
#Expects:
#$poplocID: a locus ID
#$outdirname: name of outfiledir, ends with slash
#$fname: stem of outfilename
#$loclen: sequence length
#ref to %ali: {ID}=seq
#
#prints sequences in %ali in PHYLIP format
#to an outfile named "outdir/poplocID_fname.phy"
#returns $success: everything oK: 1, else: error message
#PHYLIP format: first line: space number of sequences space number of characters
#next lines: ID three spaces sequence
#IDs are padded with spaces to length 10 or length of longest ID (whatever is longer)
sub phylip_out {
	#declare and initialize: _r means a reference
	my ($poplocID,$outdirname,$fname,$loclen,$ali_r) = @_;
	my $outfname = '';
	my $outsuff = '.phy';
	my $ntax = 0;
	my $ID = '';
	my $seq = '';
	my $success = '1';
	my $lonID = '';#longest ID
	my $lonIDlen = 0;#length of longest ID
	my $minIDlen = 10;#minimum length of ID
	my $spacer = '   ';#spacer between ID and seq
	
	#build outfilename and open outfile
	$outfname = $outdirname . $poplocID . $fname . $outsuff;
	unless(open(OUTFILE, ">", $outfname)) {
		$success = "Could not open $outfname.\n";
		return $success;
	}
	#Determine longest ID length
	$lonID = (reverse sort {length($a) <=> length($b)} keys %{$ali_r})[0];
	$lonIDlen = length($lonID);
	if ($lonIDlen < $minIDlen) {
		$lonIDlen = $minIDlen;
	}
	#print first line
	$ntax = keys %{$ali_r};
	print OUTFILE " $ntax $loclen\n";
	#print data
	for $ID (sort keys %{$ali_r}) {
		$seq = $$ali_r{$ID};
		printf OUTFILE "%-*s%s%s\n", $lonIDlen, $ID, $spacer, $seq;
	}
	close OUTFILE;
	return $success;
}
