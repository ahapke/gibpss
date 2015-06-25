#!/usr/bin/perl -w
#createstdat version 03.0 Copyright 2015 Andreas Hapke
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

#uses subroutines
#simloc
#simall
#simdup
#simfreq
#simdel
#simread

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "createstdat version 03.0 Copyright 2015 Andreas Hapke\n",
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
#Variables for user settings
my %user_settings = (
nI => '20',#number of individuals
nL => '1200',#number of loci
nR => '40',#number of reads per locus and ind
Lsl => '120',#locus sequence length
Ld => '20',#minimum distance between loci
Rl => '80',#read length (default: int(Lsl/2)
nM => '6',#number of mutation steps for creation of alleles
e => '0.001',#sequencing/polymerase error rate
N => '0',#N rate (proportion of N in reads)
i => '200',#number of indel loci, which are selected at random among the nL loci
il => '3',#indel length
d => '0',#number of duplicated loci
rc => '1',#Boolean: print reads as revcomp with prob. 0.5
o => 'testdata',#output directory
rs => ''#random number seed
);
my %defaults = (
nI => '10',#number of individuals
nL => '100',#number of loci
nR => '40',#number of reads per locus and ind
Lsl => '120',#locus sequence length
Ld => '0',#minimum distance between loci
nM => '0',#number of mutation steps for creation of alleles
e => '0.001',#sequencing/polymerase error rate
N => '0',#N rate (proportion of N in reads)
i => '0',#number of indel loci, which are selected at random among the nL loci
il => '1',#indel length
d => '0',#number of duplicated loci
rc => '0',#Boolean: print reads as revcomp with prob. 0.5
o => 'sim_dat',#output directory
rs => ''#random number seed
);
my $n_arg = 0;#number command line arguments
my $flag = '';#a flag
my $val = '';#a value
#Output settings
my $user_files_dirname = 'user_files';#Name of directory for user_files
my $indloc_dist_fname = 'indloc_dist.txt';#Name of distance settings file for indloc
my $poploc_dist_fname = 'poploc_dist.txt';#Name of distance settings file for poploc
my $repfilename = 'sim_dat_rep.txt';#Name of report file
my $popall_fname = 'ori_popall.txt';#Name of popall outfile
my $poploc_fname = 'ori_poploc.txt';#Name of poploc outfile
my $gt_fname = 'ori_gt.txt';#Name of genotypes outfile
my $indsinfiles_fname = 'inds_infiles.txt';#Name of inds infiles outfile
my $popallfas_fname = 'popall.fas';#Name of fasta file with all alleles
my $intact_fname = 'intact_reads.txt';#Name of intact reads outfile
my $slash = "\/";
#Other variables
my $seed = '';
my %popall = ();# {poplocID}{popall_ID} = popallseq
my %popallfreq = ();# {poplocID}{popall_ID} = popallfreq
					#popallfreq: count in population of nI inds
my %poploc_indel = ();# {poplocID} = 0/1 1 for has indel variation
my %poploc_dup = ();# {poplocID} = 0/1 1 for locus has a duplicate
my %gt = ();# {ind}{poplocID}=(popall_ID,popall_ID)
my $ind = '';#ID of an individual
my $poplocID = 0;#ID of a locus
my $popall_ID = 0;#ID of an allele
my $popallseq = '';#sequence of an allele
my $popallseqnm = '';#unmerged sequence of an allele
my $f_read = '';#forward read
my $r_read = '';#reverse read (rev. comp.)
my $r_start = 0;#start position of r_read (in f-strand orientation)
my $head = '';#fastq header
my $fqname = '';#name of a fastq file
my $i = 0;
#}

########################
#Take over user-settings
########################

#{
if (@ARGV) {#if user provided some
	$n_arg = @ARGV;#determine number of arguments
	for ($i = 0; $i < ($n_arg - 1); $i += 2) {
		$flag = $ARGV[$i];
		$flag =~ s/^-//;
		$val = $ARGV[$i+1];
		#if flag is defined, take over,
		if (defined $user_settings{$flag}) {
			$user_settings{$flag} = $val;
		}
	}
}
#Check user settings for plausibility
#nI integer > 0
unless(($user_settings{'nI'} =~ /^\d+$/) and ($user_settings{'nI'} > 0)) {
	$user_settings{'nI'} = $defaults{'nI'};
}
#nL integer > 0
unless(($user_settings{'nL'} =~ /^\d+$/) and ($user_settings{'nL'} > 0)) {
	$user_settings{'nL'} = $defaults{'nL'};
}
#nR integer > 0
unless(($user_settings{'nR'} =~ /^\d+$/) and ($user_settings{'nR'} > 0)) {
	$user_settings{'nR'} = $defaults{'nR'};
}
#Lsl integer >= 20
unless(($user_settings{'Lsl'} =~ /^\d+$/) and ($user_settings{'Lsl'} >= 20)) {
	$user_settings{'Lsl'} = $defaults{'Lsl'};
}
#Ld integer >=0, <= Lsl/3
unless(($user_settings{'Ld'} =~ /^\d+$/) and ($user_settings{'Ld'} >= 0) and
($user_settings{'Ld'} <= ($user_settings{'Lsl'} / 3))) {
	$user_settings{'Ld'} = $defaults{'Ld'};
}
#i integer >= 0, <= nL
unless(($user_settings{'i'} =~ /^\d+$/) and ($user_settings{'i'} >= 0) and
($user_settings{'i'} <= $user_settings{'nL'})) {
	$user_settings{'i'} = $defaults{'i'};
}
if ($user_settings{'i'} == 0) {#if indels inactive
	$user_settings{'il'} = 0;#set indel length to 0
} else {#if indels active, test indel length setting
	#il integer >= 1, <= Lsl / 2
	unless(($user_settings{'il'} =~ /^\d+$/) and ($user_settings{'il'} >= 1) and
	($user_settings{'il'} <= ($user_settings{'Lsl'} / 2))) {
		$user_settings{'il'} = $defaults{'il'};
	}
}
#Rl integer >= 10, <= (Lsl - il), default int(Lsl/2)
unless(($user_settings{'Rl'} =~ /^\d+$/) and ($user_settings{'Rl'} >= 10) and
($user_settings{'Rl'} <= ($user_settings{'Lsl'} - $user_settings{'il'}))) {
	$user_settings{'Rl'} = int($user_settings{'Lsl'} / 2);
}
#nM integer >= 0, <= Lsl
unless(($user_settings{'nM'} =~ /^\d+$/) and ($user_settings{'nM'} >= 0) and
($user_settings{'nM'} <= $user_settings{'Lsl'})) {
	$user_settings{'nM'} = $defaults{'nM'};
}
#e decimal >= 0, <= 1
unless((($user_settings{'e'} =~ /^\d$/) or ($user_settings{'e'} =~ /^\d\.\d+$/))
and ($user_settings{'e'} >= 0) and ($user_settings{'e'} <= 1)) {
	$user_settings{'e'} = $defaults{'e'};
}
#N decimal >= 0, <= 1
unless((($user_settings{'N'} =~ /^\d$/) or ($user_settings{'N'} =~ /^\d\.\d+$/))
and ($user_settings{'N'} >= 0) and ($user_settings{'N'} <= 1)) {
	$user_settings{'N'} = $defaults{'N'};
}
#d number of duplicated loci, must be integer >= 0, <= nL-i, requires nM > 0
unless(($user_settings{'d'} =~ /^\d+$/) and ($user_settings{'d'} >= 0)
and ($user_settings{'d'} <= ($user_settings{'nL'} - $user_settings{'i'}))) {
	$user_settings{'d'} = $defaults{'d'};
}
if (($user_settings{'d'} > 0) and ($user_settings{'nM'} == 0)) {
	$user_settings{'d'} = $defaults{'d'};
}
#rc must be 0 or 1
unless(($user_settings{'rc'} eq '0') or ($user_settings{'rc'} eq '1')) {
	$user_settings{'rc'} = $defaults{'rc'};
}
#o if defined, try to create or exit; if not defined, try to create default output dir or exit
if (length $user_settings{'o'} > 0) {
	#remove trailing slash or backslash if any
	if ($user_settings{'o'} =~ m/\//) {
		$user_settings{'o'} =~ s/\/$//;
	}
	elsif ($user_settings{'o'} =~ m/\\/) {
		$user_settings{'o'} =~ s/\\$//;
	}
	if (length $user_settings{'o'} == 0) {
		print "unusable outdirname, using default $defaults{'o'}\n";
		$user_settings{'o'} = $defaults{'o'};
	}
} else {
	$user_settings{'o'} = $defaults{'o'};
}
if (-d "$user_settings{'o'}")  {
	print "Directory $user_settings{'o'} already exists.\n",
	"Please rename or delete. Exiting..\n";
	exit;
} else {
	unless (mkdir "$user_settings{'o'}") {
		print "Cannot create directory $user_settings{'o'}, exiting..\n";
		exit;
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


################################
#Create directory for user files
################################

#{
unless (-d "$user_files_dirname") {
	unless(mkdir "$user_files_dirname") {
		print "Cannot create directory $user_files_dirname, exiting..";
		exit;
	}
}
#}

#######################################
#Include outdirectory into outfilenames
#######################################

#{
if ($user_settings{'o'} =~ m/\\/) {#if outdirname contains backslash
	$slash = "\\";
}
$repfilename = $user_settings{'o'} . $slash . $repfilename;
$popall_fname = $user_settings{'o'} . $slash . $popall_fname;
$poploc_fname = $user_settings{'o'} . $slash . $poploc_fname;
$gt_fname = $user_settings{'o'} . $slash . $gt_fname;
$popallfas_fname = $user_settings{'o'} . $slash . $popallfas_fname;
$intact_fname = $user_settings{'o'} . $slash . $intact_fname;
$indsinfiles_fname = $user_files_dirname . $slash . $indsinfiles_fname;
$indloc_dist_fname = $user_files_dirname . $slash . $indloc_dist_fname;
$poploc_dist_fname = $user_files_dirname . $slash . $poploc_dist_fname;
#}

######################################################
#Produce distance settings files for indloc and poploc
######################################################

#{
unless(open(ILDIST, ">", $indloc_dist_fname)) {
	print "Cannot open $indloc_dist_fname, exiting..";
	exit;
}
print ILDIST "160\t12\t15\n";
close ILDIST;
unless(open(PLDIST, ">", $poploc_dist_fname)) {
	print "Cannot open $poploc_dist_fname, exiting..";
	exit;
}
print PLDIST "160\t15\n";
close PLDIST;
#}

#######################################################
#Open report file and print settings to file and screen
#######################################################

#{
#unless(open(OUTREP, ">$repfilename")) {
#	print "Cannot open file $repfilename, exiting ...\n\n";
#	exit;
#}
#print OUTREP "SIM_DAT used settings:\n";
#print OUTREP "nI $user_settings{'nI'}\n";
#print OUTREP "nL $user_settings{'nL'}\n";
#print OUTREP "nR $user_settings{'nR'}\n";
#print OUTREP "Lsl $user_settings{'Lsl'}\n";
#print OUTREP "Ld $user_settings{'Ld'}\n";
#print OUTREP "Rl $user_settings{'Rl'}\n";
#print OUTREP "nM $user_settings{'nM'}\n";
#print OUTREP "e $user_settings{'e'}\n";
#print OUTREP "N $user_settings{'N'}\n";
#print OUTREP "i $user_settings{'i'}\n";
#print OUTREP "il $user_settings{'il'}\n";
#print OUTREP "d $user_settings{'d'}\n";
#print OUTREP "rc $user_settings{'rc'}\n";
#print OUTREP "rs $seed\n";

#print  "Used settings:\n";
#print  "nI $user_settings{'nI'}\n";
#print  "nL $user_settings{'nL'}\n";
#print  "nR $user_settings{'nR'}\n";
#print  "Lsl $user_settings{'Lsl'}\n";
#print  "Ld $user_settings{'Ld'}\n";
#print  "Rl $user_settings{'Rl'}\n";
#print  "nM $user_settings{'nM'}\n";
#print  "e $user_settings{'e'}\n";
#print  "N $user_settings{'N'}\n";
#print  "i $user_settings{'i'}\n";
#print  "il $user_settings{'il'}\n";
#print  "d $user_settings{'d'}\n";
#print  "rc $user_settings{'rc'}\n";
#print  "rs $seed\n";
#}

##############
#Open outfiles
##############

#{
#unless(open(OUTPOPALL, ">$popall_fname")) {
#	print "Cannot open file $popall_fname, exiting ..\n";
#	exit;
#}
#unless(open(OUTPOPLOC, ">$poploc_fname")) {
#	print "Cannot open file $poploc_fname, exiting..\n";
#	exit;
#}
#unless(open(OUTGT, ">$gt_fname")) {
#	print "Cannot open file $gt_fname, exiting..\n";
#	exit;
#}
unless(open(OUTINDS, ">$indsinfiles_fname")) {
	print "Cannot open file $indsinfiles_fname, exiting..\n";
	exit;
}
#unless(open(OUTFAS, ">$popallfas_fname")) {
#	print "Cannot open file popall.fas, exiting ...\n\n";
#	exit;
#}
#}

######################################################
#Simulate data,print fq files and intact reads outfile
######################################################

#{
#Create first allele of each locus, populate %popall {poplocID}{popall_ID} = popallseq
print "Simulating loci..\n";
simloc(\%popall,$user_settings{'nL'},$user_settings{'Lsl'},
$user_settings{'Ld'},$user_settings{'nM'});
#Create nM additional alleles for each locus, add to %popall
print "Simulating alleles..\n";
simall(\%popall,$user_settings{'Lsl'},$user_settings{'nM'});
#Set "locus has duplicate" to 0 for all poplocs
for $poplocID (keys %popall) {
	$poploc_dup{$poplocID} = 0;
}
if ($user_settings{'d'} > 0) {#if loci shall be duplicated
	simdup(\%popall,\%poploc_dup,$user_settings{'d'});#duplicate d loci
}
#Create random allele frequencies for alleles (including 0)
#populate %popallfreq {poplocID}{popall_ID} = popallfreq
simfreq(\%popall,\%popallfreq,$user_settings{'nI'});
#Remove alleles with frequency 0
for $poplocID (keys %popallfreq) {
	for $popall_ID (keys %{$popallfreq{$poplocID}}) {
		if ($popallfreq{$poplocID}{$popall_ID} == 0) {
			delete $popall{$poplocID}{$popall_ID};
			delete $popallfreq{$poplocID}{$popall_ID};
		}
	}
}
#Set "has indel" to 0 for all poplocs
for $poplocID (keys %popall) {
	$poploc_indel{$poplocID} = 0;
}
#If selected, create indel variation
if ($user_settings{'i'} > 0) {
	print "Simulating indel variation..\n";
	simdel(\%popall,\%poploc_indel,\%poploc_dup,$user_settings{'Lsl'},$user_settings{'i'},$user_settings{'il'});
}
#Sample genotypes
print "Simulating genotypes..\n";
simgt(\%popallfreq,$user_settings{'nI'},\%gt);

#Create reads and print to fq files, create intact reads outfile
print "Simulating reads..\n";
simread(\%popall,\%gt,$user_settings{'nR'},$user_settings{'Rl'},
$user_settings{'e'},$user_settings{'N'},$user_settings{'rc'},$intact_fname,$user_settings{'o'},$slash);
#}

########################
#Print to other outfiles
########################

#{
#Print data to popall outfile
#print OUTPOPALL "poplocID\tpopall_ID\tpopallseq_merged\tpopallseq_notmerged\tfreq\n";
#for $poplocID (sort {$a <=> $b} keys %popall) {
#	for $popall_ID (sort {$a <=> $b} keys %{$popall{$poplocID}}) {
#		$popallseq = $popall{$poplocID}{$popall_ID};
#		$f_read = substr($popallseq,0,$user_settings{'Rl'});
#		$r_start = (length $popallseq) - $user_settings{'Rl'};
#		$r_read = substr($popallseq,$r_start,$user_settings{'Rl'});
#		$popallseqnm = $f_read . $r_read;
#		print OUTPOPALL "$poplocID\t$popall_ID\t$popallseq\t$popallseqnm\t$popallfreq{$poplocID}{$popall_ID}\n";
#	}
#}
#close OUTPOPALL;
#print data to poploc outfile
#print OUTPOPLOC "poplocID\tindel\tdup\n";
#for $poplocID (sort {$a <=> $b} keys %poploc_indel) {
#	print OUTPOPLOC "$poplocID\t$poploc_indel{$poplocID}\t$poploc_dup{$poplocID}\n";
#}
#close OUTPOPLOC;
#print data to gt outfile
#print OUTGT "ind\tpoplocID\tpopall_1\tpopall_2\n";
#for $ind (sort keys %gt) {
#	for $poplocID (sort {$a <=> $b} keys %{$gt{$ind}}) {
#		print OUTGT "$ind\t$poplocID\t$gt{$ind}{$poplocID}[0]\t$gt{$ind}{$poplocID}[1]\n";
#	}
#}
#close OUTGT;
#print data to inds_infiles outfile
for $ind (sort keys %gt) {
	$fqname = $user_settings{'o'} . $slash . $ind . '.fq';
	print OUTINDS "$ind\t$fqname\n";
}
close OUTINDS;

#for $poplocID (sort {$a <=> $b} keys %popall) {
#	for $popall_ID (sort {$a <=> $b} keys %{$popall{$poplocID}}) {
#		$popallseq = $popall{$poplocID}{$popall_ID};
#		$head = '>' . $poplocID . '_' . $popall_ID;
#		print OUTFAS "$head\n";
#		print OUTFAS "$popallseq\n";
#	}
#}
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";
#print OUTREP "Analysis completed.\nRun took $run_s seconds.\n\n";
#close OUTREP;

exit;

############
#Subroutines
############

#definition of subroutine simloc
#Creates random loci (one allele per locus)
#Loci have a configurable minimum distance to each other
#all have same length and equal base compositions (as close as possible)
#Expects:
#ref to %popall {poplocID}{popall_ID} = seq		will hold IDs and sequences
#nL  number of loci
#Lsl sequence length
#Ld  minimum distance between two loci
#nM  number of additional alleles that will be created for each locus later
#
#poplocIDs start with 1 and increment by 1
#popall_IDs start with 1 and increment by nM+1
sub simloc {
	#Declare and initialize, _r means a reference
	my ($popall_r,$nL,$Lsl,$Ld,$nM) = @_;
	my $poplocID = 1;
	my $popall_ID = 1;
	my $nA = 0;#number of A in sequence
	my $nC = 0;#number of C in sequence
	my $nG = 0;#number of G in sequence
	my $nT = 0;#number of T in sequence
	my $startseq = '';#starting sequence
	my @seq = ();#current sequence as array
	my $currseq = '';#current sequence of a locus
	my %poploc = ();# {poplocID}=seq
	my $dist = 0;#number of differences between two seqs
	my $i = 0;

	#create startseq
	$nA = int($Lsl / 4);
	$nC = int($Lsl / 4);
	$nG = int($Lsl / 4);
	$nT = $Lsl - $nA - $nC -$nG;
	$startseq = 'A' x $nA;
	$startseq .= 'C' x $nC;
	$startseq .= 'G' x $nG;
	$startseq .= 'T' x $nT;
	@seq = split('', $startseq);
	@seq = shuffle @seq;

	#create loci
	$dist = $Ld;
	while ($poplocID <= $nL) {
		@seq = shuffle @seq;#randomize characters in seq
		$currseq = join('',@seq);
		#test if $currseq is sufficiently distant from all known ones
		for $poplocID (keys %poploc) {
			$dist = $Lsl - (($currseq ^ $poploc{$poplocID}) =~ tr/\0//);
			last if ($dist < $Ld);
		}
		if ($dist >= $Ld) {#if sequence was sufficiently distant
			$poploc{$poplocID} = $currseq;
			++$poplocID;
		}
	}
	#populate %popall owned by main
	for $poplocID (sort {$a <=> $b} keys %poploc) {
		$$popall_r{$poplocID}{$popall_ID} = $poploc{$poplocID};
		$popall_ID += ($nM + 1);
	}	
}

#definition of subroutine simall
#Takes a set of loci, each with one allele (seq: ACGT, all capital)
#For each locus: creates nM additional alleles
#by stepwise mutation: one mutation per mutation step, at random position
#ensures that all alleles are different
#nM: number of mutation steps
#Expects:
#ref to %popall  {poplocID}{popall_ID}=seq contains one allele per locus
#Lsl sequence length (equal for all loci)
#nM  number of mutation steps
#adds the new alleles to %popall
#poplocIDs are ascending from 1
#popall_IDs in input start with 1 and increase by nM+1
#new popall IDs increment by 1 from that of first allele of each locus
sub simall {
	#Declare and initialize, _r means a reference
	my ($popall_r,$Lsl,$nM) = @_;
	my %subst = ();# {character}=array of substition characters, e.g. {A}=(C,G,T)
	my $ori_char = '';#original character at a position
	my $sub_char = '';#substition character
	my $poplocID = 0;
	my $popall_ID = 0;
	my $popallseq = '';
	my $mutpos = 0;#position of a mutation
	my $i = 0;
	my $j = 0;
	my %alleles = ();# {seq} = 1 alleles of one locus
	my $unique = 0;# 1 if a new allele is different from all others of the locus
	
	#initialize %subst
	@{$subst{'A'}} = ('C','G','T');
	@{$subst{'C'}} = ('A','G','T');
	@{$subst{'G'}} = ('A','C','T');
	@{$subst{'T'}} = ('A','C','G');
	
	for $poplocID (sort {$a <=> $b} keys %{$popall_r}) {#loop through loci
		$popall_ID = (sort {$a <=> $b} keys %{$$popall_r{$poplocID}})[0];
		$popallseq = $$popall_r{$poplocID}{$popall_ID};#get first allele
		$alleles{$popallseq} = 1;#store in %alleles
		for ($i = 0; $i < $nM; ++$i) {#perform $nM steps
			$unique = 0;
			while ($unique == 0) {
				$mutpos = int(rand $Lsl);#select a random position
				$ori_char = substr($popallseq,$mutpos,1);
				$j = int(rand 3);#select a random number 0..2 to choose a substition character
				$sub_char = $subst{$ori_char}[$j];#get substitution character
				substr($popallseq,$mutpos,1) = $sub_char;#substitute
				#test if new allele is different from already existing ones
				unless(defined $alleles{$popallseq}) {
					$unique = 1;
					++$popall_ID;
					$$popall_r{$poplocID}{$popall_ID} = $popallseq;#store new allele
					$alleles{$popallseq} = 1;
				}
			}
		}
		%alleles = ();
	}
}

#definition of subroutine simdup
#Expects
#ref to %popall  {poplocID}{popall_ID}=seq
#ref to %poploc_dup {poplocID}=0 will set 0 to 1 if locus is duplicate
#$ndup number of loci in %popall and %poploc_dup to duplicate
#Each locus must have at least 2 alleles, guaranteed by settings control
#Duplicates loci:
#Selects ndup loci at random, removes first allele, allele becomes a new, monomorphic locus
sub simdup {
	#Declare and initialize, _r means a reference
	my ($popall_r,$poploc_dup_r,$ndup) = @_;
	my $poplocID = 0;
	my $popall_ID = 0;
	my $seq = '';
	my $nextpoplocID = 0;#next ID for a new poploc
	my @all_poplocIDs = ();#all existing poplocIDs
	my $i = 0;
	
	#determine next poplocID
	$nextpoplocID = (reverse sort {$a <=> $b} keys %{$popall_r})[0];#get highest existing ID
	++$nextpoplocID;#increment
	
	@all_poplocIDs = (sort {$a <=> $b} keys %{$popall_r});#collect all existing poplocIDs
	@all_poplocIDs = shuffle @all_poplocIDs;#randomize their order
	for ($i = 0; $i < $ndup; ++$i) {#duplicate the first ndup loci
		$poplocID = shift @all_poplocIDs;#select a locus
		$popall_ID = (sort {$a <=> $b} keys %{$$popall_r{$poplocID}})[0];#get first popall_ID
		$seq = $$popall_r{$poplocID}{$popall_ID};#get sequence
		$$popall_r{$nextpoplocID}{$popall_ID} = $seq;#create new locus with that single allele
		delete $$popall_r{$poplocID}{$popall_ID};#delete first allele from original locus
		$$poploc_dup_r{$poplocID} = 1;#flag original locus as having duplicate
		$$poploc_dup_r{$nextpoplocID} = 1;#flag new duplicate locus as having duplicate
		++$nextpoplocID;
	}
}

#definition of subroutine simfreq
#Creates random allele frequencies for each allele
#Random allele frequencies include 0
#Frequencies are counts in population of nI individuals
#If a locus has at least 2 alleles, ensures that at least
#2 alleles get frequency > 0
#Expects
#ref to %popall {poplocID}{popall_ID} = popallseq
#ref to %popallfreq {poplocID}{popall_ID} = popallfreq
#$nI  number of individuals
#populates %popallfreq
sub simfreq {
	#Declare and initialize, _r means a reference
	my ($popall_r, $popallfreq_r, $nI) = @_;
	my $poplocID = 0;
	my $popall_ID = 0;
	my @all_IDs = ();#all popall_IDs of a locus
	my $first_all_ID = 0;#ID of first allele of a locus after order randomization
	my $first_all_freq = 0;#frequency of first allele
	my $last_all_ID = 0;#ID of last allele of a locus after order randomization
	my $last_all_freq = 0;#frequency of last allele
	my $freq = 0;#frequency
	my $remcount = 0;#total remaining count for alleles that have not yet a frequency
	
	for $poplocID (sort {$a <=> $b} keys %{$popall_r}) {
		$remcount = 2 * $nI - 1;
		$last_all_freq = 1;
		@all_IDs = (sort {$a <=> $b} keys %{$$popall_r{$poplocID}});#get popall_IDs of this loc
		if (@all_IDs > 1) {#if there is more than one allele, make sure that at least 2 alleles have frequency > 0 
			@all_IDs = shuffle @all_IDs;#randomize order of alleles of this loc
			$first_all_ID = shift @all_IDs;#separate first allele
			$first_all_freq = 1;#set frequency of first allele to 1
			$remcount -= 1;#reduce remaining count accordingly
			$freq = int(rand ($remcount + 1));#determine a random frequency
			$first_all_freq += $freq;#add to frequency of first allele
			$remcount -= $freq;#reduce remaining count accordingly
			$$popallfreq_r{$poplocID}{$first_all_ID} = $first_all_freq;
		}
		$last_all_ID = pop @all_IDs;#separate last allele
		for $popall_ID (@all_IDs) {
			$freq = int(rand ($remcount + 1));
			$$popallfreq_r{$poplocID}{$popall_ID} = $freq;
			$remcount -= $freq;
		}
		$last_all_freq += $remcount;
		$$popallfreq_r{$poplocID}{$last_all_ID} = $last_all_freq;
	}	
}

#definition of sub simdel
#Creates indel variation within existing loci: one deletion in first allele
#Selects i loci with indel variation at random among loci that are not duplicates
#Takes first allele
#Selects a position at random
#Deletes il characters starting at that position
#il: indel length minimum 1, max Lsl/2
#Expects
#ref to %popall {poplocID}{popall_ID}=popallseq
#ref to %poploc_indel {poplocID} = 0/1 1 for indel
#ref to %poploc_dup {poplocID} = 0/1 1 for locus has duplicate: don't introduce indel variation
#$Lsl original sequence length (equal for all loci and alleles)
#$i number of loci with indel variation (0 <= i <= number of loci)
#$il length of the indel
#
#Modifies the data in %popall
sub simdel {
	#Declare and initialize, _r means a reference
	my ($popall_r, $poploc_indel_r,$poploc_dup_r, $Lsl, $i, $il) = @_;
	my $ipos = 0;#first position to delete
	my $poplocID = 0;
	my $popall_ID = 0;
	my @allpoplocs = ();#IDs of all poplocs
	my $j = 0;
	
	#assemble poplocIDs and randomize their order
	@allpoplocs = (sort {$a <=> $b} keys %{$poploc_indel_r});
	@allpoplocs = shuffle @allpoplocs;
	
	while($j < $i) {#for i loci
		$poplocID = shift @allpoplocs;#select poploc
		if ($$poploc_dup_r{$poplocID} == 0) {#if locus has no duplicate
			#select first allele
			$popall_ID = (sort {$a <=> $b} keys %{$$popall_r{$poplocID}})[0];
			$ipos = int(rand ($Lsl + 1 - $il));#select random start for del (consider il)
			substr($$popall_r{$poplocID}{$popall_ID},$ipos,$il) = '';#delete
			$$poploc_indel_r{$poplocID} = 1;
			++$j;
		}
	}
}

#definition of subroutine simgt
#Creates genotypes
#Expects:
#ref to %popallfreq {poplocID}{popall_ID}=freq (count of allele in nI individuals)
#$nI number of individuals
#ref to %gt {ind}{poplocID}=(popall_ID,popall_ID)
#
#Sum of freq for one poploc in %popallfreq must be equal to 2*$nI
#Populates %gt
sub simgt {
	#Declare and initialize, _r means a reference
	my ($popallfreq_r,$nI,$gt_r) = @_;
	my $poplocID = 0;
	my $popall_ID = 0;
	my $popall_ID1 = 0;
	my $popall_ID2 = 0;
	my @gt = ();
	my $freq = 0;
	my @popalls = ();
	my %inds = ();# {ind} = 1
	my $ind = '';
	my $i = 0;
	
	#Create individual IDs
	for ($i = 1; $i <= $nI; ++$i) {
		$ind = sprintf('%04d',$i);
		$ind = 'ind' . $ind;
		$inds{$ind} = 1;
	}
	#Create genotypes
	for $poplocID (sort {$a <=> $b} keys %{$popallfreq_r}) {
		for $popall_ID (sort {$a <=> $b} keys %{$$popallfreq_r{$poplocID}}) {
			$freq = $$popallfreq_r{$poplocID}{$popall_ID};
			#build array of popall_IDs, each one $freq times in array
			for ($i = 0; $i < $freq; ++$i) {
				push @popalls, $popall_ID;
			}
		}
		@popalls = shuffle @popalls;#randomize order of popall_IDs
		#draw individual genotypes for this loc
		for $ind (sort keys %inds) {
			$popall_ID1 = shift @popalls;
			$popall_ID2 = shift @popalls;
			@gt = sort {$a <=> $b} ($popall_ID1,$popall_ID2);
			@{$$gt_r{$ind}{$poplocID}} = @gt;
		}		
	}
}

#Definition of subroutine simread
#Simulates Reads and prints them to fastq files
#For each ind and locus
#Takes a genotype with two allele sequences
#Creates nR reads, with random proportions of the two alleles
#Creates f-read and r-read (in f-orientation: revcomp) and concatenates them to one read
#f- and r- reads have length Rl
#Introduces errors and N for missing data with configurable probabilities
#Expects
#ref to %popall {poplocID}{popall_ID}=popallseq
#ref to %gt {ind}{poplocID}=(popall_ID,popall_ID) popall_IDs sorted ascending
#$nR	number of reads per locus and ind
#$Rl	read length
#$e		error rate
#$N		N rate
#$rc    1: print reads as reverse complement with prob. 0.5, 0: don't
#$intact_fname name of intact_reads outfile
#$outdirname name of outfile directory for fastq files
#$slash contains slash or backslash for building of fastq outfilenames
sub simread {
	#Declare and initialize, _r means a reference
	my ($popall_r,$gt_r,$nR,$Rl,$e,$N,$rc,$intact_fname,$outdirname,$slash) = @_;
	my $poplocID = 0;
	my $popall_ID = 0;
	my $popallseq = '';
	my $ind = '';
	my $fqsuff = '.fq';#suffix for fastq outfile
	my $fqname = '';#fastq filename
	my @gt = ();#one genotype (two allele IDs)
	my @alleles = ();#both allele IDs of genotype, each nR times
	my $f_read = '';#f read
	my $r_read = '';#r read (in f orientation)
	my $r_start = 0;#start of r read in popallseq
	my %catseqs = ();# {popall_ID} = catseq
	my $ori_char = '';#original character at a position
	my $err_char = '';#erroneous character
	my %error = ();# {ori_char}=array of possible err_chars
	my $pos = 0;#position in sequence
	my $head = '';#fastq header
	my $phredline = '';#phredscore line
	my $readNo = 0;#number of read for fastq_header
	my $readintact = 0;#1 if read is intact
	my %thisgtintact = ();# {popall_ID} = number of intact reads for each allele
	my %intact = ();# {poplocID}{ind} = greater number of intact reads from both alleles
	my $nintact = 0;#greater number of intact reads from both alleles
	
	my $i = 0;
	my $j = 0;
	
	#Initialize %error
	@{$error{'A'}} = ('C','G','T');
	@{$error{'C'}} = ('A','G','T');
	@{$error{'G'}} = ('A','C','T');
	@{$error{'T'}} = ('A','C','G');
	#Build phredscore line
	$phredline = 'D' x ($Rl * 2);
	
	#open intact reads outfile and print headerline
	#unless(open(OUTINTACT, ">$intact_fname")) {
	#	print "Cannot open $intact_fname. Exiting..\n";
	#	exit;
	#}
	#print OUTINTACT "poplocID";
	#for $ind (sort keys %{$gt_r}) {
	#	print OUTINTACT "\t$ind";
	#}
	#print OUTINTACT "\n";
	
	#create concatenated f- r- reads for each popall
	for $poplocID (keys %{$popall_r}) {
		for $popall_ID (keys %{$$popall_r{$poplocID}}) {
			$popallseq = $$popall_r{$poplocID}{$popall_ID};
			$f_read = substr($popallseq,0,$Rl);
			$r_start = (length $popallseq) - $Rl;
			$r_read = substr($popallseq,$r_start,$Rl);
			$catseqs{$popall_ID} = $f_read . $r_read;
		}
	}
	
	for $ind (sort keys %{$gt_r}) {
		print "Individual $ind\n";
		#create outfilename and open outfile or die
		$fqname = $outdirname . $slash . $ind . $fqsuff;
		unless(open(OUTFQ, ">$fqname")) {
			print "Cannot open $fqname. Exiting..";
			exit;
		}
		for $poplocID (sort {$a <=> $b} keys %{$$gt_r{$ind}}) {
			@gt = @{$$gt_r{$ind}{$poplocID}};
			#sample each popall_ID of genotype nR times into @alleles
			for ($i = 0; $i < $nR; ++$i) {
				push @alleles, $gt[0];
				push @alleles, $gt[1];
			}
			@alleles = shuffle @alleles;#randomize their order
			#Create reads
			for ($i = 0; $i < $nR; ++$i) {
				$readintact = 1;
				$popallseq = $catseqs{$alleles[$i]};#get concatenated f-r read for allele
				#if reads shall be saved as revcomp with probability 0.5,do
				if (($rc == 1) and ((rand) < 0.5)) {
					$popallseq = reverse $popallseq;
					$popallseq =~ tr/ACGTacgt/TGCAtgca/;
				}
				#If error rate > 0: Create sequencing errors
				if ($e > 0) {
					for ($pos = 0; $pos < length $popallseq; ++$pos) {
						if ((rand) < $e) {#select error position with probability $e
							$readintact = 0;
							$ori_char = substr($popallseq,$pos,1);
							$j = int(rand 3);#select a random number 0..2 to choose an error character
							$err_char = $error{$ori_char}[$j];#get error character
							substr($popallseq,$pos,1) = $err_char;#substitute
						}
					}
				}
				#If N rate > 0: Create missing characters
				if ($N > 0) {
					for ($pos = 0; $pos < length $popallseq; ++$pos) {
						if ((rand) < $N) {#select N position with probability $N
							$readintact = 0;
							substr($popallseq,$pos,1) = 'N';#substitute
						}
					}
				}
				#count intact reads
				$thisgtintact{$alleles[$i]} += $readintact;
				$nintact = (reverse sort {$a <=> $b} values %thisgtintact)[0];
				$intact{$poplocID}{$ind} = $nintact;
				#build fastq header
				$readNo = $i + 1;
				$head = '@1_' . $poplocID . '_' . $alleles[$i] . '_' . $readNo . '_1';
				#print read
				print OUTFQ "$head\n$popallseq\n+\n$phredline\n";
			}			
			@alleles = ();#set back
			%thisgtintact = ();#set back and next loc
		}
		close OUTFQ;#and next ind
	}
	#print data to intact outfile
	#for $poplocID (sort {$a <=> $b} keys %intact) {
	#	print OUTINTACT "$poplocID";
	#	for $ind (sort keys %{$intact{$poplocID}}) {
	#		print OUTINTACT "\t$intact{$poplocID}{$ind}";
	#	}
	#	print OUTINTACT "\n";
	#}
	#close OUTINTACT;
}
