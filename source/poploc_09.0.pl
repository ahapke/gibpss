#!/usr/bin/perl -w
#poploc version 09.0 Copyright 2015 Andreas Hapke
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

#uses subroutines
#idpoploc
#splitvar
#fno_f_ID
#popall_pairs
#dist
#al_poploc
#pop_outfiles
#merge_poploc
#overlap1
#varpos
#consensus_IUPAC

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "poploc version 09.0 Copyright 2015 Andreas Hapke\n",
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
my %user_settings = (
d => '',#distance-filename
msl => 0,#sequence length of loci that shall be merged, 0: don't merge
mino => 10,#minimum overlap for merging
spl => 0#shortest plausible length of merged loci (will be changed to msl/2 if msl>0 and no plausible value for spl provided)
);
my %defaults = (
d => '',#distance-filename
msl => 0,#sequence length of loci that shall be merged, 0: don't merge
mino => 10,#minimum overlap for merging
spl => 0#shortest plausible length of merged loci (will be changed to msl/2 if msl>0 and no plausible value for spl provided)
);
my $n_arg = 0;#number of arguments provided by user
my $flag = '';#command flag
my $val = '';#value for command flag
my $msl = 0;#sequence length of loci that shall be merged, 0: don't merge
my $mino = 0;#minimum overlap required for merging of alleles
my $maxo = 0;#maximum overlap allowed for merging of alleles
my $individuals_fname = 'individuals.txt';#name of file with list of all individuals in database
										#out of indloc
my $allelefilesuffix = '_alleles.txt';#suffix of allele-outfiles from indloc
my $report_name = 'poploc_report.txt';#outfile: used settings
my $dist_fname = '';#distance-filename
my $inter_ind_dist = 6;#maximum distance allowed between two neighboring alleles
						#in a poplocus-network, default: 6
my %sl_dist = ();# {sl}=inter_ind_dist
my $sl_col = 2;#col of sl in infile
my $allseq_col = 4;#col of allseq in infile
my $indlocID_col = 0;#col of indlocID in infile
my $indall_ID_col = 3;#col of indall_ID in infile
#other variables
my $ind = '';#ID of one ind
my $infile = '';#name of an infile
my %ind_infile = ();# {ind}=infilename
my %ind_all_in = ();# {sl}{allseq}{ind}=@[0] indlocID, [1] indall_ID
my $sl = 0;#seqlength
my $allseq = '';#seq of an allele
my $indlocID = 0;#locus ID in an individual
my $indall_ID = 0;#allele ID in an individual
my %ind_all = ();# {sl}{popall_ID}:
				# {seq}=allseq
				# {inds}{ind}= @ of [0] indlocID [1] indall_ID
my $popall_ID = 1;# allele ID at the population level
my $poplocID = 0;#ID of a locus at the population level
my $next_poplocID = 1;#first poplocID in the next set of loci (with a given seqlength)
my %poploc = ();# {sl}{poplocID} = @ of popall_ID
my $success = 'analysis aborted';#success message from last subroutine

my @temparr1 = ();
my $tempstring1 = '';
my $i = 0;
#}

########################
#Take over user settings
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
#}

#################
#Open report file
#################

#{
unless(open(OUTREP, ">$report_name")) {
	print "Cannot open file $report_name, exiting ...\n\n";
	exit;
}
print OUTREP "POPLOC: used settings:\n";
#}

#####################################################
#Check settings for plausibility or change to default
#Determine minimum and maximum overlap for merging
#Print settings and error messages to report
#####################################################

#{
#Distance settings
$dist_fname = $user_settings{'d'};#look up name of distance settings file
if (length $dist_fname > 0) {#if user provided one
	unless(open(DISTFILE,$dist_fname)) {#try to open the file
		#Could not open it:
		print "Could not open $dist_fname.\n";
		print OUTREP "Could not open $dist_fname.\n";
	}
	while ($tempstring1 = <DISTFILE>) {#read in the file
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);#split line into array
		#take only lines with 2 entries
		#check entries for plausibility, all integers, sl >= 0, distance >= 0
		if ((@temparr1 == 2) and ($temparr1[0] =~ /^\d+$/) and ($temparr1[1] =~ /^\d+$/) and
		($temparr1[0] >= 0) and ($temparr1[1] >= 0)) {
			#take over into %sl_dist
			$sl_dist{$temparr1[0]} = $temparr1[1];
		} else {#unusable line
			print "Ignored unusable line in distance settings file: $tempstring1\n";
			print OUTREP "Ignored unusable line in distance settings file: $tempstring1\n";
		}
	}
	close DISTFILE;	
} else {#if user did not provide the name of a distance settings file
	print "You did not provide the name of a distance settings file";
	print OUTREP "You did not provide the name of a distance settings file";
}
#if you have no distance settings now
if (keys %sl_dist == 0) {
	print "I use the default distance for all sequence lengths: $inter_ind_dist\n";
	print OUTREP "I use the default distance for all sequence lengths: $inter_ind_dist\n";
	#initialize smallest possible seqlength with default distance
	$sl = 1;
	$sl_dist{$sl} = $inter_ind_dist;
}
#print distance settings to report-file
print OUTREP "Distance settings:\nseqlength distance\n";
for $sl (sort {$a <=> $b} keys %sl_dist) {
	print OUTREP "$sl $sl_dist{$sl}\n";
}
#-msl even integer >= 0
unless(($user_settings{'msl'} =~ /^\d+$/) and ($user_settings{'msl'} >= 0)
and ($user_settings{'msl'} % 2 == 0)) {
	print "-msl $user_settings{'msl'} not plausible.\n";
	print OUTREP "-msl $user_settings{'msl'} not plausible.\n";
	$user_settings{'msl'} = $defaults{'msl'};
}
#if -msl > 0, check -mino and -spl
if ($user_settings{'msl'} > 0) {
	#-mino positive integer <= (msl/2), if not: set to default 10 or smaller appropriate value
	unless(($user_settings{'mino'} =~ /^\d+$/) and ($user_settings{'mino'} > 0)
	and ($user_settings{'mino'} <= ($user_settings{'msl'} / 2))) {
		print "-mino $user_settings{'mino'} not plausible.\n";
		print OUTREP "-mino $user_settings{'mino'} not plausible.\n";
		#if length of loci to merge >= default minimum overlap * 2
		if ($user_settings{'msl'} >= ($defaults{'mino'} * 2)) {
			$user_settings{'mino'} = $defaults{'mino'};
			print "I use default: -mino $user_settings{'mino'}\n";
			print OUTREP "I use default: -mino $user_settings{'mino'}\n";
		} else {#if not (that won't really happen in a realistic analysis)
			$user_settings{'mino'} = ($user_settings{'msl'} / 2);#set minimum overlap to half of locus length
			print "Setting -msl $user_settings{'msl'} is too small for default -mino.\n",
			"I use this setting: -mino $user_settings{'mino'}\n";
			print OUTREP "Setting -msl $user_settings{'msl'} is too small for default -mino.\n",
			"I use this setting: -mino $user_settings{'mino'}\n";
		}
	}
	#-spl integer >= msl/2 <= msl, if not, set to msl/2
	unless(($user_settings{'spl'} =~ /^\d+$/) and ($user_settings{'spl'} >= ($user_settings{'msl'} / 2))
	and ($user_settings{'spl'} <= $user_settings{'msl'})) {
		print "-spl $user_settings{'spl'} not plausible.\n";
		print OUTREP "-spl $user_settings{'spl'} not plausible.\n";
		$user_settings{'spl'} = ($user_settings{'msl'} / 2);
		print "I use the smallest possible value: -spl $user_settings{'spl'}\n";
		print OUTREP "I use the smallest possible value: -spl $user_settings{'spl'}\n";
	}
	#Check if merging is possible with these settings
	if ($user_settings{'mino'} > ($user_settings{'msl'} - $user_settings{'spl'})) {#if impossible
		print "Merging is not possible with these settings:\n",
		"-msl $user_settings{'msl'} -mino $user_settings{'mino'} -spl $user_settings{'spl'}\n";
		print OUTREP "Merging is not possible with these settings:\n",
		"-msl $user_settings{'msl'} -mino $user_settings{'mino'} -spl $user_settings{'spl'}\n";
		$user_settings{'msl'} = 0;#inactivate merging
	} else {#if possible, determine variables
		$msl = $user_settings{'msl'};
		$maxo = $user_settings{'msl'} - $user_settings{'spl'};
		$mino = $user_settings{'mino'};
	}
}
if ($user_settings{'msl'} == 0) {#if merging inactive now
	print "I will not merge any sequences.\n";
	print OUTREP "I will not merge any sequences.\n";
} else {#if merging active now
	print "I will merge sequences with length $user_settings{'msl'}\n",
	"minimum overlap: $user_settings{'mino'}\n",
	"shortest allowed locus length after merging: $user_settings{'spl'}\n",
	"maximum overlap: $maxo\n";
	print OUTREP "I will merge sequences with length $user_settings{'msl'}\n",
	"minimum overlap: $user_settings{'mino'}\n",
	"shortest allowed locus length after merging: $user_settings{'spl'}\n",
	"maximum overlap: $maxo\n";
}
#}

###################
#get individual-IDs
###################

#{
unless(open(INFILE,$individuals_fname)) {
	print "Cannot open $individuals_fname. Exiting..\n";
	print OUTREP "Could not open $individuals_fname. Analysis aborted.\n";
	exit;
}
@temparr1 = <INFILE>;#put into array
close INFILE;
for $ind (@temparr1) {
	chomp $ind;
	$infile = $ind.$allelefilesuffix;#build allele-infilename
	$ind_infile{$ind} = $infile;#store in %ind_infile
}
#}

#########################################
#Read in file *_alleles.txt from each ind
#########################################

#{
while (($ind, $infile) = each %ind_infile) {#loop through inds and their infiles
	unless(open(INFILE,$infile)) {
		print "Cannot open $infile. Exiting..\n";
		print OUTREP "Could not open $infile. Analysis aborted.\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip header line
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;#remove line ending
		@temparr1 = split(/\t/,$tempstring1);#split into array
		$sl = $temparr1[$sl_col];#get seqlength out of its field (settings)
		$allseq = $temparr1[$allseq_col];#get sequence out of its field (settings)
		$indlocID = $temparr1[$indlocID_col];#get indlocID out of its field (settings)
		$indall_ID = $temparr1[$indall_ID_col];#get indall_ID out of its field (settings)
		if ($allseq ne '0') {#if an allele is defined in this line
			#add data to %ind_all_in:
			@{$ind_all_in{$sl}{$allseq}{$ind}} = ($indlocID,$indall_ID);
		}
	}
}
#}

######################################################
#Transfer input-data into number-indexed datastructure
######################################################

#{
for $sl (sort {$a <=> $b} keys %ind_all_in) {#loop through sl: key1 of %ind_all_in
	for $allseq (keys %{$ind_all_in{$sl}}) {#loop through allseq: key2 of %ind_all_in
		$ind_all{$sl}{$popall_ID}{"seq"} = $allseq;#add sequence to %ind_all
		for $ind (keys %{$ind_all_in{$sl}{$allseq}}) {#loop through ind: key3 of %ind_all_in
			@temparr1 = @{$ind_all_in{$sl}{$allseq}{$ind}};#get value array of %ind_all_in
			@{$ind_all{$sl}{$popall_ID}{"inds"}{$ind}} = @temparr1;#store in %ind_all
		}
		++$popall_ID;#increment for next popall
	}
}
%ind_all_in = ();#We don't need %ind_all_in any more.
#}

##################
#Identify pop-loci
##################

#{
#determine smallest defined sequence length in %sl_dist
$sl = (sort {$a <=> $b} keys %sl_dist)[0];
#initialize inter_ind_dist with distance for this seqlength
$inter_ind_dist = $sl_dist{$sl};

for $sl (sort {$a <=> $b} keys %ind_all) {#loop through seqlengths: key1 %ind_all
	#if current sequence length is defined in %sl_dist
	if (defined $sl_dist{$sl}) {
		#update inter_ind_dist
		$inter_ind_dist = $sl_dist{$sl};
	}
	#call sub idpoploc to identify loci with this seqlength
	#and populate %poploc
	idpoploc($sl,$inter_ind_dist,\$next_poplocID,\%ind_all,\%poploc);
}
#}

###########################################
#produce outfiles poploc.txt and popall.txt
###########################################

$success = pop_outfiles(\%poploc,\%ind_all,$msl,$mino,$maxo);

##################################
#print final information to report
##################################

#{
#if sub pop_outfiles returned positive success-message
if ($success =~ /analysis completed/i) {
	#print used individuals to report
	print OUTREP "Analysis based on data from individuals:\n";
	for $ind (sort keys %ind_infile) {
		print OUTREP "$ind\n";
	}
} else {#if no positive success-message
	print "Analysis not completed\n";
	print "$success";
	print OUTREP "Analysis not completed\n";
	print OUTREP "$success";
}
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";
print OUTREP "Run took $run_s seconds.\n";
close OUTREP;

exit;

############
#Subroutines
############

#definition of subroutine idpoploc
#identifies loci at population level
#all loci have same seqlength
#expects:
#$sl: seqlength
#$inter_ind_dist: #maximum distance allowed between two neighboring alleles
				#in a poplocus-network
#ref to $next_poplocID: 1st poplocID to use
#ref to %ind_all:	# {sl}{popall_ID}:
					# {seq}=allseq
					# {inds}{ind}= @ of [0] indlocID [1] indall_ID
#ref to %poploc: {sl}{poplocID} = @ of popall_ID
#
#populates %poploc for this $sl and updates $next_poplocID
sub idpoploc {
	#declare and initialize: _r means a reference
	my ($sl,$inter_ind_dist,$next_poplocID_r,$ind_all_r,$poploc_r) = @_;
	my @splitvar = ();#variables for splitting of alleles into nonoverlapping fragments:
					#2d array: [0]: array of startpositions [1] array of fragment-lengths						
	my %fno_f_ID = ();# {fno}{frag}=@of popall_IDs
	my %pairs = ();# {popall_ID1}{popall_ID2}=dist
	
	#determine variables for splitting of alleles into fragments
	@splitvar = splitvar($sl,$inter_ind_dist);
	#call subroutine fno_f_ID to build fragment catalog in %fno_f_ID
	fno_f_ID($sl,\@splitvar,$ind_all_r,\%fno_f_ID);
	#call subroutine popall_pairs to identify pairs of alleles (popall_IDs) with distance <= $inter_ind_dist
	#populate %pairs
	popall_pairs($sl,$inter_ind_dist,$ind_all_r,\%fno_f_ID,\%pairs);
	%fno_f_ID = ();#We don't need the fragment catalog any more.
	#call subroutine al_poploc to assemble alleles into loci
	#populate %poploc owned by main for this seqlength
	#update $next_poploc in main for next seqlength
	al_poploc($sl,$ind_all_r,$next_poplocID_r,\%pairs,$poploc_r);	
}

#definition of subroutine splitvar
#expects
#$sl: sequence length
#$hapdist: max dist between neighboring allele sequences in a locus network
#determines variables for splitting of a sequence into nonoverlapping fragments:
#$hapdist+1 fragments
#returns 2d-@splitvar:[0]: array of startpositions [1] array of fragment-lengths
#the last fragment has a greater length than the others
#if the seqlength cannot be divided into $hapdist+1 fragments
sub splitvar {
	#declare and initialize
	my ($sl, $hapdist) = @_;
	my @splitvar = ();#see above
	my $nfrag = 0;#number of fragments
	my $fragl = 0;#length of all but last fragment
	my $lastfragl = 0;#length of last fragment
	my $pos = 0;#start position of a fragment
	my $i = 0;#counter
	my @temparr = ();#temporary array	
	
	#determine number of fragments
	$nfrag = $hapdist + 1;
	#determine fragment lengths
	$fragl = int($sl / $nfrag);#lengths of all but last fragment
	$lastfragl = $sl - ($fragl * ($nfrag - 1));
	#determine start positions
	for ($i = 0; $i < $nfrag; ++$i) {
		$temparr[$i] = $pos;
		$pos += $fragl;
	}
	#add start positions to @splitvar
	push @splitvar, [@temparr];
	
	#put fragment-lengths into array
	#all but last fragment
	for ($i = 0; $i < ($nfrag - 1); ++$i) {
		$temparr[$i] = $fragl;
	}
	#last fragment:
	$temparr[$nfrag-1] = $lastfragl;
	#add fragment lengths to @splitvar
	push @splitvar, [@temparr];
		
	return @splitvar;
}

#definition of subroutine fno_f_ID
#builds a fragment catalog of alleles with a given seqlength
#expects:
#$sl: seqlength
#ref to @splitvar: #variables for splitting of alleles into nonoverlapping fragments:
					#2d array: [0]: array of startpositions [1] array of fragment-lengths
#ref to %ind_all:	# {sl}{popall_ID}: needed
					# {seq}=allseq: needed
					# {inds}{ind}= @ of [0] indlocID [1] indall_ID: not needed
#ref to %fno_f_ID: {fno}{frag}=@of popall_IDs
#populates %fno_f_ID
sub fno_f_ID {
	#declare and initialize: _r means a reference
	my ($sl,$splitvar_r,$ind_all_r,$fno_f_ID_r) = @_;
	my @splitvar = @{$splitvar_r};#copy of (outer array of) array splitvar
	my $fno = 0;#current fragment number
	my $popall_ID = 0;#ID of an allele
	my $allseq = '';#sequence of current allele
	my $frag = '';#a fragment of allele as key2 for %fno_f_ID
	my $i = 0;
	
	#loop through the popall_IDs (key2 in %ind_all
	for $popall_ID (keys %{$$ind_all_r{$sl}}) {
		$allseq = $$ind_all_r{$sl}{$popall_ID}{"seq"};#get allele-sequence
		#loop through fragment numbers
		for ($fno = 0; $fno < @{$splitvar[0]}; ++$fno) {
			#build current fragment
			$frag = substr($allseq,${$splitvar[0]}[$fno],${$splitvar[1]}[$fno]);
			#store $frag and $popall_ID in %fno_f_ID
			push @{$$fno_f_ID_r{$fno}{$frag}}, $popall_ID;
		}
	}	
}

#definition of subroutine popall_pairs
#identifies pairs of alleles (popall_IDs) with a distance <= $inter_ind_dist
#expects:
#$sl: seqlength
#$inter_ind_dist: #maximum distance allowed between two neighboring alleles
				#in a poplocus-network
#ref to %ind_all:	# {sl}{popall_ID}: needed
					# {seq}=allseq: needed
					# {inds}{ind}= @ of [0] indlocID [1] indall_ID: not needed
#ref to %fno_f_ID: {fno}{frag}=@of popall_IDs
#ref to %pairs: {popall_ID1}{popall_ID2}=dist
sub popall_pairs {
	#declare and initialize: _r means a reference
	my ($sl,$inter_ind_dist,$ind_all_r,$fno_f_ID_r,$pairs_r) = @_;
	my $fno = 0;#current fragment number
	my $frag = '';#current fragment
	my @candidates = ();#candidate popall_IDs for pw distance calculation
	my $ncandidates = 0;#number of candidates
	my $popall_ID1 = 0;#ID of first allele in a pair
	my $popall_ID2 = 0;#ID of second allele in a pair
	my $popall_seq1 = '';#sequence of first allele in a pair
	my $popall_seq2 = '';#sequence of second allele in a pair
	my $dist = 0;#distance (number of differences)
	my $i = 0;
	my $j = 0;
	
	for $fno (keys %{$fno_f_ID_r}) {#loop through fragment-numbers
		for $frag (keys %{$$fno_f_ID_r{$fno}}) {#loop through the fragments
			@candidates = @{$$fno_f_ID_r{$fno}{$frag}};#get the candidates
			$ncandidates = @candidates;#determine their number
			if ($ncandidates > 1) {#if there is more than 1 candidate
				#loop through candidates and form pairs
				for ($i = 0; $i < ($ncandidates - 1); ++$i) {
					for ($j = ($i + 1); $j < $ncandidates; ++$j) {
						#get a pair of candidate numbers
						$popall_ID1 = $candidates[$i];
						$popall_ID2 = $candidates[$j];
						#if this pair has not yet been seen
						unless ((defined $$pairs_r{$popall_ID1}) and (defined $$pairs_r{$popall_ID1}{$popall_ID2})) {
							#get the corresponding sequences
							$popall_seq1 = $$ind_all_r{$sl}{$popall_ID1}{"seq"};
							$popall_seq2 = $$ind_all_r{$sl}{$popall_ID2}{"seq"};
							#call sub dist to determine the distance between them
							$dist = dist($popall_seq1,$popall_seq2,$sl);
							if ($dist <= $inter_ind_dist) {#if distance acceptable
								$$pairs_r{$popall_ID1}{$popall_ID2} = $dist;#add the pair to %pairs
							}
						}
					}
				}				
			}
		}
	}	
}

#definition of subroutine dist
#expects 2 DNA sequences and a scalar that holds their length
#returns number of differing nucleotides between them
#naively treats N as a fifth character, but no N occur here - no problem
sub dist {
	#declare and initialize
	my ($seq1,$seq2,$sl) = @_;
	my $dist = 0;#number of different positions in the 2 seqs
	
	$dist = $sl - (($seq1 ^ $seq2) =~ tr/\0//);
	return $dist;	
}

#definition of subroutine al_poploc
#assembles alleles into loci
#expects:
#$sl: seqlength
#ref to %ind_all:	# {sl}{popall_ID}:
					# {seq}=allseq
					# {inds}{ind}= @ of [0] indlocID [1] indall_ID
#ref to $next_poplocID: first poplocID to use,
#ref to %pairs: {popall_ID1}{popall_ID2}=dist
#ref to %poploc: {sl}{poplocID} = @ of popall_ID
#populates %poploc, owned by main for the current seqlength
#after having determined the highest current poplocID, increments by 1 (as always)
#updates $next_poplocID in main for next set of loci (next seqlength)
sub al_poploc {
	#declare and initialize: _r means a reference
	my ($sl,$ind_all_r,$next_poplocID_r,$pairs_r,$poploc_r) = @_;
	my $popall_ID1 = 0;#ID of first allele in a pair
	my $popall_ID2 = 0;#ID of second allele in a pair
	my $curr_popall_ID = 0;#a current popall_ID
	my $popall_group1 = 0;#current group of first allele in a pair
	my $popall_group2 = 0;#current group of second allele in a pair
	my $groupNo = 1;#a group number
	my $poplocID = $$next_poplocID_r;#first final poplocID to use
	my %group_popall_IDs = ();# {groupNo} = @ of popallIDs
	my %popall_ID_groupNo = ();# {popall_ID} = groupNo
	my @temparr1 = ();
	
	#initialize %popall_ID_groupNo with groupNo 0 for all popall_IDs of this seqlength
	#loop through popall_IDs for this seqlength
	for $popall_ID1 (keys %{$$ind_all_r{$sl}}) {
		$popall_ID_groupNo{$popall_ID1} = 0;
	}
	
	#start assembling loci with those loci that are partners of a pair
	#loop through pairs in %pairs
	for $popall_ID1 (keys %{$pairs_r}) {
		for $popall_ID2 (keys %{$$pairs_r{$popall_ID1}}) {
			#check current groupNo of these alleles
			$popall_group1 = $popall_ID_groupNo{$popall_ID1};
			$popall_group2 = $popall_ID_groupNo{$popall_ID2};
			#if both are not yet in a group
			#make up a new group for them
			if ($popall_group1 == 0 and $popall_group2 == 0) {
				@{$group_popall_IDs{$groupNo}} = ($popall_ID1,$popall_ID2);
				$popall_ID_groupNo{$popall_ID1} = $groupNo;
				$popall_ID_groupNo{$popall_ID2} = $groupNo;
				++$groupNo;
			}
			#if allele 1 already is in a group but allele 2 is not
			elsif ($popall_group1 != 0 and $popall_group2 == 0) {
				#put allele2 into group of allele1
				push @{$group_popall_IDs{$popall_group1}}, $popall_ID2;
				$popall_ID_groupNo{$popall_ID2} = $popall_group1;
			}
			#if allele 2 is not yet in a group but allele 2 is
			elsif ($popall_group1 == 0 and $popall_group2 != 0) {
				#put allele1 into group of allele2
				push @{$group_popall_IDs{$popall_group2}}, $popall_ID1;
				$popall_ID_groupNo{$popall_ID1} = $popall_group2;
			}
			#if both are in different groups
			elsif ($popall_group1 != $popall_group2) {
				#join their groups:
				#get all members of the group allele 2 is in:
				@temparr1 = @{$group_popall_IDs{$popall_group2}};
				#loop through them and change their group-numbers to group of allele1
				for $curr_popall_ID (@temparr1) {
					$popall_ID_groupNo{$curr_popall_ID} = $popall_group1;
				}
				#add them all to the group of allele 1
				push @{$group_popall_IDs{$popall_group1}}, @temparr1;
				#delete the former group of allele 2
				delete $group_popall_IDs{$popall_group2};
			}			
		}
	}
	#Continue with alleles in %popall_ID_groupNo that have not yet a group (groupNo == 0)
	#Make up new, homozygous groups for them
	#loop through %popall_ID_groupNo
	for $popall_ID1 (keys %popall_ID_groupNo) {
		#if the allele has not yet a group
		if ($popall_ID_groupNo{$popall_ID1} == 0) {
			#assign a new group number
			$popall_ID_groupNo{$popall_ID1} = $groupNo;
			#build a corresponding group
			push @{$group_popall_IDs{$groupNo}}, $popall_ID1;
			++$groupNo;#increment for next group
		}
	}
	#all alleles are in a group now
	#populate %poploc: {sl}{poplocID} = @ of popall_ID
	#use new, consecutive poplocIDs instead of the temporary $groupNo
	#loop through %group_popall_IDs
	for $groupNo (keys %group_popall_IDs) {
		#get all members of the group
		@temparr1 = @{$group_popall_IDs{$groupNo}};
		#make up a new locus for them in %poploc
		@{$$poploc_r{$sl}{$poplocID}} = @temparr1;
		++$poplocID;
	}
	#finally, update $next_poplocID in main
	$$next_poplocID_r = $poplocID;
}

#definition of subroutine pop_outfiles
#creates outfiles poploc.txt and popall.txt
#determines output locus by locus and prints to outfiles
#expects:
#ref to %poploc: {sl}{poplocID} = @ of popall_ID
#ref to %ind_all:	# {sl}{popall_ID}:
					# {seq}=allseq
					# {inds}{ind}= @ of [0] indlocID [1] indall_ID
#$msl	length of sequences that shall be merged
#mino	minimum overlap required for merging
#maxo	maximum overlap allowed for merging
#returns $success: Analysis completed or failure message
sub pop_outfiles {
	#declare and initialize: _r means a reference
	my ($poploc_r,$ind_all_r,$msl,$mino,$maxo) = @_;
	my $poplocID = 0;#ID of a locus
	my $popall_ID = 0;#ID of an allele
	my $popall_seq = '';#seq of an allele, may be merged
	my $popall_seq_nomerg = '';#seq of an allele, original, not merged
	my $sl = 0;#seqlength
	my $cons = '';#consensus sequence of a locus
	my $nSNP = 0;#number of SNPs
	my @varpos = ();#variable positions of a locus
					#count starting with 0
	my $varpos = '';#string: variable positions separated by "/"
					#count starting with 1
	my %popall_ID_var = ();# {popall_ID} = $varchar characters at variable positions  in an allele
	my $varchar = '';#characters at variable positions in an allele
	my $n_all = 0;#number of alleles of a locus
	my $nInd = 0;#number of inds that founded a locus
	my $nIndloc = 0;#number of individual loci a locus is based on
	my $ties = 0;#1 if $nIndloc > $nInd
	my @popall_IDs = ();#all popall_IDs of a locus
	my %popall_ID_seq = ();# {popall_ID} = seq; for one locus; these sequences will be merged if appropriate
	my %popall_ID_seq_nomerg = ();# {popall_ID} = seq; unmerged version of each popall
	my $lon_over = 0;#longest possible overlap of any allele of a locus
	my $lon_com_over = 0;#longest overlap possible for all alleles of a locus
	my $merge_conflict = 0;#1: longest possible overlap of one allele is longer
							#than longest possible overlap common to all alleles of a locus, 0: not so
	my $new_sl = 0;#new sequence length after merging						
	my %locfounders = ();# {ind}{indlocID} = number of alleles; for one locus:
						#individuals that contributed, individual loci that contributed, number of their alleles in that ind
	my $ind = '';#an individual
	my @indlocdata = ();#[0] indlocID, [1] indall_ID, for one allele
	my $indlocID = 0;#ID of a locus in an individual
	my @locali = ();#alignment of one locus' alleles, 2d, d1: rows, d2 columns
	my $success = "Analysis aborted.\n";
	my $i = 0;
	my $tempstring1 = '';
	my @temparr1 = ();
	
	#open outfiles and print headerlines
	unless(open(OUTPOPLOCI, ">poploc.txt")) {
		print "Cannot open file poploc.txt. Analysis aborted.\n";
		$success = "Could not open file poploc.txt, Analysis aborted.\n";
		return $success;
	}
	unless(open(OUTPOPALL, ">popall.txt")) {
		print "Cannot open file popall.txt. Analysis aborted.\n";
		$success = "Could not open file popall.txt, Analysis aborted.\n";
		return $success;
	}
	print OUTPOPLOCI "poplocID\tsl\tcons\tnSNP\tvarpos\tn_all\tnInd\tnIndloc\tties\tlon_over\tlon_com_over\tmerge_conflict\n";
	print OUTPOPALL "poplocID\tpopall_ID\tpopall_seq\tpopallvar\tpopall_seq_notmerged\n";
	
	for $sl (sort {$a <=> $b} keys %{$poploc_r}) {#loop through %poploc: key1: sl key2: poplocID
		for $poplocID (sort {$a <=> $b} keys %{$$poploc_r{$sl}}) {
			@popall_IDs = @{$$poploc_r{$sl}{$poplocID}};#get all allele IDs of this locus
			for $popall_ID (@popall_IDs) {#loop through popall_IDs
				#get sequences and store in %popall_ID_seq, %popall_ID_seq_nomerg
				$popall_ID_seq{$popall_ID} = $$ind_all_r{$sl}{$popall_ID}{"seq"};
				$popall_ID_seq_nomerg{$popall_ID} = $$ind_all_r{$sl}{$popall_ID}{"seq"};
				#get individuals and their individual locus IDs for this locus
				for $ind (keys %{$$ind_all_r{$sl}{$popall_ID}{"inds"}}) {
					@indlocdata = @{$$ind_all_r{$sl}{$popall_ID}{"inds"}{$ind}};
					$indlocID = $indlocdata[0];
					++$locfounders{$ind}{$indlocID};
				}
			}
			if ($sl == $msl) {#if locus shall be merged
				#Call sub merge_poploc to merge if possible
				#It merges alleles in %popall_ID_seq, populates $new_sl with new seqlength after merging
				#and stores info about merging conflicts in $merge_conflict
				merge_poploc(\%popall_ID_seq,\$new_sl,\$lon_over,\$lon_com_over,\$merge_conflict,$sl,$mino,$maxo);
			} else {#no merging
				$new_sl = $sl;#new seqlength is old one
				$lon_over = 0;#no overlap possible
				$lon_com_over = 0;#no overlap possible
				$merge_conflict = 0;#no merge conflict
			}
			$n_all = @popall_IDs;#determine $n_all
			$nInd = keys %locfounders;#determine $nInd
			#determine nIndloc
			for $ind (keys %locfounders) {
				$i = keys %{$locfounders{$ind}};
				$nIndloc += $i;
			}
			#determine $ties
			if ($nIndloc > $nInd) {
				$ties = 1;
			}
			if ($n_all == 1) {#if there is only one allele in the locus
				$varpos = 'NA';#no variable positions
				$nSNP = 0;#no SNPs
				#get the single allele's ID
				for $popall_ID (keys %popall_ID_seq) {
					#set the consensus sequence of the locus to the sequence of this allele
					$cons = $popall_ID_seq{$popall_ID};
					#set the characters at variable positions of the single allele  to "consensus"
					$popall_ID_var{$popall_ID} = 'consensus';
				}				
			}
			elsif ($n_all > 1) {#if there is more than 1 allele in the locus
				#build an alignment @locali as 2d array: d1: rows, d2 columns
				#loop through the alleles
				for $popall_ID (sort {$a <=> $b} keys %popall_ID_seq) {
					$tempstring1 = $popall_ID_seq{$popall_ID};#store one allele in a string
					@temparr1 = split('',$tempstring1);#split into array
					push @locali, [@temparr1];#add to @locali
				}
				#call sub varpos to determine variable positions (count starting with 0)
				@varpos = varpos(\@locali);
				#determine variable positions-string for output
				#position-count starting with 1
				@temparr1 = @varpos;#copy variable positions
				for ($i = 0; $i < @temparr1; ++$i) {
					++$temparr1[$i];#increment each position by 1
				}
				$varpos = 'p:';#varpos-string starts with "p:"
				$varpos .= join('/',@temparr1);#join variable positions with "/" as separator and append to varpos-string
				$nSNP = @varpos;#determine $nSNP
				#determine consenus sequence
				$cons = consensus_IUPAC(\@locali,\@varpos);
				#determine characters at variable positions for all alleles
				#store in %popall_ID_var
				varchar(\@varpos,\%popall_ID_seq,\%popall_ID_var);				
			}
			#print output for this locus to poploc.txt
			print OUTPOPLOCI "$poplocID\t$new_sl\t$cons\t$nSNP\t$varpos\t$n_all\t$nInd\t$nIndloc\t$ties\t$lon_over\t$lon_com_over\t$merge_conflict\n";
			#print output for this locus to popall.txt
			#loop through popall_IDs
			for $popall_ID (sort {$a <=> $b} keys %popall_ID_seq) {
				$popall_seq = $popall_ID_seq{$popall_ID};
				$varchar = $popall_ID_var{$popall_ID};
				$popall_seq_nomerg = $popall_ID_seq_nomerg{$popall_ID};
				print OUTPOPALL "$poplocID\t$popall_ID\t$popall_seq\t$varchar\t$popall_seq_nomerg\n";
			}
			
			#set variables back
			%popall_ID_seq = ();
			%locfounders = ();
			$nIndloc = 0;
			$ties = 0;
			@locali = ();
			%popall_ID_var = ();
			#go to next locus
		}
	}
	close OUTPOPLOCI;
	close OUTPOPALL;
	#if you arrived here, give back positive success message
	$success = "analysis completed\n";
	return $success;
}

#definition of subroutine merge_poploc
#Analyzes all alleles of one poploc and merges them (without mismatch) if possible:
#Splits each allele in middle,
#determines longest possible overlap for each allele
#merges with longest overlap possible for all alleles.
#expects references to:
#ref to %popall_ID_seq: {popall_ID}=seq holds the input sequences, which will be replaced by new merged ones
#ref to $new_sl: new sequence length after merging or old sequence length if no merging
#ref to $lon_over: longest possible overlap for any allele of the locus
#ref to $lon_com_over: longest overlap possible for all alleles of the locus
#ref to $merge_conflict: 1 if $lon_over > $lon_com_over, 0 if not
#$sl: sequence length before merging
#$mino: minimum overlap required for merging
#$maxo: maximum overlap allowed for merging
sub merge_poploc {
	#Declare and initialize: _r means a reference
	my ($popall_ID_seq_r,$new_sl_r,$lon_over_r,$lon_com_over_r,
	$merge_conflict_r,$sl,$mino,$maxo) = @_;
	my $n_all = 0;#number of alleles in locus
	my $popall_ID = 0;#ID of one allele
	my $popallseq = '';#sequence of one allele
	my @overlaps_thisall = ();#(lengths of) all possible overlaps for one allele,
							#ascending, first value always 0
	my %overlaps_thisloc = ();# {overlap}=number of alleles with this overlap
	my $overlap = 0;#length of an overlap
	my $lon_over = 0;#longest possible overlap in any allele of this locus
	my $lon_com_over = 0;#longest overlap possible for all alleles of this locus
	my $frag1 = '';#first fragment of seq for merging
	my $frag2 = '';#second fragment of seq for merging
	my $frag1len = 0;#length of first fragment
	my $frag2len = 0;#length of second fragment
	my $startfrag2 = 0;#first position of second fragment
	
	$n_all = keys %{$popall_ID_seq_r};#determine number of alleles
	for $popall_ID (keys %{$popall_ID_seq_r}) {#loop through alleles
		$popallseq = $$popall_ID_seq_r{$popall_ID};#get allele sequence
		#Call sub overlap1 to determine all possible overlaps
		#into @overlaps_thisall in ascending order
		@overlaps_thisall = overlap1($popallseq,$sl,$mino,$maxo);
		for $overlap (@overlaps_thisall) {#loop through found overlaps
			++$overlaps_thisloc{$overlap};#count overlap
		}
	}
	$lon_over = (reverse sort {$a <=> $b} keys %overlaps_thisloc)[0];#get longest overlap for this locus
	#determine longest common overlap of all alleles of this locus 
	for $overlap (reverse sort {$a <=> $b} keys %overlaps_thisloc) {
		if ($overlaps_thisloc{$overlap} == $n_all) {#if all alleles have this overlap
			$lon_com_over = $overlap;
			last;
		}
	}
	$$lon_over_r = $lon_over;#report
	$$lon_com_over_r = $lon_com_over;#report
	if ($lon_over == $lon_com_over) {#no merge-conflict
		$$merge_conflict_r = 0;#report
	} else {#merge-conflict
		$$merge_conflict_r = 1;#report
	}
	if ($lon_com_over > 0) {#if merging possible, do
		$frag1len = $sl / 2;#length of first fragment to join for merging
		$frag2len = $frag1len - $lon_com_over;#length of 2nd frag. to join for merging
		$startfrag2 = $frag1len + $lon_com_over;#startpos of 2nd frag. in unmerged seq
		$$new_sl_r = $frag1len + $frag2len;#report new merged sequence length
		for $popall_ID (keys %{$popall_ID_seq_r}) {#loop through alleles
			$popallseq = $$popall_ID_seq_r{$popall_ID};#get allele sequence
			$frag1 = substr($popallseq,0,$frag1len);
			$frag2 = substr($popallseq,$startfrag2,$frag2len);
			$$popall_ID_seq_r{$popall_ID} = $frag1.$frag2;#replace unmerged seq with merged seq
		}
	} else {#merging not possible
		$$new_sl_r = $sl;#report old sequence length as new one
	}
}

#definition of subroutine overlap1 ver.02
#Expects
#$inseq		DNA sequence
#$sl		length of DNA sequence, must be even
#$mino      minimum overlap
#$maxo      maximum overlap
#
#Splits sequence into two fragments of equal length
#Determines all possible overlaps of the two fragments
#between minimum and maximum overlap
#Returns array with all overlap-lengths in ascending order
#First value in array is always 0 for no overlap
sub overlap1 {
	#declare and initialize
	my ($inseq,$sl,$mino,$maxo) = @_;
	my $fraglen = 0;#length of each fragment
	my $frag1 = '';#fragment1
	my $frag2 = '';#fragment2
	my $startfrag1 = 0;#start position of overlap in $startfrag1
	my $overlap = 0;#current tested overlap between fragments
	my @overlaps = ();#found overlap-lengths in ascending order
	
	#Build fragments
	$fraglen = $sl / 2;
	$frag1 = substr($inseq,0,$fraglen);
	$frag2 = substr($inseq,($sl - $fraglen),$fraglen);
	#Analyze
	$overlaps[0] = 0;#enter 0 as first value
	$startfrag1 = $fraglen - $mino;
	for ($overlap = $mino; $overlap <= $maxo; ++$overlap) {
		#test current overlap
		if (substr($frag1,$startfrag1,$overlap) eq substr($frag2,0,$overlap)) {
			push @overlaps, $overlap;
		}
		--$startfrag1;
	}
	return @overlaps;
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

#definition of subroutine consensus_IUPAC
#expects
#reference to an alignment as 2d-array, all seqs same length
	#allowed characters (all capital: ACGTN)
#reference to an array of variable positions
#returns a consensus-sequence as string
#consensus: variable positions are expressed as IUPAC ambiguity symbol
#no majority-rule
sub consensus_IUPAC {
	#declare and initialize
	my ($ali_r,$varpos_r) = @_;
	my @ali = @{$ali_r};#attention: 2d array, don't modify inner arrays
	my @varpos = @{$varpos_r};#the variable positions in own array
	my %IUPAC_ambi = (#IUPAC-ambiguity symbols
	AG => "R",
	CT => "Y",
	CG => "S",
	AT => "W",
	GT => "K",
	AC => "M",
	CGT => "B",
	AGT => "D",
	ACT => "H",
	ACG => "V",	
	);#what's not defined here, will be N, characters will be sorted before lookup
	my @consarr = ();#consensus sequence as array
	my $cons = '';#final consensus sequence as string
	my $pos = 0;#current position
	my %chars_h = ();#all chars at a position as string
	my $char = '';#a character
	my $chars_s = '';#all chars at a position, sorted alphabetically, as string
	my $row = 0;#a row
	my $col = 0;#a column
	
	@consarr = @{$ali[0]}; #copy first sequence
	#loop through variable positions
	for $pos (@varpos) {
		#loop through rows at that position
		for ($row = 0; $row < @ali; ++$row) {
			$char = ${$ali[$row]}[$pos];#get character
			++$chars_h{$char};#add char to %chars_h
		}
		#loop through sorted characters in %chars_h
		for $char (sort keys %chars_h) {
			$chars_s .= $char;#append to $chars_s
		}
		#if this set of characters is defined in %IUPAC_ambi
		if (defined $IUPAC_ambi{$chars_s}) {
			#replace character in consensus at this pos with ambiguity symbol
			$consarr[$pos] = $IUPAC_ambi{$chars_s};
		} else {#if not: use N
			$consarr[$pos] = 'N';
		}
		#set variables back
		%chars_h = ();
		$chars_s = '';
		#go to next position
	}
	$cons = join('',@consarr);#join consensus into string
	return $cons;
}

#definition of sub varchar
#expects:
#ref to @varpos: variable positions, sorted ascending, count starting with 0
#ref to %ID_seq: {ID}=seq seq: DNA sequence as string
#ref to %ID_varchar: {ID}=varchar varchar: characters at variable positions
#IDs are numbers
#populates %ID_varchar
sub varchar {
	#declare and initialize: _r means a reference
	my ($varpos_r, $ID_seq_r, $ID_varchar_r) = @_;
	my @varpos = @{$varpos_r};#variable positions
	my $seqID = 0;#a sequence ID
	my $seq = '';# a sequence
	my $varchar = '';# its characters at the variable positions
	my $pos = 0;#a position
	my $char = '';#a character
	
	#loop through input sequences
	for $seqID (sort {$a <=> $b} keys %{$ID_seq_r}) {
		$seq = $$ID_seq_r{$seqID};#get a sequence
		#loop through variable positions
		for $pos (@varpos) {
			$char = substr($seq,$pos,1);#get character at that position
			$varchar .= $char;#append to $varchar
		}
		$$ID_varchar_r{$seqID} = $varchar;#store result
		$varchar = '';#set back and next sequence
	}
}
