#!/usr/bin/perl -w
#export_tables version 07.0 Copyright 2015 Andreas Hapke
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
#varpos
#consensus_IUPAC

#Keyword Infilecolumns! Order of columns in an infile is important.

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "export_tables version 07.0 Copyright 2015 Andreas Hapke\n",
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
my $individuals_fname = 'individuals.txt';#out of indloc: list of all individual IDs in database
my $poplocfilename = 'poploc.txt';#poploc-outfile from poploc
my $popallfilename = 'popall.txt';#popall-outfile from poploc
my $indpopallsuffix = '_indpopall.txt';#suffix of indpopall-outfile from indpoploc
my $alleles_suffix = '_alleles.txt';#suffix of alleles-outfile from indloc
#other variables
my %poploc = ();# {poploc_ID}=(sl,cons,nSNP,varpos,n_all,nInd,nIndloc,ties)
my %popall = ();# {poploc_ID}{popall_ID}=(popall_seq,popallvar)
my $popallvar = '';# short form of an allele (chars at variable positions)
my %sel_ind = ();# {ind}=1
my %sel_loc = ();# {poplocID}=1
my %sel_gt = ();# {poplocID}{ind}=1
my $poplocID = 0;#population-locus ID
my $popall_ID = 0;#population-allele ID
my $indlocID = 0;#individual locus ID
my $indall_ID = 0;#individual allele ID
my $ind = '';#individual ID
my $infilename = '';#an infilename
my %indpoploc = ();# {poplocID}{indlocID}{indall_ID}=popall_ID
my %alldep = ();# {indall_ID}=Alldep
my $thisalldep = 0;# an allele depth
my %poplocdatasel = ();
# {poplocID}{"popall_IDs"}{popall_ID} = counter
# {poplocID}{"nInd"} = number of inds at this locus in this selection
# {poplocID}{"nIndloc"} = number of indlocs at this locus in this selection
my %gtreadsout = ();# {poplocID}{ind}{popallvar}=thisalldep
my @gt_ind = ();# one genotype of an individual
my $gt_ind_string = '';#same as string, / as separator
my @alldep_ind = ();# read-depths of one genotype of an individual
my $alldep_ind_string = '';#same as string, / as separator
my @gt_thisloc = ();# all ind-genotypes of one locus for output
my @dep_thisloc = ();# all ind-genotype-depths of one locus for output
my $sl = 0;#sequence length
my $nInd = 0;#number of individuals represented in a locus
my $nIndloc = 0;#number of individual loci represented in a locus
my $ties = 0;#1 if $nIndloc > $nInd, else 0
my $n_all = 0;#number of alleles of a locus
my @locali = ();#2d-array,alignment of a locus' alleles, d1: rows, d2: cols
my $nSNP = 0;#number of SNPs at a locus
my @varpos = ();#variable positions of a locus, count starting with 0
my @varposoutarr = ();#variable positions of a locus, count starting with 1 for output
my $varposout = '';#variable positions of a locus, count starting with 1 for output: p:6/13 or NA if none
my $cons = '';#consensus-sequence
my %countall = ();# for one loc: {popallvar} = count
my @allcounts = ();# for one loc: allele counts in reverse order
my $countmaj = 0;#count of major allele

my $tempstring1 = '';
my @temparr1 = ();
my @temparr2 = ();
my $i = 0;
#}

#########################################
#read in files needed for all individuals
#########################################

#{
print "\nEXPORT TABLES\n\n", "Loading data...\n";

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
	$poplocID = $temparr1[0];#get poplocID
	@{$poploc{$poplocID}} = @temparr1[1 .. 8];#load needed data into %poploc
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
	$poplocID = $temparr1[0];#get poplocID
	$popall_ID = $temparr1[1];#get popallID
	@{$popall{$poplocID}{$popall_ID}} = @temparr1[2 .. 3];#load needed data into %popall

}
close INFILE;

#if directory "export" exists, check for selection files and read them in
if (-d "export") {
	if (-f "export/sel_ind.txt") {#if there is an individual selection-file
		print "Loading your selection of individuals.\n";
		unless(open(INFILE,"export/sel_ind.txt")) {
			print "Cannot open export/sel_ind.txt. Exiting..\n";
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
			"I didn't produce any new file. Exiting..\n";
			exit;			
		}
	}
	if (-f "export/sel_loc.txt") {#if there is a locus selection-file
		print "Loading your selection of loci.\n";
		unless(open(INFILE,"export/sel_loc.txt")) {
			print "Cannot open export/sel_loc.txt. Exiting..\n";
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
			"I didn't produce any new file. Exiting..\n";
			exit;			
		}
	}
	if (-f "export/sel_gt.txt") {#if there is a genotype selection-file
		print "Loading your selection of genotypes.\n";
		unless(open(INFILE,"export/sel_gt.txt")) {
			print "Cannot open export/sel_gt.txt. Exiting..\n";
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
			"I didn't produce any new file. Exiting..\n";
			exit;			
		}
	}	
} else {#if directory "export" doesnt exist, create it
	unless (mkdir 'export') {
		print "Cannot create folder export. Exiting.\n";
		exit;
	}
}
#if %sel_ind is still empty, create it from file individuals.txt (using all inds)
if (keys %sel_ind == 0) {
	print "Exporting data about all individuals in the database.\n";
	unless(open(INFILE,$individuals_fname)) {
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
	print "Exporting data about all loci in the database.\n";
	for $poplocID (keys %poploc) {
		$sel_loc{$poplocID} = 1;
	}
}
#if %sel_gt is still empty, create from %sel_loc and %sel_ind
#(using all genotypes from selected inds and locs)
if (keys %sel_gt == 0) {
	print "Exporting data about all individuals and loci in the database.\n";
	for $poplocID (keys %sel_loc) {
		for $ind (keys %sel_ind) {
			$sel_gt{$poplocID}{$ind} = 1;
		}
	}
}
#}

#######################################
#read in files about single individuals
#populate individual datastructures and
#datastructures about all individuals
#######################################

#{
print "Loading data...\n";
#loop through selected individuals in %sel_ind
 for $ind (keys %sel_ind) {
	#read in file *_indpopall.txt for this ind (outfile from indpoploc)
	my $infilename = $ind.$indpopallsuffix;
	unless(open(INFILE,$infilename)) {
		print "Cannot open $infilename. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip headerline
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		#Infilecolumns!
		@temparr1 = split(/\t/,$tempstring1);
		#store in %indpoploc: {poplocID}{indlocID}{indall_ID}=popall_ID
		$indpoploc{$temparr1[2]}{$temparr1[0]}{$temparr1[1]} = $temparr1[3];
	}
	close INFILE;
	#read in file *_alleles.txt for this ind (outfile from indloc)
	$infilename = $ind.$alleles_suffix;
		unless(open(INFILE,$infilename)) {
		print "Cannot open $infilename. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip headerline
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		#Infilecolumns!
		@temparr1 = split(/\t/,$tempstring1);
		if ($temparr1[1] eq 'valid') {#if this individual locus is valid (data present)
			#store in %alldep: {indall_ID}=Alldep
			$alldep{$temparr1[3]} = $temparr1[5];
		}
	}
	close INFILE;
	#loop through selected poplocID in %sel_loc
	for $poplocID (keys %sel_loc) {
		#if this locus from this individual (this genotype) is selected
		#and if this locus has been genotyped in this individual
		if ((defined $sel_gt{$poplocID}) and (defined $sel_gt{$poplocID}{$ind})
			and (defined $indpoploc{$poplocID})) {
			++$poplocdatasel{$poplocID}{"nInd"};#count ind for this poploc
			#loop through indlocIDs of this ind for this poploc
			for $indlocID (keys %{$indpoploc{$poplocID}}) {
				++$poplocdatasel{$poplocID}{"nIndloc"};#count indloc for this poploc
				#loop through indall_IDs for this indloc
				for $indall_ID (keys %{$indpoploc{$poplocID}{$indlocID}}) {
					#get corresponding popall_ID
					$popall_ID = $indpoploc{$poplocID}{$indlocID}{$indall_ID};
					#store it in %poplocdatasel by counting it
					++$poplocdatasel{$poplocID}{"popall_IDs"}{$popall_ID};
					#look up depth of this allele in %alldep
					$thisalldep = $alldep{$indall_ID};
					#look up popallvar in %popall: {poploc_ID}{popall_ID}=(popall_seq,popallvar)
					$popallvar = $popall{$poplocID}{$popall_ID}[1];
					#store in %gtreadsout: {poplocID}{ind}{popallvar}=thisalldep
					$gtreadsout{$poplocID}{$ind}{$popallvar} = $thisalldep;
				}
			}
		}
	}
	#set variables back
	%indpoploc = ();
	%alldep = ();
	#next ind
}
#}

##########################################
#print data to genotypes.txt and reads.txt
##########################################

#{
#open outfiles and print headerlines
unless(open(OUTGT, ">export/genotypes.txt")) {
	print "Cannot open file export/genotypes.txt, exiting ...\n\n";
	exit;
}
unless(open(OUTRD, ">export/reads.txt")) {
	print "Cannot open file export/reads.txt, exiting ...\n\n";
	exit;
}
print "Printing to files export/genotypes.txt and export/reads.txt...\n";
#print headerlines
print OUTGT "poplocID";
print OUTRD "poplocID";
for $ind (sort keys %sel_ind) {
	print OUTGT "\t$ind";
	print OUTRD "\t$ind";
}
print OUTGT "\n";
print OUTRD "\n";
#print data
#loop through poplocIDs (sorted) in %gtreadsout: {poplocID}{ind}{popallvar}=thisalldep
for $poplocID (sort {$a <=> $b} keys %gtreadsout) {
	print OUTGT "$poplocID\t";
	print OUTRD "$poplocID\t";
	#loop through selected individuals in %sel_ind
	for $ind (sort keys %sel_ind) {
		#if there are data for this locus and this ind in %gtreadsout
		if (defined $gtreadsout{$poplocID}{$ind}) {
			#loop through alleles (popallvar) in %gtreadsout
			for $popallvar (sort keys %{$gtreadsout{$poplocID}{$ind}}) {
				#get depth of this allele
				$thisalldep = $gtreadsout{$poplocID}{$ind}{$popallvar};
				#add popallvar and thisalldep to arrays for this individual genotype
				push @gt_ind, $popallvar;
				push @alldep_ind, $thisalldep;			
			}
			#join alleles and read-depths of this individual genotype into strings
			$gt_ind_string = join('/',@gt_ind);
			$alldep_ind_string = join('/',@alldep_ind);
		} else {#if there are NO data for this locus and this ind in %gtreadsout
			$gt_ind_string = '-999';#missing genotype: -999
			$alldep_ind_string = '-999';#missing genotype: -999
		}
		#add the strings to an array for the whole locus
		push @gt_thisloc, $gt_ind_string;
		push @dep_thisloc, $alldep_ind_string;
		#set variables back
		@gt_ind = ();
		@alldep_ind = ();
		#next ind
	}
	#join all genotypes for this loc into string and print
	$tempstring1 = join("\t",@gt_thisloc);
	print OUTGT "$tempstring1\n";
	#join all reads for this loc into string and print
	$tempstring1 = join("\t",@dep_thisloc);
	print OUTRD "$tempstring1\n";
	#set variables back
	@gt_thisloc = ();
	@dep_thisloc = ();
	#next poploc
}
close OUTGT;
close OUTRD;
print "Saved genotypes to file export/genotypes.txt\n";
print "Saved readcounts to file export/reads.txt\n";
#}

#########################################################################
#print data about selected loci based on complete database to locdata.txt
#########################################################################

#{
#open outfile and print headerline
unless(open(OUTFILE, ">export/locdata.txt")) {
	print "Cannot open file export/locdata.txt, exiting ...\n\n";
	exit;
}
print "Printing to file export/locdata.txt...\n";
print OUTFILE "poplocID\tsl\tcons\tnSNP\tvarpos\tn_all\tnInd\tnIndloc\tties\n";
#loop through poplocIDs in %sel_loc
for $poplocID (sort {$a <=> $b} keys %sel_loc) {
	@temparr1 = @{$poploc{$poplocID}};#get value-array for this locus out of %poploc
	$tempstring1 = join("\t",@temparr1);
	print OUTFILE "$poplocID\t$tempstring1\n";
}
close OUTFILE;
print "Saved data about exported loci based on whole database to file export/locdata.txt\n";
#}

##########################################################################################
#determine data about selected loci based on current selection and print to locdatasel.txt
##########################################################################################

#{
#open outfile and print headerline
unless(open(OUTFILE, ">export/locdatasel.txt")) {
	print "Cannot open file export/locdatasel.txt, exiting ...\n\n";
	exit;
}
print "Calculating and printing to file export/locdatasel.txt...\n";
print OUTFILE "poplocID\tsl\tcons\tnSNP\tvarpos\tn_all\tnInd\tnIndloc\tties\tcount_maj\n";
#loop through poplocIDs in %sel_loc
for $poplocID (sort {$a <=> $b} keys %sel_loc) {
	$sl = $poploc{$poplocID}[0];#look up sequence-length
	$nInd = $poplocdatasel{$poplocID}{"nInd"};#number of inds represented in locus
	$nIndloc = $poplocdatasel{$poplocID}{"nIndloc"};#number of indloc represented in locus
	if ($nIndloc > $nInd) {
		$ties = 1;#poplocus combines 2 or more individual loci in at least one ind
	} else {
		$ties = 0;#not so
	}
	#get array of all popall_IDs for the locus
	@temparr1 = ();
	for $popall_ID (sort {$a <=> $b} keys %{$poplocdatasel{$poplocID}{"popall_IDs"}}) {
		push @temparr1, $popall_ID;
	}
	$n_all = @temparr1;#number of alleles
	if ($n_all == 1) {#if there is only one allele
		$nSNP = 0;#number of SNPs at the locus is 0
		$varposout = 'NA';#no variable positions
		#sequence of the single allele is consensus-sequence of the locus
		$cons = $popall{$poplocID}{$temparr1[0]}[0];
	}
	elsif ($n_all > 1) {#if there is more than one allele
		#build an alignment of the alleles as 2d-array @locali
		#loop through popall_IDs
		for $popall_ID (@temparr1) {
			$tempstring1 = $popall{$poplocID}{$popall_ID}[0];#look up sequence
			@temparr2 = split('',$tempstring1);
			push @locali, [@temparr2];
		}
		#call sub varpos to determine variable positions (count starting with 0)
		@varpos = varpos(\@locali);
		$nSNP = @varpos;#number of SNPs
		#determine variable positions-string for output
		#position-count starting with 1
		@temparr1 = @varpos;#copy variable positions
		for ($i = 0; $i < @temparr1; ++$i) {
			++$temparr1[$i];#increment each position by 1
		}
		$varposout = 'p:';#varpos-string starts with "p:"
		$varposout .= join('/',@temparr1);#join variable positions with "/" as separator and append to varpos-string
		#determine consenus sequence
		$cons = consensus_IUPAC(\@locali,\@varpos);
		
	}
	#determine count of major allele
	for $ind (keys %{$gtreadsout{$poplocID}}) {
		#if the individual is homozygous
		if (keys %{$gtreadsout{$poplocID}{$ind}} == 1) {
			for $popallvar (keys %{$gtreadsout{$poplocID}{$ind}}) {
				$countall{$popallvar} += 2;
			}
		} else {#if it has more than one allele
			for $popallvar (keys %{$gtreadsout{$poplocID}{$ind}}) {
				++$countall{$popallvar};
			}
		}
	}
	@allcounts = reverse sort {$a <=> $b} values %countall;
	%countall = ();#set back
	$countmaj = $allcounts[0];
	#print data to outfile
	print OUTFILE "$poplocID\t$sl\t$cons\t$nSNP\t$varposout\t$n_all\t$nInd\t$nIndloc\t$ties\t$countmaj\n";
	@locali = ();#set back and next locus
}
close OUTFILE;
print "Saved data about exported loci based on selection to file export/locdatasel.txt\n";
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "run took $run_s seconds.\n";

exit;

############
#Subroutines
############

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
