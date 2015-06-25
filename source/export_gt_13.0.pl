#!/usr/bin/perl -w
#export_gt version 13.0 Copyright 2015 Andreas Hapke
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

#Keyword Infilecolumns! Order of columns in an infile is important.

#Copyright message to screen
print "export_gt version 13.0 Copyright 2015 Andreas Hapke\n",
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
uS => '',#0: export alleles with potentially several SNPs
		 #2: export unlinked SNPs, include only loci with 2 alleles (2 chars at 1st variable pos.)
		 #3: export unlinked SNPs, include all variable loci (2 or more chars at 1st variable pos.)
bin => '',#0: print DNA characters to FASTA, Nexus, PHYLIP
          #1: print integer data to FASTA, Nexus, PHYLIP, coded as 012 in one sequence per ind: 0 homo A, 1 het, 2 homo B
		  #2: print binary data to FASTA, Nexus, PHYLIP as 01 in two sequences per ind
S => '',#Boolean: produce outfile for STRUCTURE
Sc => '',#name of text file with pre-genotype columns for STRUCTURE infile
Sm => '',#value for missing datapoint in STRUCTURE infile
G => '',#Boolean: produce outfile for Genepop
Gp => '',#name of text file with population assignments for Genepop
A => '',#Boolean: produce Arlequin project file
Ap => '',#name of text file with population assignments for Arlequin
Fv => '',#Boolean: produce fasta outfile with variable positions
Fmis => '',#missing data value for fasta output
Nxv => '',#Boolean: produce nexus outfile with variable positions
Nxgap => '',#missing data value for nexus output, will be coded as gap character
Pv => '',#Boolean: produce PHYLIP outfile with variable positions
Pmis => '',#missing data value for PHYLIP output
indf => '',#name of file with strings to be appended to indIDs in genepop, FASTA, Nexus, PHYLIP
rs => ''#random number seed
);
my %defaults = (
uS => '0',#0: export alleles with potentially several SNPs
		  #2: export unlinked SNPs, include only loci with 2 alleles (2 chars at 1st variable pos.)
		  #3: export unlinked SNPs, include all variable loci (2 or more chars at 1st variable pos.)
bin => '0',#0: print DNA characters to FASTA, Nexus, PHYLIP
          #1: print integer data to FASTA, Nexus, PHYLIP, coded as 012 in one sequence per ind: 0 homo A, 1 het, 2 homo B
		  #2: print binary data to FASTA, Nexus, PHYLIP as 01 in two sequences per ind
S => '1',#Boolean: produce outfile for STRUCTURE
Sc => '',#name of text file with pre-genotype columns for STRUCTURE infile
Sm => '-9',#value for missing datapoint in STRUCTURE infile
G => '1',#Boolean: produce outfile for Genepop
Gp => '',#name of text file with population assignments for Genepop
A => '1',#Boolean: produce Arlequin project file
Ap => '',#name of text file with population assignments for Arlequin
Fv => '1',#Boolean: produce fasta outfile with variable positions
Fmis => 'N',#missing data value for fasta output
Nxv => '1',#Boolean: produce nexus outfile with variable positions
Nxgap => '-',#missing data value for nexus output, will be coded as gap character
Pv => '1',#Boolean: produce PHYLIP outfile with variable positions
Pmis => 'N',#missing data value for PHYLIP output
indf => '',#name of file with strings to be appended to indIDs in genepop, FASTA, Nexus, PHYLIP
rs => ''#random number seed
);
my $n_arg = 0;#number of arguments provided by user
my $flag = '';#command flag
my $val = '';#value
my $ncols1 =0;#number of columns in first (not empty) line in a user file
my $ncols = 0;#number of columns in a line in user files
#program settings
my $settings_fname = 'export_gt_settings.txt';#outfile with settings and random number seed
my $individuals_fname = 'individuals.txt';#out of indloc: list of all individual IDs in database
my $popallfilename = 'popall.txt';#poploc-outfile from poploc
my $indpopallsuffix = '_indpopall.txt';#suffix of indpopall-outfile from indpoploc
my $indpoplocsuffix = '_indpoploc.txt';#suffix of indpoploc-outfile from indpoploc
my $indfname = '';#infilename specified by flag -indf
my $outdirname = 'export/';#name of directory for outfiles
my $struc_outname = 'in_structure.txt';#name of outfile for STRUCTURE
my $genepop_outname = 'in_genepop.txt';#name of outfile for Genepop
my $arl_outname = 'in_arlequin.arp';#name of outfile for Arlequin
my $genepopmis = 0;#missing datapoint for Genepop
my $fasvarfas_outname = 'varfas.fas';#name of fasta outfile with variable positions
my $varnx_outname = 'varnx.nex';#name of nexus outfile with variable positions
my $Ph_outname = 'varph.phy';#name of PHYLIP outfile with variable positions
my $varpos_outname = 'varpos.txt';#name of locus positions outfile for fasta variable pos
#other variables
my $seed = '';
my $infilename = '';
my %struc_cols = ();# {ind} = string: tab-separated string with pre-genotype contents for STRUCTURE format
my %struc_order = ();# {struc_orderNo} = ind; order of individuals in user-file for STRUCTURE
my $struc_orderNo = 0;#order number for individuals for STRUCTURE outfile
my %poplocallren = ();# {poplocID}{popallID}=popallren: popall renumbered from 1 for each locus
my %poplocallbin = ();# {poplocID}{popallID}= 0/1 for rarer/more frequent character at a biallelic SNP
my %poplocallvar = ();# {poplocID}{popallID}=popallvar
my $popallvar = '';#characters at variable positions in one allele (variable given all known alleles)
my @popallvarthisloc = ();#all popallvars of one locus in popall_ID order
my $pos = 0;#a position
my $char = '';#a character
my %varchar = ();# {character} = count of character
my %charnum = (
A => '1',
C => '2',
G => '3',
T => '4',
);
my %charnumbin = ();#contains charnums 0 and 1 for the two alleles of a biallelic SNP, 1 for the more frequent one
my $poplocID = 0;
my $n_invarloc = 0;#number of invariable loci in selection
my $n_multSNP = 0;#number of SNPs with more than 2 alleles that have been excluded
my $popall_ID = 0;
my $popallren = 1;#renumbered popall
my $ind = '';#individual ID
my %sel_ind = ();# {ind}=1
my %sel_loc = ();# {poplocID}=1
my %sel_gt = ();# {poplocID}{ind}=1
my %popgenepop = ();# {genepop}{ind}=1 populations for Genepop
my $genepop = 0;#population for Genepop
my %arlpops = ();# {arlpop}{ind}=1 populations for Arlequin
my $arlpop = 'Sample1';#population for Arlequin
my %indIDs = ();# {ind}=indID for Genepop, FASTA, Nexus, PHYLIP
my $indID = '';#indID for Genepop, FASTA, Nexus, PHYLIP
my $indID1 = '';#indID for seq 1 of an ind in FASTA, Nexus, PHYLIP
my $indID2 = '';#indID for seq 2 of an ind in FASTA, Nexus, PHYLIP
my %ali = ();# {indID}=seq; IDs and sequence data for FASTA, Nexus, PHYLIP
my $seq = '';#sequence
my $alimis = '87';#missing data character dummy in %ali
my $arl_nsamp = 0;#number of populations for Arlequin
my $lenidarl = 0;#length of indID field for Arlequin output
my $arlsamsize = 0;#size of a population for Arlequin
my %indpopall = ();# {ind}{poplocID}{popall_ID} = 1 for one individual
my %ex_gt = ();# {ind}{poplocID}=(popallren1,popallren2) note: all genotypes must be diallelic
my $strucmis = 0;#missing datapoint for STRUCTURE, changes later
my $start = 0;#start pos of a locus
my $end = 0;#endpos of a locus
my $loclen = 0;#length of a locus
my $ntax = 0;#number of sequences in Nexus and PHYLIP outfile
my $lonID = '';#longest ID of an individual
my $lonIDlen = '';#length of longest ID of an individual
my $PhminIDlen = 10;#minimum ID length for PHYLIP (IDs are filled up with blanks to that length
my $Pspacer = '   ';#PHYLIP outfile: spacer between ID and seqs
my $tempstring1 = '';
my $tempstring2 = '';
my @temparr1 = ();
my %temphash1 = ();
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
#-uS must be 0,2 or 3
if (length $user_settings{'uS'} > 0) {
	unless (($user_settings{'uS'} eq '0') or ($user_settings{'uS'} eq '2') or ($user_settings{'uS'} eq '3')) {
		print "Implausible setting: -uS, using default: $defaults{'uS'}.\n";
		$user_settings{'uS'} = $defaults{'uS'};
	}
} else {
	$user_settings{'uS'} = $defaults{'uS'};
}
#-bin must be 0, 1 or 2, if bin 1 or 2: automatically activate uS 2
if (length $user_settings{'bin'} > 0) {
	unless (($user_settings{'bin'} eq '0') or ($user_settings{'bin'} eq '1') or ($user_settings{'bin'} eq '2')) {
		print "Implausible setting -bin, using default: $defaults{'bin'}.\n";
		$user_settings{'bin'} = $defaults{'bin'};
	}
	if (($user_settings{'bin'} > 0) and ($user_settings{'uS'} ne '2')) {
		print "Setting -bin 1 or 2 is active, I must change -uS $user_settings{'uS'} to -uS 2\n";
		$user_settings{'uS'} = 2;
	}
} else {
	$user_settings{'bin'} = $defaults{'bin'};
}
#-S must be 0 or 1
if (length $user_settings{'S'} > 0) {
	unless (($user_settings{'S'} eq '0') or ($user_settings{'S'} eq '1')) {
		print "Implausible setting: -S, using default: $defaults{'S'}.\n";
		$user_settings{'S'} = $defaults{'S'};
	}
} else {
	$user_settings{'S'} = $defaults{'S'};
}
#-Sc: can you open the file?
if (length $user_settings{'Sc'} > 0) {
	$tempstring1 = $user_settings{'Sc'};
	unless(open(INFILE,$tempstring1)) {
		print "Cannot open file provided under -Sc: $tempstring1.\n",
		"I will ignore it.\n";
		$user_settings{'Sc'} = $defaults{'Sc'};
	}
	close INFILE;
} else {
	$user_settings{'Sc'} = $defaults{'Sc'};
}

#I don't test -Sm provide by user It depends on STRUCTURE if the choice of the user works.
unless (length $user_settings{'Sm'} > 0) {
	$user_settings{'Sm'} = $defaults{'Sm'};
}

#-G must be 0 or 1
if (length $user_settings{'G'} > 0) {
	unless (($user_settings{'G'} eq '0') or ($user_settings{'G'} eq '1')) {
		print "Implausible setting: -G, using default: $defaults{'G'}.\n";
		$user_settings{'G'} = $defaults{'G'};
	}
} else {
	$user_settings{'G'} = $defaults{'G'};
}
#-Gp: can you open the file?
if (length $user_settings{'Gp'} > 0) {
	$tempstring1 = $user_settings{'Gp'};
	unless(open(INFILE,$tempstring1)) {
		print "Cannot open file provided under -Gp: $tempstring1.\n",
		"I will ignore it.\n";
		$user_settings{'Gp'} = $defaults{'Gp'};
	}
	close INFILE;
} else {
	$user_settings{'Gp'} = $defaults{'Gp'};
}
#-A must be 0 or 1
if (length $user_settings{'A'} > 0) {
	unless (($user_settings{'A'} eq '0') or ($user_settings{'A'} eq '1')) {
		print "Implausible setting: -A, using default: $defaults{'A'}.\n";
		$user_settings{'A'} = $defaults{'A'};
	}
} else {
	$user_settings{'A'} = $defaults{'A'};
}
#-Ap: can you open the file?
if (length $user_settings{'Ap'} > 0) {
	$tempstring1 = $user_settings{'Ap'};
	unless(open(INFILE,$tempstring1)) {
		print "Cannot open file provided under -Ap: $tempstring1.\n",
		"I will ignore it.\n";
		$user_settings{'Ap'} = $defaults{'Ap'};
	}
	close INFILE;
} else {
	$user_settings{'Ap'} = $defaults{'Ap'};
}
#-Fv must be 0 or 1
if (length $user_settings{'Fv'} > 0) {
	unless (($user_settings{'Fv'} eq '0') or ($user_settings{'Fv'} eq '1')) {
		print "Implausible setting: -Fv, using default: $defaults{'Fv'}.\n";
		$user_settings{'Fv'} = $defaults{'Fv'};
	}
} else {
	$user_settings{'Fv'} = $defaults{'Fv'};
}
#-Fmis must have length 1
if (length $user_settings{'Fmis'} > 0) {
	unless (length $user_settings{'Fmis'} == 1) {
		print "-Fmis: $user_settings{'Fmis'} is unusable, using default: $defaults{'Fmis'}.\n";
		$user_settings{'Fmis'} = $defaults{'Fmis'};
	}
} else {
	$user_settings{'Fmis'} = $defaults{'Fmis'};
}
#-Nxv must be 0 or 1
if (length $user_settings{'Nxv'} > 0) {
	unless (($user_settings{'Nxv'} eq '0') or ($user_settings{'Nxv'} eq '1')) {
		print "Implausible setting: -Nxv, using default: $defaults{'Nxv'}.\n";
		$user_settings{'Nxv'} = $defaults{'Nxv'};
	}
} else {
	$user_settings{'Nxv'} = $defaults{'Nxv'};
}
#-Nxgap must have length 1
if (length $user_settings{'Nxgap'} > 0) {
	unless (length $user_settings{'Nxgap'} == 1) {
		print "-Nxgap: $user_settings{'Nxgap'} is unusable, using default: $defaults{'Nxgap'}.\n";
		$user_settings{'Nxgap'} = $defaults{'Nxgap'};
	}
} else {
	$user_settings{'Nxgap'} = $defaults{'Nxgap'};
}
#-Pv must be 0 or 1
if (length $user_settings{'Pv'} > 0) {
	unless (($user_settings{'Pv'} eq '0') or ($user_settings{'Pv'} eq '1')) {
		print "Implausible setting: -Pv, using default: $defaults{'Pv'}.\n";
		$user_settings{'Pv'} = $defaults{'Pv'};
	}
} else {
	$user_settings{'Pv'} = $defaults{'Pv'};
}
#-Pmis must have length 1
if (length $user_settings{'Pmis'} > 0) {
	unless (length $user_settings{'Pmis'} == 1) {
		print "-Pmis: $user_settings{'Pmis'} is unusable, using default: $defaults{'Pmis'}.\n";
		$user_settings{'Pmis'} = $defaults{'Pmis'};
	}
} else {
	$user_settings{'Pmis'} = $defaults{'Pmis'};
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

###################################################
#Print used settings to file export_gt_settings.txt
###################################################

#{
$settings_fname = $outdirname . $settings_fname;
unless(open(SETOUT, ">", $settings_fname)) {
	print "Cannot open file $settings_fname, exiting..\n";
	exit;
}
print SETOUT
"export_gt used these settings:\n",
"-uS $user_settings{'uS'}\n",
"-bin $user_settings{'bin'}\n",
"-S $user_settings{'S'}\n",
"-Sc $user_settings{'Sc'}\n",
"-Sm $user_settings{'Sm'}\n",
"-G $user_settings{'G'}\n",
"-Gp $user_settings{'Gp'}\n",
"-A $user_settings{'A'}\n",
"-Ap $user_settings{'Ap'}\n",
"-Fv $user_settings{'Fv'}\n",
"-Fmis $user_settings{'Fmis'}\n",
"-Nxv $user_settings{'Nxv'}\n",
"-Nxgap $user_settings{'Nxgap'}\n",
"-Pv $user_settings{'Pv'}\n",
"-Pmis $user_settings{'Pmis'}\n",
"-indf $user_settings{'indf'}\n",
"-rs $seed\n\n",
"You can use these settings to reproduce the same outfiles.\n",
"Some flags above may not have values.\n",
"These flags are undefined per default and you did not define them.\n",
"Do not use these flags when you call the program to reproduce your outfiles.\n",
"The random number seed -rs is necessary to reproduce\n",
"the random phasing of alleles in FASTA, Nexus, and PHYLIP files.\n",
"Please make sure to use the same selection of data with data_selector\n",
"to reproduce the same outfiles.\n";
close SETOUT;
#}

###################################################################
#Read in popall.txt, populate %poplocallren with 0 as allele number
#and %poplocallvar
###################################################################

#{
print "Loading data...\n";
unless(open(INFILE,$popallfilename)) {
	print "Cannot open $popallfilename. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	#Infilecolumns!
	@temparr1 = split(/\t/,$tempstring1);
	$poplocallren{$temparr1[0]}{$temparr1[1]} = 0;#store poplocID and popall_ID with value 0
	$poplocallvar{$temparr1[0]}{$temparr1[1]} = $temparr1[3];#store poplocID and popall_ID with value popallvar
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
			"I didn't produce any new file. Exiting..\n";
			exit;			
		}
		print "Loading your selection of loci...\n";
		unless(open(INFILE,"export/sel_loc.txt")) {
			print "Cannot open export/sel_loc.txt. Exiting..\n";
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
			"I didn't produce any new file. Exiting..\n";
			exit;			
		}
		print "Loading your selection of genotypes...\n";
		unless(open(INFILE,"export/sel_gt.txt")) {
			print "Cannot open export/sel_gt.txt. Exiting..\n";
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
			"I didn't produce any new file. Exiting..\n";
			exit;			
		}
	}
	else {#if not all selection files are there
		if	(-f "export/sel_ind.txt") {
			print "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_ind.txt later. I will ignore it.\n";
		}
		if	(-f "export/sel_loc.txt") {
			print "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_loc.txt later. I will ignore it.\n";
		}
		if	(-f "export/sel_gt.txt") {
			print "Selection files in directory export are incomplete.\n",
			"Please delete file export/sel_gt.txt later. I will ignore it.\n";
		}
		print "I will export all genotypes in the database.\n";
	}
} else {#if directory "export" doesn't exist, create it
	unless (mkdir 'export') {
		print "Cannot create folder export. Exiting.\n";
		exit;
	}
	print "I will export all genotypes in the database.\n";
}
#if %sel_ind is still empty, create it from file individuals.txt (using all inds)
if (keys %sel_ind == 0) {
	unless(open(INFILE,$individuals_fname)) {
		print "Cannot open $individuals_fname. Exiting..\n";
		exit;
	}
	@temparr1 = <INFILE>;
	close INFILE;
	for $ind (@temparr1) {
		chomp $ind;
		$sel_ind{$ind} = 1;#store in %sel_ind
	}
}
#if %sel_loc is still empty, create from %poplocallren (using all loci)
if (keys %sel_loc == 0) {
	for $poplocID (keys %poplocallren) {
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

#############################################################
#If defined, read in file with additional strings
#to append to indIDs in Genepop FASTA, Nexus and PHYLIP files
#if or if not, initialize %indIDs
#############################################################

#{
#initialize %indIDs with individual IDs in database
for $ind (keys %sel_ind) {
	$indIDs{$ind} = $ind;
}

if (length $user_settings{'indf'} > 0) {
	$indfname = $user_settings{'indf'};
	unless(open(INFILE, "<", $indfname)) {
		print "Cannot open $indfname, exiting..\n";
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
		} else {
			if (defined $indIDs{$temparr1[0]}) {
				$tempstring1 = join("_",@temparr1);
				$indIDs{$temparr1[0]} = $tempstring1;
			} else {
				print "File $indfname: ind $temparr1[0] is not in current selection.\n",
				"I ignore it.\n";
			}
		}
	}
	close INFILE;
}
#}

###################################################################
#If defined, read in user file with pre-genotype cols for STRUCTURE
#Check contents, if or if not: populate %struc_cols: {ind} = string
###################################################################

#{
if (length $user_settings{'Sc'} > 0) {
	$tempstring1 = $user_settings{'Sc'};
	unless(open(INFILE, $tempstring1)) {
		print "Cannot open $tempstring1. Exiting..\n";
		exit;
	}
	#read in first lines including first non-empty one
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		$tempstring1 =~ s/^s*//;#remove leading whitespace
		if (length $tempstring1 > 0) {#if the line is not empty
			@temparr1 = split(/\t/,$tempstring1);
			$ncols1 = @temparr1;#get number of columns in first line
			$struc_cols{$temparr1[0]} = $tempstring1;#store in %struc_cols
			$struc_order{$struc_orderNo} = $temparr1[0];#store ind in %struc_order
			++$struc_orderNo;
			last;
		}	
	}
	#read in remaining lines
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		$tempstring1 =~ s/^s*//;#remove leading whitespace
		if (length $tempstring1 > 0) {#if the line is not empty
			@temparr1 = split(/\t/,$tempstring1);
			$ncols = @temparr1;
			unless($ncols == $ncols1) {#if the number of columns deviates from first line: exit
				print "File $user_settings{'Sc'} has invalid format. Exiting..\n";
				exit;
			}
			$struc_cols{$temparr1[0]} = $tempstring1;#store in %struc_cols
			$struc_order{$struc_orderNo} = $temparr1[0];#store ind in %struc_order
			++$struc_orderNo;
		}
	}
	close INFILE;
	#cross-check individuals in selection (%sel_ind) with individual in %struc_cols
	for $ind (keys %sel_ind) {
		unless (defined $struc_cols{$ind}) {
			print "Individual $ind does not appear in file $user_settings{'Sc'}. Exiting..\n";
			exit;
		}
	}
	for $ind (keys %struc_cols) {
		unless (defined $sel_ind{$ind}) {
			print "Individual $ind in file $user_settings{'Sc'} is not in current selection. Exiting..\n";
			exit;
		}
	}
} else {#if user did not provide a file with pre-genotype-columns for STRUCTURE
	for $ind (sort keys %sel_ind) {
		$struc_cols{$ind} = $ind;#populate %struc_cols with individual IDs only
		$struc_order{$struc_orderNo} = $ind;#store ind in %struc_order
		++$struc_orderNo;
	}
}
#}

################################################################################
#If defined, read in user file with population assignments for Genepop and check
#if or if not: populate %popgenepop {popgenepop}{ind}=1
################################################################################

#{
if (length $user_settings{'Gp'} > 0) {
	$tempstring1 = $user_settings{'Gp'};
	unless(open(INFILE, $tempstring1)) {
		print "Cannot open $tempstring1. Exiting..\n";
		exit;
	}
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		$tempstring1 =~ s/^s*//;#remove leading whitespace
		if (length $tempstring1 > 0) {#if the line is not empty
			@temparr1 = split(/\t/,$tempstring1);
			#check if the line contains two entries, the second one an integer
			unless((@temparr1 == 2) and ($temparr1[1] =~/^\d+$/)) {
				print "File $user_settings{'Gp'}:",
				"I don't understand this line: $tempstring1. Exiting..\n";
				exit;
			}
			$popgenepop{$temparr1[1]}{$temparr1[0]} = 1;
		}
	}
	close INFILE;
	#cross-check individuals in user file with selected individuals in %sel_ind
	for $genepop (keys %popgenepop) {
		for $ind (keys %{$popgenepop{$genepop}}) {
			$temphash1{$ind} = 1;
			unless (defined $sel_ind{$ind}) {
				print "Individual $ind in file $user_settings{'Gp'} is not in current selection. Exiting..\n";
				exit;
			}
		}
	}
	for $ind (keys %sel_ind) {
		unless (defined $temphash1{$ind}) {
			print "Individual $ind is missing in file $user_settings{'Gp'}. Exiting..\n";
			exit;
		}
	}
	%temphash1 = ();
} else {#if user did not provide a file with population assignments for genepop
	for $ind (keys %sel_ind) {
		$popgenepop{'0'}{$ind} = 1;#put all individuals into one pop for Genepop
	}
}
#}

#################################################################################
#If defined, read in user file with population assignments for Arlequin and check
#if or if not: populate %arlpops {arlpop}{ind}=1
#################################################################################

#{
if (length $user_settings{'Ap'} > 0) {
	$tempstring1 = $user_settings{'Ap'};
	unless(open(INFILE, $tempstring1)) {
		print "Cannot open $tempstring1. Exiting..\n";
		exit;
	}
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		$tempstring1 =~ s/^s*//;#remove leading whitespace
		if (length $tempstring1 > 0) {#if the line is not empty
			@temparr1 = split(/\t/,$tempstring1);
			#check if the line contains two entries
			unless(@temparr1 == 2) {
				print "File $user_settings{'Ap'}:",
				"I don't understand this line: $tempstring1. Exiting..\n";
				exit;
			}
			$arlpops{$temparr1[1]}{$temparr1[0]} = 1;
		}
	}
	close INFILE;
	#cross-check individuals in user file with selected individuals in %sel_ind
	%temphash1 = ();
	for $arlpop (keys %arlpops) {
		for $ind (keys %{$arlpops{$arlpop}}) {
			$temphash1{$ind} = 1;
			unless (defined $sel_ind{$ind}) {
				print "Individual $ind in file $user_settings{'Ap'} is not in current selection. Exiting..\n";
				exit;
			}
		}
	}
	for $ind (keys %sel_ind) {
		unless (defined $temphash1{$ind}) {
			print "Individual $ind is missing in file $user_settings{'Ap'}. Exiting..\n";
			exit;
		}
	}
	%temphash1 = ();
} else {#if user did not provide a file with population assignments for Arlequin
	for $ind (keys %sel_ind) {
		$arlpops{$arlpop}{$ind} = 1;#put all individuals into one pop for Arlequin
	}
}
#}

##################################################################################
#Read in data from individuals, populate %indpopall {ind}{poplocID}{popall_ID} = 1
##################################################################################

#{
print "Loading genotype data from individuals...\n";
for $ind (keys %sel_ind) {#loop through selected individuals
	#read in file *_indpopall.txt for this ind (outfile from indpoploc)
	$infilename = $ind.$indpopallsuffix;
	unless(open(INFILE,$infilename)) {
		print "Cannot open $infilename. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip headerline
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		if ((exists $sel_gt{$temparr1[2]}) and (exists $sel_gt{$temparr1[2]}{$ind})) {#if poploc selected for this ind
			$indpopall{$ind}{$temparr1[2]}{$temparr1[3]} = 1;#store in %indpopall
		}
	}
	close INFILE;
}
#}

########################################################################
#Remove unselected loci and alleles from %poplocallren and %poplocallvar
########################################################################

#{
for $ind (keys %indpopall) {
	for $poplocID (keys %{$indpopall{$ind}}) {
		for $popall_ID (keys %{$indpopall{$ind}{$poplocID}}) {
			$poplocallren{$poplocID}{$popall_ID} = 1;
		}
	}
}
for $poplocID (keys %poplocallren) {
	for $popall_ID (keys %{$poplocallren{$poplocID}}) {
		if ($poplocallren{$poplocID}{$popall_ID} == 0) {#if this allele did not occur in any selected genotype
			delete $poplocallren{$poplocID}{$popall_ID};#delete it
			delete $poplocallvar{$poplocID}{$popall_ID};#delete it
		}
	}
}
for $poplocID (keys %poplocallren) {
	if (keys %{$poplocallren{$poplocID}} == 0) {#if no alleles left in locus
		delete $poplocallren{$poplocID};
		delete $poplocallvar{$poplocID};
	}
}
#}

#############################################################################################
#If export of genepop or Arlequin format active, verify if any locus has more than 99 alleles
#If yes, warn and inactivate genepop/Arlequin export
#############################################################################################

#{
if (($user_settings{'G'} == 1) or ($user_settings{'A'} == 1)) {
	for $poplocID (keys %poplocallren) {
		if (keys %{$poplocallren{$poplocID}} > 99) {
			print "Your selection contains at least one locus with more than 99 alleles.\n";
			if ($user_settings{'G'} == 1) {
				print "I cannot export genotypes in genepop 2-digit format.\n";
				$user_settings{'G'} = 0;
			}
			if ($user_settings{'A'} == 1) {
				print "I cannot export genotypes in Arlequin format.\n";
				$user_settings{'A'} = 0;
			}
			last;
		}
	}
}
#}

#############################################
#If -uS 2, -uS 3, -Fv 1, -Nxv 1, or Pv 1:
#Remove invariable loci and empty individuals
#############################################

#{
if (($user_settings{'uS'} == 2) or ($user_settings{'uS'} == 3)
	or ($user_settings{'Fv'} == 1) or ($user_settings{'Nxv'} == 1) or ($user_settings{'Pv'} == 1)) {
	for $poplocID (keys %poplocallren) {
		if (keys %{$poplocallren{$poplocID}} == 1) {#if there is only one allele
			delete $poplocallren{$poplocID};
			delete $poplocallvar{$poplocID};
			delete $sel_loc{$poplocID};
			delete $sel_gt{$poplocID};
			++$n_invarloc;#count locus
		}
	}
	if ($n_invarloc > 0) {#if you removed a locus
		print "I removed $n_invarloc invariable loci from your selection.\n";
		#Check for empty individuals and remove them
		%temphash1 = ();
		for $ind (keys %sel_ind) {
			$sel_ind{$ind} = 0;
			$temphash1{$ind} = 1;
		}
		for $poplocID (keys %sel_gt) {
			for $ind (keys %{$sel_gt{$poplocID}}) {
				$sel_ind{$ind} = 1;
				delete $temphash1{$ind};
			}
			if (keys %temphash1 == 0) {
				last;
			}
		}
		%temphash1 = ();
		for $ind (sort keys %sel_ind) {
			if ($sel_ind{$ind} == 0) {
				delete $sel_ind{$ind};
				delete $indpopall{$ind};
				delete $struc_cols{$ind};
				for $genepop (keys %popgenepop) {#loop through populations in %popgenepop
					if (exists $popgenepop{$genepop}{$ind}) {#if this ind is in there
						if (keys %{$popgenepop{$genepop}} == 1) {#if it is the only one
							delete $popgenepop{$genepop};#delete the population
						} else {
							delete $popgenepop{$genepop}{$ind};#delete the individual
						}
					}
				}
				for $arlpop (keys %arlpops) {#loop through populations in %arlpops
					if (exists $arlpops{$arlpop}{$ind}) {#if this ind is in there
						if (keys %{$arlpops{$arlpop}} == 1) {#if it is the only one
							delete $arlpops{$arlpop};#delete the population
						} else {
							delete $arlpops{$arlpop}{$ind};#delete the individual
						}
					}
				}
				print "I removed individual $ind because it has no remaining genotypes.\n";
			}
		}
	}
}
#}

if (($user_settings{'uS'} == 2) or ($user_settings{'uS'} == 3)) {# if unlinked SNPs is active
	###################################################################################
	#Determine 1st SNP in each locus, replace popallvar with character in %poplocallvar
	#renumber alleles in poplocallren as 1234 for ACGT
	#if -uS 2: additionally renumber alleles as 01 in poplocallbin
	###################################################################################
	
	#{
	for $poplocID (keys %poplocallvar) {#loop through loci
		for $popall_ID (sort {$a <=> $b} keys %{$poplocallvar{$poplocID}}) {#loop through alleles
			push @popallvarthisloc, $poplocallvar{$poplocID}{$popall_ID};#copy popallvar
		}
		for ($pos = 0; $pos < length $popallvarthisloc[0]; ++$pos) {#loop through positions
			%varchar = ();
			for ($i = 0; $i < @popallvarthisloc; ++$i) {#loop through alleles
				$char = substr($popallvarthisloc[$i],$pos,1);
				++$varchar{$char};
			}
			if (keys %varchar > 1) {#if position variable
				last;#stop searching
			}
		}
		if ($user_settings{'uS'} == 2) {#if only biallelic SNPS shall be included
			if (keys %varchar > 2) {#if SNP has more than 2 alleles, remove locus
				delete $sel_loc{$poplocID};
				delete $sel_gt{$poplocID};
				delete $poplocallren{$poplocID};
				delete $poplocallvar{$poplocID};
				++$n_multSNP;#count excluded locus
				#locus is still in %indpopall but that doesn't matter
			} else {#SNP has 2 alleles
				#populate %charnumbin less frequent allele: 0, more frequent allele: 1
				$i = 0;
				for $char (sort {$varchar{$a} <=> $varchar{$b}} keys %varchar) {
					$charnumbin{$char} = $i;
					++$i;
				}
				#loop through alleles again
				#replace popallvar with character at first SNP
				$i = 0;
				for $popall_ID (sort {$a <=> $b} keys %{$poplocallvar{$poplocID}}) {
					$poplocallvar{$poplocID}{$popall_ID} = substr($popallvarthisloc[$i],$pos,1);
					++$i;
					#Determine allele number as 1234 for ACGT in %poplocallren
					$poplocallren{$poplocID}{$popall_ID} = $charnum{$poplocallvar{$poplocID}{$popall_ID}};
					#Determine allele number as 0 or 1 for rarer/more frequent char in %poplocallbin
					$poplocallbin{$poplocID}{$popall_ID} = $charnumbin{$poplocallvar{$poplocID}{$popall_ID}};
				}
			}
		}
		else {#if SNPs with more than 2 alleles shall be included
			#loop through alleles again
			#replace popallvar with character at first SNP
			$i = 0;
			for $popall_ID (sort {$a <=> $b} keys %{$poplocallvar{$poplocID}}) {
				$poplocallvar{$poplocID}{$popall_ID} = substr($popallvarthisloc[$i],$pos,1);
				++$i;
				#Determine allele number as 1234 for ACGT in %poplocallren
				$poplocallren{$poplocID}{$popall_ID} = $charnum{$poplocallvar{$poplocID}{$popall_ID}};
			}
		}
		%varchar = ();#set back
		%charnumbin = ();#set back
		@popallvarthisloc = ();#set back and next loc
	}
	if ($n_multSNP > 0) {#if SNPs with more than 2 alleles have been excluded
		print "I have excluded $n_multSNP unlinked SNP loci with more than 2 alleles.\n";
		#check for empty individuals again and remove
		%temphash1 = ();
		for $ind (keys %sel_ind) {
			$sel_ind{$ind} = 0;
			$temphash1{$ind} = 1;
		}
		for $poplocID (keys %sel_gt) {
			for $ind (keys %{$sel_gt{$poplocID}}) {
				$sel_ind{$ind} = 1;
				delete $temphash1{$ind};
			}
			if (keys %temphash1 == 0) {
				last;
			}
		}
		%temphash1 = ();
		for $ind (sort keys %sel_ind) {
			if ($sel_ind{$ind} == 0) {
				delete $sel_ind{$ind};
				delete $indpopall{$ind};
				delete $struc_cols{$ind};
				for $genepop (keys %popgenepop) {#loop through populations in %popgenepop
					if (exists $popgenepop{$genepop}{$ind}) {#if this ind is in there
						if (keys %{$popgenepop{$genepop}} == 1) {#if it is the only one
							delete $popgenepop{$genepop};#delete the population
						} else {
							delete $popgenepop{$genepop}{$ind};#delete the individual
						}
					}
				}
				for $arlpop (keys %arlpops) {#loop through populations in %arlpops
					if (exists $arlpops{$arlpop}{$ind}) {#if this ind is in there
						if (keys %{$arlpops{$arlpop}} == 1) {#if it is the only one
							delete $arlpops{$arlpop};#delete the population
						} else {
							delete $arlpops{$arlpop}{$ind};#delete the individual
						}
					}
				}
				print "I removed individual $ind because it has no remaining genotypes.\n";
			}
		}
	}
	#}

} else {#if unlinked SNPs is not active
	if (($user_settings{'Fv'} == 1) or ($user_settings{'Nxv'} == 1) or ($user_settings{'Pv'} == 1)) {#if Fv, Nxv or Pv is active
		######################################################################
		#Determine all variable positions in each locus (within selected data)
		#Store in %poplocallvar (instead of original  popallvar)
		######################################################################
		
		#{
		for $poplocID (keys %poplocallvar) {#loop through loci
			for $popall_ID (sort {$a <=> $b} keys %{$poplocallvar{$poplocID}}) {#loop through alleles
				push @popallvarthisloc, $poplocallvar{$poplocID}{$popall_ID};#copy popallvar
			}
			for ($pos = 0; $pos < length $popallvarthisloc[0]; ++$pos) {#loop through positions
				for ($i = 0; $i < @popallvarthisloc; ++$i) {#loop through alleles
					$char = substr($popallvarthisloc[$i],$pos,1);
					$varchar{$char} = 1;
				}
				if (keys %varchar == 1) {#if position is invariable
					for ($i = 0; $i < @popallvarthisloc; ++$i) {#loop through alleles again
						substr($popallvarthisloc[$i],$pos,1) = "";#remove position
					}
				}
				%varchar = ();#set back and next position
			}
			#Replace old popallvar with new popallvar
			$i = 0;
			for $popall_ID (sort {$a <=> $b} keys %{$poplocallvar{$poplocID}}) {
				$poplocallvar{$poplocID}{$popall_ID} = $popallvarthisloc[$i];
				++$i;
			}
			@popallvarthisloc = ();#set back and next loc
		}
		#}
	}
	
	##################################
	#renumber alleles in %poplocallren
	##################################
	
	#{
	for $poplocID (sort {$a <=> $b} keys %poplocallren) {
		$popallren = 1;
		for $popall_ID (sort {$a <=> $b} keys %{$poplocallren{$poplocID}}) {
			$poplocallren{$poplocID}{$popall_ID} = $popallren;
			++$popallren;
		}
	}
	#}
}

################
#Populate %ex_gt
################

#{
@temparr1 = ();
for $poplocID (keys %sel_gt) {
	for $ind (keys %{$sel_gt{$poplocID}}) {
		#collect numbers of occuring alleles in this ind
		for $popall_ID (keys %{$indpopall{$ind}{$poplocID}}) {
			push @temparr1, $poplocallren{$poplocID}{$popall_ID};
		}
		#if there is one allele
		if (@temparr1 == 1) {
			@{$ex_gt{$ind}{$poplocID}} = ($temparr1[0],$temparr1[0]);
		}
		#if there are two alleles
		elsif (@temparr1 == 2) {
			@temparr1 = sort {$a <=> $b} @temparr1;#sort allele numbers ascending
			@{$ex_gt{$ind}{$poplocID}} = @temparr1;
		}
		#if there are more than two alleles
		elsif (@temparr1 > 2) {
			print "Your selection contains genotypes with more than 2 alleles. Exiting..\n";
			exit;
		}
		#else: don't know
		else {
			print "$ind $poplocID Something went wrong. Exiting..\n";
			exit;
		}
		@temparr1 = ();#set back and next ind	
	}
}
#}

#########################################
#If active, produce outfile for STRUCTURE
#########################################

#{
if ($user_settings{'S'} == 1) {
	if ($user_settings{'uS'} == 2) {#if -uS 2 is active, add uSb_ to outfilename
		$struc_outname = 'uSb_' . $struc_outname;
	}
	elsif ($user_settings{'uS'} == 3) {#if -uS 3 is active, add uS_ to outfile name
		$struc_outname = 'uS_' . $struc_outname;
	}
	$struc_outname = $outdirname . $struc_outname;
	#open outfile or die
	unless(open(OUTFILE,">$struc_outname")) {
		print "Cannot open $struc_outname. Exiting..\n";
		exit;
	}
	print "Printing data for STRUCTURE to file $struc_outname...\n";
	$strucmis = $user_settings{'Sm'};
	for $struc_orderNo (sort {$a <=> $b} keys %struc_order) {#loop through user-input-order of individuals
		$ind = $struc_order{$struc_orderNo};#get individual ID
		if (defined $struc_cols{$ind}) {
			#print first line for this ind
			print OUTFILE "$struc_cols{$ind}";#print pre-genotype columns
			for $poplocID (sort {$a <=> $b} keys %sel_gt) {#loop through selected loci, sorted
				if (defined $sel_gt{$poplocID}{$ind}) {#if genotype is selected
					print OUTFILE "\t$ex_gt{$ind}{$poplocID}[0]";#print first allele of genotype
				} else {#if not, print missing datapoint
					print OUTFILE "\t$strucmis";
				}
			}
			print OUTFILE "\n";
			#print second line for this ind
			print OUTFILE "$struc_cols{$ind}";#print pre-genotype columns
			for $poplocID (sort {$a <=> $b} keys %sel_gt) {#loop through selected loci, sorted
				if (defined $sel_gt{$poplocID}{$ind}) {#if genotype is selected
					print OUTFILE "\t$ex_gt{$ind}{$poplocID}[1]";#print second allele of genotype
				} else {#if not, print missing datapoint
					print OUTFILE "\t$strucmis";
				}
			}
			print OUTFILE "\n";
		}
	}
	close OUTFILE;
}
#}

#######################################
#If active, produce outfile for Genepop
#######################################

#{
if ($user_settings{'G'} == 1) {
	if ($user_settings{'uS'} == 2) {#if -uS 2 is active, add uSb_ to outfilename
		$genepop_outname = 'uSb_' . $genepop_outname;
	}
	elsif ($user_settings{'uS'} == 3) {#if -uS 3 is active, add uS_ to outfilename
		$genepop_outname = 'uS_' . $genepop_outname;
	}
	$genepop_outname = $outdirname . $genepop_outname;
	#open outfile or die
	unless(open(OUTFILE,">$genepop_outname")) {
		print "Cannot open $genepop_outname. Exiting..\n";
		exit;
	}
	print "Printing data for Genepop to file $genepop_outname...\n";
	print OUTFILE "Title\n";#title line for genepop
	#print locus section
	for $poplocID (sort {$a <=> $b} keys %sel_loc) {
		print OUTFILE "$poplocID\n";
	}
	for $genepop (sort {$a <=> $b} keys %popgenepop) {#loop through populations, sorted
		print OUTFILE "Pop\n";#Pop keyword for genepop
		for $ind (sort keys %{$popgenepop{$genepop}}) {#loop through individuals, sorted
			$indID = $indIDs{$ind};
			print OUTFILE "$indID ,";#individual ID comma
			for $poplocID (sort {$a <=> $b} keys %sel_gt) {#loop through selected loci, sorted
				print OUTFILE " ";#print whitespace
				if (defined $sel_gt{$poplocID}{$ind}) {#if genotype is selected
					printf OUTFILE ("%02d",$ex_gt{$ind}{$poplocID}[0]);#print first allele in 2digit format
					printf OUTFILE ("%02d",$ex_gt{$ind}{$poplocID}[1]);#print second allele in 2digit format
				} else {#if not selected: print missing datapoint in 2 digit format
					printf OUTFILE ("%02d",$genepopmis);#first missing allele
					printf OUTFILE ("%02d",$genepopmis);#second missing allele
				}
			}
			print OUTFILE "\n";
		}
	}
	close OUTFILE;
}
#}

########################################
#If active, produce outfile for Arlequin
########################################

#{
if ($user_settings{'A'} == 1) {
	if ($user_settings{'uS'} == 2) {#if -uS 2 is active, add uSb_ to outfilename
		$arl_outname = 'uSb_' . $arl_outname;
	}
	elsif ($user_settings{'uS'} == 3) {#if -uS 3 is active, add uS_ to outfilename
		$arl_outname = 'uS_' . $arl_outname;
	}
	$arl_outname = $outdirname . $arl_outname;
	#open outfile or die
	unless(open(OUTFILE,">$arl_outname")) {
		print "Cannot open $arl_outname. Exiting..\n";
		exit;
	}
	#determine some variables for output
	$arl_nsamp = keys %arlpops;
	@temparr1 = ();
	for $ind (keys %sel_ind) {
		$lenidarl = length $ind;
		push @temparr1, $lenidarl;
	}
	@temparr1 = reverse sort {$a <=> $b} @temparr1;
	$lenidarl = ($temparr1[0] + 1);
	print "Printing data for Arlequin to file $arl_outname...\n";
	print OUTFILE "[Profile]\nTitle=\"Title\"\n\n";
	print OUTFILE "NbSamples=$arl_nsamp\n";
	print OUTFILE "DataType=STANDARD\n",
	"GenotypicData=1\n",
	"LocusSeparator=WHITESPACE\n",
	"GameticPhase=0\n",
	"RecessiveData=0\n",
	"MissingData=\"?\"\n\n",
	"[Data]\n",
	"[[Samples]]\n";
	for $arlpop (sort keys %arlpops) {#loop through populations, alphanumeric sorting
		$arlsamsize = keys %{$arlpops{$arlpop}};
		print OUTFILE "SampleName=\"$arlpop\"\n";
		print OUTFILE "SampleSize=$arlsamsize\n",
		"SampleData= {\n";
		for $ind (sort keys %{$arlpops{$arlpop}}) {#loop through individuals in pop, alphanumeric sorting
			#print 1st line for this ind
			print OUTFILE "$ind";
			$tempstring1 = ' ' x ($lenidarl - length $ind);
			print OUTFILE "$tempstring1","1 ";
			for $poplocID (sort {$a <=> $b} keys %sel_gt) {#loop through loci ascending
				print OUTFILE " ";
				if (defined $sel_gt{$poplocID}{$ind}) {#if genotype is selected
					printf OUTFILE ("%02d",$ex_gt{$ind}{$poplocID}[0]);#print first allele in 2digit format
				} else {#if not selected: print missing datatype in 2 digit format
					print OUTFILE " ?";
				}
			}
			print OUTFILE "\n";
			#print 2nd line for this ind
			$tempstring1 = ' ' x $lenidarl;
			print OUTFILE "$tempstring1","  ";
			for $poplocID (sort {$a <=> $b} keys %sel_gt) {#loop through loci ascending
				print OUTFILE " ";
				if (defined $sel_gt{$poplocID}{$ind}) {#if genotype is selected
					printf OUTFILE ("%02d",$ex_gt{$ind}{$poplocID}[1]);#print second allele in 2digit format
				} else {#if not selected: print missing datatype in 2 digit format
					print OUTFILE " ?";
				}
			}
			print OUTFILE "\n";			
		}
		print OUTFILE "}\n";
	}
	#print a structure with all samples in one group
	print OUTFILE "\n[[Structure]]\n\n",
	"StructureName=\"One group for all samples\"\n",
	"NbGroups=1\n\n",
	"Group={\n";
	for $arlpop (sort keys %arlpops) {#loop through populations, alphanumeric sorting
		print OUTFILE "\"$arlpop\"\n";
	}
	print OUTFILE "}\n\n\n";
	close OUTFILE;
}
#}

###############################################################
#If FASTA, Nexus or PHYLIP is active, produce positions outfile
#with start and end positions of loci
###############################################################

#{
if (($user_settings{'Fv'} == 1) or ($user_settings{'Nxv'} == 1) or ($user_settings{'Pv'} == 1)) {
	if ($user_settings{'bin'} == 0) {#if DNA characters shall be exported
		if ($user_settings{'uS'} == 2) {#if -uS 2 is active, add uSb_ to outfilename
			$varpos_outname = 'uSb_' . $varpos_outname;
		}
		elsif ($user_settings{'uS'} == 3) {#if -uS 3 is active, add uS_ to outfilename
			$varpos_outname = 'uS_' . $varpos_outname;
		}
	}
	elsif ($user_settings{'bin'} == 1) {#if integer data shall be exported coded as 012 in one seq per ind
		$varpos_outname = 'uSbin1_' . $varpos_outname;#add uSbin1_ to outfilename
	}
	else {#if binary data shall be exported in two seqs per ind
		$varpos_outname = 'uSbin2_' . $varpos_outname;#add uSbin2_ to outfilename
	}
	$varpos_outname = $outdirname . $varpos_outname;
	unless(open(POSOUT,">$varpos_outname")) {
		print "Cannot open $varpos_outname. Exiting..\n";
		exit;
	}
	#Determine data for positions outfile and print
	print "Printing start and end positions of loci to file $varpos_outname...\n";
	print POSOUT "poplocID\tstart\tend\n";#header line
	for $poplocID (sort {$a <=> $b} keys %poplocallvar) {
		$popall_ID = (sort {$a <=> $b} keys %{$poplocallvar{$poplocID}})[0];
		$start = $end + 1;
		$end = $start - 1 + length $poplocallvar{$poplocID}{$popall_ID};
		print POSOUT "$poplocID\t$start\t$end\n";
	}
	close POSOUT;
}
#}

###################################################
#If FASTA, Nexus or PHYLIP is active: populate %ali
###################################################

#{
if (($user_settings{'Fv'} == 1) or ($user_settings{'Nxv'} == 1) or ($user_settings{'Pv'} == 1)) {
	$tempstring1 = '';
	$tempstring2 = '';
	for $ind (sort keys %sel_ind) {
		for $poplocID (sort {$a <=> $b} keys %sel_gt) {
			if (defined $sel_gt{$poplocID}{$ind}) {#if genotype selected
				#get popall_IDs
				@temparr1 = (sort {$a <=> $b} keys %{$indpopall{$ind}{$poplocID}});
				#look up popallvar for both alleles and append to out-sequences
				if (@temparr1 == 1) {#if ind is homozygous
					if ($user_settings{'bin'} eq '1') {#collect integer data in one sequence
						$val = $poplocallbin{$poplocID}{$temparr1[0]} * 2;
						$tempstring1 .= $val;
					}
					elsif ($user_settings{'bin'} eq '2') {#collect binary data in 2 sequences
						$tempstring1 .= $poplocallbin{$poplocID}{$temparr1[0]};
						$tempstring2 .= $poplocallbin{$poplocID}{$temparr1[0]};
					} else {#collect DNA characters in 2 sequences
						$tempstring1 .= $poplocallvar{$poplocID}{$temparr1[0]};
						$tempstring2 .= $poplocallvar{$poplocID}{$temparr1[0]};
					}
				}
				elsif (@temparr1 == 2) {#if ind is heterozygous
					if ($user_settings{'bin'} eq '1') {#collect integer data in one sequence
						$val = $poplocallbin{$poplocID}{$temparr1[0]} + $poplocallbin{$poplocID}{$temparr1[1]};
						$tempstring1 .= $val;
					} else {
						@temparr1 = shuffle @temparr1;#randomize order of popall_IDs
						if ($user_settings{'bin'} eq '2') {#collect binary data in 2 sequences
							$tempstring1 .= $poplocallbin{$poplocID}{$temparr1[0]};
							$tempstring2 .= $poplocallbin{$poplocID}{$temparr1[1]};
						} else {#collect DNA characters in 2 sequences
							$tempstring1 .= $poplocallvar{$poplocID}{$temparr1[0]};
							$tempstring2 .= $poplocallvar{$poplocID}{$temparr1[1]};
						}
					}
				}
				else {#don't know
					print "$ind poploc $poplocID number of alleles not 1 not 2, Exiting..\n";
					exit;
				}				
			} else {#genotype not selected
				#look up length of popallvar for this locus
				$popall_ID = (sort {$a <=> $b} keys %{$poplocallvar{$poplocID}})[0];
				$loclen = length $poplocallvar{$poplocID}{$popall_ID};
				$popallvar = $alimis x $loclen;#build string of missing data character dummies
				if ($user_settings{'bin'} eq '1') {#only one sequence
					$tempstring1 .= $popallvar;
				} else {#two sequences
					$tempstring1 .= $popallvar;
					$tempstring2 .= $popallvar;
				}
			}
		}
		#store data in %ali
		$indID = $indIDs{$ind};
		if ($user_settings{'bin'} eq '1') {#one sequence
			$ali{$indID} = $tempstring1;
		} else {#two sequences
			$indID1 = $indID . '_1';
			$indID2 = $indID . '_2';
			$ali{$indID1} = $tempstring1;
			$ali{$indID2} = $tempstring2;
		}
		$tempstring1 = '';#set back
		$tempstring2 = '';#set back
	}
}
#}


########################################################
#If active produce FASTA outfile with variable positions
########################################################

#{
if ($user_settings{'Fv'} == 1) {
	if ($user_settings{'bin'} == 0) {#if DNA characters shall be exported
		if ($user_settings{'uS'} == 2) {#if -uS 2 is active, add uSb_ to outfilename
			$fasvarfas_outname = 'uSb_' . $fasvarfas_outname;
		}
		elsif ($user_settings{'uS'} == 3) {#if -uS 3 is active, add uS_ to outfilename
			$fasvarfas_outname = 'uS_' . $fasvarfas_outname;
		}
	}
	elsif ($user_settings{'bin'} == 1) {#if integer data shall be exported coded as 012 in one seq per ind
		$fasvarfas_outname = 'uSbin1_' . $fasvarfas_outname;#add uSbin1_ to outfilename
	}
	else {#if binary data shall be exported in two seqs per ind
		$fasvarfas_outname = 'uSbin2_' . $fasvarfas_outname;# add uSbin2_ to outfilename
	}
	$fasvarfas_outname = $outdirname . $fasvarfas_outname;
	#open outfile or die
	unless(open(FASOUT,">$fasvarfas_outname")) {
		print "Cannot open $fasvarfas_outname. Exiting..\n";
		exit;
	}
	#Print data to fasta outfile
	print "Printing data in FASTA format to file $fasvarfas_outname...\n";
	for $indID (sort keys %ali) {
		$seq = $ali{$indID};
		$seq =~ s/$alimis/$user_settings{'Fmis'}/g;#replace missing char dummy with miss char
		print FASOUT ">$indID\n$seq\n";
	}
	close FASOUT;	
}
#}

#########################################################
#If active, produce Nexus outfile with variable positions
#########################################################

#{
if ($user_settings{'Nxv'} == 1) {
	if ($user_settings{'bin'} == 0) {#if DNA characters shall be exported
		if ($user_settings{'uS'} == 2) {#if -uS 2 is active, add uSb_ to outfilename
			$varnx_outname = 'uSb_' . $varnx_outname;
		}
		elsif ($user_settings{'uS'} == 3) {#if -uS 3 is active, add uS_ to outfilename
			$varnx_outname = 'uS_' . $varnx_outname;
		}
	}
	elsif ($user_settings{'bin'} == 1) {#if integer data shall be exported coded as 012 in one seq per ind
		$varnx_outname = 'uSbin1_' . $varnx_outname;#add uSbin1_ to outfilename
	}
	else {#if binary data shall be exported in two seqs per ind
		$varnx_outname = 'uSbin2_' . $varnx_outname;#add uSbin2_ to outfilename
	}
	$varnx_outname = $outdirname . $varnx_outname;
	#open outfile or die
	unless(open(NXOUT,">$varnx_outname")) {
		print "Cannot open $varnx_outname. Exiting..\n";
		exit;
	}
	#print section before data
	$ntax = keys %sel_ind;
	if ($user_settings{'bin'} == 1) {
		print NXOUT "#NEXUS\n\nBegin data;\n",
		"\tDimensions ntax=$ntax nchar=$end;\n",
		"\tFormat datatype=integerdata symbols=\"012\" gap=$user_settings{'Nxgap'};\n",
		"\tMatrix\n";
	}
	elsif ($user_settings{'bin'} == 2) {
		$ntax *= 2;
		print NXOUT "#NEXUS\n\nBegin data;\n",
		"\tDimensions ntax=$ntax nchar=$end;\n",
		"\tFormat datatype=binary symbols=\"01\" gap=$user_settings{'Nxgap'};\n",
		"\tMatrix\n";
	}
	else {
		$ntax *= 2;
		print NXOUT "#NEXUS\n\nBegin data;\n",
		"\tDimensions ntax=$ntax nchar=$end;\n",
		"\tFormat datatype=dna gap=$user_settings{'Nxgap'};\n",
		"\tMatrix\n";
	}
	#Print data to nexus outfile
	print "Printing data in Nexus format to file $varnx_outname...\n";
	for $indID (sort keys %ali) {
		$seq = $ali{$indID};
		$seq =~ s/$alimis/$user_settings{'Nxgap'}/g;#replace missing char dummy with miss char
		print NXOUT "$indID $seq\n";
	}
	#Print section after data
	print NXOUT "\t;\nEnd;\n";
	close NXOUT;	
}
#}

##########################################################
#If active, produce PHYLIP outfile with variable positions
##########################################################

#{
if ($user_settings{'Pv'} == 1) {
	if ($user_settings{'bin'} == 0) {#if DNA characters shall be exported
		if ($user_settings{'uS'} == 2) {#if -uS 2 is active, add uSb_ to outfilename
			$Ph_outname = 'uSb_' . $Ph_outname;
		}
		elsif ($user_settings{'uS'} == 3) {#if -uS 3 is active, add uS_ to outfilename
			$Ph_outname = 'uS_' . $Ph_outname;
		}
	}
	elsif ($user_settings{'bin'} == 1) {#if integer data shall be exported coded as 012 in one seq per ind
		$Ph_outname = 'uSbin1_' . $Ph_outname;#add uSbin1_ to outfilename
	}
	else {#if binary data shall be exported in two seqs per ind 
		$Ph_outname = 'uSbin2_' . $Ph_outname;#add uSbin2_ to outfilename
	}
	$Ph_outname = $outdirname . $Ph_outname;
	#open outfile or die
	unless(open(PHOUT,">$Ph_outname")) {
		print "Cannot open $Ph_outname. Exiting..\n";
		exit;
	}
	#Determine ID length, IDs will be filled up to this length with spaces
	#ID length will be that of the longest ID or 10, whatever is longer
	$lonID = (reverse sort {length($a) <=> length($b)} keys %ali)[0];
	$lonIDlen = length($lonID);
	if ($lonIDlen < $PhminIDlen) {
		$lonIDlen = $PhminIDlen;
	}
	#Print first line
	$ntax = keys %sel_ind;
	unless ($user_settings{'bin'} == 1) {
	$ntax *= 2;
	}
	print PHOUT " $ntax $end\n";
	#Print data to PHYLIP outfile
	print "Printing data in PHYLIP format to file $Ph_outname...\n";
	for $indID (sort keys %ali) {
		$seq = $ali{$indID};
		$seq =~ s/$alimis/$user_settings{'Pmis'}/g;#replace missing char dummy with miss char
		printf PHOUT "%-*s%s%s\n", $lonIDlen, $indID, $Pspacer, $seq;
	}	
	close PHOUT;	
}
#}

##############
#Final message
##############

print "Done. New files are in directory export.\n";
