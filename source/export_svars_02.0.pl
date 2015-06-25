#!/usr/bin/perl -w
#export_svars version 02.0 Copyright 2015 Andreas Hapke
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
use Math::Round qw(nearest);
use IO::File; #requires Perl 5.004 or higher
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

#Keyword Infilecolumns! Order of columns in an infile is important.

#Copyright message to screen
print "export_svars version 02.0 Copyright 2015 Andreas Hapke\n",
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
my $popallfn = 'popall.txt';#name of popall file out of poploc
my $indpopallsuff = '_indpopall.txt';#suffix of indpopall file out of indpoploc
my $allelessuff = '_alleles.txt';#suffix of alleles file out of indloc
my $svarssuff = '_svars.txt';#suffix of svars file out of indloc
my $svarsfh = '';#handle for svars infile
my $infilename = '';#an infilename
my $outdir = 'export';#directory for outfile
my $outfilesuff = '.fas';#suffix for outfile
my $outfilename = '';#the outfilename
my $ind = '';#individual ID
my $poplocID = 0;#locus ID
my $popall_ID = 0;#allele ID
my $popalldep = 0;#depth of allele in individual
my $indlocID = 0;#individual locus ID
my $indall_ID = 0;#individual allele ID
my $svarID = 0;#svar ID
my $svardep = 0;#depth of svar
my %indpopIDs = ();# {indlocID}{indall_ID}=popall_ID
my %popindIDs = ();# {popall_ID}=(indlocID,indall_ID)
my %popallseq = ();# {popall_ID}=(popalldep,popallseq)
my %indalldep = ();# {indall_ID} = depth
my %used_svars = ();# {indlocID}{indall_ID}{svarID}=(svardep,svarseq)
my %discarded_svars = ();# {indlocID}{svarID}=(svardep,svarseq)
my $header = '';#fasta header
my $outseq = '';#sequence in outfile
my $tempstring1 = '';
my @temparr1 = ();
#}

print "EXPORT_SVARS: ";

####################
#Take over arguments
####################

#{
if (@ARGV == 2) {
	$ind = $ARGV[0];
	$poplocID = $ARGV[1];
} else {
	print "I need two arguments: individual ID and poplocID.Exiting..\n";
	exit;
}
#}

print "individual $ind locus (poplocID): $poplocID\nReading data..\n";

#####################################################
#Read in *_indpopall.txt
#Populate %indpopIDs {indlocID}{indall_ID}=popall_ID
#Populate %popindIDs {popall_ID}=(indlocID,indall_ID)
#####################################################

#{
$infilename = $ind . $indpopallsuff;
unless(open(INFILE, "<", $infilename)) {
	print "Cannot open $infilename, exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	if ($temparr1[2] == $poplocID) {
		$indpopIDs{$temparr1[0]}{$temparr1[1]} = $temparr1[3];
		@{$popindIDs{$temparr1[3]}} = ($temparr1[0],$temparr1[1]);
	}
}
close INFILE;
if (keys %indpopIDs == 0) {
	print "Poploc $poplocID not genotyped in ind $ind, exiting..\n";
	exit;
}
#}

######################################################
#Read in popall.txt
#Populate %popallseq {popall_ID}=(popalldep,popallseq)
######################################################

#{
unless(open(INFILE, "<", $popallfn)) {
	print "Cannot open $popallfn, exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	if (($temparr1[0] == $poplocID) and (defined $popindIDs{$temparr1[1]})) {
		@{$popallseq{$temparr1[1]}} = (0,$temparr1[4]);
	}
}
close INFILE;
if (keys %popallseq == 0) {
	print "Poploc $poplocID does not exist, exiting..\n";
	exit;
}
#}

########################################
#Read in *_alleles.txt
#Populate %indalldep {indall_ID} = depth
########################################

#{
$infilename = $ind . $allelessuff;
unless(open(INFILE, "<", $infilename)) {
	print "Cannot open $infilename, exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	if (defined $indpopIDs{$temparr1[0]}) {
		$indalldep{$temparr1[3]} = $temparr1[5];
	}
}
close INFILE;
#}

###############################
#Transfer depth into %popallseq
###############################

#{
for $popall_ID (keys %popallseq) {
	$indall_ID = $popindIDs{$popall_ID}[1];
	$popallseq{$popall_ID}[0] = $indalldep{$indall_ID};
}
#}

#####################################################################
#Read in svars.txt
#Populate %used_svars {indlocID}{indall_ID}{svarID}=(svardep,svarseq)
#Populate %discarded_svars {indlocID}{svarID}=(svardep,svarseq)
#####################################################################

#{
$infilename = $ind . $svarssuff;
$tempstring1 = $infilename;
if (-f $infilename) {
	$svarsfh = IO::File->new("< $infilename");
} else {
	$infilename .= '.gz';#look for a gzipped svars file
	$svarsfh = new IO::Uncompress::Gunzip $infilename;
}
unless(defined $svarsfh) {
	print "Cannot open file $tempstring1 or $infilename, exiting..\n";
	exit;
}
$tempstring1 = <$svarsfh>;#skip headerline
while ($tempstring1 = <$svarsfh>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	if (defined $indpopIDs{$temparr1[0]}) {
		if ($temparr1[5] eq 'used') {
			@{$used_svars{$temparr1[0]}{$temparr1[3]}{$temparr1[4]}} = ($temparr1[7],$temparr1[6]);
		} else {
			@{$discarded_svars{$temparr1[0]}{$temparr1[4]}} = ($temparr1[7],$temparr1[6]);
		}
	}
}
close $svarsfh;
undef $svarsfh;
#}

#################################################
#Create directory export for outfile if necessary
#Open outfile
#################################################

#{
unless(-d $outdir) {
	unless(mkdir $outdir) {
		print "Cannot create directory export, exiting..\n";
		exit;
	}
}
$outfilename = $outdir . '/' . $ind . '_' . $poplocID . $outfilesuff;
unless(open(OUTFILE, ">", $outfilename)) {
	print "Cannot open $outfilename, exiting..\n";
	exit;
}
#}

######################
#Print data to outfile
######################

#{
#Print alleles
for $popall_ID (sort {$a <=> $b} keys %popallseq) {
	$popalldep = nearest(0.1,$popallseq{$popall_ID}[0]);
	$header = ">$popall_ID" . "_" . "$popalldep\n";
	$outseq = "$popallseq{$popall_ID}[1]\n";
	print OUTFILE "$header","$outseq";
}
#Print used svars
for $popall_ID (sort {$a <=> $b} keys %popallseq) {
	($indlocID,$indall_ID) = @{$popindIDs{$popall_ID}};
	for $svarID (sort {$a <=> $b} keys %{$used_svars{$indlocID}{$indall_ID}}) {
		$header = ">$popall_ID" . "_" . "$indlocID" . "_" . "$indall_ID" . "_" . "$svarID"
		. "_used_" . "$used_svars{$indlocID}{$indall_ID}{$svarID}[0]\n";
		$outseq = "$used_svars{$indlocID}{$indall_ID}{$svarID}[1]\n";
		print OUTFILE "$header","$outseq";
	}
}
#Print discarded svars
for $indlocID (sort {$a <=> $b} keys %discarded_svars) {
	for $svarID (sort {$a <=> $b} keys %{$discarded_svars{$indlocID}}) {
		$header = ">0_" . "$indlocID" . "_0_" . "$svarID" . "_discarded_"
		. "$discarded_svars{$indlocID}{$svarID}[0]\n";
		$outseq = "$discarded_svars{$indlocID}{$svarID}[1]\n";
		print OUTFILE "$header","$outseq";	
	}
}
close OUTFILE;
print "Saved data to file $outfilename\n";
exit;
#}
