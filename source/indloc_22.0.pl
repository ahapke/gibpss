#!/usr/bin/perl -w
#indloc version 22.0 Copyright 2015 Andreas Hapke
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
use Parallel::ForkManager;
use IO::File; #requires Perl 5.004 or higher
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError :constants);

#Note for the reader: Within the code of this program,
#"haplotye" refers to an svar (see documentation).

#uses subroutines
#ind
#indlen
#splitvar
#fno_f_h
#nper
#hpd
#distpwdN
#al
#rare_match_loc
#snphapcal1
#snphapcalbin
#varpos
#errbase1
#errbasebin
#hapgroups
#hapgroupsbin
#build_alleles
#build_alleles_bin
#consensus1

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "indloc version 22.0 Copyright 2015 Andreas Hapke\n",
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
#output-settings
my $individuals_fname = 'individuals.txt';#outfile: individual IDs, one per line
my $report_name = 'indloc_report.txt';#outfile: used settings
#input-settings
my %user_settings = (
i => 'indloc_infiles.txt',#name of infile-list
z => '0',#1: write svars.txt and reads.txt as gzip compressed files
d => 'none',#name of distance-file
D => '4',#hapdep
P => '1',#max_p
r => '0',#1: fastq-infiles contain only forward reads of rc-duplicates, build reverse reads as reverse complement
M => 'f',#f: frequency threshold method, b: binomial likelihood ratio method
c => '0.2',#minimum character frequency for frequency threshold method
a => '0.2',#minimum allele-frequency for frequency threshold method
e => '0.01'#sequencing error rate for binomial likelihood ratio method
);
my %defaults = (
i => 'indloc_infiles.txt',#name of infile-list
z => '0',#1: write svars.txt and reads.txt as gzip compressed files
d => 'none',#name of distance-file
D => '4',#hapdep
P => '1',#max_p
r => '0',#1: fastq-infiles contain only forward reads of rc-duplicates, build reverse reads as reverse complement
M => 'f',#f: frequency threshold method, b: binomial likelihood ratio method
c => '0.2',#minimum character frequency for frequency threshold method
a => '0.2',#minimum allele-frequency for frequency threshold method
e => '0.01'#sequencing error rate for binomial likelihood ratio method
);
#other variables
my $infilelistname = '';#each line: ind TAB fastq-infilename
my $z = 0;#1: write svars.txt and reads.txt as gzip compressed files
my %ind_infile = ();# {ind}=fastqinfilename
my @inds_inputorder = ();#individual IDs in input order
my $ind = '';#individual ID
my $fqfile = '';#fastq-infile
my $testfqfileh = '';#handle for test-open of fastq-infiles
my %dist_set = ();#distance-settings: {sl}=(hapdist,rarehapdist);
my $sl = 0;#sequence length
my $infilename = '';
my $hapdist = 4;#max distance between neighboring haplotypes in a locus
my $rarehapdist = 6;#max dist between a rare haplotype and at least one good haplotype
					#in a locus to allow for assignment of the rare haplotype to the locus
my $hapdep = 4;#min read-depth for a haplotype
my $build_rc = 0;#1: build reverse complement of each read in infile as read 2
my $method = 'f';#Method: f: frequency threshold, b: binomial likelihood ratio
my $min_cf = 0.2;#minimum character frequency in an individual locus alignment
my $min_af = 0.2;#minimum allele frequency in an individual locus alignment
my $binerr = 0.01;#sequencing error rate for binomial likelihood ratio method
my $max_p = 1;#maximum number of processes in parallel
my $n_arg = 0;#number of user-arguments
my $flag = '';#command-flag
my $val = '';#value for command-flag
my $indreportsuffix = '_indloc_report.txt';#suffix of report files from individuals
my $indreportname = '';#full name of individual report file
my $indreport = '';#contents of individual report file
my $tempstring1 = '';
my $tempstring2 = '';
my @temparr1 = ();
my %temphash1 = ();
my %temphash2 = ();
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
#check (at least some) arguments for plausibility (very limited test)
#z: must be 0 or 1
unless(($user_settings{'z'} eq '0') or ($user_settings{'z'} eq '1')) {
	$user_settings{'z'} = $defaults{'z'};
}
#hapdep: positive integer
unless(($user_settings{'D'} =~ /^\d+$/) and ($user_settings{'D'} > 0)) {
	$user_settings{'D'} = $defaults{'D'};
}
#max_p: positive integer
unless(($user_settings{'P'} =~ /^\d+$/) and ($user_settings{'P'} > 0)) {
	$user_settings{'P'} = $defaults{'P'};
}
#r: 0 or 1
unless(($user_settings{'r'} eq '0') or ($user_settings{'r'} eq '1')) {
		$user_settings{'r'} = $defaults{'r'};
}
#Method (-M) f (default): frequency threshold method, b binomial likelihood ratio method
unless(($user_settings{'M'} eq 'f') or ($user_settings{'M'} eq 'b')) {
	$user_settings{'M'} = $defaults{'M'};
}
#minimum character frequency for frequency threshold method: -1 or float between 0 and 1
#-1 inactivates error correction by sub errbase1, float between 0 and 1 is used as
#minimum frequency
unless((($user_settings{'c'} =~ /^-?\d$/) or ($user_settings{'c'} =~ /^\d\.\d+$/))
and (($user_settings{'c'} == -1) or (($user_settings{'c'} >= 0) and ($user_settings{'c'} <= 1)))) {
	$user_settings{'c'} = $defaults{'c'};
}
#minimum allele frequency for frequency threshold method: float between 0 and 1
unless((($user_settings{'a'} =~ /^\d$/) or ($user_settings{'a'} =~ /^\d\.\d+$/))
and ($user_settings{'a'} >= 0) and ($user_settings{'a'} <= 1)) {
	$user_settings{'a'} = $defaults{'a'};
}
#sequencing error rate for binomial likelihood ratio method (-e)
#must be >= 0 and <1
unless((($user_settings{'e'} =~ /^\d$/) or ($user_settings{'e'} =~ /^\d\.\d+$/))
and ($user_settings{'e'} >= 0) and ($user_settings{'e'} < 1)) {
	$user_settings{'e'} = $defaults{'e'};
}

#Take over settings into variables
$infilelistname = $user_settings{'i'};
$z = $user_settings{'z'};
$hapdep = $user_settings{'D'};
$max_p = $user_settings{'P'};
$build_rc = $user_settings{'r'};
$method = $user_settings{'M'};
$min_cf = $user_settings{'c'};
$min_af = $user_settings{'a'};
$binerr = $user_settings{'e'};
#}

############################################
#Open report-file for append, print settings 
############################################

#{
unless(open(OUTREP, ">>$report_name")) {
	print "Cannot open file $report_name, exiting ...\n\n";
	exit;
}
print OUTREP
"INDLOC: used settings:\n",
"i $user_settings{'i'}\n",
"z $user_settings{'z'}\n",
"d $user_settings{'d'}\n",
"D $user_settings{'D'}\n",
"P $user_settings{'P'}\n",
"r $user_settings{'r'}\n",
"M $user_settings{'M'}\n";
if ($user_settings{'M'} eq 'f') {
print OUTREP
"c $user_settings{'c'}\n",
"a $user_settings{'a'}\n";
} elsif ($user_settings{'M'} eq 'b') {
print OUTREP
"e $user_settings{'e'}\n";
}
#}

######################################################
#Initialize distance settings and print to report-file
######################################################

#{
#if user did not provide a distance-filename
if ($user_settings{'d'} eq 'none') {
	#initialize the smallest possible sequence length with defaults
	$sl = 1;
	@{$dist_set{$sl}}=($hapdist,$rarehapdist);
}
else {#else take over distance filename
	$infilename = $user_settings{'d'};
}
#if you can open file with the distance-filename provided by user
if(open(INFILE,$infilename)) {
	#read in
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		@temparr1 = split (/\t/,$tempstring1);
		#take only lines with 3 entries
		#check entries for plausibility, all integers, sl >= 0, distances > 0
		if ((@temparr1 == 3) and ($temparr1[0] =~ /^\d+$/) and ($temparr1[1] =~ /^\d+$/)
		and ($temparr1[2] =~ /^\d+$/) and ($temparr1[0] >= 0)
		and ($temparr1[1] > 0) and ($temparr1[2] > 0)) {
			#take over
			$sl = shift @temparr1;
			@{$dist_set{$sl}}=@temparr1;
		} else {#unusable line
			print OUTREP "Ignored unusable line in distance settings file.\n";
		}
	}
	#if no line in distance file was usable
	if (keys %dist_set == 0) {
		#initialize the smallest possible sequence length with defaults
		$sl = 1;
		@{$dist_set{$sl}}=($hapdist,$rarehapdist);
		print OUTREP "All lines in distance settings file were unusable.\n";
	}	
}
#else: you could not open the file provided by user
else {
	print "Cannot open distance file: $infilename, use default values: $hapdist $rarehapdist\n";
	print OUTREP "Could not open distance settings file, used defaults: $hapdist $rarehapdist\n";
	#initialize the smallest possible sequence length with defaults
	$sl = 1;
	@{$dist_set{$sl}}=($hapdist,$rarehapdist);
}
print OUTREP "\nDistance settings:\n";
for $sl (sort {$a <=> $b} keys %dist_set) {
	print OUTREP "$sl $dist_set{$sl}[0] $dist_set{$sl}[1]\n";
}
#}

##############################################################################
#load infilelist with names of new individuals and corresponding fastq-infiles
##############################################################################

#{
unless(open(INFILE,$infilelistname)) {
	print "Cannot open $infilelistname. Exiting..\n";
	print OUTREP "Could not open $infilelistname. Analysis aborted.\n\n";
	exit;
}
%temphash1 = ();
%temphash2 = ();
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Unless the line has 2 entries: exit
	unless (@temparr1 == 2) {
		print "I don't understand this line:\n$tempstring1\nExiting..\n";
		print OUTREP "Unusable line in file $infilelistname:\n$tempstring1\nAnalysis aborted.\n\n";
		exit;
	}
	#Warn and stop if this individual ID already appeared in infilelist
	if (defined $temphash1{$temparr1[0]}) {
		print "At least one individual ID appears more than once in file $infilelistname. Exiting..\n";
		print OUTREP "At least one individual ID appeared twice in file $infilelistname.\nAnalysis aborted.\n";
		exit;		
	}
	#Warn and stop if this fastq infile already appeared in infilelist
	if (defined $temphash2{$temparr1[1]}) {
		print "At least one fastq file appears more than once in file $infilelistname. Exiting..\n";
		print OUTREP "At least one fastq file appeared twice in file $infilelistname.\nAnalysis aborted.\n";
		exit;		
	}
	$temphash1{$temparr1[0]} = 1;#store individual ID for check
	$temphash2{$temparr1[1]} = 1;#store fastq filename for check
	#store individual ID and fastq-infilenames in %ind_infile:
	$ind_infile{$temparr1[0]} = $temparr1[1];
	push @inds_inputorder, $temparr1[0];
}
close INFILE;
%temphash1 = ();
%temphash2 = ();
#If an individuals file already exists:
#Read and make sure that no new individual has same ID as an existing one.
if (-f $individuals_fname) {
	#open it or die
	unless(open(INFILE,$individuals_fname)) {
		print "Cannot open $individuals_fname. Exiting..\n";
		print OUTREP "Could not open $individuals_fname. Analysis aborted.\n\n";
		exit;
	}
	print "Cross-checking names of new and existing individuals.\n";
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		#if existing individual has same name a new one: warn, exit
		if (exists $ind_infile{$tempstring1}){
			print "At least one new individual has the same ID as an existing one. Exiting..\n";
			print OUTREP "At least one new individual had the same ID as an existing one. Analysis aborted.\n\n";
			exit;
		}
	}
	close INFILE;
}
#}

####################################################################
#If build_rc is active, check first two fastq headers in each infile
####################################################################

#{
if ($build_rc == 1) {
	for $ind (sort keys %ind_infile) {#loop through individuals
		$fqfile = $ind_infile{$ind};
		if ($fqfile =~ /.gz$/) {#open as gzip compressed file
			$testfqfileh = new IO::Uncompress::Gunzip $fqfile;
		} else {#open as text file
			$testfqfileh = IO::File->new("< $fqfile");
		}
		unless(defined $testfqfileh) {
			print "Cannot open $fqfile, exiting..\n";
			print OUTREP "Could not open $fqfile\nanalysis aborted.\n";
			close OUTREP;
			exit;
		}
		#Check the first two fastq headers: must end with "_1"
		$tempstring1 = <$testfqfileh>;
		$tempstring2 = <$testfqfileh>;
		$tempstring2 = <$testfqfileh>;
		$tempstring2 = <$testfqfileh>;
		$tempstring2 = <$testfqfileh>;
		chomp $tempstring1;
		chomp $tempstring2;
		unless(($tempstring1 =~ m/_1$/) and ($tempstring2 =~ m/_1$/)) {
			print "Not all fastq headers in file\n",
			"$fqfile\nend with \"_1\", analysis aborted\n";
			print OUTREP "Not all fastq headers in file\n",
			"$fqfile\nend with \"_1\", analysis aborted\n";
			close OUTREP;
			close INFILE;
			exit;
		}
		close $testfqfileh;
		undef $testfqfileh;
	}
}
#}

##############
#Data analysis
##############

#{
#loop through the individuals
my $pm1 = Parallel::ForkManager->new($max_p);
for $ind (@inds_inputorder) {
	$pm1->start and next;
	$fqfile = $ind_infile{$ind};#get fastq-filename for this ind
	#call sub ind to analyze infile of this ind and print output for this ind
	ind($ind,$fqfile,$build_rc,$hapdep,\%dist_set,$method,$min_cf,$min_af,$binerr,$indreportsuffix,$z);
	$pm1->finish;
}
$pm1->wait_all_children;
#}

##############################
#Check individual report files
#remove failed individuals
#from %ind_infile
##############################

#{
print OUTREP "\nAnalysis of individuals:\n";
for $ind (sort keys %ind_infile) {#loop through individuals
	$indreportname = $ind.$indreportsuffix;#build filename
	if(open(INDREPFILE,$indreportname)) {#if you could open report file
		while ($tempstring1 = <INDREPFILE>) {#read in report file
			$indreport .= $tempstring1;#collect contents
		}
		close INDREPFILE;
		if ($indreport =~ /analysis completed/i)	{#if report reports success
			print OUTREP "$ind analysis completed.\n";#print message to main report
			print "$ind analysis completed.\n";#and to screen
		} else {#no success message in report
			print OUTREP "$ind $indreport";#print message to main report
			print "$ind analysis failed.\n";#and to screen
			delete $ind_infile{$ind};#remove ind from %ind_infile
		}
	} else {#if you could not open the report file
		print OUTREP "$ind analysis failed.\n";#print message to main report
		print "$ind analysis failed.\n";#and to screen
		delete $ind_infile{$ind};#remove ind from %ind_infile
	}
	$indreport = '';#set back and next ind
}
#}

#######################################################
#Print IDs of new individuals
#which are still in %ind_infile to file individuals.txt
#If this file already exists, append new IDs
#######################################################

#{
unless(open(OUTFILE, ">>$individuals_fname")) {
	print "Analyses are finished but I can't open $individuals_fname.\n",
	"IDs of the new individuals will be missing from that file.\n",
	"Please add them before you proceed with the next program.\n",
	"Exiting..\n";
	print OUTREP "Could not open $individuals_fname.\n",
	"IDs of the new individuals will be missing from that file.\n",
	"Please add them before you proceed with the next program or next run of this program.\n\n";
	exit;	
}
for $ind (sort keys %ind_infile) {
	print OUTFILE "$ind\n";
}
close OUTFILE;
#}


#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";
print OUTREP "Analysis completed.\nRun took $run_s seconds.\n\n";
close OUTREP;

exit;

####################################
#Subroutines
####################################

#definition of subroutine ind
#does the work for one individual
#expects
#$ind: individual ID
#$fqfile: fastq-infile for this ind
#$build_rc: 1: fastq-infile contains only forward reads of rc-duplicates
			#all headers end with _1
			#build reverse read with same header ending with _2 as reverse complement
			#0: Not so
#$hapdep: (minimum read-depth for a "good" haplotype)
#ref to %dist_set: {sl}=(hapdist,rarehapdist);
		#hapdist: max dist between two neighboring "good" haplotypes in a locus
		#rarehapdist: max dist between a rare haplotype and at least one good haplotype
		#in a locus to allow for assignment of the rare haplotype to the locus
#$method: f: frequency threshold, b: binomial likelihood ratio
#$min_cf: minimum character frequency in an individual locus alignment for frequency threshold method
#$min_af: minimum allele frequency in an individual locus alignment for frequency threshold method
#$binerr: sequencing error rate for binomial likelihood ratio method
#$indreportsuffix: suffix for name of individual report file
#$z: 1: write svars.txt and reads.txt as gzip compressed files svars.txt.gz reads.txt.gz; 0: don't
sub ind {
	#declare and initialize: _r means a reference
	my($ind,$fqfile,$build_rc,$hapdep,$dist_set_r,$method,$min_cf,$min_af,$binerr,$indreportsuffix,$z) = @_;
	#own variables
	my $indreportname = '';#name of report file for current individual
	my $fqfileh = '';#handle for fastq infile
	my $hapdist = 0;#max dist between two neighboring "good" haplotypes in a locus
	my $rarehapdist = 0;#max dist between a rare haplotype and at least one good haplotype
						#in a locus to allow for assignment of the rare haplotype to the locus
	my $seqline1 = '';#line1 of current seqentry
	my $seqline2 = '';#line2 of current seqentry
	my $seqline3 = '';#line3 of current seqentry
	my $seqline4 = '';#line4 of current seqentry
	my $rseqline2 = '';#new reverse complement sequence of current seqentry
	my $sl = 0;#length of current sequence
	my %lhd = ();# {seqlength}{haplotype}= depth (number of reads)
	my %l_loc_cat_h_f =();# {seqlength}{locusID}{cat}{haplotype}= haplotype-depth
	my $locID = 0;#ID of a locus
	my $nextlocID = 1;#locID to use for the first new locus in a new set of loci defined by its seqlength
	my $allele = '';#an allele
	my $depth = 0;#depth of an allele
	my $success = "analysis aborted\n";#will change to 'analysis completed' after completion of this ind
	
	#Open individual report file or return
	$indreportname = $ind.$indreportsuffix;#build filename
	unless(open(INDREPFILE, ">$indreportname")) {
		print "$ind: Can't open $indreportname. Analysis aborted\n";
		return;
	}	
	#open the fastq-infile or leave failure message and return
	if ($fqfile =~ /.gz$/) {#open as gzip compressed file
		$fqfileh = new IO::Uncompress::Gunzip $fqfile;
	} else {#open as text file
		$fqfileh = IO::File->new("< $fqfile");
	}
	unless(defined $fqfileh) {
		print INDREPFILE "Could not open $fqfile. Analysis aborted.\n";#print failure message to individual report file
		print "$ind: Cannot open $fqfile. Analysis aborted.\n";#and to screen
		close INDREPFILE;
		return;
	}
	print "Analyzing individual $ind\n";
	#Read in the fastq-infile
	if ($build_rc == 1) {
		#Fastq infile contains only f reads of rc-duplicates, build r-reads
		while ($seqline1 = <$fqfileh>) {
			$seqline2 = <$fqfileh>;
			$seqline3 = <$fqfileh>;
			$seqline4 = <$fqfileh>;
			chomp $seqline1;
			chomp $seqline2;
			unless ($seqline1 =~ m/_1/) {#Check if fastq-header ends with "_1"
				print "Not all fastq-headers end with \"_1\".Analysis aborted.\n";
				print INDREPFILE "Not all fastq-headers end with \"_1\".Analysis aborted.\n";
				close INDREPFILE;
				return;
			}
			$rseqline2 = reverse $seqline2;
			$rseqline2 =~ tr/ACGTacgt/TGCAtgca/;#reverse complement of current sequence
			$sl = length $seqline2;#determine length of sequences
			#store sequences in %lhd (seqlength}{haplotype}=depth and count depth
			++$lhd{$sl}{$seqline2};
			++$lhd{$sl}{$rseqline2};
		}
	} else {#Do not build r-reads
		while ($seqline1 = <$fqfileh>) {
			$seqline2 = <$fqfileh>;
			$seqline3 = <$fqfileh>;
			$seqline4 = <$fqfileh>;
			chomp $seqline2;
			$sl = length $seqline2;#determine length of sequence
			#store sequence in %lhd (seqlength}{haplotype}=depth and count depth
			++$lhd{$sl}{$seqline2};
		}
	}
	close $fqfileh;
	undef $fqfileh;
	#determine smallest defined sequence length in %dist_set owned by main
	$sl = (sort {$a <=> $b} keys %{$dist_set_r})[0];
	#initialize distances with those for this sequence length
	$hapdist = $$dist_set_r{$sl}[0];
	$rarehapdist = $$dist_set_r{$sl}[1];
	
	#loop through the sequence lengths:
	for $sl (sort {$a <=> $b} keys %lhd) {
		#if current sequence length is defined in %dist_set,
		if (defined $$dist_set_r{$sl}) {
			#update distance settings
			$hapdist = $$dist_set_r{$sl}[0];
			$rarehapdist = $$dist_set_r{$sl}[1];
		}
		#call sub indlen to analyze all sequences of this ind with this length
		#indlen will populate %l_loc_cat_h_f for this ind and this seq-length
		indlen(\%lhd,\%l_loc_cat_h_f,$sl,$hapdep,$hapdist,$rarehapdist,\$nextlocID);
	}
	if ($method eq 'f') {#use frequency threshold method:
		#call sub snphapcal1 to call SNPs and haplotypes and produce outfiles for this ind
		#sub snphapcal1 will delete haplotypes from %l_loc_cat_h_f once it doesn't need them any more
		#$succes will contain failure message or "analysis completed"
		$success = snphapcal1(\%l_loc_cat_h_f,$ind,$min_cf,$min_af,$fqfile,$build_rc,$z);
		print INDREPFILE "$success";#print failure message or "analysis completed to report file
		close INDREPFILE;#close report file and return
	}
	elsif ($method eq 'b') {#use binomial likelihood ratio method
		#call sub snphapcalbin to call SNPs and haplotypes and produce outfiles for this ind
		#sub snphapcalbin will delete haplotypes from %l_loc_cat_h_f once it doesn't need them any more
		#$succes will contain failure message or "analysis completed"
		$success = snphapcalbin(\%l_loc_cat_h_f,$ind,$binerr,$fqfile,$build_rc,$z);
		print INDREPFILE "$success";#print failure message or "analysis completed to report file
		close INDREPFILE;#close report file and return
	}
}

#definition of subroutine indlen
#analyzes all haplotypes of one length of one ind
#expects 7 scalars:
#reference to %lhd (already populated)
	# {seqlength}{haplotype}= depth
#reference to %l_loc_cat_h_f
	# {seqlength}{locusID}{cat}{haplotype}=haplotype_depth
#$sl:		seqlength
#$hapdep	min read-depth for a good haplotype
#$hapdist	max dist between neighboring haplotypes in a locus
#rarehapdist (max dist between a rare haplotype and at least one good haplotype
		#in a locus to allow for assignment of the rare haplotype to the locus)
#ref to $nextlocID: locID to start with for the first locus built
#populates %l_loc_cat_h_f{$sl}
sub indlen{
	#declare and initialize: _r means a reference
	my ($lhd_r,$l_loc_cat_h_f_r,$sl,$hapdep,$hapdist,$rarehapdist,$nextlocID_r) = @_;
	my @splitvar = ();#variables for splitting of haplotypes into nonoverlapping fragments:
						#2d array: [0]: array of startpositions [1] array of fragment-lengths						
	my %goodhaplos = ();#contains only haplotypes >= hapdep:
						# {hapNo}=(sequence,depth,groupNo)
	my %rarehaplos = ();#contains only haplotypes < hapdep:
						# {hapNo}=(sequence,depth,locusID)
	my $hapNo = 0;#number of current good haplotype
	my $haplo = '';#current haplotype
	my $dep = 0;#depth of current haplotype
	my %fno_f_h = ();# {frag_No}{fragment}= @of haplotype-numbers
					#frag_No: 1..hapdist+1: frag1, frag2..frag(hapdist+1)
	my $frag_No = 0;#a frag_No
	my $fragment = '';#a fragment
	my %hpd = ();# {haplo1_No}{haplo2_No}=distance, contains all pairs of haplotypes
					#with a distance <= $hapdist;
	my $haplo1_No = 0;
	my $haplo2_No = 0;
	my %group_hNos = ();# {groupNo}=@ of haploNo
	
	#populate %goodhaplos, groupNo initially 0
	#loop through haplotypes in %lhd{$sl}
	for $haplo (keys %{$$lhd_r{$sl}}){
		$dep = $$lhd_r{$sl}{$haplo};#determine haplotype-depth
		#if depth >= hapdep, include in %goodhaplos
		if ($dep >= $hapdep) {
			$goodhaplos{$hapNo}[0] = $haplo;#sequence at pos 0
			$goodhaplos{$hapNo}[1] = $dep;#depth at pos 1
			$goodhaplos{$hapNo}[2] = 0;#group-No at pos 2, initially 0
			++$hapNo;#increment hapNo for next haplotype
			delete $$lhd_r{$sl}{$haplo};#delete haplotype from %lhd
		}
	}	
	#call subroutine splitvar to determine variables for
	#splitting of haplotypes into nonoverlapping fragments
	@splitvar = splitvar($sl,$hapdist);
	#call subroutine fno_f_h to populate %fno_f_h
	fno_f_h(\%goodhaplos,\@splitvar,\%fno_f_h,$sl);
	#call sub hpd to populate %hpd
	hpd(\%fno_f_h,\%goodhaplos,$hapdist,$sl,\%hpd);	
	undef %fno_f_h;#undefine the fragment catalog in %fno_f_h
	#call sub al to assemble good haplotypes into loci
	#sub al will populate %l_loc_cat_h_f{$sl} owned by sub ind
	#with "good haplotypes"
	#and delete haplotypes from %goodhaplos
	al(\%goodhaplos,\%hpd,\%group_hNos,$l_loc_cat_h_f_r,$sl,$nextlocID_r);
	undef %goodhaplos;
	undef %hpd;
	undef %group_hNos;
	#populate %rarehaplos, LocusID initially 0
	#loop through haplotypes in %lhd{$sl} again
	for $haplo (keys %{$$lhd_r{$sl}}){
		#determine haplotype-depth
		$dep = $$lhd_r{$sl}{$haplo};
		if ($dep < $hapdep) {#if depth < hapdep include in %rarehaplos
			$rarehaplos{$hapNo}[0] = $haplo;#sequence at pos 0
			$rarehaplos{$hapNo}[1] = $dep;#depth at pos 1
			$rarehaplos{$hapNo}[2] = 0;#locusID at pos 2, initially 0
			++$hapNo;#increment hapNo for next haplotype
			delete $$lhd_r{$sl}{$haplo};#delete haplotype from %lhd
		}
	}	
	#call sub rare_match_loc to assign rare haplotypes to existing loci
	#sub rare_match_loc adds rare haplotypes to %l_loc_cat_h_f{$sl}
	#and deletes rare haplotypes from %rarehaplos
	rare_match_loc($rarehapdist,$l_loc_cat_h_f_r,$sl,\%rarehaplos);
}

#definition of subroutine splitvar
#expects
#$sl: sequence length
#$hapdist: max dist between neighboring haplotypes in a locus
#determines variables for splitting of a sequence into nonoverlapping fragments:
#$hapdist+1 fragments
#returns 2d-@splitvar:[0]: array of startpositions [1] array of fragment-lengths
#the last fragment has a greater length than the others
#if the seqlength cannot be divided into $hapdist+1 fragments
#uses int for rounding down- should be no problem here
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

#definition of subroutine fno_f_h
#splits haplotypes into nonoverlapping fragments
#builds a hash of hashes %fno_f_h:
# {frag_No}{fragment}= @of haplotype-numbers containing that fragment
#expects:
#ref to %goodhaplos: {hapNo}=(sequence,depth,groupNo)
#ref to @splitvar: 2d: [0]: array of startpositions [1] array of fragment-lengths
#ref to %fno_f_h: to populate it
#$sl: sequence length
sub fno_f_h {
	my($goodhaplos_r,$splitvar_r,$fno_f_h_r,$sl) = @_;
	my @splitvar = @{$splitvar_r};#copy of array splitvar
	my $fno = 0;#current fragment number
	my $hapNo = 0;#current haplotype number
	my $haplo = '';#current haplotype
	my $frag = '';#current fragment
	my @fragper = ();#all possible permutations of a fragment containing N
					#permutation: replace N with A/C/G/T
	my $i = 0;
		
	#loop through the haplotype-numbers (keys) in %goodhaplos
	for $hapNo (keys %{$goodhaplos_r}) {
		$haplo = $$goodhaplos_r{$hapNo}[0];#get sequence for this haplotype number
		if ($haplo =~ /N/) {#if the haplotype contains N:
			for ($fno = 0; $fno < @{$splitvar[0]}; ++$fno) {#loop through fragment numbers for this haplo
				#build current fragment
				$frag = substr($haplo,${$splitvar[0]}[$fno],${$splitvar[1]}[$fno]);
				if ($frag =~ /N/) {#if the fragment contains N
					#call subroutine nper to build all possible permutations of the fragment
					@fragper = nper($frag,$sl);
					for ($i = 0; $i < @fragper; ++$i) {#loop through the permutations
						#add haplotype-number to %fno_f_h
						#with key1: $fno, key2: $fragper[$i] (current permutation)
						push @{$$fno_f_h_r{$fno}{$fragper[$i]}}, $hapNo;
					}					
				} else {#if the fragment does not contain N
					#add haplotype-number to %fno_f_h with $fno as key1, $frag as key2
					push @{$$fno_f_h_r{$fno}{$frag}}, $hapNo;
				}
			}
		} else {#if haplotype does not contain N:
			for ($fno = 0; $fno < @{$splitvar[0]}; ++$fno) {#loop through fragment numbers for this haplo
				#build current fragment
				$frag = substr($haplo,${$splitvar[0]}[$fno],${$splitvar[1]}[$fno]);
				#add haplotype-number to %fno_f_h with $fno as key1, $frag as key2
				push @{$$fno_f_h_r{$fno}{$frag}}, $hapNo;			
			}
		}
	}
}

#definition of subroutine nper
#creates all possible permutations of a DNA-sequence containing at least 1 N.
#Permutations: N replaced by A/C/G/T
#expects:
#$seq	the sequence
#$seql	its length
#returns @per: array containing all permutations
#knows chars: ACGTN (all capital only)
sub nper {
	#declare and initialize
	my ($seq, $seql) = @_;
	my $Nstring = '';#consists of $seql 'N'
	my $Nmatch = '';#contains \0 at all positions where $seq has N
					#defined values at other positions
	my @per = ();#will hold the permutations
	my $nseq = 0;#number of sequences currently in @per
	my $newper = '';#one sequence that is currently permuted	
	my $pos = 0;#current sequence position
	my $i  = 0;
	
	$Nstring = "N" x $seql;#create string of seqlength x N
	$Nmatch = $seq ^ $Nstring;#see above
	$seq =~ tr/N/A/;#transliterate all N in seq by A
	$per[0] = $seq;#put into @per as first permutation
	#loop through positions
	for ($pos = 0; $pos < $seql; ++$pos) {
		#if this is an N-position
		if (substr($Nmatch,$pos,1) =~ /\0/) {
			#determine number of seqs currently in @per
			$nseq = @per;
			#loop through these seqs
			for ($i = 0; $i < $nseq; ++$i) {
				$newper = $per[$i];#copy current seq
				substr ($newper,$pos,1) = 'C';#replace current pos with C
				push @per, $newper;#put the new permutation into @per
				substr ($newper,$pos,1) = 'G';#replace current pos with G
				push @per, $newper;#put the new permutation into @per
				substr ($newper,$pos,1) = 'T';#replace current pos with T
				push @per, $newper;#put the new permutation into @per
			}
		}
	}
	return @per;
}

#definition of subroutine hpd
#identifies pairs of haplotypes with a distance <= a threshold ($hapdep)
#populates %hpd
#expects:
#ref to %fno_f_h: {frag_No}{fragment}= @of haplotype-numbers
#ref to %goodhaplos: {hapNo}=(sequence,depth,groupNo)
#$hapdist: max dist between neighboring haplotypes in a locus
#$sl: sequence length (all sequences that will be compared have same length)
#ref to %hpd: {haplo1_No}{haplo2_No}=distance; to populate it
#all haplotype-numbers in one @ in %fno_f_h are candidates for pairwise comparison
sub hpd {
	#declare and initialize
	my($fno_f_h_r,$goodhaplos_r,$hapdist,$sl,$hpd_r) = @_;
	my $fragNo = 0;#current fragment number
	my $frag = '';#current fragment
	my @candidates = ();#copy of one value array in %fno_f_h
	my $ncandidates = 0;#number of candidates
	my $hapNo1 = 0;#number of first haplotype in a pair
	my $hapNo2 = 0;#number of second haplotype in a pair
	my $haplo1 = 0;#sequence of first haplotype in a pair
	my $haplo2 = 0;#sequence of second haplotype in a pair
	my $haplo1hasN = 0;#1 if the sequence contains N, else 0
	my $haplo2hasN = 0;#1 if the sequence contains N, else 0
	my $dist = 0;#number of differences between 2 sequences
	my $i = 0;#counter
	my $j = 0;#counter
	
	#loop through fragment-numbers (key1 in %fno_f_h)
	for $fragNo (keys %{$fno_f_h_r}) {
		#loop through the fragments (key2 in %fno_f_h)
		for $frag (keys %{$$fno_f_h_r{$fragNo}}) {
			#get the candidates
			@candidates = @{$$fno_f_h_r{$fragNo}{$frag}};
			$ncandidates = @candidates;
			#if there is more than one candidate
			if ($ncandidates > 1) {
				#loop through the candidates and form pairs
				for ($i = 0; $i < ($ncandidates - 1); ++$i) {
					for ($j = ($i + 1); $j < $ncandidates; ++$j) {
						#get a pair of candidate-numbers
						$hapNo1 = $candidates[$i];
						$hapNo2 = $candidates[$j];
						#if this pair has not yet been seen
						unless ((defined $$hpd_r{$hapNo1}) and (defined $$hpd_r{$hapNo1}{$hapNo2})){
							#get the corresponding sequences
							$haplo1 = $$goodhaplos_r{$hapNo1}[0];
							$haplo2 = $$goodhaplos_r{$hapNo2}[0];
							#check each sequence if it contains N
							if ($haplo1 =~ /N/) {
								$haplo1hasN = 1;
							} else {
								$haplo1hasN = 0;
							}
							if ($haplo2 =~ /N/) {
								$haplo2hasN = 1;
							} else {
								$haplo2hasN = 0;
							}
							#call sub distpwdN to determine the distance between them
							#with pw deletion (ignorance) of N
							$dist = distpwdN($haplo1,$haplo2,$haplo1hasN,$haplo2hasN,$sl);
							#if $dist <= $hapdist, add the pair to %hpd
							#key1: hapNo1, key2: hapNo2, value: $dist
							if ($dist <= $hapdist) {
								$$hpd_r{$hapNo1}{$hapNo2} = $dist;
							}
						}
					}
				}				
			}
		}
	}
}

#definition of subroutine distpwdN version 04
#calculates number of differences between 2 DNA sequences
#with pairwise deletion of positions containing N
#knows chars: (uppercase only): ACGTN
#expects 5 scalars
#$seq1, $seq2
#$seq1hasN, $seq2hasN: 1 if the seq contains N, 0 if not
#$seql: seqlength (assumed to be equal in both seqs
#returns $dist
sub distpwdN {
	#declare and initialize
	my ($seq1,$seq2,$seq1hasN,$seq2hasN,$seql) = @_;
	my $dist = 0;#the result: number of differences
	my $Nstring = '';#a string of N of length seqlength
	my $NAseq1 = '';#seq1 with CGT replaced by AAA
	my $NAseq2 = '';#seq2 with CGT replaced by AAA
	my $changepos = '';#string that contains defined values at positions
						#where only one of the 2 seqs has N, other pos: \0
	my$pos = 0;#a position;
	
	#if none of the sequences contains N: simple algorithm
	if ($seq1hasN == 0 and $seq2hasN == 0) {
		$dist = $seql - (($seq1 ^ $seq2) =~ tr/\0//);
	} else {#else: algorithm with pw del of N
		$Nstring = "N" x $seql;#create string of seqlength x N
		#create the 2 NAseqs
		($NAseq1 = $seq1) =~ tr/CGT/AAA/;
		($NAseq2 = $seq2) =~ tr/CGT/AAA/;
		#determine N-positions in seq1 and seq2 as \0 others as defined
		#match these against each other to get $changepos
		#:pos with N in only 1 seq: defined, others \0
		$changepos = (($NAseq1 ^ $Nstring) ^ ($NAseq2 ^ $Nstring));
		#loop through $changepos, for each changepos, replace corresponding
		#character in seq1 and seq2 with N
		while ($changepos =~ /[^\0]/g) {
			$pos = sprintf "%d", pos($changepos)-1;
			substr($seq1,$pos,1) = 'N';
			substr($seq2,$pos,1) = 'N';
		}
		#calculate the distance
		$dist = $seql - (($seq1 ^ $seq2) =~ tr/\0//);
	}
	return $dist;	
}

#definition of subroutine al
#assembles all good haplotypes into loci
#expects:
#ref to %goodhaplos: {hapNo}=(sequence,depth,groupNo)
	#all good haplos
	#groupNo is still 0 for all good haplos
#ref to %hpd: {haplo1-No}{haplo2-No}=distance
	#all pairwise comparisons of all haplos that HAVE a closely related haplo
	#each pair occurs only once
#ref to %group_hNos: {groupNo}=@ of haplo-numbers, still empty
#ref to %l_loc_cat_h_f: {seqlength}{locusID}{cat}{haplotype}=haplo_dep
	#will be populated as final product
	#owned by sub ind
#$sl: sequence length: all good haplos have same length
	#needed as key1 for %l_loc_cat_h_f
#ref to $nextlocID: locID to start with for the first locus built in this set
	#owned by sub ind
sub al {
	#declare and initialize: _r means a reference
	my ($goodhaplos_r,$hpd_r,$group_hNos_r,$l_loc_cat_h_f_r,$sl,$nextlocID_r) = @_;
	my $hap1No = 0;#number of first haplo of a pair
	my $hap2No = 0;#number of second haplo of a pair
	my $hap1group = 0;#current group of haplo1
	my $hap2group = 0;#current group of haplo2
	my $groupNo = 1;#temporary groupNo
	my $locID = $$nextlocID_r;#final locus-ID, initialized with starting ID for this set
	my @temparr1 = ();
	my @temparr2 = ();
	my $i = 0;
	
	#start assembling loci with those loci that are partners of a pair
	#loop through the pairs in %hpd
	for $hap1No (keys %{$hpd_r}) {
		for $hap2No (keys %{$$hpd_r{$hap1No}}) {
			#check the current group-numbers of these haplos
			$hap1group = $$goodhaplos_r{$hap1No}[2];
			$hap2group = $$goodhaplos_r{$hap2No}[2];
			#if both are not yet in a group (group-No 0)
			if ($hap1group == 0 and $hap2group == 0) {
				#make up a new group for them
				@{$$group_hNos_r{$groupNo}} = ($hap1No,$hap2No);
				$$goodhaplos_r{$hap1No}[2] = $groupNo;
				$$goodhaplos_r{$hap2No}[2] = $groupNo;
				++$groupNo;
			}
			#if haplo1 already is in a group but haplo2 is not
			elsif ($hap1group != 0 and $hap2group == 0) {
				#put haplo 2 into group of haplo1
				push @{$$group_hNos_r{$hap1group}}, $hap2No;
				$$goodhaplos_r{$hap2No}[2] = $hap1group;
			}
			#if haplo1 is not yet in a group but haplo2 is
			elsif ($hap1group == 0 and $hap2group != 0) {
				#put haplo 1 into group of haplo2
				push @{$$group_hNos_r{$hap2group}}, $hap1No;
				$$goodhaplos_r{$hap1No}[2] = $hap2group;
			}
			#if both are in different groups
			elsif ($hap1group != $hap2group) {
				#join their groups:
				#get all members of the group hap2 is in
				@temparr1 = @{$$group_hNos_r{$hap2group}};
				#loop through them and change their group-numbers to the number of hap1
				for ($i = 0; $i < @temparr1; ++$i) {
					$$goodhaplos_r{$temparr1[$i]}[2] = $hap1group;
				}
				#add them all to the group of hap1
				push @{$$group_hNos_r{$hap1group}}, @temparr1;
				#delete the former group of hap2
				delete $$group_hNos_r{$hap2group};
			}			
		}
	}
	#Continue with haplos in %goodhaplos that have not yet a group-number now.
	#Make up new, homozygous loci for them
	#loop through %goodhaplos
	for $hap1No (keys %{$goodhaplos_r}) {
		#if the locus has not yet a groupNo (groupNo 0)
		if ($$goodhaplos_r{$hap1No}[2] == 0) {
			#give a group-number to this haplotype-number
			$$goodhaplos_r{$hap1No}[2] = $groupNo;
			#build a corresponding group for it
			$$group_hNos_r{$groupNo}[0] = $hap1No;
			++$groupNo;#increment group number for next group
		}
	}
	#All good haplos have a group-number and are in a group now.
	#Populate %l_loc_cat_h_f now, use new consecutive locID instead of the former group-numbers
	#loop through the groups in %group_nNos_r
	for $groupNo (keys %{$group_hNos_r}) {
		#get numbers of all members of the group
		@temparr1 = @{$$group_hNos_r{$groupNo}};
		#loop through the members of the group
		for ($i = 0; $i < @temparr1; ++$i) {
			#get all data about current member out of %goodhaplos
			@temparr2 = @{$$goodhaplos_r{$temparr1[$i]}};
			#data in @temparr2: [0]sequence,[1]depth[2]groupNo
			#delete this haplotype from %goodhaplos
			delete $$goodhaplos_r{$temparr1[$i]};
			#Add this haplotype to %l_loc_cat_h_f
			#key1: $sl, $key2: $locID, key3: "good" key4: $temparr2[0], value: $temparr2[1]
			$$l_loc_cat_h_f_r{$sl}{$locID}{"good"}{$temparr2[0]} = $temparr2[1];
		}
		++$locID;#increment and next group
	}
	#finally determine $nextlocID owned by sub ind for the next set of loci if any
	$$nextlocID_r = $locID;
}

#definition of subroutine rare_match_loc
#assigns rare haplotypes to loci
#stores them in %l_loc_cat_h_f under cat="rare"
#expects:
#$rarehapdist (max dist between a rare haplotype and at least one good haplotype
		#in a locus to allow for assignment of the rare haplotype to the locus)
#ref to %l_loc_cat_h_f
	# {seqlength}{locusID}{cat}{haplotype}=haplotype_depth
#$sl:		seqlength as key to %l_loc_cat_h_f
#ref to %rarehaplos: {hapNo}=(sequence,depth,locusID), locusID still 0
sub rare_match_loc {
	#declare and initialize: _r means a reference
	my ($rarehapdist,$l_loc_cat_h_f_r,$sl,$rarehaplos_r) = @_;
	my @splitvar = ();#variables for splitting of haplotypes into nonoverlapping fragments:
						#2d array: [0]: array of startpositions [1] array of fragment-lengths
	my %fragCat = ();#fragment catalog of existing loci:
					# {fragment number}{fragment}{locID}=nhap, nhap: simple counter, number of haplos of this locus
					#that have this frag, there must be some value
	my $locID = 0;#a locusID
	my $goodhap = '';#sequence of a good haplotype
	my $goodhaphasN = 0;#1 if good haplotype contains N, else 0
	my $rarehap = '';#sequence of a rare haplotype
	my $rarehapdep = 0;#depth of a rare haplotype
	my $rarehaphasN = 0;#1 if rare haplotype contains N, else 0
	my $fno = 0;#fragment number
	my $frag = '';#current fragment of a haplotype
	my @fragper = ();#all possible permutations of a fragment containing N
					#permutation: replace N with A/C/G/T
	my $hapNo = 0;#number of current haplotype
	my %candloc = ();#candidate loci for a current rare haplotype
					#key: locID, value: a scalar (number of times this locus is found in catalog)
	my $candlocID = 0;#locusID of a candidate locus
	my $dist = 0;#distance (number of differences) between two sequences
	my $matchlocID = 0;#ID of matching locus
	my $nmatchloc = 0;#number of matching loci
	my $i = 0;
		
	#call subroutine splitvar to determine variables for splitting
	@splitvar = splitvar($sl,$rarehapdist);
	
	#build the fragment catalog of existing loci
	#loop through locusIDs in %l_loc_cat_h_f
	for $locID (keys %{$$l_loc_cat_h_f_r{$sl}}) {
		#loop through the good haplotypes of this locusID
		for $goodhap (keys %{$$l_loc_cat_h_f_r{$sl}{$locID}{"good"}}) {
			if ($goodhap =~ /N/) {#if the good haplotype contains N
				#loop trough fragment numbers
				for ($fno = 0; $fno < @{$splitvar[0]}; ++$fno) {
					#build current fragment
					$frag = substr($goodhap,${$splitvar[0]}[$fno],${$splitvar[1]}[$fno]);
					if ($frag =~ /N/) {#if the fragment contains N
						#call subroutine nper to build
						#all possible permutations of the fragment
						@fragper = nper($frag,$sl);
						#loop through the permutations
						for ($i = 0; $i < @fragper; ++$i) {
							#store current permutation in the catalog
							++$fragCat{$fno}{$fragper[$i]}{$locID};
						}						
					} else {#if the fragment doesn't contain N
						#store in the catalog
						++$fragCat{$fno}{$frag}{$locID};
					}
				}
			} else {#if the good haplotype doesn't contain N
				#loop trough fragment numbers
				for ($fno = 0; $fno < @{$splitvar[0]}; ++$fno) {
					#build current fragment
					$frag = substr($goodhap,${$splitvar[0]}[$fno],${$splitvar[1]}[$fno]);
					#store in the catalog
					++$fragCat{$fno}{$frag}{$locID};
				}
			}
		}
	}
	
	##########################################
	#match rare haplotypes against the catalog
	##########################################
	
	#loop through haplotype numbers in %rarehaplos
	for $hapNo (keys %{$rarehaplos_r}) {
		#get the rare haplotype and its depths
		$rarehap = $$rarehaplos_r{$hapNo}[0];
		$rarehapdep = $$rarehaplos_r{$hapNo}[1];
		if ($rarehap =~ /N/) {#if rare haplotype contains N
			$rarehaphasN = 1;
			#loop through fragment numbers
			for ($fno = 0; $fno < @{$splitvar[0]}; ++$fno) {
				#build current fragment
				$frag = substr($rarehap,${$splitvar[0]}[$fno],${$splitvar[1]}[$fno]);
				if ($frag =~ /N/) {#if the fragment contains N
					#call subroutine nper to build
					#all possible permutations of the fragment
					@fragper = nper($frag,$sl);
					#loop through the permutations
					for ($i = 0; $i < @fragper; ++$i) {
						#if you find this fragment permutation under this fragment number in the catalog:
						if (exists $fragCat{$fno}{$fragper[$i]}) {
							#loop through the locusIDs under this fragment-permutation
							for $candlocID (keys %{$fragCat{$fno}{$fragper[$i]}}) {
								#store them in %candloc
								++$candloc{$candlocID};
							}
						}
					}					
				} else {#if the fragment doesn't contain N
					#if you find this fragment under this fragment number in the catalog:
					if (exists $fragCat{$fno}{$frag}) {
						#loop through the locusIDs under this fragment
						for $candlocID (keys %{$fragCat{$fno}{$frag}}) {
							#store them in %candloc
							++$candloc{$candlocID};
						}				
					}					
				}				
			}			
		} else {#if rare haplotype doesn't contain N
			$rarehaphasN = 0;
			#loop through fragment numbers
			for ($fno = 0; $fno < @{$splitvar[0]}; ++$fno) {
				#build current fragment
				$frag = substr($rarehap,${$splitvar[0]}[$fno],${$splitvar[1]}[$fno]);
				#if you find this fragment under this fragment number in the catalog:
				if (exists $fragCat{$fno}{$frag}) {
					#loop through the locusIDs under this fragment
					for $candlocID (keys %{$fragCat{$fno}{$frag}}) {
						#store them in %candloc
						++$candloc{$candlocID};
					}				
				}
			}
		}
		#loop through the locusIDs in %candloc
		for $candlocID (keys %candloc) {
			#loop through the haplotypes of current candidate locus
			#in %l_loc_cat_h_f
			for $goodhap (keys %{$$l_loc_cat_h_f_r{$sl}{$candlocID}{"good"}}) {
				#check if $goodhap contains N:
				if ($goodhap =~ /N/) {
					$goodhaphasN = 1;
				} else {
					$goodhaphasN = 0;
				}
				#determine distance between $rarehap and $goodhap
				$dist = distpwdN($rarehap,$goodhap,$rarehaphasN,$goodhaphasN,$sl);
				#if the distance is tolerable (<= $rarehapdist)
				if($dist <= $rarehapdist) {
					#store the current locusID as matching locus
					$matchlocID = $candlocID;
					#count a matching locus
					++$nmatchloc;
					#don't look at the other haplotypes of this locus
					last;
				}
			}
			#if the number of matching loci is now >1, stop
			if($nmatchloc > 1) {
				last;
			}
		}
		#if you found one and only one matching locus
		if ($nmatchloc == 1) {
			#store the rare haplotype and its depth in %l_loc_cat_h_f
			$$l_loc_cat_h_f_r{$sl}{$matchlocID}{"rare"}{$rarehap} = $rarehapdep;
		}
		#delete rare haplotype from %rarehaplos
		delete $$rarehaplos_r{$hapNo};
		#set variables back
		undef %candloc;
		$matchlocID = 0;
		$nmatchloc = 0;
		#go to next rare haplotype
	}	
}

#definition of subroutine snphapcal1
#calls SNPs with frequency threshold method, then haplotypes
#prints output to files
#*_loci.txt
#*_alleles.txt
#*_svars.txt
#*_reads.txt
#expects:
#ref to %l_loc_cat_h_f: {seqlength}{locusID}{cat}{haplotype}= haplotype-depth
#$ind: ID of current individual for building of outfilestem
#$min_cf: minimum character frequency in an individual locus alignment
#$min_af: minimum allele frequency in an individual locus alignment
#$fqfile: name of fastq infile for current individual
#$build_rc: 1: fastq-infile contains only forward reads of rc-duplicates
			#all headers end with _1
			#sub ind has built reverse reads with same header ending with _2 as reverse complement
			#0: Not so
#$z: 1: write svars.txt and reads.txt as gzip compressed files svars.txt.gz reads.txt.gz; 0: don't
sub snphapcal1 {
	#declare and initialize, _r means a reference
	my ($l_loc_cat_h_f_r,$ind,$min_cf,$min_af,$fqfile,$build_rc,$z) = @_;
	my $sl = 0;#seqlength
	my $locID = 0;#locus ID
	my $cat = '';#"good" or "rare"
	my $hapseq = '';#sequence of a haplotype
	my $haprd = '';#read-depth of a haplotype
	my @catinarr = ();#cat of all haplotypes in inputorder
	my @hapseqinarr = ();#seq of all haplotypes in inputorder
	my @haprdinarr = ();#read-depth of all haplotypes in inputorder
	my @hapIDinarr = ();#IDs of all haplotypes in inputorder, created here just for output
	my $hapNo = 0;
	my $hapID = 1;#first haplotype ID
	my @seqmatin = ();#d1: rows, d2: cols original sequences of a locus
	my @varposin = ();#variable positions in input-sequences
	my @seqmatcorr = ();#copy of @seqmatin, sequencing errors corrected or changed to N
	my @temparr1 = ();
	my @dischaplos = ();#numbers of discarded haplos (input-order-numbers)
	my %hapgroups = ();# {groupNo}=@ of hapNo, groups will later be collapsed to alleles
	my %hapgroup_dep = ();# {groupNo} = depth of each hapgroup
	my $locdep = 0;#total depth of a locus, including discarded reads
	my $loclostall = 0;#number of alleles a locus lost
						#through discard of ambiguous haplotypes, alleles of low frequency
	my $locpotlostmoreall = 0;#1 if the number of lost alleles could be greater than $loclostall
							#0 if not							
	my %alleles = ();# {allele_ID} = @: [0]: allele_seq [1]: depth, remainder: hapNos
	my $next_allele_ID = 1;#ID for the next allele that will be created
	my $allele_ID = 0;#an allele ID
	my @alleledata = ();#value array of one allele in %alleles
	my $n_alleles = 0;#number of alleles in a locus
	my $loc_cat = 'valid';#"lost" if the locus is lost, i.e. no valid alleles can be determined
	my $seqID = '';#a sequence ID from fastq-infile
	my $out_loci = '';#name of loci-outfile
	my $out_alleles = '';#name of alleles-outfile
	my $out_haplotypes = '';#name of svars-outfile
	my $out_reads = '';#name of reads-outfile
	my $out_haplotypes_h = '';#handle for svars-outfile
	my $out_reads_h = '';#handle for reads-outfile
	my %hdat = ();# {hapseq}= string: LocID\tLoc_cat\tseql\tAll_ID\thapID\thapcat
					#this hash contains data for the reads-outfile with haplotype sequence as key
	my $fqfileh = '';#filehandle for fastq infile
	my $seqline1 = '';#line1 of a sequence entry in fastq infile
	my $seqline2 = '';#line2 of a sequence entry in fastq infile
	my $seqline3 = '';#line3 of a sequence entry in fastq infile
	my $seqline4 = '';#line4 of a sequence entry in fastq infile
	my $rseqline1 = '';#header of new reverse compl. of seq entry
	my $rseqline2 = '';#new reverse compl. of seq entry
	my $success = "analysis aborted\n";#return value to sub ind, if succesful: "analysis completed"
									#else: failure message
	my $i = 0;
	my $j = 0;
	
	#create outfilenames
	$out_loci = $ind . '_loci.txt';
	$out_alleles = $ind . '_alleles.txt';
	if ($z == 0) {#print as plain text
		$out_haplotypes = $ind . '_svars.txt';
		$out_reads = $ind . '_reads.txt';
	} else {#print as gzip compressed files
		$out_haplotypes = $ind . '_svars.txt.gz';
		$out_reads = $ind . '_reads.txt.gz';
	}
	#open outfiles and fastq infile
	unless(open(OUTLOCI, ">$out_loci")) {
		$success = "Could not open $out_loci. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_loci. Analysis aborted.\n";#failure message to screen
		return $success;#return failure message to sub ind
	}
	unless(open(OUTALL, ">$out_alleles")) {
		$success = "Could not open $out_alleles. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_alleles. Analysis aborted.\n";#failure message to screen
		close OUTLOCI;#close outfile already open
		return $success;#return failure message to sub ind
	}
	if ($z == 0) {#print svars outfile as plain text file
		$out_haplotypes_h = IO::File->new(">$out_haplotypes");
	} else {#print as gzip compressed file
		$out_haplotypes_h = new IO::Compress::Gzip $out_haplotypes ,-Level=>Z_BEST_SPEED;
	}
	unless(defined $out_haplotypes_h) {
		$success = "Could not open $out_haplotypes. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_haplotypes. Analysis aborted.\n";#failure message to screen
		close OUTLOCI;#close outfile already open
		close OUTALL;#close outfile already open
		return $success;#return failure message to sub ind
	}
	if ($z == 0) {#print reads outfile as plain text file
		$out_reads_h = IO::File->new(">$out_reads");
	} else {#print as gzip compressed file
		$out_reads_h = new IO::Compress::Gzip $out_reads ,-Level=>Z_BEST_SPEED;
	}
	unless(defined $out_reads_h) {
		$success = "Could not open $out_reads. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_reads. Analysis aborted.\n";#failure message to screen
		close OUTLOCI;#close outfile already open
		close OUTALL;#close outfile already open
		close $out_haplotypes_h;#close outfile already open
		undef $out_haplotypes_h;
		return $success;#return failure message to sub ind
	}
	if ($fqfile =~ /.gz$/) {#open as gzip compressed file
		$fqfileh = new IO::Uncompress::Gunzip $fqfile;
	} else {#open as text file
		$fqfileh = IO::File->new("< $fqfile");
	}
	unless(defined $fqfileh) {
		$success = "Could not open $fqfile. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $fqfile. Analysis aborted.\n";
		close OUTLOCI;#close outfile already open
		close OUTALL;#close outfile already open
		close $out_haplotypes_h;#close outfile already open
		undef $out_haplotypes_h;
		close $out_reads_h;#close outfile already open
		undef $out_reads_h;
		return $success;#return failure message to sub ind
	}
	#print headerlines to outfiles
	print OUTLOCI "LocID\tLoc_cat\tseql\tLoc_dep\tn_alleles\tlost_alleles\tpotlostmore\n";
	print OUTALL "LocID\tLoc_cat\tseql\tAll_ID\tAllseq\tAlldep\n";
	print $out_haplotypes_h "LocID\tLoc_cat\tseql\tAll_ID\tsvarID\tsvarcat\tsvarseq\tsvardep\n";
	print $out_reads_h "LocID\tLoc_cat\tseql\tAll_ID\tsvarID\tsvarcat\treadID\n";
	
	#loop through sequence lengths
	for $sl (keys %{$l_loc_cat_h_f_r}) {
		#loop through locus IDs
		for $locID (keys %{$$l_loc_cat_h_f_r{$sl}}) {
			#Collect all data about current locus in
			#@catinarr, @hapseqinarr, @haprdinarr, @hapIDinarr
			#loop through good haplotypes and collect data
			$cat = 'good';
			while (($hapseq,$haprd) = each %{$$l_loc_cat_h_f_r{$sl}{$locID}{$cat}}) {
				push @catinarr, $cat;
				push @hapseqinarr, $hapseq;
				push @haprdinarr, $haprd;
				push @hapIDinarr, $hapID;
				++$hapID;#increment for next haplotype
				delete $$l_loc_cat_h_f_r{$sl}{$locID}{$cat}{$hapseq};#delete good haplotype from %l_loc_cat_h_f
			}
			$cat = 'rare';
			#if there are also rare haplotypes in current locus:
			if (defined $$l_loc_cat_h_f_r{$sl}{$locID}{$cat}) {
				#loop through rare haplotypes and collect data
				while (($hapseq,$haprd) = each %{$$l_loc_cat_h_f_r{$sl}{$locID}{$cat}}) {
					push @catinarr, $cat;
					push @hapseqinarr, $hapseq;
					push @haprdinarr, $haprd;
					push @hapIDinarr, $hapID;
					++$hapID;#increment for next haplotype
					delete $$l_loc_cat_h_f_r{$sl}{$locID}{$cat}{$hapseq};#delete rare haplotype from %l_loc_cat_h_f
				}
			}
			#calculate total depth of locus $locdep
			for ($i = 0; $i < @haprdinarr; ++$i) {
				$locdep += $haprdinarr[$i];
			}
			#if there is only one haplotype:
			if (@hapseqinarr == 1) {
				$hapseq = $hapseqinarr[0]; #get the sequence
				$haprd = $haprdinarr[0];#get its depth
				if ($hapseq =~ /N/) {#if the sequence contains N
					$loc_cat = 'lost';#record that the locus is lost
					$loclostall = 1;#one allele was lost from the locus
					$dischaplos[0] = 0;#add hapNo to dischaplos
				} else {#no N: locus valid
					#create one allele and count it
					@{$alleles{$next_allele_ID}} = ($hapseq,$haprd,0);
					$n_alleles = 1;
					++$next_allele_ID;#increment for next allele in next locus
				}			
			}
			#if there is more than 1 haplotype
			if (@hapseqinarr > 1) {
				#build @seqmatin:d1: rows, d2: cols
				for ($i = 0; $i < @hapseqinarr; ++$i) {
					@temparr1 = split('', $hapseqinarr[$i]);
					push @seqmatin, [@temparr1];
				}
				#call sub varpos to get position numbers of variable positions into @varposin
				@varposin = varpos(\@seqmatin);
				#call sub errbase1 to correct sequencing errors
				#corrected sequences into: @seqmatcorr
				@seqmatcorr = errbase1(\@seqmatin,\@varposin,\@haprdinarr,$locdep,$min_cf);
				#call sub hapgroups to build haplotype-groups into %hapgroups
				#determine their depths into %hapgroup_dep
				#collect numbers of discarded haplos in @dischaplos
				#record how many alleles got lost in that procedure($loclostall)
				#and if the number of lost alleles could be greater($locpotlostmoreall == 1)
				#haplotype groups will later be collapsed into alleles
				hapgroups(\@seqmatcorr,\%hapgroups,\%hapgroup_dep,
				\@dischaplos,$sl,\@haprdinarr,\$loclostall,\$locpotlostmoreall);
				#if there is no hapgroup:
				if (keys %hapgroups == 0) {
					$loc_cat = 'lost';#record that the locus is lost
				} else {
					#call sub build_alleles to populate %alleles and determine $n_alleles
					build_alleles(\@seqmatcorr,\@varposin,\%hapgroups,\%hapgroup_dep,
					$locdep,\%alleles,\$n_alleles,\@dischaplos,\$loclostall,\$next_allele_ID,$min_af);
					#if there is no allele after building alleles
					if (keys %alleles == 0) {
						$loc_cat = 'lost';#record that the locus is lost
					}
				}
			}
			#print output for this locus to loci-outfile
			print OUTLOCI "$locID\t$loc_cat\t$sl\t$locdep\t$n_alleles\t$loclostall\t$locpotlostmoreall\n";
			#print output for this locus to alleles-outfile
			if ($loc_cat eq 'valid') {#for a valid locus
				#loop through allele IDs of this locus
				for $allele_ID (sort {$a <=> $b} keys %alleles) {
					@alleledata = @{$alleles{$allele_ID}};#get data about this allele
					print OUTALL "$locID\t$loc_cat\t$sl\t$allele_ID\t$alleledata[0]\t$alleledata[1]\n";
				}
			} else {#for a lost locus
				print OUTALL "$locID\t$loc_cat\t$sl\t0\t0\t0\n";
			}
			#print output for this locus to svars-outfile
			if ($loc_cat eq 'valid') {#if locus is valid, i.e. contains at least one valid allele:
				#loop through allele IDs of this locus
				for $allele_ID (sort {$a <=> $b} keys %alleles) {
					@alleledata = @{$alleles{$allele_ID}};#get data about this allele
					#loop through hapNos of this allele
					for ($i = 2; $i < @alleledata; ++$i) {
						$hapNo = $alleledata[$i];
						print $out_haplotypes_h "$locID\t$loc_cat\t$sl\t$allele_ID\t",
						"$hapIDinarr[$hapNo]\tused\t$hapseqinarr[$hapNo]\t$haprdinarr[$hapNo]\n";
					}
				}
			}
			#loop through discarded haplotypes from this locus
			for ($i = 0; $i < @dischaplos; ++$i) {
				$hapNo = $dischaplos[$i];
				print $out_haplotypes_h "$locID\t$loc_cat\t$sl\t0\t",
				"$hapIDinarr[$hapNo]\tdiscarded\t$hapseqinarr[$hapNo]\t$haprdinarr[$hapNo]\n";
			}
			#collect data about this locus for reads-outfile
			#into %hdat {hapseq}=string: LocID\tLoc_cat\tseql\tAll_ID\thapID\thapcat
			if ($loc_cat eq 'valid') {#if locus is valid, i.e. contains at least one valid allele:
				#loop through allele IDs of this locus
				for $allele_ID (sort {$a <=> $b} keys %alleles) {
					@alleledata = @{$alleles{$allele_ID}};#get data about this allele
					#loop through hapNos of this allele
					for ($i = 2; $i < @alleledata; ++$i) {
						$hapNo = $alleledata[$i];
						$hapseq = $hapseqinarr[$hapNo];#get seq of haplotype
						#store data about this haplotype in %hdat
						$hdat{$hapseq} = "$locID\t$loc_cat\t$sl\t$allele_ID\t$hapIDinarr[$hapNo]\tused";
					}
				}
			}
			#loop through discarded haplotypes from this locus
			for ($i = 0; $i < @dischaplos; ++$i) {
				$hapNo = $dischaplos[$i];
				$hapseq = $hapseqinarr[$hapNo];#get seq of haplotype
				#store data about this haplotype in %hdat
				$hdat{$hapseq} = "$locID\t$loc_cat\t$sl\t$allele_ID\t$hapIDinarr[$hapNo]\tdiscarded";
			}
			#set variables back
			@catinarr = ();
			@hapseqinarr = ();
			@haprdinarr = ();
			@seqmatin = ();
			%hapgroups = ();
			@dischaplos = ();
			$locdep = 0;
			%hapgroup_dep = ();
			$loclostall = 0;
			$locpotlostmoreall = 0;
			%alleles = ();
			$n_alleles = 0;
			$loc_cat = 'valid';
			@hapIDinarr = ();
			#go to next locus
		}
	}
	#read in fastq infile, look up sequence in %hdat, and print data to outfile OUTREADS
	if ($build_rc == 1) {#fastq infile contains only f-reads of rc-duplicates, build r-reads as reverse complement
		while ($seqline1 = <$fqfileh>) {
			$seqline2 = <$fqfileh>;
			$seqline3 = <$fqfileh>;
			$seqline4 = <$fqfileh>;
			chomp $seqline1;
			chomp $seqline2;
			$rseqline1 = $seqline1;
			$rseqline1 =~ s/_1$/_2/;
			$rseqline2 = reverse $seqline2;
			$rseqline2 =~ tr/ACGTacgt/TGCAtgca/;
			#look up f-sequence in %hdat
			if (defined $hdat{$seqline2}) {
				print $out_reads_h "$hdat{$seqline2}\t$seqline1\n";#print data to OUTREADS
			}
			#look up r-sequence in %hdat
			if (defined $hdat{$rseqline2}) {
				print $out_reads_h "$hdat{$rseqline2}\t$rseqline1\n";#print data to OUTREADS
			}
		}
	} else {#Do not build r-reads
		while ($seqline1 = <$fqfileh>) {
			$seqline2 = <$fqfileh>;
			$seqline3 = <$fqfileh>;
			$seqline4 = <$fqfileh>;
			chomp $seqline1;
			chomp $seqline2;
			#look up sequence in %hdat
			if (defined $hdat{$seqline2}) {
				print $out_reads_h "$hdat{$seqline2}\t$seqline1\n";#print data to OUTREADS
			}
		}
	}
	close OUTLOCI;
	close OUTALL;
	close $out_haplotypes_h;
	undef $out_haplotypes_h;
	close $out_reads_h;
	undef $out_reads_h;
	close $fqfileh;
	undef $fqfileh;
	$success = "analysis completed\n";#success message
	return $success;#return success message
}

#definition of subroutine snphapcalbin
#calls SNPs with binomial likelihood ratio method, then haplotypes
#prints output to files
#*_loci.txt
#*_alleles.txt
#*_svars.txt
#*_reads.txt
#expects:
#ref to %l_loc_cat_h_f: {seqlength}{locusID}{cat}{haplotype}= haplotype-depth
#$ind: ID of current individual for building of outfilestem
#$binerr: sequencing error rate
#$fqfile: name of fastq infile for current individual
#$build_rc: 1: fastq-infile contains only forward reads of rc-duplicates
			#all headers end with _1
			#sub ind has built reverse reads with same header ending with _2 as reverse complement
			#0: Not so
#$z: 1: write svars.txt and reads.txt as gzip compressed files svars.txt.gz reads.txt.gz; 0: don't
sub snphapcalbin {
	#declare and initialize, _r means a reference
	my ($l_loc_cat_h_f_r,$ind,$binerr,$fqfile,$build_rc,$z) = @_;
	my $sl = 0;#seqlength
	my $locID = 0;#locus ID
	my $cat = '';#"good" or "rare"
	my $hapseq = '';#sequence of a haplotype
	my $haprd = '';#read-depth of a haplotype
	my @catinarr = ();#cat of all haplotypes in inputorder
	my @hapseqinarr = ();#seq of all haplotypes in inputorder
	my @haprdinarr = ();#read-depth of all haplotypes in inputorder
	my @hapIDinarr = ();#IDs of all haplotypes in inputorder, created here just for output
	my $hapNo = 0;
	my $hapID = 1;#first haplotype ID
	my @seqmatin = ();#d1: rows, d2: cols, original sequences of a locus
	my @varposin = ();#variable positions in input-sequences
	my @seqmatcorr = ();#copy of @seqmatin, sequencing errors corrected or changed to N
	my %majormin = ();# {major} = minor
						#major: count of most frequent nucleotide at a variable position
						#minor: minimum count of second most frequent nucleotide to accept a SNP
	my @temparr1 = ();
	my @dischaplos = ();#numbers of discarded haplos (input-order-numbers)
	my %hapgroups = ();# {groupNo}=@ of hapNo, groups will later be collapsed to alleles
	my %hapgroup_dep = ();# {groupNo} = depth of each hapgroup
	my $locdep = 0;#total depth of a locus, including discarded reads
	my $loclostall = 0;#number of alleles a locus lost
						#through discard of ambiguous haplotypes, alleles of low frequency
	my $locpotlostmoreall = 0;#1 if the number of lost alleles could be greater than $loclostall
							#0 if not							
	my %alleles = ();# {allele_ID} = @: [0]: allele_seq [1]: depth, remainder: hapNos
	my $next_allele_ID = 1;#ID for the next allele that will be created
	my $allele_ID = 0;#an allele ID
	my @alleledata = ();#value array of one allele in %alleles
	my $n_alleles = 0;#number of alleles in a locus
	my $loc_cat = 'valid';#"lost" if the locus is lost, i.e. no valid alleles can be determined
	my $seqID = '';#a sequence ID from fastq-infile
	my $out_loci = '';#name of loci-outfile
	my $out_alleles = '';#name of alleles-outfile
	my $out_haplotypes = '';#name of svars-outfile
	my $out_reads = '';#name of reads-outfile
	my $out_haplotypes_h = '';#handle for svars-outfile
	my $out_reads_h = '';#handle for reads-outfile
	my %hdat = ();# {hapseq}= string: LocID\tLoc_cat\tseql\tAll_ID\thapID\thapcat
					#this hash contains data for the reads-outfile with haplotype sequence as key
	my $fqfileh = '';#filehandle for fastq infile
	my $seqline1 = '';#line1 of a sequence entry in fastq infile
	my $seqline2 = '';#line2 of a sequence entry in fastq infile
	my $seqline3 = '';#line3 of a sequence entry in fastq infile
	my $seqline4 = '';#line4 of a sequence entry in fastq infile
	my $rseqline1 = '';#header of new reverse compl. of seq entry
	my $rseqline2 = '';#new reverse compl. of seq entry
	my $success = "analysis aborted\n";#return value to sub ind, if succesful: "analysis completed"
									#else: failure message
	my $i = 0;
	my $j = 0;
	
	#create outfilenames
	$out_loci = $ind . '_loci.txt';
	$out_alleles = $ind . '_alleles.txt';
	if ($z == 0) {#print as plain text
		$out_haplotypes = $ind . '_svars.txt';
		$out_reads = $ind . '_reads.txt';
	} else {#print as gzip compressed files
		$out_haplotypes = $ind . '_svars.txt.gz';
		$out_reads = $ind . '_reads.txt.gz';
	}
	#open outfiles and fastq infile
	unless(open(OUTLOCI, ">$out_loci")) {
		$success = "Could not open $out_loci. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_loci. Analysis aborted.\n";#failure message to screen
		return $success;#return failure message to sub ind
	}
	unless(open(OUTALL, ">$out_alleles")) {
		$success = "Could not open $out_alleles. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_alleles. Analysis aborted.\n";#failure message to screen
		close OUTLOCI;#close outfile already open
		return $success;#return failure message to sub ind
	}
	if ($z == 0) {#print haplotypes outfile as plain text file
		$out_haplotypes_h = IO::File->new(">$out_haplotypes");
	} else {#print as gzip compressed file
		$out_haplotypes_h = new IO::Compress::Gzip $out_haplotypes ,-Level=>Z_BEST_SPEED;
	}
	unless(defined $out_haplotypes_h) {
		$success = "Could not open $out_haplotypes. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_haplotypes. Analysis aborted.\n";#failure message to screen
		close OUTLOCI;#close outfile already open
		close OUTALL;#close outfile already open
		return $success;#return failure message to sub ind
	}
	if ($z == 0) {#print reads outfile as plain text file
		$out_reads_h = IO::File->new(">$out_reads");
	} else {#print as gzip compressed file
		$out_reads_h = new IO::Compress::Gzip $out_reads ,-Level=>Z_BEST_SPEED;
	}
	unless(defined $out_reads_h) {
		$success = "Could not open $out_reads. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $out_reads. Analysis aborted.\n";#failure message to screen
		close OUTLOCI;#close outfile already open
		close OUTALL;#close outfile already open
		close $out_haplotypes_h;#close outfile already open
		undef $out_haplotypes_h;
		return $success;#return failure message to sub ind
	}
	if ($fqfile =~ /.gz$/) {#open as gzip compressed file
		$fqfileh = new IO::Uncompress::Gunzip $fqfile;
	} else {#open as text file
		$fqfileh = IO::File->new("< $fqfile");
	}
	unless(defined $fqfileh) {
		$success = "Could not open $fqfile. Analysis aborted.\n";#store failure message in $success
		print "$ind: Can't open $fqfile. Analysis aborted.\n";
		close OUTLOCI;#close outfile already open
		close OUTALL;#close outfile already open
		close $out_haplotypes_h;#close outfile already open
		undef $out_haplotypes_h;
		close $out_reads_h;#close outfile already open
		undef $out_reads_h;
		return $success;#return failure message to sub ind
	}
	#print headerlines to outfiles
	print OUTLOCI "LocID\tLoc_cat\tseql\tLoc_dep\tn_alleles\tlost_alleles\tpotlostmore\n";
	print OUTALL "LocID\tLoc_cat\tseql\tAll_ID\tAllseq\tAlldep\n";
	print $out_haplotypes_h "LocID\tLoc_cat\tseql\tAll_ID\tsvarID\tsvarcat\tsvarseq\tsvardep\n";
	print $out_reads_h "LocID\tLoc_cat\tseql\tAll_ID\tsvarID\tsvarcat\treadID\n";
	
	#loop through sequence lengths
	for $sl (keys %{$l_loc_cat_h_f_r}) {
		#loop through locus IDs
		for $locID (keys %{$$l_loc_cat_h_f_r{$sl}}) {
			#Collect all data about current locus in
			#@catinarr, @hapseqinarr, @haprdinarr, @hapIDinarr
			#loop through good haplotypes and collect data
			$cat = 'good';
			while (($hapseq,$haprd) = each %{$$l_loc_cat_h_f_r{$sl}{$locID}{$cat}}) {
				push @catinarr, $cat;
				push @hapseqinarr, $hapseq;
				push @haprdinarr, $haprd;
				push @hapIDinarr, $hapID;
				++$hapID;#increment for next haplotype
				delete $$l_loc_cat_h_f_r{$sl}{$locID}{$cat}{$hapseq};#delete good haplotype from %l_loc_cat_h_f
			}
			$cat = 'rare';
			#if there are also rare haplotypes in current locus:
			if (defined $$l_loc_cat_h_f_r{$sl}{$locID}{$cat}) {
				#loop through rare haplotypes and collect data
				while (($hapseq,$haprd) = each %{$$l_loc_cat_h_f_r{$sl}{$locID}{$cat}}) {
					push @catinarr, $cat;
					push @hapseqinarr, $hapseq;
					push @haprdinarr, $haprd;
					push @hapIDinarr, $hapID;
					++$hapID;#increment for next haplotype
					delete $$l_loc_cat_h_f_r{$sl}{$locID}{$cat}{$hapseq};#delete rare haplotype from %l_loc_cat_h_f
				}
			}
			#calculate total depth of locus $locdep
			for ($i = 0; $i < @haprdinarr; ++$i) {
				$locdep += $haprdinarr[$i];
			}
			#if there is only one haplotype:
			if (@hapseqinarr == 1) {
				$hapseq = $hapseqinarr[0]; #get the sequence
				$haprd = $haprdinarr[0];#get its depth
				if ($hapseq =~ /N/) {#if the sequence contains N
					$loc_cat = 'lost';#record that the locus is lost
					$loclostall = 1;#one allele was lost from the locus
					$dischaplos[0] = 0;#add hapNo to dischaplos
				} else {#no N: locus valid
					#create one allele and count it
					@{$alleles{$next_allele_ID}} = ($hapseq,$haprd,0);
					$n_alleles = 1;
					++$next_allele_ID;#increment for next allele in next locus
				}			
			}
			#if there is more than 1 haplotype
			if (@hapseqinarr > 1) {
				#build @seqmatin:d1: rows, d2: cols
				for ($i = 0; $i < @hapseqinarr; ++$i) {
					@temparr1 = split('', $hapseqinarr[$i]);
					push @seqmatin, [@temparr1];
				}
				#call sub varpos to get position numbers of variable positions into @varposin
				@varposin = varpos(\@seqmatin);
				#call sub errbasebin to call SNPs and correct sequencing errors
				#corrected sequences into: @seqmatcorr
				@seqmatcorr = errbasebin(\@seqmatin,\@varposin,\@haprdinarr,$binerr,\%majormin);
				#call sub hapgroups to build haplotype-groups into %hapgroups
				#determine their depths into %hapgroup_dep
				#collect numbers of discarded haplos in @dischaplos
				#record how many alleles got lost in that procedure($loclostall)
				#and if the number of lost alleles could be greater($locpotlostmoreall == 1)
				#haplotype groups will later be collapsed into alleles
				hapgroupsbin(\@seqmatcorr,\%hapgroups,\%hapgroup_dep,
				\@dischaplos,$sl,\@haprdinarr,\$loclostall,\$locpotlostmoreall);
				#if there is no hapgroup:
				if (keys %hapgroups == 0) {
					$loc_cat = 'lost';#record that the locus is lost
				} else {
					#call sub build_alleles_bin to populate %alleles and determine $n_alleles
					build_alleles_bin(\@seqmatcorr,\@varposin,\%hapgroups,\%hapgroup_dep,
					\%alleles,\$n_alleles,\@dischaplos,\$loclostall,\$next_allele_ID);
					#if there is no allele after building alleles
					if (keys %alleles == 0) {
						$loc_cat = 'lost';#record that the locus is lost
					}
				}
			}
			#print output for this locus to loci-outfile
			print OUTLOCI "$locID\t$loc_cat\t$sl\t$locdep\t$n_alleles\t$loclostall\t$locpotlostmoreall\n";
			#print output for this locus to alleles-outfile
			if ($loc_cat eq 'valid') {#for a valid locus
				#loop through allele IDs of this locus
				for $allele_ID (sort {$a <=> $b} keys %alleles) {
					@alleledata = @{$alleles{$allele_ID}};#get data about this allele
					print OUTALL "$locID\t$loc_cat\t$sl\t$allele_ID\t$alleledata[0]\t$alleledata[1]\n";
				}
			} else {#for a lost locus
				print OUTALL "$locID\t$loc_cat\t$sl\t0\t0\t0\n";
			}
			#print output for this locus to svars-outfile
			if ($loc_cat eq 'valid') {#if locus is valid, i.e. contains at least one valid allele:
				#loop through allele IDs of this locus
				for $allele_ID (sort {$a <=> $b} keys %alleles) {
					@alleledata = @{$alleles{$allele_ID}};#get data about this allele
					#loop through hapNos of this allele
					for ($i = 2; $i < @alleledata; ++$i) {
						$hapNo = $alleledata[$i];
						print $out_haplotypes_h "$locID\t$loc_cat\t$sl\t$allele_ID\t",
						"$hapIDinarr[$hapNo]\tused\t$hapseqinarr[$hapNo]\t$haprdinarr[$hapNo]\n";
					}
				}
			}
			#loop through discarded haplotypes from this locus
			for ($i = 0; $i < @dischaplos; ++$i) {
				$hapNo = $dischaplos[$i];
				print $out_haplotypes_h "$locID\t$loc_cat\t$sl\t0\t",
				"$hapIDinarr[$hapNo]\tdiscarded\t$hapseqinarr[$hapNo]\t$haprdinarr[$hapNo]\n";
			}
			#collect data about this locus for reads-outfile
			#into %hdat {hapseq}=string: LocID\tLoc_cat\tseql\tAll_ID\thapID\thapcat
			if ($loc_cat eq 'valid') {#if locus is valid, i.e. contains at least one valid allele:
				#loop through allele IDs of this locus
				for $allele_ID (sort {$a <=> $b} keys %alleles) {
					@alleledata = @{$alleles{$allele_ID}};#get data about this allele
					#loop through hapNos of this allele
					for ($i = 2; $i < @alleledata; ++$i) {
						$hapNo = $alleledata[$i];
						$hapseq = $hapseqinarr[$hapNo];#get seq of haplotype
						#store data about this haplotype in %hdat
						$hdat{$hapseq} = "$locID\t$loc_cat\t$sl\t$allele_ID\t$hapIDinarr[$hapNo]\tused";
					}
				}
			}
			#loop through discarded haplotypes from this locus
			for ($i = 0; $i < @dischaplos; ++$i) {
				$hapNo = $dischaplos[$i];
				$hapseq = $hapseqinarr[$hapNo];#get seq of haplotype
				#store data about this haplotype in %hdat
				$hdat{$hapseq} = "$locID\t$loc_cat\t$sl\t$allele_ID\t$hapIDinarr[$hapNo]\tdiscarded";
			}
			#set variables back
			@catinarr = ();
			@hapseqinarr = ();
			@haprdinarr = ();
			@seqmatin = ();
			%hapgroups = ();
			@dischaplos = ();
			$locdep = 0;
			%hapgroup_dep = ();
			$loclostall = 0;
			$locpotlostmoreall = 0;
			%alleles = ();
			$n_alleles = 0;
			$loc_cat = 'valid';
			@hapIDinarr = ();
			#go to next locus
		}
	}
	#read in fastq infile, look up sequence, in %hdat, and print data to outfile OUTREADS
	if ($build_rc == 1) {#fastq infile contains only f-reads of rc-duplicates, build r-reads as reverse complement
		while ($seqline1 = <$fqfileh>) {
			$seqline2 = <$fqfileh>;
			$seqline3 = <$fqfileh>;
			$seqline4 = <$fqfileh>;
			chomp $seqline1;
			chomp $seqline2;
			$rseqline1 = $seqline1;
			$rseqline1 =~ s/_1$/_2/;
			$rseqline2 = reverse $seqline2;
			$rseqline2 =~ tr/ACGTacgt/TGCAtgca/;
			#look up f-sequence in %hdat
			if (defined $hdat{$seqline2}) {
				print $out_reads_h "$hdat{$seqline2}\t$seqline1\n";#print data to OUTREADS
			}
			#look up r-sequence in %hdat
			if (defined $hdat{$rseqline2}) {
				print $out_reads_h "$hdat{$rseqline2}\t$rseqline1\n";#print data to OUTREADS
			}
		}
	} else {#Do not build r-reads
		while ($seqline1 = <$fqfileh>) {
			$seqline2 = <$fqfileh>;
			$seqline3 = <$fqfileh>;
			$seqline4 = <$fqfileh>;
			chomp $seqline1;
			chomp $seqline2;
			#look up sequence in %hdat
			if (defined $hdat{$seqline2}) {
				print $out_reads_h "$hdat{$seqline2}\t$seqline1\n";#print data to OUTREADS
			}
		}
	}
	close OUTLOCI;
	close OUTALL;
	close $out_haplotypes_h;
	undef $out_haplotypes_h;
	close $out_reads_h;
	undef $out_reads_h;
	close $fqfileh;
	undef $fqfileh;
	$success = "analysis completed\n";#success message
	return $success;#return success message
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

#definition of subroutine errbase1 ver02
#Call SNPs with frequency threshold method
#takes an alignment of nonidentical sequences (ACGTN, no gaps) with known read depths.
#takes a minimum character frequency provided by user $min_cf (between 0 and 1)
#calculates corresponding minimum number of reads $minrd
#rounds down minimum number of reads
#characters that are rarer then $minrd are treated as error.
#N is naturally treated as error
#if in a col, there are errorbases AND only one good base,
#errorbases are changed to the good base
#if in a col, there are errorbases AND no or more than 1 good base
#errorbases are changed to N
#expects:
#ref to @inmat: input-ali: d1: rows, d2: cols
#ref to @varpos: variable positions
#ref to @rd: frequencies of the sequences (number of reads)
#$locdep: total number of reads in the locus
#$min_cf: minimum character frequency in an individual locus alignment
#
#min_cf must be between 0 and 1 or
#special value -1 which inactivates correction
#
#returns copy of @inmat with errors corrected or changed to N.
sub errbase1 {
	#declare and initialize: _r means a reference
	my ($inmat_r, $varpos_r, $rd_r, $locdep,$min_cf) = @_;
	my @inmat = @{$inmat_r};#inputali: d1: rows, d2: cols
	my @outmat = ();#outputali: d1: rows, d2: cols
	my @varpos = @{$varpos_r};#variable positions in ali
	my @rd = @{$rd_r};#read-depth of each sequence
	my $minrd = 0;#minimum number of reads for correct base
	my %charcount = ();#key: character, value total reads for this character
	my %replace = ();#key: original base, value: new base
	my @goodchar = ();#all good characters in a col
	my @errchar = ();#all error-characters in a col
	my $char = '';#a character
	my $char_rd = 0;#total read-depth for a character
	my $row = 0;
	my $col = 0;
	my @temparr = ();
	my $i = 0;
	
	#copy inputali into outputali
	#its necessary to loop through because its 2d
	for ($i = 0; $i < @inmat; ++$i) {
		@temparr = @{$inmat[$i]};
		push @outmat, [@temparr];
	}
	unless ($min_cf == -1) {#unless user inactivated character correction
		$minrd = int($locdep * $min_cf);#determine minimum number of reads for correct base
		#loop through variable positions
		for ($i = 0; $i < @varpos; ++$i) {
			$col = $varpos[$i];
			#loop through rows and count characters
			#by adding their respective read-depths
			for ($row = 0; $row < @outmat; ++$row) {
				$charcount{${$outmat[$row]}[$col]} += $rd[$row];
			}
			%replace = %charcount;#copy %charcount into %replace
			if (exists $charcount{"N"}) {#if N occurs
				delete $charcount{"N"};#delete from %charcount
				$replace{"N"} = 'N';#set replace N with N, that could change later
				push @errchar, 'N';#add N to @errchar
			}
			if (keys %charcount > 0) {#if there are other bases than N: continue
				#loop through %charcount and check frequencies
				while (($char, $char_rd) = each %charcount) {
					if ($char_rd >= $minrd) {#if frequency good: good char
						push @goodchar, $char;
					} else {#if not: error char
						push @errchar, $char;
					}
				}
				if (@errchar > 0) {#if there are error characters
					if (@goodchar == 1) {#and if there is only 1 good character
						#loop through existing characters except N (keys in %charcount)
						for $char (keys %charcount) {
							#set replace for this character to the single good one
							$replace{$char} = $goodchar[0];
						}
						#set replace for N (wether it occurs or not) to the single good char
						$replace{"N"} = $goodchar[0];					
					} else {#else: no good char or more than one
						foreach (@errchar) {#loop through error characters
							$replace{$_} = 'N';#set replace for error characters to N
						}
						if (@goodchar > 1) {#if there are good chars (now there are 0 or >1)
							foreach (@goodchar) {#loop through good characters
								$replace{$_} = $_;#set replace for good characters to the same good characters
							}
						}
					}
					#loop through rows in this col again and replace characters
					for ($row = 0; $row < @outmat; ++$row) {
						${$outmat[$row]}[$col] = $replace{${$outmat[$row]}[$col]};
					}
				}
			}
			#set variables back
			%charcount = ();
			%replace = ();
			@goodchar = ();
			@errchar = ();
			#go to next variable position
		}
	}
	return @outmat;	
}

#definition of sub errbasebin
#Calls SNPs with binomial likelihood ratio method (Glaubitz et al. 2014).
#takes an alignment of nonidentical sequences (ACGTN, no gaps) with known read depths.
#At each variable position (more than one of ACGTN):
#Counts read depths of ACGT.
#Determines minimum acceptable count of second most frequent nucleotide
#given count of most frequent nucleotide and a sequencing error rate (provided by user)
#If sequencing error rate is 0: accepts any nucleotide (ACGT).
#
#If it accepts one nucleotide: corrects all others including N to that one.
#
#If it accepts the first and second most frequent as valid:
#accepts all nucleotides with count of second most frequent one (in a case like A8 C4 G4 T1 accepts ACG).
#changes all other nucleotides to N.
#expects:
#ref to @inmat: input-ali: d1: rows, d2: cols
#ref to @varpos: variable positions
#ref to @rd: frequencies of the sequences (number of reads)
#$binerr: sequencing error rate
#ref to %majormin: {major} = minor
		#major: count of most frequent nucleotide at a variable position
		#minor: minimum count of second most frequent nucleotide to accept a SNP
#returns copy of @inmat with errors corrected or changed to N.
sub errbasebin {
	#declare and initialize: _r means a reference
	my ($inmat_r, $varpos_r, $rd_r, $binerr, $majormin_r) = @_;
	my @inmat = @{$inmat_r};#inputali: d1: rows, d2: cols
	my @outmat = ();#outputali: d1: rows, d2: cols
	my @varpos = @{$varpos_r};#variable positions in ali
	my @rd = @{$rd_r};#read-depth of each sequence
	my $minrd = 0;#minimum number of reads for correct base
	my %charcount = ();#key: character, value total reads for this character
	my @sortcounts;#read counts of ACGT at one pos sorted in descending order
	my %replace = ();#key: original base, value: new base
	my @goodchar = ();#all good characters in a col
	my @errchar = ();#all error-characters in a col
	my $char = '';#a character
	my $char_rd = 0;#total read-depth for a character
	my $row = 0;
	my $col = 0;
	my @temparr = ();
	my $i = 0;

	#copy inputali into outputali
	#its necessary to loop through because its 2d
	for ($i = 0; $i < @inmat; ++$i) {
		@temparr = @{$inmat[$i]};
		push @outmat, [@temparr];
	}
	#loop through variable positions
	for ($i = 0; $i < @varpos; ++$i) {
		$col = $varpos[$i];
		#loop through rows and count characters
		#by adding their respective read-depths
		for ($row = 0; $row < @outmat; ++$row) {
			$charcount{${$outmat[$row]}[$col]} += $rd[$row];
		}
		%replace = %charcount;#copy %charcount into %replace
		if (exists $charcount{"N"}) {#if N occurs
			delete $charcount{"N"};#delete from %charcount
			$replace{"N"} = 'N';#set replace N with N, that could change later
			push @errchar, 'N';#add N to @errchar
		}
		if (keys %charcount > 0) {#if there are other bases than N: continue
			if ((keys %charcount == 1)) {#if there is only one of ACGT
				for $char (keys %charcount) {
					push @goodchar, $char;#put it into @goodchar
				}
			}
			elsif ($binerr == 0) {#if sequencing error rate is 0
				for $char (keys %charcount) {
					push @goodchar, $char;#put all of ACGT into @goodchar
				}
			}
			else {#if there are 2 or more of ACGT and sequencing error rate > 0
				#determine acceptable ones
				#determine highest and second highest nucleotide count
				@sortcounts = reverse sort {$a <=> $b} values %charcount;
				#try to look up minrd given most frequent nucleotide
				if (defined $$majormin_r{$sortcounts[0]}){
					$minrd = $$majormin_r{$sortcounts[0]};
				} else {#if you could not find it, calculate and store in %majormin owned by sub snphapcalbin
					$minrd = $sortcounts[0] * log((1 - $binerr) / 0.5) / log(0.5 / $binerr);
					$minrd = int ($minrd + 1);#round up
					$$majormin_r{$sortcounts[0]} = $minrd;#store in %majormin
				}
				#if the second most frequent nucleotide is more frequent than minrd
				if ($sortcounts[1] > $minrd) {
					$minrd = $sortcounts[1];#elevate minrd to its frequency
					#A third or fourth nucleotide can thus only be accepted,
					#when more than one nucleotide has the second highest frequency.
				}
				#loop through %charcount and check frequencies
				while (($char, $char_rd) = each %charcount) {
					if ($char_rd >= $minrd) {#if frequency good: good char
						push @goodchar, $char;
					} else {#if not: error char
						push @errchar, $char;
					}
				}
			}
			if (@errchar > 0) {#if there are error characters
				if (@goodchar == 1) {#and if there is only 1 good character
					#loop through existing characters except N (keys in %charcount)
					for $char (keys %charcount) {
						#set replace for this character to the single good one
						$replace{$char} = $goodchar[0];
					}
					#set replace for N (wether it occurs or not) to the single good char
					$replace{"N"} = $goodchar[0];					
				} else {#else: more than one good char
					foreach (@errchar) {#loop through error characters
						$replace{$_} = 'N';#set replace for error characters to N
					}
					foreach (@goodchar) {#loop through good characters
						$replace{$_} = $_;#set replace for good characters to the same good characters
					}
				}
				#loop through rows in this col again and replace characters
				for ($row = 0; $row < @outmat; ++$row) {
					${$outmat[$row]}[$col] = $replace{${$outmat[$row]}[$col]};
				}
			}
		}
		#set variables back
		%charcount = ();
		%replace = ();
		@goodchar = ();
		@errchar = ();
		#go to next variable position
	}
	return @outmat;
}

#definition of sub hapgroups
#Used for frequency threshold method
#builds haplotype groups that will later be collapsed into alleles
#determines their depths
#in certain cases, it discards haplotypes with N (from promiscuous groups)
#if the group is not empty thereafter and does not need to be split up
#it keeps the depth of the discarded haplotypes in the group
#if the group needs to be split up into several haplotype groups thereafter
#it adds the depths of discarded haplotypes to those daughter groups
#that have a distance of 0 to the discarded haplotype
#if it adds to several groups: it splits up the depth of the discarded haplotype
#and adds parts of it to daughter-groups
#parts are proportional to daughter-groups depths
#
#expects:
#ref to @seqmat: (corrected) sequences in inputorder
					#d1: seq(row), d2: char(col)
#ref to %hapgroups: {groupNo}=@ of hapNos, still empty, will be populated
					#hapNos correspond to d1-index of @seqmatcorr
#ref to %hapgroup_dep: {groupNo} = depth
#ref to @dischaplos: #hapNos of discarded haplos, still empty, will be populated
					#hapNos correspond to d1-index of @seqmatcorr
#$sl: sequence length
#ref to @haprdinarr: #depth of each haplo, index: hapNo
#ref to $loclostall: #number of alleles lost from this locus through discard of haplos
					#here in fact a minimum number, could always be greater
#ref to $locpotlostmoreall: #1 if the number of lost alleles could be greater than $loclostall, else 0
sub hapgroups {
	#declare and initialize: _r means a reference
	my ($seqmat_r, $hapgroups_r, $hapgroup_dep_r,
	$dischaplos_r, $sl, $haprdinarr_r,
	$loclostall_r, $locpotlostmoreall_r) = @_;
	my @seqmat = @{$seqmat_r};#copy of the inseq-matrix, don't modify - is 2d!
	my @haprdinarr = @{$haprdinarr_r};#depth of each haplo, index: hapNo
	my @inseqs = ();#all sequences in @seqmat as strings, inputorder
	my $nhap = 0;#number of haplos
	my @pwd01 = ();#full pw distance matrix: 0 equal, 1 different
					#under pw ignorance of N
	my $diststring = '';#one element of @pwd01
	my %countdiststrings = ();#for counting kinds of distance strings ($pwd01[])
	my @haphasN = ();#inputorder: 0/1 for haplo contains N
	my %dist0pairs = ();# {hapNo1}{hapNo2} = 0;
	my $hap1No = 0;#number of haplotype1 in a pair
	my $hap2No = 0;#number of haplotype2 in a pair
	my $hap1group = 0;#group-number of haplotype1 in a pair
	my $hap2group = 0;#group-number of haplotype2 in a pair
	my $dist = 0;#pairwise distance
	my @hapNo_group = ();#groupNo of each haplo, index: input-order-hapNo
	my %group_hapNos = ();# {groupNo} = @of hapNos
	my %diststring_hapNos = ();# {$diststring} = @of hapNos
				#temporarily holds groups (arrays of hapNos)
				#with members that have been taken out of a promiscuous group
	my $groupNo = 1;#a group-number
	my $currgroupNo = 0;#number of group currently evaluated
	my %group_dep = ();# {groupNo} = depth
	my %diststring_dep = ();# {diststring} = depth of temporary group
							#in %diststring_hapNos
	my %dischaplos_currgroup = ();#haplos that have been discarded from current group
	my @donors = ();#hapNos of discarded haplotypes that will donate depth to a remaining haplotype
					#donors have distance 0 to a receiver (remaining haplotype)
	my $donor = 0;#hapNo of one donor
	my @notdonors = ();#hapNos of discarded haplotypes that are not donors
						#they have distance >0 to a receiver (remaining haplotype)
	my $notdonor = 0;#hapNo of one notdonor
	my $dischapNo = 0;#number of a discarded haplo
	my $disc_dep = 0;#depth of a discarded haplo
	my $receiver = 0;#hapNo of a haplotype that will receive additional depth
	my $totaldep_receivers = 0;#total depth of groups that will receive
								#additional depth from a discarded haplo
	my %receivers_dep = ();# {diststring}= dep depths of temporary groups
						#that will receive additional depth from a discarded haplo
	my $receiver_dep = 0;#depth of one temporary group that will receive additional depth
						#from a discarded haplo
	my $add_dep = 0;#depth that will be added to a receiver
	my @temparr1 = ();
	my $tempstring = '';
	my $i = 0;
	
	#build @inseqs
	for ($i = 0; $i < @seqmat; ++$i) {
		$tempstring = join('',@{$seqmat[$i]});
		push @inseqs, $tempstring;
	}	
	$nhap = @inseqs;#determine number of haplos	
	#initialize @pwd01 with 0 on diagonal and 1 in all other cells
	$tempstring = 1 x $nhap;
	for ($i = 0; $i < $nhap; ++$i) {
		push @pwd01, $tempstring;
		substr($pwd01[$i],$i,1) = 0;
	}
	#determine which haplos contain N
	for ($i = 0; $i < @inseqs; ++$i) {
		if ($inseqs[$i] =~ /N/) {
			$haphasN[$i] = 1;
		} else {
			$haphasN[$i] = 0;
		}
	}
	#determine all pairwise distances between haplos
	#loop through all pairs of seqs in inseq
	for ($hap1No = 0; $hap1No < ($nhap - 1); ++$hap1No) {
		for ($hap2No = ($hap1No + 1); $hap2No < $nhap; ++$hap2No) {
			#get the distance
			$dist = distpwdN($inseqs[$hap1No],$inseqs[$hap2No],$haphasN[$hap1No],$haphasN[$hap2No],$sl);
			if ($dist == 0) {#if the distance is 0
				#put the pair into %dist0pairs
				$dist0pairs{$hap1No}{$hap2No} = 0;
				#enter a distance of 0 in both triangles of @pwd01:
				substr($pwd01[$hap1No],$hap2No,1) = 0;
				substr($pwd01[$hap2No],$hap1No,1) = 0;
			}
		}
	}
	#initialize @hapNo_group with group 0 for all hapNos
	for ($i = 0; $i < $nhap; ++$i) {
		$hapNo_group[$i] = 0;
	}
	#build groups of haplotypes:
	#begin with haplos that are part of a pair
	#loop through the pairs (key1-key2) in %dist0pairs
	for $hap1No (keys %dist0pairs) {
		for $hap2No (keys %{$dist0pairs{$hap1No}}) {
			#check the current group-numbers of these haplos
			$hap1group = $hapNo_group[$hap1No];
			$hap2group = $hapNo_group[$hap2No];
			#if both are not yet in a group (group-No 0)
			if ($hap1group == 0 and $hap2group == 0) {
				#make up a new group for them
				@{$group_hapNos{$groupNo}} = ($hap1No,$hap2No);
				$hapNo_group[$hap1No] = $groupNo;
				$hapNo_group[$hap2No] = $groupNo;
				$group_dep{$groupNo} = $haprdinarr[$hap1No];#initialize group_dep with depth of haplo1
				$group_dep{$groupNo} += $haprdinarr[$hap2No];#add depth of haplo2 to group_dep
				++$groupNo;
			}
			#if haplo1 already is in a group but haplo2 is not
			elsif ($hap1group != 0 and $hap2group == 0) {
				#put haplo 2 into group of haplo1
				push @{$group_hapNos{$hap1group}}, $hap2No;
				$hapNo_group[$hap2No] = $hap1group;
				$group_dep{$hap1group} += $haprdinarr[$hap2No];#add depth of haplo 2 to group_dep
			}
			#if haplo1 is not yet in a group but haplo2 is
			elsif ($hap1group == 0 and $hap2group != 0) {
				#put haplo 1 into group of haplo2
				push @{$group_hapNos{$hap2group}}, $hap1No;
				$hapNo_group[$hap1No] = $hap2group;
				$group_dep{$hap2group} += $haprdinarr[$hap1No];#add depth of haplo 1 to group_dep
			}
			#if both are in different groups
			elsif ($hap1group != $hap2group) {
				#join their groups:
				#get all members of the group hap2 is in
				@temparr1 = @{$group_hapNos{$hap2group}};
				#loop through them and change their group-numbers to the one of hap1
				for ($i = 0; $i < @temparr1; ++$i) {
					$hapNo_group[$temparr1[$i]] = $hap1group;
				}
				#add them all to the group of hap1
				push @{$group_hapNos{$hap1group}}, @temparr1;
				#add the depth of hap2group to depth of hap1group
				$group_dep{$hap1group} += $group_dep{$hap2group};
				#delete the former group of hap2
				delete $group_hapNos{$hap2group};
				delete $group_dep{$hap2group};
			}			
		}
	}
	#search promiscuous groups i.e. groups in which not all pw-distances are 0
	for $currgroupNo (keys %group_hapNos) {#loop through groups in %group_hapNos
		@temparr1 = @{$group_hapNos{$currgroupNo}};#get all hapNos in that group
		for ($i = 0; $i < @temparr1; ++$i) {#loop through hapNos in that group
			#collect distance strings (@pwd01) of members in %countdiststrings
			++$countdiststrings{$pwd01[$temparr1[$i]]};
		}
		if (keys %countdiststrings > 1) {#if there is more than 1 variant
			#the group is promiscuous:
			#remove all groupmembers that have N
			$i = 0;
			while ($i < @temparr1) {#loop through the members (hapNos)
				$hap1No = $temparr1[$i];#get current number
				if ($haphasN[$hap1No] == 1) {#if current member has N
					splice (@temparr1,$i,1);#remove it from group
					$hapNo_group[$hap1No] = -1;#set its group number to -1
					push @{$dischaplos_r}, $hap1No;#put it into @dischaplos
					$dischaplos_currgroup{$hap1No} = 1;#add it to %dischaplos_currgroup
				} else {#if not increment index
					++$i;
				}
			}
			if (@temparr1 == 0) {#if the group is empty now
				#delete it from %group_hapNos and %group_dep
				delete $group_hapNos{$currgroupNo};
				delete $group_dep{$currgroupNo};
				#the locus lost at least 2 alleles:
				$$loclostall_r += 2;
				if (keys %dischaplos_currgroup > 2) {#if at least 3 haplotypes were discarded
					$$locpotlostmoreall_r = 1;#the locus could have lost more alleles
				}
			}
			elsif (@temparr1 == 1) {#if it contains 1 member (the receiver)
				@{$group_hapNos{$currgroupNo}} = @temparr1;
				$receiver = $temparr1[0];
				#determine the depth of the receiver
				$receiver_dep = $haprdinarr[$receiver];
				#set group depth down to depth of the receiver (will grow later)
				$group_dep{$currgroupNo} = $receiver_dep;
				#loop through discarded haplotypes
				#identify donors and notdonors
				for $dischapNo (keys %dischaplos_currgroup) {
					#get diststring of current discarded haplo
					$diststring = $pwd01[$dischapNo];
					#determine distance to receiver
					$dist  = substr($diststring,$temparr1[0],1);
					if ($dist == 0) {#if the distance is 0
						push @donors, $dischapNo;#add it to the donors (will donate a proportion of their depth)
					} else {#not a donor
						push @notdonors, $dischapNo;
					}
				}
				#the locus lost at least 1 allele
				$$loclostall_r += 1;
				if (keys %dischaplos_currgroup > 1) {#if at least 2 haplotypes were discarded
					$$locpotlostmoreall_r = 1;#the locus could have lost more alleles
				}
				for $donor (@donors) {#loop through the donors
					$disc_dep = $haprdinarr[$donor];#get donor's depth
					#initialize total depth of all possible receivers with that of the (single true) receiver
					$totaldep_receivers = $receiver_dep;
					for $notdonor (@notdonors) {#loop through notdonors
						#get diststring of current notdonor
						$diststring = $pwd01[$notdonor];
						#determine distance to donor
						$dist = substr($diststring,$donor,1);
						if ($dist == 0) {#if it is 0
							#add depth of notdonor to total depth of receivers
							$totaldep_receivers += $haprdinarr[$notdonor];
						}
					}
					#calculate the proportion of depth of current donor to add to depth of receiver/group
					$add_dep = $receiver_dep / $totaldep_receivers * $disc_dep;
					#add it to depth of group
					$group_dep{$currgroupNo} += $add_dep;
				}
				#set variables back
				$totaldep_receivers = 0;
				@donors = ();
				@notdonors = ();
			}
			else {#if it has more than 1 member (several receivers)
				#delete it from %group_hapNos and %group_dep
				delete $group_hapNos{$currgroupNo};
				delete $group_dep{$currgroupNo};
				#loop through the members
				for ($i = 0; $i < @temparr1; ++$i) {
					$hap1No = $temparr1[$i];#get current hapNo
					$hapNo_group[$hap1No] = -1;#set its group number to -1
					$diststring = $pwd01[$hap1No];#get its diststring
					#sort it into a temporary group in %diststring_hapNos
					push @{$diststring_hapNos{$diststring}}, $hap1No;
					#add depth of $hap1No to depth of its temporary group
					$diststring_dep{$diststring} += $haprdinarr[$hap1No];
				}				
				#loop through numbers of discarded haplos from this group
				#identify donors and notdonors
				for $dischapNo (keys %dischaplos_currgroup) {
					#get diststring of current discarded haplo
					$diststring = $pwd01[$dischapNo];
					$dist = 1;#initialize distance with 1
					#loop through receivers
					for $receiver (@temparr1) {
						#multiply $dist with distance between current discarded haplo
						#and current receiver
						$dist *= substr($diststring,$receiver,1);
					}
					#if $dist is 0 now, distance between current discarded haplo
					#and at least one receiver is 0
					#-> current discarded haplo is a donor
					if ($dist == 0) {
						push @donors, $dischapNo;
					} else {#if not, its a notdonor
						push @notdonors, $dischapNo;
					}
				}
				if (@notdonors > 0) {#if there is at least one notdonor
					$$loclostall_r += 1;#the locus lost at least one allele
					if (keys %dischaplos_currgroup > 1) {#if at least 2 haplotypes were discarded
						$$locpotlostmoreall_r = 1;#the locus potentially lost more alleles
					}
				}
				elsif (@donors > 0) {#no notdonors but at least one donor
					$$locpotlostmoreall_r = 1;#at least we know that the locus could have lost one allele or more
				}
				#Donate proportions of donors' depths to receivers'depths
				for $donor (@donors) {#loop through donors
					#determine depth of current donor
					$disc_dep = $haprdinarr[$donor];
					#determine total depth of all receivers for this donor
						$diststring = $pwd01[$donor];#get diststring of current donor
						for $receiver (@temparr1) {#loop through receivers
							#determine distance to donor
							$dist = substr($diststring,$receiver,1);
							if ($dist == 0) {#if it is 0
								#add depth of receiver to total depth of receivers
								$totaldep_receivers += $haprdinarr[$receiver];
							}
						}
						for $notdonor (@notdonors) {#loop through notdonors
							#determine distance to donor
							$dist = substr($diststring,$notdonor,1);
							if ($dist == 0) {#if it is 0
								#add depth of notdonor to total depth of receivers
								$totaldep_receivers += $haprdinarr[$notdonor];
							}
						}
					#loop through receivers again and donate
					for $receiver (@temparr1) {
						#determine distance to donor
						$dist = substr($diststring,$receiver,1);
						if ($dist == 0) {#if it is 0
							#get diststring of receiver (needed as key to its current temporary group)
							$tempstring = $pwd01[$receiver];
							#get original depth of receiver
							$receiver_dep = $haprdinarr[$receiver];
							#calculate proportion of donor's depth to donate to receiver's depth
							$add_dep = $receiver_dep / $totaldep_receivers * $disc_dep;
							#donate depth
							$diststring_dep{$tempstring} += $add_dep;							
						}
					}
					#set variables back and next donor
					$totaldep_receivers = 0;
				}
				#set variables back
				@donors = ();
				@notdonors = ();
			}#end of "if it has more than 1 member"
			#set variables back
			%dischaplos_currgroup = ();
		}#end of block "the group is promiscuous"
		%countdiststrings = ();#set back and continue searching promiscuous groups with next group
	}#end of search promiscuous groups
	#make temporary groups in %diststring_hapNos to new regular groups
	#loop through keys of %diststring_hapNos
	for $diststring (keys %diststring_hapNos) {
		#get members of current group
		@temparr1 = @{$diststring_hapNos{$diststring}};
		#loop through members of current group
		for ($i = 0; $i < @temparr1; ++$i) {
			$hapNo_group[$temparr1[$i]] = $groupNo;#give them a new groupNo
		}
		@{$group_hapNos{$groupNo}} = @temparr1;#create that new group
		#record its depth
		$group_dep{$groupNo} = $diststring_dep{$diststring};
		++$groupNo;#increment $groupNo for next new group
	}
	#make haplos that have not yet a group (and have not been discarded)
	#to new single-member groups
	#loop through @hapNo_group
	for ($i = 0; $i < @hapNo_group; ++$i) {
		if ($hapNo_group[$i] == 0) {#if groupNo is still 0
			$hapNo_group[$i] = $groupNo;#give a new groupNo
			push @{$group_hapNos{$groupNo}}, $i;#create the new group
			#record its depth
			$group_dep{$groupNo} = $haprdinarr[$i];
			++$groupNo;#increment $groupNo for next new group
		}
	}
	#populate %hapgroups (by ref, owned by sub snphapcal1)
	#and %hapgroup_dep (by ref, owned by sub snphapcal1)
	#use new consecutive group-numbers as keys
	$groupNo = 1;
	#loop through groups in %group_hapNos
	for $currgroupNo (keys %group_hapNos) {
		#get the haplo-numbers of the members
		@temparr1 = @{$group_hapNos{$currgroupNo}};
		#make them a group in %hapgroups
		@{$$hapgroups_r{$groupNo}} = @temparr1;
		#add its depth to %hapgroup_dep
		$$hapgroup_dep_r{$groupNo} = $group_dep{$currgroupNo};
		++$groupNo;#increment for next consecutive group-number
	}
}

#definition of sub hapgroupsbin
#Used for binomial likelihood ratio method
#builds haplotype groups that will later be collapsed into alleles
#promiscuous groups: not all pairwise distances are zero
#arise through connecting haplotypes containing N
#sub hapgroups_bin treats such groups according to these rules
#1. Remove all haplotypes containing N
#2. Retained haplotypes become new groups
#3. Discarded haplotypes with distance 0 to a retained one are donors
#4. Donors only donate and never receive depth
#5. Receivers are all svars that have distance 0 to a donor
#6. Receivers only receive and never donate depth
#7. Donors that have one single receiver, which is a retained haplotype
#   donate their complete depth to this receiver
#8. Depth of other donors is discarded
#
#discards all haplotypes containing N
#retained haplotypes become new groups
#those discarded haplotypes that have distance 0 to only one retained haplotype
#donate their depth to this haplotype
#depth of other discarded haplotypes is discarded
#
#expects:
#ref to @seqmat: (corrected) sequences in inputorder
					#d1: seq(row), d2: char(col)
#ref to %hapgroups: {groupNo}=@ of hapNos, still empty, will be populated
					#hapNos correspond to d1-index of @seqmatcorr
#ref to %hapgroup_dep: {groupNo} = depth
#ref to @dischaplos: #hapNos of discarded haplos, still empty, will be populated
					#hapNos correspond to d1-index of @seqmatcorr
#$sl: sequence length (because already known)
#ref to @haprdinarr: #depth of each haplo, index: hapNo
#ref to $loclostall: #number of alleles lost from this locus through discard of haplos
					#here in fact a minimum number, could always be greater
#ref to $locpotlostmoreall: #1 if the number of lost alleles could be greater than $loclostall, else 0
#uses subroutines
#distpwdN
sub hapgroupsbin {
	#declare and initialize: _r means a reference
	my ($seqmat_r, $hapgroups_r, $hapgroup_dep_r,
	$dischaplos_r, $sl, $haprdinarr_r,
	$loclostall_r, $locpotlostmoreall_r) = @_;
	my @seqmat = @{$seqmat_r};#copy of the inseq-matrix, don't modify - is 2d!
	my @haprdinarr = @{$haprdinarr_r};#depth of each haplo, index: hapNo
	my @inseqs = ();#all sequences in @seqmat as strings, inputorder
	my $nhap = 0;#number of haplos
	my @pwd01 = ();#full pw distance matrix: 0 equal, 1 different
					#under pw ignorance of N
	my $diststring = '';#one element of @pwd01
	my %countdiststrings = ();#for counting kinds of distance strings ($pwd01[])
	my @haphasN = ();#inputorder: 0/1 for haplo contains N
	my %dist0pairs = ();# {hapNo1}{hapNo2} = 0;
	my $hap1No = 0;#number of haplotype1 in a pair
	my $hap2No = 0;#number of haplotype2 in a pair
	my $hap1group = 0;#group-number of haplotype1 in a pair
	my $hap2group = 0;#group-number of haplotype2 in a pair
	my $dist = 0;#pairwise distance
	my @hapNo_group = ();#groupNo of each haplo, index: input-order-hapNo
	my %group_hapNos = ();# {groupNo} = @of hapNos
	my %diststring_hapNos = ();# {$diststring} = @of hapNos
				#temporarily holds groups (arrays of hapNos)
				#with members that have been taken out of a promiscuous group
	my $groupNo = 1;#a group-number
	my $currgroupNo = 0;#number of group currently evaluated
	my %group_dep = ();# {groupNo} = depth
	my %diststring_dep = ();# {diststring} = depth of temporary group
							#in %diststring_hapNos
	my %dischaplos_currgroup = ();#haplos that have been discarded from current group
	my @donors = ();#hapNos of discarded haplotypes that could donate depth to a retained haplotype
					#donors have distance 0 to a receiver (retained haplotype)
	my $donor = 0;#hapNo of one donor
	my @notdonors = ();#hapNos of discarded haplotypes that are not donors
						#they have distance >0 to all retained haplotypes
	my $notdonor = 0;#hapNo of one notdonor
	my $dischapNo = 0;#number of a discarded haplo
	my $receiver = 0;#hapNo of a haplotype that will receive additional depth
	my $n_receivers = 0;#number of potential receivers for a donor
	my $add_dep = 0;#depth that will be added to a receiver
	my @temparr1 = ();
	my $tempstring = '';
	my $i = 0;
	
	#build @inseqs
	for ($i = 0; $i < @seqmat; ++$i) {
		$tempstring = join('',@{$seqmat[$i]});
		push @inseqs, $tempstring;
	}	
	$nhap = @inseqs;#determine number of haplos	
	#initialize @pwd01 with 0 on diagonal and 1 in all other cells
	$tempstring = 1 x $nhap;
	for ($i = 0; $i < $nhap; ++$i) {
		push @pwd01, $tempstring;
		substr($pwd01[$i],$i,1) = 0;
	}
	#determine which haplos contain N
	for ($i = 0; $i < @inseqs; ++$i) {
		if ($inseqs[$i] =~ /N/) {
			$haphasN[$i] = 1;
		} else {
			$haphasN[$i] = 0;
		}
	}
	#determine all pairwise distances between haplos
	#loop through all pairs of seqs in inseq
	for ($hap1No = 0; $hap1No < ($nhap - 1); ++$hap1No) {
		for ($hap2No = ($hap1No + 1); $hap2No < $nhap; ++$hap2No) {
			#get the distance
			$dist = distpwdN($inseqs[$hap1No],$inseqs[$hap2No],$haphasN[$hap1No],$haphasN[$hap2No],$sl);
			if ($dist == 0) {#if the distance is 0
				#put the pair into %dist0pairs
				$dist0pairs{$hap1No}{$hap2No} = 0;
				#enter a distance of 0 in both triangles of @pwd01:
				substr($pwd01[$hap1No],$hap2No,1) = 0;
				substr($pwd01[$hap2No],$hap1No,1) = 0;
			}
		}
	}
	#initialize @hapNo_group with group 0 for all hapNos
	for ($i = 0; $i < $nhap; ++$i) {
		$hapNo_group[$i] = 0;
	}
	#build groups of haplotypes:
	#begin with haplos that are part of a pair
	#loop through the pairs (key1-key2) in %dist0pairs
	for $hap1No (keys %dist0pairs) {
		for $hap2No (keys %{$dist0pairs{$hap1No}}) {
			#check the current group-numbers of these haplos
			$hap1group = $hapNo_group[$hap1No];
			$hap2group = $hapNo_group[$hap2No];
			#if both are not yet in a group (group-No 0)
			if ($hap1group == 0 and $hap2group == 0) {
				#make up a new group for them
				@{$group_hapNos{$groupNo}} = ($hap1No,$hap2No);
				$hapNo_group[$hap1No] = $groupNo;
				$hapNo_group[$hap2No] = $groupNo;
				$group_dep{$groupNo} = $haprdinarr[$hap1No];#initialize group_dep with depth of haplo1
				$group_dep{$groupNo} += $haprdinarr[$hap2No];#add depth of haplo2 to group_dep
				++$groupNo;
			}
			#if haplo1 already is in a group but haplo2 is not
			elsif ($hap1group != 0 and $hap2group == 0) {
				#put haplo 2 into group of haplo1
				push @{$group_hapNos{$hap1group}}, $hap2No;
				$hapNo_group[$hap2No] = $hap1group;
				$group_dep{$hap1group} += $haprdinarr[$hap2No];#add depth of haplo 2 to group_dep
			}
			#if haplo1 is not yet in a group but haplo2 is
			elsif ($hap1group == 0 and $hap2group != 0) {
				#put haplo 1 into group of haplo2
				push @{$group_hapNos{$hap2group}}, $hap1No;
				$hapNo_group[$hap1No] = $hap2group;
				$group_dep{$hap2group} += $haprdinarr[$hap1No];#add depth of haplo 1 to group_dep
			}
			#if both are in different groups
			elsif ($hap1group != $hap2group) {
				#join their groups:
				#get all members of the group hap2 is in
				@temparr1 = @{$group_hapNos{$hap2group}};
				#loop through them and change their group-numbers to the one of hap1
				for ($i = 0; $i < @temparr1; ++$i) {
					$hapNo_group[$temparr1[$i]] = $hap1group;
				}
				#add them all to the group of hap1
				push @{$group_hapNos{$hap1group}}, @temparr1;
				#add the depth of hap2group to depth of hap1group
				$group_dep{$hap1group} += $group_dep{$hap2group};
				#delete the former group of hap2
				delete $group_hapNos{$hap2group};
				delete $group_dep{$hap2group};
			}			
		}
	}
	#search promiscuous groups i.e. groups in which not all pw-distances are 0
	for $currgroupNo (keys %group_hapNos) {#loop through groups in %group_hapNos
		@temparr1 = @{$group_hapNos{$currgroupNo}};#get all hapNos in that group
		for ($i = 0; $i < @temparr1; ++$i) {#loop through hapNos in that group
			#collect distance strings (@pwd01) of members in %countdiststrings
			++$countdiststrings{$pwd01[$temparr1[$i]]};
		}
		if (keys %countdiststrings > 1) {#if there is more than 1 variant
			#the group is promiscuous:
			#remove all groupmembers that have N
			$i = 0;
			while ($i < @temparr1) {#loop through the members (hapNos)
				$hap1No = $temparr1[$i];#get current number
				if ($haphasN[$hap1No] == 1) {#if current member has N
					splice (@temparr1,$i,1);#remove it from group
					$hapNo_group[$hap1No] = -1;#set its group number to -1
					push @{$dischaplos_r}, $hap1No;#put it into @dischaplos
					$dischaplos_currgroup{$hap1No} = 1;#add it to %dischaplos_currgroup
				} else {#if not increment index
					++$i;
				}
			}
			if (@temparr1 == 0) {#if the group is empty now
				#delete it from %group_hapNos and %group_dep
				delete $group_hapNos{$currgroupNo};
				delete $group_dep{$currgroupNo};
				#the locus lost at least 2 alleles:
				$$loclostall_r += 2;
				if (keys %dischaplos_currgroup > 2) {#if at least 3 haplotypes were discarded
					$$locpotlostmoreall_r = 1;#the locus could have lost more alleles
				}
			}
			elsif (@temparr1 == 1) {#if it contains 1 member (the receiver)
				@{$group_hapNos{$currgroupNo}} = @temparr1;#the group contains only the receiver
				#set group depth down to depth of the receiver (will grow later)
				$group_dep{$currgroupNo} = $haprdinarr[$temparr1[0]];
				#loop through discarded haplotypes
				#identify donors and notdonors
				for $dischapNo (keys %dischaplos_currgroup) {
					#get diststring of current discarded haplo
					$diststring = $pwd01[$dischapNo];
					#determine distance to receiver
					$dist  = substr($diststring,$temparr1[0],1);
					if ($dist == 0) {#if the distance is 0
						push @donors, $dischapNo;#add it to the donors (will donate a proportion of their depth)
					} else {#not a donor
						push @notdonors, $dischapNo;
					}
				}
				#the locus lost at least 1 allele
				$$loclostall_r += 1;
				if (keys %dischaplos_currgroup > 1) {#if at least 2 haplotypes were discarded
					$$locpotlostmoreall_r = 1;#the locus could have lost more alleles
				}
				for $donor (@donors) {#loop through the donors
					$add_dep = $haprdinarr[$donor];#set depth to donate to donor's depth
					for $notdonor (@notdonors) {#loop through notdonors
						#get diststring of current notdonor
						$diststring = $pwd01[$notdonor];
						#determine distance to donor
						$dist = substr($diststring,$donor,1);
						if ($dist == 0) {#if it is 0
							#this donor cannot donate, it has more than one potential receiver
							$add_dep = 0;#set depth to donate to 0
							last;#don't look at other notdonors
						}
					}
					#donate depth to receiver (nothing if $add_dep == 0)
					$group_dep{$currgroupNo} += $add_dep;
				}
				#set variables back
				@donors = ();
				@notdonors = ();
			}
			else {#if it has more than 1 member (several receivers)
				#delete it from %group_hapNos and %group_dep
				delete $group_hapNos{$currgroupNo};
				delete $group_dep{$currgroupNo};
				#loop through the members
				for ($i = 0; $i < @temparr1; ++$i) {
					$hap1No = $temparr1[$i];#get current hapNo
					$hapNo_group[$hap1No] = -1;#set its group number to -1
					$diststring = $pwd01[$hap1No];#get its diststring
					#sort it into a temporary group in %diststring_hapNos
					push @{$diststring_hapNos{$diststring}}, $hap1No;
					#add depth of $hap1No to depth of its temporary group
					$diststring_dep{$diststring} += $haprdinarr[$hap1No];
				}				
				#loop through numbers of discarded haplos from this group
				#identify donors and notdonors
				for $dischapNo (keys %dischaplos_currgroup) {
					#get diststring of current discarded haplo
					$diststring = $pwd01[$dischapNo];
					$dist = 1;#initialize distance with 1
					#loop through receivers
					for $receiver (@temparr1) {
						#multiply $dist with distance between current discarded haplo
						#and current receiver
						$dist *= substr($diststring,$receiver,1);
					}
					#if $dist is 0 now, distance between current discarded haplo
					#and at least one receiver is 0
					#-> current discarded haplo is a donor
					if ($dist == 0) {
						push @donors, $dischapNo;
					} else {#if not, its a notdonor
						push @notdonors, $dischapNo;
					}
				}
				if (@notdonors > 0) {#if there is at least one notdonor
					$$loclostall_r += 1;#the locus lost at least one allele
					if (keys %dischaplos_currgroup > 1) {#if at least 2 haplotypes were discarded
						$$locpotlostmoreall_r = 1;#the locus potentially lost more alleles
					}
				}
				elsif (@donors > 0) {#no notdonors but at least one donor
					$$locpotlostmoreall_r = 1;#at least we know that the locus could have lost one allele or more
				}
				#Donate depth of donors to receivers where appropriate
				for $donor (@donors) {#loop through donors
					$add_dep = $haprdinarr[$donor];#set depth to donate to donor's depth
					#Determine potential receivers for this donor
					$diststring = $pwd01[$donor];#get diststring of current donor
					for $hap1No (@temparr1) {#loop through receivers
						#determine distance to donor
						$dist = substr($diststring,$hap1No,1);
						if ($dist == 0) {#if it is 0
							++$n_receivers;#count a potential receiver
							$receiver = $hap1No;#keep this one as receiver (will receive if it was the only one)
						}						
					}
					for $notdonor (@notdonors) {#loop through notdonors
						#determine distance to donor
						$dist = substr($diststring,$notdonor,1);
						if ($dist == 0) {#if it is 0
							++$n_receivers;#count a potential receiver
						}					
					}
					if ($n_receivers == 1) {#if there is only one receiver
						$diststring_dep{$pwd01[$receiver]} += $add_dep;#donate depth to receiver
					}
					$n_receivers = 0;#set back and go to next donor
				}
				#set variables back
				@donors = ();
				@notdonors = ();
			}#end of "if it has more than 1 member"
			#set variables back
			%dischaplos_currgroup = ();
		}#end of block "the group is promiscuous"
		%countdiststrings = ();#set back and continue searching promiscuous groups with next group
	}#end of search promiscuous groups
	#make temporary groups in %diststring_hapNos to new regular groups
	#loop through keys of %diststring_hapNos
	for $diststring (keys %diststring_hapNos) {
		#get members of current group
		@temparr1 = @{$diststring_hapNos{$diststring}};
		#loop through members of current group
		for ($i = 0; $i < @temparr1; ++$i) {
			$hapNo_group[$temparr1[$i]] = $groupNo;#give them a new groupNo
		}
		@{$group_hapNos{$groupNo}} = @temparr1;#create that new group
		#record its depth
		$group_dep{$groupNo} = $diststring_dep{$diststring};
		++$groupNo;#increment $groupNo for next new group
	}
	#make haplos that have not yet a group (and have not been discarded)
	#to new single-member groups
	#loop through @hapNo_group
	for ($i = 0; $i < @hapNo_group; ++$i) {
		if ($hapNo_group[$i] == 0) {#if groupNo is still 0
			$hapNo_group[$i] = $groupNo;#give a new groupNo
			push @{$group_hapNos{$groupNo}}, $i;#create the new group
			#record its depth
			$group_dep{$groupNo} = $haprdinarr[$i];
			++$groupNo;#increment $groupNo for next new group
		}
	}
	#populate %hapgroups (by ref, owned by sub snphapcal1)
	#and %hapgroup_dep (by ref, owned by sub snphapcal1)
	#use new consecutive group-numbers as keys
	$groupNo = 1;
	#loop through groups in %group_hapNos
	for $currgroupNo (keys %group_hapNos) {
		#get the haplo-numbers of the members
		@temparr1 = @{$group_hapNos{$currgroupNo}};
		#make them a group in %hapgroups
		@{$$hapgroups_r{$groupNo}} = @temparr1;
		#add its depth to %hapgroup_dep
		$$hapgroup_dep_r{$groupNo} = $group_dep{$currgroupNo};
		++$groupNo;#increment for next consecutive group-number
	}
}

#definition of subroutine build_alleles
#Used for frequency threshold method
#expects:
#ref to @seqmat: d1: rows, d2 cols, alignment of haplotypes of a locus
#ref to @varpos: variable positions in the locus
#ref to %hapgroups: {groupNo}=@of hapNos: hapNo: index to @seqmat
#ref to %hapgroup_dep: {groupNo} = depth of each hapgroup
#$locdep: total depth of the locus, may be greater than sum of hapgroup-depths
		#includes also depths of discarded haplos
#to populate/determine:
#ref to %alleles: {allele_ID} = @: [0]: allele_seq [1]: depth, remainder: hapNos
#ref to $n_alleles: number of alleles
#ref to @dischaplos: numbers of discarded haplos
#ref to $loclostall: number of lost alleles
#ref to $next_allele_ID: allele_ID of the first allele to be build here
#$min_af: minimum allele frequency in an individual locus alignment (between 0 and 1)
##
#Collapses haplotypes into alleles.
#Calculates minimum number of reads for a good allele
#as $minrd = int($min_af * $locdep)
#Discards alleles with depth < $minrd
#Discards alleles containing N
#Moves corresponding haplotypes into @dischaplos.
sub build_alleles {
	#declare and initialize: _r means a reference
	my ($seqmat_r,$varpos_r,$hapgroups_r,$hapgroup_dep_r,
	$locdep,$alleles_r,$n_alleles_r,$dischaplos_r,
	$loclostall_r,$next_allele_ID_r,$min_af) = @_;
	my @seqmat = @{$seqmat_r};#copy of outer array of @seqmat
	my $allele_ID = $$next_allele_ID_r;#allele ID for the next allele to be build
	my $minrd = 0;#minimum read-depth corresponding to $min_af
	my $currgroup = 0;#number of current hapgroup
	my $currgroup_dep = 0;#depth of current hapgroup
	my @currgroup_hapNos = ();#haplotype-numbers of current group
	my $hapNo = 0;#number of a haplotype
	my @ali = ();#all seqs of one hapgroup: d1: rows, d2: cols
	my $allele = '';#an allele
	my $i = 0;
	my @temparr = ();
	
	$minrd = int($locdep * $min_af); #calculate $minrd
	#loop through haplotype-groups
	for $currgroup (keys %{$hapgroups_r}) {
		$currgroup_dep = $$hapgroup_dep_r{$currgroup};
		@currgroup_hapNos = @{$$hapgroups_r{$currgroup}};
		#if the depth of the haplotype is too small
		if ($currgroup_dep < $minrd) {
			++$$loclostall_r;#count a lost allele
			#add the haplotype numbers to discarded-haplotype-numbers
			push @{$dischaplos_r}, @currgroup_hapNos;
		} else {#if the depth is acceptable
			#loop through group-members and build ali
			for $hapNo (@currgroup_hapNos) {
				#get sequence of this hapNo (as array)
				@temparr = @{$seqmat[$hapNo]};
				push @ali, [@temparr];
			}
			#call sub consensus1 to collapse the ali to an allele
			$allele = consensus1(\@ali,$varpos_r);
			if ($allele =~ /N/) {#if the allele contains N
				++$$loclostall_r;#count a lost allele
				#add the haplotype numbers to discarded-haplotype-numbers
				push @{$dischaplos_r}, @currgroup_hapNos;
			} else {#if there is no N in the allele
				#store the allele, its depth, and its hapNos
				push @{$$alleles_r{$allele_ID}}, ($allele,$currgroup_dep);
				push @{$$alleles_r{$allele_ID}}, @currgroup_hapNos;
				++$$n_alleles_r;#count the allele
				++$allele_ID;#increment for next allele
			}
			@ali = (); #set back
		}
	}
	#give next allele ID back to sub snphapcal1
	$$next_allele_ID_r = $allele_ID;
}

#definition of subroutine build_alleles_bin
#Used for binomial likelihood ratio method
#expects:
#ref to @seqmat: d1: rows, d2 cols, alignment of haplotypes of a locus
#ref to @varpos: variable positions in the locus
#ref to %hapgroups: {groupNo}=@of hapNos: hapNo: index to @seqmat
#ref to %hapgroup_dep: {groupNo} = depth of each hapgroup
#to populate/determine:
#ref to %alleles: {allele_ID} = @: [0]: allele_seq [1]: depth, remainder: hapNos
#ref to $n_alleles: number of alleles
#ref to @dischaplos: numbers of discarded haplos
#ref to $loclostall: number of lost alleles
#ref to $next_allele_ID: allele_ID of the first allele to be build here
#
#Collapses haplotypes into alleles.
#Discards alleles containing N - moves corresponding haplotypes into @dischaplos
#If more than 2 alleles remain:
#Sorts alleles according to read counts in descending order.
#Accepts the first two.
#Accepts further alleles that have same read count as second one.
#Discards other alleles.
#Moves corresponding haplotypes into @dischaplos.
sub build_alleles_bin {
	#declare and initialize: _r means a reference
	my ($seqmat_r,$varpos_r,$hapgroups_r,$hapgroup_dep_r,
	$alleles_r,$n_alleles_r,$dischaplos_r,
	$loclostall_r,$next_allele_ID_r) = @_;
	my @seqmat = @{$seqmat_r};#copy of outer array of @seqmat
	my %hapgroup_dep = %{$hapgroup_dep_r};# copy of %hapgroup_dep: {groupNo} = depth
	my %hapgroup_seq = ();# {groupNo} = allele sequence
	my $allele_ID = $$next_allele_ID_r;#allele ID for the next allele to be build
	my $currgroup = 0;#number of current hapgroup
	my @currgroup_hapNos = ();#haplotype-numbers of current group
	my $hapNo = 0;#number of a haplotype
	my @ali = ();#all seqs of one hapgroup: d1: rows, d2: cols
	my $allele = '';#an allele
	my @sortdep = ();#depths of all hapgroups in descending order
	my $i = 0;
	my @temparr = ();

	#loop through haplotype-groups and determine allele sequences
	for $currgroup (keys %{$hapgroups_r}) {
		@currgroup_hapNos = @{$$hapgroups_r{$currgroup}};#get all hapNos for this group
		#loop through group-members and build ali
		for $hapNo (@currgroup_hapNos) {
			#get sequence of this hapNo (as array)
			@temparr = @{$seqmat[$hapNo]};
			push @ali, [@temparr];
		}
		#call sub consensus1 to collapse the ali to an allele sequence
		$allele = consensus1(\@ali,$varpos_r);
		if ($allele =~ /N/) {#if the allele contains N
			++$$loclostall_r;#count a lost allele
			#add the haplotype numbers to discarded-haplotype-numbers
			push @{$dischaplos_r}, @currgroup_hapNos;
			#remove the haplotype group from %hapgroup_dep
			delete $hapgroup_dep{$currgroup};
		} else {#add the allele sequence to %hapgroup_seq
			$hapgroup_seq{$currgroup} = $allele;
		}
		@ali = ();#set back
	}
	if (keys %hapgroup_seq > 2) {#if there are more than 2 remaining alleles
		#sort read counts in descending order
		@sortdep = reverse sort {$a <=> $b} values %hapgroup_dep;
		#second value is minimum acceptable depth
		#loop through remaining hapgroups/alleles
		for $currgroup (keys %hapgroup_dep) {
			#if depth too small
			if ($hapgroup_dep{$currgroup} < $sortdep[1]) {
				delete $hapgroup_dep{$currgroup};#remove hapgroup
				delete $hapgroup_seq{$currgroup};#remove hapgroup
				push @{$dischaplos_r}, @{$$hapgroups_r{$currgroup}};#add hapNos to discarded haplotype numbers
				++$$loclostall_r;#count a lost allele
			}
		}		
	}
	#populate/determine %alleles and $n_alleles owned by sub snphapcalbin
	for $currgroup (keys %hapgroup_seq) {
		push @{$$alleles_r{$allele_ID}}, ($hapgroup_seq{$currgroup},$hapgroup_dep{$currgroup});
		push @{$$alleles_r{$allele_ID}}, @{$$hapgroups_r{$currgroup}};
		++$$n_alleles_r;#count the allele
		++$allele_ID;#increment for next allele
	}
	#give next allele ID back to sub snphapcal1
	$$next_allele_ID_r = $allele_ID;
}

#definition of sub consensus1
#expects:
#ref to @ali: d1: rows, d2: cols, a sequence-alignment
#ref to @varpos: variable positions in the alignment
#returns: $allele: consensus-sequence as string
#consensus:
#consensus1 assumes that there is only 1 character (ACGT) in each col
#possibly, there is N or N and 1 character (ACGT) in a col
#in each col where there is N and 1 other character:
#consensus contains that other character, if only N: consensus contains N
sub consensus1 {
	#declare and initialize: _r means a reference
	my ($ali_r, $varpos_r) = @_;
	my @ali = @{$ali_r};#copy of outer array of @ali
	my $allele = '';#an allele as string
	my @allelearr = ();#an allele as array
	my $row = 0;#a row
	my $pos = 0;#a position
	my %chars = ();#collects characters in a col
	my $char = '';#a character
	
	@allelearr = @{$ali[0]};#copy first sequence into array
	if (@ali > 1) {#if there is more than 1 sequence
		#loop through variable positions
		for $pos (@{$varpos_r}) {
			#loop through this col in @ali
			for ($row = 0; $row < @ali; ++$row) {
				$char = ${$ali[$row]}[$pos];#get current character
				++$chars{$char};#collect in %chars
			}
			if (keys %chars > 1) {#if there is more than 1 char
				delete $chars{"N"};#delete N
			}
			for $char (keys %chars) {#get the single remaining key as char
				$allelearr[$pos] = $char;#put it into this pos in @allelearr
			}
			%chars = ();#set back and next pos
		}
	}
	$allele = join ('',@allelearr);#join the allele into a string
	return $allele;
}
