#!/usr/bin/perl -w
#phredi version 08.0 Copyright 2015 Andreas Hapke
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

#runtime parameters
my $start_time = time();
my $end_time = 0;
my $run_s = 0;

#Copyright message to screen
print "phredi version 08.0 Copyright 2015 Andreas Hapke\n",
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
#user settings and defaults
my $fqlist = 'fq_in.txt';#name of file with fastq-infile-names
my $outdir = 'phredi_out';#name of outfile directory
my $po = 33;#Phred score offset (33 for Sanger)
my $max_p = 1;#maximum number of processes in parallel
my %user_settings = (
i => 'fq_in.txt',#name of file with fastq-infile-names
o => 'phredi_out',#name of outfile directory
p => '33',#Phred score offset
P => '1',#maximum number of processes in parallel
);
my %defaults = (
i => 'fq_in.txt',#name of file with fastq-infile-names
o => 'phredi_out',#name of outfile directory
p => '33',#Phred score offset
P => '1',#maximum number of processes in parallel
);
my $n_arg = 0;#number of arguments (flags and values)
my $flag = '';#a flag
my $val = '';#a value

#other variables
my @fq_infiles = ();#names/pathes of fq_infiles
my $fq_infile = '';#name/path of one fq_infile
my $fq_infiledirname = '';#path of fq_infile
my $fq_infilename = '';#name of fq_infile
my $fq_infile_h = 0;#handle for fq_infile
my %outfiles = ();# {$fq_infile} = $outfile
my $outfile = '';#one outfilename including path
my $outfilesuffix = '_phredi.txt';
my $seql = 0;#sequence length (all seqs in one infile must have same length)
my $seqline1 = '';#seqline1 of sequence entry
my $seqline2 = '';#seqline2 of sequence entry
my $seqline3 = '';#seqline3 of sequence entry
my $seqline4 = '';#seqline4 of sequence entry
my $pos = 0;#sequence position, count starting with 0 for internal use
my $outpos = 0;#sequence position, count starting with 1 for output
my $phredsymbol = '';
my $phredscore = 0;#a phredscore
my $phredscoremin = 0;#smallest Phred score in an infile
my $phredscoremax = 0;#greatest Phred score in an infile
my %phredscores = ();#keys: all occurring phredscores in one infile
my %phredcounts = ();# {$pos}{$phredscore}= number of sequences with this Phred score
my $nseq = 0;#number of sequences in one infile
my %phredperc = ();# {$pos}{$phredscore}= percent of sequences with this Phred score
my $i = 0;
#}

########################
#Take over user settings
########################

#{
if (@ARGV) {#if user provided some
	$n_arg = @ARGV;
	for ($i = 0; $i <($n_arg - 1); $i += 2) {
		$flag = $ARGV[$i];
		$flag =~ s/^-//;
		$val = $ARGV[$i+1];
		if (defined $user_settings{$flag}) {#if flag is defined
			$user_settings{$flag} = $val;#take over with value
		}
	}
}
#check (at least some) arguments for plausibility (very limited test)
#Phred score offset: positive integer
unless (($user_settings{'p'} =~ /^\d+$/) and ($user_settings{'p'} > 0)) {
	$user_settings{'p'} = $defaults{'p'};
	print "-p must be a positive integer, using default: $defaults{'p'}\n";
}
#max_p: positive integer
unless(($user_settings{'P'} =~ /^\d+$/) and ($user_settings{'P'} > 0)) {
	$user_settings{'P'} = $defaults{'P'};
	print "P must be a positive integer, using default: $defaults{'P'}\n";
}
#take over user settings from hash into scalar variables
$fqlist = $user_settings{'i'};
$outdir = $user_settings{'o'};
$po = $user_settings{'p'};
$max_p = $user_settings{'P'};
#}

######################################
#Create outdirectory, open file fqlist
######################################

#{
#create out-directory or die
if (-d "$outdir") {#if directory already exists, die
	print "Directory $outdir already exists.\n",
	"Please delete or rename. Exiting..\n";
	exit;
} else {#if not, try to create
	unless (mkdir "$outdir") {
		print "Can't create directory $outdir. Exiting..\n";
		exit;
	}
}
unless(open(FQLIST, "$fqlist")) {
	print "Cannot open file $fqlist. Exiting..\n";
	exit;
}
#}

#################################################################
#Read in infilenames, determine outfilenames, try to open infiles
#################################################################

#{
while ($fq_infile = <FQLIST>) {
	chomp $fq_infile;
	$fq_infile =~ s/^s*//;
	if (length $fq_infile > 0) {
		#split fq_infilename into path and filename
		($fq_infiledirname,$fq_infilename) = $fq_infile =~ m{^(.*[/\\])([^/\\]+?)$};
		$outfile = $outdir . '/' . $fq_infilename . $outfilesuffix;
		push @fq_infiles, $fq_infile;
		$outfiles{$fq_infile} = $outfile;
	}
}
close FQLIST;
for $fq_infile (@fq_infiles) {
	if ($fq_infile =~ /.gz$/) {#if infile has extension .gz open as gz file
		$fq_infile_h = new IO::Uncompress::Gunzip $fq_infile
			or die "Could not open $fq_infile: $GunzipError\n";
	} else {#open as text file
		$fq_infile_h = IO::File->new();
		unless(open($fq_infile_h, $fq_infile)) {
			print "Could not open $fq_infile_h: $!\n",
			"Exiting..\n";
			exit;
		}
	}
	close $fq_infile_h;
}
undef $fq_infile_h;
#}

##############
#Data analysis
##############

#{
my $pm = Parallel::ForkManager->new($max_p);
for $fq_infile (@fq_infiles) {
	$pm->start and next;
	#message to screen
	print "Analyzing $fq_infile\n";
	#determine sequence length of first sequence in infile
	#all sequences in the files MUST have same length
	if ($fq_infile =~ /.gz$/) {#if infile has extension .gz open as gz file
		$fq_infile_h = new IO::Uncompress::Gunzip $fq_infile
			or die "Could not open $fq_infile: $GunzipError\n";
	} else {#open as text file
		$fq_infile_h = IO::File->new();
		unless(open($fq_infile_h, $fq_infile)) {
			print "Could not open $fq_infile_h: $!\n",
			"Exiting..\n";
			exit;
		}
	}
	$seqline1 = <$fq_infile_h>;
	$seqline2 = <$fq_infile_h>;
	chomp $seqline2;
	$seql = length $seqline2;
	close $fq_infile_h;
	#reopen infile
	if ($fq_infile =~ /.gz$/) {#if infile has extension .gz open as gz file
		$fq_infile_h = new IO::Uncompress::Gunzip $fq_infile
			or die "Could not open $fq_infile: $GunzipError\n";
	} else {#open as text file
		$fq_infile_h = IO::File->new();
		unless(open($fq_infile_h, $fq_infile)) {
			print "Could not open $fq_infile_h: $!\n",
			"Exiting..\n";
			exit;
		}
	}
	#open corresponding outfile
	$outfile = $outfiles{$fq_infile};
	unless(open(OUTFILE, ">$outfile")) {
		print "Cannot open $outfile. Exiting..\n";
		exit;
	}
	#analyze infile
	while ($seqline1 = <$fq_infile_h>) {
		$seqline2 = <$fq_infile_h>;
		$seqline3 = <$fq_infile_h>;
		$seqline4 = <$fq_infile_h>;#get Phred symbol line
		chomp $seqline4;
		#loop through characters in Phred symbol line
		for ($pos = 0; $pos < $seql; ++$pos) {
			$phredsymbol = substr($seqline4,$pos,1);
			$phredscore = ord($phredsymbol) - $po;#determine Phred score
			++$phredcounts{$pos}{$phredscore};#count Phred score			
		}
	}
	close $fq_infile_h;
	#determine number of sequences
	$pos = 0;
	$nseq = 0;
	for $phredscore (keys %{$phredcounts{$pos}}) {
		$nseq += $phredcounts{$pos}{$phredscore};
	}
	#exit if infile was empty
	if ($nseq == 0) {
		print "File $fq_infile did not contain sequences. Exiting..\n";
		exit;
	}
	#determine % of sequences with each phredscore at each position
	for $pos (keys %phredcounts) {
		for $phredscore (keys %{$phredcounts{$pos}}) {
			$phredperc{$pos}{$phredscore} = $phredcounts{$pos}{$phredscore} / $nseq * 100;
		}
	}
	#determine smallest and greatest occurring Phred score
	for $pos (keys %phredcounts) {
		for $phredscore (keys %{$phredcounts{$pos}}) {
			++$phredscores{$phredscore};
		}
	}
	$phredscoremin = (sort {$a <=> $b} keys %phredscores)[0];
	$phredscoremax = (reverse sort {$a <=> $b} keys %phredscores)[0];
	#print header line to outfile
	for ($phredscore = $phredscoremin; $phredscore <= $phredscoremax; ++$phredscore) {
		print OUTFILE "\tph$phredscore";
	}
	print OUTFILE "\n";
	#print data to outfile
	for $pos (sort {$a <=> $b} keys %phredperc) {
		$outpos = $pos + 1;
		print OUTFILE "pos$outpos";
		for ($phredscore = $phredscoremin; $phredscore <= $phredscoremax; ++$phredscore) {
			if (defined $phredperc{$pos}{$phredscore}) {
				print OUTFILE "\t$phredperc{$pos}{$phredscore}";
			} else {
				print OUTFILE "\t0";
			}
		}
		print OUTFILE "\n";
	}
	close OUTFILE;
	#set variables back and go to next fq-infile
	%phredcounts = ();
	%phredperc = ();
	%phredscores = ();
	$pm->finish;
}
$pm->wait_all_children;
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "run took $run_s seconds.\n";

exit;
