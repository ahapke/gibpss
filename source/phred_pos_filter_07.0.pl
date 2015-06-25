#!/usr/bin/perl -w
#phred_pos_filter version 07.0 Copyright 2015 Andreas Hapke
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
use IO::Compress::Gzip qw(gzip $GzipError :constants);
use Parallel::ForkManager;

#uses subroutines
#filter_single
#filter_pair

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "phred_pos_filter version 07.0 Copyright 2015 Andreas Hapke\n",
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
my $outdir = 'pf_out';#name of outfile directory
my $filfile = 'filter.txt';#name of filter file
my $po = 33;#Phred score offset (33 for Sanger)
my $max_p = 1;#maximum number of processes in parallel
my $z = 0;#1: write gzipped fastq outfile, 0: don't
my %user_settings = (
i => 'fq_in.txt',#name of file with fastq-infile-names
o => 'pf_out',#name of outfile directory
f => 'filter.txt',#name of filter file
p => '33',#Phred score offset
P => '1',#maximum number of processes in parallel
z => '0'#1: write gzipped fastq outfile, 0: don't
);
my %defaults = (
i => 'fq_in.txt',#name of file with fastq-infile-names
o => 'pf_out',#name of outfile directory
f => 'filter.txt',#name of filter file
p => '33',#Phred score offset
P => '1',#maximum number of processes in parallel
z => '0'#1: write gzipped fastq outfile, 0: don't
);
my $n_arg = 0;#number of arguments (flags and values)
my $flag = '';#a flag
my $val = '';#a value

#other variables
my $fileNo = 0;#Number for fq-infile or pair of fq-infiles
my %fqin = ();# {$fileNo} = @ of infilenames (1 or 2)
my @fpos = ();#checkpositions in forward (or single) reads
my @fpmin = ();#minimum Phred scores in forward (or single read) at positions in @fpos
my @rpos = ();#checkpositions in reverse reads
my @rpmin = ();#minimum Phred scores in reverse at positions in @rpos
my $finfile = '';#forward (or single) fq-infile
my $rinfile = '';#reverse fq-infile
my $rep_single = '';#report file from single infile/pair of infiles
my $i = 0;#counter
my $tempstring1 = '';
my @temparr1 = ();
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
#z: 0 or 1
unless(($user_settings{'z'} eq '0') or ($user_settings{'z'} eq '1')) {
	$user_settings{'z'} = $defaults{'z'};
	print "z must be 0 or 1, using default: $defaults{'z'}\n";
}
#take over user settings from hash into scalar variables
$fqlist = $user_settings{'i'};
$outdir = $user_settings{'o'};
$filfile = $user_settings{'f'};
$po = $user_settings{'p'};
$max_p = $user_settings{'P'};
$z = $user_settings{'z'};
#}

###################################################################
#Create outdirectory, open report-file and files fqlist and filfile
###################################################################

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
unless(open(REPALL, ">$outdir/pf_rep.txt")) {#open or die
	print "Cannot open $outdir/pf_rep.txt, exiting ...\n\n";
	exit;
}
unless(open(FQLIST, "$fqlist")) {
	print "Cannot open file $fqlist. Exiting..\n";
	exit;
}
unless(open(FILFILE, "$filfile")) {
	print "Cannot open file $filfile. Exiting..\n";
	exit;
}
#}

#######################
#Read fqlist into %fqin
#######################

#{
while ($tempstring1 = <FQLIST>) {
	chomp $tempstring1;
	$tempstring1 =~ s/^\s*//;
	if (length $tempstring1 > 0) {
		@temparr1 = split(/\s+/,$tempstring1);
		if (@temparr1 > 2) {#if there are more than 2 filenames
			print "I do not understand this line in file $fqlist:\n $tempstring1\n Exiting..\n";
			exit;
		} else {#try to open the files
			for $tempstring1 (@temparr1) {
				unless(open(INFILE, "$tempstring1")) {
					print "Cannot open file $tempstring1. Exiting..\n";
					exit;
				}
				close INFILE;
			}
			@{$fqin{$fileNo}} = @temparr1;#take over fq-infilenames
			++$fileNo;
		}		
	}
}
close FQLIST;
#}

################
#Read in filfile
################

#{
while ($tempstring1 = <FILFILE>) {
	chomp $tempstring1;
	$tempstring1 =~ s/^\s*//;
	if (length $tempstring1 > 0) {
		@temparr1 = split(/\s+/,$tempstring1);
		#check: 3 entries, position: integer >0, Phred score min: integer >0
		if ((@temparr1 == 3) and ($temparr1[1] =~ /^\d+$/) and ($temparr1[1] > 0) and
		($temparr1[2] =~ /^\d+$/) and ($temparr1[2] > 0)) {
			--$temparr1[1];#decrement position: we start counting by 0, user by 1
			if ($temparr1[0] eq 'f') {
				push @fpos, $temparr1[1];
				push @fpmin, $temparr1[2];
			}
			elsif ($temparr1[0] eq 'r') {
				push @rpos, $temparr1[1];
				push @rpmin, $temparr1[2];
			}
			else {
				print " I do not understand this line in file $filfile:\n$tempstring1\nExiting..\n";
				exit;
			}
		}
		else {
			print " I do not understand this line in file $filfile:\n$tempstring1\nExiting..\n";
			exit;
		}
	}
}
close FILFILE;
#}

##############
#Data analysis
##############

#{
my $pm1 = Parallel::ForkManager->new($max_p);
for $fileNo (sort {$a <=> $b} keys %fqin) {
	$pm1->start and next;
	if (@{$fqin{$fileNo}} == 1) {#if there is a single fq-infile
		$finfile = ${$fqin{$fileNo}}[0];
		#message to screen:
		print "Analyzing $finfile\n";
		#call subroutine filter_single
		filter_single($fileNo,$finfile,$outdir,\@fpos,\@fpmin,$po,$z);
	}
	elsif (@{$fqin{$fileNo}} == 2) {#if there is a pair of fq-infiles
		$finfile = ${$fqin{$fileNo}}[0];
		$rinfile = ${$fqin{$fileNo}}[1];
		#message to screen:
		print "Analyzing $finfile and $rinfile\n";
		#call subroutine filter_pair
		filter_pair($fileNo,$finfile,$rinfile,$outdir,\@fpos,\@fpmin,\@rpos,\@rpmin,$po,$z);
	}
	$pm1->finish;
}
$pm1->wait_all_children;
#}

############################################################################
#Produce report and delete report files from single infiles/pairs of infiles
############################################################################

#{
#print header line:
print REPALL "f_infile\tn_seq\tn_badseq\tprop_badseq\tn_goodseq\n";
for $fileNo (sort {$a <=> $b} keys %fqin) {#loop through file-Nos
	#open report file for this fileNo
	$rep_single = $outdir . '/' . $fileNo . '_rep.txt';
	unless(open(INFILE, "$rep_single")) {
		print "Cannot open report file $rep_single\n",
			"Something went wrong - don't use the results. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#copy single line of outfile
	print REPALL "$tempstring1";
	close INFILE;
	unlink $rep_single;#delete report file for this fileNo
}
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";

exit;

############
#Subroutines
############

#Definition of subroutine filter_single
#filters a single fastq-infile
#expects:
#$fileNo: number of current fastq-infile
#$finfile: name of fastq-infile
#$outdir: name of directory for outfile
#ref to @fpos: filter positions
#ref to @fpmin: minimum Phred scores for each filter position
#$po: Phred score offset
#$z: 1: write gzipped fastq outfile, 0: don't
#checks each sequence in fastq-file:
#if Phred score is < minimum at any filter position: discards sequence
#else: stores sequence to out-fq-file in directory $outdir

sub filter_single {
	#declare and initialize, _r means a reference
	my($fileNo,$finfile,$outdir,$fpos_r,$fpmin_r,$po,$z) = @_;
	my @fpos = @{$fpos_r};
	my @fpmin = @{$fpmin_r};
	my $finfiledirname = '';#dirname of finfile
	my $finfilename = '';#name of finfile
	my $foutfile = '';#name of fastq-outfile
	my $rep = '';#name of report-file
	my $nseq = 0;#number of sequences in infile
	my $nbadseq = 0;#number of bad sequences in infile
	my $propbadseq = 0;#proportion of bad sequences in infile
	my $ngoodseq = 0;#number of good sequences in infile
	my $fseqline1 = '';
	my $fseqline2 = '';
	my $fseqline3 = '';
	my $fseqline4 = '';
	my $fbad = 0;
	my $psymbol = '';#a Phred symbol
	my $pscore = 0;#a Phred score
	my $i = 0;
	my $finfileh = 0;#handle for fastq infile
	my $foutfileh = 0;#handle for fastq outfile
	
	#determine fastq-outfilename and open or die
	($finfiledirname,$finfilename) = $finfile =~ m{^(.*[/\\])([^/\\]+?)$};
	$foutfile = $finfilename;
	$foutfile =~ s/\.[^\.]*$//;#remove file-extension, if any
	$foutfile = $outdir . '/' . $foutfile;
	if ($z == 1) {
		$foutfile .= '_pf.fq.gz';
	} else {
		$foutfile .= '_pf.fq';
	}
	if ($z == 1) {#if active open fastq outfile as gzipped file
		$foutfileh = new IO::Compress::Gzip $foutfile ,-Level=>Z_BEST_SPEED
			or die "Cannot write Gzip: $GzipError\n";
	} else {#open as text file
		$foutfileh = IO::File->new();
		unless(open($foutfileh,">",$foutfile)) {
			print "Could not open $foutfile: $!\n",
			"Exiting..\n";
			exit;
		}
	}	
	#determine report-filename and open or die
	$rep = $fileNo . '_rep.txt';
	unless(open(REPORT, ">$outdir/$rep")) {
		print "Cannot open outfile $outdir/$rep. Exiting..\n";
		exit;
	}
	if ($finfile =~ /.gz$/) {#if infile has extension .gz open as gz file
		$finfileh = new IO::Uncompress::Gunzip $finfile
			or die "Could not open $finfile: $GunzipError\n";
	} else {#open as text file
		$finfileh = IO::File->new();
		unless(open($finfileh, $finfile)) {
			print "Could not open $finfile: $!\n",
			"Exiting..\n";
			exit;
		}
	}
	#Analyze fastq-infile
	while ($fseqline1 = <$finfileh>) {
		$fseqline2 = <$finfileh>;
		$fseqline3 = <$finfileh>;
		$fseqline4 = <$finfileh>;
		chomp $fseqline1;
		chomp $fseqline2;
		chomp $fseqline3;
		chomp $fseqline4;
		++$nseq;#count sequence
		#check Phred scores
		$fbad = 0;
		for ($i = 0; $i < @fpos; ++$i) {#loop through filter positions
			$psymbol = substr($fseqline4,$fpos[$i],1);
			$pscore = ord($psymbol) - $po;#determine Phred score
			if ($pscore < $fpmin[$i]) {#if below threshold
				$fbad = 1;
				last;#stop checking
			}
		}
		if ($fbad == 0) {#if sequence was good
			++$ngoodseq;
			print $foutfileh "$fseqline1\n","$fseqline2\n","$fseqline3\n","$fseqline4\n";
		} else {#if sequence was bad
			++$nbadseq;
		}		
	}
	close $finfileh;
	close $foutfileh;
	#determine proportion of bad seqs
	$propbadseq = $nbadseq / $nseq;
	#print report to file and close report file
	print REPORT "$finfile\t$nseq\t$nbadseq\t$propbadseq\t$ngoodseq\n";
	close REPORT;	
}

#Definition of subroutine filter_pair
#filters a pair of fastq-infiles
#expects:
#$fileNo: number of current fastq-infile
#$finfile: name of f-fastq-infile
#$rinfile: name of r-fastq-infile
#$outdir: name of directory for outfile
#ref to @fpos: filter positions for f-infile
#ref to @fpmin: minimum Phred scores for each filter position in f-infile
#ref to @rpos: filter positions for r-infile
#ref to @rpmin: minimum Phred scores for each filter position in r-infile
#$po: Phred score offset
#$z: 1: write gzipped fastq outfiles, 0: don't
#
#checks each sequence-pair in fastq-files:
#if Phred score is < minimum at any filter position in any seq of pair: discards pair
#else: stores sequence-pair to out-fq-files in directory $outdir

sub filter_pair {
	#declare and initialize, _r means a reference
	my($fileNo,$finfile,$rinfile,$outdir,$fpos_r,$fpmin_r,$rpos_r,$rpmin_r,$po,$z) = @_;
	my @fpos = @{$fpos_r};
	my @fpmin = @{$fpmin_r};
	my @rpos = @{$rpos_r};
	my @rpmin = @{$rpmin_r};
	my $finfiledirname = '';#dirname of finfile
	my $finfilename = '';#name of finfile
	my $rinfiledirname = '';#dirname of rinfile
	my $rinfilename = '';#name of rinfile
	my $foutfile = '';#name of f-fastq-outfile
	my $routfile = '';#name of r-fastq-outfile
	my $rep = '';#name of report-file
	my $nseq = 0;#number of sequence pairs in infiles
	my $nbadseq = 0;#number of bad sequence pairs in infiles
	my $propbadseq = 0;#proportion of bad sequence pairs in infiles
	my $ngoodseq = 0;#number of good sequence pairs in infiles
	my $fseqline1 = '';
	my $fseqline2 = '';
	my $fseqline3 = '';
	my $fseqline4 = '';
	my $rseqline1 = '';
	my $rseqline2 = '';
	my $rseqline3 = '';
	my $rseqline4 = '';
	my $fbad = 0;
	my $rbad = 0;
	my $psymbol = '';#a Phred symbol
	my $pscore = 0;#a Phred score
	my $i = 0;
	my $finfileh = 0;#handle for fastq f infile
	my $rinfileh = 0;#handle for fastq r infile
	my $foutfileh = 0;#handle for fastq f outfile
	my $routfileh = 0;#handle for fastq r outfile
	
	#determine fastq-outfilenames and open or die
	($finfiledirname,$finfilename) = $finfile =~ m{^(.*[/\\])([^/\\]+?)$};
	$foutfile = $finfilename;
	$foutfile =~ s/\.[^\.]*$//;#remove file-extension, if any
	$foutfile = $outdir . '/' . $foutfile;
	($rinfiledirname,$rinfilename) = $rinfile =~ m{^(.*[/\\])([^/\\]+?)$};
	$routfile = $rinfilename;
	$routfile =~ s/\.[^\.]*$//;#remove file-extension, if any
	$routfile = $outdir . '/' . $routfile;
	if ($z == 1) {
		$foutfile .= '_pf.fq.gz';
		$routfile .= '_pf.fq.gz';
	} else {
		$foutfile .= '_pf.fq';
		$routfile .= '_pf.fq';
	}	
	#if names are equal now, make them different:
	unless($foutfile ne $routfile) {
		$foutfile = 'R1_' . $fileNo . '_' . $foutfile;
		$routfile = 'R2_' . $fileNo . '_' . $routfile;
	}
	if ($z == 1) {#if active open fastq outfiles as gzipped files
		$foutfileh = new IO::Compress::Gzip $foutfile ,-Level=>Z_BEST_SPEED
			or die "Cannot write Gzip: $GzipError\n";
		$routfileh = new IO::Compress::Gzip $routfile ,-Level=>Z_BEST_SPEED
			or die "Cannot write Gzip: $GzipError\n";
	} else {#open as text files
		$foutfileh = IO::File->new();
		unless(open($foutfileh,">",$foutfile)) {
			print "Could not open $foutfile: $!\n",
			"Exiting..\n";
			exit;
		}
		$routfileh = IO::File->new();
		unless(open($routfileh,">",$routfile)) {
			print "Could not open $routfile: $!\n",
			"Exiting..\n";
			exit;
		}
	}
	#determine report-filename and open or die
	$rep = $fileNo . '_rep.txt';
	unless(open(REPORT, ">$outdir/$rep")) {
		print "Cannot open outfile $outdir/$rep. Exiting..\n";
		exit;
	}
	if ($finfile =~ /.gz$/) {#if f infile has extension .gz open as gz file
		$finfileh = new IO::Uncompress::Gunzip $finfile
			or die "Could not open $finfile: $GunzipError\n";
	} else {#open as text file
		$finfileh = IO::File->new();
		unless(open($finfileh, $finfile)) {
			print "Could not open $finfile: $!\n",
			"Exiting..\n";
			exit;
		}
	}
	if ($rinfile =~ /.gz$/) {#if r infile has extension .gz open as gz file
		$rinfileh = new IO::Uncompress::Gunzip $rinfile
			or die "Could not open $rinfile: $GunzipError\n";
	} else {#open as text file
		$rinfileh = IO::File->new();
		unless(open($rinfileh, $rinfile)) {
			print "Could not open $rinfile: $!\n",
			"Exiting..\n";
			exit;
		}
	}
	
	#Analyze fastq-infiles
	while ($fseqline1 = <$finfileh>) {
		$fseqline2 = <$finfileh>;
		$fseqline3 = <$finfileh>;
		$fseqline4 = <$finfileh>;
		$rseqline1 = <$rinfileh>;
		$rseqline2 = <$rinfileh>;
		$rseqline3 = <$rinfileh>;
		$rseqline4 = <$rinfileh>;
		chomp $fseqline1;
		chomp $fseqline2;
		chomp $fseqline3;
		chomp $fseqline4;
		chomp $rseqline1;
		chomp $rseqline2;
		chomp $rseqline3;
		chomp $rseqline4;
		++$nseq;#count sequence-pair
		#check Phred scores
		$fbad = 0;
		$rbad = 0;
		#if user specified filters for forward reads, check:
		if (@fpos > 0) {
			for ($i = 0; $i < @fpos; ++$i) {#loop through filter positions
				$psymbol = substr($fseqline4,$fpos[$i],1);
				$pscore = ord($psymbol) - $po;#determine Phred score
				if ($pscore < $fpmin[$i]) {#if below threshold
					$fbad = 1;
					last;#stop checking
				}
			}
		}
		#if f-seq good, AND user specified filters for reverse reads, check:
		if (($fbad == 0) and (@rpos > 0)) {
			for ($i = 0; $i < @rpos; ++$i) {#loop through filter positions
				$psymbol = substr($rseqline4,$rpos[$i],1);
				$pscore = ord($psymbol) - $po;#determine Phred score
				if ($pscore < $rpmin[$i]) {#if below threshold
					$rbad = 1;
					last;#stop checking
				}
			}
		}
		if (($fbad == 0) and ($rbad == 0)) {#if both sequences were good
			++$ngoodseq;
			print $foutfileh "$fseqline1\n","$fseqline2\n","$fseqline3\n","$fseqline4\n";
			print $routfileh "$rseqline1\n","$rseqline2\n","$rseqline3\n","$rseqline4\n";			
		} else {#if any sequence was bad
			++$nbadseq;
		}		
	}
	close $finfileh;
	close $rinfileh;
	close $foutfileh;
	close $routfileh;
	#determine proportion of bad seqs
	$propbadseq = $nbadseq / $nseq;
	#print report to file and close report file
	print REPORT "$finfile\t$nseq\t$nbadseq\t$propbadseq\t$ngoodseq\n";
	close REPORT;	
}
