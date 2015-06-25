#!/usr/bin/perl -w
#fdm version 07.0 Copyright 2015 Andreas Hapke
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
use Parallel::ForkManager;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError :constants);

#uses subroutines
#ambi_var
#tree
#ch_h
#sliding_w
#r1_c
#r1_rep
#bcr_c
#bc_rrep
#bcrep_r
#bcrep_rrep
#cl_bc
#rint_c
#ad_c_pos
#cl_r
#trunc
#rcdupa
#catsrc

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "fdm version 07.0 Copyright 2015 Andreas Hapke\n",
"This program is part of GIbPSs Copyright 2015 Andreas Hapke\n",
"GIbPSs comes with ABSOLUTELY NO WARRANTY. GIbPSs is free software,\n",
"and you are welcome to redistribute it under certain conditions.\n",
"For details see the GNU General Public License version 3.\n",
"You should have obtained a copy of the GNU General Public License\n",
"along with GIbPSs. If not, see <http://www.gnu.org/licenses/>.\n\n";

##############
#main_settings
##############

#{
my ($fastq_listname) = 'fastq_list.txt';#the list of fastq-infilenames
										#with single reads: one filename per line
										#with paired reads: each line: forward-filename whitespace reverse-filename
my ($barcode_listname) = 'barcodes.txt';#the file with the barcodes, one per line
my ($bcad_listname) = 'bc_adapters.txt';#file with several barcode adapters to search for in reverse reads
										#each line: forward-barcode whitespace barcode-adapter in reverse read (revcomp)
my ($report_filename)  = 'report.txt';#report file with analysis result counts
my ($P) = 0;#1: paired fastq-infiles with paired reads, 0: single fastq-infiles
my ($save_f) = 0;#1: save forward reads
my ($save_r) = 0;#1: save reverse reads
my ($save_d) = 0;#1: save discarded reads
my ($ch_h) = 0;#1: change headerline of all sequences to simpler format, 0: don't
my ($po) = 33;#Phred score offset (here: 33 for Sanger)
my ($wl) = 10;#window length for sliding window
my ($pt) = 10;#avg Phred score threshold in window
my ($r_reppos) = 0;#number of positions in restr_site1 that may be repaired
my ($bc_reppos) = 0;#number of positions in barcode that may be repaired
my ($compbcl) = 0;#1: compensate for different barcode lengths by truncation of sequences with short barcodes, 0: don't
my ($lonbc) = 0;#length of longest barcode (in the whole study)
my ($admis) = 0.1;#fraction of mismatch allowed when searching for adapter
					#corresponding number of mismatch pos is rounded down to nearest integer
					#e.g. 0.1 in 19 positions evaluates to 1 mismatch allowed
my ($cropia) = 0;#crop insert of adapter containing sequences
my ($ilamin) = 0;#minimum length of (cropped) insert of adapter containing seqs
				#if cropped insert would be shorter, read or pair of reads is discarded
my ($rcdupa) = 0;#adapter contaning seq: replace r-read with reverse complement of f-read or vice-versa
my ($catsrc) = 0;#paired reads without adapter:
				#append reverse complement of r-read to f-read
				#replace r-read with rev.comp. of new concatenated f-read
my ($adsf) = 0;#1:save adapter containing sequences to separate outfiles,0: don't
my ($r1r2sf) = 0;#1: save paired reads (f- and r-) to separate outfiles,0: don't
my ($z) = 0;#1: write gzip compressed fastq outfiles, 0: write plain text fastq outfiles
my ($max_p) = 0;#maximum number of processes for parallel execution
#}

####################
#forwardseq_settings
####################

#{
my ($f_trf) = 0;#1: truncate first to specified length before any other analysis or filtering is done
my ($f_restr_site1) = 'GATCC';#restriction site at beginning of sequence or after barcode
my ($f_q) = 0;#1 do sliding window quality check, 0: don't
my ($f_bc) = 0;#1: barcodes are present at beginning of seq, 0: no barcodes
my ($f_bc_rrep) = 0;#1: barcode present, repair restr_site1 if necessary, 0: don't
my ($f_bcrep_r) = 0;#1: repair barcode if necessary, 0: don't
my ($f_bcrep_rrep) = 0;#1: repair barcode and restr_site1 if necessary, 0: don't
my ($f_cl_bc) = 0;#1: clip barcodes after succesful identification of barcode+restr
my ($f_r1_c) = 0;#no barcodes present: 1: check restr_site1 at beginning of seq
				 #0: don't
my ($f_r1_rep) = 0;#no barcodes present: 1: repair restr_site1 at beginning of seq if necessary
				#0: don't
my ($f_rint_c) = 0;#1: check for internal restriction site, 0: don't
my ($f_restr_site_int) = 'GGATCC';#internal restriction site, complete sequence
my ($f_ad_c) = 1;#1: check for adapter at sequence end, 0: don't
my ($f_adseq) = 'GATCAGATCGGAAGAGCACACGTCTGAACTCCAGTCA';#adapter sequence to search for
my ($f_ao) = 0;#offset for adapter clipping
my ($f_cl_r) = 0;#1: clip restriction site1, 0: don't
my ($f_tr) = 0;#0: don't truncate, positive integer: truncate to that value from the end
#}

####################
#reverseseq_settings
####################

#{
my ($r_trf) = 0;#1: truncate first to specified length before any other analysis or filtering is done
my ($r_restr_site1) = 'GATCC';#restriction site at beginning of sequence or after barcode
my ($r_q) = 0;#1 do sliding window quality check, 0: don't
my ($r_r1_c) = 0;#no barcodes present: 1: check restr_site1 at beginning of seq
				#0: don't
my ($r_r1_rep) = 0;#no barcodes present: 1: repair restr_site1 at beginning of seq if necessary
				#0: don't
my ($r_rint_c) = 0;#1: check for internal restriction site, 0: don't
my ($r_restr_site_int) = 'GGATCC';#internal restriction site, complete sequence
my ($r_sad_c) = 0;#1: check for single adapter at sequence end, 0: don't
my ($r_sadseq) = 'GGATCAGATCGGAAGAGCACACGTCTGAACTCCAGTCA';#single adapter sequence to search for
my ($r_bcad_c) = 0;#check for one of several barcode adapters at sequence end, 0: don't
my ($r_ao) = 0;#offset for adapter clipping
my ($r_cl_r) = 0;#1: clip restriction site1, 0: don't
my ($r_tr) = 0;#0: don't truncate, positive integer: truncate to that value from the end
#}

###########################################################
#change settings according to settingsfile provided by user
###########################################################

#{
#declare a hash of references to all setting variables:
my %settingsrefs = (
fastq_listname => \$fastq_listname,
barcode_listname => \$barcode_listname,
bcad_listname => \$bcad_listname,
report_filename => \$report_filename,
compbcl => \$compbcl,
lonbc => \$lonbc,
cropia => \$cropia,
ilamin => \$ilamin,
rcdupa => \$rcdupa,
catsrc => \$catsrc,
adsf => \$adsf,
r1r2sf => \$r1r2sf,
P => \$P,
save_f => \$save_f,
save_r => \$save_r,
save_d => \$save_d,
z => \$z,
ch_h => \$ch_h,
po => \$po,
wl => \$wl,
pt => \$pt,
r_reppos => \$r_reppos,
bc_reppos => \$bc_reppos,
admis => \$admis,
max_p => \$max_p,
f_trf => \$f_trf,
f_restr_site1 => \$f_restr_site1,
f_q => \$f_q,
f_bc => \$f_bc,
f_bc_rrep => \$f_bc_rrep,
f_bcrep_r => \$f_bcrep_r,
f_bcrep_rrep => \$f_bcrep_rrep,
f_cl_bc => \$f_cl_bc,
f_r1_c => \$f_r1_c,
f_r1_rep => \$f_r1_rep,
f_rint_c => \$f_rint_c,
f_restr_site_int => \$f_restr_site_int,
f_ad_c => \$f_ad_c,
f_adseq => \$f_adseq,
f_ao => \$f_ao,
f_cl_r => \$f_cl_r,
f_tr => \$f_tr,
r_trf => \$r_trf,
r_restr_site1 => \$r_restr_site1,
r_q => \$r_q,
r_r1_c => \$r_r1_c,
r_r1_rep => \$r_r1_rep,
r_rint_c => \$r_rint_c,
r_restr_site_int => \$r_restr_site_int,
r_sad_c => \$r_sad_c,
r_sadseq => \$r_sadseq,
r_bcad_c => \$r_bcad_c,
r_ao => \$r_ao,
r_cl_r => \$r_cl_r,
r_tr => \$r_tr,
);
#if the user did not provide an argument, exit:
unless(@ARGV) {
	print "Please provide a settings-filename as argument. Exiting..\n\n";
	exit;
}
my $settingsfilename = $ARGV[0];#get the settings-filename
#Open settingsfile or die
unless(open(INFILE,$settingsfilename)) {
	print "Cannot open file $settingsfilename. Exiting..\n\n";
	exit;
}
#Read in the file
my $setline = '';
my @setlinearr = ();
while ($setline = <INFILE>) {
	chomp $setline;
	$setline =~ s/#.*//;#remove comment
	$setline =~ s/^\s*//;
	if (length $setline > 0) {#if current line is not empty now
		@setlinearr = split('\s+',$setline);
		if (@setlinearr > 2) {#if there are more than 2 entries (key and value)
			print "I don't understand this line in file $settingsfilename:\n@setlinearr\nExiting..\n\n";
			exit;
		}
		unless($settingsrefs{$setlinearr[0]}){#if the argument is unknown
			print "File $settingsfilename:\nI don't know this variable: $setlinearr[0]\nExiting..\n\n";
			exit;
		}
		#execute the setting of this variable:
		if (@setlinearr == 2) {#if there is an argument and a value
			${$settingsrefs{$setlinearr[0]}} = $setlinearr[1];
		}
		if (@setlinearr == 1) {#if no value is provided for an argument
			${$settingsrefs{$setlinearr[0]}} = '';
		}
	}
}
close INFILE;
#}

#######################################
#declare and initialize other variables
#######################################

#{
my (@fastq_f_infiles) = ();#array of fastq-infilenames for forward reads
my (@fastq_r_infiles) = ();#array of fastq-infilenames for reverse reads
my ($r_bc) = 0;#0: no barcodes present at begin of reverse seqs
my (@barcodes) = ();#barcodes in entry order
my ($n_barcodes) = 0;#number of barcodes
my (%barcodes_restr) = ();# {barcode}= @ of matching variants to $f_restr_site1 preceded by this barcode
my (@f_adseq) = ();#all matching variants of $f_adseq (forwardseq_settings)
my ($r_ad_c) = 0;#1: search for adapter in reverse sequence, 0: don't
my (@r_adseq) = ();#all matching variants of adapter to search in reverse seq
my (%bc_adapters) = ();# {barcode}= @ of matching variants of corresponding barcode-adapter
						#barcode: as it appears at beginning of forward read
						#barcode-adapter-variants: as they appear in reverse read (revcomp)
my (@f_restr_site1) = ();#all matching variants to $f_restr_site1 (forwardseq_settings)
my (@r_restr_site1) = ();#all matching variants to $r_restr_site1 (reverseseq_settings)
my ($f_restr_site1_l) = 0;#length of forward restr_site1
my ($r_restr_site1_l) = 0;#length of reverse restr_site1
my (@f_restr_site_int) = ();#all matching variants to $f_restr_site_int (forwardseq_settings)
my (@r_restr_site_int) = ();#all matching variants to $r_restr_site_int (reverseseq_settings)
my (@f_destinations) = ();#all destinations for forward seqs that may exist according to settings
my ($f_destination) = '';#destination for current f-seq
my (@r_destinations) = ();#all destinations for reverse seqs that may exist according to settings
my ($r_destination) = '';#destination for current r-seq
my ($outfolder) = '';#name of folder for final outfiles
my $outfilename = '';#an outfilename
my (@f_outfilenames) = ();#all names of final fastq-outfiles for forward seqs
my (@I_f_outfilenames) = ();#one f-outfilename for destination+sequence for each fastq_f_infile in infileorder
my ($I_f_outfilename) = '';#f-outfilename for destination+sequence for one fastq_f_infile
my ($f_outfilesuffix) = '';#suffix for foward outfilenames
my (%f_outfilehandles) = ();#keys: destinations, values: filehandles for forward-outfiles
my ($I_f_outfilehandle) = '';#handle for destination+sequence outfile from single fastq_f_infile
my ($f_outfilehandle) = '';#handle for current f_outfile
my (@r_outfilenames) = ();#all names of final fastq-outfiles for reverse seqs
my (@I_r_outfilenames) = ();#one r-outfilename for destination+sequence for each fastq_r_infile in infileorder
my ($I_r_outfilename) = '';#r-outfilename for destination+sequence for one fastq_r_infile
my ($r_outfilesuffix) = '';#suffix for reverse outfilenames
my (%r_outfilehandles) = ();#keys: destinations, values: filehandles for reverse-outfiles
my ($I_r_outfilehandle) = '';#handle for destination+sequence outfile from single fastq_r_infile
my ($r_outfilehandle) = '';#handle for current r_outfile
my ($fr_outfilesuffix) = '';#suffix for outfiles than contain forward- AND reverse reads in the same file (if r1r2sf==0)
my (@fr_outfilenames) = ();#names of outfiles than contain forward- AND reverse reads in the same file (if r1r2sf==0)
my (@fr_destinations) = ();#if r1r2sf==0: destinations of paired reads (f-,r-) that will be saved in the same outfile
my ($fseqline1) = '';#forward-seqline1
my ($fseqline2) = '';#forward-seqline2
my ($fseqline3) = '';#forward-seqline3
my ($fseqline4) = '';#forward-seqline4
my ($rseqline1) = '';#reverse-seqline1
my ($rseqline2) = '';#reverse-seqline2
my ($rseqline3) = '';#reverse-seqline3
my ($rseqline4) = '';#reverse-seqline4
my ($f_read_No) = '1';#f-reads get number one when ch_h == 1
my ($r_read_No) = '2';#r-reads get number two when ch_h == 1
my (%f_tree_arguments) = ();#arguments for decision tree for forward-seqs
my (%r_tree_arguments) = ();#arguments for decision tree for reverse-seqs
my (@f_tree_result) = ();#results from decision tree for current forward seq
						#0: goods (1 good, 0 bad)
						#1: dest (destination)
						#2: ana (analysis result)
						#3: f_override: overrides r_tree results (none, rint, ad)
my (@r_tree_result) = ();#results from decision tree for current reverse seq
						#0: goods (1 good, 0 bad)
						#1: dest (destination)
						#2: ana (analysis result)
						#3: r_override: overrides f_tree results (none, rint, ad)
my (%f_analysis_counts) = ();#counts of goods, dest and ana for forward seqs
my (@f_analysis_results) = ();#all occuring names of analysis results for f-seqs
my (%f_destination_counts) = ();#f-destination counts
my (@f_destination_results) = ();#all destinations, forward-seqs have been sent to
my (%r_analysis_counts) = ();#counts of goods, dest and ana for reverse seqs
my (@r_analysis_results) = ();#all occuring names of analysis results for r-seqs
my (%r_destination_counts) = ();#r-destination counts
my (@r_destination_results) = ();#all destinations, reverse-seqs have been sent to
my (@I_f_ana_reportfilenames) = ();#filenames for analysis counts for each f-infile
my ($I_f_ana_reportfilename) = '';#filename for analysis counts for one f infile
my (@I_f_dest_reportfilenames) = ();#filenames for destination counts for each f-infile
my ($I_f_dest_reportfilename) = '';#filename for destination counts for one f infile
my (@I_r_ana_reportfilenames) = ();#filenames for analysis counts for each r-infile
my ($I_r_ana_reportfilename) = '';#filename for analysis counts for one r infile
my (@I_r_dest_reportfilenames) = ();#filenames for destination counts for each r-infile
my ($I_r_dest_reportfilename) = '';#filename for destination counts for one r infile
my ($testfh) = 0;#holds a filehandle for open tests
my ($f_fastqin_h) = 0;#handle for forward fastq infile
my ($r_fastqin_h) = 0;#handle for reverse fastq infile
my (@tempinfile) = ();#array that temporarily holds an infile
my (@one_line_arr) = ();#array that holds one line of an infile, e.g. split by whitespace
my ($i) = 0;#counter
my ($j) = 0;#counter
my ($key1) = '';#a key
my ($value) = '';#a value
my ($tempstring) = '';#a string
#}

#####################################			
#determine variables for forward-seqs
#####################################

#{
#determine f_restr_site1_l: length of f_restr_site1
$f_restr_site1_l = length $f_restr_site1;
#determine all matching variants of $f_restr_site1 (one if without ambiguity)
@f_restr_site1 = ambi_var($f_restr_site1);
#determine all matching variants of $f_restr_site_int (one if without ambiguity)
@f_restr_site_int = ambi_var($f_restr_site_int);

#get barcodes if user provided a file with barcodes
if ($f_bc == 1) {
	unless (open(INFILE,$barcode_listname)) {
		print "Cannot open file \"$barcode_listname\".";
		exit;
	}
	#Read in the file
	@barcodes = <INFILE>;
	$n_barcodes = @barcodes;#determine number of barcodes
	for ($i = 0; $i < $n_barcodes; ++$i) {
		chomp $barcodes[$i];
	}
	close INFILE;
	#populate %barcodes_restr
	for ($i = 0; $i < @barcodes; ++$i) {
		$key1 = $barcodes[$i];
		for ($j = 0; $j < @f_restr_site1; ++$j) {
			$value = $key1.$f_restr_site1[$j];
			push @{$barcodes_restr{$key1}}, $value;
		}
	}
}
#determine all matching variants of $f_adseq (forwardseq_settings)
@f_adseq = ambi_var($f_adseq);
#}

#####################################			
#determine variables for reverse-seqs
#####################################

#{
if ($P == 1) {#if sequences are paired
	#determine r_restr_site1_l: length of r_restr_site1
	$r_restr_site1_l = length $r_restr_site1;
	#determine all matching variants of $r_restr_site1 (one if without ambiguity)
	@r_restr_site1 = ambi_var($r_restr_site1);
	#determine all matching variants of $r_restr_site_int (one if without ambiguity)
	@r_restr_site_int = ambi_var($r_restr_site_int);
	#if search for single adapter sequence active
	if ($r_sad_c == 1) {
		#determine all matching variants of $r_sadseq (reverseseq_settings)
		@r_adseq = ambi_var($r_sadseq);
	}
	#get barcodes and corresponding barcode adapters if provided
	#open the file containing them
	if ($r_bcad_c == 1) {
		unless (open(INFILE,$bcad_listname)) {
			print "Cannot open file \"$bcad_listname\".";
			exit;
		}
		@tempinfile = <INFILE>;
		close INFILE;
		#determine all matching variant of each barcode adapter (1 if without ambiguity)
		#store barcodes and corresponding barcode-adapter-variants in:
		#%bc_adapters: {bc} = @ of all matching variants of corresponding barcode
		for ($i = 0; $i < @tempinfile; ++$i) {
			chomp $tempinfile[$i];
			$tempinfile[$i] =~ s/^\s*//;
			@one_line_arr = split('\s+',$tempinfile[$i]);
			@{$bc_adapters{$one_line_arr[0]}} = ambi_var($one_line_arr[1]);
		}
	}
}
#}

############################################################
#put references to tree-arguments for forward-seqs into hash
############################################################

#{
$f_tree_arguments{"seqline1"} = \$fseqline1;
$f_tree_arguments{"seqline2"} = \$fseqline2;
$f_tree_arguments{"seqline4"} = \$fseqline4;
$f_tree_arguments{"trf"} = \$f_trf;
$f_tree_arguments{"q"} = \$f_q;
$f_tree_arguments{"bc"} = \$f_bc;
$f_tree_arguments{"bc_rrep"} = \$f_bc_rrep;
$f_tree_arguments{"bcrep_r"} = \$f_bcrep_r;
$f_tree_arguments{"bcrep_rrep"} = \$f_bcrep_rrep;
$f_tree_arguments{"cl_bc"} = \$f_cl_bc;
$f_tree_arguments{"compbcl"} = \$compbcl;
$f_tree_arguments{"lonbc"} = \$lonbc;
$f_tree_arguments{"r1_c"} = \$f_r1_c;
$f_tree_arguments{"r1_rep"} = \$f_r1_rep;
$f_tree_arguments{"rint_c"} = \$f_rint_c;
$f_tree_arguments{"ad_c"} = \$f_ad_c;
$f_tree_arguments{"cropia"} = \$cropia;
$f_tree_arguments{"ilamin"} = \$ilamin;
$f_tree_arguments{"cl_r"} = \$f_cl_r;
$f_tree_arguments{"wl"} = \$wl;
$f_tree_arguments{"pt"} = \$pt;
$f_tree_arguments{"po"} = \$po;
$f_tree_arguments{"barcodes"} = \@barcodes;
$f_tree_arguments{"n_barcodes"} = \$n_barcodes;
$f_tree_arguments{"barcodes_restr"} = \%barcodes_restr;
$f_tree_arguments{"restr_site1"} = \@f_restr_site1;
$f_tree_arguments{"restr_site1_l"} = \$f_restr_site1_l;
$f_tree_arguments{"r_reppos"} = \$r_reppos;
$f_tree_arguments{"bc_reppos"} = \$bc_reppos;
$f_tree_arguments{"restr_site_int"} = \@f_restr_site_int;
$f_tree_arguments{"adseq"} = \@f_adseq;
$f_tree_arguments{"admis"} = \$admis;
$f_tree_arguments{"ao"} = \$f_ao;
#}

############################################################
#put references to tree-arguments for reverse-seqs into hash
############################################################

#{
$r_tree_arguments{"seqline1"} = \$rseqline1;
$r_tree_arguments{"seqline2"} = \$rseqline2;
$r_tree_arguments{"seqline4"} = \$rseqline4;
$r_tree_arguments{"trf"} = \$r_trf;
$r_tree_arguments{"q"} = \$r_q;
$r_tree_arguments{"bc"} = \$r_bc;
$r_tree_arguments{"r1_c"} = \$r_r1_c;
$r_tree_arguments{"r1_rep"} = \$r_r1_rep;
$r_tree_arguments{"rint_c"} = \$r_rint_c;
$r_tree_arguments{"ad_c"} = \$r_ad_c;
$r_tree_arguments{"cropia"} = \$cropia;
$r_tree_arguments{"ilamin"} = \$ilamin;
$r_tree_arguments{"cl_r"} = \$r_cl_r;
$r_tree_arguments{"wl"} = \$wl;
$r_tree_arguments{"pt"} = \$pt;
$r_tree_arguments{"po"} = \$po;
$r_tree_arguments{"restr_site1"} = \@r_restr_site1;
$r_tree_arguments{"restr_site1_l"} = \$r_restr_site1_l;
$r_tree_arguments{"r_reppos"} = \$r_reppos;
$r_tree_arguments{"restr_site_int"} = \@r_restr_site_int;
$r_tree_arguments{"adseq"} = \@r_adseq;
$r_tree_arguments{"admis"} = \$admis;
$r_tree_arguments{"ao"} = \$r_ao;
#}

##########################
#get fastq-infilenames
#Try to open fastq-infiles
##########################

#{
#open file with fastq-infilenames and read in
unless (open(INFILE,$fastq_listname)) {
	print "Cannot open file \"$fastq_listname\".";
	exit;
}
@tempinfile = <INFILE>;
close INFILE;
for ($i = 0; $i < @tempinfile; ++$i) {
	chomp $tempinfile[$i];
	$tempinfile[$i] =~ s/^\s*//;
}
#if there are no paired infiles, put forward-infilenames into array
#Try to open forward infiles
if ($P == 0) {
	@fastq_f_infiles = @tempinfile;
	for $tempstring (@fastq_f_infiles) {
		if ($tempstring =~ /.gz$/) {#if infile has extension .gz open as gz file
			$testfh =  new IO::Uncompress::Gunzip $tempstring
				or die "Could not open $tempstring: $GunzipError\n";
		} else {#open as text file
			$testfh = IO::File->new("< $tempstring")
				or die "Could not open $tempstring: $!\n";
		}
		close $testfh;
	}
undef $testfh;
}
#if there are paired infiles, put forward and reverse infilenames into separate arrays
#Try to open forward and reverse infiles
if ($P == 1) {
	for ($i = 0; $i < @tempinfile; ++$i) {
		@one_line_arr = split('\s+',$tempinfile[$i]);
		$fastq_f_infiles[$i] = $one_line_arr[0];
		$fastq_r_infiles[$i] = $one_line_arr[1];
	}
	for $tempstring (@fastq_f_infiles) {
		if ($tempstring =~ /.gz$/) {#if infile has extension .gz open as gz file
			$testfh =  new IO::Uncompress::Gunzip $tempstring
				or die "Could not open $tempstring: $GunzipError\n";
		} else {#open as text file
			$testfh = IO::File->new("< $tempstring")
				or die "Could not open $tempstring: $!\n";
		}
		close $testfh;
	}
	for $tempstring (@fastq_r_infiles) {
		if ($tempstring =~ /.gz$/) {#if infile has extension .gz open as gz file
			$testfh =  new IO::Uncompress::Gunzip $tempstring
				or die "Could not open $tempstring: $GunzipError\n";
		} else {#open as text file
			$testfh = IO::File->new("< $tempstring")
				or die "Could not open $tempstring: $!\n";
		}
		close $testfh;
	}
undef $testfh;
}
#}

#################################
#create folder for final outfiles
#################################

#{
#determine name of outfolder
$outfolder = $settingsfilename;#take over name of settingsfile
$outfolder =~ s/\.[^\.]*$//;#remove file extension if existing
$outfolder .= '_out';
#create outfolder or die
unless (mkdir $outfolder) {
	print "Cannot create directory $outfolder. Exiting..\n\n";
	exit;
}
#}

################################################
#print used settings to file "used_settings.txt"
################################################

#{
my $printsettings = 1;
if ($printsettings == 1) {
	unless (open(OUTFILE,">$outfolder/used_settings.txt")){
		print "cannot print your settings to used_settings.txt, Exiting..";
		exit;
	}
	print OUTFILE
"fastq_listname $fastq_listname
barcode_listname $barcode_listname
bcad_listname $bcad_listname
report_filename $report_filename
P $P
save_f $save_f
save_r $save_r
save_d $save_d
z $z
ch_h $ch_h
po $po
wl $wl
pt $pt
r_reppos $r_reppos
bc_reppos $bc_reppos
compbcl $compbcl
lonbc $lonbc
admis $admis
cropia $cropia
ilamin $ilamin
rcdupa $rcdupa
catsrc $catsrc
adsf $adsf
r1r2sf $r1r2sf
max_p $max_p
f_trf $f_trf
f_restr_site1 $f_restr_site1
f_q $f_q
f_bc $f_bc
f_bc_rrep $f_bc_rrep
f_bcrep_r $f_bcrep_r
f_bcrep_rrep $f_bcrep_rrep
f_cl_bc $f_cl_bc
f_r1_c $f_r1_c
f_r1_rep $f_r1_rep
f_rint_c $f_rint_c
f_restr_site_int $f_restr_site_int
f_ad_c $f_ad_c
f_adseq $f_adseq
f_ao $f_ao
f_cl_r $f_cl_r
f_tr $f_tr
r_trf $r_trf
r_restr_site1 $r_restr_site1
r_q $r_q
r_r1_c $r_r1_c
r_r1_rep $r_r1_rep
r_rint_c $r_rint_c
r_restr_site_int $r_restr_site_int
r_sad_c $r_sad_c
r_sadseq $r_sadseq
r_bcad_c $r_bcad_c
r_ao $r_ao
r_cl_r $r_cl_r
r_tr $r_tr";
}
close OUTFILE;
#}

##################################################################################
#if f-sequences will be saved, create names of temporary outfiles for forward-seqs
##################################################################################
if ($save_f == 1) {
	#create f-outfilenames for one outfile per fastq_f_infile
	#that will hold destination plus seqentry
	for ($i = 0; $i < @fastq_f_infiles; ++$i) {
		$I_f_outfilenames[$i] = $i.'_f_out.fq';
	}
}

##################################################################################
#if r-sequences will be saved, create names of temporary outfiles for reverse-seqs
##################################################################################

if (($P == 1) and ($save_r == 1)) {
	#create r-outfilenames for one outfile per fastq_r_infile
	#that will hold destination plus seqentry
	for ($i = 0; $i < @fastq_r_infiles; ++$i) {
		$I_r_outfilenames[$i] = $i.'_r_out.fq';
	}
}

#################################################################
#create analysis and destination report filenames for each infile
#################################################################

#{
#for forward infiles
for ($i = 0; $i < @fastq_f_infiles; ++$i) {
	$I_f_ana_reportfilenames[$i] = $i.'_f_ana.txt';
	$I_f_dest_reportfilenames[$i] = $i.'_f_dest.txt';
}
#if paired infiles: for reverse infiles
if ($P == 1) {
	for ($i = 0; $i < @fastq_r_infiles; ++$i) {
		$I_r_ana_reportfilenames[$i] = $i.'_r_ana.txt';
		$I_r_dest_reportfilenames[$i] = $i.'_r_dest.txt';
	}
}
#}

##################################
#Main analysis: only forward reads
##################################

if ($P == 0) {#if sequences are not paired
	#loop through fastq-f-infiles
	my $pm1 = Parallel::ForkManager->new($max_p);
	for ($i = 0; $i < @fastq_f_infiles; ++$i) {
		$pm1->start and next;
		#open current forward-fastq-infile
		if ($fastq_f_infiles[$i] =~ /.gz$/) {#if infile has extension .gz open as gz file
			$f_fastqin_h = new IO::Uncompress::Gunzip $fastq_f_infiles[$i]
				or die "Cannot open $fastq_f_infiles[$i]: $GunzipError\n";
		} else {#open as text file
			$f_fastqin_h = IO::File->new("< $fastq_f_infiles[$i]")
				or die "Cannot open $fastq_f_infiles[$i]: $!\n";
		}
		#create outfile for destination+seq for this fastq f infile
		if ($save_f == 1) {
			$I_f_outfilename = $I_f_outfilenames[$i];
			$I_f_outfilehandle = IO::File->new(">$I_f_outfilename")
				or die "Could not open $I_f_outfilename: $!\n";
		}
		#message to screen
		print "analyzing $fastq_f_infiles[$i]\n";
		#go through current fastq-infile
		while ($fseqline1 = <$f_fastqin_h>) {
			#get the four seqlines of current forward-entry
			$fseqline2 = <$f_fastqin_h>;
			$fseqline3 = <$f_fastqin_h>;
			$fseqline4 = <$f_fastqin_h>;
			#remove line endings
			chomp $fseqline1;
			chomp $fseqline2;
			chomp $fseqline3;
			chomp $fseqline4;	
			#if active, change header
			if ($ch_h == 1) {
				ch_h(\$fseqline1,$f_read_No);
			}
			#call decision-tree to analyze and/or modify forward-seq-entry
			@f_tree_result = tree(\%f_tree_arguments);
			#store analysis-results in %f_report_counts
			$f_analysis_counts{$f_tree_result[2]} += 1;
			#determine f_destination
			$f_destination = $f_tree_result[1];
			#if adapter found in good sequence and adapter-containing sequences
			#will be saved in separate outfiles
			#append "_ad" to destination
			if (($f_tree_result[3] eq 'ad') and ($f_tree_result[1] ne 'NN') and ($adsf == 1)) {
				$f_destination .= '_ad';
			}
			#count destination
			$f_destination_counts{$f_destination} += 1;
			#if the sequence does not contain adapter
			#and if its destination is not NN
			#truncate it finally if selected by user
			if (($f_tree_result[3] ne 'ad') and ($f_tr > 0) and ($f_destination ne 'NN')) {
				trunc(\$fseqline2,\$fseqline4,\$f_tr);
			}			
			#if f-sequences will be saved, save to fastq-f-outfile for this fastq-f-infile
			if ($save_f == 1) {
				print $I_f_outfilehandle "$f_destination\n$fseqline1\n$fseqline2\n$fseqline3\n$fseqline4\n";
			}
			#go to next sequence
		}
		#print analysis- and destination counts for this infile to outfiles
		$I_f_ana_reportfilename = $I_f_ana_reportfilenames[$i];
		$f_outfilehandle = IO::File->new(">$I_f_ana_reportfilename")
			or die "Could not open $I_f_ana_reportfilename: $!\n";
		while (($key1, $value) = each %f_analysis_counts) {
			print $f_outfilehandle "$key1 $value ";
		}
		close $f_outfilehandle;
		$I_f_dest_reportfilename = $I_f_dest_reportfilenames[$i];
		$f_outfilehandle = IO::File->new(">$I_f_dest_reportfilename")
			or die "Could not open $I_f_dest_reportfilename: $!\n";
		while (($key1, $value) = each %f_destination_counts) {
			print $f_outfilehandle "$key1 $value ";
		}
		close $f_outfilehandle;
		#set counthashes back
		undef %f_analysis_counts;
		undef %f_destination_counts;
		#close current fastq-infile and seq-outfile
		close $f_fastqin_h;
		close $I_f_outfilehandle;
		$pm1->finish;
		#go to next fastq-infile
	}
	$pm1->wait_all_children;
}

############################
#Main analysis: paired reads
############################

#{
if ($P == 1) {
	#loop through fastq-infile pairs
	my $pm2 = Parallel::ForkManager->new($max_p);
	for ($i = 0; $i < @fastq_f_infiles; ++$i) {
		$pm2->start and next;
		#open current forward-fastq-infile
		if ($fastq_f_infiles[$i] =~ /.gz$/) {#if infile has extension .gz open as gz file
			$f_fastqin_h = new IO::Uncompress::Gunzip $fastq_f_infiles[$i]
				or die "Cannot open $fastq_f_infiles[$i]: $GunzipError\n";
		} else {#open as text file
			$f_fastqin_h = IO::File->new("< $fastq_f_infiles[$i]")
				or die "Cannot open $fastq_f_infiles[$i]: $!\n";
		}
		#open current reverse-fastq-infile
		if ($fastq_r_infiles[$i] =~ /.gz$/) {#if infile has extension .gz open as gz file
			$r_fastqin_h = new IO::Uncompress::Gunzip $fastq_r_infiles[$i]
				or die "Cannot open $fastq_r_infiles[$i]: $GunzipError\n";
		} else {#open as text file
			$r_fastqin_h = IO::File->new("< $fastq_r_infiles[$i]")
				or die "Cannot open $fastq_r_infiles[$i]: $!\n";
		}
		#create outfiles for destination+seq for this pair of fastq infiles
		if ($save_f == 1) {
			$I_f_outfilename = $I_f_outfilenames[$i];
			$I_f_outfilehandle = IO::File->new(">$I_f_outfilename")
				or die "Could not open $I_f_outfilename: $!\n";
		}
		if ($save_r == 1) {
			$I_r_outfilename = $I_r_outfilenames[$i];
			$I_r_outfilehandle = IO::File->new(">$I_r_outfilename")
				or die "Could not open $I_r_outfilename: $!\n";
		}
		#message to screen
		print "analyzing $fastq_f_infiles[$i] and $fastq_r_infiles[$i]\n";
		#go through current pair of fastq-infiles
		while ($fseqline1 = <$f_fastqin_h>) {
			#get the four seqlines of current forward-entry
			$fseqline2 = <$f_fastqin_h>;
			$fseqline3 = <$f_fastqin_h>;
			$fseqline4 = <$f_fastqin_h>;
			#remove line endings
			chomp $fseqline1;
			chomp $fseqline2;
			chomp $fseqline3;
			chomp $fseqline4;
			#get the four seqlines of current forward-entry
			$rseqline1 = <$r_fastqin_h>;
			$rseqline2 = <$r_fastqin_h>;
			$rseqline3 = <$r_fastqin_h>;
			$rseqline4 = <$r_fastqin_h>;
			#remove line endings
			chomp $rseqline1;
			chomp $rseqline2;
			chomp $rseqline3;
			chomp $rseqline4;
			#if active, change headers
			if ($ch_h == 1) {
				ch_h(\$fseqline1,$f_read_No);
				ch_h(\$rseqline1,$r_read_No);
			}
			#call decision tree to analyze and/or modify forward-seq-entry
			@f_tree_result = tree(\%f_tree_arguments);
			#in the following cases, don't analyze reverse-seq-entry:
			#f-seq contains internal restriction site
			if ($f_tree_result[3] eq 'rint') {
				@r_tree_result = ('0','NN','rint_pair','rint');
				$f_destination = 'NN';
				$r_destination = 'NN';
			}
			#f-seq contains adapter, insert too short
			elsif ($f_tree_result[3] eq 'ads') {
				@r_tree_result = ('0','NN','ads_pair','ads');
				$f_destination = 'NN';
				$r_destination = 'NN';
			}
			#f-seq good, contains adapter, rcdupa active
			elsif (($f_tree_result[1] ne 'NN') and ($f_tree_result[3] eq 'ad') and ($rcdupa == 1)){
				#call rcdupa with f-seqlines first because adapter has been found in f-seq
				rcdupa(\$fseqline2,\$fseqline4,\$rseqline2,\$rseqline4);
				$r_tree_result[0] = 1;#goods
				$r_tree_result[1] = $f_tree_result[1];#copy destination from f-seq
				$r_tree_result[2] = 'ad_pair';
				$r_tree_result[3] = 'ad';
				$f_destination = $f_tree_result[1];
				if ($adsf==1){#if adapter containing sequences will be saved in separate outfiles
					$f_destination .= '_ad';
				}
				$r_destination = $f_destination;
			}
			#None of these cases: analyze reverse-seq-entry
			else {
				#modify settings for analysis of reverse read
				#according to r-settings
				#f_destination ($f_tree_result[1]) and
				#f_override ($f_tree_result[3])
				#
				#If user wants to search for barcode adapter:
				if ($r_bcad_c == 1) {
					#if barcode of forward-seq is unknown:
					if (($f_tree_result[1] eq 'NN') or ($f_tree_result[1] eq 'SAMP')) {
						$r_ad_c = 0;#don't search (not possible)
					}
					#if barcode of forward-seq is known:
					else {
						#Specify all matching variants of barcode adapter corresponding
						#to barcode of forward seq and search
						@r_adseq = @{$bc_adapters{$f_tree_result[1]}};
						$r_ad_c = 1;
					}
				}
				#call decision tree to analyze and/or modify reverse-seq-entry
				@r_tree_result = tree(\%r_tree_arguments);
				#determine final destinations of forward- and reverse seq
				#r-seq contains internal restriction site
				if ($r_tree_result[3] eq 'rint') {
					$f_destination = 'NN';
					$r_destination = 'NN';
				}
				#r-seq contains adapter, insert too short
				elsif ($r_tree_result[3] eq 'ads') {
					$f_destination = 'NN';
					$r_destination = 'NN';
				}
				#r-seq is bad for other reasons
				elsif ($r_tree_result[1] eq 'NN') {
					$r_destination = 'NN';
					#f-seq is also bad
					if ($f_tree_result[1] eq 'NN') {
						$f_destination = 'NN';
					}
					#f-seq good and contains adapter (->rcdupa inactive, if active, we would not have analyzed r-seq)
					elsif ($f_tree_result[3] eq 'ad') {
						if ($adsf == 1) {#separate outfiles for adapter cont. seqs
							$f_destination = $f_tree_result[1].'_ad_sg';
						} else {#no separate outfiles for adapter cont. seqs
							$f_destination = $f_tree_result[1].'_sg';
						}
					}
					#f-seq good, no adapter
					else {
						$f_destination = $f_tree_result[1].'_sg';
					}
				}
				#r-seq is good, with adapter, rcdupa is active
				elsif (($r_tree_result[3] eq 'ad') and ($rcdupa == 1)) {
					#call rcdupa with r-seqlines first because adapter has been detected in r-seq
					rcdupa(\$rseqline2,\$rseqline4,\$fseqline2,\$fseqline4);
					#f-seq is bad
					if ($f_tree_result[1] eq 'NN') {
						$f_destination = 'SAMP';
						if ($adsf == 1) {#separate outfiles for adapter cont. seqs
							$f_destination .= '_ad';
						}
						$r_destination = $f_destination;
					}
					#f-seq good (no adapter, rcdupa is active, we would not have analyzed r-seq)
					else {
						$f_destination = $f_tree_result[1];
						if ($adsf == 1) {#separate outfiles for adapter cont. seqs
							$f_destination .= '_ad';
						}
						$r_destination = $f_destination;
					}
				}
				#r-seq is good, with adapter, rcdupa inactive
				elsif (($r_tree_result[3] eq 'ad') and ($rcdupa == 0)) {
					#f-seq is bad
					if ($f_tree_result[1] eq 'NN') {
						$f_destination = 'NN';
						if ($adsf == 1) {#separate outfiles for adapter cont. seqs
							$r_destination = 'SAMP_ad_sg';
						} else {#no separate outfiles for adapter cont. seqs
							$r_destination = 'SAMP_sg';
						}
					}
					#f-seq is good (with or without adapter)
					else {
						$f_destination = $f_tree_result[1];
						if ($adsf == 1) {#separate outfiles for adapter cont. seqs
							$f_destination .= '_ad';
						}
						$r_destination = $f_destination;
					}
				}
				#r-seq is good, no adapter
				elsif ($r_tree_result[1] ne 'NN') {
					#f-seq is bad
					if ($f_tree_result[1] eq 'NN') {
						$f_destination = 'NN';
						$r_destination = 'SAMP_sg';
					}
					#f-seq is good with adapter
					elsif ($f_tree_result[3] eq 'ad') {
						$f_destination = $f_tree_result[1];
						if ($adsf == 1) {#separate outfiles for adapter cont. seqs
							$f_destination .= '_ad';
						}
						$r_destination = $f_destination;
					}
					#f-seq is good, no adapter
					else {
						$f_destination = $f_tree_result[1];
						$r_destination = $f_destination;
					}
				}
			}
			#if none of the sequences contains adapter
			#and if sequence will not be stored in NN
			#truncate finally if selected by user
			if (($f_tree_result[3] ne 'ad') and ($r_tree_result[3] ne 'ad')) {
				if (($f_tr > 0) and ($f_destination ne 'NN')) {
					trunc(\$fseqline2,\$fseqline4,\$f_tr);
				}
				if (($r_tr > 0) and ($r_destination ne 'NN')) {
					trunc(\$rseqline2,\$rseqline4,\$r_tr);
				}
			}			
			#if active, call catsrc if none of the sequences contains adapter and none has destination NN
			if (($catsrc == 1) and ($f_tree_result[3] ne 'ad') and ($r_tree_result[3] ne 'ad') and
			($f_destination ne 'NN') and ($r_destination ne 'NN')) {
				catsrc(\$fseqline2,\$fseqline4,\$rseqline2,\$rseqline4);
			}			
			#count analysis-results and final destinations
			$f_analysis_counts{$f_tree_result[2]} += 1;
			$r_analysis_counts{$r_tree_result[2]} += 1;
			$f_destination_counts{$f_destination} += 1;
			$r_destination_counts{$r_destination} += 1;			
			#if active, save sequences to fastq-outfiles for these fastq-infiles
			if ($save_f == 1) {
				print $I_f_outfilehandle "$f_destination\n$fseqline1\n$fseqline2\n$fseqline3\n$fseqline4\n";
			}
			if ($save_r == 1) {
				print $I_r_outfilehandle "$r_destination\n$rseqline1\n$rseqline2\n$rseqline3\n$rseqline4\n";
			}
			#go to next pair of sequences			
		}
		#print analysis- and destination counts for these infiles to outfiles
		$I_f_ana_reportfilename = $I_f_ana_reportfilenames[$i];
		$f_outfilehandle = IO::File->new(">$I_f_ana_reportfilename")
			or die "Could not open $I_f_ana_reportfilename: $!\n";
		while (($key1, $value) = each %f_analysis_counts) {
			print $f_outfilehandle "$key1 $value ";
		}
		close $f_outfilehandle;
		$I_f_dest_reportfilename = $I_f_dest_reportfilenames[$i];
		$f_outfilehandle = IO::File->new(">$I_f_dest_reportfilename")
			or die "Could not open $I_f_dest_reportfilename: $!\n";
		while (($key1, $value) = each %f_destination_counts) {
			print $f_outfilehandle "$key1 $value ";
		}
		close $f_outfilehandle;
		$I_r_ana_reportfilename = $I_r_ana_reportfilenames[$i];
		$r_outfilehandle = IO::File->new(">$I_r_ana_reportfilename")
			or die "Could not open $I_r_ana_reportfilename: $!\n";
		while (($key1, $value) = each %r_analysis_counts) {
			print $r_outfilehandle "$key1 $value ";
		}
		close $r_outfilehandle;
		$I_r_dest_reportfilename = $I_r_dest_reportfilenames[$i];
		$r_outfilehandle = IO::File->new(">$I_r_dest_reportfilename")
			or die "Could not open $I_r_dest_reportfilename: $!\n";
		while (($key1, $value) = each %r_destination_counts) {
			print $r_outfilehandle "$key1 $value ";
		}
		close $r_outfilehandle;
		#set counthashes back
		undef %f_analysis_counts;
		undef %f_destination_counts;
		undef %r_analysis_counts;
		undef %r_destination_counts;
		#close current pair of fastq-infiles and seq outfiles
		close $f_fastqin_h;
		close $r_fastqin_h;
		close $I_f_outfilehandle;
		close $I_r_outfilehandle;
		$pm2->finish;
		#go to next pair of fastq-infiles
	}
	$pm2->wait_all_children;
}
#}

#Determine runtime up to now and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Analysis phase took $run_s seconds.\n";

###################################
#Determine fastq outfile extensions
###################################

#{
if ($z == 1) {#save fastq outfiles as gzip
	$f_outfilesuffix = '_R1.fq.gz';
	$r_outfilesuffix = '_R2.fq.gz';
	$fr_outfilesuffix = '_R12.fq.gz';
	
} else {#save as text
	$f_outfilesuffix = '_R1.fq';
	$r_outfilesuffix = '_R2.fq';
	$fr_outfilesuffix = '_R12.fq';
}
#}

###################################################################
#if f-sequences will be saved, open final fastq outfiles for f-seqs
#if f- and r-sequences will be saved in combined outfiles: open
###################################################################

#{
if ($save_f == 1) {#if foward sequences will be saved
	#create destinations
	if ($save_d == 1) {#if discarded sequences will be saved
		push @f_destinations, 'NN';
	}
	if ($f_bc == 0) {#no barcodes
		push @f_destinations, 'SAMP';
		if ($P == 1) {#paired sequences
			push @f_destinations, 'SAMP_sg';
		}
		if ($adsf == 1) {#adapter containing sequences in separate outfiles
			if (($f_ad_c == 1) or (($P == 1) and ($r_sad_c == 1))) {
				push @f_destinations, 'SAMP_ad';
			}
			if (($P == 1) and ($f_ad_c == 1) and ($rcdupa == 0)) {
				push @f_destinations, 'SAMP_ad_sg';
			}
		}
	}
	else {#barcodes
		for ($i = 0; $i < $n_barcodes; ++$i) {
			push @f_destinations, $barcodes[$i];
		}
		if ($P == 1) {#paired sequences
			for ($i = 0; $i < $n_barcodes; ++$i) {
				push @f_destinations, $barcodes[$i].'_sg';
			}
		}
		if ($adsf == 1) {#adapter containing sequences in separate outfiles
			if (($f_ad_c == 1) or (($P == 1) and (($r_sad_c == 1) or ($r_bcad_c == 1)))) {
				for ($i = 0; $i < $n_barcodes; ++$i) {
					push @f_destinations, $barcodes[$i].'_ad';
				}
			}
			if (($P == 1) and ($f_ad_c == 1) and ($rcdupa == 0)) {
				for ($i = 0; $i < $n_barcodes; ++$i) {
					push @f_destinations, $barcodes[$i].'_ad_sg';
				}
			}
		}
	}
	#if there are no reverse reads
	#or if only forward reads shall be saved
	#or if forward- and reverse-reads shall be saved in separate outfiles
	if (($P == 0) or (($P == 1) and ($save_r == 0)) or
	(($P == 1) and ($r1r2sf == 1) and ($save_f == 1) and ($save_r == 1))) {
		#create f-outfilenames for the concatenated outfiles
		for ($i = 0; $i < @f_destinations; ++$i) {
			$f_outfilenames[$i] = $f_destinations[$i].$f_outfilesuffix;
		}
		#create f-outfilehandles for the concatenated outfiles, open f-outfiles,
		#store f-outfilehandles in hash with destination as key
		for ($i = 0; $i < @f_destinations; ++$i) {
			$outfilename = "$outfolder/$f_outfilenames[$i]";
			if ($z == 1) {#save as gzip
				$f_outfilehandles{$f_destinations[$i]} = new IO::Compress::Gzip $outfilename ,-Level=>Z_BEST_SPEED
					or die "Could not write Gzip: $GzipError\n";
			} else {#save as text
				$f_outfilehandles{$f_destinations[$i]} = IO::File->new(">$outfilename")
					or die "Could not open $outfilename: $!\n";
			}
		}
	}
	#if there are also reverse reads and if forward- and reverse-reads
	#shall be saved in the same outfile
	elsif (($P == 1) and ($r1r2sf == 0) and ($save_f == 1) and ($save_r == 1)) {
		#create f-outfilenames for single f-outfiles (_NN and _sg) and
		#fr-outfilenames for combined outfiles (_SAMP_ and _BC_)
		for ($i = 0; $i < @f_destinations; ++$i) {
			if (($f_destinations[$i] eq 'NN') or ($f_destinations[$i] =~ /_sg$/)) {
				push @f_outfilenames, $f_destinations[$i].$f_outfilesuffix;
			} else {
				push @fr_outfilenames, $f_destinations[$i].$fr_outfilesuffix;
				push @fr_destinations, $f_destinations[$i];
			}
		}
		#create handles for combined outfiles and open
		#store handles in %f_outfilehandles with destination as key
		for ($i = 0; $i < @fr_destinations; ++$i) {
			$outfilename = "$outfolder/$fr_outfilenames[$i]";
			if ($z == 1) {#save as gzip
				$f_outfilehandles{$fr_destinations[$i]} = new IO::Compress::Gzip $outfilename ,-Level=>Z_BEST_SPEED
					or die "Could not write Gzip: $GzipError\n";
			} else {#save as text
				$f_outfilehandles{$fr_destinations[$i]} = IO::File->new(">$outfilename")
					or die "Could not open $outfilename: $!\n";
			}
		}
		#copy handles for combined outfiles into %r_outfilehandles
		%r_outfilehandles = %f_outfilehandles;
		#create handles for single f-outfiles (_NN and _sg) and open
		for ($i = 0; $i < @f_destinations; ++$i) {
			if (($f_destinations[$i] eq 'NN') or ($f_destinations[$i] =~ /_sg$/)) {
				$outfilename = "$outfolder/$f_destinations[$i]$f_outfilesuffix";
				if ($z == 1) {#save as gzip
					$f_outfilehandles{$f_destinations[$i]} = new IO::Compress::Gzip $outfilename ,-Level=>Z_BEST_SPEED
						or die "Could not write Gzip: $GzipError\n";
				} else {#save as text
					$f_outfilehandles{$f_destinations[$i]} = IO::File->new(">$outfilename")
						or die "Could not open $outfilename: $!\n";
				}
			}
		}
	}
}
#}

###################################################################
#if r-sequences will be saved, open final fastq outfiles for r-seqs
###################################################################

#{
if (($P == 1) and ($save_r == 1)) {#paired reads and r-sequences will be saved
	#create destinations
	if ($save_d == 1) {#save discarded sequences
		push @r_destinations, 'NN';
	}
	push @r_destinations, 'SAMP_sg';
	if ($f_bc == 0) {#no barcodes in f-sequences
		push @r_destinations, 'SAMP';
		if ($adsf == 1) {#adapter containing sequences in separate outfiles
			if (($f_ad_c == 1) or ($r_sad_c == 1)) {
				push @r_destinations, 'SAMP_ad';
			}
			if (($r_sad_c == 1) and ($rcdupa == 0)) {
				push @r_destinations, 'SAMP_ad_sg';
			}
		}
	}
	if ($f_bc == 1) {#barcodes in f-sequences
		for ($i = 0; $i < $n_barcodes; ++$i) {
			push @r_destinations, $barcodes[$i];
		}
		if ($adsf == 1) {#adapter containing sequences in separate outfiles
			if (($f_ad_c == 1) or ($r_sad_c == 1) or ($r_bcad_c == 1)) {
				for ($i = 0; $i < $n_barcodes; ++$i) {
					push @r_destinations, $barcodes[$i].'_ad';
				}
			}
			if (($r_sad_c == 1) and ($rcdupa == 0)) {
				push @r_destinations, 'SAMP_ad_sg';
			}
		}
	}	
	#if forward- and reverse reads shall be saved to separate outfiles
	#or if only reverse reads shall be saved
	if ((($save_f == 1) and ($r1r2sf == 1)) or ($save_f == 0)) {
		#create r-outfilenames for concatenated r-outfiles
		for ($i = 0; $i < @r_destinations; ++$i) {
			$r_outfilenames[$i] = $r_destinations[$i].$r_outfilesuffix;
		}
		#create r-outfilehandles, open r-outfiles,
		#store r-outfilehandles in hash with destination as key
		for ($i = 0; $i < @r_destinations; ++$i) {
			$outfilename = "$outfolder/$r_outfilenames[$i]";
			if ($z == 1) {#save as gzip
				$r_outfilehandles{$r_destinations[$i]} = new IO::Compress::Gzip $outfilename ,-Level=>Z_BEST_SPEED
					or die "Could not write Gzip: $GzipError\n";
			} else {#save as text
				$r_outfilehandles{$r_destinations[$i]} = IO::File->new(">$outfilename")
					or die "Could not open $outfilename: $!\n";
			}
		}
	}
	#if forward- and reverse-reads (if both not discarded) shall be combined
	#in the same outfile
	elsif(($save_f == 1) and ($r1r2sf == 0)) {
		#create r-outfilenames for single outfiles (NN and _sg)
		#create handles for single r-outfiles, open and store in %r_outfilehandles
		for ($i = 0; $i < @r_destinations; ++$i) {
			if (($r_destinations[$i] eq 'NN') or ($r_destinations[$i] =~ /_sg$/)) {
				push @r_outfilenames, $r_destinations[$i].$r_outfilesuffix;
				$outfilename = "$outfolder/$r_destinations[$i]$r_outfilesuffix";
				if ($z == 1) {#save as gzip
					$r_outfilehandles{$r_destinations[$i]} = new IO::Compress::Gzip $outfilename ,-Level=>Z_BEST_SPEED
						or die "Could not write Gzip: $GzipError\n";
				} else {#save as text
					$r_outfilehandles{$r_destinations[$i]} = IO::File->new(">$outfilename")
						or die "Could not open $outfilename: $!\n";
				}
			}
		}
	}
}
#}

###############################################################
#save names of all final fastq-outfiles to file "fastq_out.txt"
###############################################################

#{
	if (($save_f == 1) or ($save_r == 1)) {#if sequences will be saved at all
		#open the file
		unless(open(OUTFILE,">$outfolder/fastq_out.txt")) {
			print "Cannot create file fastq_out.txt. Exiting..\n\n";
			exit;
		}
		if ($save_f == 1) {#if forward reads will be saved
			#print names of fastq-outfiles from forward reads
			for ($i = 0; $i < @f_outfilenames; ++$i) {
				$tempstring = $f_outfilenames[$i];
				print OUTFILE "$tempstring\n";
			}
		}
		if (($P == 1) and ($save_r == 1)) {#if there are paired reads and r reads will be saved
			for ($i = 0; $i < @r_outfilenames; ++$i) {
				$tempstring = $r_outfilenames[$i];
				print OUTFILE "$tempstring\n";				
			}
			#if combined outfiles for forward- and reverse-reads will be saved
			if (($save_f == 1) and ($r1r2sf == 0)) {
				for ($i = 0; $i < @fr_outfilenames; ++$i) {
					$tempstring = $fr_outfilenames[$i];
					print OUTFILE "$tempstring\n";
				}
			}
		}
		close OUTFILE;
	}
#}

#############################################################
#Demultiplex and delete sequence-outfiles from single infiles
#############################################################

#{
#if sequences have been saved
if (($save_f == 1) or ($save_r == 1)) {
	#demultiplex
	#if only forward reads have been saved
	if ((($P == 0) or ($save_r == 0)) and ($save_f == 1)) {
		#forward: loop through sequence outfiles from single infiles
		for ($i = 0; $i < @I_f_outfilenames; ++$i) {
			#open current outfile to read in
			unless (open(INFILE, $I_f_outfilenames[$i])) {
				print "Cannot open $I_f_outfilenames[$i]. Exiting..\n";
				exit;
			}
			#loop through the file
			while ($f_destination = <INFILE>) {
				#get the sequence entry
				$fseqline1 = <INFILE>;
				$fseqline2 = <INFILE>;
				$fseqline3 = <INFILE>;
				$fseqline4 = <INFILE>;
				chomp $f_destination;
				#determine outfilehandle and print sequence to appropriate outfile
				if (defined $f_outfilehandles{$f_destination}) {
					$f_outfilehandle = $f_outfilehandles{$f_destination};
					print $f_outfilehandle "$fseqline1$fseqline2$fseqline3$fseqline4";
				}
			}
			close INFILE;
			unlink $I_f_outfilenames[$i];#delete the current sequence-outfile from single infile
		}
	}
	#if only reverse reads have been saved
	if (($P == 1) and ($save_f == 0) and ($save_r == 1)) {
		#reverse: loop through sequence outfiles from single infiles
		for ($i = 0; $i < @I_r_outfilenames; ++$i) {
			#open current outfile to read in
			unless (open(INFILE, $I_r_outfilenames[$i])) {
				print "Cannot open $I_r_outfilenames[$i]. Exiting..\n";
				exit;
			}
			#loop through the file
			while ($r_destination = <INFILE>) {
				#get the sequence entry
				$rseqline1 = <INFILE>;
				$rseqline2 = <INFILE>;
				$rseqline3 = <INFILE>;
				$rseqline4 = <INFILE>;
				chomp $r_destination;
				#determine outfilehandle and print sequence to appropriate outfile
				if (defined $r_outfilehandles{$r_destination}) {
					$r_outfilehandle = $r_outfilehandles{$r_destination};
					print $r_outfilehandle "$rseqline1$rseqline2$rseqline3$rseqline4";
				}
			}
			close INFILE;
			unlink $I_r_outfilenames[$i];#delete the current sequence-outfile from single infile
		}
	}
	#if there are paired reads and f- and r-reads have been saved:
	#Treat pairs of sequence outfiles from pairs of single infiles in parallel
	#This enables the following:
	#If forward- and reverse reads are combined (r1r2sf==0), forward-read and corresponding reverse-read 
	#appear consecutively (forward then reverse) in the combined outfile
	elsif (($P == 1) and ($save_f == 1) and ($save_r == 1)) {
		#loop through pairs of sequence outfiles from pairs of single infiles
		for ($i = 0; $i < @I_f_outfilenames; ++$i) {
			#open current pair of outfiles to read in
			unless (open(FINFILE, $I_f_outfilenames[$i])) {
				print "Cannot open $I_f_outfilenames[$i]. Exiting..\n";
				exit;
			}
			unless (open(RINFILE, $I_r_outfilenames[$i])) {
				print "Cannot open $I_r_outfilenames[$i]. Exiting..\n";
				exit;
			}
			#loop through the pair of files
			while ($f_destination = <FINFILE>) {#get f-destination
				#get f-sequence entry
				$fseqline1 = <FINFILE>;
				$fseqline2 = <FINFILE>;
				$fseqline3 = <FINFILE>;
				$fseqline4 = <FINFILE>;
				#get the r-destination and r-sequence entry
				$r_destination = <RINFILE>;
				$rseqline1 = <RINFILE>;
				$rseqline2 = <RINFILE>;
				$rseqline3 = <RINFILE>;
				$rseqline4 = <RINFILE>;
				chomp $f_destination;
				chomp $r_destination;
				if ((defined $f_outfilehandles{$f_destination}) and (defined $r_outfilehandles{$r_destination})) {
					#determine f-outfilehandle and print sequence to appropriate outfile
					$f_outfilehandle = $f_outfilehandles{$f_destination};
					print $f_outfilehandle "$fseqline1$fseqline2$fseqline3$fseqline4";
					#determine r-outfilehandle and print sequence to appropriate outfile
					$r_outfilehandle = $r_outfilehandles{$r_destination};
					print $r_outfilehandle "$rseqline1$rseqline2$rseqline3$rseqline4";
				}
			}
			close FINFILE;
			close RINFILE;
			unlink $I_f_outfilenames[$i];#delete the current f-sequence-outfile from single infile
			unlink $I_r_outfilenames[$i];#delete the current r-sequence-outfile from single infile
		}
	}
}
#}

######################################################################
#summarize the analysis counts from the separate files for each infile
######################################################################

#{
#loop through the files from forward reads
for ($i = 0; $i < @fastq_f_infiles; ++$i) {
	unless(open(INFILE, $I_f_ana_reportfilenames[$i])){
		print "cannot open file $I_f_ana_reportfilenames[$i].Exiting..\n\n";
		exit;
	}
	$tempstring = <INFILE>;
	close INFILE;
	chomp $tempstring;
	$tempstring =~ s/\s*$//;#delete whitespace at end
	@one_line_arr = split ('\s+',$tempstring);#explode into array
	#loop through the array and take every second element (as key..)
	for ($j = 0; $j < (@one_line_arr - 1); $j += 2) {
		$f_analysis_counts{$one_line_arr[$j]} += $one_line_arr[$j+1];#next element as value
	}
}
#loop through the files from reverse reads if existing
if ($P == 1) {
	for ($i = 0; $i < @fastq_r_infiles; ++$i) {
		unless(open(INFILE, $I_r_ana_reportfilenames[$i])){
			print "cannot open file $I_r_ana_reportfilenames[$i].Exiting..\n\n";
			exit;
		}
		$tempstring = <INFILE>;
		close INFILE;
		chomp $tempstring;
		$tempstring =~ s/\s*$//;#delete whitespace at end
		@one_line_arr = split ('\s+',$tempstring);#explode into array
		#loop through the array and take every second element (as key..)
		for ($j = 0; $j < (@one_line_arr - 1); $j += 2) {
			$r_analysis_counts{$one_line_arr[$j]} += $one_line_arr[$j+1];#next element as value
		}
	}
}
#}

#########################################################################
#summarize the destination counts from the separate files for each infile
#########################################################################

#{
#loop through the files from forward reads
for ($i = 0; $i < @fastq_f_infiles; ++$i) {
	unless(open(INFILE, $I_f_dest_reportfilenames[$i])){
		print "cannot open file $I_f_dest_reportfilenames[$i].Exiting..\n\n";
		exit;
	}
	$tempstring = <INFILE>;
	close INFILE;
	chomp $tempstring;
	$tempstring =~ s/\s*$//;#delete whitespace at end
	@one_line_arr = split ('\s+',$tempstring);#explode into array
	#loop through the array and take every second element (as key..)
	for ($j = 0; $j < (@one_line_arr - 1); $j += 2) {
		$f_destination_counts{$one_line_arr[$j]} += $one_line_arr[$j+1];#next element as value
	}
}
#loop through the files from reverse reads if existing
if ($P == 1) {
	for ($i = 0; $i < @fastq_r_infiles; ++$i) {
		unless(open(INFILE, $I_r_dest_reportfilenames[$i])){
			print "cannot open file $I_r_dest_reportfilenames[$i].Exiting..\n\n";
			exit;
		}
		$tempstring = <INFILE>;
		close INFILE;
		chomp $tempstring;
		$tempstring =~ s/\s*$//;#delete whitespace at end
		@one_line_arr = split ('\s+',$tempstring);#explode into array
		#loop through the array and take every second element (as key..)
		for ($j = 0; $j < (@one_line_arr - 1); $j += 2) {
			$r_destination_counts{$one_line_arr[$j]} += $one_line_arr[$j+1];#next element as value
		}
	}
}
#}

#####################
#print report to file
#####################

#{
#open report-outfile for summary report over all infiles
unless (open(REPORTFILE,">$outfolder/$report_filename")) {
	print "Cannot open file \"$report_filename\", exiting ...\n\n";
	exit;
}
#get array of forward analysis results
@f_analysis_results = sort keys(%f_analysis_counts);
#print analysis results to reportfile
print REPORTFILE "analysis results forward seqs:\n";
for ($i = 0; $i < @f_analysis_results; ++$i) {
	print REPORTFILE "$f_analysis_results[$i]\t$f_analysis_counts{$f_analysis_results[$i]}\n";
}
#when infiles are paired
if ($P == 1) {
	#get array of reverse analysis results
	@r_analysis_results = sort keys(%r_analysis_counts);
	#print analysis results to reportfile
	print REPORTFILE "analysis results reverse seqs:\n";
	for ($i = 0; $i < @r_analysis_results; ++$i) {
		print REPORTFILE "$r_analysis_results[$i]\t$r_analysis_counts{$r_analysis_results[$i]}\n";
	}
}

#get array of destinations forward-seqs have been sent to
@f_destination_results = sort keys(%f_destination_counts);
#print f-destination counts to reportfile
print REPORTFILE "foward-seqs destinations:\n";
for ($i = 0; $i < @f_destination_results; ++$i) {
	print REPORTFILE "$f_destination_results[$i]\t$f_destination_counts{$f_destination_results[$i]}\n";
}
#when infiles are paired
if ($P == 1) {
	#get array of destinations reverse-seqs have been sent to
	@r_destination_results = sort keys(%r_destination_counts);
	#print r-destination counts to reportfile
	print REPORTFILE "reverse-seqs destinations:\n";
	for ($i = 0; $i < @r_destination_results; ++$i) {
		print REPORTFILE "$r_destination_results[$i]\t$r_destination_counts{$r_destination_results[$i]}\n";
	}
}
#}

##################################################
#Delete all remaining outfiles from single infiles
##################################################

#{
unlink @I_f_ana_reportfilenames;
unlink @I_f_dest_reportfilenames;
unlink @I_r_ana_reportfilenames;
unlink @I_r_dest_reportfilenames;
#}


#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "run took $run_s seconds.\n";

exit;

########################################################
#subroutines
########################################################

#definition of subroutine ambi_var
#expects: a scalar: DNA sequence
#only capitals, only IUPAC symbols ACGTRYSWKMBDHVN
#returns: array of DNA sequences
#if the sequence does not contain ambiguities:
#one entry, identical with the input sequence
#else: all possible matching variants of the sequence
sub ambi_var {
	#declare and initialize
	my ($inseq) = @_;
	my $inlen = 0;#length of input sequence
	my $inchar = '';#one character of input sequence
	my $varchar = '';#one matching character for an input-ambiguity character
	my $pos = 0;#current position in input-sequence
	my @variants = ();#all possible matching variants of $inseq
	my %IUPAC = (
	R => ['A','G'],
	Y => ['C','T'],
	S => ['C','G'],
	W => ['A','T'],
	K => ['G','T'],
	M => ['A','C'],
	B => ['C','G','T'],
	D => ['A','G','T'],
	H => ['A','C','T'],
	V => ['A','C','G'],
	N => ['A','C','G','T']
	);
	my $currseq = '';#current sequence that will be permuted to variants
	my $newseq = '';#new permuted variant
	my @temparr1 = ();
	
	$inlen = length $inseq;#determine length of input sequence
	$variants[0] = $inseq;#initialize @variants with $inseq
	#loop through positions of input sequence
	for ($pos = 0; $pos < $inlen; ++$pos) {
		$inchar = substr($inseq,$pos,1);#get current character
		if (defined $IUPAC{$inchar}){#if it is an ambiguity
			#loop through sequences in @variants
			for $currseq (@variants) {
				#loop through characters that match current ambiguity
				for $varchar (@{$IUPAC{$inchar}}) {
					$newseq = $currseq;#copy current sequence in @variants into $newseq
					substr($newseq,$pos,1) = $varchar;#change current position to current matching character
					push @temparr1, $newseq;#add new sequence to temporary array
				}
			}
			@variants = @temparr1;#replace old variants with (more numerous) new variants
			@temparr1 = ();#set back and go to next pos in input sequence
		}
	}
	return @variants;
}

#definition of subroutine tree
#expects a hash with references to many variables in main
#calls other subs to do work
#returns: an array with 3 entries:
#0:goods
#1:dest
#2:ana
sub tree {
	#declare and initialize arguments and variables
	my ($argumentsref) = @_;#reference to the hash of references
	my (%arg) = %$argumentsref;
	my (@results) = ();
	my ($goods) = 1;
	my ($dest) = 0;
	my ($ana) = 'clean';
	my ($override) = 'none';
	my ($adstart) = 0;#adapter start position, count starting with 0, -1 if no adapter found
	my ($ila) = 0;#insert length of adapter containing sequences
	my ($len) = 0;#a sequence length
	
	#initialize @results
	$results[0] = $goods;
	$results[1] = $dest;
	$results[2] = $ana;
	$results[3] = $override;
		
	#truncate sequence first if active
	if (${$arg{"trf"}} > 0) {
		trunc($arg{"seqline2"},$arg{"seqline4"},$arg{"trf"});
	}
	
	#do quality check if active
	if (${$arg{"q"}} == 1) {
		#yes: call subroutine sliding_w
		$goods = sliding_w($arg{"seqline4"},$arg{"wl"},$arg{"pt"},$arg{"po"});
		#if sequence is bad, return bad result
		if ($goods == 0) {
			$results[0] = 0;
			$results[1] = 'NN';
			$results[2] = 'q_fail';
			return @results;
		}
	}
	#if barcode not present: enter no-bc-chain
	if (${$arg{"bc"}} == 0) {
		$dest = 'SAMP';
		#if active, check restriction-site1
		if (${$arg{"r1_c"}} == 1) {
			$goods = r1_c($arg{"seqline2"},$arg{"restr_site1"});
			#if sequence bad:
			if ($goods == 0) {
				#if active, try to repair restriction site1
				if (${$arg{"r1_rep"}} == 1) {
					$goods = r1_rep($arg{"seqline2"},$arg{"restr_site1"},$arg{"r_reppos"});
					#if sequence repaired, store analysis result
					if ($goods == 1) {
						$ana = 'r1_rep';
					}
					#if sequence still bad (not repaired) return bad result
					else {
						$results[0] = 0;
						$results[1] = 'NN';
						$results[2] = 'r1_rep_fail';
						return @results;
					}
				}
				#if repair not active, return bad result
				else {
					$results[0] = 0;
					$results[1] = 'NN';
					$results[2] = 'r1_c_fail';
					return @results;
				}
			}
		}
	}
	#if barcode present, enter bc-chain
	else {
		#check barcode and restriction site
		$dest = bcr_c($arg{"seqline2"},$arg{"barcodes_restr"});
		#if check failed (dest==NN), enter repair-chain
		if ($dest eq 'NN') {
			#record bad analysis result
			$ana = 'bcr_c_fail';
			#if active try to repair restriction site
			if (${$arg{"bc_rrep"}} == 1) {
				$dest = bc_rrep($arg{"seqline2"},$arg{"barcodes_restr"},$arg{"r_reppos"});
				#record result
				if ($dest ne 'NN') {
					$ana .= '_bc_rrep';#repaired
				}
				else {
					$ana .= '_bc_rrep_fail';#repair failed
				}
			}
			#if barcode still unclear, if active, try to repair barcode
			if ($dest eq 'NN' and ${$arg{"bcrep_r"}} == 1) {
				$dest = bcrep_r($arg{"seqline2"},$arg{"barcodes"},
								$arg{"restr_site1"},$arg{"bc_reppos"});
				#record result
				if ($dest ne 'NN') {
					$ana .= '_bcrep_r';#repaired
				}
				else {
					$ana .= '_bcrep_r_fail';#repair failed
				}
			}
			#if barcode still unclear, if active, try to repair barcode and restr.site
			if ($dest eq 'NN' and ${$arg{"bcrep_rrep"}} == 1) {
				$dest = bcrep_rrep($arg{"seqline2"},$arg{"barcodes"},$arg{"restr_site1"},
									$arg{"bc_reppos"},$arg{"r_reppos"});
				#record result
				if ($dest ne 'NN') {
					$ana .= '_bcrep_rrep';#repaired
				}
				else {
					$ana .= '_bcrep_rrep_fail';#repair failed
				}
			}
			#if barcode still unclear, return bad result
			if ($dest eq 'NN') {
				$results[0] = 0;
				$results[1] = 'NN';
				$results[2] = $ana;
				return @results;
			}
		}
		#if active, clip barcode
		if (${$arg{"cl_bc"}} == 1) {
			cl_bc($arg{"seqline2"},$arg{"seqline4"},$dest);
		}
		#if active, compensate for different barcode lengths by
		#truncation of sequences with short barcodes
		if (${$arg{"compbcl"}} == 1) {
			#determine new length
			$len = (length ${$arg{"seqline2"}}) + (length $dest) - ${$arg{"lonbc"}};
			#truncate
			trunc($arg{"seqline2"},$arg{"seqline4"},\$len);
		}
	}
	#if active, check for internal restriction site
	if (${$arg{"rint_c"}} == 1) {
		$goods = rint_c($arg{"seqline2"},$arg{"restr_site_int"});
		#if sequence is bad, stop and return results
		if ($goods == 0) {
			$ana .= '_rint';
			$dest = 'NN';
			$override = 'rint';
			$results[0] = 0;
			$results[1] = $dest;
			$results[2] = $ana;
			$results[3] = $override;
			return @results;
		}
	}
	#if ad_c is active, enter adapter chain
	if (${$arg{"ad_c"}} == 1) {
		$adstart = ad_c_pos($arg{"seqline2"},$arg{"adseq"},$arg{"admis"});
		#if sequence contains adapter, store info in $ana and $override
		if ($adstart != -1) {
			$ana .= '_ad';
			$override = 'ad';
			#if cropia is active, continue
			if (${$arg{"cropia"}} == 1) {
				#determine cropped insert length
				$ila = $adstart + ${$arg{"ao"}} - ${$arg{"restr_site1_l"}};
				#if it is too short, return bad result
				if ($ila < ${$arg{"ilamin"}}) {
					$ana .= 's';
					$results[0] = 0;
					$results[1] = 'NN';
					$results[2] = $ana;
					$results[3] = 'ads';
					return @results;
				}
				#still alive: determine truncation length for adapter removal
				$ila += ${$arg{"restr_site1_l"}};
				#remove adapter
				trunc($arg{"seqline2"},$arg{"seqline4"},\$ila);
				#remove leading restriction site
				$goods = cl_r($arg{"seqline2"},$arg{"seqline4"},$arg{"restr_site1"});
				#return
				$results[0] = $goods;
				$results[1] = $dest;
				$results[2] = $ana;
				$results[3] = $override;
				return @results;
			}
		}
	}
	#if active, remove leading restriction site
	if (${$arg{"cl_r"}} == 1) {
		$goods = cl_r($arg{"seqline2"},$arg{"seqline4"},$arg{"restr_site1"});
		#if clipping failed ($goods == 0), stop and return results
		if ($goods == 0) {
			$ana .= '_cl_r_fail_seq_in_NN';
			$dest = 'NN';
			$results[0] = 0;
			$results[1] = $dest;
			$results[2] = $ana;
			return @results;
		}
	}
	#return results
	$results[0] = $goods;
	$results[1] = $dest;
	$results[2] = $ana;
	$results[3] = $override;
	return @results;
}

#definition of subroutine ch_h ver. 2
#expects
#$headerref: reference to a string that holds header-line of a fastq-sequence-entry
#$read_No: number that will appear at end of new header (1 for f read, 2 for r read)
#changes the header-line to simpler format
#example:@4_2312_4042_29567_1
#explanation: @lane_sideswathtile_coordinate_coordinate_readNo
sub ch_h {
	#declare and initialize arguments and variables
	my ($headerref,$readNo) = @_;
	my $headerline = '';#my copy of the headerline
	my @header_arr = ();#the headerline in an array
	my @header_arr2 = ();#some needed fields of the headerline in a 2nd array
	my $i = 0;#counter

	$headerline = $$headerref;
	chomp $headerline;
	$headerline =~ s/^\s*//;
	$headerline =~ s/\s*$//;
	$headerline =~ s/\s+/:/g;
	@header_arr = split(':',$headerline);
	#copy the needed fields into second array
	#(lane, tile, coord, coord)
	for ($i = 3; $i < 7; ++$i) {
		push @header_arr2, $header_arr[$i];
	}
	push @header_arr2, $readNo;#add readNo
	$headerline = join('_',@header_arr2);
	substr($headerline, 0, 0, "@");#add "@ at the begin of $headerline
	$$headerref = $headerline;
}

#definition of subroutine sliding_w
#expects as arguments:
#reference to a scalar: string, sequence of Phred symbols
#reference to a scalar: window length (integer)
#reference to a scalar: Phred score threshold (integer)
#reference to a scalar: Phred score offset (integer) (e.g. 33 for Sanger)
#analyzes the Phred symbol line with a sliding window
#of length window_length
#step-length: 1
#starts from the end of the sequence and moves to the beginning
#in each window: calculates average Phred score
#if below Phred score threshold in any window:
#returns "0" (bad), else: returns: "1"(good)
sub sliding_w {
	#initialize arguments and variables
	my($seqref, $wlref, $ptref, $offsetref) = @_;
	my($seq) = $$seqref;#string: sequence of Phred symbols
	my($wl) = $$wlref;#window length
	my($pt) = $$ptref;#Phred score threshold
	my($offset) = $$offsetref;#Phred score offset
	my(@seqarr) = ();#array: Phred symbol sequence
	my($seql) = 0;#sequence length
	my($i) = 0;#counter
	my($ws) = 0;#window-sum: sum of Phred scores in window
	my($result) = 1;
	
	#add offset to Phred score threshold
	#Phred scores are ASCII-order-numbers minus offset
	#it is faster to add the offset to the threshold once
	#than to substract it from every Phred score
	$pt += $offset;
	
	#multiply $pt with $wl
	#saves calculating the window average
	#we can just compare the window sum to $pt
	$pt = $pt * $wl;

	chomp $seq;
	@seqarr = split('',$seq);
	$seql = @seqarr;
	#calculate window sum for first window (at end of sequence)
	for ($i = $seql - 1 ; $i >= $seql - $wl; --$i) {
		$ws += ord($seqarr[$i]);
	}
		
	#as long as window-sum >= threshold and sequence-beginning not reached,
	#move window through sequence (from end to beginning) and calculate window average
	$i = $seql - $wl - 1;
	while (($i >= 0) and ($ws >= $pt)) {
		$ws += ord($seqarr[$i]);#add next Phred score to windowsum
		$ws -= ord($seqarr[$i+$wl]);#subtract last phredscore from windowsum
		--$i;#decrement counter
	}
	if ($ws < $pt) {#if window-sum below threshold,
		$result = 0;#change result to 0
	}
	return $result;
}

#definition of subroutine r1_c ver.2
#expects references to
#$seq			DNA sequence
#@restr_site1	all sequence variants that match restriction site at beginning of DNA sequence
#				specified by the user (only 1 if restriction site does not contain ambiguities)
#Checks if DNA sequence begins with one of the restriction site variants.
#Returns 1 if yes, 0 if no
sub r1_c {
	#declare and initialize, _r means a reference
	my($seq_r, $restr_site1_r) = @_;
	my $restr = '';#one restriction site variant
	my $result = 0;
	
	for $restr (@{$restr_site1_r}) {#loop through restriction site variants
		if ($$seq_r =~ /^$restr/) {#if DNA starts with restriction site
			$result = 1;
			last;
		}
	}
	return $result;
}

#definition of subroutine r1_rep ver.2
#expects references to
#$seq			DNA sequence
#@restr_site1	all sequence variants that match restriction site at beginning of DNA sequence
#				specified by the user (only 1 if restriction site does not contain ambiguities)
#				variants are alphabetically sorted
#$r_reppos		number of positions that may be repaired in restriction site
#Checks if DNA begins with one of the restriction site variants.
#Repairs up to $r_reppos positions in restriction site if necessary.
#If restriction site contains ambiguities:
#Repairs to the variant with the smallest number of mismatches.
#If several variants have same smallest number of mismatches:
#Takes the one that comes first alphabetically.
#Does not repair if smallest number of mismatches > $r_reppos.
#Returns 1 if succesfully repaired, else 0.
sub r1_rep {
	#declare and initialize, _r means a reference
	my($seq_r, $restr_site1_r, $r_reppos_r) = @_;
	my $restr = '';#a restriction site variant
	my $mis = 0;#number of mismatches
	my $minmis = 0;#smallest number of mismatches
	my %misrestr = ();# {$mis} = array of tested $restr
	my $result = 0;
	
	for $restr (@{$restr_site1_r}) {#loop through restriction site variants
		#determine number of mismatches
		$mis = (length $restr) - (($restr ^ substr($$seq_r,0,length $restr)) =~ tr/\0//);
		#store result in %misrestr
		push @{$misrestr{$mis}}, $restr;
	}
	#determine smallest number of mismatches
	$minmis = (sort {$a <=> $b} keys %misrestr)[0];
	if ($minmis <= $$r_reppos_r) {#if smallest number of mismatches tolerable (repairable)
		$restr = ${$misrestr{$minmis}}[0];#get alphabetically first matching variant with minimum mismatch
		substr($$seq_r,0,length $restr) = $restr;#repair DNA
		$result = 1;#good result
	}
	return $result;
}

#definition of subroutine bcr_c ver.2
#expects references to
#$seq	a DNA sequence (string)
#%barcodes_restr: {barcode} = array of all matching variants to $f_restr_site1 specified by user
				#each preceded by this barcode
#checks if DNA sequence starts with one and only one of the barcodes
#followed by one of the matching variants to $f_restr_site1
#no mismatch allowed, N is a mismatch to ACGT
#returns: the matching barcode or 'NN' if none matches or more than one matches
sub bcr_c {
	#declare and initialize, _r means a reference
	my($seq_r, $barcodes_restr_r) = @_;
	my $bc = '';#a barcode
	my $bcr = '';#a barcode followed by a restriction site variant
	my $nmatch = 0;#number of matching barcodes
	my $result = '';#result for return
	
	#loop through the barcodes
	BC: for $bc (keys %{$barcodes_restr_r}) {
		#loop through the barcode-restriction site combinations
		BCR: for $bcr (@{$$barcodes_restr_r{$bc}}) {
			if ($$seq_r =~/^$bcr/) {#if sequence starts with this combination
				++$nmatch;#count a matching barcode
				$result = $bc;#take current barcode as result
				next BC;#continue check with next barcode
			}
		}
	}
	unless ($nmatch == 1) {#unless one and only one barcode matches
		$result = 'NN';#result is NN
	}
	return $result;
}

#definition of subroutine bc_rrep ver.2
#expects references to
#$seq	a DNA sequence (string)
#%barcodes_restr: {barcode} = array of all matching variants to $f_restr_site1 specified by user
				#each preceded by this barcode
				#variants are sorted alphabetically in the array
#$r_reppos		#number of positions that may be repaired in restriction site
#Identifies the single matching combination of a barcode and the following restriction site.
#Repairs up to $r_reppos positions in the restriction site if necessary.
#If the restriction site contains ambiguities: evaluates all possible matching variants.
#If several variants match the sequence with tolerable (repairable) mismatch:
#Takes the variant with the smallest number of mismatches.
#If there are several variants with that number of mismatches:
#Takes the one that comes alphabetically first.
#Returns the matching barcode or 'NN' if none matches or
#if several barcodes would match with repair.
#Does not repair the DNA sequence if it returns 'NN'.
sub bc_rrep {
	#declare and initialize, _r means a reference
	my($seq_r, $barcodes_restr_r, $r_reppos_r) = @_;
	my $bc = '';#a barcode
	my $matchbc = '';#matching barcode
	my $nmatchbc = 0;#number of matching barcodes
	my $bcr = '';#a barcode, followed by a restriction site variant
	my $matchbcr = '';#matching bcr
	my $bestbcr = '';#bcr with smallest number of mismatches that comes alphabetically first
	my $mis = 0;#number of mismatches
	my $minmis = 0;#smallest number of mismatches
	my %misbcr = ();# {$mis}= array of tested $bcr
	my $result = 'NN';#result for return
	
	#loop through the barcodes
	for $bc (keys %{$barcodes_restr_r}) {
		if ($$seq_r =~ /^$bc/){#if sequence begins with that barcode
			#loop through the bcr (this barcode followed by a variant of restr_site1
			for $bcr (@{$$barcodes_restr_r{$bc}}) {
				#determine number of mismatches between $bcr and beginning of sequence
				$mis = (length $bcr) - (($bcr ^ substr($$seq_r,0,length $bcr)) =~ tr/\0//);
				#store result in %misbcr
				push @{$misbcr{$mis}}, $bcr;
			}
			#determine smallest number of mismatches
			$minmis = (sort {$a <=> $b} keys %misbcr)[0];
			if ($minmis <= $$r_reppos_r) {#if number of mismatches acceptable
				++$nmatchbc;#count a matching barcode
				$matchbc = $bc;#record matching barcode
				$matchbcr = ${$misbcr{$minmis}}[0];#record matching bcr
			}
			%misbcr = ();#set back and go to next barcode
		}
	}
	if ($nmatchbc == 1) {#if tests were successful for one single barcode
		#repair DNA sequence
		substr($$seq_r,0,length $matchbcr) = $matchbcr;
		$result = $matchbc;#matching barcode is result for return
	}
	return $result;
}

#definition of subroutine bcrep_r ver.2
#expects references to
#$seq			DNA sequence
#@barcodes		all existing barcodes
#@f_restr_site1	all sequence variants that match restriction site after barcode specified
#				by the user (only 1 if restriction site does not contain ambiguities)
#$bc_reppos		number of positions that may be repaired in barcode
#Identifies the single matching combination of a barcode and the following restriction site.
#Repairs up to $bc_reppos positions in the barcode if necessary.
#If the restriction site contains ambiguities: evaluates all possible matching variants.
#Returns the matching barcode or 'NN' if none matches or
#if several barcodes would match with repair.
#Does not repair the DNA sequence if it returns 'NN'.
sub bcrep_r {
	#declare and initialize, _r means a reference
	my($seq_r, $barcodes_r, $f_restr_site1_r, $bc_reppos_r) = @_;
	my $bc = '';#a barcode
	my $matchbc = '';#matching barcode
	my $nmatchbc = 0;#number of matching barcodes
	my $restr = '';#a restriction site variant
	my $lenrestr = 0;#length of restriction site (equal for all variants)
	my $mis = 0;#number of mismatches between barcode and beginning of DNA sequence
	my $result = 'NN';
	
	#determine length of restriction site (all variants have same length)
	$lenrestr = length ${$f_restr_site1_r}[0];
	BC: for $bc (@{$barcodes_r}) {#loop through the barcodes
		#determine number of mismatches between barcode and beginning of DNA sequence
		$mis = (length $bc) - (($bc ^ substr($$seq_r,0,length $bc)) =~ tr/\0//);
		if ($mis <= $$bc_reppos_r) {#if number of mismatches tolerable (repairable)
			#loop through the restriction site variants
			RESTR: for $restr (@{$f_restr_site1_r}) {
				#if the next bases in DNA sequence match this variant exactly
				if ($restr eq substr($$seq_r,length $bc,$lenrestr)) {
					++$nmatchbc;#count matching barcode
					$matchbc = $bc;#store matching barcode
					next BC;#continue searching with next barcode
				}
			}
		}
	}
	if ($nmatchbc == 1) {#if one and only one barcode matches
		$result = $matchbc;#matching barcode is result
		#repair DNA sequence
		substr($$seq_r,0,length $matchbc) = $matchbc;
	}
	return $result;
}

#definition of subroutine bcrep_rrep ver.2
#expects references to
#$seq			DNA sequence
#@barcodes		all existing barcodes
#@f_restr_site1	all sequence variants that match restriction site after barcode specified
#				by the user (only 1 if restriction site does not contain ambiguities)
#				variants are sorted alphabetically
#$bc_reppos		number of positions that may be repaired in barcode
#$r_reppos		number of positions that may be repaired in restriction site
#Identifies the single matching combination of a barcode and the following restriction site.
#Repairs up to $bc_reppos positions in the barcode if necessary.
#Repairs up to $r_reppos positions in the restriction site if necessary.
#If the restriction site contains ambiguities:
#evaluates all possible matching variants,
#repairs to the variant with the smallest mismatch.
#If several variants have same smallest mismatch:
#takes the one that comes first alphabetically.
#Returns the matching barcode or 'NN' if none matches or
#if several barcodes would match with repair.
#Does not repair the DNA sequence if it returns 'NN'.
sub bcrep_rrep {
	#declare and initialize, _r means a reference
	my($seq_r, $barcodes_r, $f_restr_site1_r, $bc_reppos_r,$r_reppos_r) = @_;
	my $bc = '';#a barcode
	my $matchbc = '';#matching barcode
	my $nmatchbc = 0;#number of matching barcodes
	my $matchbcr = '';#matching barcode followed by selected matching restr.site variant
	my $bcmis = 0;#number of mismatches between barcode and DNA sequence
	my $restr = '';#a restriction site variant
	my $lenrestr = 0;#length of restriction site (equal for all variants)
	my $restrmis = 0;#number of mismatches between restr.site and DNA sequence
	my $minrestrmis = 0;#smallest number of mismatches between restr.site and DNA sequence
	my %restrmis_restr = ();# {$restrmis}=array of tested $restr
	my $result = 'NN';#result for return

	#determine length of restriction site (all variants have same length)
	$lenrestr = length ${$f_restr_site1_r}[0];
	for $bc (@{$barcodes_r}) {#loop through the barcodes
		#determine number of mismatches between barcode and beginning of DNA sequence
		$bcmis = (length $bc) - (($bc ^ substr($$seq_r,0,length $bc)) =~ tr/\0//);
		if ($bcmis <= $$bc_reppos_r) {#if number of mismatches tolerable (repairable)
			#loop through the restriction site variants
			for $restr (@{$f_restr_site1_r}) {
				#determine number of mismatches between
				#restriction site variant and characters in DNA after barcode
				$restrmis = ($lenrestr) - (($restr ^ substr($$seq_r,length $bc,$lenrestr)) =~ tr/\0//);
				#store result in %restrmis_restr
				push @{$restrmis_restr{$restrmis}}, $restr;
			}
			#determine smallest number of mismatches
			$minrestrmis = (sort {$a <=> $b} keys %restrmis_restr)[0];
			if ($minrestrmis <= $$r_reppos_r) {#if number of mismatches tolerable (repairable)
				++$nmatchbc;#count a matching barcode
				$matchbc = $bc;#record matching barcode
				#build seq of matching barcode followed by bestmatching (then alphabetically first)
				#restriction site variant
				$matchbcr = $bc.${$restrmis_restr{$minrestrmis}}[0];
			}
			%restrmis_restr = ();#set back
		}#go to next barcode
	}
	if ($nmatchbc == 1) {#if tests were successful for one single barcode
		#repair DNA sequence
		substr($$seq_r,0,length $matchbcr) = $matchbcr;
		$result = $matchbc;#matching barcode is result for return
	}
	return $result;
}

#definition of subroutine cl_bc
#expects:
#reference to a string (DNA)
#reference to a string (Phred symbols)
#string: barcode
#determines the length of the barcode
#removes the number of characters
#corresponding to the length of the barcode
#from the beginning of DNA and Phred symbols
sub cl_bc {
	#declare and initialize arguments and variables
	my ($seqlineref, $phredlineref,$BC) = @_;
	my ($BCl) = 0;#length of barcode
	$BCl = length $BC;
	$$seqlineref =~ s/^.{$BCl}//;#clip barcode from DNA
	$$phredlineref =~ s/^.{$BCl}//;#clip barcode from Phred symbols
}

#definition of subroutine rint_c ver.2
#expects references to
#$seqline2			DNA-sequence
#@restr_site_arr	array of internal restriction sites to search for
#returns 1 if sequence is good (no internal restriction site)
#returns 0 if sequence is not good (internal restriction site found)
sub rint_c {
	#declare and initialize, _r means a reference
	my($seqline2_r,$restr_site_arr_r) = @_;
	my $restr_site = '';#one restriction site to search for
	my $result = 1;
	
	#loop through restriction sites to search for
	for $restr_site (@{$restr_site_arr_r}) {
		if ($$seqline2_r =~ /$restr_site/) {
			$result = 0;
			return $result;
		}
	}
	return $result;
}

#definition of subroutine ad_c_pos ver.4
#expects references to 3 variables:
#$seqline2		DNA-sequence
#@adseq			array of adapter sequences to search for
#$admis			proportion of adapter sequence that may mismatch
#searches for adapter sequences in DNA sequence
#either complete adapter sequence
#or the begin of it (at least $minpos positions) at the end of the DNA sequence
#allows for mismatch, examples:
#$admis=0.1, adapter in sequence: 10 last positions: int(0.1*10) = 1 mismatch allowed
#$admis=0.1, adapter in sequence: 9 last positions: int(0.1*9) = 0 mismatch allowed
#If more than one adapter sequence (in @adseq) matches, takes the leftmost matching one.
#returns -1, when no adapter was found, else the start-position of the adapter found
#start-position: count starting with 0
sub ad_c_pos {
	#declare and initialize, _r means a reference
	my($seqline2_r,$adseq_r,$admis_r) = @_;
	#own variables
	my ($adseq) = '';#one adapter sequence
	my ($minpos) = 5;#minimum number of positions of adapter that must be there to be found
	my ($mis_tol) = 0;#number of mismatch positions tolerated
	my ($mis) = 0;#number of mismatches found
	my ($stopdiff) = 1;# $mis - $mis_tol if it is <=0: stop everything, return: adapter found
	my ($DNAstart) = 0;#current start position in DNA-sequence
	my ($DNAend) = 0;#current endposition in DNA-sequence
	my ($DNApart_l) = 0;#length of DNA-part to search in for adapter
	my (@adstartpositions) = ();#start positions of found adapters, count starting with 0
	my ($result) = -1;#return-value
	
	for $adseq (@{$adseq_r}) {#loop through adapter sequences
		#determine some variables
		$DNAend = (length $adseq) - 1;#changes later
		$DNApart_l = length $adseq;#changes later
		$mis_tol = int($$admis_r * $DNApart_l);#changes later

		#as long as stopdiff > 0: work
		#first, as long as possible, search for complete adapter sequence with mismatch
		while (($stopdiff > 0) and ($DNAend < length $$seqline2_r)) {
			#determine number of mismatches between DNApart and adseq
			$mis = $DNApart_l - ((substr($$seqline2_r,$DNAstart,$DNApart_l) ^ $adseq) =~ tr/\0//);
			$stopdiff = $mis - $mis_tol;#calculate stopdiff
			$mis = 0;#set $mis back
			++$DNAstart;
			++$DNAend;
		}
		#when you did not find an adapter up to now, DNAend is one position to big now
		#decrement it
		--$DNAend;
		#continue searching further right in DNA in smaller and smaller DNA parts
		#DNA_part_l must be one position shorter now
		--$DNApart_l;
		#recalculate number of mismatches allowed
		$mis_tol = int($$admis_r * $DNApart_l);
		#continue searching while DNA_part_l is >= $minpos
		while (($stopdiff > 0) and ($DNApart_l >= $minpos)) {
			#determine number of mismatches between DNApart and corresponding part of adseq
			$mis = $DNApart_l - ((substr($$seqline2_r,$DNAstart,$DNApart_l) ^ substr($adseq,0,$DNApart_l)) =~ tr/\0//);
			$stopdiff = $mis - $mis_tol;#calculate stopdiff
			$mis = 0;#set $mis back
			++$DNAstart;
			--$DNApart_l;
			$mis_tol = int($$admis_r * $DNApart_l);#recalculate number of mismatches allowed
		}
		#finally, if the last stopdiff was <= 0, you found an adapter,
		#last DNAstart is the start position
		if ($stopdiff <= 0) {
			push @adstartpositions, ($DNAstart - 1);#store in @adstartpositions
		}
		#set variables back
		$stopdiff = 1;
		$DNAstart = 0;
		#continue search with next adapter sequence
	}
	if (@adstartpositions > 0) {#if you found at least one adapter start position
		$result = (sort {$a <=> $b} @adstartpositions)[0];#take the smallest one as result
	}
	return $result;
}

#definition of sub cl_r ver.2
#expects references to
#$seq			DNA sequence
#$phr			Phred symbol line
#@restr_site1	all sequence variants that match restriction site at beginning of DNA sequence
#				specified by the user (only 1 if restriction site does not contain ambiguities)
#Checks if DNA begins with one of the restriction site variants.
#No mismatch allowed. N is mismatch to ACGT.
#If not: returns 0.
#If yes: removes the restriction site from DNA and phredsymbol-line
#returns 1
sub cl_r {
	#declare and initialize, _r means a reference
	my($seq_r,$phr_r,$restr_site1_r) = @_;
	my $restr = '';#a restriction site variant
	my $result = 0;
	
	for $restr (@{$restr_site1_r}) {#loop through restriction site variants
		if ($$seq_r =~ /^$restr/) {#if DNA starts with this variant
			substr($$seq_r,0,length $restr) = '';#clip from DNA
			substr($$phr_r,0,length $restr) = '';#clip from Phred symbols
			$result = 1;#good result
		}
	}
	return $result;
}

#definition of subroutine trunc
#expects:
#reference to a string (DNA)
#reference to a string (phredsymbols)
#reference to a positive integer (tr)
#truncates the two strings to the length corresponding to tr
#assumes that the two strings have the same length!
#does nothing if tr is greater than this length
sub trunc {
	my ($seqlineref, $phredlineref, $trref) = @_;
	my ($seql) = 0;#lengths of the two strings
	$seql = length $$seqlineref;
	if ($$trref < $seql) {#if tr < input-length, truncate
		$$seqlineref = substr($$seqlineref,0,$$trref);
		$$phredlineref = substr($$phredlineref,0,$$trref);
	}
}

#definition of subroutine rcdupa
#acts on paired reads,
#expects references to 4 scalars:
#DNA sequence of first read
#Phred symbols of first read
#DNA sequence of second read
#Phread symbols of second read
#replaces data of second read with reverse complement of first read
#knows characters: ACGTNacgtn
sub rcdupa {
	#declare and initialize, _r means a reference
	my ($DNA1_r,$phred1_r,$DNA2_r,$phred2_r) = @_;
	
	#replace DNA2 with revcomp of DNA1
	$$DNA2_r = reverse $$DNA1_r;
	$$DNA2_r =~ tr/ACGTacgt/TGCAtgca/;
	#replace phredscores2 with reverse of phredscores1
	$$phred2_r = reverse $$phred1_r;
}

#definition of subroutine catsrc
#acts on paired reads
#expects references to
#fDNA: DNAsequence of forward read
#fphred: Phred symbols of forward read
#rDNA: DNAsequence of reverse read
#rphred: Phred symbols of reverse read
#truncates both sequences to equal length
#appends reverse complement of reverse read to forward read
#replaces reverse read with reverse complement of new concatenated forward read
#knows characters: ACGTNacgtn
sub catsrc {
	#declare and initialize, _r means a reference
	my ($fDNA_r,$fphred_r,$rDNA_r,$rphred_r) = @_;
	my $flen = 0;#length of forward seq
	my $rlen = 0;#length of reverse seq
	my $tempstring1 = '';
	
	#truncate sequences to equal length if necessary
	$flen = length $$fDNA_r;
	$rlen = length $$rDNA_r;
	#if reverse-read is longer, truncate reverse-read
	if ($rlen > $flen) {
		$$rDNA_r = substr($$rDNA_r,0,$flen);
		$$rphred_r = substr($$rphred_r,0,$flen);
	}
	#if forward-read is longer, truncate forward-read
	elsif ($flen > $rlen) {
		$$fDNA_r = substr($$fDNA_r,0,$rlen);
		$$fphred_r = substr($$fphred_r,0,$rlen);
	}
	#append revcomp of reverse read to forward read
	$tempstring1 = reverse $$rDNA_r;
	$tempstring1 =~ tr/ACGTacgt/TGCAtgca/;
	$$fDNA_r .= $tempstring1;
	$tempstring1 = reverse $$rphred_r;
	$$fphred_r .= $tempstring1;
	#replace reverse read with revcomp of new concatenated forward read
	$$rDNA_r = reverse $$fDNA_r;
	$$rDNA_r =~ tr/ACGTacgt/TGCAtgca/;
	$$rphred_r = reverse $$fphred_r;
}
