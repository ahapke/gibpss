#!/usr/bin/perl -w
#depth_analyzer version 09.0 Copyright 2015 Andreas Hapke
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

#uses subroutine:
#slr_da

#Keyword "Infilecolumns!" Column order of an infile is important.

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#Copyright message to screen
print "depth_analyzer version 09.0 Copyright 2015 Andreas Hapke\n",
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
my $poploc_name = 'poploc.txt';#out of poploc: poplocID sl cons nSNP varpos n_all nInd nIndloc ties lon_over lon_com_over merge_conflict
my $sel_ind_name = 'export/sel_ind.txt';#out of data_selector: ind
my $sel_loc_name = 'export/sel_loc.txt';#out of data_selector: poplocID
my $sel_gt_name = 'export/sel_gt.txt';#out of data_selector: poplocID ind
my $splitloc_name = 'split_loci.txt';#out of indpoploc: poplocID ind indlocID
my $loci_suffix = '_loci.txt';#suffix of loci-outfile of indloc
my $indpoploc_suffix = '_indpoploc.txt';#suffix of indpoploc outfile of indpoploc
my $clusters_name = 'indelcheck/clusters.txt';#clusters outfile of indel_checker
my $all_loc_outname = 'depth_analysis/da_all_loc.txt';#outfile with data about all loci across seq-length ranges
#other variables
my %poploc_sl = ();# {poplocID} = sl
my %clust_poploc = ();# {clustNo}{poplocID}=sl clusters identified by indelchecker
					#sl: length of poploc
my $clustNo = 0;#number of a cluster
my %poploc_clust = ();# {poplocID} = (clustNo,clust_sl) clust_sl: length of longest poploc in cluster
my $clust_sl = 0;#length of longest poploc in a cluster: cluster seqlength
my $slrNo = 0;#consecutive number of sequence-length-range
my %slr = ();# {slrNo} = (slmin,slmax)
my $slmin = 0;#minimum sequence length
my $slmax = 0;#maximum sequence length
my $sl = 0;#sequence length
my $ind = '';#an individual ID
my $thisind = '';#current individual ID
my %sel_ind = ();# {ind}=1
my %sel_loc = ();# {poplocID}=1
my %sel_gt = ();# {poplocID}{ind}=1
my $poplocID = 0;#population-level locus ID
my %splitloc = ();# {poplocID}{ind} = indlocID
my $n_splitloc = 0;#number of split loci
my %indlocdata = ();#temp. struc. for one ind: {indlocID} = Loc_dep
my %clustsl_cl_i_d = ();# {clust_sl}{clustNo}{ind} = clust_dep
my $depth_i = 0;#depth of a locus in an individual
my %da_all_loc = ();# {poplocID}= (sl,min_scaldep,med_scaldep,max_scaldep,
					#min_dep_perc,med_dep_perc,max_dep_perc)
					#collects data for all loci across seq-length-ranges
					#for outfile da_all_loc.txt
my $cl_or_pl = 'cl';#analysis based on clusters: cl, analysis based on poplocs: pl
my $infile = '';
my $tempstring1 = '';
my @temparr1 = ();
#}

#Welcome
print "\nDEPTH_ANALYZER\n\n";

##################
#Read in sl-ranges
##################

#{
if (@ARGV) {#if user provided any argument
	$infile = shift @ARGV;
} else {#ask user for infile name
	print "I need a textfile with sequence-length ranges:\n",
	"Each row: min TAB max\n",
	"no headers, blank lines or extra spaces\n",
	"Please enter filename: ";
	$infile = <STDIN>;
	chomp $infile;
}
unless(open(INFILE,$infile)) {
	print "Can't open $infile. Exiting..\n";
	exit;
}
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#check for plausible values
	#Infilecolumns!
	unless(($temparr1[0] =~ /^\d+$/) and ($temparr1[1] =~ /^\d+$/)
		and ($temparr1[1] >= $temparr1[0]) and ($temparr1[1] > 0)) {
		print "I can't use these values: min: $temparr1[0] max: $temparr1[1].Exiting..\n";
		exit;
	}
	#store plausible values in %slr
	@{$slr{$slrNo}} = @temparr1;
	++$slrNo;
}
close INFILE;
#}

###############################################################
#If directory depth_analysis already exists: exit, else: create
###############################################################

#{
if (-d "depth_analysis") {
	print "Directory depth_analysis already exists.\n",
	"Please delete or rename. Exiting..\n";
	exit;
} else {#if not, try to create
	unless (mkdir "depth_analysis") {
		print "Can't create directory depth_analysis. Exiting..\n";
		exit;
	}
}
#}

###############################
#Read in existing preselections
###############################

#{
#if directory "export" exists, check for selection files and read them in
if (-d "export") {
	if (-f "$sel_ind_name") {#if there is an individual selection-file
		print "Loading your selection of individuals.\n";
		unless(open(INFILE,"$sel_ind_name")) {#open it
			print "Cannot open $sel_ind_name. Exiting..\n";
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
			"I can't continue. Exiting..\n";
			exit;			
		}
	}
	if (-f "$sel_loc_name") {#if there is a locus selection-file
		print "Loading your selection of loci.\n";
		unless(open(INFILE,"$sel_loc_name")) {#open it
			print "Cannot open $sel_loc_name. Exiting..\n";
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
			"I can't continue. Exiting..\n";
			exit;			
		}
	}
	if (-f "$sel_gt_name") {#if there is a genotype selection-file
		print "Loading your selection of genotypes.\n";
		unless(open(INFILE,"$sel_gt_name")) {#open it
			print "Cannot open $sel_gt_name. Exiting..\n";
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
			"I can't continue. Exiting..\n";
			exit;			
		}
	}	
}
#}

##############################################################
#If no preselections exist, create selection of whole database
##############################################################

#{
#if %sel_ind is still empty, create it from file individuals.txt (using all inds)
if (keys %sel_ind == 0) {
	print "Found no selection of individuals, select all in the database.\n";
	unless(open(INFILE,"$individuals_fname")) {#open it
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
#if %sel_loc is still empty, read in poploc.txt outfile of poploc
#populate %sel_loc
if (keys %sel_loc == 0) {
	print "Found no selection of loci, select all in the database.\n";
	unless(open(INFILE,$poploc_name)) {
		print "Cannot open $poploc_name. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#skip headerline
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		#Infilecolumns!
		@temparr1 = split(/\t/,$tempstring1);
		$poplocID = shift @temparr1;#get poplocID
		$sel_loc{$poplocID} = 1;#store in %sel_loc
	}
	close INFILE;
}
#if %sel_gt is still empty, create from %sel_loc and %sel_ind
#(using all genotypes from selected inds and locs)
if (keys %sel_gt == 0) {
	print "Found no selection of genotypes, select all in the database.\n";
	for $poplocID (keys %sel_loc) {
		for $ind (keys %sel_ind) {
			$sel_gt{$poplocID}{$ind} = 1;
		}
	}
}
#}

#####################################################
#Read in split_loci.txt and check if any are selected
#If yes: inform user that they will be ignored
#####################################################

#{
unless(open(INFILE,$splitloc_name)) {
	print "Cannot open $splitloc_name. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	$splitloc{$temparr1[0]}{$temparr1[1]} = $temparr1[2];
}
close INFILE;
#Check how may split loci are in the current selection
for $poplocID (keys %splitloc) {
	if (defined $sel_loc{$poplocID}) {
		++$n_splitloc;
	}	
}
#If you found any, exclude them from analysis, inform user
if ($n_splitloc > 0) {
	print "The current selection includes $n_splitloc split loci.\n",
	"I will exclude these loci from depth analysis.\n";
	for $poplocID (keys %splitloc) {
		delete $sel_loc{$poplocID};#delete from selected loci
		delete $sel_gt{$poplocID};#delete from selected genotypes
	}
	#check for empty individuals and exclude
	#set all individuals in %sel_ind to 0
	for $ind (keys %sel_ind) {
		$sel_ind{$ind} = 0;
	}
	#loop through selected individuals
	IND: for $ind (keys %sel_ind) {
		#loop through loci in %sel_gt
		POPLOC: for $poplocID (keys %sel_gt) {
			#loop through individuals in this locus in %sel_gt
			THISIND: for $thisind (keys %{$sel_gt{$poplocID}}) {
				$sel_ind{$thisind} = 1;#set this ind to 1 in %sel_ind
				last POPLOC if ($sel_ind{$ind} == 1);#go to next ind in outer loop once you found ind
			}			
		}
		#delete the individual if it was not found
		if ($sel_ind{$ind} == 0) {
			delete $sel_ind{$ind};
			print "No data left for individual $ind. I exclude it.\n";
		}
	}	
}
#}

######################################################
#Read in poploc.txt, populate %poploc_sl {poplocID}=sl
######################################################

#{
unless(open(INFILE,$poploc_name)) {
	print "Cannot open $poploc_name. Exiting..\n";
	exit;
}
$tempstring1 = <INFILE>;#skip headerline
while ($tempstring1 = <INFILE>) {#read in data
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	if (defined $sel_loc{$temparr1[0]}) {#if this poploc is in selection
		$poploc_sl{$temparr1[0]} = $temparr1[1];
	}
}
close INFILE;
#}

############################################################
#If exists, read in indelcheck/clusters.txt
#if or if not: populate %clust_poploc {clustNo}{poplocID}=sl
#and %poploc_clust {poplocID}=(clustNo,clust_sl)
#set back %poploc_sl
############################################################

#{
if(open(INFILE,$clusters_name)) {#if you find and can open indelcheck/clusters.txt
	print "Reading in file $clusters_name\n",
	"Depth analysis will be based on clusters of loci.\n",
	"Length of a cluster is defined as length of its longest poploc.\n";
	$tempstring1 = <INFILE>;#read to skip headerline
	while ($tempstring1 = <INFILE>) {#read in data
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		if (defined $sel_loc{$temparr1[1]}) {#if current locus is in selection
			#store in %clust_poploc with sl of this poploc
			$clust_poploc{$temparr1[0]}{$temparr1[1]} = $poploc_sl{$temparr1[1]};
		}
	}
	close INFILE;
	if (keys %clust_poploc == 0) {#if hash still empty
		print "Found no data in file $clusters_name.\n";
	} else {
		#determine length of longest poploc in each cluster and populate %poploc_clust
		@temparr1 = ();
		for $clustNo (keys %clust_poploc) {
			for $poplocID (keys %{$clust_poploc{$clustNo}}) {
				push @temparr1, $clust_poploc{$clustNo}{$poplocID};
			}
			$clust_sl = (reverse sort {$a <=> $b} @temparr1)[0];
			@temparr1 = ();
			for $poplocID (keys %{$clust_poploc{$clustNo}}) {
				@{$poploc_clust{$poplocID}} = ($clustNo,$clust_sl);
			}		
		}
		#Check if current selection of loci contains loci that did not appear
		#in the output out of indel_checker ($clusters_name)
		#If so, inform user and exit.
		for $poplocID (keys %sel_loc) {
			unless (defined $poploc_clust{$poplocID}) {
				print "Your current selection includes loci that do not appear\n",
				"in the output of indelchecker (file $clusters_name)\n",
				"I cannot perform depth analysis based on clusters.\n",
				"Please do one of the following before using depth_analyzer:\n",
				"Remove or rename directory indelcheck to perform\n",
				"depth analysis based on loci.\n",
				"If you wish to perform depth analysis based on clusters:\n",
				"Remove or rename directory export, use data_selector\n",
				"to make a new selection, which is the same that you made\n",
				"before you used indel_checker.\n",
				"See report file out of indel_checker for your previous selection.\n",
				"And, in any case: Please remove directory depth_analysis.\n",
				"Exiting..\n";
				exit;
			}
		}
	}
}
if (keys %clust_poploc == 0) {#if hash still empty
	print "Found no clusters of loci identified by indel_checker.\n",
	"Depth analysis will be based on loci.\n";
	$cl_or_pl = 'pl';
	#build one cluster for each locus, use poplocID as cluster number
	#populate %clust_poploc and %poploc_clust
	for $poplocID (keys %poploc_sl) {
		$clust_poploc{$poplocID}{$poplocID} = $poploc_sl{$poplocID};
		@{$poploc_clust{$poplocID}} = ($poplocID,$poploc_sl{$poplocID});
	}
}
%poploc_sl = ();#set back - not needed any more
#}

##########################################
#Read in data from selected individuals
#populate %clustsl_cl_i_d
#set back %poploc_clust, %sel_loc, %sel_gt
##########################################

#{
for $ind (keys %sel_ind) {#Loop through selected individuals
	#read in *_loci.txt out of indloc: indlocID Loc_cat sl Loc_dep n_alleles lost_alleles potlostmore
	$infile = $ind.$loci_suffix;#create infilename
	unless(open(INFILE,$infile)) {#open or die
		print "Can't open $infile. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#rskip headerline
	#read data and populate %indlocdata for this ind: {indlocID} = (sl,Loc_dep)
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		#Infilecolumns!
		$indlocdata{$temparr1[0]} = $temparr1[3];
	}
	close INFILE;
	#read in *_indpoploc out of indpoploc: indlocID poplocID
	$infile = $ind.$indpoploc_suffix;#create infilename
	unless(open(INFILE,$infile)) {#open or die
		print "Can't open $infile. Exiting..\n";
		exit;
	}
	$tempstring1 = <INFILE>;#read to skip headerline
	while ($tempstring1 = <INFILE>) {#read in data
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		#Infilecolumns!
		if ((defined $sel_gt{$temparr1[1]}) and (defined $sel_gt{$temparr1[1]}{$ind})) {#if poploc for this ind is included in selection
			#look up depth of this individual locus in this ind
			$depth_i = $indlocdata{$temparr1[0]};
			#look up cluster number of corresponding poploc and clust_sl (length of longest poploc in cluster)
			$clustNo = $poploc_clust{$temparr1[1]}[0];
			$clust_sl = $poploc_clust{$temparr1[1]}[1];
			#store data in %clustsl_cl_i_d, store $depth_i by addition:
			#There may be several indlocs in this ind for this poploc
			#There may be several poplocs in this cluster
			$clustsl_cl_i_d{$clust_sl}{$clustNo}{$ind} += $depth_i;			
		}
	}
	close INFILE;
	%indlocdata = ();#set back and next ind
}
%poploc_clust = ();#set back not needed any more
%sel_loc = ();#set back not needed any more
%sel_gt = ();#set back not needed any more
#}

###############################################################
#Loop through sl-ranges and call slr_da to populate %da_all_loc
###############################################################

#{
$slrNo = 0;#start with first range
for $slrNo (sort {$a <=> $b} keys %slr) {
	$slmin = $slr{$slrNo}[0];#determine min sl
	$slmax = $slr{$slrNo}[1];#determine max sl
	slr_da($slmin,$slmax,\%clustsl_cl_i_d,\%da_all_loc,\%clust_poploc,$cl_or_pl);
}
#}

############################################################
#Produce outfile over all loci across sequence length ranges
############################################################

#{
unless(open(OUTFILE, ">$all_loc_outname")) {#open or die
	print "Cannot open $all_loc_outname, exiting ...\n\n";
	exit;
}
#print headerline
print OUTFILE "poplocID\tsl\tmin_scaldep\tmed_scaldep\tmax_scaldep\t",
"min_dep_perc\tmed_dep_perc\tmax_dep_perc\n";
#loop through loci in %da_all_loc and print data
for $poplocID (sort {$a <=> $b} keys %da_all_loc) {
	$tempstring1 = join("\t",@{$da_all_loc{$poplocID}});
	print OUTFILE "$poplocID\t$tempstring1\n";
}
close OUTFILE;
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "run took $run_s seconds.\n";

exit;

############
#Subroutines
############

#definition of subroutine slr_da
#performs analyses for a given range of sequence lengths (sl)
#produces outfiles for this sl-range
#populates %da_all_loc in main
#expects
#$slmin: min sl
#$slmax: max sl
#ref to %clustsl_cl_i_d {clust_sl}{clustNo}{ind}=clust_dep
#ref to %da_all_loc: {poplocID}= (sl,min_scaldep,med_scaldep,max_scaldep,
								#min_dep_perc,med_dep_perc,max_dep_perc)
#ref to %clust_poploc: {clustNo}{poplocID}=sl (sl of poploc)
#$cl_or_pl: when analysis is based on clusters: 'cl', based on poplocs: 'pl'
sub slr_da {
	#declare and initialize, _r means a reference
	my ($slmin, $slmax, $clustsl_cl_i_d_r,$da_all_loc_r,$clust_poploc_r) = @_;
	my $dep = 0;#absolute depth of a locus in an individual
	my $scaldep = 0;#depth within this ind by interval scale variable ranging
	my $minscaldep = 0;#minimum scaldep across inds
	my $medscaldep = 0;#median of scaldep across inds
	my $maxscaldep = 0;#maximum scaldep across inds
	my $depperc = 0;#percent of loci with greater depth in this ind
	my $mindepperc = 0;#minimum depperc across inds
	my $meddepperc = 0;#median of depperc across inds
	my $maxdepperc = 0;#maximum depperc across inds
	my $mindep = 0;#smallest depth in an individual
	my $maxdep = 0;#greatest depth in an individual
	my $deprange = 0;#maxdep - mindep
	my %i_d_cl = ();# {ind}{dep}{clustNo} = 1
	my %i_d_data = ();# {ind}{dep} = (scaldep,depperc)
	my %cl_i_data = ();# {clustNo}{ind} = (dep,scaldep,depperc)
						#dep: absolute depth
						#scaldep: depth within this ind by interval scale variable ranging
						#depperc: percent of loci with greater depth in this ind
	my %cl_data = ();# {clustNo} = (minscaldep,medscaldep,maxscaldep,mindepperc,meddepperc,maxdepperc)
					#mscaldep: median of scaldep across inds
					#mdepperc: median of depperc across inds
	my $sl = 0;#sequence length
	my $clustNo = 0;#cluster number
	my $poplocID = 0;#population-level locus ID
	my $ind = '';#individual ID
	my $nclust = 0;#number of clusters of loci in an individual
	my $nclustcum = 0;#cumulative number of clusters including current dep
	my @scaldepinds = ();#scaled depths of a locus of all inds
	my @deppercinds = ();#depperc of a locus of all inds
	my $n_ind = 0;#number of individuals
	my $i = 0;
	my $outfilename = '';
	my @temparr1 = ();
	my $tempstring1 = '';
	
	###########################################
	#populate %i_d_cl: {ind}{dep}{clustNo} = 1
	###########################################
	
	#{
	$sl = $slmin;#start with the min sequence length and go through range up to slmax
	while ($sl <= $slmax) {
		for $clustNo (keys %{$$clustsl_cl_i_d_r{$sl}}) {#go through clustNos in %clustsl_cl_i_d
			for $ind (keys %{$$clustsl_cl_i_d_r{$sl}{$clustNo}}) {#go through inds for that clustNo
				$dep = $$clustsl_cl_i_d_r{$sl}{$clustNo}{$ind};#get dep for this ind and cluster
				$i_d_cl{$ind}{$dep}{$clustNo} = 1;#store in %i_d_cl
			}
		}	
		++$sl;#increment sl for next round
	}
	#}
	
	####################################################
	#populate
	#%i_d_data: {ind}{dep} = (scaldep,depperc)
	#%cl_i_data: {clustNo}{ind} = (dep,scaldep,depperc)
	####################################################
	
	#{
	#loop through individuals
	for $ind (keys %i_d_cl) {
		#determine mindep, maxdep and deprange
		$mindep = (sort {$a <=> $b} keys %{$i_d_cl{$ind}})[0];
		$maxdep = (reverse sort {$a <=> $b} keys %{$i_d_cl{$ind}})[0];
		$deprange = $maxdep - $mindep;
		#determine number of clusters $nclust in this ind
		#loop through depths in this ind
		for $dep (keys %{$i_d_cl{$ind}}) {
			$nclust += keys %{$i_d_cl{$ind}{$dep}};#count loci			
		}
		#determine scaldep and depperc for each dep in this ind
		#and for each cluster in this ind
		#loop through sorted depths in this ind
		for $dep (sort {$a <=> $b} keys %{$i_d_cl{$ind}}) {
			$scaldep = ($dep - $mindep) / $deprange;#determine scaldep
			$nclustcum += keys %{$i_d_cl{$ind}{$dep}};#update cumulative cluster count
			$depperc = ($nclust - $nclustcum) / $nclust * 100;#determine depperc
			@{$i_d_data{$ind}{$dep}} = ($scaldep,$depperc);#populate %i_d_data
			#loop through clusters with this dep in this ind
			for $clustNo (keys %{$i_d_cl{$ind}{$dep}}) {
				@{$cl_i_data{$clustNo}{$ind}} = ($dep,$scaldep,$depperc);#populate %cl_i_data
			}
		}		
		#set variables back
		$nclust = 0;
		$nclustcum = 0;
		#next ind
	}	
	#}
	
	#########################################################################
	#populate %cl_data: {clustNo} = (min med max scaldep,min med max depperc)
	#########################################################################
	
	#{
	#loop through clustNos in %cl_i_data: {clustNo}{ind} = (dep,scaldep,depperc)
	for $clustNo (keys %cl_i_data) {
		#loop through individuals and collect scaled depths and depperc for this cluster
		for $ind (keys %{$cl_i_data{$clustNo}}) {
			push @scaldepinds, $cl_i_data{$clustNo}{$ind}[1];
			push @deppercinds, $cl_i_data{$clustNo}{$ind}[2];
		}
		#sort the 2 arrays
		@scaldepinds = sort {$a <=> $b} @scaldepinds;
		@deppercinds = sort {$a <=> $b} @deppercinds;
		$n_ind = @scaldepinds;#determine n_ind for this locus
		$minscaldep = $scaldepinds[0];#determine minimum scaldep
		$maxscaldep = $scaldepinds[$n_ind-1];#determine maximum scaldep
		$mindepperc = $deppercinds[0];#determine minimum depperc
		$maxdepperc = $deppercinds[$n_ind-1];#determine maximum depperc
		if ($n_ind % 2 == 1) {#if the number of individuals is odd:
			$i = (($n_ind / 2) - 0.5);#determine index for middle value
			$medscaldep = $scaldepinds[$i];#get median scaldep
			$meddepperc = $deppercinds[$i];#get median depperc
		} else {#if it is even
			$i = (($n_ind / 2) - 1);#get index for first value to calculate median
			$medscaldep = ($scaldepinds[$i] + $scaldepinds[$i+1]) / 2;#calculate median scaldep
			$meddepperc = ($deppercinds[$i] + $deppercinds[$i+1]) / 2;#calculate median depperc
		}
		#populate %cl_data
		@{$cl_data{$clustNo}} = ($minscaldep,$medscaldep,$maxscaldep,
		$mindepperc,$meddepperc,$maxdepperc);	
		#set variables back
		@scaldepinds = ();
		@deppercinds = ();
		#next clustNo
	}	
	#}
	
	#######################################
	#produce an outfile for each individual
	#######################################
	
	#{
	#loop through individuals in %i_d_data
	for $ind (keys %i_d_data) {
		$outfilename = 'depth_analysis/'.$slmin.'_'.$slmax.'_da_'.$ind.'.txt';#create outfilename
		unless(open(OUTFILE, ">$outfilename")) {#open or die
			print "Cannot open $outfilename, exiting ...\n\n";
			exit;
		}
		print OUTFILE "dep\tscaldep\tdep_perc\n";#print headerline
		#loop through sorted depths for this ind and print data
		for $dep (sort {$a <=> $b} keys %{$i_d_data{$ind}}) {
			print OUTFILE "$dep\t$i_d_data{$ind}{$dep}[0]\t$i_d_data{$ind}{$dep}[1]\n";
		}
		close OUTFILE;
	}	
	#}
	
	############################################################
	#produce an outfile with min median max over all individuals
	############################################################
	
	#{
	$outfilename = 'depth_analysis/'.$slmin.'_'.$slmax.'_da_allind.txt';#create outfilename
	unless(open(OUTFILE, ">$outfilename")) {#open or die
		print "Cannot open $outfilename, exiting ...\n\n";
		exit;
	}
	#print headerline
	if ($cl_or_pl eq 'cl') {#if analysis is based on clusters
		print OUTFILE "clustNo\t";
	} else {#analysis is based on poplocs, cluster numbers are poplocIDs
		print OUTFILE "poplocID\t";
	}
	print OUTFILE "min_scaldep\tmed_scaldep\tmax_scaldep\t",
	"min_dep_perc\tmed_dep_perc\tmax_dep_perc\n";
	#loop through sorted clustNos in %cl_data and print data
	for $clustNo (sort {$a <=> $b} keys %cl_data) {
		@temparr1 = @{$cl_data{$clustNo}};#get data (min med max) into array
		$tempstring1 = join ("\t",@temparr1);#join with tab
		print OUTFILE "$clustNo\t$tempstring1\n";
	}
	close OUTFILE;
	#}

	#############################
	#populate %da_all_loc in main
	#############################
	
	#{
	#loop through sl and clustNo in %clustsl_cl_i_d within current seq-length range
	$sl = $slmin;#start with the min sequence length and go through range up to slmax
	while ($sl <= $slmax) {
		for $clustNo (keys %{$$clustsl_cl_i_d_r{$sl}}) {#go through clustNos in %clustsl_cl_i_d
			#get data for this cluster out of %cl_data
			@temparr1 = @{$cl_data{$clustNo}};
			unshift (@temparr1, $sl);#add sl as first element to array
			#loop through poplocIDs in this cluster
			for $poplocID (keys %{$$clust_poploc_r{$clustNo}}) {
				#store poplocID with data in @temparr1 in %da_all_loc
				@{$$da_all_loc_r{$poplocID}} = @temparr1;
			}
		}
		++$sl;#increment sl for next round
	}	
	#}	
}
