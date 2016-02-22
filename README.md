# GIbPSs

Genotype Individuals based on Paired Sequences or single reads

## Introduction
GIbPSs is a package for the analysis of GBS (genotyping-by-sequencing) and RAD (restriction site associated DNA) sequence data. It handles paired-end or single reads from protocols such as GBS, two-enzyme GBS, and ddRAD. It can also analyze forward reads from RAD libraries. Its main purpose is genotyping of individuals from diploid organisms without a reference genome for population genetic and phylogeographic analyses. GIbPSs is written in Perl 5 and runs under Linux and Windows. I have developed it for the analysis of sequence data from "classical" GBS libraries (Elshire et al. 2011). This technique uses a single restriction enzyme. Reads of a locus originate from both strands and display full, partial, or no overlap. This can lead to duplicate locus identification when a reference genome is not available. GIbPSs solves this problem through reverse complement clustering of paired sequence reads. GIbPSs can handle sequences of different lengths and analyze short restriction fragments in adapter-containing sequences. GIbPSs employs a rapid genotyping strategy that eliminates errors due to indel variation without the need for time-consuming multiple alignments. It identifies loci with indel variation with the aid of USEARCH (Edgar, R. C. 2010) and excludes them from the final dataset. GIbPSs stores genotype data in a database that offers multiple options for data exploration, quality filtering, data selection, and export. It produces various data formats for population genetic and phylogenetic programs. Examples are STRUCTURE, Genepop, Arlequin, and files in FASTA, Nexus, and Phylip format that contain DNA characters, integers, or binary data. GIbPSs can export one SNP per locus, encode alleles based on all SNPs of each locus or export full allele sequences.

## Downloads
The latest release is available [here](https://github.com/ahapke/gibpss/releases).

The downloadable archive contains a detailed documentation in directory docu. The first document, 01_GIbPSs_1_0_overview.pdf contains an overview of GIbPSs, installation instructions, a tutorial for a quick exploration of GIbPSs with simulated data, and instructions for getting started with real data. The remaining documents explain the programs of GIbPSs.

## Installation instructions
All programs of GIbPSs are written in Perl 5. I have tested them with Strawberry Perl 5.16.2 under Windows 7 and with Perl 5.18.2 under Ubuntu 14.04 LTS. Perl must be installed on your system to use GIbPSs. If you use Linux, it is probably already installed. If you use Windows, you can download Strawberry Perl at http://strawberryperl.com/. To check, which version of Perl is installed, enter the following in a command prompt: perl -v i

GIbPSs needs the following Perl modules: IO::File, IO::Uncompress::Gunzip, IO::Compress::Gzip, List::Util, Math::Round, Parallel::ForkManager. The documentation contains instructions how to see if a module is already installed and how to install a module.

USEARCH (Edgar, R.C. 2010) must be installed on your system. The program indel_checker, which is part of GIbPSs, uses it to identify loci with indel variation. USEARCH is available for Linux, Windows, and Mac OSX at http://www.drive5.com/usearch. I have tested indel_checker with USEARCH v7.0.1090.

The programs of GIbPSs themselves require no installation. Download the archive, unpack it and store it anywhere on your computer. Directory source contains the programs. Each program is a single file that you can copy to a different location and execute there. Several programs need additional text files with settings that you must prepare for your analyses. Directories user_files_Lin and user_files_Win contain example files for Linux and Windows. Under Linux, all text files read by GIbPSs must have Linux line endings. If you use Strawberry Perl under Windows, they may have Linux or Windows line endings. The program files themselves have Windows line endings but they run under Windows and Linux.

## License
GIbPSs is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GIbPSs is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The downloadable archive contains a copy of the GNU General Public License version 3 in file License.txt. You can also obtain a copy at http://www.gnu.org/licenses/.

## Citation
Please cite GIbPSs as follows:
Hapke A, Thiele D (2016) GIbPSs: a toolkit for fast and accurate analyses of genotyping-by-sequencing data without a reference genome. Molecular Ecology Resources. doi: 10.1111/1755-0998.12510

## Contact
Andreas Hapke, Johannes Gutenberg University Mainz, Institute of Anthropology, Anselm-Franz-von-Benzel-Weg 7, D-55099 Mainz, Germany, ahapke2@gmail.com.

Please contact me via e-mail when you have questions or suggestions concerning GIbPSs. Any bug reports are highly welcome as well.

## References
Edgar, R.C.(2010) Search and clustering orders of magnitudes faster than BLAST, *Bioinformatics* 26(19), 2460-2461.

Elshire, R.J., Glaubitz, J.C., Sun, Q., Poland, J.A., Kawamoto, K., Buckler, E.S., Mitchell, S.E. (2011) A robust, simple Genotyping-by-Sequencing (GBS) approach for high diversity species. *PLoS ONE* 6(5): e19379.

