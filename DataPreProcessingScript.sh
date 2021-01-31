#!/bin/bash
#*********************************************************************#
#Data Preprocessing Script by Lindsey Fenderson                       #
#Current version 1.0 - January 30, 2021                               #
#Simple script for creating fastqc reports on newly sequenced data,   #
#  (code for double checking barcodes if needed),                     #
#  & trimming adapters with AdapterRemoval2 on UNH's Premise cluster. #
#*********************************************************************#

#Process assumes that sequencing data for the same individual generated on multiple lanes in the same sequencing run has already been concatenated and that the raw sequencing files have been named informatively according to the Kovach_Lab_File_Naming_SOPs.
module purge
module load linuxbrew/colsa

##Generate FastQC reports for your sequencing data##
mkdir ./fastqcRawData
#Create a list of files to be processed
ls *fastq.gz > FilesToProcess
#Run fastqc
while read FilesToProcess; do
     fastqc -o fastqcRawData -t 24 -f fastq $FilesToProcess
done< FilesToProcess

#Create a root name list
cut -d"_" -f1-4 FilesToProcess > RootName
sort RootName | uniq > UniqRootName

Suffix1='_R1.fastq.gz'
Suffix2='_R2.fastq.gz'

while read UniqRootName; do
     File1="$UniqRootName$Suffix1"
     File2="$UniqRootName$Suffix2"
     #Optional - Double check the adapters before trimming
     #If needed, you can confirm the adaptors to search for and remove with the following command (I haven't integrated it into the script yet, but you can just run manually if you are unsure of the adapter sequence used; then just replace the consensus output with the adapters used below in the next step):
     ####AdapterRemoval --identify-adapters --file1 $File1 --file2 $File2
     #FYI the adapters used below are standard for Illumina Truseq dual-indexed libraries (adapter 1 written as 5'-3'; adapter 2 always given as the reverse complement sequence), containing the Illumina Universal adapter sequencing primer binding site sequence (read 1=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC; read 2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT), followed by the sample's unique 8-mer index, and the i7/i5 primer sequence in every Illumina adapter that binds to the flow cell (adapter 1: CGTATGCCGTCTTCTGCTTG; adapter 2: TCGGTGGTCGCCGTATCATT). Note that this also removes the read 1 3' A overhangs added during library prep (hence why there is an 'A' at the start of the read 1 adapter sequence before the primer binding site sequence). One of the other most common Illumina Truseq adapter combinations used is just a single index adapter, in which case just delete 'NNNNNNNNGTGT' from the middle of adapter 2; also note that the length of the indexing barcode has changed over time so your adapter 1 may have fewer N's (e.g., 6 (rather than the 8 used here) is also common).
     
#Trim Adapters. Parameters set to use 24 threads, the fastq read number in the input files is separated from the header by a space, trim consecutive Ns from the 5' and 3' termini, and trim consecutive stretches of low quality bases at the 5' and 3' termini (below a minquality of 25), discard reads that have 150 or more Ns (with 250bp reads, this would require a minimum of about 25-30 bp of data after trimming the adapters),  and require reads to be at least 30 bp after trimming. The adapter sequences to be trimmed from the 3' ends of the reads are given, allowing ambiguity for the individual sample index so the code will run across samples. Read pairs that overlap by a minimum of 11 bases (the default) are collapsed into a single contig. Output files are combined, meaning all of the trimmed read1 and read2 reads as well as collapsed read1/read2 reads (added to the read1 file) and singletons where one of the mate pairs was discarded (with reads discarded due to quality filters or read merging are replaced with a single 'N' and Phred score of 0). These output files are  bzip2 compressed at the max compression level (9). Remaining parameters were left as default (i.e., --mm mismatch rate of 3, --shift of 2 bp in read2 to allow for missing bases, essentially no limit on max read length, and the default --minalignmentlength of 11bp overlap required before collapsing read1 and read2 pairs.)
     AdapterRemoval --file1 $File1 --file2 $File2 --basename $UniqRootName --threads 24 --mate-separator ' ' --trimns --trimqualities --minquality 25 --maxns 150 --minlength 30 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --collapse --combined-output --bzip2
     ##bzip2 $UniqRootName.pair1.truncated $UniqRootName.pair2.truncated $UniqRootName.singleton.truncated $UniqRootName.discarded
done< UniqRootName

#Recommended next step is to start mapping pipeline
#see BWA-MEM_Alignment.sh 
