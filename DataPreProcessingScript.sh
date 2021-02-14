#!/bin/bash
#*********************************************************************#
#Data Preprocessing Script by Lindsey Fenderson                       #
#Current version 1.1.1 - February 14, 2021                            #
#Script for creating fastqc reports on newly sequenced data,          #
#  (code for double checking barcodes if needed),                     #
#  & trimming adapters with AdapterRemoval2 on UNH's Premise cluster. #
#*********************************************************************#

#Process assumes that sequencing data for the same individual generated on multiple lanes in the same sequencing run has already been concatenated and that the raw sequencing files have been named informatively according to the Kovach_Lab_File_Naming_SOPs.

#To use this program:
#   1.) Program assumes the raw reads are in fastq.gz format. If not, change the suffix as needed in lines 35, 36, 46, 47, 108, and 133; and possibly lines 48 and 49 if the fastqc output file format differs.
#   2.) If desired, change the number of max Ns permitted before a read is discarded depending on the sequencing read length of your data (e.g., I used a max of 150Ns for 250bp reads, or no more than 60% of the read could be Ns) in line 70. Confirm other parameter choices.
#   3.) Edit "Reports" and "Stats" and/or add the directory names or add any other files as needed for the raw sequencing metadata to be archived with the raw reads as relevant for your data (e.g, "report", "statsPlate1", "MD5.txt", "Summary.txt", "checkSize.xls", "Rawdata_Readme.pdf" etc.) in lines 110-111(+). Be sure to include the absolute or relative path if the files are not in the working directory.

##Expected output of this program: 
# 1.) This program will generate fastqc analysis reports on the individual R1 and R2 raw sequencing data and output them to the folder 'fastqcRawData'. For whole genome shotgun data, at this stage it is expected the fastq files will fail the AdapterContent (others?) module. (Need to list fastqc parameters used). If any of the files fail the # modules, it is suggestive of machine failure and to contact the sequencing provider to have the samples sequenced again.
# 2.) The program then quality trims adapters and low quality termini from reads and collapses overlapping read pairs using AdapterRemoval, with the parameters used described below in line 69. Output files will have the same root name as the input fastq files, with the suffix .pair#.truncated.bz2 (where # is the number 1 or 2 for reads 1 or 2 respectively).
# 3.) The program will then generate fastqc analysis reports on the individual R1 and R2 quality trimmed sequencing data generated in step 2 and output them to the folder 'fastqcTrimmedData'. At this stage it is expected the fastq files will now pass the AdapterContent (others?) module. It is also expected that the sequence length distribution module may be flagged, as there will be more variability in sequence length due to the collapsed reads. (Need to list fastqc parameters used)
# 4.) Finally, the program performs some data and file cleanup and compression to best preserve disc space; archiving all of the raw reads along with any sequencer metadata and the adapter removal settings file in $Directory_RawSequencingReads.tar.bz2, and archiving all of the fastqc output for the raw and trimmed reads in $Directory_FastQCReports.tar.bz2. These files, along with all of the pdf Adapter Removal Settings files and fastqc reports are organized into the directory "RawReads-and-SequencingPreProcessingData". The script also records the fastqc and AdapterRemoval program version numbers that were used near the top of the slurm output file.
#   ***Note 1: As per our Kovach_Lab_Data_Management_User_Guide-START_HERE SOP, it is STRONGLY recommended (required?!) that the raw sequencing reads are backed up in 2 separate locations (both in the cloud and on a physical hard drive located on campus) and that the integrity of the backups is confirmed PRIOR to running this script, since the combined raw read file archive may be too large to effectively upload or download to other backup locations (or it could get corrupted in the process) and the individual raw read files will be deleted once they are added to the archive (so you should still be able to access the raw reads if you decompress and unarchive the $Directory_RawSequencingReads.tar.bz2 file, however there is always a risk of data loss if the file gets corrupted). 

module purge
module load linuxbrew/colsa

#Record program versions:
echo -e "Program versions used in this script: \n"
fastqc -v
AdapterRemoval --version

##Generate FastQC reports for your sequencing data##
mkdir ./fastqcRawData
#Create a list of files to be processed (DOES process Undetermined raw fastq files for fastqc)
ls --ignore=Undetermined* | grep .fastq.gz > FilesToProcess
ls *.fastq.gz > AllFilesToProcess
#Run fastqc
while read AllFilesToProcess; do
     fastqc -o fastqcRawData -t 24 -f fastq $AllFilesToProcess
done< AllFilesToProcess

#Create a root name list
cut -d"_" -f1-4 AllFilesToProcess > RootName
sort RootName | uniq > UniqRootName

Suffix1='_R1.fastq.gz'
Suffix2='_R2.fastq.gz'
Suffix1A='_R1_fastqc.html'
Suffix2A='_R2_fastqc.html'
Suffix1B='_R1_fastqc.pdf'
Suffix2B='_R2_fastqc.pdf'

while read UniqRootName; do
     #Print all of the fastqc reports to pdf
     File1A="$UniqRootName$Suffix1A"
     File2A="$UniqRootName$Suffix2A"
     File1B="$UniqRootName$Suffix1B"
     File2B="$UniqRootName$Suffix2B"
     pandoc fastqcRawData/$File1A -V geometry:landscape -t latex -o $File1B
     pandoc fastqcRawData/$File2A -V geometry:landscape -t latex -o $File2B
done< UniqRootName

#Recreate root name list without the Undetermined files
cut -d"_" -f1-4 FilesToProcess > RootName
sort RootName | uniq > UniqRootName

#Run AdapterRemoval2 (does not process Undetermined fastq files)
while read UniqRootName; do
     File1="$UniqRootName$Suffix1"
     File2="$UniqRootName$Suffix2"
     #Optional - Double check the adapters before trimming
     #If needed, you can confirm the adaptors to search for and remove with the following command (I haven't integrated it into the script yet, but you can just run manually if you are unsure of the adapter sequence used; then just replace the consensus output with the adapters used below in the next step):
     ####AdapterRemoval --identify-adapters --file1 $File1 --file2 $File2
     #FYI the adapters used below are standard for Illumina Truseq dual-indexed libraries (adapter 1 written as 5'-3'; adapter 2 always given as the reverse complement sequence), containing the Illumina Universal adapter sequencing primer binding site sequence (read 1=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC; read 2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT), followed by the sample's unique 8-mer index, and the i7/i5 primer sequence in every Illumina adapter that binds to the flow cell (adapter 1: CGTATGCCGTCTTCTGCTTG; adapter 2: TCGGTGGTCGCCGTATCATT). Note that this also removes the read 1 3' A overhangs added during library prep (hence why there is an 'A' at the start of the read 1 adapter sequence before the primer binding site sequence). One of the other most common Illumina Truseq adapter combinations used is just a single index adapter, in which case just delete 'NNNNNNNNGTGT' from the middle of adapter 2; also note that the length of the indexing barcode has changed over time so your adapter 1 may have fewer N's (e.g., 6 (rather than the 8 used here) is also common).
     
#Trim Adapters. Parameters set to use 24 threads, the fastq read number in the input files is separated from the header by a space, trim consecutive Ns from the 5' and 3' termini, and trim consecutive stretches of low quality bases at the 5' and 3' termini (below a minquality of 25), discard reads that have 150 or more Ns (with 250bp reads, this would require a minimum of 40% of the read to be called),  and require reads to be at least 30 bp after trimming. The adapter sequences to be trimmed from the 3' ends of the reads are given, allowing ambiguity for the individual sample index so the code will run across samples. Read pairs that overlap by a minimum of 11 bases (the default) are collapsed into a single contig. Output files are combined, meaning all of the trimmed read1 and read2 reads as well as collapsed read1/read2 reads (added to the read1 file) and singletons where one of the mate pairs was discarded (with reads discarded due to quality filters or read merging are replaced with a single 'N' and Phred score of 0). These output files are  bzip2 compressed at the max compression level (9). Remaining parameters were left as default (i.e., --mm mismatch rate of 3, --shift of 2 bp in read2 to allow for missing bases, essentially no limit on max read length, and the default --minalignmentlength of 11bp overlap required before collapsing read1 and read2 pairs.)
     AdapterRemoval --file1 $File1 --file2 $File2 --basename $UniqRootName --threads 24 --mate-separator ' ' --trimns --trimqualities --minquality 25 --maxns 150 --minlength 30 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --collapse --combined-output --bzip2
done< UniqRootName

##Generate FastQC reports for your trimmed sequencing data##
mkdir fastqcTrimmedData

#Create a list of files to be processed (does not process Undetermined fastq files)
ls --ignore=Undetermined* | grep .truncated.bz2 > TrimmedFiles

while read TrimmedFiles; do
     fastqc -o fastqcTrimmedData -t 24 -f fastq $TrimmedFiles
done< TrimmedFiles

#Print the fastqc reports & adapter removal settings to pdf
while read UniqRootName; do
    Suffix1C='.pair1.truncated_fastqc.html'
    Suffix2C='.pair2.truncated_fastqc.html'
    Suffix1D='.pair1.truncated_fastqc.pdf'
    Suffix2D='.pair2.truncated_fastqc.pdf'
    Suffix1E='.settings'
    Suffix1F='.AdapterRemovalSettings.pdf'
    File1C="$UniqRootName$Suffix1C"
     File2C="$UniqRootName$Suffix2C"
     File1D="$UniqRootName$Suffix1D"
     File2D="$UniqRootName$Suffix2D"
     File1E="$UniqRootName$Suffix1E"
     File1F="$UniqRootName$Suffix1F"
     pandoc fastqcTrimmedData/$File1C -V geometry:landscape -t latex -o $File1D
     pandoc fastqcTrimmedData/$File2C -V geometry:landscape -t latex -o $File2D
     pandoc $File1E -t latex -o $File1F
done< UniqRootName

#Archive the Adapter Removal settings files
mkdir AdapterRemovalSettings
rsync -acv *.settings AdapterRemovalSettings

#Archive and compress the raw reads to organize the directory and save disc space
#Create a list of files to be processed (DOES include Undetermined fastq files)
ls *.fastq.gz > AllFilesToProcess
echo Reports >> AllFilesToProcess
echo Stats >> AllFilesToProcess
echo AdapterRemovalSettings >> AllFilesToProcess

#Automatically get directory name for the output archive file names
Directory=$(pwd | rev | cut -d'/' -f 1 | rev)

#Archive and compress the raw sequencing reads and reports
tar cvfj $Directory_RawSequencingReads.tar.bz2 -T AllFilesToProcess

#Archive and compress the fastqc data to organize directory and save disc space (pdf reports still available in working directory for review and easy addition to an electronic lab notebook)
rsync -acv *_R1_fastqc.pdf fastqcRawData/
rsync -acv *_R2_fastqc.pdf fastqcRawData/
rsync -acv *.truncated_fastqc.pdf fastqcTrimmedData/
tar cvfj $Directory_FastQCReports.tar.bz2 fastqcRawData fastqcTrimmedData

#Cleanup tmp files and archived files
rm FilesToProcess
rm RootName
rm UniqRootName
rm AllFilesToProcess
rm TrimmedFiles
rm -r Reports
rm -r Stats
rm *.fastq.gz
rm *.settings
rm -r fastqcRawData
rm -r fastqcTrimmedData
rm -r AdapterRemovalSettings
mkdir FastQCReports
mkdir FastQCReports/fastqcRawData
mkdir FastQCReports/fastqcTrimmedData
rsync -acv *_R1_fastqc.pdf FastQCReports/fastqcRawData/
rsync -acv *_R2_fastqc.pdf FastQCReports/fastqcRawData/
rsync -acv *.truncated_fastqc.pdf FastQCReports/fastqcTrimmedData/
rm *_R1_fastqc.pdf
rm *_R2_fastqc.pdf
rm *.truncated_fastqc.pdf
mkdir AdapterRemovalSettings
rsync -acv *.AdapterRemovalSettings.pdf AdapterRemovalSettings/
rm *.AdapterRemovalSettings.pdf
mkdir RawReads-and-SequencingPreProcessingData
rsync -acv *_FastQCReports.tar.bz2 RawReads-and-SequencingPreProcessingData/
rsync -acv *_RawSequencingReads.tar.bz2 RawReads-and-SequencingPreProcessingData/
rsync -acv AdapterRemovalSettings RawReads-and-SequencingPreProcessingData/
rsync -acv FastQCReports RawReads-and-SequencingPreProcessingData/
rm *_FastQCReports.tar.bz2
rm *_RawSequencingReads.tar.bz2
rm -r AdapterRemovalSettings
rm -r FastQCReports
#Recommended next step is to start mapping pipeline
#see BWA-MEM_Alignment.sh 
