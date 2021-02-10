#!/bin/bash
#************************************************************************#
#BWA-MEM_Alignment Script by Lindsey Fenderson                           #
#Current version 1.0.0 - February 10, 2021                               #
#Pipeline for mapping fastq reads to a reference genome with the BWA-MEM #               #
#  algorithm and additional quality control steps based on the GATK Best # 
#  Practices Data Pre-Processing for Variant Discovery Pipeline          #
#  (including samtools fixmate, index & sort; GATK Indel Realigner,      #
#  Picard Mark Duplicates, and GATK Base Quality Score Recalibration)    #
#  to generate analysis-ready bam files.                                 #
#************************************************************************#

#Process assumes that input data have already been processed (e.g., with the DataPreProcessingScript.sh script) to evaluate the sequencing output and trim adapters and low quality bases or reads.
#To use the current version of this program, do the following: 
#   1.) Set the data path of the files you wish to map in line 24. 
#       1A.) If samples were not preprocessed with the DataPreProcessingScript.sh, change the suffix of the batch of files to be mapped as needed in lines 31, 38 & 39.
#   2.) Change the path and reference genome name if required to the species you wish to map the reads to in line 49.
#   3.) Change the output path if required to the corresponding directory of the species the reads have been mapped to in line 25. 
#   4.) Change the reference genome suffix that will be appended to the output filenames to the genome version the reads were mapped to if required in line 26.

module purge
module load linuxbrew/colsa

InputDataPath="/mnt/lustre/mel/shared/GECO/Data/SparrowRawData/20200728_INVS-SP_LFe_SparrowWholeGenomeShotgun/"
OutputDataPath="/mnt/lustre/mel/shared/GECO/Data/Mapped_Data/Gallusgallus_GRCg6a/"
ReferenceGenomeSuffix="GRCg6a"
cd $OutputDataPath

#Create a list of files to be processed

ls $InputDataPath*truncated.bz2 > FilesToProcess

#Create a root name list
cut -d"/" -f10 FilesToProcess > FilesToProcess2
cut -d"_" -f1-4 FilesToProcess2 > RootName
cut -d"." -f1 RootName > RootName2
sort RootName2 | uniq > UniqRootName
Suffix1='.pair1.truncated.bz2'
Suffix2='.pair2.truncated.bz2'

#Clean up temp files
rm FilesToProcess2
rm RootName*

while read UniqRootName; do
     File1="$UniqRootName$Suffix1"
     File2="$UniqRootName$Suffix2"
#Map input files to specified reference. Parameters include using 24 threads, smart pairing (If two adjacent reads have the same name, they are considered to form a read pair. This way, paired-end and single-end reads can be mixed in a single FASTA/Q stream), the read group header line for the output SAM file (includes the sequencing instrument ID:run number on the instrument:flow cell ID from the header of the fastq files as the unique read group identifier; the root sample name of the sample that was sequenced; an identifier signifying the library batch the sequence was derived from; and the sequencing platform that generated this data), including alignments for single-end/unpaired paired-end reads so as not to discard any data, (but they will be flagged as secondary alignments), (want to use flag -C here to appending the barcode to the SAM file but is formatted incorrectly so as to lead to incorrect SAM output; need to figure out how to fix and e.g., make header in sam format with cs:Z: prefix; have removed flag for now), and mark shorter split hits as secondary for Picard tools compatibility. Then the prefix of the indexed genome is specified, followed by the R1 and R2 filenames to be mapped (specifying that they should be decompressed from bzip2 format to the standard output so the reads can be mapped by bwa) and the desired output filename. 

#Parameters left as default (mostly because I don't have any understanding of if it is worthwhile or better to tune any of these parameters, and Heng Li is a pretty smart dude!) include: the minimum seed length of 19, band with of 100,Off-diagonal dropoff of 100, re-seeding value of 1.5, discard an alignment if it occurs more than 10000 times in the genome, did not use -P flag to rescue missing hits/ignore hits that fit a proper pair, default matching score of 1, default mismatch penalty of 4, default Gap open penalty of 6, default gap extension penalty of 1, clipping penalty of 5, unpaired read pair penalty of 9, did not use -p flag as input files are not interleaved, not including alignment in output if it has a mapping quality score below 30, leaving soft clipping in place, verbosity is what it is!
     echo $File1
     Read1="$InputDataPath$File1"
     echo $Read1
     Read2="$InputDataPath$File2"
     echo $Read2
     bwa mem -t $SLURM_CPUS_PER_TASK -p -R "@RG\tID:A00873:204:HMV7NDRXX\tSM:${UniqRootName}\tLB:DiluteKappaJuly2020\tPL:ILLUMINA" -a -M /mnt/lustre/mel/shared/GECO/ReferenceGenomes/DomesticChicken/Ggallus_GRCg6a <(bzip2 -dc $Read1) <(bzip2 -dc $Read2) > "${UniqRootName}_trimmed-bwamem-GRCg6a.sam"
     #"Because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags" 
     #This command runs samtools fixmate, which fills in mate coordinates, ISIZE and mate related flags from a name-sorted alignment, and also sets the output format to bam. I opted to just use the defaults and did not include the flags to Remove secondary and unmapped reads (-r), Disable FR proper pair check (-p), Add template cigar ct tag (-c), Add ms (mate score) tags (-m) (These are used by markdup to select the best reads to keep), or Do not add a @PG line to the header of the output file (--no-PG). I suppose I would use the mate score tags if I used the samtools markdup tool instead of Picardtools MarkDuplicates...
     samtools fixmate -O bam ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sam ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
     #Most subsequent tools require coordinate-sorted (instead of read name-sorted) bam files and bam files need to be indexed to facilitate their use with various tools, hence the next 2 commands.
     samtools sort -O bam -o ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
     samtools index ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
     #Cleanup temp files
     rm ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sam
     rm ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
     ##bzip2 "${UniqRootName}_trimmed-bwamem-GRCg6a.sam"
     #bzip2 -dc 2631-22921_Acaudacuta_MillGutRI_20190531_trimmed-bwamem-GRCg6a.sam.bz2 | samtools fixmate -O bam - 2631-22921_Acaudacuta_MillGutRI_20190531_trimmed-bwamem-GRCg6a.bam
done< UniqRootName

#Run GATK Indel Realigner
#Note that this tool requires an indexed reference genome (it it doesn't already exist, this can be generated in the linuxbrew/colsa module via e.g.,: samtools faidx GCF_000002315.6_GRCg6a_genomic.fna) and a reference sequence dictionary (if it doesn't already exist, this can be generated in the linuxbrew/colsa module via e.g.,: gatk CreateSequenceDictionary -R GCF_000002315.6_GRCg6a_genomic.fna)
#Indel realigner is not supported in GATK4, as from what I can gather, the tool has been incorporated into the GATK vcf callers. I decided to include this step regardless to better permit the use of alternate variant discovery tools that don't necessarily incorporate realignment. Thus, to utilize the tool I need to switch to a conda environment so can use GATK3.8.
module purge
module load anaconda/colsa

#Use the RealignerTargetCreator tool to generate a set of sites likely to contain indels
while read UniqRootName; do
java -jar /mnt/lustre/mel/shared/Software/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /mnt/lustre/mel/shared/GECO/ReferenceGenomes/DomesticChicken/GCF_000002315.6_GRCg6a_genomic.fna \
    -I ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam \
    -o ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.intervals
#Perform local realignment around indels
java -Xmx8G -Djava.io.tmpdir=/tmp -jar /mnt/lustre/mel/shared/Software/GenomeAnalysisTK.jar \
    -T IndelRealigner --filter_bases_not_stored \
    -R /mnt/lustre/mel/shared/GECO/ReferenceGenomes/DomesticChicken/GCF_000002315.6_GRCg6a_genomic.fna \
    -targetIntervals ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.intervals \
    -I ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam \
    -o ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.indelrealigned.bam
done< UniqRootName

#Mark PCR duplicates
module purge
module load linuxbrew/colsa
#***Note that using Picardtools probably under represents the number of PCR duplicates since I am using collapsed reads. In the future should probably incorporate the paleomix rmdup_collapsed tool (This tool "Filters PCR duplicates for merged/collapsed paired-ended reads, such as those generated by AdapterRemoval with the --collapse option enabled. Unlike SAMtools rmdup or Picard MarkDuplicates, this tool identifies duplicates based on both the 5' and the 3' alignment coordinates of individual reads.") - paleomix rmdup_collapsed --remove-duplicates < sorted.bam > < out.bam >
while read UniqRootName; do
gatk MarkDuplicates -I ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.indelrealigned.bam -O ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.marked.indelrealigned.bam -M ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.indelrealigned.dup_metrics.txt
#Resort the files
#I'm unclear as to why this needs to be done again, but various pipelines I've seen online do it, so I'm following the herd off the cliff...can't hurt at any rate!
samtools sort -O bam -o ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.resorted.marked.indelrealigned.bam ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.marked.indelrealigned.bam
done< UniqRootName

#At this stage it is ideally recommended to rescale the base quality scores. This requires a known variants file. I used a chicken vcf from http://bigd.big.ac.cn/chickensd/, which was developed from 863 chicken genomes (reported in Wang et al. (2020) 863 genomes reveal the origin and domestication of chicken. Cell Research), as it makes more sense to me to use a chicken reference in this case since all 6 species were mapped to the chicken reference genome, than to use the only sparrow vcfs known to me (from Jen's 2019 EcoEvo paper for STSP only and 1 SESP, and from Jen's 2019 EvolLett paper for SAVS, NESP, SOSP and SWSP). As I understand it, this step can also be omitted.


#Cleanup temp files
rm *.intervals
rm *.sorted.bam

#Also need to fix the read group/library stuff so it extracts the info from a Yaml - TODO  "Typically your reads will be supplied to you in two files written in the FASTQ format. It is particularly important to ensure that the @RG information here is correct as this information is used by later tools. The SM field must be set to the name of the sample being processed, and LB field to the library. The resulting mapped reads will be delivered to you in a mapping format known as SAM."

#bzip2 -dc 2631-22921_Acaudacuta_MillGutRI_20190531_trimmed-bwamem-GRCg6a.sam.bz2 | samtools fixmate -O bam - 2631-22921_Acaudacuta_MillGutRI_20190531_trimmed-bwamem-GRCg6a.bam

