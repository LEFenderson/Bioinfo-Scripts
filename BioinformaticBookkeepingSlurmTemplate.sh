#!/bin/bash
#***************************************************************#
#Example slurm template by Lindsey Fenderson - January 30, 2021 #
#Simple script to record important bioinformatic metadata       #
#  in slurm output generated by a slurm scheduler process,      #
#  and for storing the slurm output informatively               #
#  for future reference.                                        #
#***************************************************************#

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --job-name SeqPreProcess
#SBATCH --mail-type=ALL
#SBATCH --mail-user=User.Email@unh.edu

#To use this template, confirm or update the above sbatch parameters to be suitable for the script you are running and update the notification email address to your own. Be sure to update the working directory of the files that will be processed in line 23 and the directory name in line 44, and update the path/script name to be run in lines 30, 33 and 40, update the input/output directory paths and names if needed in lines 43-45 and update the name of the script used if needed in the output file name in line 44.

#Set date & time stamp variable
current_time=$(date "+%Y.%m.%d-%H.%M")

#Template assumes that the slurm job was scheduled with this template from a personal directory, while the files to be processed are stored in a shared directory
cd /absolute/path/to/working/directory/with/files/to/be/processed

#Record the date and time of the analysis, as well as the script name/date/version being run (assuming that metadata is appropriately listed in the header of the script being run), and the location of the files that are being processed. Note, if the stdout of your script doesn't list the files that are being manipulated, it is recommended to add a printout of the files being processed by the script to the slurm output, e.g.: echo -e "The list of input files being processed through this script is as follows: \n" ls *.fastq.gz

echo -e "The date and run time of this script was:" $current_time
echo
echo "The script being run, from the following location, and the script's last modified time stamp:"
ls -lh /absolute/path/to/Directory/of/script/being/run/shared/Scripts/DataPreProcessingScript.sh
echo
echo -e "The header of the script run was: \n"
head /absolute/path/to/Directory/of/script/being/run/shared/Scripts/DataPreProcessingScript.sh
echo
echo -e "The script was run on the files in the directory: \n"
pwd
echo

#Run the target script
/absolute/path/to/Directory/of/script/being/run/shared/Scripts/DataPreProcessingScript.sh

#After the analysis has run, return to your directory where you wish to store your bioinformatic bookkeeping records. Copy and informatively rename the slurm output from the directory where you initiated the script and delete the original. In this example, I have in my home directory a folder called "BioinformaticScripts" where I store all my slurm scheduler scripts and just change the path to the working directory of the files to be processed. I store all of my slurm output in my directory called "BioinformaticRecords". In this example, I have just run the DataPreProcessingScript.sh, (i.e., the programs fastqc and AdapterRemoval), so I rename the slurm output informatively with the script name that was run, the directory of files which were processed, and the date-time stamp of when I conducted the analysis in case I ever have to go back to the output, or if program versions change at some point I know exactly when this was run, etc.
cd ~/BioinformaticRecords
cp ~/BioinformaticScripts/slurm-*.out DataPreProcessing-20200728_INVS-SP_LFe_SparrowWholeGenomeShotgun-slurmOutput-$current_time.yaml
rm ~/BioinformaticScripts/slurm-*.out