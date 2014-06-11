#TODO DO NOT EXECUTE YET - NOT FINISHED
#!/bin/bash
#
#Author: Hugh A. Rand (har)
#Purpose: Run set up and test of snppipline code with S. Agona data.
#Input:
#    baseDirectory
#Output:
#    various
#Use example:
#    runAgonaVirusTest.sh /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline
#History:
#   20140611-har: Started.
#Notes:
#Bugs:
#

#Process arguments
if [ -z "$1" ]; then
    echo usage: $0 baseDirectory
    exit
fi
BASEDIRECTORY=$1

#Prep work
cd $BASEDIRECTORY'/test'
mkdir testAgona #directory to work in
cd testAgona

#Fetch sequence data
#TODO set it so it only fetched the reference sequence once!!
echo -e "ERR178926\nERR178927\nERR178928\nERR178929\nERR178930\n" > sampleDirectoryNames.txt
cat sampleDirectoryNames.txt | xargs -n 1 $BASEDIRECTORY'/scripts/prepSequenceData.sh' NC_011149

#Prep up reference sequence
$BASEDIRECTORY'/scripts/prepReference.sh' NC_011149

#Prep sample sequence
#  Note: This could be run in parallel using gnu parallel (workstation) or qsub (PBS on HPC)
cat sampleDirectoryNames.txt | xargs -n 1 $BASEDIRECTORY'/snppipeline/scripts/prepSamples.sh' NC_011149
        
#Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to generate snp matrix)
ls -d -1 --color=never $PWD/samples/* > path.txt
#$BASEDIRECTORY'/scripts/runsnppipeline.py' -n 10 -d $PWD -f path.txt -r reference/NC_011149.fasta -l snplist.txt -a snpma.fasta -i True
test/test_snppipeline.py -v

