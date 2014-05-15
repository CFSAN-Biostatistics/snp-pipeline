#Example: Set up to run the snp pipeline and run it.

#Prep work
mkdir junk #make a directory to work in
cd junk    #change to it

#Set up reference sequence
/home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepReference.sh NC_011149

#Set up sample sequence
#  Note: This could be run in parallel using gnu parallel (workstation) or qsub (HPC)
echo -e "ERR178926\nERR178927\nERR178928\nERR178929\nERR178930\n" > prepInput
cat prepInput | xargs -n 1 /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepSamples.sh NC_011149
        
# Run samtools pileup in parallel and combine alignment and pileup to generate snp matrix
ls -d -1 --color=never $PWD/samples/* > path.txt
/home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/runsnppipeline.py -n 10 -d ~/mnt/biob/svn/Biostats/rand/snppipeline/test/junk/ -f path.txt -r NC_011149.fasta -l snplist.txt -a snpma.fasta -i True
/home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/runsnppipeline.py -n 10 -d ~/mnt/biob/svn/Biostats/rand/snppipeline/test/junk/ -f path.txt -r reference/NC_011149.fasta -l snplist.txt -a snpma.fasta -i True
#./snppipeline.py -n 10 -d ~/projects/snppipeline/test/testLambdaVirus/ -f path.txt -r lambda_virus.fa -l snplist.txt -a snpma.fasta -i True

