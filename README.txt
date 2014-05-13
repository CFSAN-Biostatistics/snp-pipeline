#Example: Set up to run the snp pipeline and run it.

#Prep work
mkdir #make a directory to work in
cd    #change to it

#Set up reference sequence
/home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepReference.sh NC_011149

#Set up sample sequence
#  Note: This could be run in parallel using gnu parallel (workstation) or qsub (HPC)
echo -e "ERR178926\nERR178927\nERR178928\nERR178929\nERR178930\n" > prepInput
cat prepInput | xargs -n 1 /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepSamples.sh NC_011149
        
# Run samtools pileup in parallel and combine alignment and pileup to generate snp matrix
ls -d -1 $PWD/samples/* > path.txt
./snppipeline.py -n 10 -d CHANGETHIS -f path.txt -r reference/NC_011149.fasta -l snplist.txt -a snpma.fasta -i True
#./snppipeline.py -n 10 -d ~/projects/snppipeline/test/testLambdaVirus/ -f path.txt -r lambda_virus.fa -l snplist.txt -a snpma.fasta -i True

#Prep test of running in parallel with qsub
qsub /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepSamples.sh NC_011149 ERR178926
qsub /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepSamples.sh NC_011149 ERR178927
qsub /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepSamples.sh NC_011149 ERR178928
qsub /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepSamples.sh NC_011149 ERR178929
qsub /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/scripts/prepSamples.sh NC_011149 ERR178930
