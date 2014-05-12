#Example: Set up to run the snp pipeline and run it.

#Set up reference sequence
/home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/prepReference NC_011149

#Set up sample sequence
#    Note: This could be run in parallel (thank Jamie) using gnu parallel
echo -e "ERR178926\nERR178927\nERR178928\nERR178929\nERR178930" > prepInput
cat prepInput | xargs -n 1 /home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/prep.sh NC_011149
        
# Run samtools pileup in parallel and combine alignment and pileup to generate snp matrix
./snppipeline.py -n 10 -d ~/projects/snppipeline/test/testLambdaVirus/ -f path.txt -r lambda_virus.fa -l snplist.txt -a snpma.fasta -i True

-----------------OLD STUFF FOLLOWS----------------------------
Set up to run the snp pipeline - an example:

    0. Set up some directories
        mkdir reference samples

    1. Get the Reference sequence
        ~/mnt/biob/svn/Biostats/rand/cfsanutils/scripts/fetch.py NC_011149 -e hugh.rand@fda.hhs.gov > reference/NC_011149.fasta

    2. Get the sample sequences
        fastq-dump --outdir samples/ERR178926 --split-files ERR178926
        fastq-dump --outdir samples/ERR178927 --split-files ERR178927
        fastq-dump --outdir samples/ERR178928 --split-files ERR178928
        fastq-dump --outdir samples/ERR178929 --split-files ERR178929
        fastq-dump --outdir samples/ERR178930 --split-files ERR178930

    3. Create index file for reference
        ~/software/bowtie2-2.2.2/bowtie2-build reference/NC_011149.fasta reference/NC_011149


    4. Align sequences to reference
        ~/software/bowtie2-2.2.2/bowtie2 -p 11 -q -x reference/NC_011149 -1 samples/ERR178926/ERR178926_1.fastq -2 samples/ERR178926/ERR178926_2.fastq > samples/ERR178926/ERR178926.sam
        ~/software/bowtie2-2.2.2/bowtie2 -p 11 -q -x reference/NC_011149 -1 samples/ERR178927/ERR178927_1.fastq -2 samples/ERR178927/ERR178927_2.fastq > samples/ERR178927/ERR178927.sam
        ~/software/bowtie2-2.2.2/bowtie2 -p 11 -q -x reference/NC_011149 -1 samples/ERR178928/ERR178928_1.fastq -2 samples/ERR178928/ERR178928_2.fastq > samples/ERR178928/ERR178928.sam
        ~/software/bowtie2-2.2.2/bowtie2 -p 11 -q -x reference/NC_011149 -1 samples/ERR178929/ERR178929_1.fastq -2 samples/ERR178929/ERR178929_2.fastq > samples/ERR178929/ERR178929.sam
        ~/software/bowtie2-2.2.2/bowtie2 -p 11 -q -x reference/NC_011149 -1 samples/ERR178930/ERR178930_1.fastq -2 samples/ERR178930/ERR178930_2.fastq > samples/ERR178930/ERR178930.sam

    5. Convert to bam file with only mapped positions
        samtools view -bS -F 4 -o samples/ERR178926/ERR178926.bam samples/ERR178926/ERR178926.sam
        samtools view -bS -F 4 -o samples/ERR178927/ERR178927.bam samples/ERR178927/ERR178927.sam
        samtools view -bS -F 4 -o samples/ERR178928/ERR178928.bam samples/ERR178928/ERR178928.sam
        samtools view -bS -F 4 -o samples/ERR178929/ERR178929.bam samples/ERR178929/ERR178929.sam
        samtools view -bS -F 4 -o samples/ERR178930/ERR178930.bam samples/ERR178930/ERR178930.sam

    6. Convert to a sorted bam 
        samtools sort samples/ERR178926/ERR178926.bam samples/ERR178926/ERR178926.sorted.bam
        samtools sort samples/ERR178927/ERR178927.bam samples/ERR178927/ERR178927.sorted.bam
        samtools sort samples/ERR178928/ERR178928.bam samples/ERR178928/ERR178928.sorted.bam
        samtools sort samples/ERR178929/ERR178929.bam samples/ERR178929/ERR178929.sorted.bam
        samtools sort samples/ERR178930/ERR178930.bam samples/ERR178930/ERR178930.sorted.bam

    7. Get a bcf file from the pileup and bam file
        samtools mpileup -uf reference/NC_011149.fasta samples/ERR178926/ERR178926.sorted.bam.bam | bcftools view -bvcg - > samples/ERR178926/ERR178926.bcf
        samtools mpileup -uf reference/NC_011149.fasta samples/ERR178927/ERR178927.sorted.bam.bam | bcftools view -bvcg - > samples/ERR178927/ERR178927.bcf
        samtools mpileup -uf reference/NC_011149.fasta samples/ERR178928/ERR178928.sorted.bam.bam | bcftools view -bvcg - > samples/ERR178928/ERR178928.bcf
        samtools mpileup -uf reference/NC_011149.fasta samples/ERR178929/ERR178929.sorted.bam.bam | bcftools view -bvcg - > samples/ERR178929/ERR178929.bcf
        samtools mpileup -uf reference/NC_011149.fasta samples/ERR178930/ERR178930.sorted.bam.bam | bcftools view -bvcg - > samples/ERR178930/ERR178930.bcf

    8. Convert bcf to vcf
        bcftools view samples/ERR178926/ERR178926.bcf | vcfutils.pl varFilter -D1000 > samples/ERR178926/ERR178926.vcf
        bcftools view samples/ERR178927/ERR178927.bcf | vcfutils.pl varFilter -D1000 > samples/ERR178927/ERR178927.vcf
        bcftools view samples/ERR178928/ERR178928.bcf | vcfutils.pl varFilter -D1000 > samples/ERR178928/ERR178928.vcf
        bcftools view samples/ERR178929/ERR178929.bcf | vcfutils.pl varFilter -D1000 > samples/ERR178929/ERR178929.vcf
        bcftools view samples/ERR178930/ERR178930.bcf | vcfutils.pl varFilter -D1000 > samples/ERR178930/ERR178930.vcf

    9. Run samtools pileup in parallel and combine alignment and pileup to generate snp matrix
        ./snppipeline.py -n 10 -d ~/projects/snppipeline/test/testLambdaVirus/ -f path.txt -r lambda_virus.fa -l snplist.txt -a snpma.fasta -i True

