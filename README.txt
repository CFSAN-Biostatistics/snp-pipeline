Set up to run the snp pipeline - an example:

    1. Get the Reference sequence
        NC_011149

    2. Get the sample sequences
        fastq-dump --split-files ERR178930

    0. Create index file for reference
        bowtie-build   NC_011149

    1. Align sequences to reference (bowtie2-align or just bowtie2 wrapper?)
        bowtie2-align -p 11 -q -x /Users/james.pettengill/Downloads/bowtie2-2.2.0/example/reference/lambda_virus -1 reads4_1.fq -2 reads4_2.fq > reads4.sam

    2. Convert to bam file with only mapped positions
        samtools view -bS -F 4 -o /Users/james.pettengill/Downloads/bowtie2-2.2.0/example/reads/reads1_F4.bam reads1.sam

    3. Convert to a sorted bam 
        samtools sort  reads4_F4.bam reads4_F4.sorted.bam

    4. Get a bcf file from the pileup and bam file
        samtools mpileup -uf /Users/james.pettengill/Downloads/bowtie2-2.2.0/example/reference/lambda_virus.fa reads1_F4.sorted.bam.bam | bcftools view -bvcg - > reads1_F4.bcf

    5. Convert bcf to vcf
        bcftools view reads1_F4.bcf | vcfutils.pl varFilter -D1000 > var1_F4.flt.vcf    

    6. Run samtools pileup in parallel and combine alignment and pileup to generate snp matrix
        ./snppipeline.py -n 10 -d ~/projects/snppipeline/test/testLambdaVirus/ -f path.txt -r lambda_virus.fa -l snplist.txt -a snpma.fasta -i True

