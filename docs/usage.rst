.. _usage-label:

========
Usage
========

.. highlight:: bash

The SNP Pipeline is run from the Unix command line.  The pipeline consists of a collection
of shell scripts and python scripts.


+---------------------------+--------------------------------------------------------------------+
| Script                    | | Description                                                      |
+===========================+====================================================================+
| copy_snppipeline_data.py  | | Copies supplied example data to a work directory                 |
+---------------------------+--------------------------------------------------------------------+
| prepReference.sh          | | Indexes the reference genome                                     |
+---------------------------+--------------------------------------------------------------------+
| alignSampleToReference.sh | | Aligns samples to the reference genome                           |
+---------------------------+--------------------------------------------------------------------+
| prepSamples.sh            | | Finds variants in each sample                                    |
+---------------------------+--------------------------------------------------------------------+
| create_snp_list           | | Combines the SNP positions across all samples into a single      |
|                           | | unified SNP list file                                            |
+---------------------------+--------------------------------------------------------------------+
| create_snp_pileup         | | Creates the SNP pileup file for a sample -- the pileup file at   |
|                           | | the positions where SNPs were called in any of the samples       |
+---------------------------+--------------------------------------------------------------------+
| create_snp_matrix.py      | | Creates a matrix of SNPs across all samples                      |
+---------------------------+--------------------------------------------------------------------+
| create_snp_reference_seq  | | Writes the reference sequence bases at SNP locations to          |
|                           | | a fasta file                                                     |
+---------------------------+--------------------------------------------------------------------+


Step-by-Step Example Workflow Based on Lamda Virus Test Data Provided with Code
-------------------------------------------------------------------------------

Step 1 - Gather data::

    # The SNP Pipeline distribution includes sample data organized as shown below:
    snppipeline/data/lambdaVirusInputs/reference/lambda_virus.fasta
    snppipeline/data/lambdaVirusInputs/samples/sample1/sample1_1.fastq
    snppipeline/data/lambdaVirusInputs/samples/sample1/sample1_2.fastq
    snppipeline/data/lambdaVirusInputs/samples/sample2/sample2_1.fastq
    snppipeline/data/lambdaVirusInputs/samples/sample2/sample2_2.fastq
    snppipeline/data/lambdaVirusInputs/samples/sample3/sample3_1.fastq
    snppipeline/data/lambdaVirusInputs/samples/sample3/sample3_2.fastq
    snppipeline/data/lambdaVirusInputs/samples/sample4/sample4_1.fastq
    snppipeline/data/lambdaVirusInputs/samples/sample4/sample4_2.fastq

    # Copy the supplied test data to a work area:
    cd test
    copy_snppipeline_data.py lambdaVirusInputs testLambdaVirus
    cd testLambdaVirus

Step 2 - Prep work::

    # Create files of sample directories and fastQ files:
    ls -d --color=never samples/* > sampleDirectories.txt
    rm sampleFullPathNames.txt 2>/dev/null
    cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq >> sampleFullPathNames.txt; done
    # Determine the number of CPU cores in your computer
    NUMCORES=$(grep -c ^processor /proc/cpuinfo)

Step 3 - Prep the reference::

    prepReference.sh reference/lambda_virus.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    cat sampleFullPathNames.txt | xargs --max-args=2 --max-lines=1 alignSampleToReference.sh -p $NUMCORES reference/lambda_virus.fasta

Step 5 - Prep the samples::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES prepSamples.sh reference/lambda_virus.fasta

Step 6 - Combine the SNP positions across all samples into the SNP list file::

    create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt

Step 7 - Create pileups at SNP positions for each sample::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES -I XX create_snp_pileup.py -l snplist.txt -a XX/reads.all.pileup -o XX/reads.snp.pileup

Step 8 - Create the SNP matrix::

    create_snp_matrix.py -l snplist.txt -p reads.snp.pileup -o snpma.fasta sampleDirectories.txt

Step 9 - Create the reference base sequence::

    create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta reference/lambda_virus.fasta

        
Step 10 - View and verify the results:

Upon successful completion of the pipeline, the snplist.txt file should have 165 entries.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l snplist.txt
    ls -l snpma.fasta
    ls -l referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py lambdaVirusExpectedResults expectedResults
    diff -q -s snplist.txt         expectedResults/snplist.txt
    diff -q -s snpma.fasta         expectedResults/snpma.fasta
    diff -q -s referenceSNP.fasta  expectedResults/referenceSNP.fasta



Step-by-Step Example Workflow Using S. Agona Data Downloaded from SRA
---------------------------------------------------------------------

Step 1 - Gather data::

    # The SNP Pipeline distribution includes sample data organized as shown below:
    snppipeline/data/agonaInputs/sha256sumCheck
    snppipeline/data/agonaInputs/reference/NC_011149.fasta

    # Copy the supplied test data to a work area:
    cd test
    copy_snppipeline_data.py agonaInputs testAgona
    cd testAgona
    
    # Create sample directories
    mkdir -p samples/ERR178926  samples/ERR178927  samples/ERR178928  samples/ERR178929  samples/ERR178930
    
    # Download sample data from SRA at NCBI. Note that we use the fastq-dump command from
    #   the NCBI SRA-toolkit to fetch sample sequences. There are other ways to get the data,
    #   but the SRA-toolkit is easy to install, and does a good job of downloading large
    #   files.
    fastq-dump --split-files --outdir samples/ERR178926 ERR178926
    fastq-dump --split-files --outdir samples/ERR178927 ERR178927
    fastq-dump --split-files --outdir samples/ERR178928 ERR178928
    fastq-dump --split-files --outdir samples/ERR178929 ERR178929
    fastq-dump --split-files --outdir samples/ERR178930 ERR178930
    
    # Check the data
    #   The original data was used to generate a hash as follows:
    #     sha256sum reference/*.fasta samples/*/*.fastq > sha256sumCheck
    #   The command below checks the downloaded data (and the reference sequence) against the
    #     hashes that are saved in the sha256sumCheck file using sha256sum command, which is
    #     generally available on unix systems.
    sha256sum -c sha256sumCheck

Step 2 - Prep work::

    # Create files of sample directories and fastQ files:
    ls -d --color=never samples/* > sampleDirectories.txt
    rm sampleFullPathNames.txt 2>/dev/null
    cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq >> sampleFullPathNames.txt; done
    # Determine the number of CPU cores in your computer
    NUMCORES=$(grep -c ^processor /proc/cpuinfo)

Step 3 - Prep the reference::

    prepReference.sh reference/NC_011149.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    cat sampleFullPathNames.txt | xargs --max-args=2 --max-lines=1 alignSampleToReference.sh -p $NUMCORES reference/NC_011149.fasta

Step 5 - Prep the samples::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES prepSamples.sh reference/NC_011149.fasta

Step 6 - Combine the SNP positions across all samples into the SNP list file::

    create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt

Step 7 - Create pileups at SNP positions for each sample::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES -I XX create_snp_pileup.py -l snplist.txt -a XX/reads.all.pileup -o XX/reads.snp.pileup

Step 8 - Create the SNP matrix::

    create_snp_matrix.py -l snplist.txt -p reads.snp.pileup -o snpma.fasta sampleDirectories.txt

Step 9 - Create the reference base sequence::

    create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta reference/NC_011149.fasta

        
Step 10 - View and verify the results:

Upon successful completion of the pipeline, the snplist.txt file should have 3624 entries.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l snplist.txt
    ls -l snpma.fasta
    ls -l referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py agonaExpectedResults expectedResults
    diff -q -s snplist.txt         expectedResults/snplist.txt
    diff -q -s snpma.fasta         expectedResults/snpma.fasta
    diff -q -s referenceSNP.fasta  expectedResults/referenceSNP.fasta


Step-by-Step Workflow - General Case
------------------------------------

Step 1 - Gather data:

You will need the following data:

* Reference genome
* Fastq input files for multiple samples

Organize the data into separate directories for each sample as well as the reference.  One possible
directory layout is shown below.  Note the mix of paired and unpaired samples::

    ./myProject/reference/my_reference.fasta
    ./myProject/samples/sample1/sampleA.fastq
    ./myProject/samples/sample2/sampleB.fastq
    ./myProject/samples/sample3/sampleC_1.fastq
    ./myProject/samples/sample3/sampleC_2.fastq
    ./myProject/samples/sample4/sampleD_1.fastq
    ./myProject/samples/sample4/sampleD_2.fastq

Step 2 - Prep work::

    # Optional step: Copy your input data to a safe place:
    cp -r myProject myProjectClean
    # The SNP pipeline will generate additional files into the reference and sample directories
    cd myProject
    
    # Create file of sample directories:
    ls -d --color=never samples/* > sampleDirectories.txt
    
    # get the *.fastq or *.fq files in each sample directory, possibly compresessed, on one line per sample, ready to feed to bowtie
    TMPFILE1=$(mktemp tmp.fastqs.XXXXXXXX)
    cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq* >> $TMPFILE1; echo $dir/*.fq* >> $TMPFILE1; done
    grep -v '*.fq*' $TMPFILE1 | grep -v '*.fastq*' > sampleFullPathNames.txt
    rm $TMPFILE1
    
    # Determine the number of CPU cores in your computer
    NUMCORES=$(grep -c ^processor /proc/cpuinfo)

Step 3 - Prep the reference::

    prepReference.sh reference/my_reference.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    cat sampleFullPathNames.txt | xargs --max-args=2 --max-lines=1 alignSampleToReference.sh -p $NUMCORES reference/my_reference.fasta

Step 5 - Prep the samples::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES prepSamples.sh reference/my_reference.fasta

Step 6 - Combine the SNP positions across all samples into the SNP list file::

    create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt

Step 7 - Create pileups at SNP positions for each sample::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $NUMCORES -I XX create_snp_pileup.py -l snplist.txt -a XX/reads.all.pileup -o XX/reads.snp.pileup

Step 8 - Create the SNP matrix::

    create_snp_matrix.py -l snplist.txt -p reads.snp.pileup -o snpma.fasta sampleDirectories.txt

Step 9 - Create the reference base sequence::

    # Note the .fasta file extension
    create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta reference/my_reference.fasta

Step 10 - View the results:

Upon successful completion of the pipeline, the snplist.txt identifies the SNPs in all samples.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    ls -l snplist.txt
    ls -l snpma.fasta
    ls -l referenceSNP.fasta
