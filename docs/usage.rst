.. _usage-label:

========
Usage
========

The SNP Pipeline is run from the Unix command line.  The pipeline consists of a collection
of shell scripts and python scripts:

    * prepReference.sh : indexes the reference genome
    * prepSamples.sh : finds variants in each sample
    * runsnppipeline.py : creates a matrix of SNPs across all samples

Step-by-Step Example Workflow Based on Lamda Virus Test Data Provided with Code
-------------------------------------------------------------------------------

Step 1 - Gather data::

    # The SNP Pipeline distribution includes sample data organized as shown below:
    test/testLambdaVirusClean/reference/lambda_virus.fasta
    test/testLambdaVirusClean/samples/sample1/sample1_1.fastq
    test/testLambdaVirusClean/samples/sample1/sample1_2.fastq
    test/testLambdaVirusClean/samples/sample2/sample2_1.fastq
    test/testLambdaVirusClean/samples/sample2/sample2_2.fastq
    test/testLambdaVirusClean/samples/sample3/sample3_1.fastq
    test/testLambdaVirusClean/samples/sample3/sample3_2.fastq
    test/testLambdaVirusClean/samples/sample4/sample4_1.fastq
    test/testLambdaVirusClean/samples/sample4/sample4_2.fastq

Step 2 - Prep work::

    # Copy the supplied test data to a work area:
    cd test
    cp -r testLambdaVirusClean testLambdaVirus
    cd testLambdaVirus
    # Create files of sample directories and fastQ files:
    ls -d --color=never samples/* > sampleDirectoryNames.txt
    find samples -type f | grep fastq | sort -u > sampleFullPathNames.txt

Step 3 - Prep the reference::

    prepReference.sh reference/lambda_virus

Step 4 - Prep the samples::

    # Note: This could be run in parallel using gnu parallel (workstation) or qsub (PBS on HPC)
    # Note: We use xargs parameter "-n 2" because the samples are paired
    cat sampleFullPathNames.txt | xargs -n 2 prepSamples.sh reference/lambda_virus
        
Step 5 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to
generate snp matrix)::

    runsnppipeline.py -n 10 -d ./ -f sampleDirectoryNames.txt -r reference/lambda_virus.fasta -l snplist.txt -a snpma.fasta -i True


Step-by-Step Example Workflow Based on S. Agona Data Downloaded from SRA
------------------------------------------------------------------------
TODO: do this


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
    # Create files of sample directories:
    ls -d --color=never samples/* > sampleDirectoryNames.txt

Step 3 - Prep the reference::

    # Note: do not specify the .fasta file extension here
    prepReference.sh reference/reference/my_reference

Step 4 - Prep the samples::

    # Run prepSamples once per sample, note the mix of paired and unpaired sample here
    # Note: This could be run in parallel using gnu parallel (workstation) or qsub (PBS on HPC)
    prepSamples.sh  reference/my_reference  samples/sample1/sampleA.fastq
    prepSamples.sh  reference/my_reference  samples/sample2/sampleB.fastq
    prepSamples.sh  reference/my_reference  samples/sample3/sampleC_1.fastq  samples/sample3/sampleC_2.fastq
    prepSamples.sh  reference/my_reference  samples/sample4/sampleD_1.fastq  samples/sample4/sampleD_2.fastq

Step 5 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to
generate snp matrix)::

    runsnppipeline.py -n 10 -d ./ -f sampleDirectoryNames.txt -r reference/my_reference.fasta -l snplist.txt -a snpma.fasta -i True

runsnppipeline.py Command Syntax
--------------------------------
Help for the SNP Pipeline command-line arguments can be found with the --help parameter::


    runsnppipeline.py  --help


