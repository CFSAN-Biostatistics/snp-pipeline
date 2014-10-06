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
| run_snp_pipeline.sh       | | This do-it-all script runs all the other scripts listed below,   |
|                           | | comprising all the pipeline steps                                |
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


Inputs
------

Before using the SNP Pipeline, make sure your input data is organized and
named the way the pipeline expects.  Follow these guidelines:

* No spaces in file names and directory names.

* A fasta genome reference file must exist in a separate directory.

* The samples must be organized with a separate directory for each sample.  
  Each sample directory should contain the fastq files for that sample.  
  The name of the directory should match the name of the sample.
  When using paired-end fastq files, the forward and reverse files must be 
  in the same directory.

* The script needs to know how to find all the samples.  You have two choices:

    #. You can organize all the sample directories under a common parent directory.

    #. You can have sample directories anywhere you like, but you will need to 
       create a file listing the path to all the sample directories.

* The sample fastq files must be named with one of the following file
  patterns: (\*.fastq, \*.fq, \*.fastq.gz, \*.fq.gz).  It's okay if different
  samples are named differently, but the two mate files of paired-end samples
  must be named with the same extension.

Outputs
-------

By default, the SNP Pipeline generates the following output files.  If you 
need more control over the output, you can run the pipeline one step at a time.  
See the section *Step-by-Step Workflows* below.

* snplist.txt : contains a combined list of the SNP positions across all 
  samples in a single unified SNP list file identifing the postions and sample 
  names where SNPs were called.

* reads.snp.pileup : for each sample, the pileup file at the positions where 
  SNPs were called in any of the samples.

* snpma.fasta : the SNP matrix containing the consensus base for each of 
  the samples at the positions where SNPs were called in any of the samples.  
  The matrix contains one row per sample and one column per SNP position.  
  Non-SNP positions are not included in the matrix.  The matrix is formatted 
  as a fasta file, with each sequence (all of identical length) corresponding 
  to the SNPs in the correspondingly named sequence.

* referenceSNP.fasta : a fasta file containing the reference sequence bases at
  all the SNP locations.


All-In-One SNP Pipeline Script
------------------------------

Most users should be able to run the SNP pipeline by launching a single script, 
``run_snp_pipeline.sh``.  This script is easy to use and works equally well on
your desktop workstation or on a High Performance Computing cluster.  You can 
find examples of using the script in the sections below.

If you need more flexibility, you can run the individual pipeline scripts one 
step at a time.  See the section *Step-by-Step Workflows* below.


All-In-One Workflow - Lambda Virus
----------------------------------

The SNP Pipeline software distribution includes a small Lambda Virus data set 
that can be quickly processed to verify the basic functionality of the software.

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

Step 2 - Run the SNP Pipeline::

    # Run the pipeline, specifing the locations of samples and the reference
    #
    # Specify the following options:
    #   -s : samples parent directory
    run_snp_pipeline.sh -s samples reference/lambda_virus.fasta


Step 3 - View and verify the results:

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


All-In-One Workflow - Salmonella Agona
--------------------------------------

The Salmonella Agona data set contains a small number of realistic sequences that 
can be processed in a reasonable amount of time.  Due to the large size of real
data, the sequences must be downloaded from the NCBI SRA.  Follow the instructions 
below to download and process the data set.

This workflow illustrates how to run the SNP Pipeline on a High Performance Computing 
cluster (HPC) running the Torque job queue.  If you do not have a cluster available,
you can still work through this example -- just remove the ``-Q torque`` command line 
option in step 2.

Step 1 - Gather data::

    # The SNP Pipeline distribution includes sample data organized as shown below:
    snppipeline/data/agonaInputs/sha256sumCheck
    snppipeline/data/agonaInputs/reference/NC_011149.fasta

    # Copy the supplied test data to a work area:
    mkdir testAgona
    cd testAgona
    copy_snppipeline_data.py agonaInputs cleanInputs
    cd cleanInputs
    
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
    cd ..

Step 2 - Run the SNP Pipeline::

    # Run the pipeline
    # Specify the following options:
    #   -m : mirror link the input samples and reference files
    #   -o : working directory
    #   -s : samples parent directory
    #   -Q : HPC job queue manager
    run_snp_pipeline.sh -m -Q torque -o work -s cleanInputs/samples cleanInputs/reference/NC_011149.fasta
      
Step 3 - View and verify the results:

Upon successful completion of the pipeline, the snplist.txt file should have 3624 entries.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l work/snplist.txt
    ls -l work/snpma.fasta
    ls -l work/referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py agonaExpectedResults expectedResults
    diff -q -s work/snplist.txt         expectedResults/snplist.txt
    diff -q -s work/snpma.fasta         expectedResults/snpma.fasta
    diff -q -s work/referenceSNP.fasta  expectedResults/referenceSNP.fasta

All-In-One Workflow - Listeria monocytogenes
--------------------------------------------

This Listeria monocytogene data set is based on an oubreak investigation related
to contamination in stone fruit. It only contains environmental/produce isolates,
though the full investigation contained data obtained from clinical samples as well.
Due to the large size of the data, the sequences must be downloaded from the NCBI
SRA.  The instructions below show how to create the data set and process it. 
We do the processing with the run_snp_pipeline.sh script, which does much of the
work in one step, but provides less insight into (and control of) the analysis
process.  

Step 1 - Create dataset::


    # The SNP Pipeline distribution does not include the sample data, but does
    #   include information about the sample data, as well as the reference
    #   sequence:
    snppipeline/data/listeriaInputs/sha256sumCheck
    snppipeline/data/listeriaInputs/reference/CFSAN023463.HGAP.draft.fasta
    snppipeline/data/listeriaInputs/sampleList

    # Copy the supplied test data to a work area:
    mkdir testDir
    cd testDir
    copy_snppipeline_data.py listeriaInputs cleanInputs
    cd cleanInputs
    
    # Create sample directories and download sample data from SRA at NCBI. Note that
    #   we use the fastq-dump command from the NCBI SRA-toolkit to fetch sample
    #   sequences. There are other ways to get the data, but the SRA-toolkit is
    #   easy to install, and does a good job of downloading large files.
    mkdir samples
    < sampleList xargs -I % sh -c ' mkdir samples/%; fastq-dump --split-files --outdir samples/% %;'

    # Check the data
    #   The original data was used to generate a hash as follows:
    #     sha256sum sampleList reference/*.fasta samples/*/*.fastq > sha256sumCheck
    #   The command below checks the downloaded data (and the reference sequence) against the
    #     hashes that are saved in the sha256sumCheck file using sha256sum command, which is
    #     generally available on unix systems.
    sha256sum -c sha256sumCheck

    
Step 2 - Run the SNP Pipeline::

    # Run the pipeline. In this case we show the command line options as we might use them
    #   to run on a large workstation. Depending on the amount of memory and number of cores
    #   on your workstation, there are a couple of parameters you may want/need to adjust
    #   for this analysis or other analysis work that your do. These parameters are the
    #   number of cores that are used, and the amount of memory that is used by the java
    #   virtual machine. The number of cores can be altered by changing the 'numCores' 
    #   variable in the run_snp_pipeline.sh script. The amount of memory used by the
    #   javavm can be set by using the -Xmx flag in the call to java in the prepSamples.sh
    #   script. Remember that if you have installed this code using a python virtual
    #   environment, you will need to re-run 'python setup.py develop' again, or you will
    #   be wondering why your changes are not affecting anything. 
    run_snp_pipeline.sh -o outputDirectory -s data/samples data/reference/CFSAN023463.HGAP.draft.fasta

Step 3 - View and verify the results::

Upon successful completion of the pipeline, the snplist.txt file should have 11,746
entries.  The SNP Matrix can be found in snpma.fasta.  The corresponding reference
bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l work/snplist.txt
    ls -l work/snpma.fasta
    ls -l work/referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py listeriaExpectedResults expectedResults
    diff -q -s work/snplist.txt         expectedResults/snplist.txt
    diff -q -s work/snpma.fasta         expectedResults/snpma.fasta
    diff -q -s work/referenceSNP.fasta  expectedResults/referenceSNP.fasta


Step-by-Step Workflow - Lambda Virus 
------------------------------------

The SNP Pipeline software distribution includes a small Lambda Virus data set 
that can be quickly processed to verify the basic functionality of the software.

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
    ls -d samples/* > sampleDirectories.txt
    rm sampleFullPathNames.txt 2>/dev/null
    cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq >> sampleFullPathNames.txt; done
    # Determine the number of CPU cores in your computer
    NUMCORES=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)

Step 3 - Prep the reference::

    prepReference.sh reference/lambda_virus.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    cat sampleFullPathNames.txt | xargs -n 2 -L 1 alignSampleToReference.sh -p $NUMCORES reference/lambda_virus.fasta

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



Step-by-Step Workflow - Salmonella Agona
----------------------------------------

The Salmonella Agona data set contains realistic sequences that can be processed
in a reasonable amount of time.  Due to the large size of real data, the sequences
must be downloaded from the NCBI SRA.  Follow the instructions below to download 
and process the data set.

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
    ls -d samples/* > sampleDirectories.txt
    rm sampleFullPathNames.txt 2>/dev/null
    cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq >> sampleFullPathNames.txt; done
    # Determine the number of CPU cores in your computer
    NUMCORES=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)

Step 3 - Prep the reference::

    prepReference.sh reference/NC_011149.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    cat sampleFullPathNames.txt | xargs -n 2 -L 1 alignSampleToReference.sh -p $NUMCORES reference/NC_011149.fasta

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
    ls -d samples/* > sampleDirectories.txt
    
    # get the *.fastq or *.fq files in each sample directory, possibly compresessed, on one line per sample, ready to feed to bowtie
    TMPFILE1=$(mktemp tmp.fastqs.XXXXXXXX)
    cat sampleDirectories.txt | while read dir; do echo $dir/*.fastq* >> $TMPFILE1; echo $dir/*.fq* >> $TMPFILE1; done
    grep -v '*.fq*' $TMPFILE1 | grep -v '*.fastq*' > sampleFullPathNames.txt
    rm $TMPFILE1
    
    # Determine the number of CPU cores in your computer
    NUMCORES=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)

Step 3 - Prep the reference::

    prepReference.sh reference/my_reference.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    cat sampleFullPathNames.txt | xargs -n 2 -L 1 alignSampleToReference.sh -p $NUMCORES reference/my_reference.fasta

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
