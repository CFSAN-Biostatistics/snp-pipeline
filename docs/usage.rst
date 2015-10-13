.. _usage-label:

========
Usage
========

.. highlight:: bash



The SNP Pipeline is run from the Unix command line.  The pipeline consists of a collection
of shell scripts and python scripts.


+-----------------------------+--------------------------------------------------------------------+
| Script                      | | Description                                                      |
+=============================+====================================================================+
| copy_snppipeline_data.py    | | Copies supplied example data to a work directory                 |
+-----------------------------+--------------------------------------------------------------------+
| run_snp_pipeline.sh         | | This do-it-all script runs all the other scripts listed below,   |
|                             | | comprising all the pipeline steps                                |
+-----------------------------+--------------------------------------------------------------------+
| prepReference.sh            | | Indexes the reference genome                                     |
+-----------------------------+--------------------------------------------------------------------+
| alignSampleToReference.sh   | | Aligns samples to the reference genome                           |
+-----------------------------+--------------------------------------------------------------------+
| prepSamples.sh              | | Finds variants in each sample                                    |
+-----------------------------+--------------------------------------------------------------------+
| create_snp_list.py          | | Combines the SNP positions across all samples into a single      |
|                             | | unified SNP list file                                            |
+-----------------------------+--------------------------------------------------------------------+
| create_snp_pileup.py        | | Deprecated -- this command is not used by the pipeline since     |
|                             | | v0.3.5.  Replaced by call_consensus.py                           |
|                             | |                                                                  |
|                             | | Creates the SNP pileup file for a sample -- a subset of the      |
|                             | | pileup file at only the positions where SNPs were called in any  |
|                             | | of the samples                                                   |
+-----------------------------+--------------------------------------------------------------------+
| call_consensus.py           | | Calls the consensus SNPs for each sample                         |
+-----------------------------+--------------------------------------------------------------------+
| create_snp_matrix.py        | | Creates a matrix of SNPs across all samples                      |
+-----------------------------+--------------------------------------------------------------------+
| create_snp_reference_seq.py | | Writes the reference sequence bases at SNP locations to          |
|                             | | a fasta file                                                     |
+-----------------------------+--------------------------------------------------------------------+
| collectSampleMetrics.sh     | | Collects useful coverage and variant statistics about            |
|                             | | each sample                                                      |
+-----------------------------+--------------------------------------------------------------------+
| combineSampleMetrics.sh     | | Creates a table of coverage and variant statistics for           |
|                             | | all samples                                                      |
+-----------------------------+--------------------------------------------------------------------+
| mergeVcf.sh                 | | Creates a multi-sample VCF file with the snps found in all       |
|                             | | samples                                                          |
+-----------------------------+--------------------------------------------------------------------+


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
See :ref:`step-by-step-workflows`.

* snplist.txt : contains a combined list of the SNP positions across all 
  samples in a single unified SNP list file identifing the positions and sample 
  names where SNPs were called.

* consensus.fasta : for each sample, the consensus base call at the positions 
  where SNPs were previously detected in any of the samples.

* consensus.vcf : for each sample, the VCF file of snps called, as well as 
  failed snps at the positions where SNPs were previously detected in any of 
  the samples.

* snpma.fasta : the SNP matrix containing the consensus base for each of 
  the samples at the positions where SNPs were called in any of the samples.  
  The matrix contains one row per sample and one column per SNP position.  
  Non-SNP positions are not included in the matrix.  The matrix is formatted 
  as a fasta file, with each sequence (all of identical length) corresponding 
  to the SNPs in the correspondingly named sequence.

* snpma.vcf : contains the merged multi-sample VCF file identifying the positions
  and snps for all samples.

* referenceSNP.fasta : a fasta file containing the reference sequence bases at
  all the SNP locations.

* metrics : for each sample, contains the size of the sample, number of reads, 
  alignment rate, pileup depth, and number of SNPs found.

* metrics.tsv : a tab-separated table of metrics for all samples containing 
  the size of the samples, number of reads, alignment rate, pileup depth, and 
  number of SNPs found.

.. _all-in-one-script-label:

All-In-One SNP Pipeline Script
------------------------------

Most users should be able to run the SNP Pipeline by launching a single script, 
``run_snp_pipeline.sh``.  This script is easy to use and works equally well on
your desktop workstation or on a High Performance Computing cluster.  You can 
find examples of using the script in the sections below.

If you need more flexibility, you can run the individual pipeline scripts one 
step at a time.  See :ref:`step-by-step-workflows`.

.. _logging-label:

Logging
-------

When the SNP Pipeline is launched with the ``run_snp_pipeline.sh`` script,
it generates log files for each processing step of the pipeline.  The logs for 
each pipeline run are stored in a time-stamped directory under the output directory.
If the pipeline is re-run on the same samples, the old log files are kept and
a new log directory is created for the new run.  For example, the output 
directory might look like this after two runs::

    drwx------ 2 me group 4096 Oct 17 16:37 logs-20141017.154428/
    drwx------ 2 me group 4096 Oct 17 16:38 logs-20141017.163848/
    drwx------ 2 me group 4096 Oct 17 16:37 reference/
    -rw------- 1 me group  194 Oct 17 16:38 referenceSNP.fasta
    -rw------- 1 me group  104 Oct 17 16:38 sampleDirectories.txt
    drwx------ 6 me group 4096 Oct 17 16:37 samples/
    -rw------- 1 me group 7216 Oct 17 16:38 snplist.txt
    -rw------- 1 me group  708 Oct 17 16:38 snpma.fasta

A log file is created for each step of the pipeline for each sample.  For 
performamnce reasons, the samples are sorted by size and processed largest
first.  This sorting is reflected in the naming of the log files.  The log files
are named with a suffix indicating the sample number::

    -rw------- 1 me group  1330 Oct 17 16:37 alignSamples.log-1
    -rw------- 1 me group  1330 Oct 17 16:37 alignSamples.log-2
    -rw------- 1 me group  1330 Oct 17 16:37 alignSamples.log-3
    -rw------- 1 me group 12045 Oct 17 16:37 prepReference.log
    -rw------- 1 me group  1686 Oct 17 16:37 prepSamples.log-1
    -rw------- 1 me group  1686 Oct 17 16:37 prepSamples.log-2
    -rw------- 1 me group  1686 Oct 17 16:37 prepSamples.log-3
    -rw------- 1 me group   983 Oct 17 16:37 snpList.log
    -rw------- 1 me group  1039 Oct 17 16:37 snpMatrix.log
    -rw------- 1 me group   841 Oct 17 16:37 snpPileup.log-1
    -rw------- 1 me group   841 Oct 17 16:37 snpPileup.log-2
    -rw------- 1 me group   841 Oct 17 16:37 snpPileup.log-3
    -rw------- 1 me group   806 Oct 17 16:37 snpReference.log

To determine which samples correspond to which log files, you can either grep the
log files for the sample name or inspect the sorted sampleDirectories.txt file to determine
the sequential position of the sample.  The file names are consistent regardless of whether 
the pipeline is run on a workstation or HPC cluster.

In addition to the processing log files, the log directory also contains a copy of the
configuration file used for each run -- capturing the parameters used during the run.


.. _mirrored-input-label:

Mirrored Inputs
---------------

When the SNP Pipeline is launched with the ``run_snp_pipeline.sh`` script, it has the
optional capability to create a mirrored copy of the input fasta and fastq files.  You 
might use this feature to avoid polluting the reference directory and sample directories 
with the intermediate files generated by the snp pipeline.  The mirroring function can 
either create normal copies of the files, or it can create links to the original files 
-- saving both time and disk space.  With linked files, you can easily run multiple 
experiments on the same data or different overlapping sets of samples without having 
duplicate copies of the original sample files.  See the :ref:`cmd-ref-run-snp-pipeline` 
command reference for the mirroring syntax.

The mirroring function creates a "reference" subdirectory and a "samples" subdirectory under
the main output directory.  One directory per sample is created under the "samples" directory.  
The generated intermediate files are placed into the mirrored directories, not in the original
locations of the inputs. The SNP Pipeline attempts to preserve the time stamps of the original 
files in the mirrored directories.

Keep in mind the following limitations when mirroring the inputs.

* Some file systems do not support soft (symbolic) links.  If you attempt to create a soft link
  on a file system without the capability, the operation will fail with an error message.
* Hard links cannot be used to link files across two different file systems.  The original 
  file and the link must both reside on the same file system.
* Normal file copies should always work, but the copy operation can be lengthy and the duplicate 
  files will consume extra storage space.


.. _hpc-usage-label:

High Performance Computing
--------------------------
The SNP Pipeline can be executed on a High Performamce Computing cluster.  The
Torque and Grid Engine job queue managers are supported.

Torque
~~~~~~
To run the SNP Pipeline on torque::

    run_snp_pipeline.sh -Q torque -s mySamplesDir myReference.fasta

Grid Engine
~~~~~~~~~~~
To run the SNP Pipeline on grid engine you must use a configuration file to specify
the name of your parallel environment.

Grab the default configuration file::

    copy_snppipeline_data.py configurationFile


Edit the snppipeline.conf file and make the following change::
    
    PEname="myPE" # substitute the name of your PE

Then run the pipeline with the -c and -Q command line options::    
    
    run_snp_pipeline.sh -c snppipeline.conf -Q grid -s mySamplesDir myReference.fasta

See also: :ref:`faq-performance-label`.


.. _tool-selection-label:

Tool Selection
--------------
The SNP Pipeline lets you choose either the Bowtie2 aligner or the Smalt aligner.  Your choice
of aligner, as well as the command line options for the aligner are specified in the
SNP Pipeline configuration file.

Grab the default configuration file::

    copy_snppipeline_data.py configurationFile

To run the SNP Pipeline with Bowtie2, edit ``snppipeline.conf`` with these settings::

    SnpPipeline_Aligner="bowtie2"
    Bowtie2Build_ExtraParams="" # substitute the command line options you want here
    Bowtie2Align_ExtraParams="" # substitute the command line options you want here

To run the SNP Pipeline with Smalt, edit ``snppipeline.conf`` with these settings::

    SnpPipeline_Aligner="smalt"
    SmaltIndex_ExtraParams="" # substitute the command line options you want here
    SmaltAlign_ExtraParams="" # substitute the command line options you want here

Then run the pipeline with the -c command line option::    
    
    run_snp_pipeline.sh -c snppipeline.conf -s mySamplesDir myReference.fasta
    
See also :ref:`configuration-label`.


All-In-One SNP Pipeline Workflows
---------------------------------
The sections below give detailed examples of workflows you can run with the
all-in-one run_snp_pipeline.sh script.

| :ref:`all-in-one-workflow-lambda`
| :ref:`all-in-one-workflow-agona`
| :ref:`all-in-one-workflow-listeria`
|


.. _all-in-one-workflow-lambda:

All-In-One Workflow - Lambda Virus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    ls -l snpma.vcf
    ls -l referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py lambdaVirusExpectedResults expectedResults
    diff -q -s snplist.txt         expectedResults/snplist.txt
    diff -q -s snpma.fasta         expectedResults/snpma.fasta
    diff -q -s referenceSNP.fasta  expectedResults/referenceSNP.fasta

    # View the per-sample metrics
    xdg-open metrics.tsv

.. _all-in-one-workflow-agona:

All-In-One Workflow - Salmonella Agona
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Salmonella Agona data set contains a small number of realistic sequences that 
can be processed in a reasonable amount of time.  Due to the large size of real
data, the sequences must be downloaded from the NCBI SRA.  Follow the instructions 
below to download and process the data set.

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
    #   -m : mirror the input samples and reference files
    #   -o : output directory
    #   -s : samples parent directory
    run_snp_pipeline.sh -m soft -o outputDirectory -s cleanInputs/samples cleanInputs/reference/NC_011149.fasta
      
Step 3 - View and verify the results:

Upon successful completion of the pipeline, the snplist.txt file should have 3624 entries.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l outputDirectory/snplist.txt
    ls -l outputDirectory/snpma.fasta
    ls -l outputDirectory/snpma.vcf
    ls -l outputDirectory/referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py agonaExpectedResults expectedResults
    diff -q -s outputDirectory/snplist.txt         expectedResults/snplist.txt
    diff -q -s outputDirectory/snpma.fasta         expectedResults/snpma.fasta
    diff -q -s outputDirectory/referenceSNP.fasta  expectedResults/referenceSNP.fasta

    # View the per-sample metrics
    xdg-open outputDirectory/metrics.tsv

.. _all-in-one-workflow-listeria:

All-In-One Workflow - Listeria monocytogenes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This Listeria monocytogene data set is based on an oubreak investigation related
to contamination in stone fruit. It only contains environmental/produce isolates,
though the full investigation contained data obtained from clinical samples as well.
Due to the large size of the data, the sequences must be downloaded from the NCBI
SRA.  The instructions below show how to create the data set and process it. 
We do the processing with the run_snp_pipeline.sh script, which does much of the
work in one step, but provides less insight into (and control of) the analysis
process.  

This workflow illustrates how to run the SNP Pipeline on a High Performance Computing 
cluster (HPC) running the Torque job queue manager.  If you do not have a cluster available,
you can still work through this example -- just remove the ``-Q torque`` command line 
option in step 2.

Step 1 - Create dataset::


    # The SNP Pipeline distribution does not include the sample data, but does
    #   include information about the sample data, as well as the reference
    #   sequence.  The files are organized as shown below:
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
    cd ..
    
Step 2 - Run the SNP Pipeline:

There are a couple of parameters you may need to adjust for this analysis or other analysis
work that your do. These parameters are the number of CPU cores that are used, and the 
amount of memory that is used by the java virtual machine.  Both can be set in a
configuration file you can pass to run_snp_pipeline.sh with the ``-c`` option.  
See :ref:`faq-performance-label`.

Launch the pipeline::

    # Run the pipeline. 
    # Specify the following options:
    #   -m : mirror the input samples and reference files
    #   -Q : HPC job queue manager
    #   -o : output directory
    #   -s : samples parent directory
    run_snp_pipeline.sh -m soft -Q torque -o outputDirectory -s cleanInputs/samples cleanInputs/reference/CFSAN023463.HGAP.draft.fasta

Step 3 - View and verify the results:

Upon successful completion of the pipeline, the snplist.txt file should have 11,787
entries.  The SNP Matrix can be found in snpma.fasta.  The corresponding reference
bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l outputDirectory/snplist.txt
    ls -l outputDirectory/snpma.fasta
    ls -l outputDirectory/snpma.vcf
    ls -l outputDirectory/referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py listeriaExpectedResults expectedResults
    diff -q -s outputDirectory/snplist.txt         expectedResults/snplist.txt
    diff -q -s outputDirectory/snpma.fasta         expectedResults/snpma.fasta
    diff -q -s outputDirectory/referenceSNP.fasta  expectedResults/referenceSNP.fasta

    # View the per-sample metrics
    xdg-open outputDirectory/metrics.tsv

.. _step-by-step-workflows:

Step-by-Step Workflows
----------------------

The run_snp_pipeline.sh script described above provides a simplified interface
for running all the pipeline steps from a single command.  If you need more
control over the inputs, outputs, or processing steps, you can run the pipeline 
one step at a time.

The sections below give detailed examples of workflows you can run with the
component tools of the pipeline.

| :ref:`step-by-step-workflow-lambda`
| :ref:`step-by-step-workflow-agona`
| :ref:`step-by-step-workflow-general-case`
|


.. _step-by-step-workflow-lambda:

Step-by-Step Workflow - Lambda Virus 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    numCores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)

Step 3 - Prep the reference::

    prepReference.sh reference/lambda_virus.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cat sampleFullPathNames.txt | xargs -n 2 -L 1 alignSampleToReference.sh reference/lambda_virus.fasta

Step 5 - Prep the samples::

    # Process the samples in parallel using all CPU cores
    export VarscanMpileup2snp_ExtraParams="--min-var-freq 0.90"
    cat sampleDirectories.txt | xargs -n 1 -P $numCores prepSamples.sh reference/lambda_virus.fasta

Step 6 - Combine the SNP positions across all samples into the SNP list file::

    create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt

Step 7 - Call the consensus base at SNP positions for each sample::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $numCores -I XX call_consensus.py -l snplist.txt --vcfFileName consensus.vcf -o XX/consensus.fasta XX/reads.all.pileup

Step 8 - Create the SNP matrix::

    create_snp_matrix.py -c consensus.fasta -o snpma.fasta sampleDirectories.txt

Step 9 - Create the reference base sequence::

    create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta reference/lambda_virus.fasta

Step 10 - Collect metrics for each sample::

    cat sampleDirectories.txt | xargs -n 1 -P $numCores -I XX collectSampleMetrics.sh -m snpma.fasta -o XX/metrics XX reference/lambda_virus.fasta

Step 11 - Tabulate the metrics for all samples::

    combineSampleMetrics.sh -n metrics -o metrics.tsv sampleDirectories.txt

Step 12 - Merge the VCF files for all samples into a multi-sample VCF file::

    mergeVcf.sh -n consensus.vcf -o snpma.vcf sampleDirectories.txt

Step 13 - View and verify the results:

Upon successful completion of the pipeline, the snplist.txt file should have 165 entries.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l snplist.txt
    ls -l snpma.fasta
    ls -l snpma.vcf
    ls -l referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py lambdaVirusExpectedResults expectedResults
    diff -q -s snplist.txt         expectedResults/snplist.txt
    diff -q -s snpma.fasta         expectedResults/snpma.fasta
    diff -q -s referenceSNP.fasta  expectedResults/referenceSNP.fasta

    # View the per-sample metrics
    xdg-open metrics.tsv


.. _step-by-step-workflow-agona:

Step-by-Step Workflow - Salmonella Agona
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    numCores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)

Step 3 - Prep the reference::

    prepReference.sh reference/NC_011149.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cat sampleFullPathNames.txt | xargs -n 2 -L 1 alignSampleToReference.sh reference/NC_011149.fasta

Step 5 - Prep the samples::

    # Process the samples in parallel using all CPU cores
    export VarscanMpileup2snp_ExtraParams="--min-var-freq 0.90"
    cat sampleDirectories.txt | xargs -n 1 -P $numCores prepSamples.sh reference/NC_011149.fasta

Step 6 - Combine the SNP positions across all samples into the SNP list file::

    create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt

Step 7 - Call the consensus base at SNP positions for each sample::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $numCores -I XX call_consensus.py -l snplist.txt --vcfFileName consensus.vcf -o XX/consensus.fasta XX/reads.all.pileup

Step 8 - Create the SNP matrix::

    create_snp_matrix.py -c consensus.fasta -o snpma.fasta sampleDirectories.txt

Step 9 - Create the reference base sequence::

    create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta reference/NC_011149.fasta

Step 10 - Collect metrics for each sample::

    cat sampleDirectories.txt | xargs -n 1 -P $numCores -I XX collectSampleMetrics.sh -m snpma.fasta -o XX/metrics XX reference/NC_011149.fasta

Step 11 - Tabulate the metrics for all samples::

    combineSampleMetrics.sh -n metrics -o metrics.tsv sampleDirectories.txt

Step 12 - Merge the VCF files for all samples into a multi-sample VCF file::

    mergeVcf.sh -n consensus.vcf -o snpma.vcf sampleDirectories.txt

Step 13 - View and verify the results:

Upon successful completion of the pipeline, the snplist.txt file should have 3624 entries.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    # Verify the result files were created
    ls -l snplist.txt
    ls -l snpma.fasta
    ls -l snpma.vcf
    ls -l referenceSNP.fasta

    # Verify correct results
    copy_snppipeline_data.py agonaExpectedResults expectedResults
    diff -q -s snplist.txt         expectedResults/snplist.txt
    diff -q -s snpma.fasta         expectedResults/snpma.fasta
    diff -q -s referenceSNP.fasta  expectedResults/referenceSNP.fasta

    # View the per-sample metrics
    xdg-open metrics.tsv

.. _step-by-step-workflow-general-case:

Step-by-Step Workflow - General Case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    numCores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)

Step 3 - Prep the reference::

    prepReference.sh reference/my_reference.fasta

Step 4 - Align the samples to the reference::

    # Align each sample, one at a time, using all CPU cores
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cat sampleFullPathNames.txt | xargs -n 2 -L 1 alignSampleToReference.sh reference/my_reference.fasta

Step 5 - Prep the samples::

    # Process the samples in parallel using all CPU cores
    export VarscanMpileup2snp_ExtraParams="--min-var-freq 0.90"
    cat sampleDirectories.txt | xargs -n 1 -P $numCores prepSamples.sh reference/my_reference.fasta

Step 6 - Combine the SNP positions across all samples into the SNP list file::

    create_snp_list.py -n var.flt.vcf -o snplist.txt sampleDirectories.txt

Step 7 - Call the consensus base at SNP positions for each sample::

    # Process the samples in parallel using all CPU cores
    cat sampleDirectories.txt | xargs -n 1 -P $numCores -I XX call_consensus.py -l snplist.txt --vcfFileName consensus.vcf -o XX/consensus.fasta XX/reads.all.pileup

Step 8 - Create the SNP matrix::

    create_snp_matrix.py -c consensus.fasta -o snpma.fasta sampleDirectories.txt

Step 9 - Create the reference base sequence::

    # Note the .fasta file extension
    create_snp_reference_seq.py -l snplist.txt -o referenceSNP.fasta reference/my_reference.fasta

Step 10 - Collect metrics for each sample::

    cat sampleDirectories.txt | xargs -n 1 -P $numCores -I XX collectSampleMetrics.sh -m snpma.fasta -o XX/metrics XX reference/my_reference.fasta

Step 11 - Tabulate the metrics for all samples::

    combineSampleMetrics.sh -n metrics -o metrics.tsv sampleDirectories.txt

Step 12 - Merge the VCF files for all samples into a multi-sample VCF file::

    mergeVcf.sh -n consensus.vcf -o snpma.vcf sampleDirectories.txt

Step 13 - View the results:

Upon successful completion of the pipeline, the snplist.txt identifies the SNPs in all samples.  The SNP Matrix 
can be found in snpma.fasta.  The corresponding reference bases are in the referenceSNP.fasta file::

    ls -l snplist.txt
    ls -l snpma.fasta
    ls -l snpma.vcf
    ls -l referenceSNP.fasta

    # View the per-sample metrics
    xdg-open metrics.tsv


.. _metrics-usage-label:

Metrics
-------

After creating the SNP matrix, the pipeline collects and tabulates metrics for all of the samples.  The metrics 
are first collected in one file per sample in the sample directories.  A subsequent step combines the
metrics for all the samples together into a single tab-separated file with one row per sample and one column
per metric.  The tabulated metrics file is named metrics.tsv by default.

The metrics are:

+-------------------------+------------------------------------------------------------------+
| Column                  | | Description                                                    |
+=========================+==================================================================+
| Sample                  | | The name of the directory containing the sample fastq files.   |
+-------------------------+------------------------------------------------------------------+
| Fastq Files             | | Comma separated list of fastq file names in the sample         |
|                         | | directory.                                                     |
+-------------------------+------------------------------------------------------------------+
| Fastq File Size         | | The sum of the sizes of the fastq files. This will be the      |
|                         | | compressed size if the files are compressed.                   | 
+-------------------------+------------------------------------------------------------------+
| Machine                 | | The sequencing instrument ID extracted from the compressed     |
|                         | | fastq.gz file header.  If the fastq files are not compressed,  |
|                         | | the machine ID is not captured.                                |
+-------------------------+------------------------------------------------------------------+
| Flowcell                | | The flowcell used during the sequencing run, extracted from    |
|                         | | the compressed fastq.gz file header. If the fastq files are    |
|                         | | not compressed, the flowcell is not captured.                  |
+-------------------------+------------------------------------------------------------------+
| Number of Reads         | | The number of reads in the SAM file.  When using paired fastq  |
|                         | | files, this number will be twice the number of reads reported  |
|                         | | by bowtie.                                                     |
+-------------------------+------------------------------------------------------------------+
| Percent of Reads Mapped | | The percentage of reference-aligned reads in the SAM file.     |
+-------------------------+------------------------------------------------------------------+
| Average Pileup Depth    | | The average depth of coverage in the sample pileup file.  This |
|                         | | is calculated as the sum of the depth of the pileup across all |
|                         | | pileup positions divided by the number of positions in the     |
|                         | | reference.                                                     |
+-------------------------+------------------------------------------------------------------+
| Number of SNPs          | | The number of SNPs found for this sample.  The count is        |
|                         | | computed as the number of SNP records in the VCF file          | 
|                         | | generated by the snp caller (VarScan).                         |
+-------------------------+------------------------------------------------------------------+
| Missing SNP Matrix      | | The number of positions in the SNP matrix for which a          |
| Positions               | | consensus base could not be called for this sample.  The       |
|                         | | inability to call a consensus base is caused by either a       |
|                         | | pileup file with no coverage at a SNP position, or by          |
|                         | | insufficient agreement among the pileup bases at the SNP       |
|                         | | position.  The minimum fraction of reads that must agree at a  |
|                         | | position to make a consensus call is controlled by the         |
|                         | | ``minConsFreq`` parameter.                                     |
+-------------------------+------------------------------------------------------------------+
| Warnings and Errors     | | A list of warnings or errors encountered while collecting the  |
|                         | | metrics.                                                       |
+-------------------------+------------------------------------------------------------------+




