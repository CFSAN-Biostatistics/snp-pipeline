===========================
FAQ / Troubleshooting Guide
===========================

.. highlight:: bash

Installation
------------

**Q: How can I avoid polluting my global python installation when installing the SNP pipeline?**


A: You can either use a python virtual environment or install into your user area.  To use a python virtual 
environment, see the :ref:`get-started-label` section for developers who want to contribute.  To install into 
your user area instead of installing into your global site packages, do this::

	$ pip install --user snp-pipeline

**Q: The SNP Pipeline cannot find VarScan.  How should I install it?**

A: Download the VarScan jar file from SourceForge.  Put the jar file anywhere.  You need read-access to the
jar file, but not execute-access.  The suppled shell scripts expect the CLASSPATH environment variable to 
specify the path to the VarScan jar file.  The CLASSPATH should include the filename, not just the directory.
Define it like this in your .bashrc file::

    export CLASSPATH=~/software/varscan.v2.3.6/VarScan.v2.3.6.jar:$CLASSPATH


**Q: How can I uninstall the SNP pipeline?**

A: If you installed with pip, you can uninstall from the command line::

    $ pip uninstall snp-pipeline


**Q: How can I rollback to an older version of the SNP pipeline?**

A: You can revert to an older version with these commands::

    $ pip uninstall snp-pipeline
    $ pip install --user snp-pipeline==0.3.2  # substitute the version you want here


**Q: Is there a way to install a specific release of the SNP pipeline from the github repository?**

A: Yes, you can install a release from github with this command::

    $ pip install --user https://github.com/CFSAN-Biostatistics/snp-pipeline/archive/v0.3.2  # substitute the version you want here


Running the Pipeline
--------------------

**Q: Nothing works.**

A: Make sure you have the proper dependencies on your path.  Modify your path if necessary to include bowtie 
and samtools.  See the question above about installing VarScan.  Also make sure you install Biopython.  If 
you are using virtual environments, make sure you issue the command "toggleglobalsitepackages" so Biopython 
can be found.  See the :ref:`installation-label` section of this documentation.

**Q: How can I verify the pipeline is installed and working properly?**

A: The SNP Pipeline includes two small sets of test data with result files.  You can run the pipeline against the 
test data to verify correct results.  Follow the lambda virus workflow steps here: :ref:`all-in-one-workflow-lambda`.

Upon successful completion of the pipeline, the snplist.txt file should have 165 entries.  The SNP Matrix 
can be found in snpma.fasta.  It should have the following contents::

    >sample1
    AGCACCGGGACCCACGGCGCACGCAAAGATCCGAATTGCAGGGCGTACCTGGACCCCGGT
    GACGGGGGATCGGGGACTCTTGGTGAGGAACTAAAACGAACATCCACGTTTTCATGGCGA
    CTGCTTGCCAGGTGTCAGCACATTCCCTATATCGGTGGACACGTA
    >sample2
    GGCGCTAGGAGGCAAGCCTTGGTCGTGGTTAATAGTTACAAGGCGTGCGCGTACTGCCGT
    CTCCTACTATCTCTGCCGCCTCTCGCGATCCGGACCGCAACACCAACTCTCTGGTGGCAT
    CCTCTGAATCGTCGTGAGCATCTCAATTATATATTCGTCCGCGCG
    >sample3
    GGCGCTAGGAGGTACGCTTCGCGTGTGGATCAGCGCTACGGTGCCTATGCGTGACCCGCG
    GAACTGGGTTCGCGTAAGGCAGTTCAGGTACGGCAACGTAGATCAAAGTTTAGAAACCAT
    ACTCGTAATCCGCCTGACGCTACTCATTATGTATGTGGACGCCTG
    >sample4
    GAGGTTAACTGGCTCACCTCGCGCGTGTAACAGAGTAATAGGTTGAACGCCTACCCTGGT
    GACCTGGGACGGCGGACGCCTGTTGAAGTAAGGAAACGATCCTAAGCGTCTTGATGGGAT
    CCTATTAATCGGCGCGTGCATATTCATCGGACATGTCGAGGGGTG

Note: the expected pipeline results are also included in the distribution.  To fetch the expected result files::

    copy_snppipeline_data.py lambdaVirusExpectedResults myDirectoryForExpectedResults

After verifying correct results on the lambda data set, you can follow the workflow steps for the agona data
set, :ref:`all-in-one-workflow-agona` or the listeria data set, :ref:`all-in-one-workflow-listeria`.  
To fetch the expected result files::

    copy_snppipeline_data.py agonaExpectedResults myDirectoryForExpectedResults
    copy_snppipeline_data.py listeriaExpectedResults myDirectoryForExpectedResults

**Q: My results for the included test data do not match the expected results. What is the cause?**

A: Different versions of the executable tools can generate different results.  The test data was generated with 
these versions:
	
	* bowtie2 2.2.2
	* samtools 0.1.19
	* varscan 2.3.6

**Q: How can I run the SNP Pipeline with a mix of paired and unpaired samples?**

A: Run the alignSampleToReference script once per sample with either 1 fastq file or 2 fastq files.  
For example::

    alignSampleToReference.sh  reference/NC_011149  samples/CFSAN000448/G0H235M04.RL10.fastq
    alignSampleToReference.sh  reference/NC_011149  samples/CFSAN000449/G00JH2D03.RL11.fastq
    alignSampleToReference.sh  reference/NC_011149  samples/CFSAN000450/HB4DJL101.RL1.fastq
    alignSampleToReference.sh  reference/NC_011149  samples/ERR178930/ERR178930_1.fastq  samples/ERR178930/ERR178930_2.fastq
    alignSampleToReference.sh  reference/NC_011149  samples/ERR178931/ERR178931_1.fastq  samples/ERR178931/ERR178931_2.fastq


**Q: How can I re-run some of the SNP Pipeline processing steps when I see a message that the results are already freshly built?**

A: The SNP Pipeline detects freshly built result files and does not rebuild them.  Result files are
not rebuilt when the file timestamp is newer than all of the input files.  To force a rebuild, 
specfify the ``-f`` option on the command line of any of the tools.  To re-run only some of the steps, 
you can either delete the output files for that step or touch the input files for that step.  All 
subsequent processing steps will also be re-run since their results will be out-of-date.

**Q: How does the SNP Pipeline know which processing steps should be re-run after changing the configuration file?**

A: It doesn't.  If you change the configuration file, you may want to re-run some parts of the pipeline.  The SNP 
Pipeline does not detect which parameters have changed since the last run.  You must manually intervene to cause the 
pipeline to re-run the impacted processing steps.  See the question above for guidance.

.. _faq-performance-label:

Performance
-----------

**Q: How can I control the number of concurrent processes lauched on my workstation?**

A: If you are using a HPC with a job queue manager, the pipeline will automatically run multiple
concurrent processes across multiple servers -- there are no options to control the number of 
concurrent processes.  On a workstation, the pipeline uses all available CPU cores by default
and spawns multiple concurrent processes to use all the cores.  However, you may want to 
control the number of concurrent processes.  There are three steps in the pipeline where multiple
processes are launched on a workstation.  You can control the number of processes with the
following parameters in the configuration file.  These parameters are used only by the
run_snp_pipeline.sh script::

    # Maximum concurrent prepSamples.sh processes (SAMtools and Varscan)
    MaxConcurrentPrepSamples=
    
    # Maximum concurrent create_snp_pileup.py processes
    MaxConcurrentCreateSnpPileup=

    # Maximum concurrent collectSampleMetrics.sh processes
    MaxConcurrentCollectSampleMetrics=

**Q: How can I control the number of CPU cores used by the bowtie2 aligner?**

A: By default, the SNP Pipeline will give bowtie2 all available CPU cores on a workstation and 8 CPU cores per
sample on a high performance computing cluster.  You can override the defaults with the ``-p`` bowtie2 option.  Set
the option either in the configuration file if you are running run_snp_pipeline.sh, or in the Bowtie2Align_ExtraParams 
environment variable if you are running alignSampleToReference.sh directly.  For example, to run alignments with
16 concurrent threads::

    Bowtie2Align_ExtraParams="--reorder -p 16"

On a workstation, alignments are run one at a time using multiple threads per alignment.  On a cluster with
a job queue, multiple alignments are run concurrently, each with multiple threads.

**Q: How can I control the amount of memory that is used by the VarScan java virtual machine?**

A: The amount of memory used by the java VM can be set by using the ``-Xmx`` java VM option.  Set the 
option either in the configuration file if you are running run_snp_pipeline.sh, or in the VarscanJvm_ExtraParams 
environment variable if you are running prepSamples.sh directly. For example, to set maximum java heap 
size to 3000 MB::

    VarscanJvm_ExtraParams="-Xmx3000m"

Developer Questions
-------------------

**Q: What causes "ImportError: No module named sphinx_rtd_theme" when building the documentation?**

A: The documentation uses the *Read The Docs* theme.  Install it like this::

	$ pip install --user sphinx_rtd_theme

**Q: I installed sphinx_rtd_theme, but I still get error "ImportError: No module named sphinx_rtd_theme".**

A: Try running sphinx like this::

	$ python /usr/bin/sphinx-build -b html  .  ./_build

**Q: I changed one of the shell scripts, but the changes are ignored.**

A: Reinstall the distribution.  Do this::

	$ python setup.py develop