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
jar file, but not execute-access.  Set the CLASSPATH environment variable to
specify the path to the VarScan jar file.  The CLASSPATH should include the filename, not just the directory.
Define it like this in your .bashrc file::

    export CLASSPATH=~/software/varscan.v2.3.9/VarScan.v2.3.9.jar:$CLASSPATH


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

A: Make sure you have the proper dependencies on your path.  Modify your path if necessary to include bowtie,
samtools, and bcftools.  See the question above about installing VarScan.
See the :ref:`installation-label` section of this documentation.

**Q: How can I verify the pipeline is installed and working properly?**

A: The SNP Pipeline includes sets of test data with result files.  You can run the pipeline against the
test data to verify correct results.  Follow the lambda virus workflow steps here: :ref:`all-in-one-workflow-lambda`.

Upon successful completion of the pipeline, the snplist.txt file should have 166 entries.  The SNP Matrix
can be found in snpma.fasta.  It should have the following contents::

    >sample1
    AGCACCGGGACCCACGGCGCACGCAAAGATCCGAATTGCAGGGCGTACCTGGACCCCGGT
    GACGGGGGATCGGGGACTCTTGGTCGAGGAACTAAAACGAACATCCACGTTTTCATGGCG
    ACTGCTTGCCAGGTGTCAGCACATTCCCTATATCGGTGGACACGTA
    >sample2
    GGCGCTAGGAGGCAAGCCTTGGTCGTGGTTAATAGTTACAAGGCGTGCGCGTACTGCCGT
    CTCCTACTATCTCTGCCGCCTCTCCGCGATCCGGACCGCAACACCAACTCTCTGGTGGCA
    TCCTCTGAATCGTCGTGAGCATCTCAATTATATATTCGTCCGCGCG
    >sample3
    GGCGCTAGGAGGTACGCTTCGCGTGTGGATCAGCGCTACGGTGCCTATGCGTGACCCGCG
    GAACTGGGTTCGCGTAAGGCAGTTTCAGGTACGGCAACGTAGATCAAAGT-TAGAAACCA
    TACTCGTAATCCGCCTGACGCTACTCATTATGTATGTGGACGCCTG
    >sample4
    GAGGTTAACTGGCTCACCTCGCGCGTGTAACAGAGTAATAGGTTGAACGCCTACCCTGGT
    GACCTGGGACGGCGGACGCCTGTTCGAAGTAAGGAAACGATCCTAAGCGTCTTGATGGGA
    TCCTATTAATCGGCGCGTGCATATTCATCGGACATGTCGAGGGGTG

Note: the expected pipeline results are also included in the distribution.  To fetch the expected result files::

    cfsan_snp_pipeline data lambdaVirusExpectedResults myDirectoryForExpectedResults

After verifying correct results on the lambda data set, you can follow the workflow steps for the agona data
set, :ref:`all-in-one-workflow-agona` or the listeria data set, :ref:`all-in-one-workflow-listeria`.
To fetch the expected result files::

    cfsan_snp_pipeline data agonaExpectedResults myDirectoryForExpectedResults
    cfsan_snp_pipeline data listeriaExpectedResults myDirectoryForExpectedResults

**Q: My results for the included test data do not match the expected results. What is the cause?**

A: Different versions of the executable tools can generate different results.  The test data was generated with
these versions:

	* bowtie2 2.3.4.1
	* samtools 1.8
	* varscan 2.3.9

**Q: How can I run the SNP Pipeline with a mix of paired and unpaired samples?**

A: This is handled automatically if you use the ``run`` command.  If you are running the ``map_reads`` command,
run the script once per sample with either 1 fastq file or 2 fastq files.
For example::

    cfsan_snp_pipeline map_reads  reference/NC_011149  samples/CFSAN000448/G0H235M04.RL10.fastq
    cfsan_snp_pipeline map_reads  reference/NC_011149  samples/CFSAN000449/G00JH2D03.RL11.fastq
    cfsan_snp_pipeline map_reads  reference/NC_011149  samples/CFSAN000450/HB4DJL101.RL1.fastq
    cfsan_snp_pipeline map_reads  reference/NC_011149  samples/ERR178930/ERR178930_1.fastq  samples/ERR178930/ERR178930_2.fastq
    cfsan_snp_pipeline map_reads  reference/NC_011149  samples/ERR178931/ERR178931_1.fastq  samples/ERR178931/ERR178931_2.fastq


**Q: How can I re-run some of the SNP Pipeline processing steps when I see a message that the results are already freshly built?**

A: The SNP Pipeline detects freshly built result files and does not rebuild them.  Result files are
not rebuilt when the file timestamp is newer than all of the input files.  To force a rebuild,
specify the ``-f`` option on the command line of any of the tools.  To re-run only some of the steps,
you can either delete the output files for that step or touch the input files for that step.  All
subsequent processing steps will also be re-run since their results will be out-of-date.

**Q: How does the SNP Pipeline know which processing steps should be re-run after changing the configuration file?**

A: It doesn't.  If you change the configuration file, you may want to re-run some parts of the pipeline.  The SNP
Pipeline does not detect which parameters have changed since the last run.  You must manually intervene to cause the
pipeline to re-run the impacted processing steps.  See the question above for guidance.


**Q: What do the dashes (“-“) in the snp matrix indicate?**

A: Gaps, “-“, are either missing bases (indels) or cases where there is insufficient information to make a consensus call
(coverage depth too low, or consensus base frequency too low).

**Q: Why are some snps missing from the snp matrix even when the snps were called by VarScan?**

A: Older versions of VarScan failed to generate the header section of some VCF files.  This in turn, caused the SNP Pipeline
to ignore the first snp in the VCF file.  Upgrade to a newer version VarScan.

.. _optical-dup-read-label:

**Q: Why are there no optical duplicate reads and why am I seeing the warning message "Default READ_NAME_REGEX '<optimized capture of last three ':' separated fields as numeric values>' did not match read name"?**

A: First, this is not a serious problem -- optical duplicate reads occur much less frequently than PCR amplification duplicates.
This message appears in the log file when Picard MarkDuplicates cannot identify the tile number, x-position, and y-position in
the read names in the BAM file.  Without those data elements, optical duplicates cannot be identified.  You will see this warning
only once, but it's usually a problem for every read in the file.  When downloading fastq files from NCBI with ``fastq-dump``,
you can specify the ``--origfmt`` option to format the read names in the original Illumina format when possible.  It is not always
possible because NCBI does not always store the original read names in the SRA database.

.. _faq-performance-label:

Performance
-----------

**Q: How can I limit the CPU resources consumed by the pipeline?**

A: By default, the pipeline will use all available CPU resources.  You can limit the number of CPU cores the
pipeline will use with the :ref:`MaxCpuCores-label` parameter in the configuration file.  This works on your workstation and
also on a high performance computing cluster running Grid Engine or Torque.  See also the questions below.

**Q: How can I control the number of CPU cores used by the bowtie2, smalt, samtools, and GATK?**

A: You can set the number of CPU cores per process with the :ref:`CpuCoresPerProcessOnHPC-label` and :ref:`CpuCoresPerProcessOnWorkstation-label`
parameters.  This is the recommended way to control CPU cores per process for all multi-threaded processes in the pipeline::

    CpuCoresPerProcessOnHPC = 16
    CpuCoresPerProcessOnWorkstation = 8

If you want to allocate a different number of CPU cores for different processes, you can customize
each command.  For example, you can set the number of bowtie threads with the ``-p`` option and
the number of samtools threads with the ``-@`` option.  Set the options either in the configuration
file if you are using the ``run`` command, or in the environment variables if you are running the
map_reads command directly.  ::

    Bowtie2Align_ExtraParams = "-p 20"
    SamtoolsSort_ExtraParams = "-@ 15"

You cannot use more threads than the number of allowed CPU cores.  For example, if you set bowtie to use
10 threads and you set MaxCpuCores to 8, bowtie will only get 8 threads, not 10.

**Q: How can I control the number of concurrent processes?**

A: You can set configuration parameters to indirectly control the number of processes.  First, set the :ref:`MaxCpuCores-label` parameter.
Then set the :ref:`CpuCoresPerProcessOnHPC-label` and :ref:`CpuCoresPerProcessOnWorkstation-label` parameters.  The number of concurrently
executing processes will be the total number of CPUs divided by the number of CPU cores per process.

Multiple processes are run concurrently, each with multiple threads.  The number of concurrently executing
processes depends on the :ref:`MaxCpuCores-label` parameter setting.  For example if you set bowtie to use
10 CPU cores per process and set MaxCpuCores to 20, you will have 2 concurrent bowtie2 processes each with
10 CPU threads.

**Q: How can I control the amount of memory that is used by the Picard, GATK, and VarScan java virtual machines?**

A: The amount of memory used by the java VM can be set by using the ``-Xmx`` java VM option.  Set the
option either in the configuration file if you are using the ``run`` command, or in the VarscanJvm_ExtraParams
environment variable if you are running the call_sites command directly. For example, to set maximum java heap
size to 3000 MB::

    PicardJvm_ExtraParams="-Xmx3000m"
    GatkJvm_ExtraParams="-Xmx3000m"
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
