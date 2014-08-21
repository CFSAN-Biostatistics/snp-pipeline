===========================
FAQ / Troubleshooting Guide
===========================

.. highlight:: bash

Installation
------------

**Q: How can I avoid polluting my global python installation when installing the SNP pipeline?**


A: You can either use a python virtual environment or install into your user area.  To use a python virtual 
environment, see :ref:`contributing-label`.  To install into your user area instead of installing into your 
global site packages, do this::

	$ pip install --user snp-pipeline

**Q: The SNP Pipeline cannot find VarScan.  How should I install it?**

A: Download the VarScan jar file from SourceForge.  Put the jar file anywhere.  You need read-access to the
jar file, but not execute-access.  The suppled shell scripts expect the CLASSPATH environment variable to 
specify the path to the VarScan jar file.  The CLASSPATH should include the filename, not just the directory.
Define it like this in your .bashrc file::

    export CLASSPATH=~/software/varscan.v2.3.6/VarScan.v2.3.6.jar:$CLASSPATH



Running the Pipeline
--------------------

**Q: Nothing works.**

A: Make sure you have the proper dependencies on your path.  Modify your path if necessary to include bowtie 
and samtools.  See the question above about installing VarScan.  Also make sure you install Biopython.  If 
you are using virtual environments, make sure you issue the command "toggleglobalsitepackages" so Biopython 
can be found.  See also: the :ref:`installation-label` section of this documentation.

**Q: How can I verify the pipeline is installed and working properly?**

A: The SNP Pipeline includes a small set of test data with result files.  You can run the pipeline against the 
test data to verify correct results.  Follow the lambda virus workflow steps in the :ref:`usage-label` section.

Upon successful completion of the pipeline, the snplist.txt file should have 165 entries.  The SNP Matrix 
can be found in snpma.fasta.  It should have the following contents::

    >sample1 <unknown description>
    AGCACCGGGACCCACGGCGCACGCAAAGATCCGAATTGCAGGGCGTACCTGGACCCCGGT
    GACGGGGGATCGGGGACTCTTGGTGAGGAACTAAAACGAACATCCACGTTTTCATGGCGA
    CTGCTTGCCAGGTGTCAGCACATTCCCTATATCGGTGGACACGTA
    >sample2 <unknown description>
    GGCGCTAGGAGGCAAGCCTTGGTCGTGGTTAATAGTTACAAGGCGTGCGCGTACTGCCGT
    CTCCTACTATCTCTGCCGCCTCTCGCGATCCGGACCGCAACACCAACTCTCTGGTGGCAT
    CCTCTGAATCGTCGTGAGCATCTCAATTATATATTCGTCCGCGCG
    >sample3 <unknown description>
    GGCGCTAGGAGGTACGCTTCGCGTGTGGATCAGCGCTACGGTGCCTATGCGTGACCCGCG
    GAACTGGGTTCGCGTAAGGCAGTTCAGGTACGGCAACGTAGATCAAAGTTTAGAAACCAT
    ACTCGTAATCCGCCTGACGCTACTCATTATGTATGTGGACGCCTG
    >sample4 <unknown description>
    GAGGTTAACTGGCTCACCTCGCGCGTGTAACAGAGTAATAGGTTGAACGCCTACCCTGGT
    GACCTGGGACGGCGGACGCCTGTTGAAGTAAGGAAACGATCCTAAGCGTCTTGATGGGAT
    CCTATTAATCGGCGCGTGCATATTCATCGGACATGTCGAGGGGTG

Note: the expected pipeline results are also included in the distribution.  To fetch the expected result files::

    copy_snppipeline_data.py lambdaVirusExpectedResults myDirectoryForExpectedResults

**Q: My results for the included test data do not match the expected results. What is the cause?**

A: Different versions of the executable tools can generate different results.  The test data was generated with 
these versions:
	
	* bowtie2 2.2.2
	* samtools 0.1.19
	* varscan 2.3.6

**Q: How can I run the SNP Pipeline with a mix of paired and unpaired samples?**

A: Run the alignSampleToReference script once per sample with either 1 fastq file or 2 fastq files.  See 
the following example::

    alignSampleToReference.sh  8  reference/NC_011149  samples/CFSAN000448/G0H235M04.RL10.fastq
    alignSampleToReference.sh  8  reference/NC_011149  samples/CFSAN000449/G00JH2D03.RL11.fastq
    alignSampleToReference.sh  8  reference/NC_011149  samples/CFSAN000450/HB4DJL101.RL1.fastq
    alignSampleToReference.sh  8  reference/NC_011149  samples/ERR178930/ERR178930_1.fastq  samples/ERR178930/ERR178930_2.fastq
    alignSampleToReference.sh  8  reference/NC_011149  samples/ERR178931/ERR178931_1.fastq  samples/ERR178931/ERR178931_2.fastq



Developer Questions
-------------------

**Q: What causes "ImportError: No module named sphinx_rtd_theme" when building the documentation?**

A: The documentation uses the *Read The Docs* theme.  Install it like this::

	$ pip install sphinx_rtd_theme

**Q: I installed sphinx_rtd_theme, but I still get error "ImportError: No module named sphinx_rtd_theme".**

A: Try running sphinx like this::

	$ python /usr/bin/sphinx-build -b html  .  ./_build

**Q: I changed one of the shell scripts, but the changes are ignored.**

A: Reinstall the distribution.  Do this::

	$ python setup.py develop