.. _reproducible-label:

====================
Reproducible Results
====================

This page is a work in progress.

We strive to make the SNP Pipeline results fully reproducible -- not just the
final SNP matrix, but each intermediate file as well.  Reproducible results
help us test and debug the pipeline and also facilitate collaborative efforts
between researchers.

Version Control
---------------
The SNP Pipeline source code is available on GitHub so anyone can download
and reproduce our results.  Each new release is tagged with a version identifier
and also released to the Python Package Index for easy installation.

Control Test Set
----------------
* We should probably include a control data set for which the SNPs are known.
* We have the lambda virus data set, but we don't know the correct list of snps.

Parameters
----------
* Debating whether to include this section...
* Parameter values can affect results.
* Some parameters are exposed, some are hard-coded.
* We recommend changing one parameter at a time in a controlled manner to determine the effects of the parameter change.
* We recommend accurately recording the parameter values used for any important results, ideally in a script or configuration file.

Concurrency
-----------
The SNP Pipeline takes advantage of multiple CPU cores to run portions of the processing 
in parallel.  However, concurrency can lead to non-deterministic behavior and different 
results when the pipeline is run repeatedly.  The pipeline addresses known concurrency issues
with bowtie and samtools.

**Bowtie**

The SNP Pipeline uses multiple CPU cores during the bowtie alignment.  Unless told otherwise, 
when bowtie runs multiple concurrent threads, it generates output records in the SAM file in 
non-deterministic order.  The consequence of this is the SAM files and Pileup 
files can differ between runs.  This may appear as two adjacent read-bases swapped in the 
pileup files.

To workaround this problem, the pipeline uses the ``--reorder`` bowtie command line option.  
The reorder option causes bowtie to generate output records in the same order as the
reads in the input file.  This is discussed in the bowtie documentation here:
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#performance-options

**SAMtools**

The SNP Pipeline runs multiple samtools processes concurrently to generate pileups for
each sample.  When the samtools pileup process runs, it checks for the existence of the
reference faidx file, ``*.fai``.  If the faidx file does not exist, samtools creates it 
automatically.  However, multiple samtools mpileup processes can interfere with each 
other by attempting to create the file at the same time.  This interference causes 
incorrect pipeline results.

To workaround this problem, the pipeline explicitly creates the faidx file by running 
``samtools faidx`` on the reference before running the mpileup processes.  This prevents errors 
later when multiple samtools mpileup processes run concurrently.
  

Software Versions
-----------------
Different versions of the executable tools can generate different results.  This may be important
if you intend to compare the results between runs.

**Bowtie**

Bowtie 2.1.0 and 2.2.2 produce functionally identical SAM files.  The only difference
observed is a header record in bowtie 2.2.2 documenting the program version and command line parameters.

**SAMtools**

Different versions of SAMtools can produce different pileup files which can subsequently impact
the snp list and snp matrix.  SAMtools version 0.1.19 has the capability to exclude read bases 
with low quality. By default, it excludes read bases with quality below 13 (95% accuracy).
SAMtools mpileup version 0.1.18 does not filter read bases by quality and there is no option to 
change this behavior.

On one of our data sets with 116 samples, we observed these results:

* 36030 snps found when pileups generated with SAMtools 0.1.18
* 38154 snps found when pileups generated with SAMtools 0.1.19

