.. _reproducible-label:

================================
Correct and Reproducible Results
================================


It is our goal to make the SNP Pipeline results both fully
reproducible and as correct as current scientific understanding
allows. As part of this effort, we document here problems we have
found that have affected correctness. We also detail how we have built
this software so the results are as reproducible as possible. In addition,
we are building, collecting, and collating data sets that we use to
assess both correctness and reproducibility of our software. As a
project that is made possible only by very recent developments in
science and technology, our efforts to ensure correctness and
reproducibility are an ongoing effort. The publications we have
produced in an effort to ensure scientific correctness are listed as
references at the bottom of this document. This document will continue
to evolve as we improve our process and as scientific advances occur.

Reproducible Results
====================

We have made the SNP Pipeline results fully reproducible -- not just the
final SNP matrix, but each intermediate file as well.  Reproducible results
help us test and debug the pipeline and also facilitate collaborative efforts
between researchers.

Public Availability
-------------------
The SNP Pipeline source code is available on GitHub so anyone can download our
source code. The GitHub repository contains some data sets that can be used to
reproduce selected results. We also provide information on how to obtain and
verify other data sets that we have used. (These data sets are large, so we
do not provide them directly.)

Version Control
---------------
We use git internally for code development to ensure that we have control over
our source code and can identify which version of code was used to produce any
particular result. We tag/version commits of the code that we consider production
releases and use them for the majority of our internal analyses. We also release
each of these tagged versions to GitHub and to the Python Package Index for easy
installation.

Parameters
----------
The SNP pipeline behavior depends on the setting of a number of parameters that
determine the behavior of various software packages that the pipeline uses. These
parameters affect both the correctness and reproducibility of results. We have set
all the parameters so that the results are reproducible. This entails
setting seeds for all random number dependent processes, as well as specific
choices for other parameters that can affect such behavior as the order of the
results. We discuss these aspects of ensuring reproducibility in more detail in
other portions of this document.

The pipeline depends on some fairly complex software packages, and these packages have
large numbers of parameters. We do not expose all of these parameters, but only those
we have found it useful to modify in our work. Those wishing to further modify the
software behavior will have to adjust the code to meet their needs.

We recommend recording the parameter values used for any important
results, ideally in a script or configuration file that is under
version control. Future versions of this software will probably
contain both a configuration file that sets all parameters, as well as
an output file documenting the run.


Concurrency
-----------
The SNP Pipeline takes advantage of multiple CPU cores to run portions of the
processing in parallel.  However, concurrency can lead to non-deterministic behavior
and different results when the pipeline is run repeatedly.  The pipeline addresses
known concurrency issues with bowtie and samtools.

**Bowtie**

The SNP Pipeline uses multiple CPU cores during the bowtie alignment.  Unless told
otherwise, when bowtie runs multiple concurrent threads, it generates output records
in the SAM file in non-deterministic order.  The consequence of this is the SAM
files and Pileup files can differ between runs.  This may appear as two adjacent
read-bases swapped in the pileup files.

To work around this problem, the pipeline uses the ``--reorder`` bowtie command line
option. The reorder option causes bowtie to generate output records in the same order
as the reads in the input file.  This is discussed in the bowtie documentation here:
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#performance-options

**SAMtools**

The SNP Pipeline runs multiple samtools processes concurrently to generate
pileups for each sample.  When the samtools pileup process runs, it checks for the
existence of the reference faidx file, ``*.fai``.  If the faidx file does not exist,
samtools creates it automatically.  However, multiple samtools mpileup processes can
interfere with each other by attempting to create the file at the same time.  This
interference causes incorrect pipeline results.

To work around this problem, the pipeline explicitly creates the faidx
file by running ``samtools faidx`` on the reference before running the
mpileup processes.  This prevents errors later when multiple samtools
mpileup processes run concurrently.


Software Versions
-----------------

Different versions of the software packages this pipeline uses can
generate different results.  This is important to be aware of if you
end up comparing the results between runs. We share our observations
from the versions of SAMtools and Bowtie that we have used below.

**Bowtie**

Bowtie 2.1.0 and 2.2.2 produce functionally identical SAM files.  The
only difference we observed is a header record in bowtie 2.2.2
documenting the program version and command line parameters.

**SAMtools**

SAMtools mpileup version 0.1.18 and version 0.1.19 differ in their
default behavior. Version 0.1.19 can filter out bases with low
quality, and by default, it excludes bases with quality score below 13
(95% accuracy). Version 0.1.18 does not have this capability, and thus
different versions of SAMtools mpileup when run with the default
parameters can produce different pileup files which can impact the snp
list and snp matrix.

On one of our data sets with 116 samples, we observed these results:

* 36030 snps found when pileups generated with SAMtools 0.1.18
* 38154 snps found when pileups generated with SAMtools 0.1.19

Correct Results
===============

As we have constructed our pipeline, we have found problems both in
our own software and in the various packages we use. To this point we
have found two problems worth mentioning here.

**SAMtools first and last base SNPs not called**

TODO

**SAMtools snp pileup difference from genome-wide pileup**

An important processing step in the SNP Pipeline is creation of a pileup 
file per sample containing read pileups at the positions where snps were called 
in *any* of the samples.  This pileup file should be a subset of the genome-wide
pileup for each sample.  However, the SAMtools software does not generate 
pileup records exactly matching the genome-wide pileup when given a list of 
positions for which the pileup should be generated.  The differences are 
particularly evident at the first few snp positions and cause missing
values in the SNP matrix.  To work around this problem, the SNP Pipeline 
internally extracts the desired pileup records from the genome-wide pileup.

Test Data Sets
==============

We have created/curated a number of data sets for use in testing both the
reproducibility and correctness of the pipeline. In the following sections
we briefly describe theses data sets.

Lambda Virus
------------

This data set was built using the bowtie2 example, and intended to be a small
test case and example that will run quickly and verify the basic functionality
of the code.

Salmonella Agona
----------------

This data set was designed to contain realistic sequences, but not very many
of them, so that it could be run in a reasonable amount of time. The data must
be downloaded from the NCBI due to its large size. We provide a file of hashes
that can easily be used to verify that the data downloaded matches the data
originally used to produce our results. (Use sha256sum at the unix command line.)

Listeria monocytogenes
----------------------

*Coming soon in a future release*

This is designed to be a realistic-sized data set based on an outbreak
of L. m.  in Roos cheese. The data must be downloaded from the NCBI
due to its large size. We provide a file of hashes that can easily be
used to verify that the data downloaded matches the data originally
used to produce our results. (Use sha256sum at the unix command line.)

Synthetic data sets
-------------------

*Coming soon in a future release*

We are currently creating synthetic data sets based on simulating
various evolutionary scenarios. The simulations are designed to be
similar to what we would expect in the types of organisms we study
(food-borne pathogens), with error structure appropriate for the
platforms we use to do sequencing.

