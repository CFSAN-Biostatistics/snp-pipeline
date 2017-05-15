.. :changelog:

History
-------

0.8.0 (2017-05-09) - `docs <http://snp-pipeline.readthedocs.io/en/0.8-branch/history.html>`_ 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Changes Impacting Backwards Compatibility:**

* Changed the collectSampleMetrics script to only accept input files in the sample directory,
  not in arbitrary locations.
* Changed the combineSampleMetrics script to write to metrics.tsv by default, not stdout.
* Leading zeros are stripped from Miseq flowcell identifiers in the metrics files.
* Added a dependency on Picard.  You need to install Picard and change your CLASSPATH.
  See :ref:`installation-label`.
* Removed the unused create_snp_pileup.py script.

**Bug Fixes:**

* Fixed the machine and flow cell reporting in the metrics file when the fastq read names are not
  in the original Illumina format.
* Fixed the calculation of average pileup depth in the metrics file.  The formula previously
  included whitespace characters when calculating the length of the reference.  The correct
  average depth is slightly deeper than previously calculated.

**Other Changes:**

* Sweeping changes under the hood replacing most shell scripts with equivalent python code.
  Repackaged the SNP Pipeline as a single executable with multiple sub-commands.  The old scripts
  still exist for backwards compatibility and are rewritten as one-liners calling the new
  replacement commands.  The main executable program is called :ref:`cmd-ref-cfsan-snp-pipeline`.
* Added the capability to remove duplicate reads from BAM files prior to creating the pileup and
  calling snps.  See :ref:`remove-duplicate-reads-label`.  This change introduces a dependency on
  ``Picard`` and will require changing your CLASSPATH.  See :ref:`installation-label`. You can
  disable this step and keep the duplicate reads by configuring ``SnpPipeline_RemoveDuplicateReads=false``
  in the configuration file.
* Added a new metric to count the number of duplicate reads in each sample.
* Capture read-group metadata in the SAM/BAM files during the read mapping step.
* Added a new configuration parameter, ``BcftoolsMerge_ExtraParams`` to allow customizing the
  snpma.vcf files created when merging the consensus VCF files.  See :ref:`configuration-label`.
* Removed the hard-coded wall-clock run-time limits for Torque and Sun Grid Engine jobs.  Added
  default limits (12 hours) to the configuration file.  You can change the runtime limits for
  all SNP Pipeline job steps with the ``Torque_QsubExtraParams`` or ``GridEngine_QsubExtraParams``
  configuration parameters.
* Log the SNP Pipeline version in the header of all the log files.
* Changed the composition of the included Salmonella Agona data set to remove the excessively large
  sample ERR178930 and include a more diverse set of isolates from different geographic locations,
  different environmental sources, and different types of sequencing instruments.


0.7.0 (2016-11-30) - `docs <http://snp-pipeline.readthedocs.io/en/0.7-branch/history.html>`_ 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Added a new script to the pipeline: ``snp_filter.py`` removes snps from the ends of contigs and
  from regions where the snp density is abnormally high.  This is an important change to the
  pipeline with additional processing and new output files.  See :ref:`snp-filtering-label`.
* NOTE: You cannot re-use an old configuration file when running SNP Pipeline version 0.7.0.  You
  must create a new configuration file.  See :ref:`configuration-label`.
* Fixed compatibility with bcftools 1.2 and higher.
* Updated the result files in the included data sets with the results obtained using bcftools v1.3.1
  and bowtie2 v2.2.9.  Note: upgrading from bowtie 2.2.2 to 2.2.9 did not change the snp matrix
  on any of the included datasets.


0.6.1 (2016-05-23)
~~~~~~~~~~~~~~~~~~

* Fixed compatibility with SAMtools 1.3.
* Changed the expected results data sets to match the results obtained using SAMtools
  version 1.3.1.  Starting with SAMtools version 1.0, the samtools mpileup command implemented
  a feature to avoid double counting the read depth when the two ends of a paired-end read
  overlap.  If you use this feature of SAMtools, the pileup depth will be noticably reduced.
  You can still count the overlapping read sections twice by using SAMtools v0.1.19 or by using
  a configuration file specifying the ``-x`` option in ``SamtoolsMpileup_ExtraParams``.
* Removed the obsolete ``reads.snp.pileup`` files from the included results data sets.

0.6.0 (2016-04-11)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Fixed compatibility with the newly released PyVCF 0.6.8 package.

**Other Changes:**

* A new configuration parameter, ``SnpPipeline_MaxSnps``, controls the maximum number of snps
  allowed for each sample.  Samples with excessive snps exceeding this limit are excluded
  from the snp list and snp matrix.
  See :ref:`excessive-snps-label`.
* A new column in the metrics.tsv file, ``Excluded_Sample``, indicates when a sample has been
  excluded from the snp matrix.  This column is normally blank.
* Added a new script to the pipeline: ``calculate_snp_distances.py`` computes the SNP distances between
  all pairs of samples. The SNP distances are written to the output files ``snp_distance_pairwise.tsv``
  and ``snp_distance_matrix.tsv``.
* Changed Sun Grid Engine execution to use array-slot dependency where possible, resulting
  in less idle time waiting for job steps to complete.


0.5.2 (2016-03-07)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* An empty snplist.txt file should not cause errors when creating the referenceSNP.fasta.
* An empty snplist.txt file should not preclude re-running subsequent steps of the pipeline.
* When configured to ignore single-sample errors, a missing var.flt.vcf file should not
  preclude rebuilding the snplist.txt file during a pipeline re-run.
* The metrics file did not properly capture the total number of snps per sample. See below for the details.

**Other Changes:**

* Capture separate metrics counting phase 1 snps (varscan) and phase 2 snps (consensus). Previously, the
  metrics only included phase 1 snps.  This changes the contents of both the ``metrics`` and ``metrics.tsv``
  files. The metrics file now contains a new tag ``phase1Snps``.  The old tag ``snps`` now correctly counts
  the total number of snps. The metrics.tsv file now has separate column headers for phase 1 snps and
  phase 2 snps.  Any code that parses those files may need modifications to work properly with v0.5.2.
* Added the ``Average Insert Size`` metric.
* The metrics.tsv column headings now contain underscores instead of spaces for better interoperability
  with some downstream analysis tools. Column headings with spaces can be generated by specifing the
  combineSampleMetrics.sh ``-s`` option in the configuration file.
* Remove the dependence on the snp matrix when collecting sample metrics.
* Improve the speed of metrics calculation when rerunning the pipeline.  Reuse the previously computed metrics
  when recalculation would be slow.


0.5.1 (2016-02-19)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Do not shutdown the pipeline when the generated snplist is empty when there are no snps.
* Do not attempt to merge VCF files when there are fewer than two VCF files to merge.

**Other Changes:**

* Added the ``vcfFailedSnpGt`` option to the call_consensus.py script to control how the VCF file GT data
  element is emitted when the snp is failed because of depth, allele frequency, or some other filter.  If
  not specified, the GT element will contain a dot.  Prior to this release, the behavior was to emit the
  ALT allele index.  The old behavior can be retained by setting ``--vcfFailedSnpGt 1``
* Changed the setup to require PyVCF version 0.6.7 or higher.  It will automatically upgrade if necessary.
* Added error checking after running SamTools and VarScan to detect missing, empty, or erroneous output files.


0.5.0 (2016-01-19)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Changed VCF file generator to not emit multiple alleles when the reference base is lowercase.

**Other Changes:**

* Trap errors, shutdown the pipeline, and prevent execution of subsequent steps when earlier processing
  steps fail. A summary of errors is written to the ``error.log`` file.
  See :ref:`error-handling-label`.
* Check for the necessary software tools (bowtie, samtools, etc.) on the path at the start of each
  pipeline run.
* Check for missing or empty input files at the start of each processing step.
* Added two new parameters, ``GridEngine_QsubExtraParams`` and ``Torque_QsubExtraParams``, to the
  configuration file to pass options to qsub when running the SNP Pipeline on an HPC computing cluster.
  Among other things, you can control which queue the snp-pipeline will use when executing on an HPC
  with multiple queues.  See :ref:`configuration-label`.
* Removed the "job." prefix to shorten job names when running on an HPC.
* Changed the vcf file generator to emit reference bases in uppercase.  Added the ``vcfPreserveRefCase``
  flag to the call_consensus.py script to cause the vcf file generator to emit each reference base in
  uppercase/lowercase as it appears in the original reference sequence file.  If not specified, the
  reference bases are emitted in uppercase.  Prior to this release, the behavior was to always preserve the
  original case.
* Added support for Python 3.3, 3.4, 3.5.
* Implemented a regression test suite for the bash shell scripts, using the shUnit2 package.


0.4.1 (2015-10-30)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Fixed a Python 2.6 incompatibility with the new consensus caller.

**Other Changes:**

* Added Tox support for automatically testing installation and execution with multiple Python versions.


0.4.0 (2015-10-22)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* When run on Grid Engine with the default settings, bowtie2 was consuming all available CPU cores
  per node while scheduled with Grid to use only 8 cores. On a lightly loaded cluster, this bug made
  the pipeline run faster, but when the cluster was full or nearly full, it would cause contention
  for available CPU resources and cause jobs to run more slowly.  Changed to use only 8 CPU cores
  by default.
* The consensus snp caller miscounted the number of reference bases when the pileup record
  contained the ^ symbol marking the start of a read segment followed by a dot or comma.  In this
  situation, the dot or comma should not be counted as reference bases.


**Other Changes:**

* Added support for the Smalt aligner.  You can choose either bowtie2 or smalt in the configuration file.
  A new parameter in the configuration file, ``SnpPipeline_Aligner``, selects the aligner to use.
  Two additional configuration parameters, ``SmaltIndex_ExtraParams`` and ``SmaltAlign_ExtraParams``
  can be configured with any Smalt command line options.  See :ref:`tool-selection-label`.  The
  default aligner is still bowtie2.
* Split the create_snp_matrix.py script into two pieces.  The new script, call_consensus.py, is a redesigned
  consensus caller which is run in parallel to call snps for multiple samples concurrently.  The
  create_snp_matrix.py script simply merges the consensus calls for all samples into a multi-fasta file.
* The new consensus caller has the following adjustable parameters.
  See the :ref:`cmd-ref-call-consensus` command reference.

  * ``minBaseQual`` : Mimimum base quality score to count a read.
  * ``minConsFreq`` : Minimum consensus frequency.
  * ``minConsStrdDpth`` : Minimum consensus-supporting strand depth.
  * ``minConsStrdBias``: Strand bias.
* Added the capability to generate VCF files.  By default, a file named consensus.vcf is generated
  by the consensus caller for each sample, and the merged multi-sample VCF file is called snpma.vcf.
  This capability introduces a new dependency on bgzip, tabix, and bcftools.  You can disable VCF file
  generation by removing the ``--vcfFileName`` option in the configuration file. Also, be aware the
  contents of the VCF files may change in future versions of the SNP Pipeline.
* Added configuration parameters ``Torque_StripJobArraySuffix`` and ``GridEngine_StripJobArraySuffix`` to
  improve compatibility with some HPC environments where array job id suffix stripping is
  incompatible with qsub.
* Renamed the configuration parameter ``PEname`` to ``GridEngine_PEname``.

0.3.4 (2015-06-25)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* The referenceSNP.fasta file was missing newlines between sequences when the reference fasta file
  contained multiple sequences.  In addition, each sequence was written as a single long string of
  characters.  Changed to emit a valid fasta file.  Updated the expected result files for the
  datasets included with the distribution accordingly.
* Changed the run_snp_pipeline.sh script to allow blank lines in the file of sample directories
  when called with the -S option.
* Changed the run_snp_pipeline.sh script to allow trailing slashes in the file of sample directories
  when called with the -S option.
* Do not print system environment information when the user only requests command line help.
* Fixed the broken pypi downloads per month badge on the readme page.

**Other Changes:**

* Changed the default configuration file to specify the ``-X 1000`` option to the bowtie2 aligner.  This
  parameter is the maximum inter-mate distance (as measured from the furthest extremes of the mates)
  for valid concordant paired-end alignments.  Previously this value was not explicitly set and
  defaulted to 500.  As a result of this change, the generated SAM files may have a different number
  of mapped reads, the pileup files may have different depth, and the number of snps called may change.
* We now recommend using VarScan version 2.3.9 or later.  We discoved VarScan v2.3.6 was occasionally
  omitting the header section of the generated VCF files.  This in turn, caused the SNP Pipeline
  to miss the first snp in the VCF file.  This is not a SNP Pipeline code change, only a
  documentation and procedural change.
* Updated the result files in the included data sets with the results obtained using VarScan v2.3.9
  and the Bowtie -X 1000 option.
* Log the Java classpath to help determine which version of VarScan is executed.
* Changed the python unit tests to execute the non-python processes in a temporary directory instead
  of assuming the processes were already run in the test directory.



0.3.3 (2015-04-14)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Improve HPC qsub submission speed throttling to avoid errors with the HPC job scheduler when
  submitting large and small jobs.  Dynamically adjust the delays between HPC array job submission so
  small datasets have small delays and large datasets have large delays between qsub submissions.
* Process the sample directories in order by size, largest first, considering only the size of fastq
  files and ignoring all other files.  Previously non-fastq files were affecting the processing order.
* Fixed divide-by-zero error in create_snp_matrix when no snps are detected.
* Don't skip the last sample when run_snp_pipeline is started with the -S option and the file of
  sample directories is not terminated with a newline.
* Gracefully exit run_snp_pipeline with error messages when run with -S option and any of the sample
  directories in the sample directory file is missing, empty, or does not contain fastq files.
* Gracefully exit run_snp_pipeline with an error message when run with -s option and the samples directory
  is empty or contains no subdirectories with fastq files.
* Fixed the sun grid engine "undefined" task id reported in non-array job log files.

**Other Changes:**

* Sample Metrics.  The pipeline generates a table of sample metrics capturing various alignment, coverage, and snp statistics per sample.
  See :ref:`metrics-usage-label`.
* Explicitly expose the ``minConsFreq`` parameter in the supplied default configuration file to make it easier to adjust.
* Updated the FAQ with instructions to install to an older version.



0.3.2 (2015-01-14)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Fixed (again) a Python 2.6 incompatibility with formatting syntax when printing the available RAM.
  This affected the shell scripts (prepReference.sh, alignSampleToReference.sh, prepSamples.sh).
* Improved installation in a Python 2.6 environment.  Added several Python packages to the automatic
  setup script.

**Other Changes:**

* Added support for the Grid Engine job queue manager.  See :ref:`hpc-usage-label`.
* Added a configurable parameter, ``minConsFreq``, to the create_snp_matrix.py script.  This parameter specifies
  the mimimum fraction of reads that must agree at a position to make a consensus call.  Prior to version
  0.3.2, the snp pipeline required that a majority (more than half) of the reads must agree to make
  a snp call.  In version 0.3.2, the default behavior requires at least 60% of reads must
  agree to make a consensus call.
* Changed the included snp matrix files for the agona and listeria data sets to match the new results
  obtained by setting minConsFreq=0.6.  The lambda virus results were not impacted by this change.
* Revised the Installation instructions with more detailed step-by-step procedures.
* Added a Dockerfile for automated docker builds.  This feature is still experimental.


0.3.1 (2014-10-27)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Fixed a Python 2.6 incompatibility with formatting syntax when printing the available RAM.
  Also added the Python version to the log files.


0.3.0 (2014-10-22)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Fixed some Mac OSX incompatibilities.
* Fixed a bug in copy_snppipeline_data.py that caused copy failure when the destination
  directory did not exist.
* Fixed alignSampleToReference.sh to properly handle unpaired gzipped fastq files.

**Installation Changes:**

* There is a new dependency on the python psutil package.  When you install the SNP Pipeline,
  pip will attempt to install the psutil package automatically.  If it fails, you may need to
  manually install the python-dev package.  In Ubuntu, ``sudo apt-get install python-dev``


**Other Changes:**

*Note a possible loss of backward compatibilty for existing workflows using
alignSampleToReference.sh and prepSamples.sh*


* All-in-one script: Added a new script, run_snp_pipeline.sh, to run the entire pipeline either on
  a workstation or on a High Performance Computing cluster with the Torque job
  queue manager.  See :ref:`all-in-one-script-label`.
* Logging: The run_snp_pipeline.sh script adds consistent logging functionality for
  workstation and HPC runs.  The logs for each pipeline run are stored in a
  time-stamped directory under the output directory.  See :ref:`logging-label`.
* Timestamp checking: Changed the python scripts (create_snp_list.py, create_snp_pileup.py, create_snp_matrix.py, create_snp_reference.py)
  to skip processing steps when result files already exist and are newer than the input
  files.  If you modify an upstream file, any dependent downstream files will be rebuilt.
  You can force processing regardless of file timestamps with the ``-f`` option.
  Similar functionality for the shell scripts was previously implemented in release 0.2.0.
* Mirrored input files: The run_snp_pipeline.sh script has the capability to make a mirrored copy
  of the input reference and samples to avoid polluting a clean repository.  You have the
  choice to create copies, soft links, or hard links.  See :ref:`mirrored-input-label`.
* Configuration file: Added the capability to customize the behavior of the SNP Pipeline by specifying parameters
  either in a configuration file, or in environment variables.  You can create a configuration
  file with default values pre-set by executing ``copy_snppipeline_data.py configurationFile``
  from the command line.  Pass the configuration file to the run_snp_pipeline.sh script with
  the ``-c`` option.  Alternatively, environment variables matching the names of the
  parameters in the configuration file can be manually set (be sure to export the variables).
  When the run_snp_pipeline.sh script is run, it copies the configuration file for the run into
  the log directory for the run. See :ref:`configuration-label`.
* Removed the ``-p INT`` command line option, to specify the number of cpu cores, from the
  alignSampleToReference.sh script.  You can now control the number of cpu cores used by bowtie2
  with the ``-p INT`` option either in the configuration file when running run_snp_pipeline.sh, or
  in the ``Bowtie2Align_ExtraParams`` environment variable when running alignSampleToReference.sh
  directly. If not specified, it defaults to 8 cpu cores on a HPC cluster, or all cpu cores on
  a workstation.
* Removed the ``--min-var-freq 0.90`` varscan mpileup2snp option from the prepSamples.sh script.
  This parameter is now specified in the ``VarscanMpileup2snp_ExtraParams`` environment variable
  or in the configuration file.
* Listeria monocytogenes data set: Added a Listeria monocytogenes data set.  Updated the usage instructions, illustrating
  how to download the Listeria samples from NCBI and how to run the SNP Pipeline on the
  Listeria data set.  The distribution includes the expected result files for the Listeria
  data set.  Note that due to the large file sizes, the Listeria expected results data set
  does not contain all the intermediate output files.
* Added a command reference page to the documentation.  See :ref:`cmd-ref-label`.


0.2.1 (2014-09-24)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* Version 0.2.0 was missing the Agona data files in the Python distribution.  The
  GitHub repo was fine.  The missing files only impacted PyPi.  Add the Agona
  data files to the Python distribution file list.


0.2.0 (2014-09-17)
~~~~~~~~~~~~~~~~~~

**Changes Impacting Results:**

* Previously, the pipeline executed SAMtools mpileup twice -- the first pileup across
  the whole genome, and the second pileup restricted to those positions where snps
  were identified by varscan in *any* of the samples.  This release removes the
  second SAMtools pileup, and generates the snp pileup file by simply extracting a
  subset of the pileup records from the genome-wide pileup at the positions where
  variants were found in *any* sample.  The consequence of this change is faster run
  times, but also an improvement to the results -- there will be fewer missing
  values in the snp matrix.
* Changed the the supplied lambda virus expected results data set to match the
  results obtained with the pipeline enhancements in this release and now using SAMtools
  version 0.1.19.  SAMtools mpileup version 0.1.19 excludes read bases with low quality.
  As a reminder, the expected results files are fetched with the copy_snppipeline_data.py
  script.
* Removed the "<unknown description>" from the snp matrix fasta file.

**Other Changes:**

*Note the loss of backward compatibilty for existing workflows using prepReference.sh,
alignSampleToReference.sh, prepSamples.sh, create_snp_matrix.py*

* Split the create_snp_matrix script into 4 smaller scripts to simplify the code
  and improve performance when processing many samples in parallel.  Refer to the
  :ref:`usage-label` section for the revised step-by-step usage instructions. The
  rewritten python scripts emit their version number, arguments, run timestamps,
  and other diagnostic information to stdout.
* Changed the default name of the reads.pileup file to reads.snp.pileup.  You can
  override this on the command line of the create_snp_pileup.py script.
* Added the referenceSNP.fasta file to the supplied lambda virus expected results
  data set.
* Updated the usage instructions, illustrating how to download the Agona samples from
  NCBI and how to run the SNP Pipeline on the Agona data set.
* Updated the supplied expected result files for the Agona data set.  Note that due to
  the large file sizes, the Agona expected results data set does not contain all
  the intermediate output files.
* Improved the online help (usage) for all scripts.
* The copy_snppipeline_data.py script handles existing destination directories more
  sensibly now.  The example data is copied into the destination directory if the directory
  already exists.  Otherwise the destination directory is created and the example data
  files are copied there.
* Changed the alignSampleToReference.sh script to specify the number of CPU cores with
  the -p flag, rather than a positional argument.  By default, all CPU cores are
  utilized during the alignment.
* Changed the shell scripts (prepReference.sh, alignSampleToReference.sh, prepSamples.sh)
  to expect the full file name of the reference including the fasta extension, if any.
* Changed the shell scripts (prepReference.sh, alignSampleToReference.sh, prepSamples.sh)
  to skip processing steps when result files already exist and are newer than the input
  files.  If you modify an upstream file, any dependent downstream files will be rebuilt.
  You can force processing regardless of file timestamps with the ``-f`` option.
* Changed the name of the sorted bam file to reads.sorted.bam.
* Changed the general-case usage instructions to handle a variety of fastq file
  extensions (\*.fastq\* and \*.fq\*).


0.1.1 (2014-07-28)
~~~~~~~~~~~~~~~~~~

**Bug fixes:**

* The snp list, snp matrix, and referenceSNP files were incorrectly sorted by
  position alphabetically, not numerically.
* The SNP Pipeline produced slightly different pileups each time we ran the pipeline.
  Often we noticed two adjacent read-bases swapped in the pileup files.  This was
  caused by utilizing multiple CPU cores during the bowtie alignment.  The output
  records in the SAM file were written in non-deterministic order when bowtie ran
  with multiple concurrent threads.  Fixed by adding the ``--reorder`` option to the
  bowtie alignment command line.
* The snp list was written to the wrong file path when the main working directory
  was not specified with a trailing slash.

**Other Changes:**

*Note the loss of backward compatibilty for existing workflows using prepSamples.sh*

* Moved the bowtie alignment to a new script, alignSampleToReference.sh, for
  better control of CPU core utilization when running in HPC environment.
* Changed the prepSamples.sh calling convention to take the sample directory,
  not the sample files.
* prepSamples.sh uses the CLASSPATH environment variable to locate VarScan.jar.
* Changed prepReference.sh to run ``samtools faidx`` on the reference.  This
  prevents errors later when multiple samtools mpileup processes run concurrently.
  When the faidx file does not already exist, multiple samtools mpileup processes
  could interfere with each other by attempting to create it at the same time.
* Added the intermediate lambda virus result files (\*.sam, \*.pileup, \*.vcf) to the
  distribution to help test the installation and functionality.
* Changed the usage instructions to make use of all CPU cores.
* Log the executed commands (bowtie, samtools, varscan) with all options to stdout.

0.1.0 (2014-07-03)
~~~~~~~~~~~~~~~~~~

* Basic functionality implemented.
* Lambda virus tests created and pass.
* S. Agona tests created -- UNDER DEVELOPMENT
* Installs properly from PyPI.
* Documentation available at ReadTheDocs.
