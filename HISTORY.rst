.. :changelog:

History
-------

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
