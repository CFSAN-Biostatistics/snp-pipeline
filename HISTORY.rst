.. :changelog:

History
-------

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
