.. _configuration-label:

=============
Configuration
=============

You can customize the behavior of the SNP Pipeline by configuring parameters.
Each step in the pipeline has a corresponding parameter allowing you to set one
or more options for each tool in the pipeline.

Parameters can be configured either in a configuration file if you are using the
``run_snp_pipeline.sh`` script, or in environment variables.

The pipeline comes with a default configuration file.  When you run the run_snp_pipeline.sh
script without specifying a configuration file, it automatically uses the
default supplied configuration file.

To get a copy of the default configuration file, run the following command.  This
will create a file called ``snppipeline.conf``::

    cfsan_snp_pipeline data configurationFile

To customize the pipeline behavior, edit the configuration file and pass the file to
the run_snp_pipeline.sh script::

    run_snp_pipeline.sh -c snppipeline.conf ...

When the run_snp_pipeline.sh script runs, it copies the configuration file to the
log directory for the run, capturing the configuration used for each run.

If you decide not to use the run_snp_pipeline.sh script, you can still customize the
behavior of the pipeline tools, but you will need to set (and export) environment
variables with the same names as the parameters in the configuration file.

The available configuration parameters are described below.

StopOnSampleError
-----------------
Controls whether the pipeline exits upon detecting errors affecting only a single
sample.  The pipeline will always stop upon detecting global errors affecting all
samples.

**Default**:

    When this parameter is not set to a value, the pipeline will stop upon detecting
    single sample errors.  If you want the pipeline to continue, you must explicitly set
    this parameter false.

**Example**::

    StopOnSampleError=false


MaxCpuCores
-----------
Controls the number of CPU cores used by the SNP Pipeline.  This parameter controls
CPU resources on both your workstation and on an HPC cluster.  By limiting the
number of CPU cores used, you are also limiting the number of concurrently executing
processes.

**Default**:

    When this parameter is not set to a value, the pipeline will launch multiple concurrent
    processes using all available CPU cores.

**Example**::

    MaxCpuCores=2


MaxSnps
-------
Controls the maximum number of snps allowed for each sample. Any sample with excessive snps exceeding
this limit will be excluded from the snp list, snp matrix, and snpma.vcf file. When set to -1, this
parameter is disabled.

**Default**:

    Do not leave this parameter unset.  To disable the excessive snp filtering and include all samples
    regardless of the number of snps, set the parameter to -1

**Example**::

    MaxSnps=1000


SnpPipeline_Aligner
-------------------
Controls which reference-based aligner is used to map reads to the reference genome.
The choices are ``bowtie2`` or ``smalt``.

**Default**:

    When this parameter is not set to a value, the pipeline will use the bowtie2 aligner.

**Example**::

    SnpPipeline_Aligner="smalt"


Bowtie2Build_ExtraParams
------------------------

Specifies options passed to the bowtie2 indexer.  Any of the bowtie2-build options
can be specified.

**Default**: none

**Example**::

    Bowtie2Build_ExtraParams="--offrate 3"


SmaltIndex_ExtraParams
------------------------

Specifies options passed to the smalt indexer.  Any of the smalt index options
can be specified.

**Default**: none

**Example**::

    SmaltIndex_ExtraParams="-k 20 -s 1"


SamtoolsFaidx_ExtraParams
-------------------------

Specifies options passed to the SAMtools faidx indexer.  Any of the SAMtools faidx options
can be specified.

**Default**: none

**Example**::

    SamtoolsFaidx_ExtraParams=""


Bowtie2Align_ExtraParams
------------------------

Specifies options passed to the bowtie2 aligner.  Any of the bowtie2 aligner options
can be specified.

**Default**:

|   If you do not specify the ``-p`` option, it defaults to 8 threads on an HPC or all cpu cores otherwise.
|      There is no way to completely suppress the -p option.
|   If Bowtie2Align_ExtraParams is not set to any value, the ``--reorder`` option is enabled by default.
|      Any value, even a single space, will suppress this default option.
|

**Parameter Notes**:

| ``-p``        : bowtie2 uses the specified number of parallel search threads
| ``--reorder`` : generate output records in the same order as the reads in the input file
| ``-X``        : maximum inter-mate fragment length for valid concordant paired-end alignments
|

**Example**::

    Bowtie2Align_ExtraParams="--reorder -p 16 -X 1000"


SmaltAlign_ExtraParams
----------------------

Specifies options passed to the smalt mapper.  Any of the smalt map options
can be specified.

**Default**:

|   If you do not specify the ``-n`` option, it defaults to 8 threads on an HPC or all cpu cores otherwise.
|      There is no way to completely suppress the -n option.
|   If SmaltAlign_ExtraParams is not set to any value, the ``-O`` option is enabled by default.
|      Any value, even a single space, will suppress this default option.
|

**Parameter Notes**:

| ``-n`` : number of parallel alignment threads
| ``-O`` : generate output records in the same order as the reads in the input file
| ``-i`` : maximum insert size for paired-end reads
| ``-r`` : random number seed, if seed < 0 reads with multiple best mappings are reported as 'not mapped'
| ``-y`` : filters output alignments by a threshold in the number of exactly matching nucleotides
|

**Example**::

    SmaltAlign_ExtraParams="-O -i 1000 -r 1"


.. _SamtoolsSamFilter-ExtraParams-label:

SamtoolsSamFilter_ExtraParams
-----------------------------
Specifies options passed to the SAMtools view tool when filtering the SAM file.
Any of the SAMtools view options can be specified.

**Default**:

| If SamtoolsSamFilter_ExtraParams is not set, the "-F 4" option is enabled by default.
|    Any value, even a single space, will suppress the -F option.
|

**Parameter Notes**:

``-F 4``
  Exclude unmapped reads.
``-q threshold``
  Exclude reads with map quality below threshold.

**Example**::

    SamtoolsSamFilter_ExtraParams="-F 4 -q 30"


SamtoolsSort_ExtraParams
------------------------
Specifies options passed to the SAMtools sort tool when sorting the BAM file.
Any of the SAMtools sort options can be specified.

**Default**: None

**Example**::

    SamtoolsSort_ExtraParams=""


RemoveDuplicateReads
--------------------
Controls whether the pipeline removes duplicate reads prior to creating the pileup
and calling snps.

**Default**:

    When this parameter is not set to a value, the pipeline removes duplicate reads.

**Example**::

    RemoveDuplicateReads=false


PicardMarkDuplicates_ExtraParams
--------------------------------
Specifies options passed to the Picard MarkDuplicates tool when removing duplicate reads.

**Default**: None

**Example**::

    PicardMarkDuplicates_ExtraParams="DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH"


PicardJvm_ExtraParams
---------------------
Specifies options passed to the Picard Java Virtual Machine.
Any of the JVM options can be specified.

**Default**: None

**Parameter Notes**:

| ``-Xmx300m``  : use 300 MB memory (modify as needed)
|

**Example**::

    PicardJvm_ExtraParams="-Xmx300m"


SamtoolsMpileup_ExtraParams
---------------------------
Specifies options passed to the SAMtools mpileup tool.
Any of the SAMtools mpileup options can be specified.

**Default**: None

**Parameter Notes**:

| ``-q``    : minimum mapping quality for an alignment to be used
| ``-Q``    : minimum base quality for a base to be considered
| ``-x``    : disable read-pair overlap detection
|

**Example**::

    SamtoolsMpileup_ExtraParams="-q 0 -Q 13"


VarscanMpileup2snp_ExtraParams
------------------------------
Specifies options passed to the Varscan mpileup2snp tool.
Any of the Varscan mpileup2snp options can be specified.

**Default**: None

**Parameter Notes**:

| ``--min-avg-qual`` : minimum base quality at a position to count a read
| ``--min-var-freq`` : minimum variant allele frequency threshold
|

**Example**::

    VarscanMpileup2snp_ExtraParams="--min-avg-qual 15 --min-var-freq 0.90"


VarscanJvm_ExtraParams
----------------------
Specifies options passed to the Varscan Java Virtual Machine.
Any of the JVM options can be specified.

**Default**: None

**Parameter Notes**:

| ``-Xmx300m``  : use 300 MB memory (modify as needed)
|

**Example**::

    VarscanJvm_ExtraParams="-Xmx300m"


.. _FilterRegions-ExtraParams-label:

FilterRegions_ExtraParams
------------------------------
Specifies options passed to the filter_regions command.

**Default**: None

**Parameter Notes**:

``--edge_length``
  The length of the edge regions in a contig, in which all SNPs will be removed.
``--window_size``
  The length of the window in which the number of SNPs should be no more than max_num_snp.
``--max_snp``
  The maximum number of SNPs allowed in a window.
``--out_group``
    Relative or absolute path to the file indicating outgroup samples, one sample ID per line.

You can filter snps more than once by specifying multiple window sizes and max snps.  For example "-m 3 2 -w 1000 100" will filter more than 3 snps in 1000 bases and also more than 2 snps in 100 bases.

**Example**::

    FilterRegions_ExtraParams="--edge_length 500 --window_size 1000 125 15 --max_snp 3 2 1 --out_group /path/to/outgroupSamples.txt"


MergeSites_ExtraParams
-------------------------
Specifies options passed to the merge_sites command.

**Default**: None

**Example**::

    MergeSites_ExtraParams="--verbose 1"


CallConsensus_ExtraParams
-------------------------
Specifies options passed to the call_consensus command.

**Default**: None

**Parameter Notes**:

``--minBaseQual``
    Mimimum base quality score to count a read. All other snp filters take effect after the low-quality reads
    are discarded.
``--minConsFreq``
    Consensus frequency. Mimimum fraction of high-quality reads supporting the consensus to make a call.
``--minConsStrdDpth``
    Consensus strand depth. Minimum number of high-quality reads supporting the consensus which must be present
    on both the forward and reverse strands to make a call
``--minConsStrdBias``
    Strand bias. Minimum fraction of the high-quality consensus-supporting reads which must be present on both
    the forward and reverse strands to make a call. The numerator of this fraction is the number of high-quality
    consensus-supporting reads on one strand. The denominator is the total number of high-quality
    consensus-supporting reads on both strands combined.
``--vcfFileName``
    VCF Output file name. If specified, a VCF file with this file name will be created in the same directory
    as the consensus fasta file for this sample.
``--vcfAllPos``
    Flag to cause VCF file generation at all positions, not just the snp positions. This has no effect on the
    consensus fasta file, it only affects the VCF file. This capability is intended primarily as a diagnostic
    tool and enabling this flag will greatly increase execution time.
``--vcfPreserveRefCase``
    Flag to cause the VCF file generator to emit each reference base in uppercase/lowercase as it appears in the
    reference sequence file.  If not specified, the reference bases are emitted in uppercase.

**Example**::

    CallConsensus_ExtraParams="--verbose 1 --minBaseQual 15 --vcfFileName consensus.vcf"


SnpMatrix_ExtraParams
---------------------------
Specifies options passed to the snp_matrix command.

**Default**: None

**Example**::

    SnpMatrix_ExtraParams="--verbose 1"


SnpReference_ExtraParams
---------------------------------
Specifies options passed to the snp_reference command.

**Default**: None

**Example**::

    SnpReference_ExtraParams="--verbose 1"


MergeVcfs_ExtraParams
---------------------
Specifies options passed to the merge_vcfs command.

**Default**: none

**Example**::

    MergeVcfs_ExtraParams="-n sample.vcf"


BcftoolsMerge_ExtraParams
-------------------------
Specifies options passed to the bcftools merge tool.

**Default**:

    When this parameter is not set to a value, the pipeline uses the settings:
    ``--merge all --info-rules NS:sum``.  Any value, even a single space, will
    suppress the default settings.

**Parameter Notes**:

``--merge``
    Controls the creation of multiallelic records.
        - none   = no new multiallelics, output multiple records instead
        - snps   = allow multiallelic SNP records
        - indels = allow multiallelic indel records
        - both   = both SNP and indel records can be multiallelic
        - all    = SNP records can be merged with indel records
        - id     = merge by ID
``--filter-logic``
    Controls the content of the filter data element.
        - x = set the output record filter to PASS if any of the inputs pass
        - \+ = set the output record filter to PASS when all of the inputs pass
``--info-rules``
    Rules for merging INFO fields (scalars or vectors) or - to disable the default rules. METHOD is one of
    sum, avg, min, max, join. Default is DP:sum,DP4:sum if these fields exist in the input files. Fields
    with no specified rule will take the value from the first input file.

**Example**::

    BcftoolsMerge_ExtraParams="--merge all --info-rules NS:sum"


CollectSampleMetrics_ExtraParams
--------------------------------
Specifies options passed to the collect_metrics command.

**Default**: none

**Example**::

    CollectSampleMetrics_ExtraParams="-v consensus.vcf"


CombineSampleMetrics_ExtraParams
--------------------------------
Specifies options passed to the combine_metrics command.

**Default**: none

**Parameter Notes**:

| ``-s``  : Emit column headings with spaces instead of underscores
|

**Example**::

    CombineSampleMetrics_ExtraParams="-s"


Torque_StripJobArraySuffix
--------------------------
Controls stripping the suffix from the job id when specifying Torque job array dependencies.
It may be necessary to change this parameter if run_snp_pipeline.sh fails with an illegal qsub
dependency error.

**Example**::

    Torque_StripJobArraySuffix=false


GridEngine_StripJobArraySuffix
------------------------------
Controls stripping the suffix from the job id when specifying Grid Engine job array dependencies.
It may be necessary to change this parameter if run_snp_pipeline.sh fails with an illegal qsub
dependency error.

**Example**::

    GridEngine_StripJobArraySuffix=true


GridEngine_PEname
-----------------
Specifies the name of the Grid Engine parallel environment.  This is only needed when running
the SNP Pipeline on a High Performance Computing cluster with the Grid Engine job manager.
Contact your HPC system administrator to determine the name of your parallel environment.
Note: the name of this parameter was PEname in releases prior to 0.4.0.

**Example**::

    GridEngine_PEname="mpi"


GridEngine_QsubExtraParams
--------------------------
Specifies extra options passed to qsub when running the SNP Pipeline on the Grid Engine job scheduler.

**Default**: None

**Example**::

    GridEngine_QsubExtraParams="-q bigmem.q -l h_rt=12:00:00"


Torque_QsubExtraParams
--------------------------
Specifies extra options passed to qsub when running the SNP Pipeline on the Torque job scheduler.

**Default**: None

**Example**::

    Torque_QsubExtraParams="-l pmem=16gb -l walltime=12:00:00"
