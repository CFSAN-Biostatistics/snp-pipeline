#!/usr/bin/env python2.7

from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pprint
import utils
import sys
import time
import platform
import psutil
import locale
import pileup as pileup
import vcf_writer as vcf_writer
from __init__ import __version__


verbose_print  = lambda *a, **k: None
verbose_pprint = lambda *a, **k: None

def set_logging_verbosity(options_dict):
    """Enable or disable logging.

    Args:
        verbose : Verbosity value, any value greater than 0 enables logging
    """
    global verbose_print
    global verbose_pprint
    verbose_print  = print         if options_dict['verbose'] > 0 else lambda *a, **k: None
    verbose_pprint = pprint.pprint if options_dict['verbose'] > 0 else lambda *a, **k: None


def timestamp():
    """Return a timestamp string."""
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

def program_name():
    """Return the basename of the python script being executed."""
    return os.path.basename(sys.argv[0])

def command_line_short():
    """Return the command line string without the full path to the program."""
    return "%s %s" % (program_name(), " ".join(sys.argv[1:]))

def command_line_long():
    """Return the command line string with the full path to the program."""
    return " ".join(sys.argv)

def print_log_header():
    """Print a standardized header for the log with starting conditions."""
    verbose_print("# Command           : %s" % command_line_long())
    verbose_print("# Working Directory : %s" % os.getcwd())
    pbs_jobid = os.environ.get("PBS_JOBID")
    sge_jobid = os.environ.get("JOB_ID")
    sge_task_id = os.environ.get("SGE_TASK_ID")
    if sge_task_id == "undefined":
        sge_task_id = None
    if pbs_jobid:
        verbose_print("# Job ID            : %s" % pbs_jobid)
    elif sge_jobid and sge_task_id:
        verbose_print("# Job ID            : %s[%s]" % (sge_jobid, sge_task_id))
    elif sge_jobid:
        verbose_print("# Job ID            : %s" % sge_jobid)

    verbose_print("# Hostname          : %s" % platform.node())
    locale.setlocale(locale.LC_ALL, '')
    ram_mbytes = psutil.virtual_memory().total / 1024 / 1024
    ram_str = locale.format("%d", ram_mbytes, grouping=True)
    verbose_print("# RAM               : %s MB" % ram_str)
    verbose_print("# Python Version    : %s" % sys.version.replace("\n", " "))
    verbose_print("")


def print_arguments(options_dict):
    """Print the program options.

    Inputs:
        options_dict : Dictionary of program arguments
    """
    verbose_print("Options:")
    for key in options_dict.keys():
        verbose_print("    %s=%s" % (key, options_dict[key]))


def create_snp_list(options_dict):
    """Create SNP list file

    Description:
    Create the SNP list -- the list of positions where variants were found
    and the corresponding list of samples having a variant at each position. 
    This function expects, or creates '(*)', the following files arranged 
    in the following way:
            sampleDirectories.txt
            samples
                sample_name_one/var.flt.vcf
                ...
            snplist.txt (*)

    The files are used as follows:
        1. The sampleDirectories.txt input file contains a list of the paths to 
           the sample directories.
        2. The var.flt.vcf variant input files are used to construct the 
           SNP position list.
        3. The snplist.txt output file contains the union of the SNP positions 
           and sample names extracted from all the var.flt.vcf files.

    The sampleDirectories.txt and var.flt.vcf files are created outside of 
    this function. The package documentation provides an example of creating 
    these files based on the lambda_virus sequence that is used as one test 
    for this package.

    Args:
        sampleDirsFile: File path (not just file name) of file containing paths 
            to directories containing var.flt.vcf file for each sequence.
        vcfFileName: File name of the VCF files which must exist in each of the
            sample directories
        snpListFile: File path (not just file name) of text format list 
            of SNP positions

    Raises:

    Examples:
    options_dict = {'sampleDirsFile':'sampleDirectories.txt',
                    'vcfFileName':'var.flt.vcf'
                    'snpListFile':'snplist.txt',
                   }
    create_snp_list(options_dict)
    """
    print_log_header()
    verbose_print("# %s %s" % (timestamp(), command_line_short()))
    verbose_print("# %s version %s" % (program_name(), __version__))
    print_arguments(options_dict)

    #==========================================================================
    # Prep work
    # Note use of filter on list_of_sample_directories to remove blank lines.
    #==========================================================================
    sample_directories_list_filename = options_dict['sampleDirsFile']
    list_of_sample_directories = [line.rstrip() for line in open(sample_directories_list_filename, "r")]
    list_of_sample_directories = sorted(filter(None, list_of_sample_directories))

    #==========================================================================
    # Read in all vcf files and process into dict of SNPs passing various
    # criteria. Do this for each sample. Write to file.
    #==========================================================================
    snp_list_file_path = options_dict['snpListFile']
    vcf_file_name = options_dict['vcfFileName']
    list_of_vcf_files = [os.path.join(dir, vcf_file_name) for dir in list_of_sample_directories]

    if options_dict['forceFlag'] or utils.target_needs_rebuild(list_of_vcf_files, snp_list_file_path):
        snp_dict = utils.convert_vcf_files_to_snp_dict(list_of_vcf_files)
        verbose_print('Found %d snp positions across %d sample vcf files.' % (len(snp_dict), len(list_of_vcf_files)))
        utils.write_list_of_snps(snp_list_file_path, snp_dict)
        verbose_print("")
    else:
        verbose_print("SNP list %s has already been freshly built.  Use the -f option to force a rebuild." % snp_list_file_path)
    verbose_print("# %s %s finished" % (timestamp(), program_name()))


def create_snp_pileup(options_dict):
    """Create the SNP pileup file for a sample.

    Description:
    Create the SNP pileup file for a sample -- the pileup file restricted to 
    only positions where variants were found in any sample.
    This function expects, or creates '(*)', the following files arranged 
    in the following way:
            snplist.txt
            samples
                sample_name_one/reads.all.pileup
                sample_name_one/reads.snp.pileup (*)
                ...

    The files are used as follows:
        1. The snplist.txt input file contains the list of SNP positions 
           extracted from the var.flt.vcf file.
        2. The reads.all.pileup input file is the genome-wide pileup file 
           for this sample.
        3. The reads.snp.pileup output file is the pileup file for this sample,
           restricted to only positions where variants were found in any 
           sample.

    The snplist.txt and reads.all.pileup files are created outside of this 
    function. The package documentation provides an example of creating these
    files based on the lambda_virus sequence that is used as one test for 
    this package.

    Args:
        snpListFile: File path (not just file name) of text format list 
            of SNP positions across all samples
        allPileupFile: File path (not just file name) of the whole-genome
            pileup file fot this sample
        snpPileupFile: File path (not just file name) of the snp pileup file

    Raises:

    Examples:
    options_dict = {'snpListFile':'snplist.txt',
                    'allPileupFile':'samples/SRR555888/reads.all.pileup'
                    'snpPileupFile':'samples/SRR555888/reads.snp.pileup'
                   }
    create_snp_pileup(options_dict)
    """
    print_log_header()
    verbose_print("# %s %s" % (timestamp(), command_line_short()))
    verbose_print("# %s version %s" % (program_name(), __version__))
    print_arguments(options_dict)

    snp_list_file_path = options_dict['snpListFile']
    all_pileup_file_path = options_dict['allPileupFile']
    snp_pileup_file_path = options_dict['snpPileupFile']

    source_files = [snp_list_file_path, all_pileup_file_path]
    if options_dict['forceFlag'] or utils.target_needs_rebuild(source_files, snp_pileup_file_path):
        # Create a pileup file with a subset of the whole-genome pileup restricted
        # to locations with SNPs only.
        snp_list = utils.read_snp_position_list(snp_list_file_path)
        utils.create_snp_pileup(all_pileup_file_path, snp_pileup_file_path, set(snp_list))
        verbose_print("")
    else:
        verbose_print("SNP pileup %s has already been freshly built.  Use the -f option to force a rebuild." % snp_pileup_file_path)
    verbose_print("# %s %s finished" % (timestamp(), program_name()))


def call_consensus(options_dict):
    """Call the consensus base for a sample

    Call the consensus base for a sample at the positions where SNPs were found
    in any of the samples.
    This function expects, or creates '(*)', the following
        files arranged in the following way:
            snplist.txt
            samples
                sample_name_one/reads.all.pileup
                sample_name_one/consensus.fasta (*)

    The files are used as follows:
        1. The snplist.txt input file contains the list of SNP positions 
           extracted from all the var.flt.vcf files combined.
        2. The reads.all.pileup input file is a pileups at all positions
           used to determine the nucleotide base at each SNP position.
        3. The consensus.fasta output file contains the SNP calls for each 
           sequence, arranged as a fasta file with one sequence per sample.

    The snplist.txt, and reads.snp.pileup are created outside of this function.
       The package documentation provides an example 
        of creating these files based on the lambda_virus sequence that is used 
        as one test for this package.

    Args:
        forceFlag : boolean 
            flag to force processing even when result file already exists and
            is newer than inputs
        snpListFile : str
            File path (not just file name) of text format list of SNP positions
        allPileupFile : str
            Relative or absolute path to the genome-wide pileup file for this 
            sample
        consensusFile : str
            Output file. Relative or absolute path to the consensus fasta file
            for this sample.
        minBaseQual : int
            Mimimum base quality score to count a read. All other snp filters
            take effect after the low-quality reads are discarded.
        minConsFreq : float
            Consensus frequency. Mimimum fraction of high-quality reads
            supporting the consensus to make a call.
        minConsStrdDpth : int
            Consensus strand depth. Minimum number of high-quality reads 
            supporting the consensus which must be present on both the
            forward and reverse strands to make a call.
        minConsStrdBias : float
            Strand bias. Minimum fraction of the high-quality 
            consensus-supporting reads which must be present on both the 
            forward and reverse strands to make a call. The numerator of this
            fraction is the number of high-quality consensus-supporting reads
            on one strand.  The denominator is the total number of high-quality
            consensus-supporting reads on both strands combined.

    Raises:

    Examples:
    options_dict = {'snpListFile':'snplist.txt',
                    'allPileupFile':'reads.all.pileup',
                    'consensusFile':'consensus.fasta',
                    'minBaseQual':15,
                    'minConsFreq':0.6,
                    'minConsStrdDpth':4
                    'minConsStrdBias':0.10
                   }
    call_consensus(options_dict)
    """
    print_log_header()
    verbose_print("# %s %s" % (timestamp(), command_line_short()))
    verbose_print("# %s version %s" % (program_name(), __version__))
    print_arguments(options_dict)

    snp_list_file_path = options_dict['snpListFile']
    all_pileup_file_path = options_dict['allPileupFile']
    sample_directory = os.path.dirname(os.path.abspath(all_pileup_file_path))
    sample_name = os.path.basename(sample_directory)
    consensus_file_path = options_dict['consensusFile']
    consensus_file_dir = os.path.dirname(os.path.abspath(consensus_file_path))
    vcf_file_name = options_dict['vcfFileName']
    vcf_file_path = os.path.join(consensus_file_dir, vcf_file_name) if vcf_file_name else None

    # Check if the result is already fresh
    source_files = [snp_list_file_path, all_pileup_file_path]
    if not options_dict['forceFlag'] and not utils.target_needs_rebuild(source_files, consensus_file_path):
        verbose_print("Consensus call file %s has already been freshly built.  Use the -f option to force a rebuild." % consensus_file_path)
        verbose_print("# %s %s finished" % (timestamp(), program_name()))
        return

    # Load the list of which positions to called
    snp_list = utils.read_snp_position_list(snp_list_file_path)
    snplist_length = len(snp_list)
    verbose_print("snp position list length = %d" % snplist_length)

    # Call consensus. Write results to file.
    position_consensus_base_dict = dict()

    caller = pileup.ConsensusCaller(options_dict['minConsFreq'], 
                                    options_dict['minConsStrdDpth'], 
                                    options_dict['minConsStrdBias'])
    snp_positions = set(snp_list)
    parse_positions = None if options_dict['vcfAllPos'] else snp_positions
    pileup_reader = pileup.Reader(all_pileup_file_path, 
                                  options_dict['minBaseQual'], 
                                  parse_positions)
    if vcf_file_name:
        writer = vcf_writer.SingleSampleWriter(vcf_file_path)
        filters = caller.get_filter_descriptions()
        writer.write_header(sample_name, filters, options_dict['vcfRefName'])
    for pileup_record in pileup_reader:
        chrom = pileup_record.chrom
        pos = pileup_record.position
        consensus_base, fail_reasons = caller.call_consensus(pileup_record)
        if (chrom, pos) in snp_positions:
            if fail_reasons:
                position_consensus_base_dict[(chrom, pos)] = '-'
            else:
                position_consensus_base_dict[(chrom, pos)] = consensus_base

        if vcf_file_name:
            writer.write_from_pileup(pileup_record, fail_reasons)
    if vcf_file_name:
        writer.close()

    verbose_print("called consensus positions = %i" % (len(position_consensus_base_dict)))

    consensus_list = [position_consensus_base_dict.get(key, '-') for key in snp_list]
    consensus_str = ''.join(consensus_list)
    snp_seq_record = SeqRecord(Seq(consensus_str), id=sample_name, description="")

    # Write the consensus calls to a fasta file
    with open(consensus_file_path, "w") as fasta_file_object:
        SeqIO.write([snp_seq_record], fasta_file_object, "fasta")

    verbose_print("")
    verbose_print("# %s %s finished" % (timestamp(), program_name()))


def create_snp_matrix(options_dict):
    """Create SNP matrix

    Description:
    Create the SNP matrix containing the consensus base for each of the samples 
    at the positions where SNPs were found in any of the samples.  The matrix 
    contains one row per sample and one column per SNP position.  Non-SNP 
    positions are not included in the matrix.
    This function expects, or creates '(*)', the following
        files arranged in the following way:
            sampleDirectories.txt
            samples
                sample_name_one/consensus.fasta
                ...
            snpma.fasta (*)

    The files are used as follows:
        1. The sampleDirectories.txt input file contains a list of the paths to 
           the sample directories.
        2. The consensus.fasta input files are previously called consensus
           for each sample to construct the SNP matrix fasta file.
        3. The snpma.fasta output file contains the SNP calls for each 
           sequence, arranged as a multi-fasta file with one sequence per 
           sample.

    The sampleDirectories.txt, and consensus.fasta are created outside of this
        function. The package documentation provides an example of creating 
        these files based on the lambda_virus sequence that is used as one 
        test for this package.

    Args:
        sampleDirsFile : str
            File path (not just file name) of file containing paths 
            to directories containing consensus.fasta file for each sequence.
        snpListFile : str
            File path (not just file name) of text format list of SNP positions
        consFileName : str
            File name of the previously called consensus fasta files which must
            exist in each of the sample directories
        snpmaFile : str
            File path (not just file name) of the output snp matrix, formatted
            as a fasta file, with each sequence (all of identical length) 
            corresponding to the SNPs in the correspondingly named sequence.

    Raises:

    Examples:
    options_dict = {'sampleDirsFile':'sampleDirectories.txt',
                    'consFileName':'consensus.fasta',
                    'snpmaFile':'snpma.fasta',
                    'minConsFreq':0.6,
                   }
    create_snp_matrix(options_dict)
    """
    print_log_header()
    verbose_print("# %s %s" % (timestamp(), command_line_short()))
    verbose_print("# %s version %s" % (program_name(), __version__))
    print_arguments(options_dict)

    #==========================================================================
    # Prep work
    # Note use of filter on list_of_sample_directories to remove blank lines.
    #==========================================================================
    sample_directories_list_filename = options_dict['sampleDirsFile']
    list_of_sample_directories = [line.rstrip() for line in open(sample_directories_list_filename, "r")]
    list_of_sample_directories = sorted(filter(None, list_of_sample_directories))

    #==========================================================================
    # Check if the result is already fresh
    #==========================================================================
    snpma_file_path = options_dict['snpmaFile']
    source_files = []
    if not options_dict['forceFlag']:
        for sample_directory in list_of_sample_directories:
            consensus_file_path = os.path.join(sample_directory, options_dict['consFileName'])
            source_files.append(consensus_file_path)
        if not utils.target_needs_rebuild(source_files, snpma_file_path):
            verbose_print("SNP matrix %s has already been freshly built.  Use the -f option to force a rebuild." % snpma_file_path)
            verbose_print("# %s %s finished" % (timestamp(), program_name()))
            return

    #==========================================================================
    #   Create snp matrix. Write results to file.
    #==========================================================================
    with open(snpma_file_path, "w") as output_file:
        for sample_directory in list_of_sample_directories:
            consensus_file_path = os.path.join(sample_directory, options_dict['consFileName'])
            verbose_print("Merging " + consensus_file_path)
            with open(consensus_file_path, "r") as input_file:
                for line in input_file:
                    output_file.write(line)

    verbose_print("")
    verbose_print("# %s %s finished" % (timestamp(), program_name()))



def create_snp_reference_seq(options_dict):
    """Write reference sequence bases at SNP locations to a fasta file.

    Description:
    Write reference sequence bases at SNP locations to a fasta file. 
    This function expects, or creates '(*)', the following files:
            reference.fasta
            snplist.txt
            referenceSNP.fasta (*)

    The files are used as follows:
        1. The reference.fasta input file contains the whole-genome reference 
           bases.
        2. The snplist.txt input file contains the list of SNP positions across
           all the samples.
        2. The referenceSNP.fasta output file contains the reference bases at 
           the identified SNP locations.

    The snplist.txt file is created outside of this function.  The package 
        documentation provides an example of creating this file based on the
        lambda_virus sequence that is used as one test for this package.

    Args:
        referenceFile: File path (not just file name) for reference sequence 
            (in fasta format
        snpListFile: File path (not just file name) of text format list of SNP
            positions
        snpRefFile: File path (not just file name) for the SNP reference 
            sequence file.

    Raises:

    Examples:
    options_dict = {'referenceFile':'reference.fasta',
                    'snpListFile':'snplist.txt',
                    'snpRefFile':'referenceSNP.fasta'
                   }
    create_snp_reference_seq(options_dict)
    """
    print_log_header()
    verbose_print("# %s %s" % (timestamp(), command_line_short()))
    verbose_print("# %s version %s" % (program_name(), __version__))
    print_arguments(options_dict)

    #==========================================================================
    #    Write reference sequence bases at SNP locations to a fasta file.
    #==========================================================================
    reference_file = options_dict['referenceFile']
    snp_list_file_path = options_dict['snpListFile']
    snp_ref_seq_path = options_dict['snpRefFile']

    source_files = [reference_file, snp_list_file_path]
    if options_dict['forceFlag'] or utils.target_needs_rebuild(source_files, snp_ref_seq_path):
        utils.write_reference_snp_file(reference_file, snp_list_file_path, snp_ref_seq_path)
        verbose_print("")
    else:
        verbose_print("SNP reference sequence %s has already been freshly built.  Use the -f option to force a rebuild." % snp_ref_seq_path)

    verbose_print("# %s %s finished" % (timestamp(), program_name()))


