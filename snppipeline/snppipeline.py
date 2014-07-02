#!/usr/bin/env python2.7

from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing
import os
import pprint
import utils


#TODO use os.path.join consistently through code for path creation

def run_snp_pipeline(options_dict):
    """Create SNP matrix

    Description:
    Create a SNP matrix. This function expects, or creates '(*)', the following
        files arranged in the following way:
            mainPath
                reference.fasta
                sampleDirectoryNames.txt
                snplist.txt (*)
                snpma.fasta (*)
            samples
                sample_name_one/reads.bam
                sample_name_one/reads.pileup (*)
                sample_name_one/var.flt.vcf
                ...

    The files are used as follows:
        1. The reference.fasta file is used for alignment of the sequence
            data in the running of the samtools pileup.
        2. The sampleDirectoryNames.txt file contains a list of the paths to the sample
            directories.
        3. The snplist.txt file contains the list of SNPs extracted from the
            var.flt.vcf file.
        4. The snpma.fasta file contains the SNP calls for each sequence,
            arranged as a fasta file with one sequence per sample.
        5. The reads.bam file is used by the samtools pileup program to
            generate the pileup files for each sequence.
        6. The reads.pileup files are used to determine the nucleotide base at
            each SNP position for each sample to construct the SNP fasta file.
        7. The variant file var.flt.vcf is used to construct the SNP position
            list.

    The samtool pileups are run in parallel using the python multiprocessing
        package.

    The vcf and bam files are created outside of this function. The package
        documentation provides an example of creating these files based on the
        lambda_virus sequence that is used as one test for this package.

    Args:
        maxThread: (15) Max number of concurrent threads.
        mainPath:  (no default) Directory containing all input files. Output
            files will also be written here.
        Reference: (no default) File name for reference sequence (in fasta
            format) for mapping.
        pathFileName: ("sampleDirectoryNames.txt") Name of file containing full paths to
            directories containing information for each sequence.
        snplistFileName: Text format list of SNP positions in samples.
        snpmaFileName: File name for snp matrix, formatted as a fasta file,
            with each sequence (all of identical length) corresponding to the
            SNPs in the correspondingly named sequence.
        bamFileName: #TODO
        pileupFileName: Name for pileup files. One is generated for each
            sample, and placed in the corresponding directory for each sample.

    Raises:

    Examples:
    args_dict = {'maxThread':2,
                 'mainPath':'',
                 'Reference':'',
                 'pathFileName':'sampleDirectoryNames.txt',
                 'snplistFileName':'snplist.txt',
                 'snpmaFileName':'snpma.fa',
                 'bamFileName':'reads.bam',
                 'pileupFileName':'reads.pileup',
                }
    run_snp_pipeline(options_dict)
    """

    #==========================================================================
    #Prep work
    #Note use of filter on list_of_sample_directories to remove blank lines.
    #==========================================================================

    verbose        = options_dict['verbose'] > 0
    verbose_print  = print         if verbose else lambda *a, **k: None
    verbose_pprint = pprint.pprint if verbose else lambda *a, **k: None

    verbose_print("Options:")
    verbose_print(options_dict)

    sample_directories_list_filename = os.path.join(options_dict['mainPath'], options_dict['pathFileName'])
    list_of_sample_directories = [line.rstrip() for line in open(sample_directories_list_filename, "r")]
    list_of_sample_directories = filter(None, list_of_sample_directories)

    #==========================================================================
    #read in all vcf files and process into dict of SNPs passing various
    #  criteria. Do this for each sample. Write to file.
    #==========================================================================

    snp_dict = utils.convert_vcf_files_to_snp_dict(list_of_sample_directories, options_dict)
    snp_list_file_path = options_dict['mainPath'] + options_dict['snplistFileName'] # TODO Fix this
    utils.write_list_of_snps(snp_list_file_path, snp_dict)

    #==========================================================================
    # Generate Pileups of samples (in parallel)
    # Note the rather peculiar use of the pool. I.e., we use:
    #   result_many = pool.map_async(utils.pileup_wrapper,
    #                                parameter_list).get(9999999)
    #  And not the maybe more sensible:
    #   result_many = pool.map(utils.pileup_wrapper, parameter_list)
    #  This is to ensure that a keyboard interupt (Ctrl-C) is properly handled.
    #  The approach used is from:
    #   http://stackoverflow.com/questions/1408356/keyboard-interrupts-with
    #   -pythons-multiprocessing-pool
    #  There are alternate approaches, one of which is in:
    #   http://noswap.com/blog/python-multiprocessing-keyboardinterrupt
    #==========================================================================

    #create a list of tuples containing values need for pileup code (as passed
    #  via pileup code wrapper)
    parameter_list = zip(list_of_sample_directories,
                         len(list_of_sample_directories)*[options_dict])

    verbose_print("Starting Pileups.")

    pool        = multiprocessing.Pool(processes=options_dict['maxThread'])
    result_many = pool.map_async(utils.pileup_wrapper,
                                 parameter_list).get(9999999)

    #verbose_pprint(result_many)
    verbose_print("Pileups are finished.")

    #==========================================================================
    #   Create snp matrix. Write results to file.
    #==========================================================================

    snp_sequence_records_list = []

    for sample_directory in list_of_sample_directories:

        sample_name                  = sample_directory.split(os.sep)[-1]
        pileup_file_path             = os.path.join(sample_directory, options_dict['pileupFileName'])
        position_consensus_base_dict = utils.create_consensus_dict(pileup_file_path)

        snp_seq_string = ""
        #TODO - is there more than one chromosome?  If so, why is there only one snp_seq_string?
        for key in sorted(snp_dict.iterkeys()):  #TODO - Why sorted?
            chrom, pos = key.split()
            # TODO how large is this dictionary?  Do we need a more efficient storage?
            if position_consensus_base_dict.has_key(chrom + ":" + pos):
                snp_seq_string += position_consensus_base_dict[chrom + ":" + pos]
            else:
                snp_seq_string += "-"

        snp_seq_record = SeqRecord(Seq(snp_seq_string), id=sample_name)
        snp_sequence_records_list.append(snp_seq_record)

    #Write bases for snps for each sequence to a fasta file
    snpma_file_path = os.path.join(options_dict['mainPath'], options_dict['snpmaFileName'])
    with open(snpma_file_path, "w") as fasta_file_object:
        SeqIO.write(snp_sequence_records_list, fasta_file_object, "fasta")

    #==========================================================================
    #    Write reference sequence bases at SNP locations to a fasta file.
    #==========================================================================
    if options_dict['includeReference']:
        snp_list_file_path       = os.path.join(options_dict['mainPath'], options_dict['snplistFileName'])
        reference_file_path      = os.path.join(options_dict['mainPath'], options_dict['Reference'])
        snp_reference_file_path  = os.path.join(options_dict['mainPath'], "referenceSNP.fasta")
        utils.write_reference_snp_file(reference_file_path, snp_list_file_path,
                                       snp_reference_file_path)


