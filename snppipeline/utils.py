#!/usr/bin/env python2.7

from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import operator
import os
import pileup
import pprint
import re
import sys
import vcf


#==============================================================================
#Prep work
#==============================================================================

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



def target_needs_rebuild(source_files, target_file):
    """Determine if a target file needs a fresh rebuild, i.e. the target does
    not exist or its modification time is older than any of its source files.

    Args:
        source_files : relative or absolute path to a list of files
        target_file : relative or absolute path to target file
    """
    if not os.path.isfile(target_file):
        return True;

    if os.path.getsize(target_file) == 0:
        return True;

    target_timestamp = os.stat(target_file).st_mtime

    for source_file in source_files:
        source_timestamp = os.stat(source_file).st_mtime
        if source_timestamp > target_timestamp:
            return True

    return False



#==============================================================================
#Define functions
#==============================================================================

def create_snp_pileup(all_pileup_file_path, snp_pileup_file_path, snp_set):
    """Create a subset pileup with SNP locations only.

    Given a whole-genome pileup file, create a new pileup file with a subset
    of the records at the SNP locations only.

    Args:
        all_pileup_file_path: path to a whole-genome pileup file
        all_pileup_file_path: path to a snp pileup file to be created
        snp_set: set of (CHROM, POS) tuples identifying the locations with SNPs
    """
    with open(all_pileup_file_path, "r") as all_pileup_file_object:
        with open(snp_pileup_file_path, "w") as snp_pileup_file_object:
            for pileup_line in all_pileup_file_object:
                current_line_data = pileup_line.rstrip().split()
                seq_id, pos = current_line_data[:2]
                key = (seq_id, int(pos))
                if key in snp_set:
                    snp_pileup_file_object.write(pileup_line)




def write_list_of_snps(file_path, snp_dict):
    """Write out list of snps for all samples to a single file.

    Args:
        file_path : path to snplist file to be written
        snp_dict  : dictionary with key = tuple(CHROM, POS) -> value = list[count, sampleName1, sampleName2, ..., sampleNameN]

    Returns:
        Nothing
    """

    with open(file_path, "w") as snp_list_file_object:
        for key in sorted(snp_dict.iterkeys()):
            snp_list_file_object.write(key[0] + "\t" + str(key[1]))
            values = snp_dict[key]
            for value in values:
                snp_list_file_object.write("\t" + str(value))
            snp_list_file_object.write("\n")


def read_snp_position_list(snp_list_file_path):
    """Read list of snp positions across all samples from the snplist.txt.

    Args:
        snp_list_file_path : path to snplist file to be written

    Returns:
        snp_list  : sorted list of tuple(str(CHROM), int(POS))
    """

    snp_list = list()
    with open(snp_list_file_path, "r") as snp_list_file_object:
        for line in snp_list_file_object:
            chrom, pos = line.split()[0:2]
            snp_list.append((chrom, int(pos)))
    return snp_list


def write_reference_snp_file(reference_file_path, snp_list_file_path,
                             snp_reference_file_path):
    """Write out the snp fasta file for the reference.fasta using the snp
    position file ( snplist.txt).
    """
    #TODO finish documentation
    #TODO actual code is more general than stated. Fix this.

    position_list = [line.split()[0:2] for line in open(snp_list_file_path, "r")]
    match_dict    = SeqIO.to_dict(SeqIO.parse(reference_file_path, "fasta"))

    with open(snp_reference_file_path, "w") as snp_reference_file_object:
        for ordered_id in sorted(match_dict.keys()):
            ref_str = ""
            for chrom_id, pos in position_list:
                if chrom_id == ordered_id:
                    ref_str += match_dict[ordered_id][int(pos) - 1].upper()
            record = SeqRecord(Seq(ref_str), id=ordered_id, description="")
            SeqIO.write([record], snp_reference_file_object, "fasta")


def convert_vcf_files_to_snp_dict(sample_vcf_file_list):
    """convert list of vcf files to a single dict of quality SNPs.

    Args:
        sample_vcf_file_list : list of relative or absolute paths to the sample VCF files

    Returns:
        snp_dict  : dictionary with key = (CHROM, POS) -> value = [count, sampleName1, sampleName2, ..., sampleNameN]

    """

    snp_dict = dict()

    for vcf_file_path in sample_vcf_file_list:

        if not os.path.isfile(vcf_file_path):
            verbose_print("ERROR: Missing VCF file %s" % vcf_file_path)
            continue
        if os.path.getsize(vcf_file_path) == 0:
            verbose_print("ERROR: Empty VCF file %s" % vcf_file_path)
            continue

        verbose_print("Processing VCF file %s" % vcf_file_path)

        sample_name = os.path.basename(os.path.dirname(vcf_file_path))

        with open(vcf_file_path, 'r') as vcf_file_object:
            vcf_reader = vcf.Reader(vcf_file_object)
            for vcf_data_line in vcf_reader:
                key = (vcf_data_line.CHROM, vcf_data_line.POS)
                if not snp_dict.has_key(key):
                    record = [1]
                    snp_dict[key] = record
                else:
                    record = snp_dict[key]
                    record[0] += 1
                record.append(sample_name)

    return snp_dict

