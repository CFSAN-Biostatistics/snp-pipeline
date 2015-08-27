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

def get_consensus_base_from_pileup(record, min_cons_freq, min_cons_strand_depth, min_cons_strand_freq):
    """
    Call the consensus base for each SNP position with the specified thresholds.

    Parameters
    ----------
    record : pileup.record
        Parsed pileup record
    min_cons_freq : float
        Mimimum fraction of the high-quality reads supporting the consensus to
        make a consensus call (0.0 - 1.0).  The numerator of this fraction
        is the number of high-quality consensus supporting reads.  The
        denominator of this fraction is the total number of high-quality reads.
    min_cons_strand_depth : int
        Minimum number of high-quality reads supporting the consensus which must
        be present on both forward and reverse strands separately to make a
        call.
    min_cons_strand_freq : float
        Minimum fraction of the high-quality consensus supporting reads which
        must be present on both forward and reverse strands separately to make
        a call (0.0 - 0.5).  The numerator of this fraction is the number of
        high-quality consensus supporting reads on one strand at a time.  The
        denominator of this fraction is the number of high-quality consensus
        supporting reads.  

    Returns
    -------
    consensus_base : str
        Consensus base.

    Examples
    --------
    >>> r = pileup.Record(['ID', 42, 'G', 14, 'aaaaAAAA...,,,', '00001111222333'], 15)
    >>> get_consensus_base_from_pileup(r, 0.5, 0, 0.0)
    'A'
    >>> get_consensus_base_from_pileup(r, 0.6, 0, 0.0)
    '-'
    >>> get_consensus_base_from_pileup(r, 0.0, 5, 0.0)
    '-'
    >>> r = pileup.Record(['ID', 42, 'G', 14, 'aAAAAAAA...,,,', '00001111222333'], 15)
    >>> get_consensus_base_from_pileup(r, 0.0, 0, 0.13)
    '-'
    >>> r = pileup.Record(['ID', 42, 'G', 14, 'aaaAAA....,,,,', '00011122223333'], 15)
    >>> get_consensus_base_from_pileup(r, 0.0, 0, 0.0)
    'G'
    >>> r = pileup.Record(['ID', 42, 'g', 14, 'aaaAAA....,,,,', '00011122223333'], 15)
    >>> get_consensus_base_from_pileup(r, 0.0, 0, 0.0)
    'g'
    >>> r = pileup.Record(['ID', 42, 'g', 0], 15)
    >>> get_consensus_base_from_pileup(r, 0.0, 0, 0.0)
    '-'
    """
    consensus_base = record.most_common_base
    if consensus_base == None:
        return "-"
    good_depth = record.good_depth
    good_cons_depth = record.base_good_depth[consensus_base]
    fwd_good_cons_depth = record.forward_base_good_depth[consensus_base]
    rev_good_cons_depth = record.reverse_base_good_depth[consensus_base]

    # Filter: allele minimum frequency
    if good_cons_depth < (good_depth * min_cons_freq):
        return "-"
    # Filter: allele minimum depth on each strand
    if fwd_good_cons_depth < min_cons_strand_depth:
        return "-"
    if rev_good_cons_depth < min_cons_strand_depth:
        return "-"
    # Filter: strand bias
    if fwd_good_cons_depth < (good_cons_depth * min_cons_strand_freq):
        return "-"
    if rev_good_cons_depth < (good_cons_depth * min_cons_strand_freq):
        return "-"

    # Keep the reference lowercase if it was lowercase in the pileup
    if consensus_base == record.reference_base.upper():
        consensus_base = record.reference_base 

    return consensus_base


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


def create_consensus_dict(pileup_file_path, min_base_qual, min_cons_freq, 
                          min_cons_strand_depth, min_cons_strand_bias, 
                          chrom_position_set=None):
    """Create a consensus snp call dict based on the information in a pileup
    file.

    Given a path to a pileup file, create a dict based on the pileup
    data in that file. The dict contains the consensus base calls for
    positions in the reference sequence.

    Parameters
    ----------
    pileup_file_path : str
        Full path to a pileup file.
    min_base_qual : int
        Mimimum base quality score to count a read. All other snp filters
        take effect after the low-quality reads are discarded.
    min_cons_freq : float
        Mimimum fraction of reads that must agree (0.0 - 1.0)
    min_cons_freq : float
        Consensus frequency. Mimimum fraction of high-quality reads
        supporting the consensus to make a call.
    min_cons_strand_depth : int
        Consensus strand depth. Minimum number of high-quality reads 
        supporting the consensus which must be present on both the
        forward and reverse strands to make a call.
    min_cons_strand_bias : float
        Strand bias. Minimum fraction of the high-quality 
        consensus-supporting reads which must be present on both the 
        forward and reverse strands to make a call. The numerator of this
        fraction is the number of high-quality consensus-supporting reads
        on one strand.  The denominator is the total number of high-quality
        consensus-supporting reads on both strands combined.
    chrom_position_set : set of (str, int), optional
        Tuples of (chromosome name, position) identifying the positions to
        be called in the pileup file.  If not specified, all positions will be
        called.

    Returns
    -------
    snp_dict : dict
        A dict mapping a key formed from a tuple(chromosome name, position) to a
        consensus sequence base. For example:
        {('gi|9626243|ref|NC_001416.1|', 46842): 'G',
         ('gi|9626243|ref|NC_001416.1|', 47425): 'T',
         ('gi|9626243|ref|NC_001416.1|', 47893): 'A'}

    Examples
    --------
    >>> consensus_dict = create_consensus_dict("snppipeline/data/lambdaVirusExpectedResults/samples/sample2/reads.snp.pileup", 0, 0.5, 0, 0)
    >>> consensus_dict[('gi|9626243|ref|NC_001416.1|',3678)]
    'T'
    >>> consensus_dict[('gi|9626243|ref|NC_001416.1|',40984)]
    'A'
    """

    position_value_dict = dict()
    reader = pileup.Reader(pileup_file_path, min_base_qual, chrom_position_set)
    for record in reader:
        chrom = record.chrom
        pos = record.position
        position_value_dict[(chrom, pos)] = get_consensus_base_from_pileup(record, min_cons_freq, min_cons_strand_depth, min_cons_strand_bias)

    return position_value_dict


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

