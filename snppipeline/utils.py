#!/usr/bin/env python2.7

from __future__ import print_function
from Bio import SeqIO
import operator
import os
import pprint
import re
import sets
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

    target_timestamp = os.stat(target_file).st_mtime

    for source_file in source_files:
        source_timestamp = os.stat(source_file).st_mtime
        if source_timestamp > target_timestamp:
            return True

    return False



#==============================================================================
#Define functions
#==============================================================================


def get_consensus_base_from_pileup(base, length, data):
    """Call the base for each SNP position

    Description:
    Calls the base based on the pileup data for a SNP position with a given
        length cutoff and with a given reference base.

    Args:
        base: Reference base.
        length: length cutoff.
        data: information from alignment in pileup format.

    Returns:
        consensus_base: Consensus base from alignment.

    Raises:

TODO: all these examples return the value of the 1st argument, can we get some other examples?

    Examples:

    >>> print(get_consensus_base_from_pileup('T',10,',.....,,.,.,...,,,.,..A'))
    T
    >>> print(get_consensus_base_from_pileup('A',10,',.$.....,,.,.,...,,,.,..^+.'))
    A
    >>> print(get_consensus_base_from_pileup('C',10,',.....,,.,.,...,,,.,..A'))
    C
    >>> print(get_consensus_base_from_pileup('T',10,',.$....,,.,.,...,,,.,...'))
    T
    >>> print(get_consensus_base_from_pileup('G',10,',$....,,.,.,...,,,.,...^l.'))
    G
    >>> print(get_consensus_base_from_pileup('G',10,'...T,,.,.,...,,,.,....'))
    G
    >>> print(get_consensus_base_from_pileup('T',10,'....,,.,.,.C.,,,.,..G.'))
    T
    >>> print(get_consensus_base_from_pileup('T',10,'....,,.,.,...,,,.,....^k.'))
    T
    >>> print(get_consensus_base_from_pileup('A',10,'A..T,,.,.,...,,,.,.....'))
    A
    """

    consensus_base = ""
    #TODO consider using a Counter
    base_count_dict = {'.,': 0, 'A': 0, 'C': 0, 'G': 0, 'N': 0, 'T': 0}

    upper_lower_base_set = sets.Set(['A','a','C','c','T','t','G','g','N','n'])
    i = 0
    while i < len(data):
        char = data[i]
        if char in upper_lower_base_set:
            base_count_dict[char.upper()] += 1
        elif char == '.' or char == ',':
            base_count_dict['.,'] += 1
        elif char == '+' or char == '-':
            count_str = ""
            count = 1
            while re.match("\d", data[i + 1]):
                count_str += data[i + 1]
                i += 1
            if count_str != "":
                count = int(count_str)
            i += count
        elif char == '^':
            if data[i + 1] != "." and data[i + 1] != ',':
                i += 1
        i += 1

    consensus_base = max(base_count_dict.iteritems(), key=operator.itemgetter(1))[0]
    if base_count_dict[consensus_base] <= (length/2):
        consensus_base = "-"
    elif consensus_base == ".,":
        consensus_base = base

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


def create_consensus_dict(pileup_file_path):
    """Create a dict based on the information in a pileup file.

    Given a path to a pileup file, create a dict based on the pileup
        data in that file. The dict contains the consensus
        base calls for positions in the reference sequence.

    Args:
        pileup_file_path: full path to a pileup file.

    Returns:
        A dict mapping a key formed from a tuple(sequence, position) to a
        consensus sequence base. For example:

        {('gi|9626243|ref|NC_001416.1|', 46842): 'G',
         ('gi|9626243|ref|NC_001416.1|', 47425): 'T',
         ('gi|9626243|ref|NC_001416.1|', 47893): 'A'}

    Raises:

    Examples:
    >>> consensus_dict = create_consensus_dict("snppipeline/data/lambdaVirusExpectedResults/samples/sample2/reads.snp.pileup")
    >>> consensus_dict[('gi|9626243|ref|NC_001416.1|',3678)]
    'T'
    >>> consensus_dict[('gi|9626243|ref|NC_001416.1|',40984)]
    'A'
    """

    position_value_dict = dict()

    with open(pileup_file_path, "r") as pileup_file_object:
        for pileup_line in pileup_file_object:
            current_line_data = pileup_line.rstrip().split()
            if len(current_line_data) >= 5:   #don't process lines without 5 pieces of information or more
                seq_id, pos, ref_base, depth, bases_string = current_line_data[:5]
                pos = int(pos)
                depth = int(depth)
                position_value_dict[(seq_id, pos)] = get_consensus_base_from_pileup(ref_base, depth, bases_string)

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
            snp_reference_file_object.write(">" + ordered_id + "\n")
            for chrom_id, pos in position_list:
                if chrom_id == ordered_id:
                    snp_reference_file_object.write(match_dict[ordered_id][int(pos)-1].upper())


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

