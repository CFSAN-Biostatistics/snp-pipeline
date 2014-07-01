#!/usr/bin/env python2.7

from __future__ import print_function
from Bio import SeqIO
import operator
import os
import pprint
import re
import sets
import subprocess
import sys
import vcf


#==============================================================================
#Prep work
#Note use of filter on list_of_sample_directories to remove blank lines.
#==============================================================================

verbose        = False
verbose_print  = print         if verbose else lambda *a, **k: None
verbose_pprint = pprint.pprint if verbose else lambda *a, **k: None


#==============================================================================
#Define functions
#==============================================================================

def pileup_wrapper(args):
    """Wraps pileup to use multiple arguments with multiprocessing package.
    """
    return pileup(*args)


def pileup(sample_dir_path, options_dict):
    """Run samtools to generate pileup for one sample.

    Description:
    Generate pileup file, using snplist file and the reference fasta file.

    Args:
        sample_dir_path: Full path to directory for a sample.
        options_dict: Dictionary containing command-line options as per
            documention for run_snp_pipeline.
    """

    verbose_print('Generating pileup file ' + options_dict['pileupFileName'] +
                  ' in '+sample_dir_path)
    pileup_file_path  = os.path.abspath(os.path.join(sample_dir_path, options_dict['pileupFileName']))

    # SAMtools ignores the snplist when the path is relative, make it absolute here:
    snplist_file_path = os.path.abspath(os.path.join(options_dict['mainPath'], options_dict['snplistFileName']))

    #TODO - allow for use of previously done pileup via command line argument?
    if os.path.isfile(pileup_file_path):
        verbose_print('Removing old pileup file '+pileup_file_path)
        os.remove(pileup_file_path)

    # SAMtools fails when the reference file is "./reference/lambdaVirus.fasta", make it absolute here:
    reference_file_path = os.path.abspath(os.path.join(options_dict['mainPath'], options_dict['Reference']))
    bam_file_path = os.path.abspath(os.path.join(sample_dir_path, options_dict['bamFileName']))
    command_line = (
        'samtools mpileup ' +
            '-l ' + snplist_file_path +
            ' -f ' + reference_file_path + 
            ' ' + bam_file_path +
            ' > ' + pileup_file_path
    )
    
    verbose_print('Executing: '+command_line)
    #TODO review return values and clean up this next bit of code.
    #  see for details: https://docs.python.org/2/library/subprocess.html#replacing-os-system
    try:
        return_code = subprocess.call(command_line, cwd=sample_dir_path, shell=True)
        if return_code < 0:
            sys.stderr.write("Child was terminated by signal: " + str(return_code))
        else:
            verbose_print("Child returned: " + str(return_code))
    except OSError as e:
        sys.stderr.write("Execution failed: " + str(e))


    verbose_print("samtools return code: ", return_code)
    if not os.path.isfile(pileup_file_path):
        print('Pileup file not created: '+pileup_file_path)

    verbose_print('pileup function exit')

    return pileup_file_path

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

    #TODO move the char.upper out of the loop.
    #TODO move the set construction out of the loop, and we only need 4 elements in the set
    i = 0
    while i < len(data):
        char = data[i]
        if char in sets.Set(['A','a','C','c','T','t','G','g','N','n']):
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
    if base_count_dict[consensus_base] <= (int(length)/2):
        consensus_base = "-"
    elif consensus_base == ".,":
        consensus_base = base

    return consensus_base


def create_consensus_dict(pileup_file_path):
    """Create a dict based on the information in a pileup file.

    Given a path to a pileup file, create a dict based on the pileup
        data in that file. The dict contains the consensus
        base calls for positions in the reference sequence.

    Args:
        pileup_file_path: full path to a pileup file.

    Returns:
        A dict mapping a key formed from a sequence and position to a
        consensus sequence base. For example:

        {'gi|9626243|ref|NC_001416.1|:46842': 'G',
         'gi|9626243|ref|NC_001416.1|:47425': 'T',
         'gi|9626243|ref|NC_001416.1|:47893': 'A'}

    Raises:

    Examples:
    >>> consensus_dict = create_consensus_dict("snppipeline/data/lambdaVirusExpectedResults/samples/sample2/reads.pileup")
    >>> consensus_dict['gi|9626243|ref|NC_001416.1|:3678']
    'T'
    >>> consensus_dict['gi|9626243|ref|NC_001416.1|:40984']
    'A'
    """

    position_value_dict = dict()

    with open(pileup_file_path, "r") as pileup_file_object:
        for pileup_line in pileup_file_object:
            current_line_data = pileup_line.rstrip().split()
            #TODO should the test below be >= 5 ???
            if len(current_line_data) > 5:   #don't process lines without 5 pieces of information or more
                #TODO seq_id, pos, ref_base, depth, bases_string = current_line_data[:5]
                #TODO position_value_dict[seq_id + ":" + pos] = get_consensus_base_from_pileup(ref_base, depth, bases_string)
                position_value_dict[current_line_data[0] + ":" + current_line_data[1]] = get_consensus_base_from_pileup(current_line_data[2], current_line_data[3], current_line_data[4])

    return position_value_dict


def write_list_of_snps(file_path, snp_dict):
    """Write out list of snps for all samples to a single file.
    """
    #TODO finish documentation

    with open(file_path, "w") as snp_list_file_object:
        for key in sorted(snp_dict.iterkeys()):
            snp_list_file_object.write(key)
            values = snp_dict[key]
            for value in values:
                snp_list_file_object.write("\t" + str(value))
            snp_list_file_object.write("\n")


def write_reference_snp_file(reference_file_path, snp_list_file_path,
                             snp_reference_file_path):
    """Write out the snp fasta file for the reference.fasta using the snp
    position file ( snplist.txt).
    """
    #TODO finish documentation
    #TODO actual code is more general than stated. Fix this.

    position_list = [line.split() for line in open(snp_list_file_path, "r")]
    match_dict    = SeqIO.to_dict(SeqIO.parse(reference_file_path, "fasta"))

    with open(snp_reference_file_path, "w") as snp_reference_file_object:
        for ordered_id in sorted(match_dict.keys()):
            snp_reference_file_object.write(">" + ordered_id + "\n")
            for position in position_list:
                chrom_id, PosID = position[0:2]
                if chrom_id == ordered_id:
                    snp_reference_file_object.write(match_dict[ordered_id][int(PosID)-1].upper())


def convert_vcf_files_to_snp_dict(list_of_sample_directories, options_dict):
    """convert list of vcf files to a single dict of quality SNPs.

    We use several criteria from options_dict to evaluate SNPs.

    Note use of get to cleanly handle case of missing key w/o exception.
    """

    snp_dict = dict()

    for sample_directory in list_of_sample_directories:

        #TODO os.path.split(sample_directory) ???
        sample_name = sample_directory.split(os.sep)[-1]

        with open(os.path.join(sample_directory, "var.flt.vcf"), 'r') as vcf_file_object:
            vcf_reader = vcf.Reader(vcf_file_object)
            for vcf_data_line in vcf_reader:
                if not snp_dict.has_key(vcf_data_line.CHROM + "\t" + str(vcf_data_line.POS)):
                    record = [1]
                    record.append(sample_name)
                    snp_dict[vcf_data_line.CHROM + "\t" + str(vcf_data_line.POS)] = record
                else:
                    record = snp_dict[vcf_data_line.CHROM + "\t" + str(vcf_data_line.POS)]
                    record[0] += 1
                    record.append(sample_name)

    return snp_dict

