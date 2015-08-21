#!/usr/bin/env python2.7

"""
Pileup file parser.  Typical workflow would follow a pattern like this:

import pileup as pileup
reader = pileup.Reader(file, min_base_quality, chrom_position_set=None)
for record in reader:
    chrom = record.chrom
    pos = record.position
    ref = record.reference_base
    alt = record.most_common_base  # how to handle ties?
    depth = record.raw_depth
    good_depth = record.good_depth
    good_ref_depth = record.base_good_depth[ref]
    good_alt_depth = record.base_good_depth[alt]
    fwd_good_ref_depth = record.forward_base_good_depth[ref]
    rev_good_ref_depth = record.reverse_base_good_depth[ref]
    fwd_good_alt_depth = record.forward_base_good_depth[alt]
    rev_good_alt_depth = record.reverse_base_good_depth[alt]
"""

from __future__ import print_function
from collections import Counter
import re


#==============================================================================
#Define module regular expressions
#==============================================================================
_re_read_segment_start = re.compile(r"\^.")
_re_indel_count = re.compile(r"[+-](\d+)")


#==============================================================================
#Define functions
#==============================================================================

class Record(object):
    def __init__(self, arg, min_base_quality):
        """
        Parse and analyze a pileup record.

        Parameters
        ----------
        arg : str or list
            If str, a line of text from a pileup file.
            If list, a split line of text from a pileup file.
        min_base_quality : int
            Minimum base quality phred score for a base to be counted as valid 
            read depth.

        Attributes
        ----------
        chrom : str
            Chromosome.
        position : int
            One-based position of this record.
        reference_base : str
            The reference base for the chromosome at this position.  
            The upper/lower case of the base in the pileup file is preserved.
        most_common_base : str or None
            The most commonly occurring base after discarding low-quality bases.
            Always uppercase if there is any good depth, otherwise None.
        raw_depth : int
            The total depth of reads regardless of base quality.
        good_depth : int
            The total depth of all reads where the base quality meets or
            exceeds the min_base_quality.
        base_good_depth : Counter
            Count of depth per base (A,C,G,T,N,-) on both forward and reverse 
            strands combined where the base quality meets or exceeds the 
            min_base_quality.
        forward_base_good_depth : Counter
            Count of depth per base (A,C,G,T,N,-) on forward strand only where the
            base quality meets or exceeds the min_base_quality.
        reverse_base_good_depth : Counter
            Count of depth per base (A,C,G,T,N,-) on reverse strand only where the
            base quality meets or exceeds the min_base_quality.

        Examples
        -------
        >>> r = Record(['NC_011149.1', 42, 'G', 9, 'aaAaA+6TAAGAG..+5AAGAG.,', '21G1G-111'], 15)
        >>> r.chrom
        'NC_011149.1'
        >>> r.position
        42
        >>> r.reference_base
        'G'
        >>> r.most_common_base
        'A'
        >>> r.raw_depth
        9
        >>> r.good_depth
        8
        >>> r.base_good_depth['A'], r.base_good_depth['T'], r.base_good_depth['G'], r.base_good_depth['C'],
        (5, 0, 3, 0)
        >>> fwd = r.forward_base_good_depth
        >>> fwd['A'], fwd['T'], fwd['G'], fwd['C'],
        (2, 0, 2, 0)
        >>> rev = r.reverse_base_good_depth
        >>> rev['A'], rev['T'], rev['G'], rev['C'],
        (3, 0, 1, 0)
        >>>
        >>>
        >>> r = Record('NC_011149.1\t42\tG\t9\taaAaA+6TAAGAG..+5AAGAG.,\t21G1G-111', 15)
        >>> r.chrom
        'NC_011149.1'
        >>> r.position
        42
        >>> r.reference_base
        'G'
        >>> r.most_common_base
        'A'
        >>> r.raw_depth
        9
        >>> r.good_depth
        8
        >>> r.base_good_depth['A'], r.base_good_depth['T'], r.base_good_depth['G'], r.base_good_depth['C'],
        (5, 0, 3, 0)
        >>> fwd = r.forward_base_good_depth
        >>> fwd['A'], fwd['T'], fwd['G'], fwd['C'],
        (2, 0, 2, 0)
        >>> rev = r.reverse_base_good_depth
        >>> rev['A'], rev['T'], rev['G'], rev['C'],
        (3, 0, 1, 0)
        >>>
        >>>
        >>>
        >>> r = Record(['ID', 628640, 'A', 20, '**.,,.,.............', '22E?;9HF;H8EDGHHI?GH'], 15)
        >>> r.base_good_depth['-']
        2
        >>> r = Record(['gi|197247352|ref|NC_011149.1|', '4663812', 'T', '0'], 15)
        >>> r.good_depth
        0
        >>> print(r.most_common_base)
        None
        """
        if isinstance(arg, str):
            self._init_from_line(arg, min_base_quality)
        elif isinstance(arg, list):
            self._init_from_split_line(arg, min_base_quality)
        else:
            raise TypeError("Pileup record can be constructed from a string or list.  %s is not supported." % str(type(arg)))

    def _init_from_line(self, line, min_base_quality):
        """
        Method for internal use to parse and analyze a line of text from a 
        pileup file.

        Parameters
        ----------
        line : str
            A line of text from a pileup file.
        min_base_quality : int
            Minimum base quality phred score for a base to be counted as valid 
            read depth.
        """
        split_line = line.rstrip().split()
        self._init_from_split_line(split_line, min_base_quality)

    def _init_from_split_line(self, split_line, min_base_quality):
        """
        Method for internal use to parse and analyze a previously split line 
        from a pileup file.

        Parameters
        ----------
        split_line : list
            Split line of text from a pileup file.
        min_base_quality : int
            Minimum base quality phred score for a base to be counted as valid 
            read depth.
        """
        self.chrom = split_line[0]
        self.position = int(split_line[1])
        self.reference_base = split_line[2]
        self.raw_depth = int(split_line[3])
        if self.raw_depth == 0 or len(split_line) < 5:
            self.good_depth = 0
            self.base_good_depth = Counter()
            self.forward_base_good_depth = Counter()
            self.reverse_base_good_depth = Counter()
            self.most_common_base = None
            return

        bases_str = split_line[4]
        quality_str = split_line[5]
        bases_str = self._strip_unwanted_base_patterns(bases_str)
        #if len(bases_str) != len(quality_str):
        #    print("len(bases_str) != len(quality_str)")
        #    print("len(bases_str) =", len(bases_str), "len(quality_str) =", len(quality_str))
        #    print("pos", self.position)
        #    print("    orig bases =", split_line[4])
        #    print("stripped bases =", bases_str)
        #    print("quality string =", quality_str)

        # Discard low quality bases
        quality_values = [ord(c) - 33 for c in quality_str]
        bases_list = [b for b, q in zip(bases_str, quality_values) if q >= min_base_quality]
        bases_str = ''.join(bases_list)

        self.good_depth = len(bases_list)

        # Substitute - for * (deletion supported by this read)
        bases_str = bases_str.replace('*', '-')

        # Substitute reference bases markers with letter
        bases_str = bases_str.replace('.', self.reference_base.upper())
        bases_str = bases_str.replace(',', self.reference_base.lower())

        # Get counts of bases regardless of strand
        self.base_good_depth = Counter(bases_str.upper())
        if self.good_depth < 1:
            self.most_common_base = None
        else:
            self.most_common_base = self.base_good_depth.most_common(1)[0][0]

        # Get counts of bases on each strand
        forward_bases_str = ''.join([c for c in bases_str if c <= 'Z'])
        reverse_bases_str = ''.join([c for c in bases_str if c >= 'a'])
        self.forward_base_good_depth = Counter(forward_bases_str)
        self.reverse_base_good_depth = Counter(reverse_bases_str.upper())


    @staticmethod # Doesn't use self
    def _strip_unwanted_base_patterns(bases_str):
        """
        Method for internal use to remove unwanted special patterns from
        the string of pileup bases.

        Parameters
        ----------
        bases_str : str
            Raw string of pileup bases.

        Returns
        ------
        bases_str : str
            Modified string of pileup bases without unwanted patterns.

        Examples
        --------
        >>> print(Record._strip_unwanted_base_patterns(".,.actg,,,"))
        .,.actg,,,
        >>> print(Record._strip_unwanted_base_patterns("^K.,.^Fa,,,^K"))
        .,.a,,,
        >>> print(Record._strip_unwanted_base_patterns("$.,.$*$*,,,*"))
        .,.**,,,*
        >>> print(Record._strip_unwanted_base_patterns(".,.+10AAAAAAAAAAa,,,"))
        .,.a,,,
        >>> print(Record._strip_unwanted_base_patterns("+2TT.,.+10AAAAAAAAAAa,,,+2GC"))
        .,.a,,,
        >>> print(Record._strip_unwanted_base_patterns(".,.-10AAAAAAAAAAa,,,"))
        .,.a,,,
        >>> print(Record._strip_unwanted_base_patterns("-2TT.,.-10AAAAAAAAAAa,,,-2GC"))
        .,.a,,,
        >>> print(Record._strip_unwanted_base_patterns("^Kc-2TT..$a+10AAAAAAAAAAa,,*,-2GC"))
        c..aa,,*,
        """
        # Remove all read segment start markers
        bases_str, _ = _re_read_segment_start.subn('', bases_str)

        # Remove all indel markers backwards from end, to preserve indexes
        matches = [m for m in _re_indel_count.finditer(bases_str)]
        for match in reversed(matches):
            num_bases = int(match.group(1))
            drop_start_idx = match.start()
            drop_end_idx = match.end() + num_bases
            bases_str = bases_str[:drop_start_idx] + bases_str[drop_end_idx:]

        # Remove all end of read segment merkers
        bases_str = bases_str.replace('$', '')

        return bases_str


class Reader(object):
    def __init__(self, file_path, min_base_quality, chrom_position_set=None):
        """
        Contruct a reader object with the capability to read and parse lines 
        from a pileup file.

        Parameters
        ----------
        file_path : str
            Path to the pileup file to be parsed.
        min_base_quality : int
            Minimum base quality phred score for a base to be counted as valid 
            read depth.
        chrom_position_set : set of (str, int), optional
            Tuples of (chromosome name, position) identifying the positions to
            be parsed in the file.  If not specified, all positions will be
            parsed.
        """
        self.file_path = file_path
        self.min_base_quality = min_base_quality
        self.chrom_position_set = chrom_position_set
        # open and close the file to make sure it works
        f = open(file_path)
        f.close

    def __iter__(self):
        """
        Read lines from the pileup file and convert them to parsed records.

        Returns
        -------
        record : Record
            Returns a parsed pileup record each time this method is called.
        """
        with open(self.file_path, "r") as f:
            if self.chrom_position_set == None:
                for line in f:
                    record = Record(line, self.min_base_quality)
                    yield record
            else:
                for line in f:
                    split_line = line.rstrip().split()
                    chrom, pos = split_line[:2]
                    key = (chrom, int(pos))
                    if key in self.chrom_position_set:
                        record = Record(split_line, self.min_base_quality)
                        yield record

