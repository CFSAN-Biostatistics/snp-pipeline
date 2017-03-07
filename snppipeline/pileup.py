"""
Pileup file parser.  Typical workflow would follow a pattern like this:

import pileup as pileup
caller = pileup.ConsensusCaller(min_freq, min_strand_depth, min_strand_bias)
reader = pileup.Reader(file, min_base_quality, chrom_position_set=None)
for record in reader:
    chrom = record.chrom
    pos = record.position
    ref = record.reference_base
    most_common = record.most_common_base
    alt, fail_reasons = caller.call_consensus(record)
    depth = record.raw_depth
    good_depth = record.good_depth
    fwd_good_depth = record.forward_good_depth
    rev_good_depth = record.reverse_good_depth
    good_ref_depth = record.base_good_depth[ref]
    good_alt_depth = record.base_good_depth[alt]
    fwd_good_ref_depth = record.forward_base_good_depth[ref]
    rev_good_ref_depth = record.reverse_base_good_depth[ref]
    fwd_good_alt_depth = record.forward_base_good_depth[alt]
    rev_good_alt_depth = record.reverse_base_good_depth[alt]
"""

from __future__ import print_function
try:
    from collections import Counter
except ImportError:
    from counter import Counter
import re


#==============================================================================
#Define module regular expressions
#==============================================================================
_re_read_segment_start = re.compile(r"\^.")
_re_indel_count = re.compile(r"[+-](\d+)")


#==============================================================================
#Define classes and functions
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
        forward_good_depth : int
            The total depth of all reads on the forward strand where the base
            quality meets or exceeds the min_base_quality.
        reverse_good_depth : int
            The total depth of all reads on the reverse strand where the base
            quality meets or exceeds the min_base_quality.
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
        >>> r.forward_good_depth
        4
        >>> r.reverse_good_depth
        4
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
        >>> r.forward_good_depth
        4
        >>> r.reverse_good_depth
        4
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
        >>> r.forward_good_depth
        0
        >>> r.reverse_good_depth
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
            self.forward_good_depth = 0
            self.reverse_good_depth = 0
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
        self.forward_good_depth = len(forward_bases_str)
        self.reverse_good_depth = len(reverse_bases_str)
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
        f.close()

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


class ConsensusCaller(object):
    def __init__(self, min_cons_freq, min_cons_strand_depth,
                 min_cons_strand_bias):
        """
        Construct a consensus base caller object with various filter settings
        for subsequent base calling.

        Parameters
        ----------
        min_cons_freq : float
            Mimimum fraction of the high-quality reads supporting the consensus
            to make a consensus call (0.5 - 1.0).  The numerator of this
            fraction is the number of high-quality consensus supporting reads.
            The denominator of this fraction is the total number of high-quality
            reads.
        min_cons_strand_depth : int
            Minimum number of high-quality reads supporting the consensus which
            must be present on both forward and reverse strands separately to
            make a call.
        min_cons_strand_bias : float
            Minimum fraction of the high-quality consensus supporting reads which
            must be present on both forward and reverse strands separately to
            make a call (0.0 - 0.5).  The numerator of this fraction is the
            number of high-quality consensus supporting reads on one strand at
            a time.  The denominator of this fraction is the number of high-
            quality consensus supporting reads.
        """
        self.min_cons_freq = min_cons_freq
        self.min_cons_strand_depth = min_cons_strand_depth
        self.min_cons_strand_bias = min_cons_strand_bias

        self.fail_raw_depth = "RawDpth"
        self.fail_freq = "VarFreq" +  str(int(100 * min_cons_freq))
        self.fail_strand_depth = "StrDpth" + str(min_cons_strand_depth)
        self.fail_strand_bias = "StrBias" +  str(int(100 * min_cons_strand_bias))


    def get_filter_descriptions(self):
        """
        Return a list of tuples of filters and corresponding descriptions.
        This method is intended for use when writing a VCF file header.

        Returns
        -------
        filter : str
            Short string identifying the filter.
        description : str
            Long description of the filter
        """
        return [(self.fail_raw_depth, "No read depth" ),
                (self.fail_freq, "Variant base frequency below %.2f" % self.min_cons_freq),
                (self.fail_strand_depth, "Less than %i variant-supporing reads on at least one strand" % self.min_cons_strand_depth),
                (self.fail_strand_bias, "Fraction of variant supporting reads below %.2f on one strand" % self.min_cons_strand_bias),
               ]

    def call_consensus(self, record):
        """
        Call the consensus base with a list of failed filters for a given
        pileup record using the filters previously configured in the
        constructor.  Calling code should always check the list of
        failed filters before emitting the consensus base.

        Parameters
        ----------
        record : Record
            Parsed pileup record

        Returns
        -------
        consensus_base : str
            Consensus base or '-' if the most common base cannot be determined
        failed_filters : list or str or None
            List of failed filters or None if all filters passed

        Examples
        --------
        >>> r = Record(['ID', 42, 'G', 14, 'aaaaAAAA...,,,', '00001111222333'], 15)
        >>> caller = ConsensusCaller(0.5, 0, 0.0) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('A', None)
        >>> caller = ConsensusCaller(0.6, 0, 0.0) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('A', ['VarFreq60'])
        >>> caller = ConsensusCaller(0.0, 5, 0.0) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('A', ['StrDpth5'])
        >>> r = Record(['ID', 42, 'G', 14, 'aAAAAAAA...,,,', '00001111222333'], 15)
        >>> caller = ConsensusCaller(0.0, 0, 0.2) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('A', ['StrBias20'])
        >>> r = Record(['ID', 42, 'G', 14, 'aaaAAAAA...,,,', '00001111222333'], 15)
        >>> caller = ConsensusCaller(0.0, 4, 0.4) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('A', ['StrDpth4', 'StrBias40'])
        >>> r = Record(['ID', 42, 'G', 14, 'aaaAAA....,,,,', '00011122223333'], 15)
        >>> caller = ConsensusCaller(0.0, 0, 0.0) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('G', None)
        >>> r = Record(['ID', 42, 'g', 14, 'aaaAAA....,,,,', '00011122223333'], 15)
        >>> caller = ConsensusCaller(0.0, 0, 0.0) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('g', None)
        >>> r = Record(['ID', 42, 'g', 0], 15) # zero depth
        >>> caller = ConsensusCaller(0.0, 5, 0.0) # freq, strand_depth, strand_bias
        >>> caller.call_consensus(r)
        ('-', ['RawDpth'])
        """
        consensus_base = record.most_common_base
        if consensus_base == None:
            failed_filters = [self.fail_raw_depth]
            return ('-', failed_filters)

        failed_filters = None

        good_depth = record.good_depth
        good_cons_depth = record.base_good_depth[consensus_base]
        fwd_good_cons_depth = record.forward_base_good_depth[consensus_base]
        rev_good_cons_depth = record.reverse_base_good_depth[consensus_base]

        # Filter: allele minimum frequency
        if good_cons_depth < (good_depth * self.min_cons_freq):
            failed_filters = failed_filters or []
            failed_filters.append(self.fail_freq)

        # Filter: allele minimum depth on each strand
        if fwd_good_cons_depth < self.min_cons_strand_depth or \
           rev_good_cons_depth < self.min_cons_strand_depth:
            failed_filters = failed_filters or []
            failed_filters.append(self.fail_strand_depth)

        # Filter: strand bias
        min_strand_bias_depth = good_cons_depth * self.min_cons_strand_bias
        if fwd_good_cons_depth < min_strand_bias_depth or \
           rev_good_cons_depth < min_strand_bias_depth:
            failed_filters = failed_filters or []
            failed_filters.append(self.fail_strand_bias)

        # Keep the reference lowercase if it was lowercase in the pileup
        if consensus_base == record.reference_base.upper():
            consensus_base = record.reference_base

        return (consensus_base, failed_filters)


