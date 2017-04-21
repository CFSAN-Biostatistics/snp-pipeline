"""
Utility functions for fastq files.
"""

from __future__ import print_function

import collections
import glob
import gzip
import os
import re


def list_fastq_files(directory):
    """Return a list of fastq files in a directory.

    The files match any of the suffixes: .fastq, .fq, .fastq.gz, fq.gz.

    Parameters
    ----------
    directory : str
        Directory path.

    Returns
    -------
    files : list
        Sorted list of fastq files in the directory.

    Examples
    --------
    # Setup tests
    >>> import tempfile
    >>> import shutil
    >>> temp_dir = tempfile.mkdtemp(prefix="tmp.fastq.test", dir="./")
    >>> f = open(os.path.join(temp_dir, "aaa_2.fastq"), 'w'); f.close()
    >>> f = open(os.path.join(temp_dir, "aaa_1.fastq"), 'w'); f.close()
    >>> f = open(os.path.join(temp_dir, "bbb.fastq.gz"), 'w'); f.close()
    >>> f = open(os.path.join(temp_dir, "ccc.fq"), 'w'); f.close()
    >>> f = open(os.path.join(temp_dir, "ddd.fq.gz"), 'w'); f.close()
    >>> f = open(os.path.join(temp_dir, "eee.notfastq"), 'w'); f.close()
    >>> f = open(os.path.join(temp_dir, "fff.fastq.not"), 'w'); f.close()
    >>> files = list_fastq_files(temp_dir)
    >>> [os.path.basename(file) for file in files]
    ['aaa_1.fastq', 'aaa_2.fastq', 'bbb.fastq.gz', 'ccc.fq', 'ddd.fq.gz']

    # Clean up
    >>> shutil.rmtree(temp_dir)
    """
    fastq_files = list()

    for suffix in ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]:
        fastq_glob_pattern = os.path.join(directory, suffix)
        fastq_files.extend(glob.glob(fastq_glob_pattern))

    fastq_files.sort()
    return fastq_files


# Mapping of Illumina flowcell last 4 characters to instrument type
ILLUMINA_FLOWCELL_INSTRUMENT_TYPE_DICT = {
    "AAXX" : "Genome Analyzer",
    "ACXX" : "HiSeq",   # HiSeq High-Output v3
    "ADXX" : "HiSeq",   # HiSeq RR v1
    "AFXX" : "NextSeq", # NextSeq Mid-Output
    "AGXX" : "NextSeq", # NextSeq High-Output
    "ALXX" : "HiSeqX",
    "AMXX" : "HiSeq",   # HiSeq RR v2
    "ANXX" : "HiSeq",   # HiSeq High-Output v4
    "BBXX" : "HiSeq",   # See SRR3619967
    "BCXX" : "HiSeq",   # HiSeq v1.5 or HiSeq RR v2
    "BGXY" : "NextSeq", # NextSeq High-Output
}

def flowcell_to_instrument_type(flow_cell):
    """Given a flow cell id, determine the sequencing instrument type.

    Parameters
    ----------
    flow_cell : str
        Flowcell identifier

    Returns
    -------
    instrument_type : str or None
        A string identifying the instrument type if recognized, otherwise None.

    Examples
    --------
    >>> flowcell_to_instrument_type("fcAAXX")
    'Genome Analyzer'
    >>> flowcell_to_instrument_type("fcALXX")
    'HiSeqX'
    >>> flowcell_to_instrument_type("fc") is None
    True

    References
    ----------
    https://www.biostars.org/p/198143/
    """
    last4 = flow_cell[-4:]
    return ILLUMINA_FLOWCELL_INSTRUMENT_TYPE_DICT.get(last4)


# Regular expressions used to determine the Illumina instrument type from an instrument name.
ILLUMINA_MISEQ_PATTERN = "((HWI-)?M[0-9]{5}(R|L1)?)$"
ILLUMINA_HISEQ_PATTERN = "((HWI-)?(([DJK][0-9]{5})|(ST[0-9]{3,4})))$"
ILLUMINA_NEXTSEQ_PATTERN = "(NS[0-9]{6})$"
ILLUMINA_MISEQ_REGEX = re.compile(ILLUMINA_MISEQ_PATTERN)
ILLUMINA_HISEQ_REGEX = re.compile(ILLUMINA_HISEQ_PATTERN)
ILLUMINA_NEXTSEQ_REGEX = re.compile(ILLUMINA_NEXTSEQ_PATTERN)

def instrument_name_to_instrument_type(instrument_name):
    """Given an instrument_name, determine the sequencing instrument type.

    Parameters
    ----------
    instrument_name : str
        Instrument name

    Returns
    -------
    instrument_type : str or None
        A string identifying the instrument type if recognized, otherwise None.

    Examples
    --------
    @HWI-Mxxxx or @Mxxxx - MiSeq
    @HWUSI - GAIIx
    @HWI-Dxxxx - HiSeq 2000/2500
    @Jxxxxx - HiSeq 3000/4000
    @Kxxxxx - HiSeq 3000/4000
    @Nxxxx - NextSeq 500/550
    @NSxxxxxx - NextSeq

    >>> instrument_name_to_instrument_type(None) is None
    True
    >>> instrument_name_to_instrument_type("") is None
    True
    >>> instrument_name_to_instrument_type("HWI-M00229")
    'MiSeq'
    >>> instrument_name_to_instrument_type("M00229")
    'MiSeq'
    >>> instrument_name_to_instrument_type("HWI-M00229R")
    'MiSeq'
    >>> instrument_name_to_instrument_type("M00229R")
    'MiSeq'
    >>> instrument_name_to_instrument_type("M00229L1")
    'MiSeq'
    >>> instrument_name_to_instrument_type("M00229L2") is None
    True
    >>> instrument_name_to_instrument_type("M00229Z") is None
    True
    >>> instrument_name_to_instrument_type("HWI-ST1029")
    'HiSeq'
    >>> instrument_name_to_instrument_type("HWI-ST741")
    'HiSeq'
    >>> instrument_name_to_instrument_type("ST741")
    'HiSeq'
    >>> instrument_name_to_instrument_type("NS500287")
    'NextSeq'
    >>> instrument_name_to_instrument_type("HWUSIxxx")
    'GAIIx'
    >>> instrument_name_to_instrument_type("Unknown") is None
    True

    References
    ----------
    https://www.biostars.org/p/198143/
    """
    if not instrument_name:
        return None
    if ILLUMINA_MISEQ_REGEX.match(instrument_name):
        return "MiSeq"
    if ILLUMINA_HISEQ_REGEX.match(instrument_name):
        return "HiSeq"
    if ILLUMINA_NEXTSEQ_REGEX.match(instrument_name):
        return "NextSeq"
    if instrument_name.startswith("HWUSI"):
        return "GAIIx"
    return None


# Regular expression used to parse Illumina fastq sequence id lines.
# This is an abbreviated regular expression, there may be more tags after these.
ILLUMINA_INSTRUMENT_PATTERN = "([A-Z][A-Z0-9\-]*)"   # Instrument name : starts with A-Z and contains A-Z,0-9,'-'
ILLUMINA_RUN_PATTERN = "([0-9]+)" # Run number : sequence of digits
ILLUMINA_FLOWCELL_PATTERN = "([a-zA-Z0-9\-]*)" # Flowcell : contains A-Z,0-9,'-'
ILLUMINA_LANE_PATTERN = "([0-9]{1,2})" # Lane number : one or two digits
ILLUMINA_TILE_PATTERN = "([0-9]+)" # Tile number : digits
ILLUMINA_XYPOS_PATTERN = "([0-9]+)" # X pos or y pos : digits

ILLUMINA_FASTQ_SEQ_ID_ENDS_WITH_PATTERN = ILLUMINA_FLOWCELL_PATTERN + "[:_]" + \
                                          ILLUMINA_LANE_PATTERN + "[:_]" + \
                                          ILLUMINA_TILE_PATTERN + "[:_]" + \
                                          ILLUMINA_XYPOS_PATTERN + "[:_]" + \
                                          ILLUMINA_XYPOS_PATTERN

# @<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>#<index sequence>
ILLUMINA_FASTQ_SEQ_ID_REGEX1 = re.compile("@" + ILLUMINA_FASTQ_SEQ_ID_ENDS_WITH_PATTERN)

# @<sample id> <flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>#<index sequence>
ILLUMINA_FASTQ_SEQ_ID_REGEX2 = re.compile("@[SE]RR[A-Z0-9\-.]+[ _]" + ILLUMINA_FASTQ_SEQ_ID_ENDS_WITH_PATTERN)

# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
ILLUMINA_FASTQ_SEQ_ID_REGEX3 = re.compile("@" + ILLUMINA_INSTRUMENT_PATTERN + "[:_]" + ILLUMINA_RUN_PATTERN + "[:_]" + ILLUMINA_FASTQ_SEQ_ID_ENDS_WITH_PATTERN)

# @<sample id> <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
ILLUMINA_FASTQ_SEQ_ID_REGEX4 = re.compile("@[SE]RR[A-Z0-9\-.]+[ _]" + ILLUMINA_INSTRUMENT_PATTERN + "[:_]" + ILLUMINA_RUN_PATTERN + "[:_]" + ILLUMINA_FASTQ_SEQ_ID_ENDS_WITH_PATTERN)

# Named tuple to contain fastq metadata
FastqSeqTags = collections.namedtuple("FastqSeqTags", "platform instrument_type instrument run flow_cell lane")

def parse_seqid_line(seqid_line):
    """Examine a fastq sequence id line and extract various metadata tags.

    Looks at a sequence id line of a fastq file to extract the instrument,
    flowcell, and lane if possible.  There is one of these lines preceeding
    every read in a fastq file and the information is highly redundant.

    Illumina fastq files have a variety of sequence id formats:
    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index sequence>
    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>       <read>:<is filtered>:<control number>:<index sequence>
    @<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>#<index sequence>

    Parameters
    ----------
    seqid_line : str
        A sequence identifier line from a fastq file.

    Returns
    -------
    FastqSeqTags named-tuple with the following named elements or None if the line cannot be parsed
        platform : str
            Currently, only Illumina data is supported and this will always be "illumina" if the line can be parsed
        instrument_type : str or None
            The instrument type can _sometimes_ be determined from the flowcell or instrument name.
            This is not well tested.
        instrument : str
            Instrument name
        run : str
            Run number
        flow_cell : str
            Flowcell with leading zeros removed
        lane : str
            Lane number

    Examples
    --------
    # Not Illumina format
    >>> parse_seqid_line("@SRR1206159_1/1") is None
    True

    # fastq-dump default settings
    >>> parse_seqid_line("@SRR498276.1 HWI-M00229:9:000000000-A1474:1:1:15012:1874 length=151")
    FastqSeqTags(platform='illumina', instrument_type='MiSeq', instrument='HWI-M00229', run='9', flow_cell='A1474', lane='1')

    # something other than zeros in flowcell prefix
    >>> parse_seqid_line("@SRR498276.1 HWI-M00229:9:000100000-A1474:1:1:15012:1874 length=151")
    FastqSeqTags(platform='illumina', instrument_type='MiSeq', instrument='HWI-M00229', run='9', flow_cell='000100000-A1474', lane='1')

    # Lee's fastq-dump defline-seq format
    >>> parse_seqid_line("@SRR498423_HWI-M00229:7:000000000-A0WG8:1:1:12203:2225/1")
    FastqSeqTags(platform='illumina', instrument_type='MiSeq', instrument='HWI-M00229', run='7', flow_cell='A0WG8', lane='1')

    # HiSeq
    >>> parse_seqid_line("@HWI-ST741:189:C0GU5ACXX:8:1101:1219:1953 1:N:0:")
    FastqSeqTags(platform='illumina', instrument_type='HiSeq', instrument='HWI-ST741', run='189', flow_cell='C0GU5ACXX', lane='8')

    # NextSeq
    >>> parse_seqid_line("@NS500287:189:FLOW0000:8:1101:1219:1953 1:N:0:")
    FastqSeqTags(platform='illumina', instrument_type='NextSeq', instrument='NS500287', run='189', flow_cell='FLOW0000', lane='8')

    # ERR178930 pattern with "HWI-" instrument prefix
    >>> parse_seqid_line('@ERR178930.1 HWI-ST322_0214_"AC0HTNACXX":8:1101:1555:2158#ATCACG length=101')
    FastqSeqTags(platform='illumina', instrument_type='HiSeq', instrument='HWI-ST322', run='0214', flow_cell='AC0HTNACXX', lane='8')

    # ERR178930 pattern with MiSeq instrument
    >>> parse_seqid_line('@ERR178930.1 M01234_0214_"000100000-A1474":8:1101:1555:2158#ATCACG length=101')
    FastqSeqTags(platform='illumina', instrument_type='MiSeq', instrument='M01234', run='0214', flow_cell='000100000-A1474', lane='8')

    # GAIIx
    >>> parse_seqid_line("@HWUSI:189:0000-FLOW:8:1101:1219:1953 1:N:0:")
    FastqSeqTags(platform='illumina', instrument_type='GAIIx', instrument='HWUSI', run='189', flow_cell='FLOW', lane='8')

    # Unrecognized instrument type
    >>> parse_seqid_line("@UNKNOWN:189:FLOW0000:8:1101:1219:1953 1:N:0:")
    FastqSeqTags(platform='illumina', instrument_type=None, instrument='UNKNOWN', run='189', flow_cell='FLOW0000', lane='8')

    # MISEQ literal
    >>> parse_seqid_line("@MISEQ:6:000000000-A1445:1:1:16976:1440 2:N:0:CGTACTAGTAGATCGC SEQ000001383")
    FastqSeqTags(platform='illumina', instrument_type=None, instrument='MISEQ', run='6', flow_cell='A1445', lane='1')

    # Hiseq starting with flowcell, no instrument, no run
    >>> parse_seqid_line("@FCC3NWVACXX:3:1101:1161:2200#AACCGAGAA/2")
    FastqSeqTags(platform='illumina', instrument_type='HiSeq', instrument=None, run=None, flow_cell='FCC3NWVACXX', lane='3')

    # Hiseq starting with sample id, then flowcell, no instrument, no run
    >>> parse_seqid_line("@SRR1840614.1 FCC1KPRACXX:1:1101:1291:2172 length=200")
    FastqSeqTags(platform='illumina', instrument_type='HiSeq', instrument=None, run=None, flow_cell='FCC1KPRACXX', lane='1')

    # Lowercase flowcell
    >>> parse_seqid_line("@SRR1166969.1 HWI-ST406:204:d1cywacxx:7:1101:1292:1941 length=100")
    FastqSeqTags(platform='illumina', instrument_type='HiSeq', instrument='HWI-ST406', run='204', flow_cell='d1cywacxx', lane='7')

    References
    ----------
    https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    https://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFormat.htm
    """
    seqid_line = seqid_line.replace('"', '') # strip quotes

    # Does the sequence id line look like one of the Illumina formats?
    for regex in [ILLUMINA_FASTQ_SEQ_ID_REGEX1, ILLUMINA_FASTQ_SEQ_ID_REGEX2]:
        match = regex.search(seqid_line)
        if match:
            instrument = None
            run = None
            flow_cell = match.group(1)
            lane = match.group(2)
            break
    if match is None:
        for regex in [ILLUMINA_FASTQ_SEQ_ID_REGEX3, ILLUMINA_FASTQ_SEQ_ID_REGEX4]:
            match = regex.search(seqid_line)
            if match:
                instrument = match.group(1)
                run = match.group(2)
                flow_cell = match.group(3)
                lane = match.group(4)
                break
    if match is None:
        return None

    # Strip leading zeros and minus from flowcell if present
    flow_cell_parts = flow_cell.split('-')
    if len(flow_cell_parts[0].strip('0')) == 0:  # everything before - is 0
        flow_cell = flow_cell_parts[-1]

    # Lookup the instrument type from the flowcell or instrument name
    instrument_type = flowcell_to_instrument_type(flow_cell)
    if not instrument_type:
        instrument_type = instrument_name_to_instrument_type(instrument)

    return FastqSeqTags("illumina", instrument_type, instrument, run, flow_cell, lane)


def extract_metadata_tags(fastq_path):
    """Examine a fastq file and extract various metadata tags.

    Looks at the first line of the fastq file to extract the instrument,
    flowcell, and lane if possible.

    Parameters
    ----------
    fastq_path : str
        Path to the fastq file.

    Returns
    -------
    FastqSeqTags named-tuple with the following named elements or None if the line cannot be parsed
        platform : str
            Currently, only Illumina data is supported and this will always be "illumina" if the line can be parsed
        instrument_type : str or None
            The instrument type can _sometimes_ be determined from the flowcell or instrument name.
            This is not well tested.
        instrument : str
            Instrument name
        run : str
            Run number
        flow_cell : str
            Flowcell with leading zeros removed
        lane : str
            Lane number

    Examples
    --------
    # Setup tests
    >>> from tempfile import NamedTemporaryFile

    # Uncompressed file
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> num_bytes = f.write("@SRR498276.1 HWI-M00229:9:000000000-A1474:1:1:15012:1874 length=151\\n")
    >>> f.close()
    >>> extract_metadata_tags(filepath)
    FastqSeqTags(platform='illumina', instrument_type='MiSeq', instrument='HWI-M00229', run='9', flow_cell='A1474', lane='1')

    # Compressed file
    >>> os.rename(filepath, filepath + ".gz")
    >>> filepath = filepath + ".gz"
    >>> gf = gzip.open(filepath, "wt")
    >>> num_bytes = gf.write("@SRR498276.1 HWI-ST1029:398:000000000-A1444:2:1:15012:1874 length=151\\n")
    >>> gf.close()
    >>> extract_metadata_tags(filepath)
    FastqSeqTags(platform='illumina', instrument_type='HiSeq', instrument='HWI-ST1029', run='398', flow_cell='A1444', lane='2')

    >>> os.unlink(filepath)

    References
    ----------
    https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    """
    # Case 1: file is gzipped
    if fastq_path.endswith(".gz"):
        line = ""
        f = gzip.open(fastq_path, "rt")
        try:
            line = f.readline()
        finally:
            f.close()
        return parse_seqid_line(line)

    # Case 2: file is NOT gzipped
    with open(fastq_path, "r") as f:
        line = ""
        line = f.readline()
    return parse_seqid_line(line)


# Named tuple to contain read group tags parsed from fastq metadata
ReadGroupTags = collections.namedtuple("ReadGroupTags", "ID SM LB PL PU")

def construct_read_group_tags(fastq_path, sample_name):
    """Examine a fastq file and construct read group tags from the
    flowcell and lane if possible.

    A read group is a set of reads originating from a separate library
    and generated from a single run of a sequencing instrument.  This function
    assumes all the reads in the fastq file are from the same library.

    This function only looks at the first line of the fastq file and assumes
    all the reads in the file have the same read group.  When sequencing is
    spread over multiple lanes (NextSeq), the reads should be kept in separate
    fastq files until after the read groups are identified.

    Parameters
    ----------
    fastq_path : str
        Path to the fastq file.
    sample_name : str
        The name of the sample sequenced in this read group.

    Returns
    -------
    tags : ReadGroupTags or None
        ID : str
            An identifier for the source of reads within the fastq file formed
            from the flowcell and lane.  The read group id is {flowcell}.{lane}".
        SM : str
            This will always be the value of the sample_name argument.
        LB : str
            DNA preparation library identifier.  This function assumes all the reads
            in the fastq file are from the same library.  Always set to "1".
        PL : str
            Platform/technology used to produce the read.  Valid values are:
            illumina, solid, ls454, helicos and pacbio.
        PU : str
            Platform Unit is {flowcell}.{lane}.{sample}.

    Examples
    --------
    # Setup tests
    >>> from tempfile import NamedTemporaryFile

    # Sequence id line contains flowcell and lane
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> num_bytes = f.write("@SRR498276.1 HWI-M00229:9:000000000-A1474:1:1:15012:1874 length=151\\n")
    >>> f.close()
    >>> construct_read_group_tags(filepath, 'sampleA')
    ReadGroupTags(ID='A1474.1', SM='sampleA', LB='1', PL='illumina', PU='A1474.1.sampleA')
    >>> os.unlink(filepath)

    # Sequence id line missing flowcell and lane
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> num_bytes = f.write("@SRR498276.1 length=151\\n")
    >>> f.close()
    >>> construct_read_group_tags(filepath, 'sampleA') is None
    True
    >>> os.unlink(filepath)

    References
    ----------
    http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    """
    tags = extract_metadata_tags(fastq_path)
    if tags is None:
        return None

    id = tags.flow_cell + '.' + tags.lane
    sm = sample_name
    lb = '1'
    pl = tags.platform
    pu = id + '.' + sm
    return ReadGroupTags(id, sm, lb, pl, pu)
