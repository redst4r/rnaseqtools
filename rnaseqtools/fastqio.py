import gzip
import itertools
"""
here's a fast iterator over the sequences, doesnt do the checking done by skbio
"""


def read_fastq_seqs_bare(filepath):
    """returns the bare 4 lines of the fastq entry, no postprocesing (i.e. newline stripping)"""
    with gzip.open(filepath, 'rt') as fh:
        for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
            if any([line is None for line in (seq_header, seq, qual_header, qual)]):
                raise Exception(
                    "Number of lines in FASTQ file must be multiple of four "
                    "(i.e., each record must be exactly four lines long).")

            if not seq_header.startswith('@'):
                raise Exception("Invalid FASTQ sequence header: %r" % seq_header)
            if qual_header != '+\n':
                raise Exception("Invalid FASTQ quality header: %r" % qual_header)
            if qual == '\n':
                if not seq =='\n':  # sometimes theres just no sequence, hence also quality socre misses
                    raise Exception("FASTQ record is missing quality scores.")
            yield seq_header, seq, qual_header, qual


def read_fastq_seqs(filepath):
    """
    reading the fastq entries, stripping newlines, returning only read_id and read_sequence
    """
    for seq_header, seq, qual_header, qual in read_fastq_seqs_bare(filepath):
        read_id = seq_header[1:].split(" ")[0]
        yield read_id, seq.rstrip('\n')


def read_fastq_seqs_multiple_lanes(fastq_files: list):
    """
    seemlessly concatenates the reads from multiple fastq reads;

    fastq_files: a list of fastq filenames that will all be concatenated

    Return:
        generator of:
        Read_Id, Sequence
    """
    return itertools.chain.from_iterable(fastq_files)
