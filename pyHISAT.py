import h5py
import anndata
import numpy as np
import pandas as pd
import plumbum
import pandas as pd



HISAT_BINARY = '/home/michi/hisat2-2.1.0/hisat2'
# HISAT_INDEX = '/home/michi/hisat_indices/grch38/genome'
HISAT_INDEX = '/run/media/michi/618937cd-1ca7-4a88-a2c6-d3e54c1ba4b9/qiagen_pipeline_data/hisat_index/grch38/genome'

def hisat_single(R12_pairs, output_bam, cores=1):
    """
    returns the string for a HISAT call
    """
    # using plumbum (could construct the string here easily, but need plumbum somewhere else too)
    pb_cmd = hisat_single_plumbum(R12_pairs, output_bam, cores)
    return str(pb_cmd)

def hisat_single_plumbum(R12_pairs, output_bam:str, cores=1):
    """
    creates a plumbum command of HISAT
    """
    assert isinstance(output_bam, str)

    hisat = plumbum.local[HISAT_BINARY]
    samtools = plumbum.local['/home/michi/miniconda3_newtry/envs/mkl/bin//samtools']


    r2_files = [r2 for _, r2 in R12_pairs]
    input_args = ",".join(r2_files)

    hisat_bound = hisat['--threads', f'{cores}',
                        '-x', f'{HISAT_INDEX}',
                        '-U', f'{input_args}']

    samtools_convert = samtools['view',
                                '-bS', '-']

    samtools_sort = samtools['sort',
                             '-']

    cmd = hisat_bound | samtools_convert | samtools_sort > output_bam
    return cmd


GTF_FILE = '/home/michi/postdoc_seattle/HL60_SC_vs_bulk/Homo_sapiens.GRCh38.94.gtf.gz'
FEATURECOUNTS_BINARY = '/home/michi/subread-1.6.4-Linux-x86_64/bin/featureCounts'
def featureCounts(bamfile, output_file, cores=1):

    # CMD = f'{FEATURECOUNTS_BINARY}  -a {GTF_FILE}  -o {output_file} {bamfile}'
    pb_cmd = featureCounts_plumbum(bamfile, output_file, cores)
    CMD = str(pb_cmd)
    return CMD

def featureCounts_plumbum(bamfile, output_file, cores=1):
    assert isinstance(bamfile, str)
    PB = plumbum.local[FEATURECOUNTS_BINARY]
    PB_bound = PB['-a', GTF_FILE,
                  '-o', output_file,
                  '-T', cores,
                  bamfile]
    return PB_bound


def load_featureCounts(fc_file):
    df = pd.read_csv(fc_file, sep='\t', skiprows=1)
    df = df.set_index('Geneid')
    df = df.drop(['Chr', 'Start', 'End', 'Strand'],axis=1)
    count_col = [c for c in df.columns if c.endswith('.bam')]
    assert len(count_col) == 1
    count_col= count_col[0]

    df = df.rename({count_col: 'counts'}, axis=1)
    return df
