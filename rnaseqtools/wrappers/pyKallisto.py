"""
Python wrapper around some kallisto calls

also to load kallist odata
"""
import h5py
import anndata
import numpy as np
import pandas as pd


KALLISTO_INDEX = '/home/michi/postdoc_seattle/HL60_SC_vs_bulk/Homo_sapiens.GRCh38.cdna.all.idx.gz'
KALLISTO_INDEX = '/home/michi/Homo_sapiens.GRCh38.cdna.all.kallisto.idx'

GENOMEBAM_GTF  = '/home/michi/postdoc_seattle/HL60_SC_vs_bulk/Homo_sapiens.GRCh38.94.gtf'


def create_kallisto_quant_call(R12_pairs, output_dir, n_bootstraps, cores=6, genomebam=False):
    """
    kallisto qunatification with paired end reads
    """
    assert isinstance(R12_pairs, list)

    input_args = ""  # alternating R1 and R2 pairs (to aggregate multiple lanes)
    for r1, r2 in R12_pairs:
        input_args = input_args+f' {r1} {r2}'

    CMD = f"kallisto quant -b {n_bootstraps} -i {KALLISTO_INDEX} -o {output_dir} -t {cores} {input_args}"
    if genomebam:
        CMD = f'{CMD} --genomebam --gtf {GENOMEBAM_GTF}'

    return CMD


def create_kallisto_pseudo_call(R12_pairs, output_dir, cores=6):
    assert isinstance(R12_pairs, list)

    input_args = ""  # alternating R1 and R2 pairs (to aggregate multiple lanes)
    for r1, r2 in R12_pairs:
        input_args = input_args+f' {r1} {r2}'

    CMD = f"kallisto pseudo -i {KALLISTO_INDEX} -o {output_dir} -t {cores} {input_args}"
    return CMD


def create_kallisto_pseudo_call_10x(R12_pairs, output_dir, cores=6, frag_len_mean=200, frag_len_sd=20):
    """
    this one uses only the R2 file, which contains the biological reads from
    the 10x experiments (R1 only has barcoding data)
    """
    input_args = ""  # alternating R1 and R2 pairs (to aggregate multiple lanes)
    for _, r2file in R12_pairs:
        input_args = input_args+f' {r2file}'

    CMD = f"kallisto pseudo -i {KALLISTO_INDEX} -o {output_dir} -t {cores} {input_args} --single -l {frag_len_mean} -s {frag_len_sd}"
    return CMD


def create_kallisto_quant_call10_bulk(lane_dirs, output_dir, n_bootstraps, cores=6, frag_len_mean=200, frag_len_sd=20, genomebam=False):
    """
    kallisto qunatification with single end reads: R1 is IGNORED!!!!
    """
    input_args = ""  # alternating R1 and R2 pairs (to aggregate multiple lanes)
    for _, r2file in lane_dirs:
        input_args = input_args+f' {r2file}'

    CMD = f"kallisto quant -b {n_bootstraps} -i {KALLISTO_INDEX} -o {output_dir} -t {cores} {input_args} --single -l {frag_len_mean} -s {frag_len_sd}"

    if genomebam:
        CMD = f'{CMD} --genomebam --gtf {GENOMEBAM_GTF}'

    return CMD



def kallist_bootstrap_2_adata(fname):
    """
    turns a bootstrapped kallisto quantification into a scnpy object,
    where each bootstrap is considered a separate datapoint/sample!
    """
    with h5py.File(fname, 'r') as qq:
        n_bootstraps = len(qq['bootstrap'])
        X = np.stack([qq['bootstrap'][f'bs{i}'][:] for i in range(n_bootstraps)])
        obs = pd.DataFrame()
        obs['bootstrap'] = np.arange(n_bootstraps)
        var = pd.DataFrame()
        var['length'] = qq['aux/lengths'][:]
        var['target_id'] = [_.decode() for _ in qq['aux/ids'][:]]
        var['eff_length'] = qq['aux/eff_lengths'][:]
        adata = anndata.AnnData(X, obs=obs,var=var.set_index('target_id'))
    return adata


def kallisto_tsv_2_adata(fname, samplename, units):
    import scanpy as sc
    """
    kallisto quant outputs a single tsv containing a single sample
    but it has 4 rows: `target_id	length	eff_length	est_counts	tpm`
    """
    assert units in ['tpm','est_counts']
    df = pd.read_csv(fname, sep='\t')
    var = df[['target_id', 'length', 'eff_length']].set_index('target_id')

#     units = 'est_counts'
    X = df[['target_id', units]].set_index('target_id').values.T
    obs = pd.DataFrame([{'samplename': samplename}])
    return sc.AnnData(X, obs=obs, var=var)
