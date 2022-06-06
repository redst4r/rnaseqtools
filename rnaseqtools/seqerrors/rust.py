from rnaseqtools.seqerrors.transcript_errors import estimate_error_rate
import pandas as pd
import os

PREFIX = 'cd /home/mstrasse/TB4/rust_code/rust_shadow; RUSTUP_HOME=/home/mstrasse/TB4/rust_installation/.rustup CARGO_HOME=/home/mstrasse/TB4/rust_installation/.cargo'

def mk_rustcall_cb(fastq_glob, whitelist_file, outfile, topn):
    s1 = f"{PREFIX} cargo run --release -- --output {outfile} cb -w {whitelist_file} --ntop {topn}  {fastq_glob}"
    # print(s1)
    os.system(s1)

def mk_rustcall_cbumi(fastq_glob, whitelist_file, outfile, topn):
    s1 = f"{PREFIX} cargo run --release --  --output {outfile} cb-umi-sketch -w {whitelist_file} --ntop {topn}  {fastq_glob}"
    # print(s1)
    os.system(s1)


def mk_rustcall_cbumi_full(fastq_glob, whitelist_file, outfile, topn):
    s1 = f"{PREFIX} cargo run --release -- --output {outfile} cb-umi-exact -w {whitelist_file} --ntop {topn}   {fastq_glob}"
    # print(s1)
    os.system(s1)

def mk_rustcall_cbumi_cell(busfile, outfile, nmax, aggr: bool):
    cmd = 'cb-umi-cell-aggr' if aggr else 'cb-umi-cell'
    s1 = f"{PREFIX} cargo run --release -- --output {outfile} {cmd} --nmax {nmax} {busfile}"
    # print(s1)
    os.system(s1)

def rust_output_to_error_estimate(df_rust):
    n_bases = len([_ for _ in df_rust.columns if _.startswith('position_')])
    df_error2 = []
    for i in range(n_bases):
        col = f'position_{i}'
        _df = df_rust[['n_real', col]].rename({col:'n_shadow'}, axis=1)
        s = estimate_error_rate(_df)
        s['position'] = i
        df_error2.append(s)
    df_error2 = pd.DataFrame(df_error2)
    return df_error2


def rust_read_cb_results(samplename, DIRECTORY):
    # reading the precompuated error estimates in the directory
    results_df = pd.read_csv(f'{DIRECTORY}/{samplename}/cb.csv')
    df_error  = rust_output_to_error_estimate(results_df)
    df_error['samplename'] = samplename

def rust_read_umi_results(samplename, DIRECTORY):
    # reading the precompuated error estimates in the directory
    results_df = pd.read_csv(f'{DIRECTORY}/{samplename}/cb.csv')
    df_error  = rust_output_to_error_estimate(results_df)
    df_error['samplename'] = samplename


def load_rust_cb_cell_aggregate(samplename, remove_singletons, DIRECTORY):
    """
    loading the error estimates for UMIs; based on the aggregated error estimates (rustfastq cb-umi-cell-aggr)
    """
    df_raw = pd.read_csv(f'{DIRECTORY}/{samplename}/umi.csv')
    df_raw = df_raw.rename({f'position_{i}_sum': f'position_{i}' for i in range(28)}, axis=1)
    df_raw = df_raw.rename({'n_real_sum': 'n_real'}, axis=1)

    if remove_singletons:
        df_raw['n_real'] = df_raw['n_real'] - df_raw['singeltons']
        df_raw['n_total'] = df_raw['n_total'] - df_raw['singeltons']

    df_error = beta_binomial_error_estimates(df_raw)
    df_error['samplename'] = samplename
    
    return df_error

from scipy.stats import beta
def beta_binomial_error_estimates(df_raw):
    """
    given the table of shadow and real molecules, estimate the
    sequencing errors using a beta-binomial model
    """
    pos_cols = [_ for _ in df_raw.columns if _.startswith('position')]
    df_raw['n_shadow'] = df_raw[pos_cols].sum(1)
    df_raw['n_total'] = df_raw['n_real'] + df_raw['n_shadow']

    df_error = []
    for i in range(len(pos_cols)):
        col = f'position_{i}'
        pseudo_count = 1
        rv = beta(a=pseudo_count+df_raw[col].sum(), b=pseudo_count+df_raw['n_real'].sum())
        df_error.append({'error_binom': rv.mean(), 'error_binom_c95': rv.ppf(0.95), 'error_binom_c5': rv.ppf(0.05), 'position': i})
    df_error = pd.DataFrame(df_error)

    return df_error
