from rnaseqtools.seqerrors.transcript_errors import estimate_error_rate
import pandas as pd
import os
def mk_rustcall_cb(fastq_glob, whitelist_file, outfile, topn):
    s1 = f"cd /home/mstrasse/TB4/rust_code/rust_shadow; RUSTUP_HOME=/home/mstrasse/TB4/rust_installation/.rustup CARGO_HOME=/home/mstrasse/TB4/rust_installation/.cargo cargo run --release -- -w {whitelist_file} --ntop {topn} --output {outfile} --command cb {fastq_glob}"
    print(s1)
    os.system(s1)

def mk_rustcall_cbumi(fastq_glob, whitelist_file, outfile, topn):
    s1 = f"cd /home/mstrasse/TB4/rust_code/rust_shadow; RUSTUP_HOME=/home/mstrasse/TB4/rust_installation/.rustup CARGO_HOME=/home/mstrasse/TB4/rust_installation/.cargo cargo run --release -- -w {whitelist_file} --ntop {topn} --output {outfile} --command cb_umi_sketch {fastq_glob}"
    print(s1)
    os.system(s1)


def mk_rustcall_cbumi_full(fastq_glob, whitelist_file, outfile, topn):
    s1 = f"cd /home/mstrasse/TB4/rust_code/rust_shadow; RUSTUP_HOME=/home/mstrasse/TB4/rust_installation/.rustup CARGO_HOME=/home/mstrasse/TB4/rust_installation/.cargo cargo run --release -- -w {whitelist_file} --ntop {topn} --output {outfile} --command cb_umi_exact {fastq_glob}"
    print(s1)
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
