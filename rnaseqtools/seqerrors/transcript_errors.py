
"""
# Based on shadows
**important**: For plotting purposes, dont filter the n_shadow==0 observations. This will lead to bias:
- assume two datasets with the exact same (small) error rate
- dataset A has alot more datapoints then dataset B
- since the error rate is small A will have very discrete and coarse counter, B will look more continuos
- filtering out the 0 will affect A more. Similar datapoint in B will have very small errorrate (which gets rounded to 0 in A) but will persist here!

**NOTE** as similar thing happens when taking the log. `plotnine` tends to filter these observation before calculating boxplots eg
"""
import pysam
import collections
import pandas as pd
import numpy as np
import pysam
import seaborn as sns
import collections
import tqdm
import matplotlib.pyplot as plt
import plotnine as pn
import statsmodels.api as sm
from statsmodels.formula.api import ols
import warnings


def _get_1BP_mutants(seq:str, position:int):
    "return the three 1BP mutations of a sequence at the given position"
    for base in ['A', 'C', 'G', 'T']:
        new = seq[:position] + base + seq[position+1:]
        if seq != new:
            yield new


def _get_1BP_2BP_mutants(seq):
    for i in range(len(seq)):
        for bp1_mut in _get_1BP_mutants(seq, i):
            yield bp1_mut
            for j in range(i+1, len(seq)):
                for bp2_mut in _get_1BP_mutants(bp1_mut, j):
                    if bp2_mut != seq:
                        yield bp2_mut

"""
from scipy.special import binom
s = 'AAAAAAAAAAAA'
N = len(s)
assert list(_get_1BP_2BP_mutants(s)) == 3 * 3 * binom(N,2) + 3 * binom(N,1)
"""


def _get_number_of_shadows(cb, counter_some_error):
    """
    for a given sequence, check the total reads that its 1BP neighbours have
    we do this position dependent too: how many reads where CB and read differ by 1BP in the first base

    return a dict with 16 positions, with the number of 1BP neighbours at each position
    """
    shadow_counts_per_position = {}
    for i in range(len(cb)):
        shadow_counts_per_position[i] = 0
        for mutant in _get_1BP_mutants(cb, position=i):
            if mutant in counter_some_error:
                shadow_counts_per_position[i] += counter_some_error[mutant]
    return shadow_counts_per_position


def estimate_error_rate_shadows(counter_no_error, counter_some_error):
    """
    estimate the error rate based on cell barcodes (the whitelist makes it easy to distinguish correct reads from errors)
    using the "shadows"

    :param counter_no_error: prepared from `create_error_per_cb_counters()`
    """
    most_common_CBS = get_most_common_true_sequences(counter_no_error, topN=1000)
    print(f'most common {len(most_common_CBS)}/1000')

    df = []
    df_pos = []
    N = 16  # bases in the barcode

    for cb in tqdm.tqdm(most_common_CBS):
        freq = counter_no_error[cb]
        assert len(cb)==N

        shadow_counts_per_position = _get_number_of_shadows(cb, counter_some_error)
        n_shadows_total = sum(shadow_counts_per_position.values())
        df.append({'n_shadows': n_shadows_total, 'n_real': freq, 'cb': cb})

        for pos, n_shad in shadow_counts_per_position.items():
            df_pos.append({'n_shadows': n_shad, 'n_real': freq, 'cb': cb, 'position': pos})
    df = pd.DataFrame(df)
    df_pos = pd.DataFrame(df_pos)

    df['error_rate'] = df['n_shadows'] / (df['n_real'] + df['n_shadows'] )
    df_pos['error_rate'] = df_pos['n_shadows'] / (df_pos['n_real'] + df_pos['n_shadows'] )

    return df, df_pos


def read_filter(read):
    """
    filtering the R2 for uniquely mapped, no gaps, no N
    """
    if read.is_duplicate:
        return False

    if read.is_secondary or read.is_supplementary or read.is_unmapped:
        return False

    if not read.has_tag('CR') or not read.has_tag('UR'):  # UR and CR are the RAW uncorrected sequences
        return False

    length = len(read.query_sequence)
    # all bases have to be aligned without gaps
    if read.cigarstring != f'{length}M':
        return False
    seq = read.query_sequence

    if 'N' in seq:
        return False

    mismatches = read.get_tag('nM')
    # just to make sure that this read REALLY maps here
    if mismatches > 1:
        return False

    # only uniquely mapped reads (10x uses STARalign and thats 255 for unique mapped)
    if read.mapping_quality != 255:
        return False
    return True


def estimate_shadow_3prime(bamfile, region):
    """
    counts occurances of all reads in the given region.
    also takes care of reverse reads (flipping them)
    """
    read_counter = collections.defaultdict(int)

    for read in tqdm.tqdm(pysam.AlignmentFile(bamfile).fetch(region=region)):
        if read_filter(read) is False:
            continue

        seq = read.query_sequence
        if read.is_reverse:
            seq = seq[::-1]
        read_counter[seq] += 1

    return read_counter


def get_most_common_true_sequences(read_counter, topN:int):
    """
    get the most abundant sequences, but also make sure that shadows dont sneak in.
    e.g. a VERY abundant true sequence might be ~100000reads, and 1% (1000)
    will result in shadows. these shadows might end up in the top100 itself
    """
    most_common = set()
    flagged_shadows = set()  # all sequences that we've seen as a shadow of a more abundant molecule
    for seq, freq in tqdm.tqdm(collections.Counter(read_counter).most_common(topN), desc='finding most common seqs'):

        # if its already flagged as shadow, skip it
        if seq in flagged_shadows:
            continue

        most_common.add(seq)

        # flag all 1BP or 2BP mutants
        for mutant in _get_1BP_2BP_mutants(seq):
            flagged_shadows.add(mutant)

        # # flag 1BP mutatnts
        # for i in range(len(seq)):
        #     for mutant in _get_1BP_mutants(seq, position=i):
        #         flagged_shadows.add(mutant)
    return most_common


def read_counter_to_df(read_counter, topN):
    """
    turns a read counter (which has real and shadow molecules) into
    a dataframe, each row being a true barcode (based on frequency, topN
    barcodes will be used),  and the number of real reads
    and shadows annotated
    """
    df_seq = []
    most_common = get_most_common_true_sequences(read_counter, topN)
    print(len(most_common))
    for seq in most_common:
        shadow_per_position = _get_number_of_shadows(seq, read_counter)
        for pos, n_shadow in shadow_per_position.items():
            df_seq.append({'seq': seq, 'n_shadow': n_shadow,
                           'n_real': read_counter[seq], 'position': pos})

    df_seq = pd.DataFrame(df_seq)
    return df_seq


def read_counter_to_subsitution_table(read_counter, topN):
    """
    estimate real and shadow molecules, and keep track of which mutations/errors occured!
    """
    df_substitutions = []
    most_common = get_most_common_true_sequences(read_counter, topN)
    print(len(most_common))
    for seq in most_common:
        for pos in range(len(seq)):  # check all mutation at base-position

            mut_frequencies = {'A': 0, 'C': 0, 'G': 0, 'T':0}
            for mutant in _get_1BP_mutants(seq, pos):
                if mutant in read_counter:  # we found a existing one mutant
                    base_mut = mutant[pos]
                    mut_frequencies[base_mut] += read_counter[mutant]

            base_true = seq[pos]
            freq_mut = sum(mut_frequencies.values())
            mut_frequencies[base_true] = read_counter[seq]  # record the number of correct basepairs too!
            r = {'base_true': base_true,
                 'position': pos, 'frequency_true': read_counter[seq], 'frequency_mut': freq_mut,
                 'frequency_A': mut_frequencies['A'],
                 'frequency_C': mut_frequencies['C'],
                 'frequency_G': mut_frequencies['G'],
                 'frequency_T': mut_frequencies['T'],
                 'seq': seq}
            df_substitutions.append(r)
    df_substitutions = pd.DataFrame(df_substitutions)
    return df_substitutions


def estimate_error_rate(df):
    """
    estimate the error rate based on real and shadow counts.
    Dataframe is created via `read_counter_to_df()`
    """
    assert 'n_shadow' in df.columns and 'n_real' in df.columns
    # linear regression
    model = ols('n_shadow ~ -1 + n_real', df)
    res = model.fit()
    beta = res.params['n_real']
    confidence_interval = res.conf_int().loc['n_real'].values
    error_linear = beta / (1 + beta)
    error_c5 = confidence_interval[0] / (1+confidence_interval[0])
    error_c95 = confidence_interval[1] / (1+confidence_interval[1])
    # simple ratio
    naive = np.mean(df.n_shadow / (df.n_shadow+df.n_real))

    # binomail "regression": actually just estimating the paramter of a binomial,
    if np.all(df['n_shadow'] == 0):
        theta_binom = 0
        theta_binom_c5 = 0
        theta_binom_c95 = 0
    else:
        y = df[['n_shadow', 'n_real']].values
        r = np.ones(len(y))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            glm = sm.GLM(y, r, family=sm.families.Binomial(sm.families.links.log()))
            res = glm.fit()
            assert len(res.params) == 1
            theta_binom = np.exp(res.params[0])
            theta_binom_c5, theta_binom_c95 = np.exp(res.conf_int()[0])

    # robust linear reg
    if np.all(df['n_shadow'] == 0):
        robust_err = 0
    else:
        rmod = sm.RLM(df['n_shadow'], df['n_real'])
        rres = rmod.fit()
        robust_beta = rres.params[0]
        robust_err = robust_beta / (1 + robust_beta)
    return {'error_linear': error_linear, 'error_linear_c5': error_c5, 'error_linear_c95': error_c95,
            'error_naive': naive,
            'error_binom': theta_binom, 'error_binom_c5': theta_binom_c5, 'theta_binom_c95': theta_binom_c95,
            'error_robust': robust_err}
