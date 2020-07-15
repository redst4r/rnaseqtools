import pandas as pd
import numpy as np
import os
import pysam
import tqdm
import pybktree
import collections
import pandas as pd
import toolz

"""
now the shadow regresion thing:
- we know the "true" CBs from the whitelist
- for each true CB: look at the number of correct reads vs incorrect:
  these are CBs that have one subsitution
- it will take forever to do it on all CBs, but take the top1000 frequent ones
We later might use the UMIs t figure out whats PCR error or seq
"""


def _load_whitelist(fname):
    "loads a whitelist of cellbarcdes as a set"
    with open(fname, 'r') as fh:
        whitelist = set(_.strip() for _ in fh.readlines())
    return whitelist


def hamming_distance(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguements '''

    dist = sum([not a == b for a, b in zip(first, second)])
    return dist


"""
to speed thigns up:
- instead of gathering all CBs in the entire file, just look at the whitelisted ones
- build a bktree on the whitelist
- for an observed CB: check if its in the whitelisted
    - if yes: increase error_free coutn
    - if no: check the bktree for similar cells: if theres a single one: increase the error count for that one
"""


def create_error_per_cb_counters(bamfile, whitelist):
    """
    for each barcode, check if its in the whitelist, and count the occurances
    """
    bam = pysam.AlignmentFile(bamfile)
    counter_no_error = collections.defaultdict(int)
    counter_some_error = collections.defaultdict(int)

    # counter = 0
    # counter_max = 100_000_000
    for read in tqdm.tqdm(bam.fetch()):
        if read.has_tag('CR') and read.has_tag('UR'):  # UR and CR are the RAW uncorrected sequences
            CB = read.get_tag('CR')
            # UMI = read.get_tag('UR')

            if CB in whitelist:
                counter_no_error[CB] += 1
            else:
                counter_some_error[CB] += 1
    return counter_no_error, counter_some_error


if __name__ == '__main__':

    whitelist = _load_whitelist('/home/mstrasse/resources/3M-february-2018.txt')
    # bktree_whitelist = pybktree.BKTree(hamming_distance, tqdm.tqdm(whitelist))

    samples = [
        'novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs',
        'dnbseqg400.V300039753.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
        'dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
        'dnbseqg400.V300026370_88A.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
        'dnbseqg400.V300039753.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
        'novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs',
        'dnbseqg400.V300035342.L06B_2-718732_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
        'dnbseqg400.V300035342.L06C_2-718733_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
        'dnbseqg400.V300035355.L06C_2-718733_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs'
    ]

    res_dict = {}
    for s in samples:
        bamfile = f'/home/mstrasse/mount/download/{s}/possorted_genome_bam.bam'

        counter_no_error, counter_some_error = create_error_per_cb_counters(bamfile, whitelist)
        res_dict[s] = (counter_no_error, counter_some_error)


    import pickle
    with open('/home/mstrasse/seqerror.pkl', 'wb') as fh:
        pickle.dump(res_dict ,fh)

    for s in samples:
        counter_no_error, counter_some_error = res_dict[s]
        U1 = sum(counter_some_error.values())
        U0 = sum(counter_no_error.values())

        print(f'{s}: Error rate: {U1 / (U1 + U0):.4f}')

    import pickle
    with open('/home/mstrasse/seqerror.pkl', 'rb') as fh:
        res_dict = pickle.load(fh)



    """
    now the shadow regresion thing:
    - for each true CB: look at the number of correct reads vs incorrect: these are CBs that have one subsitution
    - it will take forever to do it on all CBs, but take the top1000 frequent ones
    """


    def estimate_error_rate_shadows(counter_no_error, counter_some_error):
        most_common_CBS = collections.Counter(counter_no_error).most_common(1000)
        shadows = {}
        real = {}
        df = []
        for cb_freq in tqdm.tqdm(most_common_CBS):
            cb, freq = cb_freq

            shadow_counts  = []
            for i in range(16):
                for base in ['A','C','G','T']:
                    new = cb[:i] + base + cb[i+1:]
                    if cb == new:
                        continue
                    if new in counter_some_error:
                        shadow_counts.append(counter_some_error[new])
            shadows[cb] = sum(shadow_counts)
            real[cb] = freq
            df.append({'n_shadows': shadows[cb], 'n_real': freq, 'cb': cb})
        df = pd.DataFrame(df)
        df['error_rate'] = df['n_shadows'] / df['n_real']
        # error_rates = [shadows[cb] / real[cb] for cb, freq in most_common_CBS]

        return df


    dfs = toolz.valmap(lambda x: estimate_error_rate_shadows(x[0], x[1]), res_dict)


    """
    the above its the per read error rate:  #correct reads / #wrong reads

    But how does that translate to a per-base error:
    - to get 16BP correct, we need (1-p_err)^16
    - a single error is 16 * err * (1-err)^15

    %correct = (1-p_err)^16

    c**(1/16) = 1-p
    p = 1-  c**(1/16)

    """
    right = U0 / (U1 + U0)
    1- right**(1/16)
