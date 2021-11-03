import collections
import pybktree
import toolz
import pybustools.pybustools as pb
from pybustools.butterfly import make_ec2gene_dict
import tqdm
import pandas as pd

"""
estimating the sequencing error rates in the UMIs

given as kallisto-busfile, check each cell for UMIs that map to the same gene
and are only 1BP appart. These are most likely sequencing errors
"""


def hamming_distance(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguements '''

    dist = sum([not a == b for a, b in zip(first, second)])
    return dist

def groupby_gene_and_umi(bus_records, ec2gene_dict):
    """
    group a set of record by 1) gene 2) UMI and record the mulitplicity (r.COUNT)

    returns a dict[gene -> dict[ umi->count ]]

    """
#     gene2umilist = collections.defaultdict(set)
    gene2umilist = collections.defaultdict(collections.Counter)

    assert [r.CB == bus_records[0].CB for r in bus_records], "must share same CB"

    for r in bus_records:
        g = ec2gene_dict[r.EC]
        if len(g) == 1:
            # g = g.pop()  # THIS IS BAD, modifies the ec2gene_dict
            g = list(g)[0]
        else:
            continue # multimapped

#         gene2umilist[g].add(r.UMI)
        gene2umilist[g][r.UMI] += r.COUNT

    return dict(gene2umilist)


def group_UMIs_into_cliques(umis_from_gene, cb,gene, topN=100):
    """
    greedy assignment of UMIs into cliques.

    To find out which UMIs are all derived from the same molecule,
    we look at the x most abundant UMIs (in terms of reads).
    any UMI that is 1BP away and less abundant will be assigned to taht
    clique (and removed from the pool).
    """

    already_processed = set()
    cliques = []
    max_errors = 1  # maximum erros to be considered in grouing the UMIs
    ptree = pybktree.BKTree(hamming_distance, items=umis_from_gene.keys())

    for umi, count in umis_from_gene.most_common(topN):
        if umi in already_processed:
            continue

        already_processed.add(umi)

        new_clique = collections.defaultdict(int)
        # find all sequences that are closeby
        for dist, umi2 in ptree.find(umi, n=max_errors):  # note: this also pull in the d=0 (i.e. original), which is good!
            umi_freq = umis_from_gene[umi2]
            new_clique[dist] += umi_freq
            already_processed.add(umi2)

        n = {f'errors_{i}': new_clique[i] for i in range(max_errors+1)}
        n['UMI'] = umi
        n['gene'] = gene
        n['cb'] = cb
        cliques.append(n)
    return cliques


def UMI_errors(busfolder, t2g, max_entries=1_000_000):
    """
    iterate over busfile by cell, groupby gene/UMI and check for UMIs that are 1BP appart
    """
    B = pb.Bus(busfolder)
    ec2gene_dict = make_ec2gene_dict(B, t2g)
    res = []
    counter = 0

    # perpare the stream
    I = B.iterate_cells()
    I = toolz.filter(lambda cb_recordlist: len(cb_recordlist[1]) > 100, I)
    I = toolz.random_sample(0.1, I)

    for cb, bus_records in tqdm.tqdm(I):

        if len(bus_records) < 100:
            assert 1 == 0, "should never happens, filtered above"
            continue

        counter += 1
        s = groupby_gene_and_umi(bus_records, ec2gene_dict)

        for gene, umi_set in s.items():
    ## this is too expensive, takes 99% CPU
    #         t = pd.DataFrame(group_UMIs_into_cliques(umi_set))
    #         t['gene'] = gene
    #         t['CB'] = cb

            t = group_UMIs_into_cliques(umi_set, cb, gene)
            res.extend(t)

        if counter % 100 == 0:
            print(counter, len(res))

        if len(res) > max_entries:
            break
    res = pd.DataFrame(res)

    return res
