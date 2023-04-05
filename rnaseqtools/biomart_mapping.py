import warnings
import pathlib
import io
import os
import collections
import tqdm
import numpy as np
import pandas as pd
from bioservices import biomart

HOST = 'uswest.ensembl.org'

"""
helper functions to create the biomart dataframe
"""

def _biomart_df_postprocess(df):
    df['chromosome_name'] = df['chromosome_name'].apply(lambda x: str(x))
    return df


def print_attributes():
    """
    thats all the atribures that we can get from biomart for a gene/transcript
    """
    s = biomart.BioMart(host=HOST)
    return list(s.attributes('hsapiens_gene_ensembl').keys())


def biomart_query_transcripts(transcript_ids, batchsize=10000, verbose=False):
    return biomart_query(transcript_ids, 'ensembl_transcript_id', batchsize, verbose)


def biomart_query_genes(gene_ids, batchsize=10000, verbose=False):
    return biomart_query(gene_ids, 'ensembl_gene_id', batchsize, verbose)


import hashlib
def md5(fname: str):
    """
    calcualte the m5sum of a file
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def biomart_query_all(verbose=False, extra_fields=None, force_download=False):
    """
    pulls down all entries from BIOMART for Human: symbol, trasncript, gene, length, type
    """

    THE_FILE = pathlib.Path(__file__).parent / 'biomart_all.csv.gz'

    assert  md5(str(THE_FILE)) == "34f414e3384c0e19294de6be6c947090", "biomart file changed on disk!"
    if not force_download and os.path.exists(THE_FILE):
        return _biomart_df_postprocess(pd.read_csv(THE_FILE, index_col=0))

    raise ValueError("THE BIOMART FILE SHOULD BE INCLUDED ALREADY!!")

    s = biomart.BioMart(host=HOST)
    s.new_query()
    s.add_dataset_to_xml('hsapiens_gene_ensembl')

    # what we want to get back
    # s.add_attribute_to_xml('entrezgene')

    fields = [
        'hgnc_symbol',
        'ensembl_gene_id',
        'ensembl_gene_id_version',
        'transcript_length',
        'ensembl_transcript_id',
        'ensembl_transcript_id_version',
        'transcript_biotype',
        'chromosome_name',
        'start_position',
        'end_position',
        'external_synonym',
    ]

    if extra_fields:
        fields.extend(extra_fields)

    for f in fields:
        s.add_attribute_to_xml(f)

    xml = s.get_xml()

    if verbose:
        print(xml)

    res = s.query(xml)

    df = pd.read_csv(io.StringIO(res), sep='\t', header=None)
    df.columns = fields
    df = df.drop_duplicates()

    df = _biomart_df_postprocess(df)

    df.to_csv(THE_FILE)

    return df


def batch(some_iterable, batchsize):
    """
    splits the iterable into a couple of chunks of size n
    handy for iterating over batches
    :param some_iterable:  iterable to be chunked/batched
    :param batchsize: batchSize
    :return: gnerator over iterables
    """
    assert isinstance(some_iterable, collections.Iterable)  # TODO this does not guard against np.arrays as they are also iterable (over single elements)
    l = len(some_iterable)
    for ndx in range(0, l, batchsize):
        yield some_iterable[ndx:min(ndx + batchsize, l)]


def biomart_query_new(attributes:list, filter:str, filter_values:list, batchsize=10000, verbose=False):
    """
    query biomart similar to the R-function getBM()

    :param attributes: which attributes to retrieve from biomart
    :param filter: field which to filter on
    :param filter_values: return only entries where filter \in filter_values
    """
    s = biomart.BioMart(host=HOST)

    batch_results = []
    n_batches = int(np.ceil(len(filter_values)/ batchsize))
    for id_batch in tqdm.tqdm(batch(filter_values, batchsize=batchsize), total=n_batches):

        s.new_query()
        s.add_dataset_to_xml('hsapiens_gene_ensembl')

        # what we want to get back
        for a in attributes:
            s.add_attribute_to_xml(a)

        # the query should be comma separated
        # better make sure theres no whitespace
        query = ",".join([_.strip() for _ in id_batch])

        s.add_filter_to_xml(filter, query)
        xml = s.get_xml()

        if verbose:
            print(xml)

        res = s.query(xml)

        df = pd.read_csv(io.StringIO(res), sep='\t', header=None)
        df.columns = attributes
        df = df.drop_duplicates()
        batch_results.append(df)

    df = pd.concat(batch_results, axis=0).drop_duplicates()

    return df


def biomart_query(id_list:list, id_type:str, batchsize=10000, verbose=False):
    """
    THIS IS SOME LEGACY FUNTION AS OF 2021/01/14. JUST CALLS biomart_query_new
    which has the proper naming of function arguments (as R's getBM())

    pulls down the symbol/entrez/transcriptlength for each ensebml-rtanscript_id given
    """

    warnings.warn('This function is deprecated. Use biomart_query_new', DeprecationWarning)
    assert id_type in ['ensembl_gene_id', 'ensembl_transcript_id']

    attributes = [
        'hgnc_symbol',
        'ensembl_gene_id',
        'transcript_length',
        'ensembl_transcript_id',
        'transcript_biotype',
        'entrezgene',
    ]
    return biomart_query_new(attributes, filter=id_type, filter_values=id_list, batchsize=batchsize)


import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import pandas as pd

def probeids_to_symbols(probeids):
    """
    convert microarray probe ids into gene names
    """
    importr("biomaRt")
    ro.r(
        """
        require("biomaRt")
        mart <- useMart("ENSEMBL_MART_ENSEMBL")
        mart <- useDataset("hsapiens_gene_ensembl", mart)
        """
    )
    r_probe_ids = ro.StrVector(probeids)

    # move into R
    with localconverter(ro.default_converter):
        ro.globalenv['probe_ids'] = r_probe_ids
        ro.r(
        """
        annotLookup <- getBM(
          mart=mart,
          attributes=c(
            "affy_hg_u133_plus_2",
            "ensembl_gene_id",
            "gene_biotype",
            "external_gene_name"),
          filter = "affy_hg_u133_plus_2",
          values = probe_ids, uniqueRows=TRUE)
            """
        )
    with localconverter(ro.default_converter + pandas2ri.converter):
        df = ro.r('annotLookup')
    return df
