from bioservices import biomart
import io
import collections
import tqdm
import os
import pandas as pd
# datasets = s.datasets("ENSEMBL_MART_ENSEMBL")

# filters = s.filters('hsapiens_gene_ensembl')
# filters['']

# s.attributes('hsapiens_gene_ensembl')




def print_attribtues():
    """
    thats all the atribures that we can get from biomart for a gene/transcript
    """
    s = biomart.BioMart(host='uswest.ensembl.org')
    list(s.attributes('hsapiens_gene_ensembl').keys())


def biomart_query_transcripts(transcript_ids, batchsize=10000, verbose=False):
    return biomart_query(transcript_ids, 'ensembl_transcript_id', batchsize, verbose)


def biomart_query_genes(gene_ids, batchsize=10000, verbose=False):
    return biomart_query(gene_ids, 'ensembl_gene_id', batchsize, verbose)


def biomart_query_all(verbose=False):
    """
    pulls down all entries from BIOMART for Human: symbol, trasncript, gene, length, type
    """

    s = biomart.BioMart(host='uswest.ensembl.org')
    s.new_query()
    s.add_dataset_to_xml('hsapiens_gene_ensembl')

    # what we want to get back
    # s.add_attribute_to_xml('entrezgene')
    s.add_attribute_to_xml('hgnc_symbol')
    s.add_attribute_to_xml('ensembl_gene_id')
    s.add_attribute_to_xml('transcript_length')
    s.add_attribute_to_xml('ensembl_transcript_id')
    s.add_attribute_to_xml('transcript_biotype')
    s.add_attribute_to_xml('entrezgene')


    xml = s.get_xml()

    if verbose:
        print(xml)

    res = s.query(xml)

    df = pd.read_csv(io.StringIO(res), sep='\t', header=None)
    df.columns=['hgnc_symbol', 'ensembl_gene_id', 'transcript_length', 'ensembl_transcript_id', 'transcript_biotype', 'entrezgene']
    df = df.drop_duplicates()
    return df

def biomart_query(id_list:list, id_type:str, batchsize=10000, verbose=False):
    """
    pulls down the symbol/entrez/transcriptlength for each ensebml-rtanscript_id given
    """

    assert id_type in ['ensembl_gene_id', 'ensembl_transcript_id']

    s = biomart.BioMart(host='uswest.ensembl.org')

    "we have to batch the call if the list of queries is long"
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


    batch_results = []
    n_batches = int(np.ceil(len(id_list)/ batchsize))
    for id_batch in tqdm.tqdm(batch(id_list, batchsize=batchsize), total=n_batches):

        s.new_query()
        s.add_dataset_to_xml('hsapiens_gene_ensembl')

        # what we want to get back
        # s.add_attribute_to_xml('entrezgene')
        s.add_attribute_to_xml('hgnc_symbol')
        s.add_attribute_to_xml('ensembl_gene_id')
        s.add_attribute_to_xml('transcript_length')
        s.add_attribute_to_xml('ensembl_transcript_id')
        s.add_attribute_to_xml('transcript_biotype')
        s.add_attribute_to_xml('entrezgene')

        # the query should be comma separated
        # better make sure theres no whitespace
        query = ",".join([_.strip() for _ in id_batch])

        s.add_filter_to_xml(id_type, query)
        xml = s.get_xml()

        if verbose:
            print(xml)

        res = s.query(xml)

        df = pd.read_csv(io.StringIO(res), sep='\t', header=None)
        df.columns=['hgnc_symbol', 'ensembl_gene_id', 'transcript_length', 'ensembl_transcript_id', 'transcript_biotype', 'entrezgene']
        df = df.drop_duplicates()
        batch_results.append(df)

    df = pd.concat(batch_results, axis=0).drop_duplicates()

    return df
