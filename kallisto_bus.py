import os
from pathlib import Path
import gzip
import tqdm
import hashlib
import numpy as np
import pandas as pd

from kallisto_em import KallistoEM
"""
from different sources on the kallisto github page

- https://github.com/paulranum11/kallisto_pseudo_to_expressionMatrix/blob/master/prep_TCC_matrix.py

"""

class KallistoPseudoParser():

    def __init__(self, transcriptome_fasta):
        self.transcriptome_fasta = transcriptome_fasta  # the transript IDs are just numbers (order as in the fasta file used to create the index)
        self.ecID_2_transcript = None
        self.transcript_length = None

    def create_ecID_2_transcript(self, ):
        """
        gets the mapping (int) -> TranscriptName  for the fasta file used to geenrate the kallisto index
        which we need to associate EC with the respective trasncrtips (the transcripts are encoded as int)

        also extract hte length of each transcrtipy from the fasta

        see https://github.com/paulranum11/kallisto_pseudo_to_expressionMatrix/blob/master/prep_TCC_matrix.py
        """
        print("Loading index fasta..")
        print("Extracting transcript IDs..")
        transcriptID_Dict = {}
        # with open(self.transcriptome_fasta, "r") as infile:
        with gzip.open(self.transcriptome_fasta, "r") as infile:

            id_ct = -1

            transcript_length_dict = {}
            for line in tqdm.tqdm(infile):
                line = line.decode()
                if (">" in line):
                    split1 = line.split(" ")
                    split2 = split1[0].split(">")
                    transcriptID = split2[1]
                    id_ct = id_ct + 1
                    #print(str(id_ct) + " " + transcriptID)
                    transcriptID_Dict[id_ct] = transcriptID
                    transcript_length_dict[id_ct] = 0
                else:
                    transcript_length_dict[id_ct] += len(line.strip())

        self.ecID_2_transcript = transcriptID_Dict
        self.transcript_length = transcript_length_dict

    def parse(self, dir):
        """
            dir: location of the `kallisto pseudo` output
            (folder with pseudoalignments.ec, pseudoalignments.tsv etc)
        """
        assert self.ecID_2_transcript, "run create_ecID_2_transcript() first"

        dir = Path(dir)
        assert os.path.exists(dir / 'matrix.cells' )
        assert os.path.exists(dir / 'matrix.ec' )
        assert os.path.exists(dir / 'matrix.tcc.mtx' )
        assert os.path.exists(dir / 'pseudoalignments.ec' )
        assert os.path.exists(dir / 'pseudoalignments.tsv' )


        # reading the kallisto output
        ec_dict = parse_ec(dir / 'pseudoalignments.ec' )
        ec_counts = pd.read_csv(dir / 'pseudoalignments.tsv', index_col=0, sep='\t', header=None, names=['ec_id', 'counts'])

        "problem is that the ECs (just ints) are not consistant across different experiments"
        "ec_dict is currently (int)->list(int)"
        "first replace the list(int) by list(ensembl IDS)"
        "then replace the ec_id by a hash thats representative of the EC-members"

        ec_table = []

        self.ec2tr = dict()
        for ec_id, transcript_ids in ec_dict.items():
            # translate the ids into ensembl

            ensembl_transcripts = sorted([self.ecID_2_transcript[_] for _ in transcript_ids])
            ec_size = len(ensembl_transcripts)
            ensembl_transcripts_string = "_".join(ensembl_transcripts)  # join them into a single string for easy hasinh/storage

            the_hash = hashlib.sha1(ensembl_transcripts_string.encode()).hexdigest()

            ec_table.append({'ec_hash': the_hash,
                             'ec_id': ec_id,
                             'transcripts': ensembl_transcripts_string,
                              'ec_size': ec_size})
            self.ec2tr[the_hash] = set(ensembl_transcripts)
        ec_table = pd.DataFrame(ec_table)


        assert len(ec_table.ec_hash.unique()) == len(ec_table.ec_hash), "Hash collision"


        # now, in the count matrix, swap out the meaningless IDs for the hash
        new_counts = ec_counts.merge(ec_table, left_index=True, right_on='ec_id', how='left')

        # cleanup
        new_counts = new_counts.drop(['ec_id','transcripts', 'ec_size'],axis=1).set_index('ec_hash')
        return new_counts, ec_table

    def tr2ec(self):

        assert self.ec2tr
        import collections
        self.tr2ec_dict = collections.defaultdict(set)
        for ec, trset in tqdm.tqdm(self.ec2tr.items()):
            for t in trset:
                self.tr2ec_dict[t].add(ec)

    def em_quantify(self, count_df, ec_df, dir, n_iter, mode='vectorized'):

        # swap the ec_hash for the ec_id in the coutns table
        X = count_df.merge(ec_df, how='left', left_index=True, right_on='ec_hash', validate='1:1').set_index('ec_id')[['counts']].values
        transcript_length = np.array([self.transcript_length[_] for _ in range(len(self.transcript_length))]) # from dict to list
        ec_dict = parse_ec(f'{dir}/pseudoalignments.ec' )

        EM = KallistoEM(X, ec_dict, transcript_length)

        rho, logp = EM.em_iterator(n_iter, rho_init=None, reportfreq = 20)

        # turn the transcript id->ensebml into a list
        transcripts = [self.ecID_2_transcript[_] for _ in range(len(self.ecID_2_transcript))]
        rho_final = pd.DataFrame(rho[-1,:], index= transcripts, columns=['rho'])

        return rho_final, rho, logp


def transcript_to_gene():
    "see https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_hgmm_100_python/10x_hgmm_100.ipynb"
    if not (os.path.isfile('Homo_sapiens.GRCh38.94.gtf.gz')):
        os.system("curl -O ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz")
    else: print('Human GTC already downloaded!')
    os.system('gunzip -v ./Homo_sapiens.GRCh38.94.gtf.gz')



"""
xcreates the transcripts_2_gene.csv

if called like this


    with open('./Mus_musculus.GRCm38.94.gtf') as file:
        r = create_transcript_list(file, use_name = True, use_version = True)
    with open('mouse_transcript_to_gene.tsv', "w+") as output:
        print_output(output, r, use_name = True)
    print('Created mouse_transcript_to_gene.tsv file')


"""
def create_transcript_list(input, use_name = True, use_version = True):
    r = {}
    for line in input:
        if len(line) == 0 or line[0] == '#':
            continue
        l = line.strip().split('\t')
        if l[2] == 'transcript':
            info = l[8]
            d = {}
            for x in info.split('; '):
                x = x.strip()
                p = x.find(' ')
                if p == -1:
                    continue
                k = x[:p]
                p = x.find('"',p)
                p2 = x.find('"',p+1)
                v = x[p+1:p2]
                d[k] = v


            if 'transcript_id' not in d or 'gene_id' not in d:
                continue

            tid = d['transcript_id']
            gid = d['gene_id']
            if use_version:
                if 'transcript_version' not in d or 'gene_version' not in d:
                    continue

                tid += '.' + d['transcript_version']
                gid += '.' + d['gene_version']
            gname = None
            if use_name:
                if 'gene_name' not in d:
                    continue
                gname = d['gene_name']

            if tid in r:
                continue

            r[tid] = (gid, gname)
    return r



def print_output(output, r, use_name = True):
    for tid in r:
        if use_name:
            output.write("%s\t%s\t%s\n"%(tid, r[tid][0], r[tid][1]))
        else:
            output.write("%s\t%s\n"%(tid, r[tid][0]))



def parse_t2g(fname):
    "dict of transcript->gene"
    tr2g = {}
    trlist = []
    with open(fname) as f:
        for line in f:
            l = line.split()
            tr2g[l[0]] = l[1]
            trlist.append(l[0])
    return tr2g, trlist

def parse_ec(fname):
    "dict of EC->[transcript id +`````````````````````   ]\\][POIUYT]]"
    ecs = {}
    with open(fname) as f:
        for line in f:
            l = line.split()
            ec = int(l[0])
            trs = [int(x) for x in l[1].split(',')]
            ecs[ec] = trs
    return ecs

def parse_bus(t2g_file, ec_file):
    gene_min = 200
    gene_max = 10000

    #setup working directory
    import os
    os.chdir("./bus_output/")

    from subprocess import call
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import sys, collections

    tr2g, trlist = parse_t2g(t2g_file)
    genes = list(set(tr2g[t] for t in tr2g))

    # load equivalence classes
    ecs = parse_ec(ec_file)

    def ec2g(ec):
        # looks up the list of genes for a EC
        if ec in ecs:
            return list(set(tr2g[trlist[t]] for t in ecs[ec]))
        else:
            return []

    cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
    pbar=None
    pumi=None
    with open('./output.sorted.txt') as f:
        gs = set()
        for line in f:
            l = line.split()
            barcode,umi,ec,count = line.split()
            ec = int(ec)

            if barcode == pbar:
                # same barcode
                if umi == pumi:
                    # same UMI, let's update with intersection of genelist
                    gl = ec2g(ec)
                    gs.intersection_update(gl)
                else:
                    # new UMI, process the previous gene set
                    for g in gs:
                        cell_gene[barcode][g] += 1.0/len(gs)
                    # record new umi, reset gene set
                    pumi = umi
                    gs = set(ec2g(ec))
            else:
                # work with previous gene list
                for g in gs:
                    cell_gene[pbar][g] += 1.0/len(gs)

                if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                    del cell_gene[pbar]

                pbar = barcode
                pumi = umi

                gs = set(ec2g(ec))

        for g in gs:
            cell_gene[pbar][g] += 1.0/len(gs)

        if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
            del cell_gene[pbar]

    barcode_hist = collections.defaultdict(int)
    for barcode in cell_gene:
        cg = cell_gene[barcode]
        s = len([cg[g] for g in cg])
        barcode_hist[barcode] += s

    #Output a gene count histogram
    bcv = [x for b,x in barcode_hist.items() if x > gene_min and x < gene_max]
    plt.switch_backend('agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(bcv,bins=100)
    ax.set_title("Histogram")
    plt.xlabel("number of genes detected")
    plt.ylabel("number of barcodes")
    fig.savefig('gene_hist.png')

    outfile = './matrix.mtx'

    gene_to_id = dict((g,i+1) for i,g in enumerate(genes))
    barcodes_to_use = [b for b,x in barcode_hist.items() if x > gene_min and x < gene_max]

    num_entries = 0
    for barcode in barcodes_to_use:
        num_entries += len([x for x in cell_gene[barcode].values() if x>0])

    with open(outfile, 'w') as of:
        of.write('%%MatrixMarket matrix coordinate real general\n%\n')
        #number of genes
        of.write("%d %d %d\n"%(len(genes), len(barcodes_to_use), round(num_entries)))
        bcid = 0
        for barcode in barcodes_to_use:
            bcid += 1
            cg = cell_gene[barcode]
            gl = [(gene_to_id[g],cg[g]) for g in cg if cg[g] > 0]
            gl.sort()
            for x in gl:
                of.write("%d %d %f\n"%(x[0],bcid,x[1]))

    gene_names = {}
    with open("../t2g.txt") as f:
        f.readline()
        for line in f:
            t,g,gn = line.split()
            gene_names[g] = gn

    id_to_genes = dict((i,g) for (g,i) in gene_to_id.items())
    gl = []
    for i in range(1,len(genes)+1):
        g = id_to_genes[i]
        gid = g
    #    gid = g[:g.find('.')]
        if gid in gene_names:
            gn = gene_names[gid]
        else:
            gn = ''
        gl.append((g,gn))

    with open('./genes.tsv','w') as of:
        for g,gn in gl:
            of.write("%s\t%s\n"%(g,gn))

    with open('./barcodes.tsv','w') as of:
        of.write('\n'.join(x + '' for x in barcodes_to_use))
        of.write('\n')
