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
raise ValueError('this is deprecated, checkout pybustools!')
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
