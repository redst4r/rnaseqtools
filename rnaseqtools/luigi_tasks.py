import luigi
from luigi.contrib.external_program import ExternalProgramTask
from pathlib import Path
import sys
sys.path.append('/home/michi/postdoc_seattle/rnaseqtools')
from pyHISAT import hisat_single_plumbum, featureCounts_plumbum
import plumbum

class HISAT_Single_Aligner(luigi.Task):
    """

    """
    # resources = {'n_cores': 1}  # not really tied to the cores, its my way of assuring no more than one aligner runs at the same time
    # config_file = luigi.Parameter()
    outdir = luigi.Parameter()
    r12pairs = luigi.ListParameter()

    def output(self):
        fh = Path(self.outdir) / 'hisat_aligned.bam'
        return luigi.LocalTarget(str(fh))

    def requires(self):
        pass

    def run(self):
        output_bam = self.output().path
        hisat_cmd = hisat_single_plumbum(self.r12pairs, output_bam=Path(output_bam), cores=1)
        hisat_cmd()  #


class BAMIndexer(luigi.Task):
    outdir = luigi.Parameter()
    r12pairs = luigi.ListParameter()

    def output(self):
        fh = Path(self.outdir) / 'hisat_aligned.bam.bai'
        return luigi.LocalTarget(str(fh))

    def requires(self):
        yield HISAT_Single_Aligner(self.outdir, self.r12pairs)

    def run(self):
        output_bai = self.output().path

        samtools = plumbum.local['/home/michi/miniconda3_newtry/envs/mkl/bin//samtools']

        cmd = samtools['index', '-o', output_bai]
        cmd()

class featureCount_Quantifier(luigi.Task):
    """
    """
    outdir = luigi.Parameter()
    r12pairs = luigi.ListParameter()

    def run(self):
        bamfile = self.input()
        output_file = self.output().path
        fc_cmd = featureCounts_plumbum(bamfile=bamfile,
                                       output_file=Path(output_file))
        fc_cmd()


    def output(self):
        fh = Path(self.outdir) / 'hisat_featureCount.txt'
        return luigi.LocalTarget(str(fh))

    def requires(self):
        yield HISAT_Single_Aligner(self.outdir, self.r12pairs)
