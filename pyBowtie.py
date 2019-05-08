



def hisat_single_plumbum(R12_pairs, output_bam:str, cores=1):


    bowtieCMD = '/home/michi/bowtie2-2.3.4.1-linux-x86_64/bowtie2'
    unaligned_reads = output_bam + '_unaligned.fastq'
    bowtie_bam = f'{output_bam}.bowtie.bam'
    assert mode in ['r1', 'r2', 'paired']

    r1files, r2files = zip(*r12_pairs)

    r1args = ",".join(r1files)
    r2args = ",".join(r2files)
