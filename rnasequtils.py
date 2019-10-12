import os
import gzip
import tempfile
import tqdm
import pathlib
import pysam

def Phred_symbol2errorprob(symbol:str):
    "Phred ascii reprensentation to probability of error"
    Q = ord(symbol) - 33
    P = 10**(-Q/10)
    # logP = -Q/10
    return P

def Phred2symbol(phred:str):
    "phred score to ascii"
    return str(chr(phred+33))

def merge_zip_files(zipfiles:list, targetzip):
    """
    weird naming thing: the file inside the targetzip will get the same name as the archive (seemsl ike a gzip proglem)
    workaround: targetzip must end with .gz. we create a zipfile without the .gz extension (file inside also no extension)
    and then we rename the zipfile
    """
    assert isinstance(zipfiles, list) or  isinstance(zipfiles, tuple), "input must be as list/tuple"
    assert targetzip.endswith('.gz')

    targetzip = targetzip.replace('.gz' ,'')

    # with gzip.open(targetzip, 'wb') as zip_target:
    with gzip.GzipFile(filename=targetzip, mode='w', compresslevel=1) as zip_target:
        for zf in tqdm.tqdm(zipfiles):
            with gzip.open(zf, 'rb') as zip1:
                 zip_target.writelines(zip1.readlines())  # TODO bad: reasds the entire file at once!!

    os.rename(targetzip, targetzip+'.gz')

# def merge_zip_files_fast(zipfiles:list, targetzip):
# """
# ACTUALLY NOT FASTER AT ALL
# """
#
#     """
#     instead of unpacking and writing into target on the fly,
#     first, unpack all (at the expense of disk-space), merge the uncompressed
#     file, and compress the whole thing
#     """
#     assert isinstance(zipfiles, list) or  isinstance(zipfiles, tuple), "input must be as list/tuple"
#     assert targetzip.endswith('.gz')
#
#     targetzip = targetzip.replace('.gz' ,'')
#     _, tmp_name = tempfile.mkstemp()
#
#     # decomp all zipfiles into the uncompressed single
#     with open(tmp_name, 'wb') as tmp_fh:
#         for zf in tqdm.tqdm(zipfiles, desc='Decompressing singles'):
#             with gzip.open(zf, 'r') as zip1:
#                  tmp_fh.writelines(zip1.readlines())  # TODO bad: reasds the entire file at once!!
#
#     # zip the big single file
#     with gzip.GzipFile(filename=targetzip, mode='w', compresslevel=1) as zip_target:
#     # with gzip.open(targetzip, 'wb') as zip_target:
#         with open(tmp_name, 'rb') as tmp_fh:
#             zip_target.writelines(tmp_fh)
#
#     os.remove(tmp_name)
#
#     os.rename(targetzip, targetzip+'.gz')

def bamfile_index(bamfile):
    os.system(f'samtools index {bamfile}')


def bamfile_sort_index(bamfilename, outfile=None, add_flags=''):

    f = pathlib.Path(bamfilename)
    basename = f.stem
    folder = str(f.parent)
    fullpath_and_name = str(bamfilename)

    # we sort inplace, hence the sort first goes into a separate file
    #

    sorted_name = f'{folder}/{basename}.sorted.bam' if outfile is None else outfile
    os.system(f'samtools sort {add_flags} {fullpath_and_name} > {sorted_name}')
    bamfile_index(sorted_name)
    return sorted_name


def rename_contigs(bamfile_name, outname, contig_rename_dict):
    """
    sometimes the contig names get messed up: chr1 vs 1 etc
    """
    with pysam.AlignmentFile(bamfile_name, "rb") as bamfile:
        # fix the header
        newheader = bamfile.header.to_dict().copy()
        for contig_dict in newheader['SQ']:
            if contig_dict['SN'] in contig_rename_dict.keys():
                contig_dict['SN'] = contig_rename_dict[contig_dict['SN']]
        newheader = pysam.AlignmentHeader.from_dict(newheader)
        #fix the reads
        with pysam.AlignmentFile(outname, "wb", header=newheader) as out:
            for read in tqdm.tqdm(bamfile.fetch(until_eof=True)):
                newread = read.to_dict().copy()
                if newread['ref_name'] in contig_rename_dict.keys():
                    newread['ref_name'] = contig_rename_dict[newread['ref_name']]
                    newread = pysam.AlignedSegment.from_dict(newread, newheader)
                    out.write(newread)




def main():
    import pathlib
    from collections import defaultdict
    D = pathlib.Path('/home/michi/osiris_mount/134338207_10xrun')

    exp_dict = defaultdict(list)

    for d in D.iterdir():
        if d.is_dir() and not f.name == 'merged_lanes':
            R1, R2 = None, None
            for g in d.iterdir():
                if g.name.endswith('_R1_001.fastq.gz'):
                    R1 = g
                if g.name.endswith('_R2_001.fastq.gz'):
                    R2 = g
            exp, lane_id = d.name.split('_L00')
            if R1 and R2:
                exp_dict[exp].append((R1, R2))
            else:
                raise ValueError(f'R1/R2 pair not found for {f}')

    outdir = '/home/michi/osiris_mount/134338207_10xrun/'

    args = []
    for exp, R1R2_list in exp_dict.items():
        R1list, R2list = zip(*R1R2_list)
        r1out = f'{outdir}/merged_lanes/{exp}_R1_001.fastq.gz'
        r2out = f'{outdir}/merged_lanes/{exp}_R2_001.fastq.gz'

        if not os.path.exists(r1out):
            args.append((R1list, r1out))
            # merge_zip_files(R1list, r1out)
        if not os.path.exists(r2out):
            args.append((R2list, r2out))
            # merge_zip_files(R2list, r2out)

    import multiprocessing as mp
    [merge_zip_files(rlist, rout) for rlist, rout in args]

    with mp.Pool(4) as pool:
        pool.starmap(merge_zip_files, args)

if __name__ == '__main__':
    main()
