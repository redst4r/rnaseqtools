import struct
from functools import partial
import gmpy2


def read_binary_bus(fname):
    """
    iterating over a binary busfile,
    yielding CB/UMI/EC/Counts/Flag
    """
    def decode_int_to_ACGT(the_int, seq_len):


        # the follwing line is a major performance HOG:
        # seq_decode = np.base_repr(the_int, 4) # seq_decode is a str: '10012031...'
        seq_decode = gmpy2.digits(the_int, 4)

        # note that it cannot recognized leading 0 (adenines) hence we pad 0 until we have 16BP
        seq_decode_pad = seq_decode.zfill(seq_len)
        seq_str = seq_decode_pad.replace('0', 'A').replace('1', 'C').replace('2', 'G').replace('3', 'T')
        return seq_str

    with open(fname, 'rb') as fh:
        # read the header
        # Magic (4bytes)
        # version int
        # CB len
        # umi len
        # freetext header len
        header = fh.read(20)
        magic, version, cb_len, umi_len, tlen = struct.unpack('4sIIII', header)
        assert magic == b'BUS\x00', "MAGIC doesnt match, wrong filetype??"

        # read the free header
        free_header = struct.unpack(f'{tlen}s', fh.read(tlen))

        print(f'Bustools {version}, CB length {cb_len}, UMI length {umi_len}')
        print(f'Free header  {free_header}')

        # read the bus entries
        """
        we could do a
        while True:
            bus = fh.read(32)
            if bus == b'':
                break

        but the iter(..., stop_element) is alot nice
        """
        for bus_entry in iter(partial(fh.read, 32), b''):
            cb, umi, ec, count, flags, pad = struct.unpack('QQiIII', bus_entry)
            assert pad == 0
            cb_str = decode_int_to_ACGT(cb, cb_len)
            umi_str = decode_int_to_ACGT(umi, umi_len)
            yield cb_str, umi_str, ec, count, flags


def read_text_bus(fname):
    """
    iterating over a plaintext busfile,
    yielding CB/UMI/EC/Counts/Flag
    """
    with open(fname, 'r') as fh:
        for line in fh:
            cb, umi, ec, count = read_text_bus_entry(line)
            flag = 0  # TODO Flag fixed!!
            yield cb, umi, ec, count, flag


def read_text_bus_entry(line):
    """
    read a plaintext bus
    """
    cb, umi, ec, count = line.split()
    ec = int(ec)
    count = int(count)
    return cb, umi, ec, count


def read_matrix_ec(fname):
    D = {}
    with open(fname, 'r') as fh:
        for line in fh:
            ec, transcript_list = line.split()
            ec = int(ec)
            transcripts = [int(_) for _ in transcript_list.split(',')]
            D[ec] = transcripts
    return D


def iterate_cells_of_busfile(fname, is_binary=True):
    """
    runs over the !!!SORTED!!! busfile, collecting all entries for a single CB
    and yield it as `cb,info_list`
    """
    if is_binary:
        bus_iterator = read_binary_bus(fname)
    else:
        bus_iterator = read_text_bus(fname)

    # get the first entry to get started
    cb, umi, ec, count, flag = next(bus_iterator)
    current_cell = cb
    current_info = [(umi, ec, count, flag)]

    for cb, umi, ec, count, flag in bus_iterator:
        if cb != current_cell:
            # we're finished with one cells, yield it and start the next
            yield current_cell, current_info

            # reset for the next cell
            # process results and reset
            current_cell = cb
            current_info = [(umi, ec, count, flag)]
        else:
            current_info.append((umi, ec, count, flag))

    # emitting the final cell
    yield current_cell, current_info
    return


def iterate_cells_of_busfile_old(fname):
    """
    runs over the !!!SORTED!!! busfile, collecting all entries for a single CB
    and yield it as `cb,info_list`
    """
    with open(fname, 'r') as fh:
        # first line to get going
        line = fh.readline()
        cb, umi, ec, count = read_text_bus_entry(line)

        current_cell = cb
        current_info = [(umi, ec, count)]

        for line in fh:
            cb, umi, ec, count = read_text_bus_entry(line)

            if cb != current_cell:

                # we're finished with one cells, yield it and start the next
                yield current_cell, current_info

                # reset for the next cell
                # process results and reset
                current_cell = cb
                current_info = [(umi, ec, count)]

            else:
                current_info.append((umi, ec, count))

        # emitting the final cell
        yield current_cell, current_info


def iterate_bus_cells_pairs(fname1, fname2, is_binary=True):
    """
    bustfiles must be sorted!!

    iterate over the files until we have a matchign pair of cells
    collect the info on both cells and advance to the next pair

    the pairing is tricky, given two CBs cb1, cb2:
    - if cb1 > cb2 : advance the 2nd iterator (since cb1 is in the future of that 2nd iterator)
    - if cb2 > cb1 : advance the 1st iterator (since cb2 is in the future of that 1nd iterator)
    """
    I1 = iterate_cells_of_busfile(fname1, is_binary)
    I2 = iterate_cells_of_busfile(fname2, is_binary)

    try:
        # get it started outside the loop
        cb1, info1 = next(I1)
        cb2, info2 = next(I2)
        while True:
            if cb1 == cb2:
                yield cb1, info1, info2
                # advancing both iterators
                cb1, info1 = next(I1)
                cb2, info2 = next(I2)
            elif cb1 > cb2:
                # get the next cell in I2
                cb2, info2 = next(I2)
            elif cb2 > cb1:
                cb1, info1 = next(I1)
            else:
                raise ValueError('cant happen')
    # the next() will throw an exception if one generator runs out
    # thats the signal that we're done with the pairs
    except StopIteration:
        print('One iterator finished!')
        return


if __name__ == '__main__':

    ## some speedtest
    folder = '/run/media/michi/42506642-b470-4238-be14-bb0c303b3682/cruk/bustools_testing/E14B_mgi_bus_sorted'
    text_bus = f'{folder}/output.corrected.sort.txt.bus'
    bin_bus = f'{folder}/output.corrected.sort.bus'

    I_plain = read_text_bus(text_bus)
    I_bin = read_binary_bus(bin_bus)

    import tqdm
    for _ in tqdm.tqdm(I_plain):
        pass

    for _ in tqdm.tqdm(I_bin):
        pass
