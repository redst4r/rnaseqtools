def hamming_distance(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguements '''

    dist = sum([not a == b for a, b in zip(first, second)])
    return dist


def _load_whitelist(fname):
    "loads a whitelist of cellbarcdes as a set"
    with open(fname, 'r') as fh:
        whitelist = set(_.strip() for _ in fh.readlines())
    return whitelist
