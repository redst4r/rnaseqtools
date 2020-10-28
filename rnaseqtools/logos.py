import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import tqdm

BASE_TO_INT = {'A':0,
               'T':1,
               'C':2,
               'G':3,
               'X':4,
               'N':4}

def sequence_logo(sequences):

    LOGO = np.zeros((5,20))
    for seq in tqdm.tqdm(sequences):
        for i, bp in enumerate(seq):
            LOGO[BASE_TO_INT[bp],i] += 1

    # entropy per position
    #pseudocounts
    LOGO = LOGO[:4,:]
    # LOGO = LOGO+1e-16
    probs = LOGO/ LOGO.sum(0, keepdims=True)
    
    # entropy, but take care of the probs==0 (0 *log(0) := 0)
    _q = probs * np.log2(probs)
    _q[probs==0] = 0
    entropy = -np.sum(_q, 0)
    en = (1/np.log(2)) * (4-1)/(2*len(sequences))
    information_content = np.log2(4) * (entropy + en)

    heigh = probs * information_content
    return heigh, probs

def array_to_logo_format(array):
    scores = []
    bases_only_dict = {b:i for b,i in BASE_TO_INT.items() if b in ['A','T','G','C']}
    for i in range(array.shape[1]):
        tmp = [(base, array[j,i]) for base, j in bases_only_dict.items()]
        scores.append(tmp)
    return scores


def plot_logo(sequences):
    heigh, prob2 = sequence_logo(sequences)
    return _plot_logo_from_scores(array_to_logo_format(heigh))

def _plot_logo_from_scores(all_scores):


    fp = FontProperties(family="Arial", weight="bold")
    globscale = 1.35
    LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
                "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
                "A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
                "C" : TextPath((-0.366, 0), "C", size=1, prop=fp) }
    COLOR_SCHEME = {'G': 'orange',
                    'A': 'red',
                    'C': 'blue',
                    'T': 'darkgreen'}

    def letterAt(letter, x, y, yscale=1, ax=None):
        text = LETTERS[letter]

        t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
            mpl.transforms.Affine2D().translate(x,y) + ax.transData
        p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
        if ax != None:
            ax.add_artist(p)
        return p


    fig, ax = plt.subplots(figsize=(10,3))

    x = 1
    maxi = 0
    for scores in all_scores:
        y = 0
        for base, score in scores:
            letterAt(base, x,y, score, ax)
            y += score
        x += 1
        maxi = max(maxi, y)

    plt.xticks(range(1,x))
    plt.xlim((0, x))
    plt.ylim((0, maxi))
    plt.tight_layout()
    return fig, ax


if __name__ == '__main__':

    heigh, prob = sequence_logo(['AAAA', 'TTAA'])
    plot_logo()
