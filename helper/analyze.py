import collections
import csv
import random
from operator import mod, sub

from matplotlib import pyplot as plt
from util.func_utils import (apply, filt, get, getitems, identity, nest, part,
                             pipe, transpose)
from util.stats_utils import mean, minmax, relerr, zscore

from helper.clean import VALID_KEYS, safenum
from helper.retrieve import TEMPLATE_NCBI_GENE


# CONSTANTS


PALETTE_TABLEAU = (
    "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
    "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"
)


# FUNCTIONS


def palette(n):
    """Return corresponding color for a given index (mod div into range)

    :param n: Integer
    :returns: Pyplot named color
    """
    return nest(
        part(mod, len(PALETTE_TABLEAU)),
        get(PALETTE_TABLEAU, iskey=False)
    )(n)


def convert(dct):
    """Convert a dict of str values to numeric where possible

    :param dct: dict with strs from cleaned data
    :returns: dict with numeric data
    """
    return dict(zip(dct.keys(), pipe(safenum)(dct.values())))


def todatadict(data):
    """Reorganize tuple of dicts into single dict where keys are defined by ID
    and NAME, values are numerical data corresponding to RPKM values

    :param data: tuple of dicts from loaded cleaned data
    :returns: dict containing reorganized data
    """
    full = nest(get(0), dict.keys, tuple)(data)
    k = nest(filt(nest(VALID_KEYS.__contains__, lambda x: not x)), tuple)(full)
    dct = {
        (x["ID"], x["NAME"]): pipe(get(x, iskey=False))(k)
        for x in data
    }

    # remember the keys associated with the RPKM values
    dct[None] = k
    return dct


def load(fname):
    """Load from cleaned data

    :param fname: Filename of cleaned data
    :returns: dict containing reorganized data
    """
    with open(fname) as f:
        return nest(
            csv.DictReader,
            pipe(convert),
            part(todatadict)
        )(f)


def error(a, b):
    """Calculate SD between two numerical containers

    :param a: tuple of RPKM values
    :param b: tuple of RPKM values
    :returns: Standard deviation
    """
    return nest(
        apply(nest(sub, part(pow, 2))),
        sum,
        part(pow, 0.5)
    )(a, b)


def kstep(means, values):
    """Assign values (RPKM) to the mean with minimal error

    :param means: K-means
    :param values: tuples of RPKM values
    :returns: Sum within group error, assignments, and original means
    """
    sse = 0
    labels = tuple()
    for x in values:
        err, i = min((error(x, y), i) for i, y in enumerate(means))
        sse += err
        labels += i,
    return sse, labels, means


def mstep(labels, values):
    """Generate new k-means given labeled RPKM values

    :param labels: K-mean assignments
    :param values: tuples of RPKM values
    :returns: New k-means
    """
    dct = collections.defaultdict(list)
    for i, x in zip(labels, values):
        dct[i].append(x)

    return pipe(nest(transpose, pipe(mean)))(dct.values())


def plot_kmeans(fname, means, labels, values, xlabs, scaling):
    """Pyplot rendering of k-means with all original RPKM profiles as well

    :param fname: Output filename (e.g. "plot.png")
    :param means: K-means
    :param labels: K-mean assignments
    :param values: tuples of RPKM values
    :param xlabs: Labels for the X-axis (tissue names, typically)
    :param scaling: Identity, zscore, or min-max scaling
    :returns: None
    """
    for i, x in zip(labels, values):
        plt.plot(x, color=palette(i), linewidth=0.1)

    miny = float("inf")
    maxy = -float("inf")
    for i, x in enumerate(means):
        miny = min(miny, *x)
        maxy = max(maxy, *x)
        plt.plot(x, color="white", linewidth=4)
        plt.plot(x, color=palette(i), linewidth=2, label=i)

    d = (maxy - miny) * 0.05

    plt.title(f"K-Means Clustering")
    if scaling != "identity":
        plt.ylabel(f"{scaling.upper()}(RPKM)")
    else:
        plt.ylabel("RPKM")
    plt.ylim((miny - d, maxy + d))
    plt.xticks(range(len(x)), xlabs, rotation=90)
    plt.legend(labelspacing=0.2, fontsize=8)
    plt.tight_layout()
    plt.savefig(fname, dpi=600)
    plt.clf()


def analyze(inp_fname, out_fname, plt_fname=None, k=3, thresh_init=1e-2, thresh_conv=1e-4, scaling="identity"):
    """Given cleaned data (output from `clean` function), run k-means
    clustering analysis, outputting textual and optionally plot output

    :param inp_fname: Input filename of raw data
    :param out_fname: Output filename for k-means data
    :param plt_fname: (optional) Output filename for k-means plot
    :param k: (optional) Number of k-means
    :param thresh_init: (optional) Threshold for initiation
    :param thresh_conv: (optional) Threshold for convergence
    :param scaling: (optional) Identity, zscore, or min-max scaling
    :returns: None
    """
    scaling_funcs = {
        "identity": identity,
        "zscore": zscore,
        "minmax": minmax,
    }

    dct = load(inp_fname)
    keys = nest(dict.keys, filt(bool), tuple)(dct)
    rpkm = nest(part(getitems, dct), pipe(scaling_funcs[scaling]))(keys)

    sse = 1
    sse_last = 0
    while relerr(sse, sse_last) > thresh_init:
        sse_last = sse
        sse, labels, means = kstep(random.sample(rpkm, k), rpkm)

    sse = 1
    sse_last = 0
    while relerr(sse, sse_last) > thresh_conv:
        sse_last = sse
        sse, labels, _ = kstep(means, rpkm)
        means = mstep(labels, rpkm)

    if plt_fname is not None:
        plot_kmeans(plt_fname, means, labels, rpkm, dct[None], scaling)

    organs = pipe(lambda x: x.index(max(x)))(transpose(means))
    index = sorted(range(len(organs)), key=organs.__getitem__)

    temp = collections.defaultdict(list)
    for i in index:
        temp[organs[i]].append(dct[None][i])

    op = nest(sorted, "|".join, repr)

    with open(out_fname, "w") as f:
        f.write(f"ID,NAME,GROUP,ASSIGNMENT,URL\n")
        for (id_, name), group in sorted(zip(keys, labels), key=nest(reversed, tuple)):
            f.write(
                f"{id_},{name},{group},{op(temp[group])},{TEMPLATE_NCBI_GENE.format(id_)}\n")
