import collections
import csv
import random
from operator import add, inv, mod, neg, sub
from helper.retrieve import TEMPLATE_NCBI_GENE

from matplotlib import pyplot as plt

from helper.clean import safenum, VALID_KEYS
from util.func_utils import (apply, filt, get, getitems, identity, nest, part, pipe,
                        transpose)
from util.stats_utils import geometric_mean, harmonic_mean, mean, median, minmax, relerr, zscore


# CONSTANTS


PALETTE_TABLEAU = (
    "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
    "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"
)


# FUNCTIONS


tab_color = nest(
    part(mod, len(PALETTE_TABLEAU)),
    get(PALETTE_TABLEAU, iskey=False)
)


def convert(dct):
    return dict(zip(dct.keys(), pipe(safenum)(dct.values())))


def todatadict(data):
    full = nest(get(0), dict.keys, tuple)(data)
    k = nest(filt(nest(VALID_KEYS.__contains__, lambda x: not x)), tuple)(full)
    dct = {
        (x["ID"], x["NAME"]): pipe(get(x, iskey=False))(k)
        for x in data
    }
    dct[None] = k
    return dct


def load(fname):
    with open(fname) as f:
        return nest(
            csv.DictReader,
            pipe(convert),
            part(todatadict)
        )(f)


error = nest(
    apply(nest(sub, part(pow, 2))),
    sum,
    part(pow, 0.5)
)


def kstep(means, values):
    sse = 0
    labels = tuple()
    for x in values:
        err, i = min((error(x, y), i) for i, y in enumerate(means))
        sse += err
        labels += i,
    return sse, labels, means


def mstep(labels, values):
    dct = collections.defaultdict(list)
    for i, x in zip(labels, values):
        dct[i].append(x)

    return pipe(nest(transpose, pipe(mean)))(dct.values())


def plot_kmeans(fname, means, labels, values, xlabs, scaling):
    for i, x in zip(labels, values):
        plt.plot(x, color=tab_color(i), linewidth=0.1)

    miny = float("inf")
    maxy = -float("inf")
    for i, x in enumerate(means):
        miny = min(miny, *x)
        maxy = max(maxy, *x)
        plt.plot(x, color="white", linewidth=4)
        plt.plot(x, color=tab_color(i), linewidth=2, label=i)

    d = (maxy - miny) * 0.05

    plt.title(f"K-Means Clustering")
    if scaling != "identity":
        plt.ylabel(f"{scaling.upper()}(RPKM)")
    else:
        plt.ylabel("RPKM")
    plt.ylim((miny - d, maxy + d))
    plt.xticks(range(len(x)), xlabs, rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=600)
    plt.clf()


def analyze(inp_fname, out_fname, plt_fname=None, k=3, thresh_init=1e-2, thresh_conv=1e-4, scaling="identity"):
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
            f.write(f"{id_},{name},{group},{op(temp[group])},{TEMPLATE_NCBI_GENE.format(id_)}\n")


# SCRIPT


if __name__ == "__main__":
    name_clean = "data/clean.csv"
    name_kcsv = "data/kmeans.csv"
    name_kpng = "data/kmeans.png"

    analyze(name_clean, name_kcsv, name_kpng, k=3)
