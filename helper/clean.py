import csv
from operator import sub

from util.func_utils import apply, catch, display, get, identity, nest, part, pipe
from helper.retrieve import writecsv
from util.stats_utils import geometric_mean, mean, zscore


# CONSTANTS


VALID_KEYS = ("ID", "NAME")


# FUNCTIONS


def numerical(x):
    x = float(x)
    if x == int(x):
        x = int(x)
    return x


safenum = catch(numerical, ValueError, identity)


def prune(dct):
    return {
        k.rstrip(".full_rpkm"): safenum(v) for k, v in dct.items()
        if k in VALID_KEYS or k.endswith("full_rpkm")
    }


def fields(data):
    full = nest(get(0), dict.keys, tuple)(data)
    keys = nest(set, part(sub, set(VALID_KEYS)), sorted, tuple)(full)
    vals = pipe(lambda x: pipe(get(x, iskey=False))(keys))(data)
    trans = nest(zip, tuple)(*vals)
    ind = nest(len, range)(trans)
    op = get(pipe(mean)(trans), iskey=False)
    srt = nest(part(sorted, key=op, reverse=True), tuple)(ind)
    fields = VALID_KEYS + pipe(get(keys, iskey=False))(srt)
    return fields


def clean(inp_fname, out_fname):
    with open(inp_fname) as f:
        nest(
            csv.DictReader,
            pipe(prune),
            lambda x: writecsv(x, out_fname, fields(x))
        )(f)


# TEST SCRIPT


if __name__ == "__main__":
    name_raw = "output/data.csv"
    name_clean = "output/clean.csv"

    clean(name_raw, name_clean)