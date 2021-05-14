import csv
from operator import sub

from util.func_utils import catch, get, identity, nest, part, pipe
from util.stats_utils import mean

from helper.retrieve import writecsv


# CONSTANTS


VALID_KEYS = ("ID", "NAME")


# FUNCTIONS


def numerical(x):
    """Convert to float, then reduce to int if equal

    :param x: Any, taken by float constructor
    :returns: Float or int (if equivalent)
    """
    x = float(x)
    if x == int(x):
        x = int(x)
    return x


def safenum(x):
    """Convert to numerical, if possible; otherwise, return original value

    :param x: Any
    :returns: Float or int (if equivalent)
    """
    return catch(numerical, ValueError, identity)(x)


def prune(dct):
    """Extracts "full_rpkm" from raw data

    :param dct: dict loaded from raw file
    :returns: dict without extra informaiton
    """
    return {
        k.rstrip(".full_rpkm"): safenum(v) for k, v in dct.items()
        if k in VALID_KEYS or k.endswith("full_rpkm")
    }


def sorted_fields(data):
    """Sort the fieldnames by mean RPKM

    :param data: Pruned data, tuple of dicts
    :returns: tuple of sorted fieldnames
    """
    full = nest(get(0), dict.keys, tuple)(data)
    keys = nest(set, part(sub, set(VALID_KEYS)), sorted, tuple)(full)
    vals = pipe(lambda x: pipe(get(x, iskey=False))(keys))(data)
    trans = nest(zip, tuple)(*vals)
    ind = nest(len, range)(trans)
    op = get(pipe(mean)(trans), iskey=False)
    srt = nest(part(sorted, key=op, reverse=True), tuple)(ind)
    sorted_fields = VALID_KEYS + pipe(get(keys, iskey=False))(srt)
    return sorted_fields


def clean(inp_fname, out_fname):
    """Given raw data (output from `retrieve` function), remove the unnecessary
    information and sort the fieldnames by mean RPKM

    :param inp_fname: Input filename of raw data
    :param out_fname: Output filename for cleaned data
    :returns: tuple of dicts
    """
    with open(inp_fname) as f:
        return nest(
            csv.DictReader,
            pipe(prune),
            lambda x: writecsv(x, out_fname, sorted_fields(x))
        )(f)
