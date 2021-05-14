import csv
import os
import re

import requests
from util.func_utils import (apply, chunker, display, filt, get, log, nest,
                             part, pipe, reduce, ripper)


# CONSTANTS


CHUNK_SIZE = 100
PATTERN_GENE_NAME = re.compile(r"<span.+?class=\"gn\".+?>(.+?)</span>")
PATTERN_TISSUES_DATA = re.compile(r"var\stissues_data\s=\s([\s\S]+?);")
TEMPLATE_NCBI_GENE = "https://www.ncbi.nlm.nih.gov/gene/{}"


# FUNCTIONS


def logcount(x):
    """Display the current chunk, in-line pipe

    :param x: Any, will not be used
    :returns: The original `x` variable
    """
    return display(log("Now retrieving chunk {}...", "count"))(x)


def statusok(resp):
    """Is the response status valid?

    :param resp: A requests.Response object
    :returns: Boolean
    """
    return resp.status_code == 200


def okgets(urls):
    """Multi-threaded requests.get, only returning valid response objects

    :param urls: A container of str URLs
    :returns: A tuple of requests.Response objects
    """
    return nest(
        ripper(requests.get),
        filt(statusok),
        tuple
    )(urls)


def flatten_dict(dct, old_k=str()):
    """Take nested dictionaries, combine it into one dictionary

    :param dct: A dict
    :param old_k: Previous key, for reduction
    :returns: A dict, no nested dicts within
    """
    out = dict()
    for k, v in dct.items():
        new_k = "{}.{}".format(old_k, k)
        if isinstance(v, dict):
            out = {**out, **flatten_dict(v, new_k)}
        else:
            out[new_k.strip(".")] = v
    return out


def writecsv(data, fname, fields=None):
    """Write out a container of dicts

    :param data: A container (tuple, list, etc.) of dicts
    :param fname: Filename to write to
    :param fields: (optional) A custom ordering of fieldnames
    :returns: The original `data` variable
    """
    if not data:
        raise ValueError("writecsv: no data inputted")

    if fields is None:
        fields = sorted(data[0].keys())

    filepath = os.path.realpath(fname)
    folder = os.path.dirname(filepath)
    if not os.path.exists(folder):
        os.makedirs(folder)

    with open(filepath, "w", newline=str()) as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(data)

    return data


def procpage(resp):
    """Process page; regex extraction of ENCODE RNA-Seq data

    :param resp: requests.Response object, from NCBI Gene
    :returns: dict of parsed information
    """
    html = resp.text

    match = PATTERN_TISSUES_DATA.search(html)
    if match is None:
        return dict()

    dct = flatten_dict(eval(match.group(1)))
    dct = {k.lower(): v for k, v in dct.items()}

    dct["ID"] = resp.url.split("/")[-1]
    dct["NAME"] = PATTERN_GENE_NAME.search(html).group(1)

    return dct


def procchunk(chunk):
    """Process many pages

    :param chunk: Container of requests.Response objects
    :returns: tuple of dicts
    """
    return nest(
        okgets,
        apply(procpage),
        filt(bool),
        tuple
    )(chunk)


def procfull(out_fname):
    """Process many chunks

    :param out_fname: Output filename for writing
    :returns: tuple of dicts
    """
    return nest(
        chunker(CHUNK_SIZE),
        apply(nest(logcount, procchunk)),
        reduce(tuple.__add__),
        part(sorted, key=get("NAME")),
        part(writecsv, out_fname),
    )


def retrieve(inp_fname, out_fname, include=bool):
    """Given an search export from NCBI Gene, retrieve RNA-Seq data on those
    GeneIDs, and then write out the data to an output file; the `include`
    argument is a function that will filter out particular rows of the input
    file

    :param inp_fname: Input filename of NCBI Gene search export
    :param out_fname: Output filename for writing raw data
    :param include: Function to filter out input rows
    :returns: tuple of dicts
    """
    with open(inp_fname) as f:
        ids = nest(
            part(csv.DictReader, delimiter="\t"),
            filt(include),
            apply(nest(get("GeneID"), int)),
            sorted,
            tuple
        )(f)

    urls = pipe(TEMPLATE_NCBI_GENE.format)(ids)
    n = (1 + len(urls) // CHUNK_SIZE)

    print(f"A total of {n} chunks will be retrieved.")
    return procfull(out_fname)(urls)
