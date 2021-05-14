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
    return display(log("Now retrieving chunk {}...", "count"))(x)


def statusok(resp):
    return resp.status_code == 200


def okgets(x):
    return nest(
        ripper(requests.get),
        filt(statusok),
        tuple
    )(x)


def flatten_dict(x, old_k=str()):
    out = dict()
    for k, v in x.items():
        new_k = "{}.{}".format(old_k, k)
        if isinstance(v, dict):
            out = {**out, **flatten_dict(v, new_k)}
        else:
            out[new_k.strip(".")] = v
    return out


def writecsv(data, fname, fields=None):
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
    html = resp.text

    match = PATTERN_TISSUES_DATA.search(html)
    if match is None:
        return dict()

    data = flatten_dict(eval(match.group(1)))
    data = {k.lower(): v for k, v in data.items()}

    data["ID"] = resp.url.split("/")[-1]
    data["NAME"] = PATTERN_GENE_NAME.search(html).group(1)

    return data


def procchunk(x):
    return nest(
        okgets,
        apply(procpage),
        filt(bool),
        tuple
    )(x)


def procfull(fname):
    return nest(
        chunker(CHUNK_SIZE),
        apply(nest(logcount, procchunk)),
        reduce(tuple.__add__),
        part(sorted, key=get("NAME")),
        part(writecsv, fname),
    )


def retrieve(inp_fname, out_fname, include=bool):
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
    procfull(out_fname)(urls)
