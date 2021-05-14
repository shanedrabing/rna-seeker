__author__ = "Shane Drabing"
__license__ = "MIT"
__version__ = "0.0.1"
__email__ = "shane.drabing@gmail.com"


from helper.analyze import analyze
from helper.clean import clean
from helper.retrieve import retrieve


if __name__ == "__main__":
    # example `include` functions

    def iszfp(x):
        return x["Symbol"].startswith("Zfp")

    def ispred(x):
        return x["Symbol"].startswith("Gm")

    def iscancer(x):
        return "cancer" in x["description"]

    # filenames
    input = "data/mus_musculus_exp.txt"
    raw = "data/raw.csv"
    cleaned = "data/clean.csv"
    kcsv = "data/kmeans.csv"
    kpng = "data/kmeans.png"

    # workflow
    retrieve(input, raw, include=iszfp)
    clean(raw, cleaned)
    analyze(cleaned, kcsv, kpng, k=6, scaling="zscore")
