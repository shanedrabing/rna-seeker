__author__ = "Shane Drabing"
__license__ = "MIT"
__version__ = "0.1.0"
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

    # analysis filenames
    k3csv = "data/k3_rpkm.csv"
    k3jpg = "data/k3_rpkm.jpg"
    k6csv = "data/k6_zscore.csv"
    k6jpg = "data/k6_zscore.jpg"

    # workflow
    retrieve(input, raw, include=iszfp)
    clean(raw, cleaned)
    analyze(cleaned, k3csv, k3jpg, k=3)
    analyze(cleaned, k6csv, k6jpg, k=6, scaling="zscore")
