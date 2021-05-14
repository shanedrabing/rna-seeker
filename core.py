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
    name_input = "data/mus_musculus_exp.txt"
    name_raw = "data/raw.csv"
    name_clean = "data/clean.csv"
    name_kcsv = "data/kmeans.csv"
    name_kpng = "data/kmeans.png"

    # workflow
    retrieve(name_input, name_raw, include=iszfp)
    clean(name_raw, name_clean)
    analyze(name_clean, name_kcsv, name_kpng, k=6, scaling="zscore")
