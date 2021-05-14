from analyze import analyze
from clean import clean
from retrieve import retrieve


def iszfp(x):
    return x["Symbol"].startswith("Zfp")


def ispred(x):
    return x["Symbol"].startswith("Gm")


if __name__ == "__main__":
    name_input = "input/mus_musculus_exp.txt"
    name_raw = "output/data.csv"
    name_clean = "output/clean.csv"
    name_kcsv = "output/kmeans.csv"
    name_kpng = "output/kmeans.png"

    f = lambda x: "cancer" in x["description"]
    # retrieve(name_input, name_raw, include=f)
    clean(name_raw, name_clean)
    analyze(name_clean, name_kcsv, name_kpng, k=2)
