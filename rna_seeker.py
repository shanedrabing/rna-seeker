from helper.analyze import analyze
from helper.clean import clean
from helper.retrieve import retrieve


if __name__ == "__main__":
    def iszfp(x):
        return x["Symbol"].startswith("Zfp")

    def ispred(x):
        return x["Symbol"].startswith("Gm")

    def iscancer(x):
        return "cancer" in x["description"]

    name_input = "input/mus_musculus_exp.txt"
    name_raw = "output/data.csv"
    name_clean = "output/clean.csv"
    name_kcsv = "output/kmeans.csv"
    name_kpng = "output/kmeans.png"

    retrieve(name_input, name_raw, include=iscancer)
    clean(name_raw, name_clean)
    analyze(name_clean, name_kcsv, name_kpng, k=2)
