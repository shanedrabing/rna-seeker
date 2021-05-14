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

    name_input = "data/mus_musculus_exp.txt"
    name_raw = "data/data.csv"
    name_clean = "data/clean.csv"
    name_kcsv = "data/kmeans.csv"
    name_kpng = "data/kmeans.png"

    retrieve(name_input, "text.csv", include=lambda x: x["Symbol"].startswith("Brca"))
    # clean(name_raw, name_clean)
    # analyze(name_clean, name_kcsv, name_kpng, k=2)
