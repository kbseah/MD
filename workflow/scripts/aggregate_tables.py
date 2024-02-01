import sys
import json
import os.path
import pandas as pd


if __name__ == "__main__":
    try:
        filelist = snakemake.input
        outfile = snakemake.output[0]
    except:
        filelist = sys.argv[1:]
        outfile = sys.stdout
    df_list = []
    for fn in filelist:
        bn = os.path.basename(fn)
        (tool, prefix, suffix) = bn.split(".")
        with open(fn, "r") as fh:
            dat = json.load(fh)
        df = pd.DataFrame.from_dict(dat, orient="index", columns=["value"]).reset_index(
            names=["metric"]
        )
        df["tool"] = tool
        df["sample"] = prefix
        df_list.append(df)
    df_concat = pd.concat(df_list)
    df_concat.to_csv(outfile, sep="\t", index=False)
