import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### Config processing #####

sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),


##### Target rules #####
rule all:
    input:"results/aggregated_gwas_data.tsv",
          "results/gwas_report.html"

##### Modules #####

include: "rules/gwas_merging.smk"