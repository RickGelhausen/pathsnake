import glob
from pathlib import Path
from snakemake.utils import min_version

min_version("7.18.2")


files = [Path(x) for x in glob.glob(f"{config['diffexFolder']}/*.xlsx")]

CONTRASTS = []
for x in files:
    CONTRASTS.append(x.stem.split("_")[0])


include: "rules/dependencies.smk"

input_type = config["diffexType"]
if input_type == "deltaTE":
    include: "rules/prep_deltaTE.smk"
    include: "rules/functional_deltaTE.smk"

elif input_type == "DESeq2":
    include: "rules/prep_DESeq2.smk"
    include: "rules/functional_DESeq2.smk"


output = []
output.append(expand("functional_analysis/{contrast}/RNA/ORA/tables/results_BP_up.tsv", contrast=CONTRASTS))
output.append(expand("functional_analysis/{contrast}/RNA/ORA/tables/results_MF_up.tsv", contrast=CONTRASTS))
output.append(expand("functional_analysis/{contrast}/RNA/ORA/tables/results_CC_up.tsv", contrast=CONTRASTS))

if input_type == "deltaTE":
    output.append(expand("functional_analysis/{contrast}/RIBO/ORA/tables/results_BP_up.tsv", contrast=CONTRASTS))
    output.append(expand("functional_analysis/{contrast}/RIBO/ORA/tables/results_MF_up.tsv", contrast=CONTRASTS))
    output.append(expand("functional_analysis/{contrast}/RIBO/ORA/tables/results_CC_up.tsv", contrast=CONTRASTS))


rule all:
    input:
        output

onsuccess:
    print("Done, no error")
