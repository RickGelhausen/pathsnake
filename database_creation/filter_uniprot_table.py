import pandas as pd
import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description='Extract the useful information from the xlsx file')
parser.add_argument("-t", "--tsv", action="store", dest="tsv", required=True, type=str, help= "Input tsv file")
parser.add_argument("-o", "--output", action="store", dest="output_folder", required=True, type=Path, help= "Output folder.")
args = parser.parse_args()

# Read in the data
df = pd.read_csv(args.tsv, sep='\t')

# Filter the data
df = df[df['Entry Name'].str.contains('METMA')]
df.insert(0, "ID", range(1, 1+len(df)))

def parse_goterms(go_terms):
    out_list = []
    if (go_terms == "nan"):
        return out_list

    out_list = go_terms.split("; ")
    return out_list

go_df = df[["ID", "Gene Names", "Gene Names (ordered locus)", "Gene Ontology IDs", "KEGG", "Entry"]]
go_df = go_df.astype({"ID": int, "Gene Names": str, "Gene Names (ordered locus)": str, "Gene Ontology IDs": str, "KEGG": str, "Entry": str})


rows_go_all = []
rows_refseq = []
rows_kegg = []
rows_uniprot = []
for row in go_df.itertuples(index=False, name=None):
    id = row[0]
    kegg_id = row[4]
    gid = row[1].strip().replace("'","").replace("; ", " ").replace(" ", ":")

    locus_tag = row[2]
    if ";" in locus_tag:
        locus_tags = list(filter(None,locus_tag.split(";")))
        for locus_tag in locus_tags:
            rows_refseq.append((gid, locus_tag))
    elif locus_tag not in ["", "nan"]:
        locus_tag = locus_tag.split(" ")
        rows_refseq.append((gid, locus_tag[0]))
    else:
        print("Missing refseq id: ", gid)

    if ";" in kegg_id:
        kegg_ids = list(filter(None,kegg_id.split(";")))
        for kegg_id in kegg_ids:
            rows_kegg.append((gid, kegg_id))
    elif kegg_id not in ["", "nan"]:
        rows_kegg.append((gid, kegg_id))
    else:
        print("Missing kegg id: ", gid)

    rows_uniprot.append((gid, row[5]))
    rows_go_all.extend([(gid, term, "IEA") for term in parse_goterms(row[3])])


with open(args.output_folder / "go.csv", "w") as f:
    f.write("GID,GO,EVIDENCE\n")
    for row in rows_go_all:
        f.write(f"{row[0]},{row[1]},{row[2]}\n")

with open(args.output_folder / "refseq.csv", "w") as f:
    f.write("GID,SYMBOL\n")
    for row in rows_refseq:
        f.write(f"{row[0]},{row[1]}\n")

with open(args.output_folder / "ko.csv", "w") as f:
    f.write("GID,KEGG\n")
    for row in rows_kegg:
        f.write(f"{row[0]},{row[1]}\n")

with open(args.output_folder / "uniprot.csv", "w") as f:
    f.write("GID,UNIPROT\n")
    for row in rows_uniprot:
        f.write(f"{row[0]},{row[1]}\n")

