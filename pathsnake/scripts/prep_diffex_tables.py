#!/usr/bin/env python
import argparse
import pandas as pd

from pathlib import Path


def parse_deltaTE(diffex_input, id_column, output_folder):
    ribo_df = diffex_input[[id_column, "RIBO_log2FoldChange", "RIBO_padj"]]
    rna_df = diffex_input[[id_column, "RNA_log2FoldChange", "RNA_padj"]]

    ribo_df.columns = ["Locus_tag", "log2FC", "padj"]
    rna_df.columns = ["Locus_tag", "log2FC", "padj"]

    ribo_df.dropna(inplace=True)
    rna_df.dropna(inplace=True)

    ribo_df.sort_values(by="log2FC", inplace=True, ascending=False)
    rna_df.sort_values(by="log2FC", inplace=True, ascending=False)

    ribo_df.to_csv(output_folder / "ribo_diffex.tsv", sep="\t", index=False)
    rna_df.to_csv(output_folder / "rna_diffex.tsv", sep="\t", index=False)

def parse_DESeq2(diffex_input, id_column, output_folder):
    rna_df = diffex_input[[id_column, "log2FC", "pvalue_adjusted"]]
    rna_df.columns = ["Locus_tag", "log2FC", "padj"]

    rna_df.dropna(inplace=True)
    rna_df.sort_values(by="log2FC", inplace=True, ascending=False)

    rna_df.to_csv(output_folder / "rna_diffex.tsv", sep="\t", index=False)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Extract the useful information from the xlsx file')
    parser.add_argument("-x", "--xlsx", action="store", dest="xlsx", required=True, type=Path, help= "Input xlsx file")
    parser.add_argument("-t", "--type", action="store", dest="input_type", required=True, type=str, help= "deltaTE or DESeq2")
    parser.add_argument("-i", "--id_column", action="store", dest="id_column", required=True, type=str, help= "Column containing the IDs to be used.")
    parser.add_argument("-o", "--output", action="store", dest="output_folder", required=True, type=Path, help= "Output folder.")

    args = parser.parse_args()
    diffex_input = pd.read_excel(args.xlsx, sheet_name=0)

    if args.input_type == "deltaTE":
        parse_deltaTE(diffex_input, args.id_column, args.output_folder)

    elif args.input_type == "DESeq2":
        parse_DESeq2(diffex_input, args.id_column, args.output_folder)

    else:
        raise ValueError("input_type must be either deltaTE or DESeq2")

if __name__ == '__main__':
    main()