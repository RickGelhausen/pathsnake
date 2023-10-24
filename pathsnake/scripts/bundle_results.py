#!/usr/bin/env python
import argparse
import pandas as pd

from pathlib import Path

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Bundle the results of the functional analysis into a single excel file.")
    parser.add_argument("-i", "--input_tables", nargs="+", dest="input_tables", required=True, help= "Input xlsx file")
    parser.add_argument("-m", "--method", action="store", dest="method", required=True, type=str, help= "ORA or GSEA")
    parser.add_argument("-o", "--output_xlsx", action="store", dest="output_xlsx", required=True, type=Path, help= "Output folder.")
    args = parser.parse_args()

    # read the input tables
    input_tables = [Path(x) for x in args.input_tables]

    ontology = ["BP", "MF", "CC"]
    regulation = ["up", "down"]

    if args.method == "ORA":
        ORA_dict = {}
        for file in input_tables:
            tmp = file.stem.split(".")[0].split("_")
            ont = tmp[-2]
            reg = tmp[-1]
            if (ont, reg) not in ORA_dict:
                ORA_dict[(ont, reg)] = pd.read_csv(file, sep="\t")
            else:
                print("ERROR: duplicate key")

        # create the output excel file
        writer = pd.ExcelWriter(args.output_xlsx, engine="xlsxwriter")
        for ont in ontology:
            for reg in regulation:
                cur_df = ORA_dict[(ont, reg)]
                cur_df.to_excel(writer, sheet_name=f"{ont}_{reg}", index=False, freeze_panes=(1, 0))

                for column in cur_df.columns:
                    column_width = max(cur_df[column].astype(str).map(len).max(), len(column)) + 2
                    col_idx = cur_df.columns.get_loc(column)
                    writer.sheets[f"{ont}_{reg}"].set_column(col_idx, col_idx, column_width)
        writer.close()

    elif args.method == "GSEA":
        GSEA_dict = {}
        for file in input_tables:
            tmp = file.stem.split(".")[0].split("_")
            ont = tmp[-1]
            if ont not in GSEA_dict:
                GSEA_dict[ont] = pd.read_csv(file, sep="\t")
            else:
                print("ERROR: duplicate key")

        # create the output excel file
        writer = pd.ExcelWriter(args.output_xlsx, engine="xlsxwriter")
        for ont in ontology:
            cur_df = GSEA_dict[ont]
            cur_df.to_excel(writer, sheet_name=f"{ont}", index=False, freeze_panes=(1, 0))

            for column in cur_df.columns:
                column_width = max(cur_df[column].astype(str).map(len).max(), len(column)) + 2
                col_idx = cur_df.columns.get_loc(column)
                writer.sheets[f"{ont}"].set_column(col_idx, col_idx, column_width)
        writer.close()



if __name__ == '__main__':
    main()