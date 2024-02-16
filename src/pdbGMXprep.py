import os
import sys
import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb

# Extract the summary information from the .pka file
def extract_summary(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    start_line = "SUMMARY OF THIS PREDICTION"
    summary_started = False
    summary_lines = []

    # Extract summary lines from the .pka file
    for line in lines:
        if line.strip() == start_line:
            summary_started = True
            continue
        if "----" in line.strip():
            summary_started = False
        if summary_started:
            summary_lines.append(line.strip().split())

    # Create a dataframe with the summary data
    keys = ["resn", "resid", "chainID", "pKa", "model-pKa"]
    values = summary_lines[1:]
    pka_dict = {key: [value[i] for value in values] for i, key in enumerate(keys)}
    pka_DF = pd.DataFrame.from_dict(pka_dict)

    # Convert the resid column to numeric
    pka_DF['resid'] = pd.to_numeric(pka_DF['resid'], errors='coerce')
    # Remove rows with non-numeric resid values
    pka_DF = pka_DF.dropna(subset=['resid'])
    # Reset the index
    pka_DF.reset_index(drop=True, inplace=True)

    return pka_DF

# Determine the mutated residue name based on the pKa and pH
def mutResName(pka_DF, pH):
    residue_prot_maps = {"ASP": ("ASP", "ASP"), "GLU": ("GLU", "GLU"),
                         "LYS": ("LSN", "LYS"), "ARG": ("ARG", "ARG"),
                         "HIS": ("HIS", "HSP"), "HSD": ("HSD", "HSP"),
                         "HSE": ("HSE", "HSP")}
    pka_DF["mut_res"] = None
    for index, row in pka_DF.iterrows():
        if row["resn"] in residue_prot_maps.keys():
            if float(row["pKa"]) > pH:
                pka_DF.at[index, "mut_res"] = residue_prot_maps[row["resn"]][1]
            else:
                pka_DF.at[index, "mut_res"] = residue_prot_maps[row["resn"]][0]

class mutateRes():
    def __init__(self, df, pka_DF):
        self.df = df
        self.pka_DF = pka_DF

    def renameTERs(self):
        self.df['atom_name'].replace('OT1', 'O  ', inplace=True)
        self.df['atom_name'].replace('OT2', 'OXT', inplace=True)

    def renameCYSs(self):
        self.df['atom_name'].replace('1CB', 'CB ', inplace=True)
        self.df['atom_name'].replace('1SG', 'SG ', inplace=True)
           
    def replaceTER(self):
        ter_rows = self.df['residue_name'] == 'TER'
        ter_df = self.df[ter_rows]

        for index, row in ter_df.iterrows():
            res_num, chain = row['residue_number'], row['chain_id']
            mask = (self.df['residue_number'] == res_num) & (self.df['chain_id'] == chain) & (self.df['residue_name'] != 'TER')
        
            # If a matching residue is found, replace 'TER'
            if self.df[mask].shape[0] > 0:
                new_res_name = self.df[mask]['residue_name'].iloc[0]
                self.df.loc[index, 'residue_name'] = new_res_name

    def mutateResidues(self):
        resn_list = ["GLUP", "ASPP", "ARGN"]
        for resn in resn_list:
            res_mask = self.df['residue_name'] == resn[1:]
            unique_res_groups = self.df[res_mask].groupby(['residue_number', 'chain_id']).size().reset_index().rename(columns={0:'count'})

            for _, group in unique_res_groups.iterrows():
                res_num, chain = group['residue_number'], group['chain_id']
                mask = (self.df['residue_number'] == res_num) & (self.df['chain_id'] == chain) & (self.df['residue_name'] == resn[1:])
                if len(self.df[mask]['residue_name'].unique()) == 1:
                    new_res_name = resn[:-1]
                    self.df.loc[mask, 'residue_name'] = new_res_name
                    self.df.loc[mask, 'blank_2'] = " "
                    self.df.loc[mask, 'alt_loc'] = " "
                    
    def replaceDISU(self):
        disu_rows = self.df['residue_name'] == 'ISU'
        unique_res_nums = self.df.loc[disu_rows, 'residue_number'].unique()
        for res_num in unique_res_nums:
            mask = (self.df['residue_number'] == res_num) & (self.df['residue_name'] == 'ISU')
            if len(self.df.loc[mask, 'residue_name'].unique()) == 1:
                self.df.loc[disu_rows & (self.df['residue_number'] == res_num), 'residue_name'] = "CYS"
                self.df.loc[disu_rows & (self.df['residue_number'] == res_num), 'alt_loc'] = " "

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--pkafile', type=str, required=True, help='Path to the .pka file')
    parser.add_argument('--pdb', type=str, required=True, help='Path to the input .pdb file')
    parser.add_argument('--mutatedPDB', type=str, required=True, help='Path to the output (mutated) .pdb file')
    parser.add_argument('--pH', type=float, required=True, help='pH value used for mutation')
    args = parser.parse_args()

    pka_DF = extract_summary(args.pkafile)
    mutResName(pka_DF, args.pH)

    ppdb = PandasPdb()
    ppdb.read_pdb(args.pdb)
    mutResObj = mutateRes(ppdb.df['ATOM'], pka_DF)
    mutResObj.renameTERs()
    mutResObj.renameCYSs()
    mutResObj.replaceTER()
    mutResObj.replaceDISU()
    mutResObj.mutateResidues()

    ppdb.to_pdb(path=args.mutatedPDB, records=None, gz=False, append_newline=True)


