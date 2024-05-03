#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:50:07 2024

@author: alansu
"""

'''
Opens {FbID}.SLAC.UniProt.bed files and parses as a .tsv. This is then
converted into a pandas dataframe and merged with solvent accessbility and
packing density outputs from the DSSP_output.py and WCN_output.py scripts
'''

# %% Initialize

import pandas as pd
from DSSP_output import get_dssp, make_dssp_df 
from WCN_output import get_wcn, make_wcn_df 
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression

fit_model = LinearRegression()
sns.set_theme(style="darkgrid")



file_dnds = './FBgn0000015/FBgn0000015.SLAC.UniProt.bed'
file_pdb = './FBgn0000015/refprot/FBgn0000015.pdb'
name_pdb = 'ADBD_DROME'

#file_dnds = './FBgn0000018/FBgn0000018.SLAC.UniProt.bed'
#file_pdb = './FBgn0000018/refprot/FBgn0000018.pdb'
#name_pdb = 'Q9VKM4_DROME'


# %% Parse

frame = []

labels = ['UniProt_ID',
          'Unknown',
          'AA_pos',
          'Codon_pos',
          'E[S]',
          'E[N]',
          'DS',
          'DN',
          '8',
          '9',
          '10',
          '11',
          '12',
          '13',
          '14',]

with open(file_dnds) as read:
    for line in read:
        entry = line.split() # Split by tabs
        data = entry[-1].split(":") # Split data entry by :
        data = [data[0]] + data[1].split(",") # Split data entry by ,
        
        entry = entry[:-1] + data # Stitch together
        
        if len(entry) != len(labels):# In case of missing data
            raise ValueError
        
        frame.append(entry)

df_dnds = pd.DataFrame(frame, columns = labels)

# %% Get DSSP data
sander, wilke, miller = get_dssp(name_pdb, file_pdb)
df_dssp = make_dssp_df(sander, wilke, miller)

# %% Get WCN data
output = get_wcn(name_pdb, file_pdb)
df_wcn = make_wcn_df(output)

# %% Clean and merge
'''
Merges df_dnds with df_dssp and df_wcn. df_dnds will have fewer entries as not
all residues have coverage. df_dssp and df_wcn will be matched to df_dnds by
amino acid position
'''

# Clean dssp by pulling only RASA_Wilke and renaming index
df_dssp_sub = df_dssp[['DSSP_Index', 'AA', 'Sec_Struct', 'RASA_Wilke']]
df_dssp_sub.rename(columns = {'DSSP_Index' : 'AA_pos', 'AA' : 'DSSP_AA'}, inplace = True)
df_dssp_sub['AA_pos'] = df_dssp_sub['AA_pos'].astype(int)

# Clean wcn by pulling only wcn_ca, wcn_sc and renaming index
df_wcn_sub = df_wcn[['pdb_aa', 'pdb_position', 'wcn_ca', 'wcn_sc']]
df_wcn_sub.rename(columns = {'pdb_position' : 'AA_pos', 'pdb_AA' : 'WCN_AA'}, inplace = True)
df_wcn_sub['AA_pos'] = df_wcn_sub['AA_pos'].astype(int)

# Clean dnds but dropping irrelevant entries
df_dnds_sub = df_dnds[['UniProt_ID', 'Codon_pos', 'AA_pos', 'E[S]', 'E[N]', 'DS', 'DN']] 
df_dnds_sub['AA_pos'] = df_dnds_sub[['AA_pos']].astype(int) 
df_dnds_sub[['E[S]', 'E[N]', 'DS', 'DN']] = df_dnds_sub[['E[S]', 'E[N]', 'DS', 'DN']].astype(float)

# Merge
df = pd.merge(df_dnds_sub, df_dssp_sub, how = 'left', on = 'AA_pos')
df = pd.merge(df, df_wcn_sub, how = 'left', on = 'AA_pos')
df['dN/dS'] = (df['DN']/df['E[N]'])/(df['DS']/df['E[S]'])
        
# %% dN/dS manhattan plot
df.plot('AA_pos','dN/dS')

# %% Hist for RASA
sns.histplot(data = df, x = 'RASA_Wilke')

# %% Hist for wcn_ca
sns.histplot(data = df, x = 'wcn_ca')

# %% Hist for wcn_sc
sns.histplot(data = df, x = 'wcn_sc')

# %% Scatter for wcn_ca
sns.scatterplot(data = df, x = 'wcn_ca', y = 'dN/dS')

# %% Scatter for wcn_sc
sns.scatterplot(data = df, x = 'wcn_sc', y = 'dN/dS')

# %% Scatter for RASA_Wilke
sns.scatterplot(data = df, x = 'RASA_Wilke', y = 'dN/dS')

# %% Correlations
df_sub = df[['dN/dS','wcn_sc','wcn_ca','RASA_Wilke']]
corr_matrix = df_sub.corr()['dN/dS']
print(corr_matrix)
