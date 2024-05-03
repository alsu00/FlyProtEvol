#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:34:40 2024

@author: alansu
"""

# %% Initialize
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from matplotlib import pyplot as plt
import pandas as pd

p = PDBParser()

file = './FBgn0000015/refprot/FBgn0000015.pdb'
name = 'ADBD_DROME'

 # %% Getting all dssp methods
def get_dssp(name, file):
    structure = p.get_structure(name, file)
    model = structure[0]
    dssp_sander = DSSP(model, file, dssp = 'mkdssp', acc_array = 'Sander')
    dssp_wilke = DSSP(model, file, dssp = 'mkdssp', acc_array = 'Wilke')
    dssp_miller = DSSP(model, file, dssp = 'mkdssp', acc_array = 'Miller')
    return dssp_sander, dssp_wilke, dssp_miller

# %% Stitch into a dataframe
def make_dssp_df(dssp_sander, dssp_wilke, dssp_miller):
    key_list = dssp_sander.keys() # List of tuples that will serve as keys for dssp objects

    frame = [[] for i in range(16)]

    for key in key_list:
        for i in range(14): # No of entries in DSSP object
            frame[i].append(dssp_sander[key][i])
        frame[14].append(dssp_wilke[key][3]) # Relative ASA in Wilke object
        frame[15].append(dssp_miller[key][3]) # Relative ASA in Miller object
    
    labels = ['DSSP_Index',
              'AA',
              'Sec_Struct',
              'RASA_Sander',
              'Phi',
              'Psi',
              'NH->O_1_relidx',
              'NH–>O_1_energy',
              'O–>NH_1_relidx',
              'O–>NH_1_energy',
              'NH–>O_2_relidx',
              'NH–>O_2_energy',
              'O–>NH_2_relidx',
              'O–>NH_2_energy',
              'RASA_Wilke',
              'RASA_Miller']

    frame_dict = {}

    for i in range(16):
        frame_dict[labels[i]] = frame[i]

    df = pd.DataFrame(frame_dict)

    return df

# %% Run
sander, wilke, miller = get_dssp(name, file)
df_ADBD = make_dssp_df(sander, wilke, miller)

# %% Check correlation of different ASA normalization methods
df_ADBD.plot("RASA_Sander","RASA_Wilke",'scatter')
df_ADBD.plot("RASA_Sander","RASA_Miller",'scatter')
df_ADBD.plot("RASA_Wilke","RASA_Miller",'scatter')
df_ADBD.hist("RASA_Wilke")
df_ADBD.hist('RASA_Sander')
df_ADBD.hist('RASA_Miller')