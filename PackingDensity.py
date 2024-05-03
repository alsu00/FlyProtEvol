#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:57:36 2024

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
structure = p.get_structure(name, file)
model = structure[0]

