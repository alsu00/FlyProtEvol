#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 16:27:36 2024

@author: alansu
"""
import os
import pandas as pd
from pnps_calc import get_poly_vals

# Check files are there
def file_check(dir_path = os.getcwd()):
    
    fb_id = sorted([filename for filename in os.listdir(dir_path) if 'FBgn' in filename])
    
    frame = []
    
    labels = ['FBgn_id',
              'synonymous.Poly?',
              'missense.Poly?',
              'codonStats?',
              'features?',
              'SLAC?',
              'pdb?']
    
    # Checking all important files and outputting 1/0.
    for entry in fb_id:
        syn = int(os.path.isfile(dir_path + '/' + entry + '/' 
                               + entry + '.synonymous.Poly.UniProt.bed'))
        mis = int(os.path.isfile(dir_path + '/' + entry + '/' 
                               + entry + '.missense.Poly.UniProt.bed'))
        cstat = int(os.path.isfile(dir_path + '/' + entry + '/' 
                               + entry + '.codonStats.UniProt.bed'))
        features = int(os.path.isfile(dir_path + '/' + entry + '/' 
                               + entry + '.features.UniProt.bed'))
        dnds = int(os.path.isfile(dir_path + '/' + entry + '/' 
                               + entry + '.SLAC.UniProt.bed'))
        pdb = int(os.path.isfile(dir_path + '/' + entry + '/' + 'refprot' + '/'
                               + entry + '.pdb'))
        
        frame.append([entry, syn, mis, cstat, features, dnds, pdb])
    
    df = pd.DataFrame(frame, columns = labels)
    
    return df
    
# Aggregate pnps:
def aggregate_pnps(dir_path = os.getcwd()):
    
    frame = []
    
    check_df = file_check(dir_path) # Generate file checking dataframe for all genes
    check_df['pNpS_calc'] = 0 # Set indicator to 'uncalculated', so zero, by default
    
    
    for entry in check_df.FBgn_id:
        
        file_count = check_df.loc[check_df.FBgn_id == entry, 'synonymous.Poly?'].item() \
            + check_df.loc[check_df.FBgn_id == entry, 'missense.Poly?'].item() \
            + check_df.loc[check_df.FBgn_id == entry, 'codonStats?'].item()
            
        if file_count == 3:
            entry_data = get_poly_vals(entry, dir_path)
            if len(entry_data) > 0: # Adding an additional checker for empty files
                check_df.loc[check_df.FBgn_id == entry, 'pNpS_calc'] = 1
            frame += entry_data
            
        else:
            continue
    
    # Calculate additional values and make df
    
    labels = ['UniProt_ID', 
          'AA_pos', 
          'Codon_index',
          'Species_count',
          'E[N]', 
          'E[S]', 
          'PN', 
          'PS']
    
    df = pd.DataFrame(frame, columns = labels)
    df = df.astype({'AA_pos': int,
                    'Species_count': int,
                    'E[N]': float, 
                    'E[S]': float, 
                    'PN': float,
                    'PS': float})
    
    df['pN'] = df['PN'] / df['E[N]']
    df['pS'] = df['PS'] / df['E[S]']
    df['pN/pS'] = df.pN / df.pS
    
    return df, check_df

# Run and save pnps and checker dataframes as csv.
def test():
    df, check_df = aggregate_pnps()
    df.to_csv('pNpS_aggegate.csv')
    check_df.to_csv('file_checker.csv')
    
if __name__ == '__main__':
    pass
    
    
        
    
    
    
        
        
        

