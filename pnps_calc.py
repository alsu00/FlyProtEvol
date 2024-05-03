import pandas as pd


# Parse variant list file into a dataframe
def var_parser(file_path):
    frame = []

    labels = ['UniProt_ID',
              'AA_pos',
              'AA_pos_end',
              'Data']

    with open(file_path) as read:
        for line in read:
            entry = line.split() # Split by tabs
            
            if len(entry) != len(labels):# In case of missing data
                print('Missing data')
                raise ValueError
                
            frame.append(entry)

    df = pd.DataFrame(frame, columns = labels)
    df = df.astype({"AA_pos": int, "AA_pos_end": int})
    
    return df

# Parse codonstat list file into a dataframe
def cstat_parser(file_path):
    frame = []
    
    labels = ['UniProt_ID',
              'AA_pos',
              'AA_pos_end',
              'Species',
              'Codon_index',
              'AA',
              'Codon',
              'E[N]',
              'E[S]']
    
    with open(file_path) as read:
        for line in read:
            entry = line.split() # Split by tabs
            data = entry[-1].split(":") # Split data entry by :
            data = [data[0]] + data[1].split(",") # Split data entry by ,
            
            entry = entry[:-1] + data # Stitch together

            if len(entry) != len(labels):# In case of missing data
                print('Missing data')
                raise ValueError
        
            frame.append(entry)

    df = pd.DataFrame(frame, columns = labels)
    df = df.astype({'AA_pos': int, 'AA_pos_end': int, 'E[N]': float, 'E[S]': float})
                   
    return df

# Calculate pN, pS
'''
We ignore AA_pos = 1 in codonStat, missense and synonymous files. 
This is presumed to be a start codon position in the alignment. Note that this
function outputs a df containing calculated pN/pS values if iterate=False by
default. If iterating over this function, it is convenient to keep the data in
list format and not convert to a dataframe, so if iterate=TRUE, the function 
outputs just the PN, PS, E[N], E[S] values as a nested list.
'''

def calc_pnps(mis_df, syn_df, cstat_df, iterate = False):
    frame = []
    
    labels = ['UniProt_ID', 
              'AA_pos', 
              'Codon_index',
              'Species_count',
              'E[N]', 
              'E[S]', 
              'PN', 
              'PS']
    
    # Catches cases where there is no data
    if len(cstat_df) == 0:
        
        if iterate == True:
            return frame
        
        labels = ['UniProt_ID', 
                  'AA_pos', 
                  'Codon_index',
                  'Species_count',
                  'E[N]', 
                  'E[S]', 
                  'PN', 
                  'PS',
                  'pN',
                  'pS',
                  'pN/pS']
        df = pd.DataFrame(columns = labels)
        
        return df 
    
    # If there is data, count variants
    for i in range(1, max(cstat_df.AA_pos)+1): # Skips first pos
    
        # Get observed missense (PN) and synonymous (PS)
        PN = len(mis_df[mis_df.AA_pos == i])
        PS = len(syn_df[syn_df.AA_pos == i])
        
        # Get expected missense (E_N) and synonymous (E_S) by summing across species
        cstat = cstat_df[cstat_df.AA_pos == i]
        E_N = sum(cstat['E[N]'])
        E_S = sum(cstat['E[S]'])
        
        # Check if this is a masked position (cstat missing)
        if len(cstat) == 0:
            continue
        
        # Stitch and add to frame
        entry = [cstat.UniProt_ID[:1].item(),
                 i,
                 cstat.Codon_index[:1].item(),
                 len(cstat),
                 E_N,
                 E_S,
                 PN,
                 PS]
        frame.append(entry)
    
    # Exit here and return frame if iterate = True
    if iterate == True:
        return frame
    
    # If returning dataframe for a single, calculate additional values
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
    
    return df

    
# Get pnps for a single fb_id
def get_single_pnps_csv(fb_id, dir_path, output_filepath):
    
    syn_path = dir_path + '/' + fb_id + '/' + fb_id + '.synonymous.Poly.UniProt.bed'
    mis_path = dir_path + '/' + fb_id + '/' + fb_id + '.missense.Poly.UniProt.bed'
    cstat_path = dir_path + '/' + fb_id + '/' + fb_id + '.codonStats.UniProt.bed'
    
    syn_df, mis_df, cstat_df = var_parser(syn_path), var_parser(mis_path), cstat_parser(cstat_path)
    df = calc_pnps(mis_df, syn_df, cstat_df)
    
    df.to_csv(output_filepath + '/' + fb_id + '_pnps.csv')

    return df

# Get PN, PS, E[N], E[S]
'''
This function will feed into the aggregate function to calculate pnps
for all genes in a directory
'''
def get_poly_vals(fb_id, dir_path):
    
    syn_path = dir_path + '/' + fb_id + '/' + fb_id + '.synonymous.Poly.UniProt.bed'
    mis_path = dir_path + '/' + fb_id + '/' + fb_id + '.missense.Poly.UniProt.bed'
    cstat_path = dir_path + '/' + fb_id + '/' + fb_id + '.codonStats.UniProt.bed'
    
    syn_df, mis_df, cstat_df = var_parser(syn_path), var_parser(mis_path), cstat_parser(cstat_path)
    frame = calc_pnps(mis_df, syn_df, cstat_df, iterate = True)
    
    return frame
    