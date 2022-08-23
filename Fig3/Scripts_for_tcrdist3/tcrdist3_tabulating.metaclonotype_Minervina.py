#################Meta-Clonotype Tabulation#######################################

#Apply Meta-Clonotype on Bulk TCR repertoire data
#We download a raw ImmunoSEQ file.
#Format it for use in tcrdist3.
#Search it with one of the meta-clonotypes file we made from above.

#Import library
import sys
import multiprocessing
import pandas as pd
import numpy as np
import os
from tcrdist.paths import path_to_base
from tcrdist.vdjtools_funcs import import_vdjtools
from tcrdist.repertoire import TCRrep
from tcrdist.setup_tests import download_and_extract_zip_file
from tcrdist.breadth import get_safe_chunk
from tcrdist.join import join_by_dist
from tcrdist.tabulate import tabulate
import re

#determine cpu core number
ncpus = min(multiprocessing.cpu_count(), 6)

#Import CDJtools format
#Input file must have columns ['count', 'freq', 'cdr3aa', 'v', 'j','cdr3nt'].
args = sys.argv
file = args[1]

df_bulk = import_vdjtools(   vdj_tools_file = file ,
                    chain = 'beta',
                    organism = 'human',
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    validate = True)
assert np.all(df_bulk.columns == ['count', 'freq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'cdr3_b_nucseq','valid_v', 'valid_j', 'valid_cdr3'])
#Determine rank of clonotyes to avoid collapsing different clones with same aa sequence
df_bulk = df_bulk.sort_values('count').reset_index(drop = True)
df_bulk['rank'] = df_bulk.index.to_list()

#Transfer original data into TCRrep format
tr_bulk = TCRrep(cell_df = df_bulk,
                 organism = 'human',
                 chains = ['beta'],
                 db_file = 'alphabeta_gammadelta_db.tsv',
                 compute_distances = False)

#load files containing Metaclonotypes
df_search = pd.read_csv("mydata/Metaclonotype_query/Minervina_metaclonotypes.tsv", sep = "\t") #Change query file names!!!
tr_search = TCRrep(cell_df = df_search,
                   organism = 'human',
                   chains = ['beta'],
                   db_file = 'alphabeta_gammadelta_db.tsv',
                   compute_distances = False)
	
	
#TCRs conforming to these meta-clonotypes in our dataset. 
chunk_size = get_safe_chunk(tr_search.clone_df.shape[0], tr_bulk.clone_df.shape[0])
tr_search.compute_sparse_rect_distances(
    df = tr_search.clone_df,
    df2 = tr_bulk.clone_df,
    radius = 36,
    chunk_size = chunk_size)

df_join = join_by_dist(
    how = 'inner',
    csrmat = tr_search.rw_beta,
    left_df = tr_search.clone_df,
    right_df = tr_bulk.clone_df,
    left_cols  = tr_search.clone_df.columns.to_list(),
    right_cols = tr_bulk.clone_df.columns.to_list(),
    left_suffix = '_search',
    right_suffix = '_bulk',
    max_n= 1000,
    radius = 36)

#check whether clones are witihin metaclonotype radius
df_join['RADIUS'] = df_join.apply(lambda x: x['dist'] <= x['radius_search'], axis = 1)
import re
#check whether clones posess metaclonotype motif
df_join['MOTIF'] = df_join.apply(lambda x: re.search(string = x['cdr3_b_aa_bulk'], pattern = x['regex_search']) is not None, axis = 1)

df_join['RADIUSANDMOTIF'] =  df_join['RADIUS'] & df_join['MOTIF']
df_join['unique_clones'] = 1

#Combine search result (to ensure that a clone belongs to only one metaclonotype)
output = df_join.query('RADIUSANDMOTIF').\
    sort_values('dist', ascending = True).\
    groupby(['rank_bulk']).\
    head(1)

#Calculate total frequency and total number of clones per ORF
df_join.query('RADIUSANDMOTIF').\
    sort_values('dist', ascending = True).\
    groupby(['rank_bulk']).\
    head(1).\
    groupby('protein_search').\
    sum()[['freq_bulk', 'unique_clones']]

#Output
fileout = file.replace('.pool.aaVJ.table.txt', '.minervina.metaclonotype.csv')
fileout = fileout.replace('Pool/', 'tcrdist3_output/')
output.to_csv(fileout)

#The results, which can be saved as .tsv