#!/usr/bin/env python
# coding: utf-8

from tcrdist.repertoire import TCRrep
from tcrdist.centers import calc_radii
from tcrdist.background import get_stratified_gene_usage_frequency
from tcrdist.sample import _default_tcrsampler_olga_mouse_beta
import pandas as pd
import os

def tcrdist_radii(df, cores):
  ts = _default_tcrsampler_olga_mouse_beta()
  ts = get_stratified_gene_usage_frequency(ts = ts, replace = True)
    
  df = pd.read_csv(input_path, sep='\t')
  df['v_b_gene'] = df['bestVGene'] + '*01'
  df['j_b_gene'] = df['bestJGene'] + '*01'
  df['cdr3_b_aa'] = df['cdr3aa']
    
  tr = TCRrep(
    cell_df = df[['count', 'freq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene']], 
    chains = ['beta'],
    organism = 'mouse',
    compute_distances=True)

  df_vj_background = tr.synthesize_vj_matched_background(ts = ts, chain = 'beta')

  trb = TCRrep(
    cell_df = df_vj_background.copy(),
    organism = 'mouse', 
    chains = ['beta'], 
    db_file = 'alphabeta_gammadelta_db.tsv',
    compute_distances = False
    )
    
  tr.cpus = cores
  radii, thresholds, ecdfs = calc_radii(
    tr = tr, 
    tr_bkgd = trb,
    chain = 'beta',
    ctrl_bkgd = 10**-5
    )
  tr.clone_df['radius'] = radii
  tr.clone_df
