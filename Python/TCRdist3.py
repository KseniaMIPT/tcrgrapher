#!/usr/bin/env python
# coding: utf-8

from tcrdist.repertoire import TCRrep
from tcrdist.centers import calc_radii
from tcrdist.background import get_stratified_gene_usage_frequency
from tcrdist.sample import _default_sampler_olga
import pandas as pd
import os

def tcrdist_radii(df, cores, organism, chain):
  
  ts = _default_sampler_olga(organism = organism, chain = chain)
  ts = get_stratified_gene_usage_frequency(ts = ts, replace = True)
  
  if chain == 'beta':
    df['v_b_gene'] = df['bestVGene'] + '*01'
    df['j_b_gene'] = df['bestJGene'] + '*01'
    df['cdr3_b_aa'] = df['cdr3aa']
    cell_df = df[['count', 'freq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene']]
  elif chain == 'alpha':
    df['v_a_gene'] = df['bestVGene'] + '*01'
    df['j_a_gene'] = df['bestJGene'] + '*01'
    df['cdr3_a_aa'] = df['cdr3aa']
    cell_df = df[['count', 'freq', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene']]
    
  tr = TCRrep(
    cell_df = cell_df,
    chains = [chain],
    organism = organism,
    compute_distances=True)

  df_vj_background = tr.synthesize_vj_matched_background(ts = ts, chain = chain)

  trb = TCRrep(
    cell_df = df_vj_background.copy(),
    organism = organism, 
    chains = [chain], 
    db_file = 'alphabeta_gammadelta_db.tsv',
    compute_distances = False
    )
    
  tr.cpus = cores
  radii, thresholds, ecdfs = calc_radii(
    tr = tr, 
    tr_bkgd = trb,
    chain = chain,
    ctrl_bkgd = 10**-5
    )
  tr.clone_df['radius'] = radii
  tr.clone_df
