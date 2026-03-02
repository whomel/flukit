#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import pandas as pd
import re
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from io import StringIO
from .utils import run_command
from importlib_resources import files

def select_gene(input_sequences, output_segment, gene_pattern):
    """Extract specific gene from sequences"""
    cmd = f"seqkit grep -n -r -p {gene_pattern} {input_sequences} -o {output_segment}"
    run_command(cmd, "Selecting gene sequences_"+ f"{input_sequences}")

def align_sequences(input_seq, output_align, gene_pattern, lineage):
    """Align sequences using augur"""
    reference = files('flukit').joinpath(f"config_tree/{lineage}/{gene_pattern}.fasta")
    concat = f"awk 1 {reference} {input_seq} > concat.fasta"
    run_command(concat, "Adding reference_"+ f"{input_seq}")

    reference_name = {
    "h3n2"    : "A/Texas/50/2012",
    "h1n1" : "6B.1_A/Michigan/45/2015", 
    "vic"     : "B/Brisbane/60/2008",
    "yam"     : "Y3_B/Phuket/3073/2013"}

    cmd = f"augur align -s concat.fasta -o {output_align} --nthreads 5 --reference-name {reference_name[lineage]}"
    run_command(cmd, "Aligning sequences_"+ f"{input_seq}")

def build_tree(input_align, output_tree):
    """Build phylogenetic tree"""
    cmd = f"augur tree -a {input_align} -o {output_tree} --nthreads 5"
    run_command(cmd, "Building phylogenetic tree_"+ f"{input_align}")

def plot_tree_R(tree_file, lineage, gene_pattern):
    """Extract specific gene from sequences"""
    script_plot = files('flukit').joinpath(f"resource/plot_tree.R")
    cmd = f"Rscript {script_plot} {tree_file} {lineage} {gene_pattern}"
    run_command(cmd, "Plotting tree and find the representative virus_"+ f"{tree_file}")

def annotate_virus(table):
    path = files('flukit').joinpath(f"config_tree/representative_virus.csv")

    right_df = pd.read_csv(path)

    # Convert both columns to string type
    table.index = table.index.astype(str)  # if left_index=True means index
    right_df['nearest_reference'] = right_df['nearest_reference'].astype(str)

    right_df = right_df.set_index('nearest_reference')
    
    # Concatenate along columns (axis=1)
    virus_table = pd.concat([table, right_df], axis=1)
    
    # Handle the .6 suffix condition if 'sample' and 'NA_clade' columns exist
    if 'sample' in virus_table.columns and 'NA_clade' in virus_table.columns:
        virus_table.loc[virus_table['sample'].str.endswith('.6'), 'NA_clade'] = np.nan
    
    return virus_table



def check_pattern_gene(filename, pattern):
    """
    Check if gene pattern exists
    
    Parameters:
    -----------
    filename : str
        Path to file
    pattern : str
        Regex pattern to search for
    
    Returns:
    --------
    bool : True if pattern found, False otherwise
    """
    if not os.path.exists(filename):
        print(f"Warning: File '{filename}' not found")
        return False
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                if re.search(pattern, line):
                    return True
        return False
    except Exception as e:
        print(f"Error reading file: {e}")
        return False
