from .utils import run_command
import pandas as pd
from importlib_resources import files
import os

config_path = files('flukit').joinpath('config')

def get_clade_data(lineages):
    """
        Loading HA clades from the nextclade API
    """
    cmd = f"nextclade3 dataset get -n flu_{lineages}_ha --output-dir nextclade/{lineages}_ha"
    run_command(cmd, "Loading HA clades from nextclade")

def run_nextclade3(lineages, input_sequences, output):
    """
        Extract clades informations from each sequence; default are short-clade, and ha subclade_columns
        Current available for HA segment only.
    """
    cmd = f"nextclade3 run -j 5 -D nextclade/{lineages}_ha {input_sequences} --quiet --output-tsv {output}"
    run_command(cmd, "Annotating clades_"+ f"{input_sequences}")

    df_clade = pd.read_csv(output, sep='\t')

    if lineages == 'vic':
        df_clade['short-clade'] = df_clade['clade']
    filename = os.path.splitext(os.path.basename(input_sequences))[0]
    df_clade_filtered = df_clade[["seqName", 'subclade', 'short-clade']]
    
    return df_clade_filtered
