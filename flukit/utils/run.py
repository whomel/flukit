import typer
import glob
import pandas as pd
from rich import print
from pathlib import Path
from pandas import DataFrame
from rich.progress import track

from .variants import get_ref, get_snps, get_ha_snps, get_pa_snps, get_na_mp_snps, get_na_snps, get_mp_snps
from .align_frames import align
from .utils import write_temp_fasta, read_meta, remove_files, combine_pdfs
from .clades import run_nextclade, update_dataset
from .rename import find_fasta, rename_fasta, write_sequences, write_meta, fuzee_get
from .update_clade import get_clade_data, run_nextclade3
from .phylogenetics import select_gene, align_sequences, build_tree, plot_tree_R, check_pattern_gene, annotate_virus

def call_variants(sequences: dict, lineage: str) -> tuple[DataFrame, list]:
    results = DataFrame.from_dict(
        data = {"ha_aa":[],"NA":[],"MP":[],"PA":[], "anti-viral":[],"vacc_ref":[]}, 
        dtype='str'
            )

    ha_records = []
    for record in track(sequences, description="Processing..."):
        gene = sequences[record].gene
        try:
            seqAA, refAA = align(lineage = lineage, input_record = sequences[record])

            if len(seqAA) == 0:
                print(f"[bold yellow]{record} failed to translate. Out of frame.[/bold yellow]")
                continue            
        except Exception as error:
            print(f"[bold yellow]Error: {error} - {record} [/bold yellow]")
            break

        try:
            antiviral_str = ""
            if gene in ['PB2', 'PB1', 'NS', 'NP']:
                pass
            elif gene == 'NA':
                results.at[record, gene] = get_na_snps(seqAA, gene, lineage)
                antiviral_str = get_na_mp_snps(seqAA, refAA, gene, lineage)
            elif gene == 'MP':
                results.at[record, gene] = get_mp_snps(seqAA, gene, lineage)
                antiviral_str =  get_na_mp_snps(seqAA, refAA, gene, lineage)
            elif gene == 'PA':
                results.at[record, gene] = get_pa_snps(seqAA, refAA, lineage)
                antiviral_str = get_pa_snps(seqAA, refAA, lineage)
            elif gene == 'HA':
                results.at[record, 'ha_aa'] = get_ha_snps(seqAA, refAA)
                results.at[record, 'vacc_ref'] = get_ref(lineage)
                ha_records.append(sequences[record])
            results.at[record, 'anti-viral'] = antiviral_str
        except Exception as error:
            print(f"[bold yellow]Issue with calling variants on sample {sequences[record].id}. Error: {error}[/bold yellow]")
            pass   
    
    return(results, ha_records)

# clade calling
def call_clades(
    sequences: str, 
    lineage: str) -> DataFrame:
    '''
    Call clades with nextclade. 
    The list of `sequences` are written with clades from nextclade. 
    '''
    if lineage == 'h1n1':
        lineage = 'h1n1pdm'
    else: lineage = lineage

    get_clade_data(lineage)
    
    output = f'{lineage}_clade.tsv'

    tsv = run_nextclade3(lineage, sequences, output)
    return(tsv)


## work in progress - master function for calling directly
def findrename(
    input_dir: Path, 
    input_meta: Path, 
    output_dir: Path, 
    split_by: str,
    batch_num: str,
    rename: bool,
    ) -> None:
    '''
    Find and rename fasta files
    rename if split_by is gene or multi else do not rename output as is
    '''
    
    if batch_num and input_meta:
        raise typer.BadParameter(f"specify either but not both: batch_num/input_meta")
    if batch_num: # not implemented
        meta = fuzee_get(batch_num)
    else:
        meta = read_meta(input_meta)

    seq_num = list(meta['Seq No'])
    sequences, matched, unmatched = find_fasta(seq_num=seq_num, input_dir=input_dir)

    if rename:
        if split_by == "multi":
            sequences_renamed = rename_fasta(
                sequences=sequences, 
                meta_data=meta, 
                add_gene=True, 
                add_month=True, 
                add_passage=True
            )
        if split_by == "gene":
            sequences_renamed = rename_fasta(
                sequences=sequences, 
                meta_data=meta, 
                add_gene=False, 
                add_month=True, 
                add_passage=True
                )
    else:
        sequences_renamed = sequences

    meta_matched = meta[ meta["Seq No"].isin(matched) ]
    meta_unmatched = meta[ meta["Seq No"].isin(unmatched) ]

    write_sequences(
        sequences=sequences_renamed, 
        output=output_dir, 
        split_by=split_by
        )
    write_meta(
        meta=meta_matched, 
        output=output_dir / "meta_matched.tsv", 
        split_by='multi'
        )
    if not meta_unmatched.empty:
        write_meta(
            meta=meta_unmatched, 
            output=output_dir / "meta_unmatched.tsv", 
            split_by='multi'
            )
    
#phylogenetic analysis
def tree(input_sequences, lineages):
    """
        Function to perform the phylogenetic analysis using mfft via augur tree
        The representative virus were selected from tree traversal that matched to reference
    """
    #checking the header to get the number of gene found; then for loop over the list of available genes
    segment_object = {}
    gene_list = ['.1$','.2$','.3$','.4$','.5$','.6$','.7$','.8$']
    gene_annotate = {'.1$': 'PB2','.2$': 'PB1','.3$': 'PA','.4$': 'HA','.5$': 'NP','.6$': 'NA','.7$': 'MP','.8$': 'NS'}

    ref_name = get_ref(lineages)

    file_path = Path(input_sequences)
    input_filenames = file_path.stem

    for items in gene_list:
        if check_pattern_gene(input_sequences, items):
            
            gene_name = gene_annotate[items]

            segment_filenames = f'{input_filenames}_{gene_name}'
            
            select_gene(input_sequences, f'{segment_filenames}.fasta', items)
            align_sequences(f'{segment_filenames}.fasta',f'{segment_filenames}.aligned.fasta', gene_name, lineages)
            build_tree(f'{segment_filenames}.aligned.fasta', f'{segment_filenames}.tree')
            plot_tree_R(f'{segment_filenames}.tree', lineages, {gene_name})
            
            remove_files('.', r'.*\.tree$')
            remove_files('.', r'.*\.log$')
            remove_files('.', r'.*\.aligned.fasta$')
            remove_files('.', r'.*\.aligned.fasta.insertions.csv$')
            remove_files('.', f'{segment_filenames}.fasta$')
        else:
             continue

    csv_reader = (pd.read_csv(f) for f in glob.glob("*_representative_virus.csv"))
    combined_df = pd.concat(csv_reader, ignore_index=True)
    
    combined_df.to_csv(f"{input_filenames}_representative.csv", index=False)
    combine_pdfs('.', f'{input_filenames}_tree.pdf')
    remove_files('.', r'.*_representative_virus.csv$')
    remove_files('.', r'.*_(HA|NA|PB1|PB2|NP|MP|NS|PA)\.pdf$')
    remove_files('.', 'concat.fasta')

    return combined_df
