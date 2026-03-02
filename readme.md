# Flukit - simple variant caller for influenza

update version: 3.0.0
Used interally at WHOFLUCC. Not recommended for external usage.

### install

1. clone this repo
2. install using pip

```
cd flukit && python -m pip install .
```
3. Install `nextclade`:

```
conda install -c bioconda nextclade
```

please update the nextclade version into 3.14.0 to be compatible with the current updates.

### usage

Run `flukit --help` to see detailed instructions

- input
	- fasta, file with headers ending in `.1`... `.4` which represent the gene. For example: `MySample.4` is the HA gene of MySample.
	- lineage, one of h1n1, h3n2, vic
	- the batch name used as prefix for output files
	- the output path

Example:

```
flukit tree -s flu.fasta -l h3n2 -b batch150 -o ~/Desktop/
```
command variant is not available with current updates.

- output
	- result.tsv file in the following format:
	```
	seqno - the fasta header
	ha_aa - HA mutations called against vacc_ref
	na - H275Y mutation
	mp - S31N mutation
	pa - I38X mutation
	antiviral_resistant - list of all antiviral found
	vacc_ref - called against ancestral strains
 	HA_clade - short clade name derived from nextclade
 	HA_subclade - subclade annotated by nextclade
	```
 	- representative.tsv file in the following format:
    ```
    nearest_reference - representative virus to each sample
    sample - the fasta header
    ```
    - phylogenetics.pdf contains the phylogenetics of each gene segments with the backbone.
      blue color represents the sample and the red color represents the reference.

Output is formatted specifically for internal database.
