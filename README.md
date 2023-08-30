# sc_ti_atlas

A repository containing all the code, data and information required to reproduce the computational analyses in Krzak et al. 2023 including but not limited to:

* Quality control, clustering and annotation
* Differential gene expression
* Gene set enrichment analysis
* Non-negative matrix factorisation
* Heritability partionining

## Installation and data

To install, simply use this command:
```
git clone --recursive git@github.com:andersonlab/sc_ti_atlas.git
```

Note: This command will recursively clone other GitHub repositories which are nested within this one. The overarching structure is:
```
 ┣ 📦nf_scrna_qc
 ┣ 📦cNMF-snakemake
 ┣ 📦sc_heritability_analysis
 ┃ ┣ 📦CELLECT
 ┃ ┃ ┣ 📦LDSC
 ┗ 📦sc_nf_diffexpression
```
The count matrices and necessary metadata can be found here:

If you would like to skip the resource intensive QC and clustering steps, the auto-annoted discovery and replication AnnData objects can be found here:
10.5281/zenodo.8301000


Once you have the data the following commands can be run to reproduce the key analyses found in our paper.

```
conda activate sc_ti_atlas
```

## Quality control and clustering

## Differential gene expression analysis

## Non-negative matrix factorisation

## Heritability partitioning
