# sc_ti_atlas

A repository containing all the code, data and information required to reproduce the computational analyses in Krzak et al. 2022 including but not limited to:

* Quality control, clustering and annotation
* Differential gene expression
* Gene set enrichment analysis
* Heritability partionining
* Velocity analysis

## Installation and data

To install, simply use this command:
```
git clone --recursive git@github.com:andersonlab/sc_ti_atlas.git
```

Note: This command will recursively clone other GitHub repositories which are nested within this one. The overarching structure is:

 ┣ 📦nf_scrna_qc
 ┣ 📦sc_heritability_analysis
 ┃ ┣ 📦CELLECT
 ┃ ┃ ┣ 📦LDSC
 ┗ 📦nf_sc_dge

The count matrices and necessary metadata can be found here:

If you would like to skip the resource intensive QC and clustering steps, an annotated AnnData file can be found here:


## Quality control and clustering

## Differential gene expression analysis

## Heritability partitioning