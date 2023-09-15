# sc_ti_atlas

A repository containing all the code, data and information required to reproduce the computational analyses in Krzak et al. 2023 including but not limited to:

* Quality control, clustering and annotation
* Differential gene expression
* Gene set enrichment analysis
* Non-negative matrix factorisation
* Heritability partionining

## Installation and data

Firstly, we highly reccomend installing mamba. It can be used like-for-like in place of conda but is much quicker. To install mamba, run the following commands:

```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

Next, git lfs is required to download the large files in one of the submodules. To install git lfs, run the following commands:
```
mamba install -c conda-forge git-lfs
```
Then clone the repository:

```
git clone --recurse-submodules git@github.com:andersonlab/sc_ti_atlas.git
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

```
mkdir data; cd data
wget https://zenodo.org/record/830100/files/discovery.h5ad?download=1 -O discovery_cohort.h5ad
wget https://zenodo.org/record/830100/files/replication.h5ad?download=1 -O replication_cohort.h5ad
```
Next up you will need to install the requisite packages. The easiest way to do this is to use the following commands to create a conda environment with all the necessary packages:

```
mamba env create -f env/conda.yml
mamba activate sc_ti_atlas
pip install cellex
```

Once you have the data and packages installed the following commands can be run to reproduce the key analyses found in our paper.

## Specifically expressed genes

```
cd sc_heritability_analysis
./src/run_CELLEX.py \
--h5_anndata ../data/discovery_cohort.h5ad \
--output_file ../out/$DATASET \
--annotation_columns predicted_celltype_machine \
--verbose True
```

## Differential gene expression analysis

```
mamba env create -f env/dge.yml
mamba activate sc_ti_dge
```

```
cd sc_nf_diffexpression
export CUR_DIR=$(pwd)
export SC_TI_OUTDIR="$CUR_DIR/out"
mkdir -p "${SC_TI_OUTDIR}"
nextflow run \
    "main.nf" \
     -profile "lsf" \
     --file_anndata "$CUR_DIR/../data/discovery_cohort.h5ad" \
     --output_dir "${SC_TI_OUTDIR}" \
     -params-file "$CUR_DIR/../configs/dge_config.yml" \
     -resume

```

## Non-negative matrix factorisation

## Heritability partitioning

```
cd sc_heritability_analysis
./src/run_CELLEX.py \
--h5_anndata ../data/discovery_cohort.h5ad \
--output_file ../out/$DATASET \
--annotation_columns predicted_celltype_machine \
--verbose True
```