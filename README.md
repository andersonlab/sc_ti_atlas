# sc_ti_atlas

A repository containing all the code, data and information required to reproduce the computational analyses in Krzak et al. 2023 including but not limited to:

* Quality control, clustering and annotation
* Differential gene expression
* Gene set enrichment analysis
* Non-negative matrix factorisation
* bility partionining

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
 ┣ 📦sc_h2_analysis
 ┃ ┣ 📦CELLECT
 ┃ ┃ ┣ 📦LDSC
 ┗ 📦sc_nf_diffexpression
```
The count matrices and necessary metadata can be found here:

If you would like to skip the resource intensive QC and clustering steps, the auto-annoted discovery and replication AnnData objects can be found here:
https://zenodo.org/record/8379692 

```
mkdir data; cd data
wget https://zenodo.org/record/8379692/files/Discovery_cohort.h5ad?download=1 -O discovery_cohort.h5ad
wget https://zenodo.org/record/8379692/files/Replication_cohort.h5ad?download=1 -O replication_cohort.h5ad
```

Once you have our data the following commands can be run to reproduce the key analyses found in our paper.

## Specifically expressed genes

```
mamba env create -f env/heritability_env.yml
mamba activate sc_ti_heritability
```

```
cd sc_h2_analysis
DATASET="discovery"
# DATASET="replication"
./src/run_CELLEX.py \
--h5_anndata ../data/${DATASET}_cohort.h5ad \
--output_file ../out/$DATASET \
--annotation_columns predicted_celltype_machine \
--verbose True
```

## Differential gene expression analysis

```
mamba env create -f env/dge_env.yml
mamba activate sc_ti_dge
```

```
cd sc_nf_diffexpression
DATASET="discovery"
# DATASET="replication"

export CUR_DIR=$(pwd)
export SC_TI_OUTDIR="$CUR_DIR/out/discovery"
mkdir -p "${SC_TI_OUTDIR}"
nextflow run \
    "main.nf" \
     -profile "lsf" \
     --file_anndata "$(pwd)/../data/${DATASET}_cohort.h5ad" \
     --output_dir "${SC_TI_OUTDIR}" \
     -params-file "$CUR_DIR/../configs/dge_config.yml" \
     -resume

```

## Non-negative matrix factorisation

## Heritability partitioning

Once you have run the above three analyses the results can be converted into a format that can used by the CELLECT pipeline. This can be done by running the following commands:

```
cd sc_h2_analysis
Rscript src/convert_MAST.R --outpath ../out/ --dge_file ../sc_nf_diffexpression/out/differential_expression/disease_status/disease_status_dge.tsv.gz
bash src/convert_cNMF.py

```

The CELLEX marker genes are already in a format acceptable by CELLECT. The next step is to run CELLECT on the discovery and replication cohorts. To do this the following sumstats will need to be downloaded and munged:

This can be done by running the following commands:



```
cd CELLECT
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile ../../configs/bility_config.yml
```
