# sc_ti_atlas

A repository containing all the code, data and information required to reproduce the computational analyses in Krzak et al. 2023 including but not limited to:

* Quality control, clustering and annotation
* Differential gene expression
* Gene set enrichment analysis
* Non-negative matrix factorisation
* Heratibility partionining

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
PENDING

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

## Differential abundance analysis

```

## Heritability partitioning

Once you have run the above three analyses the results can be converted into a format that can used by the CELLECT pipeline. This can be done by running the following commands:
```
mamba env create -f env/heritability_env.yml
mamba activate sc_ti_heritability
```

```
cd sc_h2_analysis
Rscript src/convert_MAST.R --outpath ../out/ --dge_file ../sc_nf_diffexpression/out/differential_expression/disease_status/disease_status_dge.tsv.gz
bash src/convert_cNMF.py

```

The CELLEX marker genes are already in a format acceptable by CELLECT. The next step is to run CELLECT on the discovery and replication cohorts. To do this the following GWAS sumstats will need to be downloaded and munged:
* [BMI GWAS Yengo et al. (HMG, 2018)](https://academic.oup.com/hmg/article/27/20/3641/5067845)
* [Educational Attainment GWAS from Lee et al. (Nat. Gen., 2018)](https://www.nature.com/articles/s41588-018-0147-3)
* [CD GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)
* [UC GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)
* [IBD GWAS from de Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)

This can be done by running the following commands:
```
mkdir gwas
wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz -P gwas/
wget https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004131/ibd_build37_59957_20161107.txt.gz -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004132/cd_build37_40266_20161107.txt.gz -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/uc_build37_45975_20161107.txt.gz -P gwas/
```

Then to munge the sumstats into a format that can be used by CELLECT, run the following commands:
```

python CELLECT/ldsc/mtag_munge.py \
--sumstats example/GWAS_EA_excl23andMe.txt \
--merge-alleles data/ldsc/w_hm3.snplist \
--n-value 766345 \
--keep-pval \
--p PVAL \
--out example/EA3_Lee2018


python CELLECT/ldsc/mtag_munge.py \
--sumstats example/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz \
--a1 Tested_Allele \
--a2 Other_Allele \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p PVAL \
--out example/BMI_Yengo2018


python CELLECT/ldsc/mtag_munge.py \
--sumstats example/cd_build37_40266_20161107.txt.gz \
--a1 A1 \
--a2 A2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--n-value 20064 \
--keep-pval \
--p P \
--out example/CD_DeLange2017

python CELLECT/ldsc/mtag_munge.py \
--sumstats example/uc_build37_45975_20161107.txt.gz \
--a1 A1 \
--a2 A2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--n-value 27046 \
--keep-pval \
--p P \
--out example/UC_DeLange2017
```

Then to run the heritability partitioning you can use the below command:
```
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile ../../configs/bility_config.yml
```
