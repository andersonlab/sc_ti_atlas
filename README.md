# IBDverse terminal ileum atlas and analyses

A repository containing the code, data and information required to reproduce the main analyses in Krzak, Alegbe, Taylor, Jones et al. 2024 including but not limited to:

* Quality control, clustering and annotation
* Differential gene expression
* Gene set enrichment analysis
* Heritability partitioning

## Installation and data

Firstly, we use a combination of pip and singularity to manage the software environments. To install pip, follow the instructions on the [pip installation page](https://pip.pypa.io/en/stable/installation/). And to install singularity, follow the instructions on the [singularity installation page](https://sylabs.io/guides/3.0/user-guide/installation.html).

Next, *git lfs is required* to download the large files in one of the submodules. To install git lfs, visit the [git lfs installation page](https://git-lfs.com/) and follow the instructions for your system.

After installing git lfs, you can clone this repository with the following command:

``` shell
git clone --recurse-submodules git@github.com:andersonlab/sc_ti_atlas.git
```

Note: This command will recursively clone other GitHub repositories which are nested within this one. The overarching structure is:

```
📦sc_ti_atlas
 ┣ 📦sc_nextflow
 ┣ 📦CELLECT
 ┃ ┃ ┣ 📦LDSC
 ┗ 📦sc_nf_diffexpression
```

Finally to ensure you have the required the dependencies you can use the singularity containers provided in our Zenodo archive: https://zenodo.org/record/14276773. Download the `containers.tar.gz` file and extract it into the root of this repository.

``` shell
wget https://zenodo.org/record/14276773/files/containers.tar.gz?download=1 -O containers.tar.gz
tar -xvzf containers.tar.gz
```

## Start from either a) Initial quality control and clustering

The raw data for this project is accessible via the EGA under the following accession IDs: .... If you would like to start from the raw data you can use the sc_nextflow pipeline to quality control, process and cluster the data into ileal cell types. Note: this step is resource intensive and may take a long time to run, also it is not entirely deterministic so the results may vary slightly each time it is run.

``` shell
singularity exec -B $PWD:/data containers/sc_nextflow.sif
bash scripts/run-qc_cluster.sh
```

Once the atlas is generated, the following command can be used to generate a CellTypist model for the ileal cell types. The model can then be used to auto-annotate the pre-QC data.

```shell
singularity exec -B $PWD:/data containers/sc_python.sif
python ./train_celltypist.py \
    --h5 data/ti-atlas.h5ad \
    --annotation cluster \
    --gene_symbols \
    --variable_genes_exclude data/data-variable_gene_filter.tsv
```

## Or start from b) Annotated anndata

If you would like to skip the resource intensive QC and clustering steps, the annotated AnnData objects can be found here:
https://zenodo.org/records/14276773

``` shell
mkdir -p data/anndata
wget https://zenodo.org/record/14276773/files/ti_atlas_cohort-clustered.h5ad?download=1 -O data/anndata/ti_atlas.h5ad
wget https://zenodo.org/record/14276773/files/ti_full_cohort-auto.h5ad?download=1 -O data/anndata/ti_full_cohort.h5ad
wget https://zenodo.org/record/14276773/files/organoid-auto.h5ad?download=1 -O data/anndata/organoid.h5ad
```

Once you have our data the following commands can be ran to reproduce the key analyses found in our paper.

## Generate the various cohorts and subsets for the analyses

Split the ti_full_cohort.h5ad into discovery and replication cohorts as well as the subsets for the DGE analysis.

``` shell
singularity exec -B $PWD:/data containers/sc_python.sif

# For the TI cohorts
cohorts=("Discovery" "Replication" "Full")
excludes=("CD-inflamed" "CD-uninflamed" "Healthy-uninflamed" "CD-All")
tissue="TI"
for COHORT in "${cohorts[@]}"; do
    for EXCLUDE in "${excludes[@]}"; do
    python -u subset_anndata.py \
        --adata data/anndata/ti_full_cohort.h5ad \
        --cohort $COHORT \
        --exclusions $EXCLUDE \
        --metadata "data/metadata/individual_metadata.tsv" \
        --tissue $tissue \
        --out_dir data/anndata/ 
    done
done

# For the organoid cohort
python -u subset_anndata.py \
        --adata data/anndata/organoid_cohort.h5ad \
        --cohort organoid \
        --exclusions $EXCLUDE \
        --metadata data/metadata/organoid_metadata.tsv \
        --tissue organoid \
        --out_dir data/anndata/
```

## Specifically expressed genes

The following commands were used to calculate the specifically expressed genes that were used as marker genes for our cell types as well as input for the heritability enrichment analysis.


``` shell
singularity exec -B $PWD:/data containers/sc_python.sif

cd sc_h2_analysis
for DATASET in discovery replication full; do
    ./scripts/run_CELLEX.py \
    --h5_anndata ../data/ti_${DATASET}_cohort.h5ad \
    --output_file ../out/$DATASET \
    --annotation_columns label_machine \
    --verbose True
done
```

## Differential gene expression and gene set enrichment analysis

The following commands were used to run the differential gene expression analyssis followed by geneset enrichment analysis to determine pathways.

For the TI DGE the following cohorts and exclusions were used:

``` shell
singularity exec -B $PWD:/data containers/dge.sif

# Define cohorts and parameters
cohorts=("Discovery" "Replication" "Full")
export INITIAL_DIR=$(realpath .)
export REPO_MODULE="${INITIAL_DIR}/sc_nf_diffexpression"

# Define which individuals to exclude when performing the analysis
excludes_params=(
    "CD-inflamed params__TI-fr005-case-control.yml"
    "CD-uninflamed params__TI-fr005-case-control.yml"
    "Healthy-uninflamed params__TI-fr005-cds_only.yml"
    "CD-All params__TI-fr005-healthy_only.yml"
)

# Loop through each cohort and parameter combination
for cohort in "${cohorts[@]}"; do
    for pair in "${excludes_params[@]}"; do
        export COHORT="$cohort"
        export EXCLUDE=$(echo "$pair" | cut -d' ' -f1)
        export PARAMS_FILE=$(echo "$pair" | cut -d' ' -f2)
        export OUTPUT_DIR="${INITIAL_DIR}/results/ti-cohort_${COHORT}-exclude_${EXCLUDE}/"
        
        mkdir -p ${OUTPUT_DIR}
        cd ${OUTPUT_DIR}
        echo "Running ${COHORT} cohort with ${EXCLUDE} excluded..."

        export NXF_HOME=$(pwd)
        export NXF_WORK="${NXF_HOME}/nextflow_work"
        export NXF_TEMP="${NXF_HOME}/nextflow_temp"

        nextflow run '${REPO_MODULE}/main.nf' -profile 'lsf' --file_anndata 'anndata/ti-${COHORT}_cohort-exclude_${EXCLUDE}.h5ad' --output_dir '${OUTPUT_DIR}' -params-file '${INITIAL_DIR}/configs/dge/${PARAMS_FILE}' -with-report -with-trace -with-timeline -with-dag flowchart.png -resume

        cd "${INITIAL_DIR}"
    done
done
```

For the organoid DGE the following cohorts and exclusions were used:

``` shell
singularity exec -B $PWD:/data containers/dge.sif

# Define cohorts and parameters
cohorts=("Organoid")
export INITIAL_DIR=$(realpath .)
export REPO_MODULE="${INITIAL_DIR}/sc_nf_diffexpression"
export exclude="Stimulated"
# export exclude="None"
export OUTPUT_DIR="${INITIAL_DIR}/results_organoid/${FREEZE}-cohort_organoid-exclude_${exclude}/"
export PARAMS_FILE="data/config/params__organoid-fr005-case-control.yml"
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}
echo "Running Organoid cohort with ${exclude} excluded..."
nextflow run ${REPO_MODULE}/main.nf -profile lsf \
 --file_anndata data/anndata/organoid_cohort-exclude_${exclude}.h5ad 
 --output_dir ${OUTPUT_DIR}\
 -params-file ${PARAMS_FILE}\
 -with-report -with-trace -with-timeline\
 -with-dag flowchart.png -resume
```

## Differential abundance analysis

To run the differetial abundance analysis we ran the following script:

```shell
singularity exec -B $PWD:/data containers/sc_R.sif

Rscript scripts/differential_abundance.R
```


## Heritability partitioning


The CELLEX marker genes are already in a format acceptable by CELLECT. The next step is to run CELLECT on the discovery, replication and full cohorts. To do this the following GWAS sumstats will need to be downloaded and munged:
* [BMI GWAS Yengo et al. (HMG, 2018)](https://academic.oup.com/hmg/article/27/20/3641/5067845)
* [Educational Attainment GWAS from Lee et al. (Nat. Gen., 2018)](https://www.nature.com/articles/s41588-018-0147-3)
* [IBD GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)
* [CD GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)
* [UC GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)
* [IBD GWAS from de Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)


This can be done by running the following commands:

``` shell
mkdir gwas
wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz -P gwas/
wget https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004131/ibd_build37_59957_20161107.txt.gz -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004132/cd_build37_40266_20161107.txt.gz -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/uc_build37_45975_20161107.txt.gz -P gwas/
```

Then to munge the sumstats into a format that can be used by CELLECT, run the following commands:

``` shell
singularity exec -B $PWD:/data containers/sc_python.sif

python CELLECT/ldsc/mtag_munge.py \
--sumstats gwas/GWAS_EA_excl23andMe.txt \
--merge-alleles data/ldsc/w_hm3.snplist \
--n-value 766345 \
--keep-pval \
--p PVAL \
--out gwas/EA3_Lee2018

python CELLECT/ldsc/mtag_munge.py \
--sumstats gwas/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz \
--a1 Tested_Allele \
--a2 Other_Allele \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p PVAL \
--out gwas/BMI_Yengo2018

python CELLECT/ldsc/mtag_munge.py \
--sumstats gwas/ibd_build37_59957_20161107.txt.gz \
--a1 Allele1 \
--a2 Allele2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p P.value
--out gwas/IBD_DeLange2017

python CELLECT/ldsc/mtag_munge.py \
--sumstats gwas/cd_build37_40266_20161107.txt.gz \
--a1 Allele1 \
--a2 Allele2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p P.value \
--out gwas/CD_DeLange2017

python CELLECT/ldsc/mtag_munge.py \
--sumstats gwas/uc_build37_45975_20161107.txt.gz \
--a1 Allele1 \
--a2 Allele2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p P.value \
--out gwas/UC_DeLange2017
```
Then to run the heritability partitioning you can use the below command:

``` shell
singularity exec -B $PWD:/data containers/sc_python.sif
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile ../../configs/h2_config.yml
```
