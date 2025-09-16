# ibdverse_terminal_ileum

A repository containing the code, data and information required to reproduce the main analyses in Krzak, Alegbe, Taylor, Jones et al. 2025 including but not limited to:

* Quality control, clustering and annotation
* Differential gene expression
* Gene set enrichment analysis
* Heritability partitioning

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
 ┣ 📦sc_h2_analysis
 ┃ ┣ 📦CELLECT
 ┃ ┃ ┣ 📦LDSC
 ┗ 📦sc_nf_diffexpression
```

## Start from either a) Initial quality control and clustering

The raw data for this project is accessible via the EGA under the following accession IDs: .... If you would like to start from the raw data you can use the nf_scrna_qc pipeline to quality control, process and cluster the data into ileal cell types. Note: this step is resource intensive and may take a long time to run, also it is not enttirely deterministic so the results may vary slightly each time it is run.

```

Once the atlas is generated, the following command can be used to generate a CellTypist model for the ileal cell types. The model can then be used to auto-annotate the pre-QC data.

``` 

## Or start from b) Annotated anndata

If you would like to skip the resource intensive QC and clustering steps, the annotated AnnData objects can be found here:
https://zenodo.org/records/14276773

```
mkdir data; cd data
wget https://zenodo.org/record/14276773/files/ti_atlas_cohort-clustered.h5ad?download=1 -O ti_atlas.h5ad
wget https://zenodo.org/record/14276773/files/ti_full_cohort-auto.h5ad?download=1 -O ti_full_cohort.h5ad
wget https://zenodo.org/record/14276773/files/organoid-auto.h5ad?download=1 -O organoid.h5ad
```

Once you have our data the following commands can be ran to reproduce the key analyses found in our paper.

## Generate the various cohorts and subsets for the analyses

## Specifically expressed genes

The following commands were used to calculate the specifically expressed genes that were used as marker genes for our cell types as well as input for the heritability enrichment analysis.

```
mamba env create -f env/heritability_env.yml
mamba activate sc_ti_heritability
```

``` bash
cd sc_h2_analysis
for DATASET in discovery replication full; do
    ./src/run_CELLEX.py \
    --h5_anndata ../data/ti_${DATASET}_cohort.h5ad \
    --output_file ../out/$DATASET \
    --annotation_columns label_machine \
    --verbose True
done
```

## Differential gene expression analysis

The following commands were used to run the differential gene expression analyssis followed by geneset enrichment analysis to determine pathways.

```
mamba env create -f env/dge_env.yml
mamba activate sc_ti_dge
```

For the TI cohorts:
```
# Define cohorts and parameters
cohorts=("Discovery" "Replication" "Full")
export FREEZE="freeze_005"
export INITIAL_DIR=$(realpath .)
export STUDY_DIR="${INITIAL_DIR}/sc_nf_diffexpression"
export REPO_MODULE="${STUDY_DIR}/sc_nf_diffexpression"

excludes_params=(
    "CD-inflamed params__TI-fr005-case-control.yml"
    "CD-uninflamed params__TI-fr005-case-control.yml"
    "Healthy-uninflamed params__TI-fr005-cds_only.yml"
    "All-uninflamed params__TI-fr005-inflamed_cd_only.yml"
    "CD-All params__TI-fr005-healthy_only.yml"
)

# Loop through each cohort and parameter combination
for cohort in "${cohorts[@]}"; do
    for pair in "${excludes_params[@]}"; do
        export COHORT="$cohort"
        export EXCLUDE=$(echo "$pair" | cut -d' ' -f1)
        export PARAMS_FILE=$(echo "$pair" | cut -d' ' -f2)
        export OUTPUT_DIR="${INITIAL_DIR}/results/${FREEZE}-cohort_${COHORT}-exclude_${EXCLUDE}/"
        
        mkdir -p ${OUTPUT_DIR}
        cd ${OUTPUT_DIR}
        echo "Running ${COHORT} cohort with ${EXCLUDE} excluded..."

        export NXF_HOME=$(pwd)
        export NXF_WORK="${NXF_HOME}/nextflow_work"
        export NXF_TEMP="${NXF_HOME}/nextflow_temp"
        
        # Submit as LSF job
        bsub -q oversubscribed -n 2 -M 10000 -R "select[mem>10000] rusage[mem=10000]" -o ${OUTPUT_DIR}/%J.out -e ${OUTPUT_DIR}/%J.err "nextflow run '${REPO_MODULE}/main.nf' -profile 'lsf' --file_anndata '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_analysis/gut/TI/anndata/ti-${COHORT}_cohort-exclude_${EXCLUDE}.h5ad' --output_dir '${OUTPUT_DIR}' -params-file '${STUDY_DIR}/${PARAMS_FILE}' -with-report -with-trace -with-timeline -with-dag flowchart.png -resume"

        cd "${INITIAL_DIR}"
    done
done
```

## Heritability partitioning


```
mamba env create -f env/heritability_env.yml
mamba activate sc_ti_heritability
```

The CELLEX marker genes are already in a format acceptable by CELLECT. The next step is to run CELLECT on the discovery, replication and full cohorts. To do this the following GWAS sumstats will need to be downloaded and munged:
* [BMI GWAS Yengo et al. (HMG, 2018)](https://academic.oup.com/hmg/article/27/20/3641/5067845)
* [Educational Attainment GWAS from Lee et al. (Nat. Gen., 2018)](https://www.nature.com/articles/s41588-018-0147-3)
* [IBD GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)
* [CD GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)
* [UC GWAS from De Lange et al. (Nat. Gen., 2017)](https://www.nature.com/articles/ng.3760)


This can be done by running the following commands:
```
mkdir gwas
wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz -P gwas/
wget https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004131/ibd_build37_59957_20161107.txt.gz -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004132/cd_build37_40266_20161107.txt.gz -P gwas/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/uc_build37_45975_20161107.txt.gz -P gwas/
```

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
--out gwas/IBD

python CELLECT/ldsc/mtag_munge.py \
--sumstats gwas/cd_build37_40266_20161107.txt.gz \
--a1 Allele1 \
--a2 Allele2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p P.value \
--out gwas/CD

python CELLECT/ldsc/mtag_munge.py \
--sumstats gwas/uc_build37_45975_20161107.txt.gz \
--a1 Allele1 \
--a2 Allele2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p P.value \
--out gwas/UC

Then to run the heritability partitioning you can use the below command:
```
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile ../../configs/h2_config.yml
```
