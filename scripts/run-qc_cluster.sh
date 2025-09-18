#!/bin/sh

# Remove old logs
# rm -r *html;
# rm .nextflow.log*;

export REPO_MODULE="./sc_nextflow/pipelines/0025-qc_cluster"
export STUDY_DIR="./configs/qc_cluster"


# Nextflow settings
export NXF_OPTS="-Xms25G -Xmx25G"
# Uncomment this if get strange bus errors
# export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality
export NXF_HOME=$(pwd)
export NXF_WORK="${NXF_HOME}/nextflow_work"
export NXF_TEMP="${NXF_HOME}/nextflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/nextflow_conda"
export NXF_SINGULARITY_CACHEDIR="/software/team152/nextflow/cache_singularity"

# One may need to unset R_LIBS to get the conda install
unset R_LIBS
unset R_LIBS_USER

# Check the nextflow work dir size and give warning if it is big
# First make the dir if it does not exist so the below code does not error out
[ -d ${NXF_WORK} ] || mkdir ${NXF_WORK}
NXF_WORK_SIZE=$(du -B 1 --max-depth=0 ${NXF_WORK} | awk -F ' ' '{print $1}')
# 10GB = 10737418240 bytes
# 75GB = 80530636800 bytes
if [ ${NXF_WORK_SIZE} -gt 80530636800 ]; then
    echo "WARNING: NXF_WORK dir is >75Gb. Consider running command:"
    echo -e "rm -r ${NXF_WORK}\n"
fi


/software/hgi/installs/nextflow_install/21_10_6/nextflow run \
    "${REPO_MODULE}/main.nf" \
    -profile "lsf" \
    --file_paths_10x "data/file_paths_atlas70_samples.tsv" \
    --file_metadata "data/sample_metadata.tsv" \
    --file_sample_qc "configs/qc_cluster/params-sample_qc.yml" \
    --genes_exclude_hvg "data/data-variable_gene_filter.tsv" \
    --genes_score "data/data-gene_scores.tsv" \
    --output_dir "$(pwd)/nf_results" \
    -params-file "configs/qc_cluster/params-analysis-cellbender_fpr0pt1-parameter_sweep-finalfreezeparams.yml" \
    --run_multiplet \
    -with-singularity "/software/team152/nextflow/cache_singularity/sc_qc_cluster_latest.img" \
    -with-report \
    -with-trace \
    -with-timeline \
    -with-dag flowchart.png \
    -resume

