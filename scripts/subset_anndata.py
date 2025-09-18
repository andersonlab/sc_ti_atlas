Ωimport pandas as pd
import numpy as np
import scanpy as sc
import argparse

# # Parse arguments
parser = argparse.ArgumentParser(description='Subset anndata object by cohort and save.')
parser.add_argument('--adata', type=str, help='Path to anndata object.', required=True)
parser.add_argument('--cohort', type=str, help='Cohort to subset by.', required=True)
parser.add_argument('--exclusions', type=str, help='Comma-separated list of subsets of individuals should be excluded', required=True)
parser.add_argument('--tissue', type=str, help='Tissue shorthand: Either r, ti, or blood', required=True)
parser.add_argument('--metadata', type=str, help='The metadata tsv to be joined to the anndata object', required=True)
parser.add_argument('--out_dir', type=str, help='Directory to save anndata object.', required=True)

args = parser.parse_args()
adata_path = args.adata
cohort = args.cohort
exclude_list = args.exclusions.split(',')
metadata_path = args.metadata
tissue = args.tissue
out_dir = args.out_dir

# Read data
print('Reading data...')
adata = sc.read_h5ad(adata_path)


# Set counts as a layer
print('Checking counts are counts...')
print(np.sum(adata.X[1]))
adata.layers['counts'] = adata.X.copy()

# Add sangar sample id column
adata.obs['sanger_sample_id'] = adata.obs['Exp']

print(f'Full dataset has {adata.shape[0]} cells and {adata.obs["sanger_sample_id"].nunique()} individuals.')

metadata = pd.read_table(metadata_path, sep='\t')

# Slice the metadata to include only samples where tissue matches the filename
metadata = metadata[metadata['biopsy_type'] == tissue]

adata.obs = adata.obs.reset_index().merge(metadata,how='left', left_on='sanger_sample_id', right_on='sanger_sample_id').set_index('index')
adata.obs.index = adata.obs.index.astype(str)
adata.obs.index.name = None

adata = adata[adata.obs['pct_counts_gene_group__mito_protein'] < 80]
    

# Subset by cohort
adata.obs.index = adata.obs.index.astype(str)
if cohort == 'Full':
    out_adata = adata[adata.obs['cohort'].isin(['Discovery', 'Replication'])]
else:
    out_adata = adata[adata.obs['cohort'] == cohort]

del adata

# Save the data
print(f'Saving {out_adata.shape[0]} cells and {out_adata.obs["sanger_sample_id"].nunique()} samples divided into {out_adata.obs["disease_status"].value_counts(dropna=False)} disease statuses and {out_adata.obs["ses_inflamed"].value_counts(dropna=False)} inflammation statuses.')
# Save h5ad
out_adata.write_h5ad(f'{out_dir}/{tissue}-{cohort}_cohort.h5ad',compression='gzip')
out_adata.obs.index = out_adata.obs.index.astype(str)
# Save obs only as csv
out_adata.obs.to_csv(f'{out_dir}/{tissue}-{cohort}_obs.csv')

if tissue == 'ti':
    for data_class in exclude_list:
        disease = data_class.split('-')[0]
        inflammation = data_class.split('-')[1]

        # Slicing
        if disease == 'All':
            exclude_out_adata = out_adata[out_adata.obs['ses_inflamed'] != inflammation]
        elif inflammation == 'All':
            exclude_out_adata = out_adata[out_adata.obs['disease_status'] != disease]
        else:
            exclude_out_adata = out_adata[~((out_adata.obs['disease_status'] == disease) & (out_adata.obs['ses_inflamed'] == inflammation))]
        
        print(f'Saving without {data_class}. Dataset to be saved has {exclude_out_adata.shape[0]} cells and {exclude_out_adata.obs["sanger_sample_id"].nunique()} samples divided into {exclude_out_adata.obs["disease_status"].value_counts(dropna=False)} disease statuses and {exclude_out_adata.obs["ses_inflamed"].value_counts(dropna=False)} inflammation statuses.')
        exclude_out_adata.write_h5ad(f'{out_dir}/{tissue}-{cohort}_cohort-exclude_{data_class}.h5ad',compression='gzip')
        del exclude_out_adata
if tissue == 'organoid':
    for data_class in exclude_list:
        exclude_out_adata = out_adata[out_adata.obs['stimulation'] != data_class]
        print(f'Saving without {data_class}. Dataset to be saved has {exclude_out_adata.shape[0]} cells and {exclude_out_adata.obs["sanger_sample_id"].nunique()} samples divided into {exclude_out_adata.obs["disease_status"].value_counts(dropna=False)} disease statuses and {exclude_out_adata.obs["stimulation"].value_counts(dropna=False)} stimulations.')
        exclude_out_adata.write_h5ad(f'{out_dir}/{tissue}-{cohort}_cohort-exclude_{data_class}.h5ad',compression='gzip')
        del exclude_out_adata
print('Done!')