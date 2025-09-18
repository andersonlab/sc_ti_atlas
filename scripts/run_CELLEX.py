#!/usr/bin/env python

__author__ = 'Moritz Przybilla, Tobi Alegbe'
__date__ = '2020-10-27'
__version__ = '0.0.1'

import argparse
import os
import random
import numpy as np
import scipy as sp
from distutils.version import LooseVersion
import cellex
import pandas as pd
import scanpy as sc
import csv
import time
from datetime import timedelta

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def run_CELLEX(
    adata,
    annotation_columns,
    output_file,
    verbose
):
    """Use CELLEX to prioritize marker genes for clusters.

    Parameters
    ----------
    adata : AnnData
        Input AnnData file.
    annotation_columns : list of string
        List of columns in the adata.obs that contain relevant cluster annotations.
    output_file : string
        Basename of output_file, will have -normalized_pca.h5ad appended to it.
    verbose : string
        Specify whether summary should be printed or not using True or False.

    Returns
    -------
    output_file : string
        output_file
    """

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Check all user-specified annotation columns are in the dataset
    for col in annotation_columns:
        if col in adata.obs.columns:
            if verbose:
                print('Found provided annotation column {}'.format(col))
        else:
            raise Exception('Provided column name {} is NOT present in h5 AnnData file obs'.format(col))

    # create a dataframe from the count matrix
    if verbose:
        print('Creating pandas dataframe')
    data = pd.DataFrame(adata.layers['counts'].toarray(
    ), index=adata.obs_names, columns=adata.var_names).transpose()



    for col in annotation_columns:
        
        if verbose:
            print('Processing annotations in {}'.format(col))

        # establish the metadata file from the annData
        # assign the column with the cluster id
        metadata = pd.DataFrame(
            data={col: adata.obs[col]}, index=adata.obs_names)
        
        # Slice for cells which do not have an annotation
        meta_not_na = metadata[col].notnull()
        if np.sum(meta_not_na) < len(metadata):
            if verbose:
                print('Dropping {} cells from {} as they are NaNs'.format(len(metadata) - np.sum(meta_not_na),col))
            metadata = metadata.loc[meta_not_na]
            data = pd.DataFrame(adata.layers['counts'].toarray(
                    ), index=adata.obs_names, columns=adata.var_names).loc[meta_not_na].transpose()

        # Sanitise annotation names, remove special characters
        num_annots = len(np.unique(metadata[col]))
        metadata[col] = metadata[col].str.replace(r'\s', '_')
        metadata[col] = metadata[col].str.replace(r'[^0-9a-zA-Z+\-\(\)_]+', '')
        
        # Check sanitation hasn't merged any annotations
        if num_annots != len(np.unique(metadata[col])):
            raise Exception(
                'Annotation names for {} could be automatically sanitised,\
                     please remove spaces and special characters manually.'.format(col))


        # check if the data shape matches up
        if (data.shape[1] == len(metadata)):

            # The ESObject encapsulates the core features of CELLEX. We set "verbose=True" to get some progress reports.
            # The computations may take a while depending on the data and available computational power.
            if verbose:
                print('Running CELLEX on {}'.format(col))
            eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
            eso.compute(verbose=True)

        else:
            raise Exception(
                'The number of cells provided in the CELLEX data {} does not match the number of cells in the metadata {}.'.format(data.shape[1],len(metadata)))

        # save the CELLEX results

        df = eso.results['esmu']
        df.to_csv(
            "{}.{}.esmu.csv.gz".format(output_file, col),
            compression='gzip',
            index=True,
            header=True
        )

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read anndata object. Run CELLEX and save the expression specificity
            object to file.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-annot', '--annotation_columns',
        action='store',
        dest='annot',
        required=True,
        help='Comma separated obs columns in h5 AnnData containing annotations that CELLEX should consider\
             e.g. major_cell_types,subtype,t_cell_subcluster'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='sc_CELLEX_data',
        help='Basename of output files.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-vb', '--verbose',
        action='store',
        dest='vb',
        default='True',
        help='Specify if the output should be printed.'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.

    # Load the AnnData file
    if options.vb:
        print('Reading in AnnData')
    adata = sc.read_h5ad(filename=options.h5)

    # Split the column names into a list
    annot_cols = options.annot.split(',')

    start_time = time.time()
    _ = run_CELLEX(
        adata,
        annotation_columns=annot_cols,
        output_file=options.of,
        verbose=options.vb
    )
    execution_summary = "Analysis execution time [{}]:\t{}".format(
        "run_CELLEX.py",
        str(timedelta(seconds=time.time()-start_time))
    )

    if (options.vb == True):
        print(execution_summary)


if __name__ == '__main__':
    main()