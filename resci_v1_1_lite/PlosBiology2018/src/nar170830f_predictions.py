import glob
import os

import numpy as np
import pandas as pd

import resci_inout as inout


def make_base(project_base, sub_name, ref_genes, u_features, df_target):
    """
    Make prediction bases; Note that it will remove
    rows with NaN as a feature, or as prediction target:

    Input:
        project_base    str, name of series of analysis
        sub_name        str, will be used as folder
                        (should be specific to dataset)
        ref_genes       list-like; reference genes that shall be filtered for
        u_features      either a dataframe with features or a dictionary
                        of many dataframes (note: always: index is gene_ncbi)
        df_target       dataframe with target to predict as sole column


    """

    p = inout.get_internal_path(os.path.join(
        project_base,
        sub_name,
        'input'
    ))

    if df_target.shape[1] != 1:
        raise EnvironmentError(
            'df_target must have exactly one column')
    else:
        df_target.columns = ['prediction_target']

    if isinstance(u_features, dict):
        df = pd.concat(u_features, join='outer', axis=1)
    else:
        df = u_features

    if isinstance(df.columns, pd.core.index.MultiIndex):
        df.columns = df.columns.droplevel(level=0)

    df = df.loc[ref_genes, :].dropna()

    master = pd.merge(
        df,
        df_target,
        left_index=True,
        right_index=True,
        how='inner'
    )

    master = master.dropna()

    if any(master.index.duplicated()):
        raise EnvironmentError(
            'At least one index is duplicated')

    features = master.drop('prediction_target', axis=1)
    target = master.loc[:, ['prediction_target']]
    _export(features, p, 'features')
    _export(target, p, 'target')


def load_predicitions(path_to_predictions, min_required=0, max_to_load=np.inf):
    """
    Aggregates predictions
    """

    def pooling_fun(x):
        return np.nanmedian(x)

    if path_to_predictions.endswith('/'):
        path_to_predictions = path_to_predictions[:-1]

    if path_to_predictions.endswith('\\'):
        path_to_predictions = path_to_predictions[:-1]

    glob_query = os.path.join(
        path_to_predictions,
        'prediction*')

    prediction_files = glob.glob(glob_query)

    if len(prediction_files) == 0:
        raise EnvironmentError(
            'Did not find a file matching {} '.format(
                glob_query))

    agg = list()
    loaded = 0
    for p in prediction_files:

        try:
            df = pd.read_csv(p)
            if df.shape[0] > 0:
                if loaded < max_to_load:
                    process = True
                    loaded = loaded + 1
                else:
                    process = False
            else:
                process = False
        except:
            process = False

        if process:   # some files are empty, likely related to cluster
            df = df.set_index('gene_ncbi', verify_integrity=False)
            f = df.index.duplicated()
            df = df.loc[~f, :]

            agg.append(df)
    dfj = pd.concat(agg, axis=1, join='outer')

    project_base, _ = os.path.split(path_to_predictions)
    p_target = os.path.join(project_base, 'input', 'target.csv.gz')
    dfo = pd.read_csv(p_target).set_index('gene_ncbi')

    dfj_pooled = dfj.apply(pooling_fun, axis=1)
    # sometimes new identifier are created, liekly related to cluster
    together = pd.concat([dfo, dfj_pooled], axis=1, join='inner')

    together = together.rename(columns={
        'prediction_target': 'target',
        0: 'predicted'
    })

    if loaded < min_required:
        print('{} only contains {} files'.format(
            path_to_predictions,
            loaded))

    return together


def pool_target_and_predicitions(path_to_predictions):
    """
    Creates table with all preditions and saves in prediction
    folder als pooled_target_and_predition.csv.gz
    """

    if path_to_predictions.endswith('/'):
        path_to_predictions = path_to_predictions[:-1]

    if path_to_predictions.endswith('\\'):
        path_to_predictions = path_to_predictions[:-1]

    glob_query = os.path.join(
        path_to_predictions,
        'prediction*')

    prediction_files = glob.glob(glob_query)

    if len(prediction_files) == 0:
        raise EnvironmentError(
            'Did not find a file matching {} '.format(
                glob_query))

    agg = list()
    for p in prediction_files:

        try:
            df = pd.read_csv(p)
            if df.shape[0] > 0:
                process = True
            else:
                process = False
        except:
            process = False

        if process:   # some files are empty, likely related to cluster
            df = df.set_index('gene_ncbi', verify_integrity=False)
            f = df.index.duplicated()
            df = df.loc[~f, :]

            agg.append(df)
    dfj = pd.concat(agg, axis=1, join='outer')
    dfj.columns = ['prediction_{}'.format(x) for x in np.arange(
        int(dfj.shape[1]))]

    project_base, _ = os.path.split(path_to_predictions)
    p_target = os.path.join(project_base, 'input', 'target.csv.gz')
    dfo = pd.read_csv(p_target).set_index('gene_ncbi')
    dfo.columns = ['target']

    # sometimes new identifier are created, liekly related to cluster
    together = pd.concat([dfo, dfj], axis=1, join='inner')

    p = os.path.join(
        path_to_predictions,
        'pooled_target_and_prediciton.csv.gz')
    together.to_csv(p, compression='gzip', index=True)


def _export(df, p, name):
    pp = os.path.join(p, '{}.csv.gz'.format(name))
    inout.ensure_presence_of_directory(pp)
    if os.path.exists(pp):
        raise EnvironmentError('File already exists: {}'.format(p))
    else:
        df.to_csv(pp, compression='gzip')
