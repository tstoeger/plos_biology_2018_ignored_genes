import re

from natsort import natsorted
import numpy as np
import pandas as pd

from access_mixed_data import genealacart

"""
A series of functions that mimic natural langauge and aim
to faciliate the formatting of feature sets for predictions
"""


def index_by_gene_ncbi(df):
    df = df.set_index('gene_ncbi', verify_integrity=True)
    return df


def filter_by_reference_index(df, ref_index):
    f = df.index.isin(ref_index)
    df = df.loc[f, :]
    return df


def filter_for_genealacart(df):
    ei = genealacart.load_genealacart_dataset(
        'ExternalIdentifiers')
    entrez_in_genealacart = ei['EntrezGene_x'].unique()
    f = df.index.isin(entrez_in_genealacart)
    df = df.loc[f, :]
    return df


def count_nonzero_in_row(df, name):
    h = df.apply(lambda x: np.count_nonzero(x), axis=1)
    h = h.to_frame(name)
    return h


def concatenate_missing_to_left(df_left, df_right):
    f = ~df_right.index.isin(df_left.index)
    r = df_right.loc[f, :]
    df = pd.concat([df_left, r], join='outer', verify_integrity=True)
    return df


def enforce_reference_index(df, ref_index):
    df = df.loc[ref_index, :]
    return df


def mask_swiss_trembl_specificity_label(df):
    df.columns = [
        re.sub('swissprot|trembl', 'swiss_or_trembl', x
               ) for x in df.columns]
    return df


def prefix_column_names(df, prefix):
    df.columns = ['{}__{}'.format(prefix, x) for x in df.columns]
    return df


def pivot_gene_ncbi_to_predictor(df, col_name):
    df.loc[:, 'is_present'] = True
    df = df.pivot(
        index='gene_ncbi',
        columns=col_name,
        values='is_present')
    df = df.fillna(False)
    return df


def remove_invariant_columns(df):
    f = df.apply(lambda x: len(set(x)) > 1)
    df = df.loc[:, f]
    return df


def replace_less_than_one_by_nan_and_log_transform(df):
    f = df < 1
    df[f] = np.nan
    df = df.apply(lambda x: np.log10(x))
    return df


def replace_less_than_0p1_by_nan_and_log_transform(df):
    f = df < 0.1
    df[f] = np.nan
    df = df.apply(lambda x: np.log10(x))
    return df


def replace_nan_by_false(df):
    df = df.fillna(False)
    return df


def replace_nan_by_minus_one(df):
    df = df.fillna(-1)
    return df


def replace_nan_by_minus_two(df):
    df = df.fillna(-2)
    return df


def replace_nan_by_zero(df):
    df = df.fillna(0)
    return df


def replace_nan_by_plus_one(df):
    df = df.fillna(+1)
    return df


def replace_natsorted_first_column_by_integers(df):
    u = df.iloc[:, 0].unique()
    u = natsorted(u)
    y = dict()
    for it, orig in enumerate(u):
        y[orig] = it
    t = df.iloc[:, 0].apply(lambda x: y[x])
    df.iloc[:, 0] = t
    return df


def require_at_least_ten_non_zero_records(df):
    c = (df > 0).sum()
    f = c >= 10
    df = df.loc[:, f]
    return df


def require_at_least_twenty_non_zero_records(df):
    c = (df > 0).sum()
    f = c >= 10
    df = df.loc[:, f]
    return df


def require_at_least_fivehundred_non_zero_records(df):
    c = (df > 0).sum()
    f = c >= 500
    df = df.loc[:, f]
    return df


def swiss_complemented_with_trembl(df_swiss, df_trembl):
    df_left = df_swiss
    df_right = df_trembl
    functions = [index_by_gene_ncbi, mask_swiss_trembl_specificity_label]
    df_left = _run_functions_sequentially(df_left, functions)
    df_right = _run_functions_sequentially(df_right, functions)
    df = concatenate_missing_to_left(df_left, df_right)
    return df


def _run_functions_sequentially(df, functions):
    for func in functions:
        df = func(df)
    return df
