import glob
import os

import pandas as pd

from access_science_shared import inout

# This module contains tools for extracting specific datasets from
# genealacart (which is anticipated in manual data)

# set version(s) of genealacart that shall be used
version_170217 = 'genealacart/170217/GeneALaCart-431-170217*'
reference_version = version_170217


def load_genealacart_dataset(dataset, version=reference_version):
    """
    - Loads genealacart dataset
    - removes ambiguous mappings between
        symbols and entrez gene (Note: which in the 170217 version of
        genelacart only corresponds to the remval of two genes)
    - introdcue gene_ncbi for common mapping to rest of science of
        biology

    Input:

        dataset         str, name of dataset, e.g.: 'Enhancers' note that a
                        given genealacart datast can contain several smaller
                        ones; corresponds to names of files within the zip'ed
                        download of genealacarrt
        version         str, glob-able pattern to folders of indivdual batches
                        within a geisen_manual data path
                        By default: set to latest reviewed version
    """

    df_identifiers = _load_batches(version, 'ExternalIdentifiers')
    allowed_symbols = _get_unambiguous_symbols(df_identifiers)

    df = _load_batches(version, dataset)

    f = df_identifiers['Symbol'].isin(allowed_symbols)
    df_identifiers = df_identifiers.loc[f, :]

    df = pd.merge(
        df_identifiers[['Symbol', 'EntrezGene']],
        df,
        left_on='Symbol',
        right_on='Symbol',
        how='left')

    # Enforce standard scinece of sciene terminology for main ID
    df = df.rename(columns={
        'EntrezGene': 'gene_ncbi'
    })

    # Safety given heterogeous dynamic input data sets
    df = df.drop_duplicates()

    return df


"""
    GENERAL SUPPORT FUNCTIONS

        _load_batches
        _get_unambiguous_symbols
"""


def _load_batches(path_pattern, dataset):
        p_scheme = inout.get_path(
            'geisen_manual',
            os.path.join(path_pattern))
        agg = []
        for p in glob.glob(p_scheme):
            p = os.path.join(p, '{}.txt'.format(dataset))
            agg.append(pd.read_table(p))
        df = pd.concat(agg)
        df = df.drop_duplicates()

        # place symbols lower-case as genealacart does not appear
        # to distinguish internally
        n = ['InputTerm', 'Symbol']
        for t in n:
            df.loc[:, t] = df.loc[:, t].str.lower()

        return df


def _get_unambiguous_symbols(df_identifiers):
        """
        Obtains the entries which have an unambiguous
        Input Term, Symbol and Entrez Gene. Will return
        unambigous Symols (which are used internally by
        GeneALaCart)

        Note that input symbols can be duplicate, if the same
        input maps to several Symbols in Genealacart

        Input:
            df_identifiers  DataFrame with identifiers containing
                            ['InputTerm', 'Symbol', 'EntrezGene']
                            e.g.: df_identifiers = _load_batches(
                                version, 'ExternalIdentifiers')

        Output:
            unambiguous_symbols  list; all symbols that can
                                    be mapped unambiguously to
                                    gene_ncbi within genealacart
        """

        df = df_identifiers

        f = df['EntrezGene'].notnull()
        dff = df.loc[f, :]

        terms_that_should_be_unique = ['InputTerm', 'Symbol', 'EntrezGene']

        dff = dff.loc[:, terms_that_should_be_unique]
        dff = dff.drop_duplicates()

        ambiguous = dict()
        for t in terms_that_should_be_unique:
            c = dff.loc[:, t].value_counts()
            c = list(c[c > 1].index)
            ambiguous[t] = c

        for k, v in ambiguous.items():
            num_ambiguous = len(v)
            if num_ambiguous > 0:
                print(k, 'has', num_ambiguous, 'ambiguous entries')
                f = dff.loc[:, k].isin(v)
                dff = dff.loc[~f, :]

        unambiguous_symbols = dff.loc[:, 'Symbol'].unique()

        return unambiguous_symbols



