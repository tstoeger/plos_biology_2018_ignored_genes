import os

import pandas as pd

from access_science_shared import inout


def gene_info(taxon_id=None, usecols=None):
    """
    Obtains gene info for individual genes or taxa

    Input:
        taxon_id    NCBI taxon ID (integer), or 'all'
        usecols     Columns to load

    Output:
        gene_info   df, Dataframe with NIH gene Info

    """

    # Input specific loading functions
    def load_all_taxa(usecols=usecols):
        p = inout.get_path(
            'geisen',
            'ncbi/gene_info/gene_info_full.gz')

        if os.path.exists(p):
            df = pd.read_csv(p, usecols=usecols)
        else:
            df = None   # Extra Safety
            raise EnvironmentError(
                'Did not find gene info for all taxa')
        return df

    def load_taxon(taxon_id, usecols):
        p = inout.get_path(
            'geisen',
            'ncbi/gene_info/gene_info_taxon_{}.gz'.format(taxon_id))

        if os.path.exists(p):
            df = pd.read_csv(p, usecols=usecols)
        else:
            df = None   # Extra Safety
            raise EnvironmentError(
                'Did not find gene info for taxon {}'.format(taxon_id))
        return df

    # Implement input specific behavior of gene_info
    if taxon_id == 'all':
        df = load_all_taxa(usecols)
    elif isinstance(taxon_id, int):
        df = load_taxon(taxon_id, usecols)
    else:
        raise EnvironmentError(
            'Did not recognize format of taxon_id')

    return df


def taxon_name(taxon_id):
    """
    Obtains name of taxon

    Note that the original reference data would also
    contain synonymous names, and the class of the name
    (e.g.: whether it is a scientific name)

    Input:
        taxon_id    int; optionally set to 'all' to load
                        table with full taxonomy information

    Output:
        taxon_name  str
    """

    p = inout.get_path(
        'geisen',
        'ncbi/taxon_names.h5')

    if taxon_id == 'all':

        df = pd.read_hdf(p, 'table')
        df = df.set_index('taxon_ncbi', verify_integrity=True)
        output = df

    else:

        q = 'taxon_ncbi=={}'.format(taxon_id)
        df = pd.read_hdf(p, 'table', where=q)

        v = df['taxon_name'].values
        if len(v) == 0:
            print(
                'Could not find name of taxon {}'.format(
                    taxon_id))
            name = 'not found taxon; id is {}'.format(taxon_id)
        elif len(v) > 1:
            raise EnvironmentError(
                'Found multipole records for taxon {}'.format(
                    taxon_id))
        elif len(v) == 1:
            name = str(v[0])
        else:
            raise EnvironmentError(
                'Some error in code. This condition should never be' +
                'triggered. Please investigate code of this function.')

        output = name

    return output
