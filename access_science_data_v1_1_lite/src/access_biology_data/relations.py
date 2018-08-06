import numpy as np
import pandas as pd

from access_science_shared import inout


def biogrid(taxon_id):
    """
    Imports BioGrid database

    Input:
        taxon_id    int

    Output:
        df_biogrid
    """

    p = inout.get_path(
        'geisen',
        'biogrid/biogrid_homologous_taxon_{}.csv.gz'.format(
            int(taxon_id)))
    df = pd.read_csv(p)
    return df


def bioplex2():
    """
    Imports the bioplex2 dataset as published by Huttlin et. al,
    2017, Nature

    """

    p = inout.get_path(
        'publications',
        'huttlin2017/nature22366-s7.xlsx'
    )

    c = ['bioplex2_id', 'gene_ncbi']
    df = pd.read_excel(p).rename(
        columns={
            'Cluster Number': 'bioplex2_id',
            'GeneID': 'gene_ncbi',
            'Symbol': 'symbol_ambiguous'
        })[c].drop_duplicates().sort_values(c).reset_index(drop=True)

    return df


def rolland_trans_interactions(complement=False):
    """
    Interactions from Rolland et al. (for human); Will remove
    self-interactions; note that the interactions are
    bi-directional and that none both must be trated equally

    complement:     default: False will no add complement
                    (intearctors a and b swapped)

    """

    p = inout.get_path(
        'geisen',
        'papers/rolland_2014/rolland_table_binary_interactions.csv.gz'
    )

    df = pd.read_csv(p)
    df = df[['entrez_gene_ida', 'entrez_gene_idb']]
    df.columns = ['a', 'b']

    if complement:
        df = pd.concat(
            [df, df.rename(
                columns={'a': 'b', 'b': 'a'})], axis=0, ignore_index=True, )

    df = df.rename(
        columns={'a': 'entrez_gene_ida',
                 'b': 'entrez_gene_idb'})

    f = df['entrez_gene_ida'] == df['entrez_gene_idb']
    df = df.loc[~f, :]
    df = df.drop_duplicates()

    return df


def homologene():
    """
    Homologene

    Output:
        dataframe   Linkage table between homologous genes of
                    ~20 model taxa
    """

    p = inout.get_path(
        'geisen',
        'ncbi/homologene.gz')
    df = pd.read_csv(p)

    df['taxon_ncbi'] = df['taxon_ncbi'].astype(float)
    df['taxon_ncbi'] = df['taxon_ncbi'].astype(int)

    return df


def relate_to_homologene(
        taxon_id,
        df_orig,
        extend_column_names=True,
        aggregate=True):
    """
    Relates dataframe with data from one organism
    to data of another organism (taxon_id)

    Input:
        taxon_id    int  taxon_id of interest (new taxon)
        df_orig     dataframe  original taxon
        extend_column_names     If True the original taxon will
                                be indicated in column names.
                                If False, column names remain
                                will be same as in df_orig
        aggregate   default: True will aggregate over multipel homologss
                    if False: will remove duplicates.

    Output:
        df_related  dataframe; data mapped to new taxon

    """

    df_h = homologene()

    # values of cis taxon (taxon of original data)
    df_c = pd.merge(
        df_h,
        df_orig,
        left_on='gene_ncbi',
        right_on='gene_ncbi',
        how='inner')

    taxa = df_c['taxon_ncbi'].unique()
    if len(taxa) > 1:
        raise EnvironmentError('Input data set contained muliple taxa.')
    elif len(taxa) == 0:
        raise EnvironmentError('Input data contained no mappable taxon.')
    else:
        cis_taxon = taxa[0]

    df_c = df_c.drop(['gene_ncbi'], axis=1)
    df_c = df_c.drop(['taxon_ncbi'], axis=1)

    # values mapped to trans traxon (taxon of interest)
    df_t = df_h[df_h['taxon_ncbi'] == taxon_id].drop(['taxon_ncbi'], axis=1)
    df_t = pd.merge(
        df_t,
        df_c,
        left_on='homologene_group',
        right_on='homologene_group',
        how='inner')

    df_t = df_t.drop('homologene_group', axis=1)
    if aggregate:
        df_t = df_t.groupby('gene_ncbi').agg(np.median)
    else:
        df_t = df_t.drop_duplicates().set_index('gene_ncbi')

    if extend_column_names:
        df_t.columns = ['{} (via taxon {})'.format(
            x, cis_taxon) for x in df_t.columns]
    df_related = df_t.reset_index()

    return df_related
