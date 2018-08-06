import os

import numpy as np
import pandas as pd

from access_science_shared import inout


def biosystems(taxon_id=None, databases='all'):
    """
    Biosystems databases, which is an aggregation of several pathway
    databases (of which some may not be defined for every taxon)

    Source: NIH

    Input:
        taxon_id    int
        databases   default is 'all', which will load
                        'biocyc', 'kegg', 'lipid_maps',
                        'pathway_interaction_database',
                        'reactome', 'wikipathways'
                    otherwise: supply a list with pathways of interest
                        e.g.: ['kegg','reactome']

    Output:
        dataframe
    """

    if taxon_id is None:
        raise EnvironmentError('No taxon specified')

    def _load_biosystem_geisen_v1_1_upwards(p):
        df = pd.read_csv(p)
        return df

    # Load datasets
    if databases == 'all':
        databases = [
            'biocyc', 'kegg', 'lipid_maps', 'pathway_interaction_database',
            'reactome', 'wikipathways']
    elif isinstance(databases, str):
        databases = [databases]
    else:
        raise EnvironmentError('databases are not defined.')

    agg = list()
    for d in databases:
        d = d.lower()  # file names are all lower

        p = inout.get_path(
            'geisen',
            'ncbi/biosystems/' +
            'biosystem_{}/biosystem_{}_taxon_{}.csv.gz'.format(
                d, d, taxon_id))

        if os.path.exists(p):
            df = _load_biosystem_geisen_v1_1_upwards(p)

            p = inout.get_path(
                'geisen',
                'ncbi/biosystems/biosystem_' +
                '{}/biosystem_{}_accession_to_term.csv.gz'.format(
                    d, d))
            df_name = pd.read_csv(p)

            # patch with names
            df = pd.merge(
                df, df_name,
                left_on='accession',
                right_on='accession',
                how='left')
        else:
            print('Biosystem {} is not defined for taxon {}'.format(
                d, taxon_id))

        agg.append(df)

    df = pd.concat(agg)
    return df


def disease_genealacart(taxon_id=9606, add_absenece=True):
    """
    Unified diseases
    Source: Genealacart

    Input:
        taxon_id    int
        add_absence bool; default is True; add genes for which there
                    is no disease entry
    Output:
        dataframe
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'genealacart/genealacart_diseases_kind.gz')
    df = pd.read_csv(p)

    df = df.set_index('gene_ncbi')
    df.columns = ['unified_disease__{}'.format(x) for x in df.columns]
    df = df.rename(
        columns={'unified_disease__disease_name': 'unified_disease'})
    df = df.reset_index()

    if add_absenece:
        p = inout.get_path(
            'geisen',
            'genealacart/genealacart_diseases_amount.gz')
        df_c = pd.read_csv(p)

        extra = np.setdiff1d(
            np.array(df_c['gene_ncbi'].unique()),
            np.array(df['gene_ncbi'].unique()))

        df_p = pd.DataFrame(index=extra, columns=['unified_disease'])
        df_p.loc[:, 'unified_disease'] = 'No known disease'
        df_p.index.name = 'gene_ncbi'
        df_p = df_p.reset_index()
        df = pd.concat([df, df_p], axis=0)

    return df


def generif(taxon_id):
    """
    Loads gene RIFs

    Input:
        taxon_id    list of taxa, or int of taxon id, or 'all'

    Output:
        generifs    dataframe

    """

    p = inout.get_path(
        'geisen',
        'ncbi/generifs_basic.gz')

    df = pd.read_csv(p, low_memory=False)

    if isinstance(taxon_id, int):
        taxon_id = [taxon_id]

    if taxon_id != 'all':
        f = df.loc[:, 'taxon_ncbi'].isin(taxon_id)
        df = df.loc[f, :]

    return df


def go(taxon_id,
        category=None,
        negating_support=None,
        any_negating_support=None,
        temporary_evidence=None,
        unmapped_evidence=None):
    """
    Loads GO annotation
    Source: NIH (GO mapped to genes)

    Input:
        taxon_id   int
        category   list, values of ['Function', 'Process', 'Component']
        negating_support   List of Bool, for single term negating support
        any_negating_support   List of Bool, spreads negation
        temporary_evidence  List of Bool, temporary (12 months limit)
        unmapped_evidence   List of Bool, wheater evidence is unmapped

    Output:
        dataframe

    """

    # Initialize defaults
    if category is None:
        category = ['Function', 'Process', 'Component']
    if negating_support is None:
        negating_support = [True, False]
    if any_negating_support is None:
        any_negating_support = [True, False]
    if temporary_evidence is None:
        temporary_evidence = [True, False]
    if unmapped_evidence is None:
        unmapped_evidence = [True, False]

    def _load_go_v1(p, taxon_id):

        print("""
            Found gene2go.gz . This is indicative of usage of
            an outdated version of the geisen datasets.

            Using legacy fallback. Note that this option might
            become deactivated in a future release.
            """)

        df = pd.read_csv(p)

        if isinstance(taxon_id, int):
            taxon_id = [taxon_id]

        if taxon_id != 'all':
            f = df.loc[:, 'taxon_ncbi'].isin(taxon_id)
            df = df.loc[f, :]
        return df

    def _load_go_v1_1(taxon_id):

        p = inout.get_path(
            'geisen', 'ncbi/gene2go/gene2go_taxon_{}.csv.gz'.format(
                int(taxon_id)))
        df = pd.read_csv(p)

        p = inout.get_path(
            'geisen', 'ncbi/gene2go/go_id_to_term.csv.gz')
        df_label = pd.read_csv(p)

        df = pd.merge(
            df, df_label, left_on='GO_ID', right_on='GO_ID', how='left')
        return df

    # Load datasets
    p = inout.get_path('geisen', 'ncbi/gene2go.gz')
    if os.path.exists(p):
        df = _load_go_v1(p, taxon_id)
    else:
        df = _load_go_v1_1(taxon_id)

    # Filter datasets
    f = \
        (df.loc[:, 'Category'].isin(category)) & \
        (df.loc[:, 'Negating support'].isin(negating_support)) & \
        (df.loc[:, 'Any negating support'].isin(any_negating_support)) & \
        (df.loc[
            :,
            'Temporary Evidence (max. 12m)'].isin(temporary_evidence)) & \
        (df.loc[:, 'Unmapped Evidence'].isin(unmapped_evidence))

    df = df.loc[f, :]

    return df


def human_phenotype_genealacart(taxon_id=9606, add_absenece=True):
    """
    Human phenotypes
    Source: Human Phenotype Ontology through Genealacart

    Input:
        taxon_id    int
        add_absence bool; default is True; add genes for which there
                    is no phenotype entry
    Output:
        dataframe
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'genealacart/genealacart_phenotype_ontology_kind.gz')
    df = pd.read_csv(p)

    df = df.set_index('gene_ncbi')
    df = df.reset_index()

    if add_absenece:
        p = inout.get_path(
            'geisen',
            'genealacart/genealacart_phenotype_ontology_amount.gz')
        df_c = pd.read_csv(p)

        extra = np.setdiff1d(
            np.array(df_c['gene_ncbi'].unique()),
            np.array(df['gene_ncbi'].unique()))

        df_p = pd.DataFrame(
            index=extra,
            columns=[
                'human_phenotype_genealacart: human_phenotype_id',
                'human_phenotype_genealacart: human_phenotype_name'
            ])
        df_p.loc[
            :,
            'human_phenotype_genealacart: human_phenotype_name'] = \
            'No known human phenotype'
        df_p.index.name = 'gene_ncbi'
        df_p = df_p.reset_index()
        df = pd.concat([df, df_p], axis=0)

    return df


def interpro(taxon_id=None, databases='all', add_names=True, usecols=None):
    """
    Loads interpro from NIH. Note that for the same interpro there can be
    multiple entries (if there is support in multiple databases)

    Input:
        taxon_id    int, or 'all'
        databases   str, or list of strs, or all; options are:
                    cd, G3, MF, PD, PF, PI, PR, PS, PT, SF, SM, SS, TI
        add_names  default: True: add name of interpro domain
        usecols     optional, list of columns to be loaded

    Output:
        interpro    df

    """

    def load_interpro(taxon_id, database, usecols):

        p = inout.get_path(
            'geisen',
            'interpro/interpro_{}.h5'.format(database))

        if os.path.exists(p) is False:
            raise EnvironmentError(
                'Did not find interpro database {}'.format(
                    database))

        if taxon_id == 'all':
            df = pd.read_hdf(p, 'table')
        elif isinstance(taxon_id, int):
            q = 'taxon_ncbi=={}'.format(taxon_id)
            df = pd.read_hdf(p, 'table', where=q)
        else:
            raise EnvironmentError(
                'Did not recognize format of taxon_id')

        if usecols is not None:
            df = df.loc[:, usecols]

        return df

    #         DBCODE   |    name
    # ------------------------------------------------
    # cd          Conserved Domain
    # G3          CATH Superfamily
    # MF          Hamap
    # PD          ProDom
    # PF          Pfam
    # PI          Protein Information Resource
    # PR          PRINTS
    # PS          Prosite
    # PT          Panther
    # SF          Structure Function Linkage Database
    # SM          SMART
    # SS          SUPERFAMILY
    # TI          TIGR

    if databases == 'all':
        databases = [
            'cd', 'G3', 'MF', 'PD', 'PF', 'PI',
            'PR', 'PS', 'PT', 'SF', 'SM', 'SS',
            'TI']
    elif isinstance(databases, str):
        databases = [databases]

    agg = list()
    for d in databases:
        df = load_interpro(taxon_id, d, usecols)
        agg.append(df)
    df = pd.concat(agg)

    if add_names:
        p = inout.get_path(
            'geisen',
            'interpro/interpro_names.csv.gz')
        df_labels = pd.read_csv(p)

        df = pd.merge(
            df,
            df_labels,
            left_on='interpro_id',
            right_on='interpro_id',
            how='left')

    return df


def omim_genealacart(taxon_id=9606, add_absenece=True):
    """
    Mendelian diseases
    Source: OMIM through Genealacart

    Note that the source OMIM database, which is quite messy, would
    in principle allow a higher level of stratication then Genealacart

    Input:
        taxon_id    int
        add_absence bool; default is True; add genes for which there
                    is no disease entry
    Output:
        dataframe
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'genealacart/genealacart_omim_kind.gz')
    df = pd.read_csv(p)

    df = df.set_index('gene_ncbi')
    df.columns = ['omim_disease__{}'.format(x) for x in df.columns]
    df = df.rename(columns={'omim_disease__disease_name': 'omim_disease'})
    df = df.reset_index()

    if add_absenece:
        p = inout.get_path(
            'geisen',
            'genealacart/genealacart_omim_amount.gz')
        df_c = pd.read_csv(p)

        extra = np.setdiff1d(
            np.array(df_c['gene_ncbi'].unique()),
            np.array(df['gene_ncbi'].unique()))

        df_p = pd.DataFrame(index=extra, columns=['omim_disease'])
        df_p.loc[:, 'omim_disease'] = 'No entry in OMIM'
        df_p.index.name = 'gene_ncbi'
        df_p = df_p.reset_index()
        df = pd.concat([df, df_p], axis=0)

    return df


def orphanet_genealacart(taxon_id=9606, add_absenece=True):
    """
    Orphan diseases
    Source: Orphanet through Genealacart

    Input:
        taxon_id    int
        add_absence bool; default is True; add genes for which there
                    is no disease entry
    Output:
        dataframe
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'genealacart/genealacart_orphanet_kind.gz')
    df = pd.read_csv(p)

    df = df.set_index('gene_ncbi')
    df.columns = ['orphanet__disease__{}'.format(x) for x in df.columns]
    df = df.rename(
        columns={'orphanet__disease__disease_name': 'orphanet_disease'})
    df = df.reset_index()

    if add_absenece:
        p = inout.get_path(
            'geisen',
            'genealacart/genealacart_orphanet_amount.gz')
        df_c = pd.read_csv(p)

        extra = np.setdiff1d(
            np.array(df_c['gene_ncbi'].unique()),
            np.array(df['gene_ncbi'].unique()))

        df_p = pd.DataFrame(index=extra, columns=['orphanet_disease'])
        df_p.loc[:, 'orphanet_disease'] = 'No entry in Orphanet'
        df_p.index.name = 'gene_ncbi'
        df_p = df_p.reset_index()
        df = pd.concat([df, df_p], axis=0)

    return df
