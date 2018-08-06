import os

import pandas as pd

from access_science_shared import inout, utils, mapper


def external_links():
    p = inout.get_path(
        'drugbank',
        'external_drug_links/drug_links.csv')
    df = pd.read_csv(p)

    df = df.rename(
        columns={'DrugBank ID': 'drug_drugbank'}).set_index(
        'drug_drugbank', verify_integrity=True)

    if df['Name'].str.lower().value_counts().max() > 1:
        raise ValueError('Drug names are not unambigous')

    return df


def genes_2_drugs_and_status(taxon_id, target_class):
    """
    Loads drugbank IDs, and status for individual genes of
    taxon_id

    Input:
        taxon_id        int  taxon ID
        target_class    str  'pharmacologically_active' or 'all'

    """

    # Define code to filter drugbank
    dictionary_of_taxa = {
        9606: 'Human'
    }

    # neglect separation by drug class (e.g: small molecule)
    considered_status = [
        'approved',
        'experimental',
        'illicit',
        'investigational',
        'nutrazeutical',
        'withdrawn'
    ]

    agg = []
    for status in considered_status:

        p = inout.get_path(
            'drugbank',
            'protein_identifiers/drug_target_identifiers/{}/{}.csv'.format(
                status, target_class))
        df = pd.read_csv(p)

        f = df['Species'] == dictionary_of_taxa[taxon_id]
        df = df.loc[f, ['UniProt ID', 'Drug IDs']]
        df.loc[:, 'status'] = status
        agg.append(df)

    df = pd.concat(agg, axis=0)
    df = utils.split_text_to_multiple_rows(
        df, 'Drug IDs', ';')
    df = df.rename(
        columns={
            'UniProt ID': 'protein_uniprot',
            'Drug IDs': 'drug_drugbank'}
    ).drop_duplicates()

    p_mapper = mapper._get_geisen_path('uniprot/uniprot_id_mapper.h5')
    if not os.path.exists(p_mapper):
        raise EnvironmentError(
            'uniprot_protein_2_gene_ncbi() requires uniprot_id_mapper')
    ma = pd.read_hdf(
        p_mapper,
        'table',
        columns=['protein_uniprot', 'gene_ncbi'],
        where='taxon_ncbi={}'.format(taxon_id))

    df = pd.merge(df, ma)[
        ['gene_ncbi', 'drug_drugbank', 'status']].drop_duplicates()

    df = df.sort_values(
        ['gene_ncbi', 'status', 'drug_drugbank']).reset_index(drop=True)

    return df
