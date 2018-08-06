import pandas as pd

from access_literature_data import medline

import json

import resci_inout as inout


def gene_pmid_research_gene_pmid_research():

    p = inout.get_internal_path(
        (
            '171106f_citation_network/research_unambiguous/'
            'ut_citation_network.json')
    )

    with open(p) as data_file:
        references_from_disk = json.load(data_file)

    manager = medline.select_medline_wos_records(
        columns_sql='''
            medline.pubmed_id,
            ut2pmid.ut AS wos_id''',
        taxon_id='all',
        kind='research',
        unambiguous=True
    )

    allowed_wos_id = set(manager['wos_id'])

    references_from_disk = {k: set(v) for k, v in references_from_disk.items()}

    # only consider those publications defined in manager
    references_from_disk = {k: v.intersection(
        allowed_wos_id) for k, v in references_from_disk.items()}

    references_from_disk = {k: list(v)
                            for k, v in references_from_disk.items()}

    agg_values = []
    agg_keys = []
    for k, v in references_from_disk.items():
        if any(v):
            agg_values.append(v)
            agg_keys.append([k] * len(v))

    def flatten_list(x):
        return [item for sublist in x for item in sublist]

    agg_values = flatten_list(agg_values)
    agg_keys = flatten_list(agg_keys)

    df = pd.DataFrame(
        index=agg_keys,
        data=agg_values,
        columns=['citing']).reset_index().rename(
        columns={'index': 'cited'})

    def update_with_pmid(df, field):
        df = pd.merge(df, manager, left_on=field,
                      right_on='wos_id', how='left')
        df = df.drop(field, 1)
        df = df.drop('wos_id', 1)
        df = df.rename(columns={'pubmed_id': field})
        return df

    df = update_with_pmid(df, 'cited')
    df = update_with_pmid(df, 'citing')

    df = df.drop_duplicates().reset_index(drop=True)

    return df
