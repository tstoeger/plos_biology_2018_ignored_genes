from access_literature_data import medline

import re

import numpy as np
import pandas as pd


def count_citations(
        df_query, earliest_year_after_publication=1, years_to_include=3):

    """
    Counts number of citations following year of publication. If the last
    year that shall be included, exceeds the last complete year
    (namely 2015), citations will not be computed (NaN will be returned)

    Input:
        df  DataFrame, containing 'pubate_year' and 'citation_year_'s
        earliest_year_after_publication     int First year after publication
                                            that will be included in counting
                                            citations
        years_to_include    int  Amount of years that shall be included
                                    when counting citations (first year will
                                    be the one indicated by
                                    earliest_year_after_publication)

    Output:
        df_citations_over_span  DataFrame with Citations over requested
                                    time span
    """

    # It appears that 2016 might isn't fully indexed in WoS yet
    last_complete_year = 2015

    lead_pattern = 'citation_year_'
    f = [x.startswith(lead_pattern) for x in df_query.columns]

    df_citations = df_query.loc[:, f]
    df_pubdate_year = df_query.loc[:, ['pubdate_year']]

    df_citations.columns = [
        int(re.findall('{}(.*)$'.format(lead_pattern), x)[0]
            ) for x in df_citations.columns]

    df_o = pd.DataFrame(index=df_citations.index, columns=['citations'])

    # all_available_years = set(df_citations.columns)

    for y in df_citations.columns:
        earliest_year = y + earliest_year_after_publication
        last_year = earliest_year + years_to_include

        f = df_pubdate_year['pubdate_year'] == y

        ideal_years_to_query = set(np.arange(earliest_year, last_year))
        # years_to_query = ideal_years_to_query.intersection(all_available_years)
        # years_to_query = sorted(list(years_to_query))
        years_to_query = sorted(list(ideal_years_to_query))

        if max(ideal_years_to_query) <= last_complete_year:
            # if max(years_to_query) <= df_pubdate_year['pubdate_year'].max():
            #     if max(years_to_query) <= last_complete_year:
            dff = df_citations.loc[f, years_to_query].sum(axis=1)
            df_o.loc[dff.index, 'citations'] = dff.values
            #     else:
            #         df_o.loc[f, 'citations'] = np.nan
            # else:
            #     df_o.loc[f, 'citations'] = np.nan
        else:
            df_o.loc[f, 'citations'] = np.nan

    df_o['citations'] = df_o['citations'].astype(float)

    return df_o


def load_medline_wos_of_taxon(taxon_id, years_to_include, kind, unambiguous):

    df_m = medline.select_medline_wos_records(
        columns_sql='''
                medline.pubmed_id,
                medline.pubdate_year,
                medline.amount_of_authors,
                medline.j_name_s,
                pubmed2taxon.taxon_ncbi''',
        years_range='all',
        taxon_id=taxon_id,
        kind=kind,
        unambiguous=unambiguous)

    df_c = count_citations(
        df_m,
        earliest_year_after_publication=1,
        years_to_include=years_to_include)

    df = pd.merge(
        df_m[['pubmed_id', 'pubdate_year', 'amount_of_authors']],
        df_c,
        left_index=True,
        right_index=True,
        how='inner')

    df = df.dropna()

    return df


def add_citations(df_m, years_to_include):
    """
    Similar to the routine included in load_medline_wos_of_taxon,
    however there is no filtering for columns
    """

    df_c = count_citations(
        df_m,
        earliest_year_after_publication=1,
        years_to_include=years_to_include)

    df = pd.merge(
        df_m,
        df_c,
        left_index=True,
        right_index=True,
        how='inner')

    df = df.dropna()

    return df


def filter_for_papers_with_reference_genes(taxon_id, df, reference_genes):

    gene2pubmed = medline.gene2pubmed(taxon_id, ['gene_ncbi', 'pubmed_id'])
    gene2pubmed = gene2pubmed[gene2pubmed['gene_ncbi'].isin(reference_genes)]

    f = df['pubmed_id'].isin(gene2pubmed['pubmed_id'])
    df = df.loc[f, :]

    return df


def add_yearly_citation_rank(df):

    df['yearly_citation_rank'] = df[
        ['pubdate_year', 'citations']].groupby(['pubdate_year']).rank(pct=True)

    return df


def add_team_scale(df):

    f = df['amount_of_authors'] == 1
    if any(f):
        df.loc[f, 'team_scale'] = 'single'

    f = df['amount_of_authors'] == 2
    if any(f):
        df.loc[f, 'team_scale'] = 'pair'

    f = df['amount_of_authors'] > 2
    if any(f):
        df.loc[f, 'team_scale'] = 'team'

    f = df['amount_of_authors'] == -1  # undefined, in practice: consortia
    if any(f):
        df.loc[f, 'team_scale'] = 'team'

    return df


def load_shared_gene2pubmed(taxon_id, reference_df, reference_genes):

    gene2pubmed = medline.gene2pubmed(taxon_id, ['gene_ncbi', 'pubmed_id'])
    gene2pubmed = gene2pubmed[gene2pubmed['gene_ncbi'].isin(reference_genes)]
    f = gene2pubmed['pubmed_id'].isin(reference_df['pubmed_id'])
    gene2pubmed = gene2pubmed.loc[f, :]

    return gene2pubmed
