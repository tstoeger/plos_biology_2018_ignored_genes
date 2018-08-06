import pandas as pd

from access_literature_data import medline


def get_year_of_discovery(taxon_id, ref_genes):
    """
    Returns earliest years within research papers covered
    in MedLine

    Input:
        taxon_id
        ref_genes

    Output:
        genes_earliest_years    df with first_year and first_solo_year

    """

    ref_gene2pubmed = medline.gene2pubmed(
        taxon_id, paper_kind='research', ref_genes=ref_genes)

    df_m = medline.select_medline_records(
        columns_sql='''
                medline.pubmed_id,
                medline.pubdate_year''',
        taxon_id=taxon_id,
        kind='research')

    df_m = df_m[df_m['pubmed_id'].isin(ref_gene2pubmed['pubmed_id'])]
    columns_to_use = ['pubmed_id', 'pubdate_year']
    df_m = df_m.loc[:, columns_to_use].drop_duplicates()

    genes_per_paper = ref_gene2pubmed['pubmed_id'].value_counts(
    ).to_frame('genes')
    df_m = pd.merge(df_m, genes_per_paper, left_on='pubmed_id',
                    right_index=True, how='inner')
    df_m.loc[:, 'taxon_ncbi'] = taxon_id

    # add genes to medline
    master = pd.merge(
        df_m,
        medline.gene2pubmed(
            taxon_id=taxon_id,
            paper_kind='research',
            ref_genes=ref_genes),
        left_on=['taxon_ncbi', 'pubmed_id'],
        right_on=['taxon_ncbi', 'pubmed_id'],
        how='inner').drop_duplicates()

    # get initial years
    is_single_gene_paper = master['genes'] == 1
    genes_earliest_years = pd.merge(
        master.loc[
            :,
            ['gene_ncbi', 'pubdate_year']].groupby(
                'gene_ncbi').agg(min).reset_index().rename(
                    columns={'pubdate_year': 'first_year'}),
        master.loc[
            is_single_gene_paper,
            ['gene_ncbi', 'pubdate_year']].groupby(
                'gene_ncbi').agg(min).reset_index().rename(
                    columns={'pubdate_year': 'first_solo_year'}),
        left_on='gene_ncbi',
        right_on='gene_ncbi',
        how='outer'
    )

    f = master['genes'] == 1
    genes_earliest_years = pd.merge(
        master.loc[:, ['gene_ncbi', 'pubdate_year']].groupby(
            'gene_ncbi').agg(
            min).rename(columns={'pubdate_year': 'first_year'}),
        master.loc[f, ['gene_ncbi', 'pubdate_year']].groupby(
            'gene_ncbi').agg(
            min).rename(columns={'pubdate_year': 'first_solo_year'}),
        left_index=True,
        right_index=True,
        how='outer'
    )

    genes_earliest_years = genes_earliest_years.loc[ref_genes, :]

    return genes_earliest_years
