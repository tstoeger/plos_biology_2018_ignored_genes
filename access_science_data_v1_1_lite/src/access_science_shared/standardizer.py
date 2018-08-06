from access_biology_data import meta
from access_literature_data import medline


def filter_by_paper_kind(df, paper_kind):
    """
    Filters a dataframe for pubmed_id that also
    appear in paper_kind, which allows filtering of medline
    by potentially diverse options; prsently supported:
    'research'  --> research artciles of diverse kinds


    Input:
        df          dataframe, containing 'pubmed_id' in column
        paper_kind  input option for get_pubmed_ids_passing_medline_filter
                        e.g.: research; will filter according to
                        MedLine
    """

    if 'pubmed_id' not in df.columns:
        raise EnvironmentError(
            'Failed to find pubmed_id in provided input df')
    else:
        r = medline.select_medline_records(
            columns_sql='''
                medline.pubmed_id''',
            taxon_id=None,
            kind='research')
        r = r['pubmed_id'].values

        df = df[df.loc[:, 'pubmed_id'].isin(r)]
    return df


def reference_genes(taxon_id, ref_code):
    """
    Obtains a list of reference genes

    Input:
        taxon_id    int
        ref_code    str;  if it contains
                        l   -> at least one medline paper
                        o   -> official nomenclature require
                        p   -> protein-coding only

    Output:
        ref_genes   sorted list of gene identifiers
    """

    df = meta.gene_info(taxon_id)

    if df.shape[0] == 0:
        raise EnvironmentError("""
            Did not find gene info for taxon {}""".format(
            int(taxon_id)))

    if 'l' in ref_code:
        genes_in_medline = medline.gene2pubmed(taxon_id, ['gene_ncbi'])
        f = df.loc[:, 'gene_ncbi'].isin(genes_in_medline['gene_ncbi'])
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes with at least one paper,
                no gene is left.""")

    if 'o' in ref_code:  # official nomeclature
        f = df.loc[:, 'Nomenclature_status'] == 'O'
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes with official nomeclature,
                no gene is left.""")

    if 'p' in ref_code:  # protein-coding
        f = df.loc[:, 'type_of_gene'] == 'protein-coding'
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for protein-coding, no gene is
                left.""")

    if 'r' in ref_code:
        genes_in_medline = medline.gene2pubmed(
            taxon_id, ['pubmed_id', 'gene_ncbi'], paper_kind='research')
        f = df.loc[:, 'gene_ncbi'].isin(genes_in_medline['gene_ncbi'])
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes with at least one research paper,
                no gene is left.""")

    ref_genes = sorted(df.loc[:, 'gene_ncbi'].values)

    return ref_genes
