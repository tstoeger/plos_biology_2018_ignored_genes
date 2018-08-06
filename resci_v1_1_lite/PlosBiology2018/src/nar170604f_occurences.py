import pandas as pd

from access_literature_data import medline


def count_papers_and_attention(ref_genes, ref_gene2pubmed):
    """
    Counts the number of articles in ref_gene2pubmed
    Uses the full gene2pubmed (not only the one of ref_gene2pubmed)
    to compute attention (thus also considerting genes that lie
    outisde of the present filtering)

    Input:
        ref_genes       list of reference genes
        ref_gene2pubmed dataframe with gene_ncbi and pubmed_id

    """

    ref_gene2pubmed = ref_gene2pubmed[
        ref_gene2pubmed['gene_ncbi'].isin(ref_genes)]

    full_gene2pubmed = medline.gene2pubmed(taxon_id='all')
    full_gene2pubmed = full_gene2pubmed[
        full_gene2pubmed['pubmed_id'].isin(ref_gene2pubmed['pubmed_id'])]

    h = full_gene2pubmed[['pubmed_id', 'gene_ncbi']].groupby('pubmed_id').agg(
        lambda x: 1 / len(x)).rename(
            columns={'gene_ncbi': 'attention'}).reset_index()

    master = pd.merge(ref_gene2pubmed, h, how='inner')

    papers = master['gene_ncbi'].value_counts().to_frame('papers')
    papers.index.name = 'gene_ncbi'
    papers.loc[:, 'attention'] = master[
        ['gene_ncbi', 'attention']].groupby('gene_ncbi').agg(sum)

    papers = papers.loc[ref_genes, :]
    papers = papers.fillna(0)

    return papers
