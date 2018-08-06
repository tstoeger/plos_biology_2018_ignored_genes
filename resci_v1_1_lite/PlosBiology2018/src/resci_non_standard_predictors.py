import numpy as np
import pandas as pd

from copy import deepcopy

from access_biology_data import meta, relations
from access_literature_data import medline
from access_science_shared import standardizer, inout

import nar170902_discovery_year as nar_year
import nar170604f_occurences as nar_occurences
import nar170604f_occurences as nar_attention


def literature_of_homologenes(taxon_id, categ):
    """
    Obtains research literature of protein-coding
    homologous genes of other organisms for the organism
    specified by taxon_id; output will be nan if no
    homolog; note that literature is not normalized
    or transformed

    Input:
        taxon_id    int; taxon of interest
        categ       str; literature category of interest
                        e.g.: 'papers' or 'attention'
    Output:
        df          literature statistics for other
                    organsisms (mean, median, sum)
    """

    taxon_id_cis = deepcopy(taxon_id)
    fame_of_interest = categ

    hg = relations.homologene()

    if taxon_id_cis not in hg['taxon_ncbi'].unique():
        raise ValueError(
            'The taxon {} is absent from homolgene.'.format(
                taxon_id_cis))

    trans_taxa = set(
        hg['taxon_ncbi'].unique()
    ) - set([taxon_id_cis])

    agg = []

    for taxon_id in sorted(list(trans_taxa)):
        taxon_id = int(taxon_id)
        ref_genes = standardizer.reference_genes(
            taxon_id, 'p')  # protein-coding
        ref_gene2pubmed = medline.gene2pubmed(
            taxon_id, ['pubmed_id', 'gene_ncbi'], paper_kind='research')
        papers = nar_occurences.count_papers_and_attention(
            ref_genes, ref_gene2pubmed)

        fame_on_hg_on_taxon = pd.merge(
            papers[[fame_of_interest]].reset_index(),
            hg,
            how='inner')

        agg.append(fame_on_hg_on_taxon)

    fame_of_homologs = pd.concat(agg, axis=0)

    stats_to_run = {
        'median': np.median,
        'sum': np.sum,
        'mean': np.mean,
    }

    agg = []
    for norm_name, norm_method in stats_to_run.items():
        df = fame_of_homologs[
            ['taxon_ncbi', 'homologene_group', fame_of_interest]].groupby(
            ['taxon_ncbi', 'homologene_group']).agg(
            norm_method)

        df = df.reset_index(
        ).pivot(
            index='homologene_group',
            columns='taxon_ncbi',
            values=fame_of_interest)

        df.columns = [
            'homolog_{}_{}_from_{}'.format(
                fame_of_interest,
                norm_name,
                meta.taxon_name(x).replace(' ', '_')
            ) for x in df.columns]

        agg.append(df)

    literature_of_homologs = pd.concat(agg, axis=1)

    cis_hg = hg[
        hg['taxon_ncbi'] == taxon_id_cis][
        ['homologene_group', 'gene_ncbi']]

    df = pd.merge(
        cis_hg,
        literature_of_homologs.reset_index(),
        how='left').drop(
            'homologene_group', axis=1).set_index(
        'gene_ncbi', verify_integrity=True).sort_index(
    ).reset_index()

    return df


def literature_of_rolland_2014_interactors(taxon_id, categ):
    """
    Obtains research literature of protein-coding
    homologous genes of interctors according to Rolland et al. 2014;
    output will be nan if no interactor
    reported, but tested in Rolland et al; note that literature
    is not transformed

    Input:
        taxon_id    int; taxon of interest
        categ       str; literature category of interest
                        e.g.: 'papers' or 'attention'
    Output:
        df          literature statistics for
                        interactors (mean, median, sum)
    """

    if taxon_id != 9606:
        raise EnvironmentError(
            'Only taxon 9606 (human) supported.')

    p = inout.get_path(
        'geisen',
        'papers/rolland_2014/rolland_table_binary_interactions.csv.gz'
    )

    df = pd.read_csv(p)
    df = df[['entrez_gene_ida', 'entrez_gene_idb']]
    df.columns = ['a', 'b']

    df = pd.concat(
        [df, df.rename(
            columns={'a': 'b', 'b': 'a'})], axis=0, ignore_index=True, )

    f = df['a'] == df['b']
    df = df.loc[~f, :]

    ref_genes = df['a'].unique()
    gene2pubmed = medline.gene2pubmed(taxon_id=9606, paper_kind='research')

    co = nar_occurences.count_papers_and_attention(ref_genes, gene2pubmed)
    co = co[[categ]].reset_index()

    jo = pd.merge(
        df,
        co,
        left_on='b',
        right_on='gene_ncbi',
        how='left')

    jo = jo[['a', categ]].rename(columns={'a': 'gene_ncbi'})

    stats_to_run = {
        'median': np.median,
        'sum': np.sum,
        'mean': np.mean,
    }

    agg = []
    for norm_name, norm_method in stats_to_run.items():
        df = jo.groupby('gene_ncbi').agg(norm_method)

        df.columns = [
            'interactor_{}_{}'.format(
                categ,
                norm_name,
            ) for x in df.columns]
        agg.append(df)

    df = pd.concat(agg, axis=1)

    p = inout.get_path(
        'geisen',
        'papers/rolland_2014/rolland_considered_genes.csv.gz'
    )

    genes_in_rolland = pd.read_csv(p)
    v = list(set(genes_in_rolland['rolland_considered_genes'].values))

    literature_of_interactor = df.loc[v, :]
    literature_of_interactor = literature_of_interactor.reset_index()

    literature_of_interactor = literature_of_interactor.sort_values(
        'gene_ncbi').reset_index(drop=True)

    return literature_of_interactor


def own_rpo_log_fame(taxon_id):
    """
    Returns log10 transformed amount of research publciations
    and fame for rpo (research, protein-coding, official)
    genes of taxon_id

    Input:
        taxon_id    int     taxon_id

    Output:
        df          DataFrame with log10 transformed papers
                        and attention
    """
    ref_genes = standardizer.reference_genes(
        taxon_id,
        'rpo')
    gene2pubmed = medline.gene2pubmed(
        taxon_id,
        paper_kind='research',
        ref_genes=ref_genes)

    fame = nar_attention.count_papers_and_attention(
        ref_genes, gene2pubmed)
    fame = fame.applymap(lambda x: np.log10(x))
    fame.columns = ['own_log_{}'.format(x) for x in fame.columns]
    fame = fame.reset_index()

    return fame


def year_of_discovery_in_homologenes(taxon_id, categ):
    """
    Obtains year of initial discovery in other homologenes;
    Note that here it only keeps information about genes that had
    been repoted earlier in other taxon

    input:
        taxon_id    int   cis taxon of interst
        categ       str   either first_solo_year or first_year
    output:
        df          year of discoveries, NaN either means no homolog
                        or discovery after - or in same year as - taxon_id
    """

    # internal settings
    only_keep_preceding_years = True

    hg_cis = _get_disovery_in_homologene(
        taxon_id,
        categ,
        only_keep_preceding_years)
    return hg_cis


def year_of_description_in_homologenes(taxon_id, categ):
    """
    Obtains year of initial discovery in other homologenes;
    Note that here it only keeps information about genes that had
    been repoted earlier in other taxon

    input:
        taxon_id    int   cis taxon of interst
        categ       str   either first_solo_year or first_year
    output:
        df          year of discoveries, NaN either means no homolog
                        or discovery after - or in same year as - taxon_id
    """

    # internal settings
    only_keep_preceding_years = False

    hg_cis = _get_disovery_in_homologene(
        taxon_id,
        categ,
        only_keep_preceding_years)
    return hg_cis


def year_of_description(taxon_id, categ):
    """
    Obtains year of description in a given taxon;
    Note that here it only keeps information about genes that had
    been repoted earlier in other taxon

    input:
        taxon_id    int   cis taxon of interst
        categ       str   either first_solo_year or first_year
    output:
        df          year of discoveries, NaN either means no homolog
                        or discovery after - or in same year as - taxon_id
    """

    ref_genes = standardizer.reference_genes(taxon_id, '')
    df = nar_year.get_year_of_discovery(taxon_id, ref_genes)

    df = df.loc[:, [categ]].dropna()
    df = df.reset_index()
    return df


def _get_disovery_in_homologene(taxon_id, categ, only_keep_preceding_years):

    # load homologene
    hg = relations.homologene()

    # get years of publications
    agg = list()
    for taxon_trans in hg['taxon_ncbi'].unique():
        taxon_trans = int(taxon_trans)
        ref_genes = standardizer.reference_genes(
            taxon_id=taxon_trans, ref_code='p')
        gene2pubmed = medline.gene2pubmed(
            taxon_id=taxon_trans,
            paper_kind='research',
            ref_genes=ref_genes)

        df_m = medline.select_medline_records(
            columns_sql='''
                    medline.pubmed_id,
                    medline.pubdate_year''',
            taxon_id=taxon_trans,
            kind='research')

        df_m = df_m[df_m['pubmed_id'].isin(gene2pubmed['pubmed_id'])]
        columns_to_use = ['pubmed_id', 'pubdate_year']
        df_m = df_m.loc[:, columns_to_use].drop_duplicates()

        genes_per_paper = gene2pubmed['pubmed_id'].value_counts(
        ).to_frame('genes')
        df_m = pd.merge(df_m, genes_per_paper, left_on='pubmed_id',
                        right_index=True, how='inner')
        df_m.loc[:, 'taxon_ncbi'] = taxon_trans
        agg.append(df_m)

    # add genes to medline
    master = pd.merge(
        pd.concat(agg, axis=0),
        medline.gene2pubmed(
            taxon_id='all',
            paper_kind='research'),
        left_on=['taxon_ncbi', 'pubmed_id'],
        right_on=['taxon_ncbi', 'pubmed_id'],
        how='inner').drop_duplicates()

    # get initial years
    d = master
    is_single_gene_paper = d['genes'] == 1
    genes_earliest_years = pd.merge(
        d.loc[
            :,
            ['gene_ncbi', 'pubdate_year']].groupby(
                'gene_ncbi').agg(min).reset_index().rename(
                    columns={'pubdate_year': 'first_year'}),
        d.loc[
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

    # create table with years of discovery
    hgp = pd.merge(hg, genes_earliest_years, left_on='gene_ncbi',
                   right_index=True, how='left')

    phg = hgp[['homologene_group', 'taxon_ncbi', categ]].groupby(
        ['homologene_group', 'taxon_ncbi']).agg(min).reset_index().pivot(
            index='homologene_group', columns='taxon_ncbi', values=categ)

    hg_cis = hg.loc[hg['taxon_ncbi'] == taxon_id,
                    ['homologene_group', 'gene_ncbi']]
    hg_cis = pd.merge(hg_cis, phg.reset_index()).drop(
        'homologene_group', axis=1)
    hg_cis = hg_cis.set_index('gene_ncbi', verify_integrity=True)

    # format output
    if only_keep_preceding_years:
        for c in hg_cis.columns:
            f = hg_cis.loc[:, c] >= hg_cis.loc[:, taxon_id]
            hg_cis.loc[f, c] = np.nan

    hg_cis = hg_cis.drop(taxon_id, axis=1)

    if only_keep_preceding_years:
        hg_cis.columns = ['preceding_{}_in_{}'.format(
            categ, meta.taxon_name(x)) for x in hg_cis.columns]
    else:
        hg_cis.columns = ['description_{}_in_{}'.format(
            categ, meta.taxon_name(x)) for x in hg_cis.columns]

    hg_cis = hg_cis.sort_index()
    hg_cis = hg_cis.reset_index()
    return hg_cis
