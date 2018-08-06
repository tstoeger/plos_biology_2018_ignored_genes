import glob
import math
import os

import numpy as np
import pandas as pd
import seaborn as sns

from access_biology_data import annotation, gwas_studies, meta, phenotype_collections, properties, relations
from access_literature_data import medline
from access_mixed_data import genealacart
from access_science_shared import standardizer, utils

from access_science_shared import inout

import nar170604f_occurences as nar_attention
import nar170605f_funding as nar_funding
import nar170823f_prediction_datasets as pred


import resci_inout as rinout
import resci_tools as ret

from natsort import natsorted


"""
Formats heterogeneous data for a resource-type display
cl      cleaned output for visualization
dd      detailed data
ge      genes considered in dataset

"""


def get_ref_genes():
    ref_genes = standardizer.reference_genes(
        taxon_id=9606,
        ref_code='orp')
    return ref_genes


def get_publications():
    ref_genes = get_ref_genes()
    ref_gene2pubmed = medline.gene2pubmed(
        taxon_id=9606,
        paper_kind='research',
        ref_genes=ref_genes)

    papers = nar_attention.count_papers_and_attention(
        ref_genes,
        ref_gene2pubmed)

    df_m = medline.select_medline_records(
        columns_sql='''
            medline.pubmed_id,
            medline.pubdate_year''',
        taxon_id=9606,
        kind='research')

    current_gene2pubmed = ref_gene2pubmed[ref_gene2pubmed['pubmed_id'].isin(
        df_m[df_m['pubdate_year'] == 2015]['pubmed_id'])]

    current_papers = nar_attention.count_papers_and_attention(
        ref_genes,
        current_gene2pubmed)

    current_papers = current_papers.rename(columns={
        'attention': 'attention_2015'
    })

    papers = pd.concat([papers, current_papers[['attention_2015']]], axis=1)

    return papers


def LoF():
    he = properties.allelepool_lek_2016_aberration(taxon_id=9606)
    he = he[he['gene_ncbi'].isin(get_ref_genes())]

    dd = he[['gene_ncbi', 'Population variability Lek pLI']
            ].set_index('gene_ncbi')

    cl = dd.copy()
    cl.loc[:, 'extremly_LoF_intolerant'] = cl[
        'Population variability Lek pLI'] >= 0.9
    cl = cl[['extremly_LoF_intolerant']]

    ge = sorted(list(cl.index))

    return cl, dd, ge


def DUF():

    he = annotation.interpro(
        taxon_id=9606)[
        ['gene_ncbi', 'interpro_name']].drop_duplicates()
    he = he[he['gene_ncbi'].isin(get_ref_genes())]

    ge = sorted(
        list(he['gene_ncbi'].unique()))

    he = he.copy()

    f = he['interpro_name'].str.contains('Domain of unknown function')
    he.loc[:, 'unknown_function'] = f

    he = he.loc[f, :].copy()
    he.loc[:, 'name_of_unknown'] = he['interpro_name'].str.extract(
        'Domain of unknown function ([A-Z0-9]*)', expand=False)

    f = he['name_of_unknown'].apply(lambda x: len(x) == 0)
    he.loc[f, 'name_of_unknown'] = 'unnamed'

    dd = he[['gene_ncbi', 'name_of_unknown', 'unknown_function']].pivot(
        index='gene_ncbi',
        values='unknown_function',
        columns='name_of_unknown')

    dd = dd.reindex(ge)
    dd = dd.fillna(False)
    dd = dd.loc[:, natsorted(dd.columns)]

    cl = dd.any(1).to_frame('has_duf')

    return cl, dd, ge


def extreme_swissprot():

    he = properties.aminoacids_swissprot(taxon_id=9606)
    he = he.set_index('gene_ncbi')
    he.columns = [x[len('Aminoacids_swissprot: '):] for x in he.columns]

    he_a = he.copy()

    he = properties.seg_swissprot(taxon_id=9606)
    he = he.set_index('gene_ncbi')
    he.columns = [x[len('SEG_swissprot: '):] for x in he.columns]

    he_b = he.copy()

    he_b = he_b[
        ['Fraction_of_protein', 'Fraction_of_protein_by_longest']].rename(
        columns={
            'Fraction_of_protein': 'SEG_fraction',
            'Fraction_of_protein_by_longest': 'SEG_fraction_longest',
        })

    he = properties.radar_swissprot(taxon_id=9606)
    he = he.set_index('gene_ncbi')

    he = he.divide(he_a[
        'amount_measured_amino_acids'], axis='rows').rename(columns={
            'Radar_swissprot_best: length': 'RADAR_fraction',
        })[['RADAR_fraction']]

    he_c = he.copy()

    he = pd.concat([he_a, he_b, he_c], axis=1)

    he.loc[:, 'minus_isoelectric_point'] = -he['isoelectric_point']

    c = ['acidic',
         'amount_measured_amino_acids', 'aromatic', 'basic', 'charged',
         'gravy_ignoring_O_and_U', 'helix_affine',
         'hydrophobic', 'isoelectric_point', 'molecular_weight', 'polar',
         'polar_uncharged', 'sheet_affine', 'turn_affine', 'SEG_fraction',
         'SEG_fraction_longest', 'RADAR_fraction', 'minus_isoelectric_point']

    he = he.loc[:, c]

    he = he[he['amount_measured_amino_acids'] >= 100]

    he = he[he.index.isin(get_ref_genes())]

    dd = he.copy()
    ge = sorted(he.index)

    cl = he.copy()
    cl = dd.rank(pct=True) >= 0.99
    cl['extreme_swissprot'] = cl.any(axis=1)

    return cl, dd, ge


def gtx():

    p = rinout.get_internal_path(
        (
            '170830f_differential_gene_expression_from_gtx/'
            'taxon9606/gtx/df_genes.csv.gz'
        )
    )

    df = pd.read_csv(p)

    gi = _get_gene_ncbi_2_ensembl()

    df = pd.merge(df, gi).drop('gene_ensembl', axis=1)
    df = df[df['gene_ncbi'].isin(get_ref_genes())]

    df.loc[:, 'differential'] = df.loc[:, 'p-value'] < 0.0001
    he = df[['gene_ncbi', 'differential']].groupby('gene_ncbi').agg(
        np.mean).rename(columns={'differential': 'gtx_fraction'})

    ge = natsorted(he.index)

    h = np.mean(df['differential'])
    f = he['gtx_fraction'] == 0
    he.loc[~f, 'gtx_fold'] = np.log2(he.loc[~f, 'gtx_fraction'] / h)
    he.loc[f, 'gtx_fold'] = -np.inf

    dd = he.copy()

    f = dd['gtx_fraction'].rank(pct=True) >= 0.99

    dd.loc[:, 'extreme_gtx'] = f

    cl = dd.copy()

    f = cl['gtx_fold'] > 2
    cl.loc[f, 'gtx_fold'] = 2

    f = cl['gtx_fold'] < -2
    cl.loc[f, 'gtx_fold'] = -2

    return cl, dd, ge


def omim_disease():

    he = annotation.omim_genealacart(
        taxon_id=9606,
        add_absenece=True).drop_duplicates()

    he = he[he['gene_ncbi'].isin(get_ref_genes())]
    ge = he['gene_ncbi'].unique()

    f = he['omim_disease'] == 'No entry in OMIM'

    he = he.loc[~f, :]

    he.loc[:, 'has'] = True

    dd = he.pivot(index='gene_ncbi', columns='omim_disease', values='has')
    dd = dd.reindex(ge)
    dd = dd.fillna(False)
    dd = dd.loc[:, natsorted(dd.columns)]

    cl = dd.any(axis=1).to_frame(
        'has_omim_disease')

    return cl, dd, ge


def human_phenotype():

    he = annotation.human_phenotype_genealacart(
        taxon_id=9606,
        add_absenece=True).drop_duplicates()

    he = he[he['gene_ncbi'].isin(get_ref_genes())]
    ge = he['gene_ncbi'].unique()

    he = he.rename(columns={
        'human_phenotype_genealacart: human_phenotype_name': 'human_phenotype'
    })

    f = he['human_phenotype'] == 'No known human phenotype'

    he = he.loc[~f, :]

    he.loc[:, 'has'] = True

    dd = he.pivot(index='gene_ncbi', columns='human_phenotype', values='has')
    dd = dd.reindex(ge)
    dd = dd.fillna(False)
    dd = dd.loc[:, natsorted(dd.columns)]

    cl = dd.any(axis=1).to_frame(
        'has_human_phenotype')

    return cl, dd, ge


def unified_disease():

    he = annotation.disease_genealacart(
        taxon_id=9606,
        add_absenece=True).drop_duplicates()

    he = he[he['gene_ncbi'].isin(get_ref_genes())]
    ge = he['gene_ncbi'].unique()

    f = he['unified_disease'] == 'No known disease'

    he = he.loc[~f, :]

    he.loc[:, 'has'] = True

    dd = he.pivot(index='gene_ncbi', columns='unified_disease', values='has')
    dd = dd.reindex(ge)
    dd = dd.fillna(False)
    dd = dd.loc[:, natsorted(dd.columns)]

    cl = dd.any(axis=1).to_frame(
        'has_unified_disease')

    return cl, dd, ge


def orphan_disease():

    he = annotation.orphanet_genealacart(
        taxon_id=9606,
        add_absenece=True).drop_duplicates()

    he = he[he['gene_ncbi'].isin(get_ref_genes())]
    ge = he['gene_ncbi'].unique()

    f = he['orphanet_disease'] == 'No entry in Orphanet'

    he = he.loc[~f, :]

    he.loc[:, 'has'] = True

    dd = he.pivot(index='gene_ncbi', columns='orphanet_disease', values='has')
    dd = dd.reindex(ge)
    dd = dd.fillna(False)
    dd = dd.loc[:, natsorted(dd.columns)]

    cl = dd.any(axis=1).to_frame(
        'has_orphan_disease')

    return cl, dd, ge


def rare_go():

    he = annotation.go(
        taxon_id=9606,
        category=['Process'],
        negating_support=[False],
        any_negating_support=[False],
        temporary_evidence=[False, True],
        unmapped_evidence=[False]
    )

    he = he[he['gene_ncbi'].isin(get_ref_genes())]
    ge = he['gene_ncbi'].unique()

    he = he[['gene_ncbi', 'GO_term']].drop_duplicates()

    go_counts = he['GO_term'].value_counts()

    he.loc[:, 'has'] = he['GO_term'].isin(go_counts[go_counts <= 3].index)

    he = he[he['has']]

    dd = he.pivot(index='gene_ncbi', columns='GO_term', values='has')
    dd = dd.reindex(ge)
    dd = dd.fillna(False)
    dd = dd.loc[:, natsorted(dd.columns)]

    cl = dd.any(axis=1).to_frame(
        'has_rare_go')

    return cl, dd, ge


def fame_rank():
    p = get_publications()
    p = p[p.index.isin(get_ref_genes())]
    dd = p[['attention']].rank(pct=True, ascending=True)
    dd = dd.rename(columns={'attention': 'attention_rank'})
    cl = dd.copy()
    ge = cl.index
    return cl, dd, ge


def signal_peptide():

    signal_peptide = pd.merge(
        properties.signalp_swissprot(taxon_id=9606)[
            ['gene_ncbi', 'SignalP_swissprot: cleaved']],
        properties.signalp_trembl(taxon_id=9606)[
            ['gene_ncbi', 'SignalP_trembl: cleaved']],
        how='outer')

    signal_peptide[
        'SignalP_swissprot: cleaved'] = signal_peptide[
        'SignalP_swissprot: cleaved'] > 0
    signal_peptide[
        'SignalP_trembl: cleaved'] = signal_peptide[
        'SignalP_trembl: cleaved'] > 0

    he = signal_peptide.set_index(
        'gene_ncbi').max(1).to_frame('signal_peptide')

    he = he[he.index.isin(get_ref_genes())]

    cl = he.copy()
    dd = he.copy()
    ge = natsorted(he.index)

    return cl, dd, ge


def rnai_phenotypes():

    he = phenotype_collections.genome_rnai(taxon_id=9606)
    f = he['phenotype'].str.contains('shRNA abundance')
    he = he[~he['pubmed_id'].isin(he.loc[f, 'pubmed_id'])].copy()

    he = he[~he['gene_ncbi'].str.contains(',')]
    he = he[he['gene_ncbi'] != '']
    he['gene_ncbi'] = he['gene_ncbi'].astype(float)

    he = he[he['gene_ncbi'].isin(get_ref_genes())]

    he = he.copy()

    he.loc[:, 'has_phenotype'] = he.loc[:, 'phenotype'] != 'none'

    g = he[['gene_ncbi', 'has_phenotype']].groupby(['gene_ncbi'])

    d = pd.concat([
        g.agg(np.mean),
        g.size().rename('studies')
    ], axis=1).reset_index()

    minimally_required_studies = 20
    d = d[d['studies'] >= minimally_required_studies]

    minimally_required_studies = 20
    d = d[d['studies'] >= minimally_required_studies]

    d = d[['gene_ncbi', 'has_phenotype']].rename(
        columns={'has_phenotype': 'fraction_rnai_phenotype'})

    d = d.set_index('gene_ncbi')

    dd = d.copy()

    d.loc[:, 'rnai_frequent'] = d.loc[:, 'fraction_rnai_phenotype'] > 0.3
    cl = d[['rnai_frequent']]

    ge = sorted(cl.index)

    return cl, dd, ge


def rare_compounds():
    he = genealacart.load_genealacart_dataset('UnifiedCompounds')

    he = he[he['gene_ncbi'].isin(get_ref_genes())]

    ge = he['gene_ncbi'].unique()

    he = he[he['CasNumbers'].notnull()]
    he = he[['gene_ncbi', 'CasNumbers']].drop_duplicates()

    c = he['CasNumbers'].value_counts()

    g = he[he['CasNumbers'].isin(c[c <= 10].index)]

    g = g['gene_ncbi'].value_counts().to_frame(
        'rare_compounds').rename_axis('gene_ncbi').reindex(ge)

    g = g.fillna(0).sort_index()

    dd = g.copy()

    cl = g > 0

    return cl, dd, ge


def fame_in_bioplex():

    he = relations.bioplex2()
    papers = get_publications()

    he = he[he['gene_ncbi'].isin(get_ref_genes())]

    ge = he['gene_ncbi'].unique()

    h = pd.merge(
        he,
        papers[['attention']].reset_index())

    hh = pd.merge(
        h,
        h,
        left_on='bioplex2_id',
        right_on='bioplex2_id',
        suffixes=('_cis', '_trans')
    )

    b = hh[hh['gene_ncbi_cis'] != hh['gene_ncbi_trans']
           ][['gene_ncbi_cis', 'attention_trans']]

    g = b.groupby(['gene_ncbi_cis'])

    dd = pd.concat([
        g.agg(max).rename(
            columns={'attention_trans': 'max_attention_others_in_complex'}),
        g.agg(min).rename(
            columns={'attention_trans': 'min_attention_others_in_complex'}),
        g.agg(np.mean).rename(
            columns={'attention_trans': 'mean_attention_others_in_complex'}),
    ], axis=1).rename_axis('gene_ncbi')

    cl = dd.copy()

    cl['max_attention_others_in_complex'] = cl[
        'max_attention_others_in_complex'] > 100
    cl['min_attention_others_in_complex'] = cl[
        'min_attention_others_in_complex'] < 1
    cl['mean_attention_others_in_complex'] = cl[
        'mean_attention_others_in_complex'] < 1

    cl = cl.rename(columns={
        'max_attention_others_in_complex': 'bp2_with_studied',
        'min_attention_others_in_complex': 'bp2_with_unstudied',
        'mean_attention_others_in_complex': 'bp2_mean_unstudied'
    })

    return cl, dd, ge


def frequent_gwas():

    ebi_gwas = gwas_studies.ebi_gwas()

    f = ebi_gwas['MAPPED_GENE'].str.contains('[;,-]') == True
    gwas = ebi_gwas.loc[
        ~f,
        ['MAPPED_GENE', 'DISEASE/TRAIT', 'PVALUE_MLOG', 'pubmed_id']].rename(
        columns={
            'MAPPED_GENE': 'symbol_ambiguous',
            'DISEASE/TRAIT': 'trait',
            'PVALUE_MLOG': 'log_pvalue'
        }
    )

    gwas = pd.merge(
        gwas,
        meta.gene_info(taxon_id=9606, usecols=[
                       'symbol_ncbi', 'gene_ncbi']),
        left_on='symbol_ambiguous',
        right_on='symbol_ncbi',
        how='inner'
    ).drop('symbol_ambiguous', axis=1).drop('symbol_ncbi', axis=1)

    gwas = gwas[gwas['gene_ncbi'].isin(get_ref_genes())]

    ge = sorted(gwas['gene_ncbi'].unique())

    gwas = gwas.sort_values('log_pvalue', ascending=False)
    gwas = gwas.drop_duplicates(
        ['trait', 'pubmed_id', 'gene_ncbi'],
        keep='first')

    studies_per_phenotype = gwas[
        ['pubmed_id', 'trait']].drop_duplicates()[
        'trait'].value_counts()

    required_studies = 10
    important_gwas = gwas.loc[
        (
            gwas['trait'].isin(
                studies_per_phenotype[
                    studies_per_phenotype >= required_studies].index)), :

    ][['pubmed_id', 'trait', 'gene_ncbi']].drop_duplicates()

    he = pd.merge(
        important_gwas.groupby(
            ['trait', 'gene_ncbi']).size(
        ).reset_index().rename(columns={0: 'records'}),
        studies_per_phenotype.to_frame(
            'studies').reset_index().rename(columns={'index': 'trait'}))

    he.loc[:, 'fraction_of_gwas_studies'] = he['records'] / he['studies']

    dd = he.pivot(
        index='gene_ncbi',
        columns='trait',
        values='fraction_of_gwas_studies'
    )

    dd.columns = ['gwas_prominent_{}'.format(x) for x in dd.columns]

    dd = dd.reindex(ge)
    dd = dd.fillna(0)

    cl = dd > 0.2
    cl.loc[:, 'any_prominent_gwas'] = cl.any(axis=1)

    return cl, dd, ge


def any_gwas():

    ebi_gwas = gwas_studies.ebi_gwas()

    f = ebi_gwas['MAPPED_GENE'].str.contains('[;,-]') == True
    gwas = ebi_gwas.loc[
        ~f,
        ['MAPPED_GENE', 'DISEASE/TRAIT', 'PVALUE_MLOG', 'pubmed_id']].rename(
        columns={
            'MAPPED_GENE': 'symbol_ambiguous',
            'DISEASE/TRAIT': 'trait',
            'PVALUE_MLOG': 'log_pvalue'
        }
    )

    gwas = pd.merge(
        gwas,
        meta.gene_info(taxon_id=9606, usecols=[
                       'symbol_ncbi', 'gene_ncbi']),
        left_on='symbol_ambiguous',
        right_on='symbol_ncbi',
        how='inner'
    ).drop('symbol_ambiguous', axis=1).drop('symbol_ncbi', axis=1)

    gwas = gwas[gwas['gene_ncbi'].isin(get_ref_genes())]

    ge = sorted(gwas['gene_ncbi'].unique())

    gwas = gwas.sort_values('log_pvalue', ascending=False)
    gwas = gwas.drop_duplicates(
        ['trait', 'pubmed_id', 'gene_ncbi'],
        keep='first')

    studies_per_phenotype = gwas[
        ['pubmed_id', 'trait']].drop_duplicates()[
        'trait'].value_counts()

    required_studies = 1
    important_gwas = gwas.loc[
        (
            gwas['trait'].isin(
                studies_per_phenotype[
                    studies_per_phenotype >= required_studies].index)), :

    ][['pubmed_id', 'trait', 'gene_ncbi']].drop_duplicates()

    he = pd.merge(
        important_gwas.groupby(
            ['trait', 'gene_ncbi']).size(
        ).reset_index().rename(columns={0: 'records'}),
        studies_per_phenotype.to_frame(
            'studies').reset_index().rename(columns={'index': 'trait'}))

    he.loc[:, 'fraction_of_any_gwas_studies'] = he['records'] / he['studies']

    dd = he.pivot(
        index='gene_ncbi',
        columns='trait',
        values='fraction_of_any_gwas_studies'
    )

    dd.columns = ['gwas_any_{}'.format(x) for x in dd.columns]

    dd = dd.reindex(ge)
    dd = dd.fillna(0)

    cl = dd > 0.0
    cl.loc[:, 'any_gwas'] = cl.any(axis=1)

    return cl, dd, ge


def challenged_proteins():

    p = inout.get_path('publications', 'ezkurdia2014/ddu309supp_tables1.xlsx')

    df_all = pd.read_excel(p, sheet_name='All G12 genes')

    df_tagged = pd.read_excel(p, sheet_name='Possible non-coding set')

    gi = _get_gene_ncbi_2_ensembl()

    if (df_tagged['ENSEMBL'].isin(df_all['ENSEMBL'])).all():
        ge = sorted(gi[gi['gene_ensembl'].isin(
            df_all['ENSEMBL'])]['gene_ncbi'])

    df_suspicious = gi.copy()

    df_suspicious = df_suspicious[df_suspicious['gene_ensembl'].isin(
        df_all['ENSEMBL'])]

    df_suspicious.loc[
        :,
        'ezkurdia_challenged'] = df_suspicious.loc[:, 'gene_ensembl'].isin(
            df_tagged['ENSEMBL'])

    dd = df_suspicious[['gene_ncbi', 'ezkurdia_challenged']
                       ].set_index('gene_ncbi')
    cl = dd.copy()

    return cl, dd, ge


def detection_in_cells():
    he = properties.transcript_abundance_uhlen_2015_cells()
    he = he[he['gene_ncbi'].isin(get_ref_genes())]
    he = he.set_index('gene_ncbi')

    dd = pd.concat([
        ((he >= 0).sum(axis=1) / he.shape[1]).to_frame('fraction_of_cells_1'),
        ((he >= 1).sum(axis=1) / he.shape[1]).to_frame('fraction_of_cells_10'),
        ((he >= 2).sum(axis=1) / he.shape[1]).to_frame('fraction_of_cells_100')
    ], 1)

    cl = dd.copy()
    ge = sorted(dd.index)

    return cl, dd, ge


def detection_in_tissues():
    he = properties.transcript_abundance_uhlen_2015_tissues()
    he = he[he['gene_ncbi'].isin(get_ref_genes())]
    he = he.set_index('gene_ncbi')

    dd = pd.concat([
        ((he >= 0).sum(axis=1) / he.shape[1]
         ).to_frame('fraction_of_tissues_1'),
        ((he >= 1).sum(axis=1) / he.shape[1]
         ).to_frame('fraction_of_tissues_10'),
        ((he >= 2).sum(axis=1) / he.shape[1]
         ).to_frame('fraction_of_tissues_100'),
    ], 1)

    cl = dd.copy()
    ge = sorted(dd.index)

    return cl, dd, ge


def biogrid_western_blot():
    he = relations.biogrid(taxon_id=9606)

    genes = pd.Series(list(
        set(he['Entrez Gene Interactor A']).union(
            set(he['Entrez Gene Interactor B']))))
    ge = natsorted(genes[genes.isin(get_ref_genes())])

    he = he[he['Experimental System'] == 'Affinity Capture-Western']

    dd = pd.DataFrame(index=ge)

    dd.loc[:, 'biogrid_western_blot'] = dd.index.isin(
        he['Entrez Gene Interactor B'])

    cl = dd.copy()

    return cl, dd, ge


def fame_of_homologs():
    gene2pubmed_research = medline.gene2pubmed(
        taxon_id='all', paper_kind='research')

    value_of_pubmed_id = gene2pubmed_research[
        'pubmed_id'].value_counts().to_frame().reset_index().rename(
        columns={'index': 'pubmed_id', 'pubmed_id': 'value'})

    value_of_pubmed_id['value'] = 1 / value_of_pubmed_id['value']

    gene2pubmed_research = pd.merge(gene2pubmed_research, value_of_pubmed_id)
    extended_attention = gene2pubmed_research[[
        'gene_ncbi', 'value']].groupby('gene_ncbi').agg(sum)

    hg = relations.homologene()
    hg_attention = pd.merge(hg, extended_attention.reset_index(), how='left')
    hg_attention['value'] = hg_attention['value'].fillna(0)
    hg_max_attention = hg_attention[
        ['homologene_group', 'taxon_ncbi', 'value']].groupby(
        ['homologene_group', 'taxon_ncbi']).agg(max).reset_index()

    hg_max_attention = hg_max_attention[hg_max_attention['taxon_ncbi'] != 9606]

    he = pd.merge(
        hg_max_attention,
        hg[hg['taxon_ncbi'] == 9606][['homologene_group', 'gene_ncbi']]
    )

    he = he[he['gene_ncbi'].isin(get_ref_genes())]

    he = he[['taxon_ncbi', 'gene_ncbi', 'value']].groupby(
        ['gene_ncbi', 'taxon_ncbi']).agg(max).reset_index()

    dd = he.pivot(index='gene_ncbi', columns='taxon_ncbi', values='value')

    dd.columns = [meta.taxon_name(x) for x in dd.columns]
    dd.columns = ['studied_{}'.format(x) for x in dd.columns]

    cl = dd > 1

    f = dd.isnull()
    cl[f] = np.nan

    ge = sorted(dd.index)

    return cl, dd, ge


def presence_of_homologs():

    hg = relations.homologene()

    hg = pd.merge(
        hg[hg['taxon_ncbi'] == 9606],
        hg[['taxon_ncbi', 'homologene_group']],
        left_on='homologene_group',
        right_on='homologene_group',
        suffixes=('_cis', '_trans')
    )

    f = hg['taxon_ncbi_cis'] == hg['taxon_ncbi_trans']
    hg = hg.loc[~f, :]

    hg = hg.copy()

    hg.loc[:, 'present'] = True

    hg = hg[['gene_ncbi', 'taxon_ncbi_trans', 'present']].drop_duplicates()

    dd = hg.pivot(
        index='gene_ncbi',
        columns='taxon_ncbi_trans',
        values='present').fillna(False)

    dd.columns = [meta.taxon_name(x) for x in dd.columns]
    dd.columns = ['presence_{}'.format(x) for x in dd.columns]

    cl = dd.copy()

    ge = natsorted(dd.index)

    return cl, dd, ge


def _get_gene_ncbi_2_ensembl():

    gi = meta.gene_info(taxon_id=9606)

    f = gi['dbXrefs'].str.contains('Ensembl:ENSG[0-9]*')
    gi.loc[f, 'gene_ensembl'] = gi.loc[f, 'dbXrefs'].str.extract(
        'Ensembl:(ENSG[0-9]*)', expand=False)
    gi = gi[['gene_ncbi', 'gene_ensembl']].drop_duplicates()
    gi = gi.drop_duplicates('gene_ensembl', keep=False)
    gi = gi[gi['gene_ncbi'].isin(get_ref_genes())]

    return gi


def load_layout(rotation_degrees=0):

    p = rinout.get_internal_path(
        os.path.join(
            '171014f_visualize_metrics_on_characteristics_tsne_and_features/',
            'genes_coordinates.csv'))

    tsne_frame = pd.read_csv(p).set_index('gene_ncbi')

    if rotation_degrees != 0:
        agg = []
        for x in zip(tsne_frame.loc[:, 'x'], tsne_frame.loc[:, 'y']):
            d = ret.rotate((0, 0), x, rotation_degrees)
            agg.append(d)

        tsne_frame.loc[:, 'x'] = [x[0] for x in agg]
        tsne_frame.loc[:, 'y'] = [x[1] for x in agg]

    return tsne_frame


def load_group_annotation(annotation_code='180318'):

    # Add manually selected groups of genes to data to plot

    if annotation_code == 'pre_180318':
        gg = rinout.get_internal_path(
            (
                '171014f_visualize_metrics_on_characteristics_tsne'
                '_and_features/gene_list*'))

        agg = []
        for g in glob.glob(gg):
            df = pd.read_table(g, names=['gene_ncbi'])
            _, fn = os.path.split(g)
            df.loc[:, 'list'] = fn
            agg.append(df)
        df_manual_label = pd.concat(agg, axis=0)
        mislabelled_genes = [
            192670,
            154810,
            140710,
            170690,
            157570,
            374860,
            133690,
            353140,
            387700,
            319100,
            404550,
            283870,
            129450,
            121260,
            131870,
            123720,
            132720,
            203430,
            284110,
            219790,
            130540,
            158830,
            347730,
            124540,
            163590,
            283600,
            118460
        ]
        f = df_manual_label['gene_ncbi'].isin(mislabelled_genes)
        df_manual_label = df_manual_label.loc[~f, :]
        df_manual_label = df_manual_label.drop_duplicates(
            subset='gene_ncbi', keep=False).copy()
        df_manual_label.loc[:, 'list_code'] = df_manual_label.loc[
            :, 'list'].copy(
        ).str.extract('gene_list([0-9]*).*', expand=False).astype(float)

    elif annotation_code == '180318':
        p = rinout.get_internal_path(
            '180317_complete_cluster_maker/mark_180318.csv')
        df_manual_label = pd.read_csv(p).rename(columns={
            'group': 'list_code'
        })

    else:
        raise ValueError('Does not support provdied annotation_code.')

    df_manual_label = df_manual_label[df_manual_label['gene_ncbi'].isin(
        get_ref_genes())]

    return df_manual_label


def pi_transition():

    p = rinout.get_internal_path(
        (
            '180311_cache_pi_transition_for_genes/'
            '180311_cache_pi_transition_for_genes.csv')

    )

    pool = pd.read_csv(p, low_memory=False)

    pubmed_year_pi = pool[
        ['pubmed_id', 'pubdate_year', 'will_be_pi', 'genes']].copy()

    tolerated_genes_per_publication = 10
    pubmed_year_pi = pubmed_year_pi[
        pubmed_year_pi['genes'] <= tolerated_genes_per_publication]

    human_gene2pubmed = medline.gene2pubmed(
        taxon_id=9606,
        paper_kind='research',
        ref_genes=get_ref_genes())[['gene_ncbi', 'pubmed_id']]

    ma = pd.merge(human_gene2pubmed, pubmed_year_pi)

    # m = ma[['gene_ncbi', 'pubdate_year', 'will_be_pi']].groupby(
    #     ['gene_ncbi', 'pubdate_year']).agg(np.mean).reset_index()

    # av = m[[
    #     'pubdate_year', 'will_be_pi']].groupby(
    #         'pubdate_year').agg(np.mean).reset_index().rename(columns={
    #             'will_be_pi': 'per_year_occurence_will_be_pi'
    #         })

    # n = pd.merge(m, av)

    # f = n['pubdate_year'].isin(range(2010, 2011))

    # nn = n.loc[f, :].copy()

    # nn['above'] = nn['will_be_pi'] > (nn['per_year_occurence_will_be_pi']*2)

    # r = nn[['gene_ncbi', 'above']].groupby('gene_ncbi').agg(np.mean)

    # r = r.rename(columns={'above': 'recent_above_average'})

    # dd = r.copy()
    # ge = natsorted(r.index)
    # cl = r > 0.9

    # y = 2010

    m = ma[ma['pubdate_year'].isin(range(2010, 2016))]

    c = m['gene_ncbi'].value_counts()

    m = m[m['gene_ncbi'].isin(c[c >= 10].index)]

    dd = m[['gene_ncbi', 'will_be_pi']].groupby('gene_ncbi').agg(np.mean)

    a = np.log2(dd / dd['will_be_pi'].mean())

    # a[a > 2] = 2
    # a[a < -2] = -2

    dd = a.copy()

    cl = dd.copy() > 1    

    ge = cl.index

    return cl, dd, ge


def pi_transition_no_limit_on_number_of_studies():

    p = rinout.get_internal_path(
        (
            '180311_cache_pi_transition_for_genes/'
            '180311_cache_pi_transition_for_genes.csv')

    )

    pool = pd.read_csv(p, low_memory=False)

    pubmed_year_pi = pool[
        ['pubmed_id', 'pubdate_year', 'will_be_pi', 'genes']].copy()

    tolerated_genes_per_publication = 10
    pubmed_year_pi = pubmed_year_pi[
        pubmed_year_pi['genes'] <= tolerated_genes_per_publication]

    human_gene2pubmed = medline.gene2pubmed(
        taxon_id=9606,
        paper_kind='research',
        ref_genes=get_ref_genes())[['gene_ncbi', 'pubmed_id']]

    ma = pd.merge(human_gene2pubmed, pubmed_year_pi)

    # m = ma[['gene_ncbi', 'pubdate_year', 'will_be_pi']].groupby(
    #     ['gene_ncbi', 'pubdate_year']).agg(np.mean).reset_index()

    # av = m[[
    #     'pubdate_year', 'will_be_pi']].groupby(
    #         'pubdate_year').agg(np.mean).reset_index().rename(columns={
    #             'will_be_pi': 'per_year_occurence_will_be_pi'
    #         })

    # n = pd.merge(m, av)

    # f = n['pubdate_year'].isin(range(2010, 2011))

    # nn = n.loc[f, :].copy()

    # nn['above'] = nn['will_be_pi'] > (nn['per_year_occurence_will_be_pi']*2)

    # r = nn[['gene_ncbi', 'above']].groupby('gene_ncbi').agg(np.mean)

    # r = r.rename(columns={'above': 'recent_above_average'})

    # dd = r.copy()
    # ge = natsorted(r.index)
    # cl = r > 0.9

    # y = 2010

    m = ma[ma['pubdate_year'].isin(range(2010, 2016))]

    c = m['gene_ncbi'].value_counts()

    m = m[m['gene_ncbi'].isin(c[c >= 1].index)]

    dd = m[['gene_ncbi', 'will_be_pi']].groupby('gene_ncbi').agg(np.mean)

    a = np.log2(dd / dd['will_be_pi'].mean())

    # a[a > 2] = 2
    # a[a < -2] = -2

    dd = a.copy()

    cl = dd.copy() > 3

    ge = cl.index

    return cl, dd, ge


def supporting_nih_institutes():

    df_prj_core, _, df_papers = nar_funding.get_paper_funding_through_nih()

    ref_genes = get_ref_genes()
    ref_gene2pubmed = medline.gene2pubmed(
        taxon_id=9606,
        paper_kind='research',
        ref_genes=ref_genes)

    f = df_papers['pubdate_year'].isin(range(2015, 2016))

    d = pd.merge(
        ref_gene2pubmed,
        pd.merge(
            df_prj_core[['project_num', 'ADMINISTERING_IC']].drop_duplicates(),
            df_papers.loc[f, ['project_num', 'pubmed_id']].drop_duplicates()))

    d = d[d['gene_ncbi'].isin(ref_genes)]

    dc = d[['gene_ncbi', 'ADMINISTERING_IC']].drop_duplicates()

    dc = dc['gene_ncbi'].value_counts()

    dd = dc.to_frame('recently_supporting_institutes').rename_axis('gene_ncbi')

    total_institutes = len(set(d['ADMINISTERING_IC']))
    dd = dd / total_institutes

    cl = dd >= 0.1

    ge = natsorted(d['gene_ncbi'].unique())

    return cl, dd, ge


def get_clustered_zscored_features():

    features_to_use = [
        'Population variability Lek mis_z',
        'SignalP_swiss_or_trembl: cmax',
        'uhlen_2015_cells_log10fpkm: appendices_4b',
        'Population variability Lek lof_z',
        'uhlen_2015_cells_log10fpkm: liver_c',
        'Genbank__gene: SumACGT',
        'uhlen_2015_cells_log10fpkm: brain_3c',
        'uhlen_2015_fraction_detection_tissues',
        'Population variability Lek pNull',
        'Protein Itzhak Itzhak2016_Contribution to cell protein mass [ppm]',
        'uhlen_2015_cells_log10fpkm: adrenal_4d',
        'Aminoacids_swiss_or_trembl: gravy_ignoring_O_and_U',
        'Wang2015: KBM7 CS', 'Aminoacids_swiss_or_trembl: basic',
        'Genbank_validated_RNA: full_SumACGT']

    ref_genes = get_ref_genes()

    bio = pred.retreive_biophysics(
        ref_genes=ref_genes,
        taxon_id=9606)
    exp = pred.retreive_human_experiments(
        ref_genes=ref_genes,
        taxon_id=9606)

    bioexp = list(bio.values()) + list(exp.values())

    features = pd.concat(bioexp, axis=1)

    features = features.loc[
        ref_genes,
        features_to_use]

    features = features.dropna().astype(float)

    def normfun(x):
        m = np.nanmean(x)
        s = np.nanstd(x)
        a = (x - m) / s
        return a
    n_features = features.apply(lambda x: normfun(x), 0)

    cl = sns.clustermap(n_features, method='ward',
                        vmin=-3, vmax=3, cmap='PiYG')

    clustered_features = n_features.iloc[
        cl.dendrogram_row.reordered_ind,
        cl.dendrogram_col.reordered_ind]

    return clustered_features
