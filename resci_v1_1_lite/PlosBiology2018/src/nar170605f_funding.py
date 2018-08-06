import numpy as np
import pandas as pd


from access_literature_data import medline
from access_economic_data import inflation, nih


import nar170604f_occurences as nar_attention
import ana170508f_human_citations as ana


from access_economic_data import nih
from access_literature_data import medline, wos
from access_science_shared import standardizer


def get_paper_funding_through_nih():
    """
    Gets NIH funding information for papers. The common column
    in all outputs is the NIH project number. In addition
    fiscal years, and publication years are given to potentially
    test time sensitivity.

    Data processing:
        - ignore applicaitons without defined core project
        - ignore applications with wildcard in project number
            (000000 in number)
        - exclude applications to sub-projects (note: which would
            still have the same core project number)
        - Filter for presence in medline with mappable year
        - Performs Inflation adjustment, according to Urban
            consumer index, and a reference year of 2016

    Output:
        df_prj_core         core / NIH institute, funding project
        df_prj_costs        budget of projects, by year
        df_papers           paper

    """

    df_year = medline.sql_medline_wos("""
        SELECT DISTINCT
            medline.pubmed_id, medline.pubdate_year
        FROM medline
    """)

    df_papers = nih.publnk()
    df_papers = pd.merge(
        df_papers, df_year,
        how='inner')   # our medline only contains the subset with genes

    cols_to_use = [
        'APPLICATION_ID',
        'CORE_PROJECT_NUM',
        'ADMINISTERING_IC',
        'SUBPROJECT_ID',
        'FY',
        'TOTAL_COST',
        'TOTAL_COST_SUB_PROJECT']
    df_prj = nih.prj(cols_to_use)

    # Ignore projects with a no-defiend project numbmer
    f = (df_prj['CORE_PROJECT_NUM'].isnull()
         ) | (df_prj['CORE_PROJECT_NUM'] == 'nan')
    df_prj = df_prj[~f]

    # ignore generic project IDs, that would seem prone to overestimating
    # the budget per gene
    f = df_prj['CORE_PROJECT_NUM'].str.contains('000000$', regex=True)
    df_prj = df_prj[~f]

    # keep only main projects; note that subprojects still use the same
    # core project id, which also is the id used by mapping to pubmed

    f = (df_prj['TOTAL_COST'].notnull()) & (
        df_prj['SUBPROJECT_ID'].isnull()) & (
        df_prj['TOTAL_COST_SUB_PROJECT'].isnull())
    df_prj = df_prj[f]

    # remove columns consisting only of NaN, or no
    f = df_prj.notnull().sum() != 0
    df_prj = df_prj.loc[:, f]

    # separeate administration from costs
    df_prj_core = df_prj[
        ['ADMINISTERING_IC', 'CORE_PROJECT_NUM', 'FY']].drop_duplicates()

    df_prj_costs = df_prj[
        ['CORE_PROJECT_NUM', 'FY', 'TOTAL_COST']].drop_duplicates()
    df_prj_costs = df_prj_costs.groupby(['CORE_PROJECT_NUM', 'FY']).agg(sum)
    df_prj_costs = df_prj_costs.reset_index()

    # perform inflation correction on costs:
    ref_year = 2016
    cpi_adjustement = inflation.consumer_price_adjustment_factor(ref_year)

    df_prj_costs = pd.merge(
        df_prj_costs,
        cpi_adjustement,
        left_on='FY',
        right_on='Year',
        how='left')

    df_prj_costs['TOTAL_COST'] = df_prj_costs[
        'TOTAL_COST'] * df_prj_costs['adjustment_factor']
    df_prj_costs = df_prj_costs.drop('adjustment_factor', axis=1)
    df_prj_costs = df_prj_costs.drop('Year', axis=1)

    # note that 10% some papers are actually mapped to proper grant
    # manual inspection suggests that some have wrong grant category within
    # NIH's list, but appear correct in medline; but cutting the first
    # 3 characters would also lead to ambiguity,... ; one example of mismatch
    # in NIH's description and Medline's is pmid 22078318, which is mapped
    # to an R55 in NIH, but to an R56 in medline (with specific code otherwise
    # being identical)

    # Filter for mutual presence
    f = df_papers['PROJECT_NUMBER'].isin(df_prj_costs['CORE_PROJECT_NUM'])
    df_papers = df_papers[f]

    f = df_prj_costs['CORE_PROJECT_NUM'].isin(df_papers['PROJECT_NUMBER'])
    df_prj_costs = df_prj_costs[f]

    f = df_prj_core['CORE_PROJECT_NUM'].isin(df_papers['PROJECT_NUMBER'])
    df_prj_core = df_prj_core[f]

    renamer = {
        'CORE_PROJECT_NUM': 'project_num',
        'PROJECT_NUMBER': 'project_num',
        'TOTAL_COST': 'budget'}

    def rename(df):
        return df.rename(columns=renamer)

    df_papers = rename(df_papers)
    df_prj_costs = rename(df_prj_costs)
    df_prj_core = rename(df_prj_core)

    return df_prj_core, df_prj_costs, df_papers


def replace_administration_code_by_institute_code(ser):
    """
    NIH grants can also be processed by individual organizations that
    are part of an NIH institute. This function will replace
    division and  office codes by the institutional code.institutional

    Codes have been imported manually from:
    https://www.nlm.nih.gov/bsd/grant_acronym.html

    Input:
        ser         pandas series with NIH 2 letter codes of administrator

    Output:
        ser         pandas series with NIH 2 letter codes of institute

    """
    ser = ser.copy()

    # National Cancer Institute:
    o = ['CO', 'BC', 'CN', 'CB', 'CP', 'CM', 'PC', 'SC']
    f = ser.isin(o)
    ser[f] = 'CA'

    # National Heart, Lung, and Blood Institute
    o = ['HV', 'HB', 'HR', 'HI', 'HO', 'HC']
    f = ser.isin(o)
    ser[f] = 'HL'

    return ser


def get_extended_funding_info(taxon_id, earliest_year, latest_year):
    """
    Function to standardize queries on budget, creates estimate
    of budget per gene

    """

    # INITIALIZATION ###

    # MedLine
    ref_genes = standardizer.reference_genes(taxon_id, 'rpo')
    gene2pubmed = medline.gene2pubmed(
        taxon_id, paper_kind='research', ref_genes=ref_genes)

    df_m = medline.select_medline_wos_records(
        columns_sql='''
                medline.pubmed_id,
                medline.pubdate_year,
                medline.amount_of_authors,
                medline.j_name_s''',
        years_range='all',
        taxon_id=taxon_id,
        kind='research',
        unambiguous=True)

    df_m = df_m[df_m['amount_of_authors'] > 0]   # exclude consortia paper (-1)
    df_m = df_m[['pubmed_id', 'pubdate_year', 'amount_of_authors', 'j_name_s']]

    df_m = df_m[df_m['pubdate_year'] >= earliest_year]
    df_m = df_m[df_m['pubdate_year'] <= latest_year]

    # <========== use later for filtering all!!!!
    _pubmed_articles_in_medline_time_span = set(df_m['pubmed_id'])

    # NIH Exporter
    df_prj_core, df_prj_budget, df_nih_papers = get_paper_funding_through_nih()
    df_nih_papers = df_nih_papers.loc[:, [
        'project_num', 'pubmed_id']]  # skip publication year

    df_prj_core = df_prj_core[df_prj_core['FY'] >= earliest_year]
    df_prj_core = df_prj_core[df_prj_core['FY'] <= latest_year]

    df_prj_budget = df_prj_budget.loc[:, [
        'project_num', 'budget']]   # skip fiscal year
    df_prj_budget = df_prj_budget.groupby('project_num').agg(sum)
    df_prj_budget = df_prj_budget.reset_index()

    # Estimations of costs for non-covered papers ###

    papers_in_nih = len(
        set(df_nih_papers['pubmed_id']).intersection(set(df_m['pubmed_id'])))
    papers_in_medline = len(set(df_m['pubmed_id']))
    multiplier_nih2medline = papers_in_medline / papers_in_nih
    print('Multiplier:', multiplier_nih2medline)

    # Synchronization ###

    # PubMed
    lis = [set(df_nih_papers['pubmed_id']), set(
        df_m['pubmed_id']), set(gene2pubmed['pubmed_id'])]
    pubmed_in_all = set.intersection(*lis)
    print('Amount of MedLine articles:', len(pubmed_in_all))

    gene2pubmed = gene2pubmed[gene2pubmed['pubmed_id'].isin(pubmed_in_all)]
    df_m = df_m[df_m['pubmed_id'].isin(pubmed_in_all)]
    df_nih_papers = df_nih_papers[df_nih_papers['pubmed_id'].isin(
        pubmed_in_all)]

    # Projects
    lis = [set(df_prj_core['project_num']), set(
        df_prj_budget['project_num']), set(df_nih_papers['project_num'])]
    project_in_all = set.intersection(*lis)

    df_prj_core = df_prj_core[df_prj_core['project_num'].isin(project_in_all)]
    df_prj_budget = df_prj_budget[df_prj_budget['project_num'].isin(
        project_in_all)]
    df_nih_papers = df_nih_papers[df_nih_papers['project_num'].isin(
        project_in_all)]

    # Resources per paper per gene

    # amount of publications per project
    papers_per_project = df_nih_papers['project_num'].value_counts()
    # overall budget per project
    budget_per_project = df_prj_budget.set_index('project_num')['budget']
    # budget per paper for each project
    budget_per_paper_per_project = budget_per_project.div(papers_per_project).to_frame(
        'budget_per_paper_per_project').reset_index().rename(columns={'index': 'project_num'})

    budget_per_pubmed_id = pd.merge(budget_per_paper_per_project, df_nih_papers)[
        ['pubmed_id', 'budget_per_paper_per_project']].groupby('pubmed_id').agg(sum).reset_index()

    attention_per_paper = (
        1 / gene2pubmed['pubmed_id'].value_counts()).to_frame('attention_per_gene').reset_index()
    attention_per_paper = attention_per_paper.rename(
        columns={'index': 'pubmed_id'})

    gene2pubmed_plus = pd.merge(gene2pubmed, budget_per_pubmed_id)
    gene2pubmed_plus = pd.merge(gene2pubmed_plus, attention_per_paper)
    gene2pubmed_plus = gene2pubmed_plus.rename(
        columns={'budget_per_paper_per_project': 'budget_for_paper', 'attention_per_gene': 'attention'})
    gene2pubmed_plus.loc[:, 'papers'] = 1

    gene2pubmed_plus['budget_for_attention'] = gene2pubmed_plus['attention'] * \
        gene2pubmed_plus['budget_for_paper']

    master = gene2pubmed_plus[
        ['gene_ncbi', 'budget_for_attention', 'attention', 'papers', 'budget_for_paper']].groupby('gene_ncbi').agg(sum)

    master['budget_by_attention'] = master['budget_for_attention'] / \
        master['attention']
    master['budget_by_papers'] = master['budget_for_paper'] / master['papers']

    gene2pubmed_full = medline.gene2pubmed(
        taxon_id, paper_kind='research', ref_genes=ref_genes)
    gene2pubmed_full = gene2pubmed_full[gene2pubmed_full['pubmed_id'].isin(
        _pubmed_articles_in_medline_time_span)]

    fame_full = nar_attention.count_papers_and_attention(
        ref_genes, gene2pubmed_full)

    n = fame_full.columns
    fame_full.columns = ['full_' + x for x in fame_full.columns]
    master = pd.merge(master.reset_index(), fame_full.reset_index())

    nih_publnk = nih.publnk().drop_duplicates()
    gene2pubmed_all_nih = gene2pubmed_full[gene2pubmed_full['pubmed_id'].isin(
        nih_publnk['pubmed_id'])]

    fame_all_nih = nar_attention.count_papers_and_attention(
        ref_genes, gene2pubmed_all_nih)
    n = fame_all_nih.columns
    fame_all_nih.columns = ['all_nih_' + x for x in fame_all_nih.columns]
    master = pd.merge(master, fame_all_nih.reset_index())

    for x in n:
        master.loc[:, 'non_nih_' + x] = master.loc[:,
                                                   'full_' + x] - master.loc[:, 'all_nih_' + x]

    master = master.set_index('gene_ncbi')

    return master
