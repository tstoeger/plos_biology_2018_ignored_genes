import os

import pandas as pd

from access_science_shared import inout, queries, standardizer

"""
Functions for accessing MedLine:

gene2pubmed
select_medline_records
select_medline_wos_records
sql_medline_wos             direct access at sql level

"""


"""
    MAIN FUNCTIONS
"""


def gene2pubmed(taxon_id=None, usecols=None, paper_kind=None, ref_genes=None):
    """
    Loads gene2pubmed from NIH; Will only return
    non-duplicated data in casse that columns are
    specified.

    Input:
        taxon_id    int, or 'all'
        usecols     optional, list of columns to be loaded
        paper_kind  optional; filter for  articles, e.g:
                        'research' to filter for papers in
                        Medline, where meta data suggests
                        that it is a research paper
        ref_genes   optional; filter for genes in ref_genes

    Output:
        gene2pubmed df

    """

    p = inout.get_path(
        'geisen',
        'ncbi/gene2pubmed.h5')

    if os.path.exists(p) is False:
        raise EnvironmentError(
            'Did not find gene2pubmed')

    def load_all_taxa():
        df = pd.read_hdf(p, 'table')
        return df

    def load_taxon(usecols):
        q = 'taxon_ncbi=={}'.format(taxon_id)
        df = pd.read_hdf(p, 'table', where=q)
        return df

    # Implement input specific behavior of gene_info
    if taxon_id == 'all':
        df = load_all_taxa()
    elif isinstance(taxon_id, int):
        df = load_taxon(taxon_id)
    else:
        raise EnvironmentError(
            'Did not recognize format of taxon_id')

    if usecols is not None:
        df = df.loc[:, usecols]
        df = df.drop_duplicates()

    if paper_kind is not None:
        df = standardizer.filter_by_paper_kind(df, paper_kind)

    if ref_genes is not None:
        df = df[df['gene_ncbi'].isin(ref_genes)]

    return df


def select_medline_records(
        columns_sql='''
            medline.pubmed_id,
            medline.pubdate_year,
            medline.amount_of_authors,
            medline.j_name_s''',
        taxon_id=None,
        kind=None):
    """
    Queries Medline; not that this function does not imply any cross-
    dependency or mappability with web of sciece; for altearntive see:
    select_confident_medline_wos_research_records()

    Input:
        columns_sql     string about columns that should be loaded
        taxon_id        int, optional; ncbi_taxon id that should be
                            loaded; leave as None (or 'all') to load all taxa
        kind            defines default filter, e.g.: 'research' to filter
                            for high-conficence resarch publicatoins

    Output:
        df              DataFrame corresponding to query, with sorted
                            pubmed_id as index

    """

    if taxon_id == 'all':
        taxon_id = None

    if taxon_id is not None:
        taxon_sql = 'AND pubmed2taxon.taxon_ncbi = {}'.format(
            int(taxon_id))
    else:
        taxon_sql = ''

    sql = '''
    SELECT DISTINCT
        {}
    FROM pubmed2taxon, medline
        WHERE
                pubmed2taxon.pubmed_id = medline.pubmed_id
            {}
            AND
                {}
    ORDER BY medline.pubmed_id ASC;
    '''.format(
        columns_sql,
        taxon_sql,
        _prefixed_medline_wos_filter(kind))

    df = sql_medline_wos(sql)

    return df


def select_medline_wos_records(
        columns_sql='''
            medline.pubmed_id,
            medline.pubdate_year,
            medline.amount_of_authors,
            medline.j_name_s''',
        years_range=None,
        taxon_id=None,
        kind=None,
        unambiguous=True):
    """
    Queries Medline-WoS. Note that the dependency on mutual presence in
    MedLine and WoS might lead to a slightly reduced number of papers
    compared to queries through select_medline_records() which only queries
    MedLine

    Input:
        columns_sql     string about columns that should be loaded
        years_range     tuple or 'all, optional; range indicating years
                            for which citations; e.g.: (1970, 2017);
                            If 'all is specified, all years will be loaded'
                            (note: in addition 'all' will load 'index' and
                            'ut')
        taxon_id        int, optional; ncbi_taxon id that should be
                            loaded; leave as None (or 'all') to load all taxa
        kind            str; defines filter for publications, e.g.:
                            'research'
        unambiguous     optional; set to True to require unambigous mapping
                            between MedLine and Web of Science, and mapping
                            score above 95

    Output:
        df              DataFrame corresponding to query, with sorted
                            pubmed_id as index

    """

    def _get_citation_years(tup):
        """
        tuple to define range

        Input:
            (start_year, first_excluded year)
        Output:
            concatanation of strings corresponding to column names
            with years
        """
        vec = range(tup[0], tup[1])
        sql_sub_string = ''
        for v in vec:
            label = 'citations.citation_year_{}, '.format(int(v))
            sql_sub_string = sql_sub_string + label
        sql_sub_string = sql_sub_string[:-2]
        return sql_sub_string

    if years_range is 'all':
        # use generic filter (which unfortunately also returns ut and index)
        columns_sql = columns_sql + ', citations.*'
    elif years_range is not None:
        columns_sql = columns_sql + ',' + _get_citation_years(years_range)

    if taxon_id == 'all':
        taxon_id = None

    if taxon_id is not None:
        taxon_sql = 'AND pubmed2taxon.taxon_ncbi = {}'.format(
            int(taxon_id))
    else:
        taxon_sql = ''

    if unambiguous:
        unambiguous_sql = 'AND ut2pmid.score >= 95 ' + \
                          'AND ut2pmid.ambiguous_pmid_2_ut = 0'
    else:
        unambiguous_sql = ''

    sql = '''
    SELECT DISTINCT
        {}
    FROM pubmed2taxon, medline, ut2pmid, citations
        WHERE
                pubmed2taxon.pubmed_id = medline.pubmed_id
            AND
                pubmed2taxon.pubmed_id = ut2pmid.pubmed_id
            AND
                medline.pubmed_id = ut2pmid.pubmed_id
            AND
                ut2pmid.ut = citations.ut
            {}
            AND
                {}
            {}
    ORDER BY medline.pubmed_id ASC;
    '''.format(
        columns_sql,
        taxon_sql,
        _prefixed_medline_wos_filter(kind),
        unambiguous_sql)

    df = sql_medline_wos(sql)

    return df


def sql_medline_wos(sql):
    """
    Queries medline_wos database with the provided SQL.
    Assumes that the database is static, and will
    thus use local on-disk caching to speed
    up large queries (e.g.: on all human publications)

    By default the cacheing function will delete cached
    results which are older than one week

    Input:
        sql     str
    Output:
        df      dataframe as queried from database
    """

    df = queries.sql_w_cache(
        sql,
        p_database=os.path.join(
            inout.get_path('medline_wos', 'db/medline_wos.sqlite')),
        refresh=False)

    return df


"""
    SUPPORT FUNCTIONS
"""


def _prefixed_medline_wos_filter(paper_kind):
    """
    Gets premade filters for medline_wos queries

    Possible Inputs:
        paper_kind      Research articles, no review,
                            no paper that appears
                            in a journal with the majority of
                            reviews

    Output:
        sql_string      with filtering options

    """

    filters = dict()


#     research:    MedLine Publication class has
#     to be one of mulitple research classes (e.g.: clinical trial or twin
#     study or journal article,...) and must not be a review; Along these
#     lines publications from Review journals (more than 50% of publications
#     are reviews) will be ignored even if the research class indicates
#     othererwise (This is motivated by manual lookup of articles of these
#     journals, which tend to be very small hypothesis / opinion articles,
#     that are very different from classical research articles, but might
#     affect composite readouts as they usuaully only have one author).
#     Will always return pubmed_id (from MedLine) as index
    filters['research'] = """
            (
                medline.atype_case_reports = 1
            OR
                medline.atype_classical_article = 1
            OR
                medline.atype_clinical_trial = 1
            OR
                medline.atype_clinical_trial_phase_i = 1
            OR
                medline.atype_clinical_trial_phase_ii = 1
            OR
                medline.atype_clinical_trial_phase_iii = 1
            OR
                medline.atype_clinical_trial_phase_iv = 1
            OR
                medline.atype_comparative_study = 1
            OR
                medline.atype_historical_article = 1
            OR
                medline.atype_journal_article = 1
            OR
                medline.atype_meta_analysis = 1
            OR
                medline.atype_multicenter_study = 1
            OR
                medline.atype_observational_study = 1
            OR
                medline.atype_randomized_controlled_trial = 1
            OR
                medline.atype_twin_study = 1
            OR
                medline.atype_validation_studies = 1
            )
        AND
            medline.atype_review = 0
        AND
            medline.j_name_s IN
                (
                SELECT j_name_s
                FROM type_ratio
                WHERE type_ratio.atype_review < 0.5
                )
        """

    return filters[paper_kind]


"""
    LEGACY FUNCTIONS
"""


def get_pubmed_ids_passing_medline_filter(paper_kind):
    """
    Gets list of pubmed ids, as present in MedLine, and
    fullfilling the criteria as specified in paper_kind

    Input:
        paper_kind      string; defining filter; e.g.:
                            'research'

    Output:
        pubmed_ids      list of pubmed_ids tha
    """

    print("""
        Using get_pubmed_ids_passing_medline_filter().
        This function will become deprecated in future.
        Use select_medline_records instead.""")

    sql = '''
        SELECT DISTINCT
            medline.pubmed_id
        FROM medline
            WHERE
                {}
        ORDER BY medline.pubmed_id ASC;
        '''.format(_prefixed_medline_wos_filter(paper_kind))

    pubmed_ids = sql_medline_wos(sql)['pubmed_id']

    return pubmed_ids


# def select_confident_medline_wos_research_records(
#         columns_sql='''
#             medline.pubdate_year,
#             medline.amount_of_authors,
#             medline.j_name_s''',
#         years_range=None,
#         taxon_id=None):

#     """
#     Queries Medline-WoS with default filters for high confidenc research
#     records; More specifically: unambiguous matching of MedLine and WoS
#     with a matching score of at least 95; MedLine Publication class has
#     to be one of mulitple research classes (e.g.: clinical trial or twin
#     study or journal article,...) and must not be a review; Along these
#     lines publications from Review journals (more than 50% of publications
#     are reviews) will be ignored even if the research class indicates
#     othererwise (This is motivated by manual lookup of articles of these
#     journals, which tend to be very small hypothesis / opinion articles,
#     that are very different from classical research articles, but might
#     affect composite readouts as they usuaully only have one author).
#     Will always return pubmed_id (from MedLine) as index

#     Input:
#         columns_sql     string about columns that should be loaded
#         years_range     tuple or 'all, optional; range indicating years
#                             for which citations; e.g.: (17970, 2017);
#                             If 'all is specified, all years will be loaded'
#                             (note: in addition 'all' will load 'index' and
#                             'ut')
#         taxon_id        int, optional; ncbi_taxon id that should be
#                             loaded; leave as None to load all taxa

#     Output:
#         df              DataFrame corresponding to query, with sorted
#                             pubmed_id as index

#     """

#     print("""
#         Using select_confident_medline_wos_research_records().
#         This function will become deprecated in future.
#         Use select_medline_wos_records instead.""")

#     df = select_medline_wos_records(
#         columns_sql,
#         years_range,
#         taxon_id,
#         kind='research',
#         unambiguous=True)

#     return df
