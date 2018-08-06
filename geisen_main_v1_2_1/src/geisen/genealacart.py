import glob
import os

import pandas as pd

import geisen.inout as io


# This module contains tools for extracting specific datasets from
# genealacart (which is anticipated in manual data)

# set version(s) of genealacart that shall be used
version_170217 = 'out/genealacart/170217/GeneALaCart-431-170217*'
reference_version = version_170217


def export_selected_genealacart_datasets(patch_absent=False):
    """
    Will load selected datasetes from genealacard and export in a format
    that is consistent with the science of biology project

    Requirements:
        geisen_manual       with genealacart

    Input:
        patch_absent    optional; default: False; If True,
                            absent files will be added (e.g.:
                            if novel features of GeneCards should
                            be extracted)
    """

    p_out = io.get_output_path('genealacart')
    io.ensure_presence_of_directory(p_out)
    if io.check_number_of_files_in_directory(p_out, 'gz') > 0:
        raise EnvironmentError('Output directory needs to be empty')

    def export(df, name):
        o = os.path.join(p_out, 'genealacart_{}.gz'.format(name))

        if patch_absent:
            if not os.path.exists(o):
                df.to_csv(o, index=True, compression='gzip')
                print('Added absent file {}'.format(o))
        else:
            io.ensure_absence_of_file(o)
            df.to_csv(o, index=True, compression='gzip')

    def add_counts_for_absent_reference_genes(df):
        d = pd.merge(
            reference_genes,
            df,
            left_on='gene_ncbi',
            right_index=True,
            how='left')
        d = d.fillna(0)
        d = d.set_index('gene_ncbi')
        d = d.astype(int)

        return d

    # Reference genes: all genes that are in genealacart, and
    # unambiguously map to gene_ncbi gene IDs
    reference_genes = load_genealacart_dataset('ExternalIdentifiers')
    reference_genes = reference_genes[['EntrezGene_x']]
    reference_genes = reference_genes.rename(columns={
        'EntrezGene_x': 'gene_ncbi'})

    print('Start processing ENCODE')
    amount_of_enhancers, tf_by_gene = _get_encode()
    export(amount_of_enhancers, 'encode_amount_of_tfs')
    export(tf_by_gene, 'encode_tfs_by_gene')

    print('Start processing Promoters (ENSRs)')
    amount_of_tfs, tf_by_gene = _get_promoters()
    export(amount_of_enhancers, 'promoters_amount_of_tfs')
    export(tf_by_gene, 'promoters_tfs_by_gene')

    print('Start processing intolerance')
    df_gdi, df_rvis = _get_intolerance()
    export(df_gdi, 'intolerance_gdi')
    export(df_rvis, 'intolerance_rvis')

    print('Start processing selected disease databases')
    dbs = ['DISEASES', 'Orphanet', 'OMIM']

    for disease in dbs:
        amount_of_diseases, df_stack_diseases = _get_disease(disease)
        amount_of_diseases = add_counts_for_absent_reference_genes(
            amount_of_diseases)
        export(amount_of_diseases, '{}_amount'.format(disease.lower()))
        export(df_stack_diseases, '{}_kind'.format(disease.lower()))

    print('Start processing human phenotypes')
    amount_of_phenotypes, df_stack_phenotype = _get_human_phenotype_ontology()
    amount_of_phenotypes = add_counts_for_absent_reference_genes(
        amount_of_phenotypes)
    export(amount_of_phenotypes, 'phenotype_ontology_amount')
    export(df_stack_phenotype, 'phenotype_ontology_kind')

    print('Start processing GIFTS score')
    gifts = _get_gifts()
    export(gifts, 'annotation_range_gifts')


def load_genealacart_dataset(dataset, version=reference_version):
    """
    - Loads genealacart dataset
    - removes ambiguous mappings between
        symbols and entrez gene (Note: which in the 170217 version of
        genelacart only corresponds to the remval of two genes)
    - introdcue gene_ncbi for common mapping to rest of science of
        biology

    Input:

        dataset         str, name of dataset, e.g.: 'Enhancers' note that a
                        given genealacart datast can contain several smaller
                        ones; corresponds to names of files within the zip'ed
                        download of genealacarrt
        version         str, glob-able pattern to folders of indivdual batches
                        within a geisen_manual data path
                        By default: set to latest reviewed version
    """

    df_identifiers = _load_batches(version, 'ExternalIdentifiers')
    allowed_symbols = _get_unambiguous_symbols(df_identifiers)

    df = _load_batches(version, dataset)

    f = df_identifiers['Symbol'].isin(allowed_symbols)
    df_identifiers = df_identifiers.loc[f, :]

    df = pd.merge(
        df_identifiers[['Symbol', 'EntrezGene']],
        df,
        left_on='Symbol',
        right_on='Symbol',
        how='left')

    # Enforce standard scinece of sciene terminology for main ID
    df = df.rename(columns={
        'EntrezGene': 'gene_ncbi'
    })

    # Safety given heterogeous dynamic input data sets
    df = df.drop_duplicates()

    return df


"""
    SUBFUNCTIONS TO ISOLATE PARTS OF GENEALACART DATASET

        Note that this list is not comprehensive, and genealacarts contains
        even more datasets. The list is about properties, which can not
        easily be obtained from other sources, or are nicer post-processed
        in genealacart

        encode          ENCODE enhancers and transcription factors in them
        intolerance     GDI and RVIS
        disease         Disease databases, call with specification (e.g. OMIM)
        human_phenotype_ontology    human phenotypes mapped to genes
        promoters       Transcription factors within promoter elements
                           (appears forked from EBI)

"""


def _get_encode():
    """
    Loads enhancers from ENCODE. And places them in a format indexable
    by gene_ncbi

    Output:
        amount_of_enhancers
        tf_by_gene             occurence matrix with individual TF as columns
    """

    df = load_genealacart_dataset('Enhancers')

    f = (df['Sources'].str.contains('ENCODE')) & \
        (df['Sources'].notnull())

    df = df[f]

    df = df.loc[:, ['gene_ncbi', 'TFBSs', 'TSSdistance']].drop_duplicates()
    df = df.loc[:, ['gene_ncbi', 'TFBSs']]

    amount_of_enhancers = df[
        'gene_ncbi'].value_counts().to_frame(
            'encode_genealacart_amount_of_enhancers')
    amount_of_enhancers.index.name = 'gene_ncbi'

    df = df.loc[df['TFBSs'].notnull(), :]  # no annotated TFs
    dfn = _split_text_to_multiple_rows(df, 'TFBSs', '\|\|')

    tf_by_gene = dfn.groupby(['gene_ncbi', 'TFBSs']).size()
    tf_by_gene = tf_by_gene.reset_index()
    tf_by_gene = tf_by_gene.pivot(index='gene_ncbi', columns='TFBSs', values=0)
    tf_by_gene = tf_by_gene.fillna(0)
    tf_by_gene.columns = ['encode_genealacart__TFs_in_enhancers: {}'.format(
        j) for j in tf_by_gene.columns]

    return amount_of_enhancers, tf_by_gene


def _get_gifts():
    """
    Obtain GIFT score, which indicates fraction of databases,
    in which a gene appears. Used by Genealacart as
    a proxy of the general annotation density for genes.

    Output:
        df_gifts       GIFTS score
    """

    df_gifts = load_genealacart_dataset('Gene')
    df_gifts = df_gifts[['Gifts']].rename(columns={'Gifts': 'gifts'})
    return df_gifts


def _get_promoters():
    """
    Loads transcription factors within promoters.
    And places them in a format indexableby gene_ncbi.
    The source datbases appears Ensemble's regulatory elements, and
    thus this function only considers entries with an ENSR ID

    Output:
        amount_of_promoters
        tf_by_gene             occurence matrix with individual TF as columns
    """

    df = load_genealacart_dataset('Promoters')

    f = df['ENSRs'].notnull()
    df = df[f]

    df = df.loc[:, ['gene_ncbi', 'TFBSs', 'TSSdistance']].drop_duplicates()
    df = df.loc[:, ['gene_ncbi', 'TFBSs']]

    amount_of_promoters = df['gene_ncbi'].value_counts().to_frame(
        'promoters_genealacart_amount_of_promoters')
    amount_of_promoters.index.name = 'gene_ncbi'

    df = df.loc[df['TFBSs'].notnull(), :]   # no annotated TF
    dfn = _split_text_to_multiple_rows(df, 'TFBSs', '\|\|')

    tf_by_gene = dfn.groupby(['gene_ncbi', 'TFBSs']).size()
    tf_by_gene = tf_by_gene.reset_index()
    tf_by_gene = tf_by_gene.pivot(index='gene_ncbi', columns='TFBSs', values=0)
    tf_by_gene = tf_by_gene.fillna(0)

    tf_by_gene.columns = ['promoters_genealacart__TFs_in_promoters: {}'.format(
        j) for j in tf_by_gene.columns]

    return amount_of_promoters, tf_by_gene


def _get_intolerance():
    """
    Loads intolerance metrics. And places them in a format indexable
    by gene_ncbi

    Output:
        df_gdi      gdi score
        df_rvis     rvis score
    """

    df = load_genealacart_dataset('Intolerance')

    df_gdi = df[['gene_ncbi', 'GDI']].dropna()
    df_gdi = df_gdi.set_index('gene_ncbi')
    df_gdi.loc[:, 'gdi_quantity1_genealacart'] = df_gdi.loc[
        :, 'GDI'].str.extract('^(.*)\:', expand=False)
    df_gdi.loc[:, 'gdi_quantity2_genealacart'] = df_gdi.loc[
        :, 'GDI'].str.extract('\:(.*)$', expand=False)
    df_gdi.loc[:, 'gdi_quantity1_genealacart'] = df_gdi.loc[
        :, 'gdi_quantity1_genealacart'].astype(float)
    df_gdi.loc[:, 'gdi_quantity2_genealacart'] = df_gdi.loc[
        :, 'gdi_quantity2_genealacart'].astype(float)
    df_gdi = df_gdi.drop('GDI', axis=1)

    df_rvis = df[['gene_ncbi', 'RVIS']].dropna()
    df_rvis = df_rvis.rename(columns={'RVIS': 'rvis_genealacart'})
    df_rvis = df_rvis.set_index('gene_ncbi')

    return df_gdi, df_rvis


def _get_disease(source):
    """
    Loads disease metrics. And places them in a format indexable
    by gene_ncbi. Will not consider proximity of diseases: Any
    disease with a unique name will be considered as a separate entry

    Input:
        source    str, source of disease information that
                    will be extracted from malacards
                    e.g.: 'DISEASES' or 'Orphanet' or 'OMIM'

    Output:
        amount_of_diseases      dataframe, with counts of diseases
        df_stack_diseases       dataframe, in stacked format
    """

    df = load_genealacart_dataset('MalaCardsDisorders')

    df['Name'].value_counts().to_frame(
        'Counts').to_csv('tmp_diseases.csv', index=True)

    f = (df['Sources'].str.contains(source)) & (df['Sources'].notnull())

    df = df.loc[f, :]
    df = df.loc[:, ['gene_ncbi', 'Name']].dropna().drop_duplicates()
    df = df.rename(columns={'Name': 'disease_name'})
    df_stack_diseases = df.set_index('gene_ncbi')

    amount_of_diseases = df['gene_ncbi'].value_counts().to_frame(
        '{}_malacard_amount_of_diseases'.format(source.lower()))
    amount_of_diseases.index.name = 'gene_ncbi'

    return amount_of_diseases, df_stack_diseases


def _get_human_phenotype_ontology():
    """
    Loads human phenotype ontology. And places them in a format indexable
    by gene_ncbi. Will not consider proximity of diseases: Any
    disease with a unique name will be considered as a separate entry

    Output:
        amount_of_phenotypes    dataframe, with counts of human phenotypes
        df_stack_phenotype      dataframe, in stacked format
                                columns: human phenotype ID and name
    """

    df = load_genealacart_dataset('HumanPhenotypeOntology')

    dfn = df.loc[:, ['gene_ncbi', 'AlternativeId', 'Name']].drop_duplicates()
    dfn = dfn.rename(columns={
        'AlternativeId': 'human_phenotype_id',
        'Name': 'human_phenotype_name'})
    dfn = dfn.dropna(subset=['human_phenotype_id'])
    dfn = dfn.set_index('gene_ncbi')
    dfn.columns = ['human_phenotype_genealacart: {}'.format(
        j) for j in dfn.columns]

    df_stack_phenotype = dfn

    nn = dfn.reset_index()
    amount_of_phenotypes = nn['gene_ncbi'].value_counts().to_frame(
        'malacard_amount_of_phenotypes')
    amount_of_phenotypes.index.name = 'gene_ncbi'

    return amount_of_phenotypes, df_stack_phenotype


"""
    GENERAL SUPPORT FUNCTIONS

        _load_batches
        _get_unambiguous_symbols
        _split_text_to_multiple_rows
"""


def _load_batches(path_pattern, dataset):
        p_scheme = io.get_geisen_manual_data_path(path_pattern)
        agg = []
        for p in glob.glob(p_scheme):
            p = os.path.join(p, '{}.txt'.format(dataset))
            agg.append(pd.read_table(p))
        df = pd.concat(agg)
        df = df.drop_duplicates()

        # place symbols lower-case as genealacart does not appear
        # to distinguish internally
        n = ['InputTerm', 'Symbol']
        for t in n:
            df.loc[:, t] = df.loc[:, t].str.lower()

        return df


def _get_unambiguous_symbols(df_identifiers):
        """
        Obtains the entries which have an unambiguous
        Input Term, Symbol and Entrez Gene. Will return
        unambigous Symols (which are used internally by
        GeneALaCart)

        Note that input symbols can be duplicate, if the same
        input maps to several Symbols in Genealacart

        Input:
            df_identifiers  DataFrame with identifiers containing
                            ['InputTerm', 'Symbol', 'EntrezGene']
                            e.g.: df_identifiers = _load_batches(
                                version, 'ExternalIdentifiers')

        Output:
            unambiguous_symbols  list; all symbols that can
                                    be mapped unambiguously to
                                    gene_ncbi within genealacart
        """

        df = df_identifiers

        f = df['EntrezGene'].notnull()
        dff = df.loc[f, :]

        terms_that_should_be_unique = ['InputTerm', 'Symbol', 'EntrezGene']

        dff = dff.loc[:, terms_that_should_be_unique]
        dff = dff.drop_duplicates()

        ambiguous = dict()
        for t in terms_that_should_be_unique:
            c = dff.loc[:, t].value_counts()
            c = list(c[c > 1].index)
            ambiguous[t] = c

        for k, v in ambiguous.items():
            num_ambiguous = len(v)
            if num_ambiguous > 0:
                print(k, 'has', num_ambiguous, 'ambiguous entries')
                f = dff.loc[:, k].isin(v)
                dff = dff.loc[~f, :]

        unambiguous_symbols = dff.loc[:, 'Symbol'].unique()

        return unambiguous_symbols


def _split_text_to_multiple_rows(df, column, separator):
    """
    Creates database where entries within one column are separated
    to multiple rows

    Input:
        df          DataFrame
        column      str; column that needs separation
        separator   regulatar expression

    Output:
        df_stacked  DataFrame, with records split to separate rows
    """

    # Separate rows, that should be processed (to save overhead)
    f = df[column].str.contains(separator)
    dff_s = df.loc[f, :]    # _s  separate
    dff_ns = df.loc[~f, :]  # _ns not separate
    orig_column_order = df.columns  # track original order of columns

    # Separate records and place them in separate rows
    s = dff_s[column].str.split(separator).apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = column
    del dff_s[column]
    dff_s = dff_s.join(s)

    # Recreate DataFrame in original input format
    dff_s = dff_s.reindex_axis(orig_column_order, axis=1)
    df = pd.concat([dff_s, dff_ns])
    df_stacked = df

    return df_stacked
