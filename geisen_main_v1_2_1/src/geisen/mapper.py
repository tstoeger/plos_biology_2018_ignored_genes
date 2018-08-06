import os

import numpy as np
import pandas as pd

import geisen.inout as io


def gene_ensembl_2_gene_ncbi_unambiguously(df, taxon_id):
    """
    Maps ensembl gene ID to NCBI (Entrez) gene IDs. Will only
    consider unambiguos 1:1 mappings of ensembl and entrez gene IDs.

    Although ncbi and ensg have a working project on creating
    a uniform mapping for mouse and humans, the mapping is not
    necessarily unambiguous; Although there are different mappers
    available from different organisations, different
    organizations have different mapping schemes.

    The present mapper will use NIH's gene_info. Note that
    this can be different from EBI Biomart (which appers to map
    by overlap of any sequence)

    If the mapping of genes is not 1:1 within gene_info, those
    genes will be ignored.

    Note that for some taxa ensembl does not carry unique identifiers,
    but external databases, which are also listed as other
    databases in NIH (e.g.: flybase or wormbase IDs). -> when moving
    to additional taxa, one may need to implement taxon specific
    exteranal references (that would also be used by ensembl)

    Furhter note: in contrast to "Science of Biology v0.1" this
    function uses NIH's gene_info rather NIH's gene2ensembl, as the
    former covers more taxa

    Input:
        df      dataframe with gene_ensembl

    Output:
        dfm     dataframe with gene_ncbi as index

    """

    id_name = 'gene_ensembl'  # Science of Biology nomeclature
    is_column, is_index = _check_for_presence(
        df,
        id_name,
        require_presence=True)

    # Construct Mapper from gene_info
    p_mapper = io.get_output_path(
        'ncbi/gene_info/gene_info_taxon_{}.gz'.format(
            taxon_id))
    if not os.path.exists(p_mapper):
        raise EnvironmentError(
            'gene_ensembl_2_gene_ncbi_unambiguously()'
            ' requires gene_info')
    m = pd.read_csv(
        p_mapper,
        usecols=[
            'gene_ncbi',
            'dbXrefs']).drop_duplicates()
    m = m.set_index('gene_ncbi')

    # get taxon specific pattern for extracting ensembl id
    # for some taxa, ensembl has interited from other databases
    if taxon_id in [6239]:
        p = 'WormBase:([A-Za-z0-9]*)'
    elif taxon_id in [7227]:
        p = 'FLYBASE:([A-Za-z0-9]*)'
    else:  # Default
        p = 'Ensembl:([A-Z0-9]*)'

    m = m.loc[:, 'dbXrefs'].str.extractall(p)
    m = m.rename(columns={0: 'gene_ensembl'})
    m = m.reset_index()

    mapper = m[['gene_ncbi', 'gene_ensembl']]

    # Tidy mapper: only consider unambiguous ones
    forbidden_ncbi = _get_duplicates(mapper['gene_ncbi'])
    forbidden_ensg = _get_duplicates(mapper['gene_ensembl'])
    f = (
        (~mapper['gene_ncbi'].isin(forbidden_ncbi)) &
        (~mapper['gene_ensembl'].isin(forbidden_ensg))
    )
    mapper = mapper.loc[f, :]

    if is_index:
        df = df.reset_index()

    dfm = pd.merge(
        df,
        mapper,
        left_on=id_name,
        right_on=id_name,
        how='inner')
    dfm = dfm.drop(id_name, axis=1)

    dfm = dfm.set_index('gene_ncbi')

    return dfm


def locustag_2_gene_ncbi_unambiguously(df, taxon_id):
    """
    Maps locus tag to NCBI (Entrez) gene IDs. Will only
    consider unambiguos 1:1 mappings.

    Input:
        df      dataframe with LocusTag

    Output:
        dfm     dataframe with gene_ncbi as index

    """

    id_name = 'LocusTag'
    is_column, is_index = _check_for_presence(
        df,
        id_name,
        require_presence=True)

    # Construct Mapper from gene_info
    p_mapper = io.get_output_path(
        'ncbi/gene_info/gene_info_taxon_{}.gz'.format(
            taxon_id))
    if not os.path.exists(p_mapper):
        raise EnvironmentError(
            'locustag_2_gene_ncbi_unambiguously()'
            ' requires gene_info')
    mapper = pd.read_csv(
        p_mapper,
        usecols=[
            'gene_ncbi',
            'LocusTag']).drop_duplicates()

    # Tidy mapper: only consider unambiguous ones
    forbidden_ncbi = _get_duplicates(mapper['gene_ncbi'])
    forbidden_locus = _get_duplicates(mapper['LocusTag'])
    f = (
        (~mapper['gene_ncbi'].isin(forbidden_ncbi)) &
        (~mapper['LocusTag'].isin(forbidden_locus))
    )
    mapper = mapper.loc[f, :]

    if is_index:
        df = df.reset_index()

    dfm = pd.merge(
        df,
        mapper,
        left_on=id_name,
        right_on=id_name,
        how='inner')
    dfm = dfm.drop(id_name, axis=1)

    dfm = dfm.set_index('gene_ncbi')

    return dfm


def ncbi_taxon_2_uniprot_taxon(taxon_id_ncbi):
    """
    Retreives a manually curated lookup between ncbi taxonomy IDs
    and taxon names as used by unprot (swissprot and trembl).

    Note that mapping might not always be unambiguous, for instance:
    uniprot contains several subtypes for HIV, whereas the taxon id
    from ncbi refers to a collection of strains. -> In such cases
    the curator (Thomas Stoeger) would attempt find an inclusive name
    within uniprot

    Input:
        taxon_id_ncbi       int; ncbi taxonomy id (e.g.: 9606 for homo sapiens)

    Output:
        taxon_uniprot       str; matching taxonomy description used in uniprot

    """

    dict_species = {
        9606: 'Homo sapiens',
        10090: 'Mus musculus',
        10116: 'Rattus norvegicus',
        7227: 'Drosophila melanogaster',
        511145: 'Escherichia coli (strain K12)',
        559292: 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)',
        3702: 'Arabidopsis thaliana',
        7955: 'Danio rerio',
        9913: 'Bos taurus',
        11676: 'Human immunodeficiency virus',  # includes several subtypes
        9031: 'Gallus gallus',
        6239: 'Caenorhabditis elegans',
        9823: 'Sus scrofa',
        8355: 'Xenopus laevis',
        284812: 'Schizosaccharomyces pombe (strain 972 / ATCC 24843)',
        386585: 'Escherichia coli O139:H28 (strain E24377A / ETEC)',
        9615: 'Canis lupus familiaris',
        224308: 'Bacillus subtilis (strain 168)',
        9986: 'Oryctolagus cuniculus',
        83332: 'Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)',
    }

    if taxon_id_ncbi in list(dict_species.keys()):
        taxon_uniprot = dict_species[taxon_id_ncbi]
    else:
        raise ValueError(
            'Did not find lookup for uniprot taxon ID'
            'Possibly the taxon has not been included in the lookup.'
            'Lookup has been gernearted by manual curation')

    return taxon_uniprot


def rna_ensembl_2_gene_ncbi(df, how):
    """
    Maps ensembl transcript ID to NCBI (Entrez) gene IDs.

    The present mapper will use gene2ensembl from NIH. Note that
    this can be different from EBI Biomart (which appers to map
    by overlap of any sequence)

    Input:
        df      dataframe with rna_ensembl
        how     str, method for aggregation (e.g.: median)

    Output:
        dfm     dataframe with gene_ncbi as index

    """

    id_name = 'rna_ensembl'  # Science of Biology nomeclature
    is_column, is_index = _check_for_presence(
        df,
        id_name,
        require_presence=True)

    p_mapper = io.get_output_path('ncbi/gene2ensembl.gz')
    if not os.path.exists(p_mapper):
        raise EnvironmentError(
            'rna_ensembl_2_gene_ncbi() requires gene2ensembl')
    mapper = pd.read_csv(
        p_mapper,
        usecols=[
            'gene_ncbi',
            'rna_ensembl']).drop_duplicates()

    if is_index:
        df = df.reset_index()

    dfm = pd.merge(
        df,
        mapper,
        left_on=id_name,
        right_on=id_name,
        how='inner')
    dfm = dfm.drop(id_name, axis=1)

    df_fused = _group_aggregate_to_gene_ncbi(dfm, how)

    return df_fused


def symbol_2_gene_ncbi(df, taxon_id, how):
    """
    - Mappes a dataframe with gene symbols IDs to gene_ncbi
    - Places gene_ncbi as the index
    - Only returns genes that could be mapped (inner join)
    - Aggregates according to how (e.g.: median)

    Input:
        df    DataFrame, with symbol_ncbi (or as fallback: symbol_ambiguous)
        taxon_id  int with ncbi taxonomy ID; Required as the same symbols
                    are often used for homologs of different taxa
        how  str, eg.: median, (or substitute to ignore aggregation)

    """

    id_name = 'symbol_ncbi'  # Science of Biology nomeclature
    is_column, is_index = _check_for_presence(
        df,
        id_name,
        require_presence=False)

    if (not(is_column)) & (not(is_index)):
            id_name = 'symbol_ambiguous'  # Science of Biology fall-back
            is_column, is_index = _check_for_presence(
                df,
                id_name,
                require_presence=True)   # throw error if no match

    p_mapper = io.get_output_path(
        'ncbi/gene_info/gene_info_taxon_{}.gz'.format(taxon_id))

    if not os.path.exists(p_mapper):
        raise EnvironmentError(
            'symbol_2_gene_ncbi() requires taxon specific gene_info')
    mapper = pd.read_csv(
        p_mapper,
        usecols=['gene_ncbi', 'symbol_ncbi'])

    if is_index:
        df = df.reset_index()

    dfm = pd.merge(
        df,
        mapper,
        left_on=id_name,
        right_on='symbol_ncbi',
        how='inner')
    dfm = dfm.drop(id_name, axis=1)

    if id_name != 'symbol_ncbi':
        dfm = dfm.drop('symbol_ncbi', axis=1)

    if how is not 'substitute':
        df_fused = _group_aggregate_to_gene_ncbi(dfm, how)
    else:
        df_fused = dfm
        if is_index:
            df_fused = df_fused.set_index(
                'gene_ncbi',
                verify_integrity=True)

    return df_fused


def uniprot_protein_2_gene_ncbi(df, how):
    """
    - Mappes a dataframe with uniprot_protein IDs to gene_ncbi
    - Places gene_ncbi as the index
    - Only returns genes that could be mapped (inner join)
    - Aggregates according to how (e.g.: median)

    Input:
        df    DataFrame, with protein_uniprot
        how  str, eg.: median

    """

    id_name = 'protein_uniprot'  # Science of Biology nomeclature
    is_column, is_index = _check_for_presence(
        df,
        id_name,
        require_presence=True)

    p_mapper = io.get_output_path('uniprot/uniprot_id_mapper.h5')
    if not os.path.exists(p_mapper):
        raise EnvironmentError(
            'uniprot_protein_2_gene_ncbi() requires uniprot_id_mapper')
    mapper = pd.read_hdf(
        p_mapper,
        'table',
        columns=['protein_uniprot', 'gene_ncbi'])

    if is_index:
        df = df.reset_index()

    dfm = pd.merge(
        df,
        mapper,
        left_on='protein_uniprot',
        right_on='protein_uniprot',
        how='inner')
    dfm = dfm.drop(id_name, axis=1)

    df_fused = _group_aggregate_to_gene_ncbi(dfm, how)

    return df_fused


def _group_aggregate_to_gene_ncbi(df, how):
    """
    Groups DataFrame and returns dataframe grouped by
    gene_ncbi and aggregated accordign to how

    Input:
        df      DataFrame, with gene_ncbi
        how    Method for aggregating, e.g.: median, nanmedian

    """

    if how == 'median':
        def fn(x):
            return np.median(x)
    elif how == 'nanmedian':
        def fn(x):
            return np.nanmedian(x)
    elif how == 'any':
        def fn(x):
            return any(x)
    else:
        raise ValueError(
            'Given value for how not supported.'
            'Presently supported options are:'
            'median, nanmedian, any')

    g = df.groupby('gene_ncbi')
    df_mapped = g.agg(fn)

    if (df_mapped.shape[1] + 1) != (df.shape[1]):  # +1 for index
        raise ValueError(
            'Some columns appear lost while aggregating to gene_ncbi'
            'Potential reason: processing function not suitable for'
            'format of at least one column.')

    return df_mapped


def _check_for_presence(df, id_name, require_presence=False):
    """
    Checks for presence of id_name in index or column of
    a given dataframe df. Will throw erros if the id is amiguous
    (appears in index and column) or absent

    Input:
        df      DataFrame, possibly with one column or index as identifier
        id_name str; name of id that should be present
        require_presence    bool; default: False; will throw error,
                            if id_name is absent

    Output:
        is_column   bool; True if id_name occurs in column
        is_index    bool; True if id_name occurs in index
    """

    is_column = False
    is_index = False

    if id_name in df.columns:
        is_column = True

    if id_name == df.index.name:
        is_index = True

    if require_presence:  # optional input
        if (not is_column) & (not is_index):
            raise ValueError('No {} defined in dataframe to be mapped.')

    if (is_column) & (is_index):
        raise ValueError('Found {} in index and in column.')

    return is_column, is_index


def _get_duplicates(ser):
    """
    obtains the duplicate values within a series

    Input:
        ser     Series

    Output:
        duplicates  list; of values which occur more than once
    """

    c = ser.value_counts()
    duplicates = list(c[c > 1].index)
    return duplicates
