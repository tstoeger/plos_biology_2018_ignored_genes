import glob
import os
import re

import pandas as pd


from access_biology_data import meta
from access_science_shared import inout, utils


def abundance(dataset):
    """
    Loads abundance measurement from GTX gene expression atlas.
    Will try to find the sample annotation from a configuration file
    placed within the dataset folder

    Input:
        dataset, str, code of dataset: e.g: E-MTAB-2770

    Output:
        df, count tables

    """

    p_experiment = inout.get_path('gtx_atlas', dataset)
    if not os.path.exists(p_experiment):
        if not os.path.exists(inout.get_path('gtx_atlas')):
            raise EnvironmentError(
                'Could not find gtx atlas')
        else:
            raise EnvironmentError(
                'Could not find gtx dataset {}'.format(
                    dataset))

    fi = '{}.tsv'.format(dataset)
    p_abundance = os.path.join(
        p_experiment, fi)
    if not os.path.exists(p_abundance):
        raise EnvironmentError(
            'Could not find the anticipated abundance table {}'.format(fi))
    else:
        df = pd.read_table(p_abundance)

    fi = '{}-configuration.xml'.format(dataset)
    p_sample_annotation = os.path.join(
        p_experiment, fi)
    if not os.path.exists(p_sample_annotation):
        raise EnvironmentError(
            'Could not find the anticipated sampe annotation {}'.format(fi))

    pattern = '^.*assay_group id=\"(.*)\" label=\"(.*)\">.*$'
    group_to_annotation = dict()
    with open(p_sample_annotation) as configuration:
        for li in configuration:
            if '<assay_group id' in li:
                catch = re.search(pattern, li)
                group_key = catch.group(1)
                label = catch.group(2)
                group_to_annotation[group_key] = label
    df = df.rename(columns=group_to_annotation)

    df = df.drop('Gene Name', axis=1)
    df = df.rename(columns={'Gene ID': 'gene_ensembl'})
    df = df.set_index('gene_ensembl')

    def get_unique_counts(x):
        unique_counts = list(set(map(float, x.split(','))))
        return unique_counts

    df_uniques = df.applymap(lambda x: get_unique_counts(x))

    df_first = df_uniques.applymap(lambda x: x[0])
    df_last = df_uniques.applymap(lambda x: x[-1])

    if any((df_first - df_last).apply(any)):
        raise ValueError(
            'Simplifiying assumpiton not met' +
            ' Need to implement slower more sophisticate pooling')
    else:
        df = df_first

    # Add GTX atlas data-set specific name in front of column name
    df.columns = ['{}_{}'.format(dataset, x) for x in df.columns]

    return df


def differential_gene_expression_of_taxon(taxon_id):
    """
    Obtains differential gene expession experiments of one taxon;

    will use comparison_key as a unique key for a
    given combination of experiemnt, analysis, and group comparison

    If genes are covered multiple times (e.g.: by microarray), the
    most significantly differentially expressed gene
    will be used.

    Input:
        taxon_id    int
    Output:
        df_genes    Dataframe with:
                    - log2foldchange
                    - p-value
                    - comparison_key
        df_go       Dataframe with:
                    - p adj (non-dir.)
                    - effect.size
        df_interpro Dataframe with:
                    - p adj (non-dir.)
                    - effect.size
        df_comparison_key Dataframe with:
                    - experiment
                    - analysis
                    - comparison
    """

    print((
        'Starting differential_gene_expression_of_taxon. '
        'Note that depending on taxon this might require ~30GB '
        'of RAM. Moreover the process might take ~20 minutes.'))

    def get_canonical_ensembl_names(taxon_id):
        beginners = {
            7227: 'FLYBASE:',
            6239: 'WormBase:',
            3702: 'Araport:',
            511145: 'EcoGene',
            559292: 'SGD'
        }
        if taxon_id in beginners.keys():
            beg = beginners[taxon_id]
        else:
            beg = 'Ensembl:'

        gene_info = meta.gene_info(taxon_id, usecols=['gene_ncbi', 'dbXrefs'])
        gene_info = utils.split_text_to_multiple_rows(
            gene_info, 'dbXrefs', '\|')
        f = gene_info['dbXrefs'].str.startswith(beg)
        gene_info = gene_info.loc[f, :].copy()
        gene_info[
            'dbXrefs'] = gene_info.loc[
                :, 'dbXrefs'].str.replace(r'^' + beg, '',)
        gene_info = gene_info.drop_duplicates()
        canonical_ensembl_names = set(gene_info['dbXrefs'].values)
        return canonical_ensembl_names

    canonical_names = get_canonical_ensembl_names(taxon_id)

    paths = glob.glob(os.path.join(
        inout.get_path('gtx_atlas'),
        'E*'))
    agg_genes = []
    current_comparison_key = 0
    agg_comparison_key = []
    agg_go = []
    agg_interpro = []

    for pe in paths:
        experiment = os.path.split(pe)[1]
        pa = glob.glob(os.path.join(
            pe, '{}[-_]*analytics.tsv'.format(experiment)))
        if len(pa) > 0:
            for paa in pa:   # sometimes there are several analysis

                df = pd.read_table(paa, low_memory=False)
                df = df[df['Gene ID'].isin(canonical_names)]

                if df.shape[0] > 0:

                    _, current_analysis = os.path.split(paa)

                    dh = pd.DataFrame(data=df.columns, columns=['name'])
                    comparisons = dh[dh['name'].str.endswith(
                        'log2foldchange')].loc[
                            :, 'name'].str.extract(
                                '^(.*)\.', expand=False).values

                    for co in comparisons:

                        agg_comparison_key.append({
                            'comparison_key': current_comparison_key,
                            'experiment': experiment,
                            'comparison': co,
                            'analysis': current_analysis
                        })

                        usecols = [
                            'Gene ID',
                            '{}.log2foldchange'.format(co),
                            '{}.p-value'.format(co),
                        ]
                        dff = df.loc[:, usecols]
                        dff = dff.rename(columns={
                            '{}.log2foldchange'.format(co): 'log2foldchange',
                            '{}.p-value'.format(co): 'p-value',
                        })
                        # sort to obtain the most significant value
                        dff['p-value'] = dff['p-value'].fillna(1)
                        dff = dff.sort_values('p-value', ascending=True)
                        f = dff.duplicated(subset='Gene ID', keep='first')
                        dff = dff.loc[~f, :]
                        dff.loc[:, 'comparison_key'] = current_comparison_key
                        agg_genes.append(dff)

                        p_go = os.path.join(
                            pe, '{}.{}.go.gsea.tsv'.format(
                                experiment,
                                co))

                        if os.path.exists(p_go):
                            try:
                                dff = pd.read_table(p_go, usecols=[
                                    'Term', 'p adj (non-dir.)', 'effect.size'],
                                    low_memory=False)
                                dff.loc[
                                    :, 'comparison_key'] = current_comparison_key
                                agg_go.append(dff)
                            except:
                                dummy = pd.read_table(p_go)
                                if dummy.shape[0] > 0:
                                    raise ValueError(
                                        'Unexpected format of {}'.format(
                                            p_go))

                        p_interpro = os.path.join(
                            pe, '{}.{}.interpro.gsea.tsv'.format(
                                experiment,
                                co))

                        if os.path.exists(p_interpro):
                            try:
                                dff = pd.read_table(p_interpro, usecols=[
                                    'Term', 'p adj (non-dir.)', 'effect.size'],
                                    low_memory=False)
                                dff.loc[
                                    :, 'comparison_key'] = current_comparison_key
                                agg_interpro.append(dff)
                            except:
                                dummy = pd.read_table(p_interpro)
                                if dummy.shape[0] > 0:
                                    raise ValueError(
                                        'Unexpected format of {}'.format(
                                            p_interpro))                                

                        current_comparison_key += 1

    # organize results of differential gene expression analysis
    df_genes = pd.concat(agg_genes, axis=0, ignore_index=True).rename(
        columns={'Gene ID': 'gene_ensembl'})
    df_genes = df_genes.loc[:, [
        'gene_ensembl',
        'log2foldchange',
        'p-value',
        'comparison_key']]
    del agg_genes

    require_canonical_id = True   # some samples carry names (diff. processed!)

    df_go = pd.concat(agg_go, axis=0, ignore_index=True).rename(
        columns={'Gene ID': 'gene_ensembl'})
    df_go = df_go.loc[:, [
        'Term',
        'p adj (non-dir.)',
        'effect.size',
        'comparison_key']].rename(columns={
            'Term': 'GO_ID'
        })
    if require_canonical_id:
        f = df_go['GO_ID'].str.contains('^GO:')
        df_go = df_go.loc[f, :]

    del agg_go

    df_interpro = pd.concat(agg_interpro, axis=0, ignore_index=True).rename(
        columns={'Gene ID': 'gene_ensembl'})
    df_interpro = df_interpro.loc[:, [
        'Term',
        'p adj (non-dir.)',
        'effect.size',
        'comparison_key']].rename(columns={
            'Term': 'interpro_id'
        })
    if require_canonical_id:
        f = df_interpro['interpro_id'].str.contains('^IPR[0-9]')
        df_interpro = df_interpro.loc[f, :]

    del agg_interpro

    df_comparison_key = pd.DataFrame(
        agg_comparison_key)
    df_comparison_key = df_comparison_key.loc[:, [
        'comparison_key',
        'experiment',
        'analysis',
        'comparison']]

    return df_genes, df_go, df_interpro, df_comparison_key
