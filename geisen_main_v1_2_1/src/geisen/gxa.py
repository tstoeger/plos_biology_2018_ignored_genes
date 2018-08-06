import glob
import os
import re

import pandas as pd

import geisen.inout as io
from geisen import mapper
from geisen.prepare import _save_orig_and_ncbi_gene_mapped_tables


def matt_antalek_170222():
    """
    Matt Antalek (Rick Morimoto lab)
    downloaded on 170222 tissue data of several
    model organisms; Used cutoff was 0, and when a filter would
    be required by the web-interface he chose reasonable
    representative ones
    """

    # manually curated condition codes:
    # dictionary with extension as key, and entries
    # - taxon_id
    # - if qualifier: [taxon_id, qualifier]
    condition_codes = {
        'rattus_norvegicus_female': [10116, 'female'],
        'rattus_norvegicus_male': [10116, 'male'],
        'ovis_aries_texel': [9940, 'texel'],
        'ovis_aries_female': [9940, 'female'],
        'ovis_aries_male': [9940, 'male'],
        'mus_musculus': 10090,
        'bos_taurus': 9913,
        'gallus_gallus': 9031,
        'macaca_mulatta': 9544,
        'homo_sapiens': 9606,
        'pabio_anubis': 9555,  # olive baboon
        'monodelphis_domestica': 13616,
        'xenopus_tropicalis': 8364,
        'anolis_carolinesis': 28377,
    }

    p_dir_in = io.get_geisen_manual_data_path(
        'out/'
        'ebi_expression_manual/'
        'matt_antalek_170222/'
        'E-*.tsv')  # filter for correct files

    p_out = io.get_output_path('gxa/matt_antalek_170222')
    io.ensure_presence_of_directory(p_out)

    files = glob.glob(p_dir_in)

    for p in files:

        df = pd.read_table(p, header=3)
        df = df.rename(columns={'Gene ID': 'gene_ensembl'})
        df = df.drop('Gene Name', axis=1)

        def add_GXA_to_label(x):    # introduced in geisen v1_1
            if not x.startswith('gene'):
                x = 'GXA_' + x
            return x
        df.columns = [add_GXA_to_label(y) for y in df.columns]

        _, fname = os.path.split(p)

        matched = re.findall('^(.*)-[0-9].*-results_(.*)\.tsv', fname)

        if len(matched) != 1:
            raise ValueError('Unexpected format. Check parsing pattern.')

        experiment = matched[0][0]

        k = matched[0][1]
        meta = condition_codes[k]

        if isinstance(meta, list):
            taxon_id = meta[0]
            condition = meta[1]
            v = '{}-taxon_id-{}-{}'.format(experiment, taxon_id, condition)
        elif isinstance(meta, int):
            taxon_id = meta
            v = '{}-taxon_id-{}'.format(experiment, taxon_id)
        else:
            raise ValueError('Unexpected format. Check condition_codes.')

        taxa_without_nih_ensembl = [
            8364]

        if taxon_id not in taxa_without_nih_ensembl:

            # If NIH has corresponding ensembl for ncbi gene IDs,
            # save original, and ncbi_gene mapped

            df_entrez = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
                df, taxon_id)

            _save_orig_and_ncbi_gene_mapped_tables(
                p_dir=p_out,
                filebase=v,
                df_orig=df,
                df_ncbi=df_entrez)

        else:  # for some taxa NIH does not have mapping to ensembl

            df.to_csv(
                os.path.join(p_out, '{}_orig.csv.gz'.format(v)),
                compression='gzip',
                index=True)
