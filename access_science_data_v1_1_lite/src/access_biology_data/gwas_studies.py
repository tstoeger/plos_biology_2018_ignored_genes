import pandas as pd


# from access_biology_data import meta
from access_science_shared import inout


def ebi_gwas():
    """
    Will load EBI gwas table

    Output:
        df     GWAS table
    """

    # import GTEX dataset
    p = inout.get_path(
        'ebi_gwas', 'full.tsv')

    df = pd.read_table(p, low_memory=False)

    df = df.rename(columns={
        'PUBMEDID': 'pubmed_id'
    })

    df['year'] = df['DATE'].copy().str.extract(
        '^([0-9]{4})-*', expand=False)

    df = df.set_index(['MAPPED_GENE', 'DISEASE/TRAIT']).reset_index()

    return df
