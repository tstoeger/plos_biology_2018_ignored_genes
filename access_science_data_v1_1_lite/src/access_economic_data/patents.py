import os
import pandas as pd

from access_science_shared import inout


def rosenfeld_2013(cols_to_use=None):
    """
    Will load gene-linked patent data as described by
    Rosenfeld et al. ; Note that most patents are filed
    for sequences, rather than genes, and are ambiguous,
    and including related genes, on purpose. See
    publication of Rosenfeld et al. on that
    socio-economical problem, and for specific
    cutoff / mapping scheme applied by them (Note:
    which sounds like a very reasonable compromise)

    Note: original datset might have some excess
    records that could not be unambiguously mapped
    to gene_ncbi (Entrez gene IDs)

    Output:
    linkage table between patent and gene_ncbi

    """

    p_folder = inout.get_path(
        'geisen',
        os.path.join(
            'papers',
            'rosenfeld_2013',
            'rosenfeld_2013_patents_ncbi_gene.csv.gz'))

    df = pd.read_csv(
        p_folder,
        usecols=['patent', 'gene_ncbi']) # skip index column
    df = df.drop_duplicates()

    return df
