import glob
import os

import numpy as np
import pandas as pd

from access_biology_data import gene_mapper, meta, relations
from access_science_shared import inout, mapper, utils

# Obtain properties of genes:
# - every function only uses taxon_id as input
# - all functions will return a table with one column being gene_ncbi


def allelepool_lek_2016_aberration(taxon_id=9606):
    """
    Variability in human populations (allele diversity): aberrations
    from anticipation;
    Original source data: Lek et al. 2015

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'lek2016_aberration_ordnum_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/lek_2016/{}'.format(n))

    df = pd.read_csv(p)

    def add_label(x):
        if not x.startswith('gene_ncbi'):
            x = 'Population variability Lek ' + x
        return x
    df.columns = [add_label(x) for x in df.columns]
    return df


def allelepool_lek_2016_anticipation(taxon_id=9606):
    """
    Anticipated variability in human populations (allele diversity)
    Original source data: Lek et al. 2015

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'lek2016_anticipation_ordnum_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/lek_2016/{}'.format(n))

    df = pd.read_csv(p)

    def add_label(x):
        if not x.startswith('gene_ncbi'):
            x = 'Population variability Lek ' + x
        return x
    df.columns = [add_label(x) for x in df.columns]
    return df


def aminoacids_swissprot(taxon_id):
    """
    Amino acids properties derived for proteins in the the swissprot
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          dataframe
    """

    p = inout.get_path(
        'geisen',
        'aminoacids/aminoacids_swissprot_{}_gene_id.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def aminoacids_trembl(taxon_id):
    """
    Amino acids properties derived for proteins in the the trembl
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          dataframe
    """

    p = inout.get_path(
        'geisen',
        'aminoacids/aminoacids_trembl_{}_gene_id.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def chromosomes(taxon_id):
    """
    Chromosomes as provided by NCBI. Note that this can also list
    non-canonical chromosme names that can convey furtber information,
    such as - or X|Y or Un or '10|19|3'

    Input:
        taxon_id    int

    Output:
        df          dataframe with categorial data
    """

    def _load_chromosomes_geisen_v1(p, taxon_id):
        print("""
            Found chromosomes.h5 . This is indicative of usage of
            an outdated version of the geisen datasets.

            Using legacy fallback. Note that this option might
            become deactivated in a future release.
            """)
        q = 'taxon_ncbi=={}'.format(taxon_id)
        df = pd.read_hdf(p, 'table', where=q)
        col_to_use = ['gene_ncbi', 'chromosome']
        df = df.loc[:, col_to_use]
        return df

    def _load_chromosomes_geisen_v1_1(taxon_id):
        p = inout.get_path(
            'geisen',
            'ncbi/chromosomes/chromosomes_taxon_{}.csv.gz'.format(
                int(taxon_id)))
        df = pd.read_csv(p, usecols=['gene_ncbi', 'chromosome'])
        return df

    p = inout.get_path(
        'geisen',
        'ncbi/chromosomes.h5')

    if os.path.exists(p):
        df = _load_chromosomes_geisen_v1(p, taxon_id)
    else:
        df = _load_chromosomes_geisen_v1_1(taxon_id)

    return df


def compartment_itzhak_2016_cytoplasmic(taxon_id=9606):
    """
    Cytoplasmic location (note complements with non-cytoplasmic
    localization to 100); Note that the amount of recoreds is higher
    than the amount of records for other predictions / features
    of the original underlying publication (the finely resolved
    localization)

    Input:
        taxon_id    int safety

    Ouput:
        df          nominal data
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    def _load_cytoplasmic_v1():
        df = compartment_itzhak_2016_global_scores(taxon_id)
        col_to_use = [
            'gene_ncbi',
            'Localization Itzhak Cytosolic Pool']
        df = df.loc[:, col_to_use]
        return df

    p = inout.get_path(
        'geisen',
        'papers/itzhak_2016/' +
        'itzhak2016_localization_cytoplasm_ncbi_gene.csv.gz')

    if os.path.exists(p):    # geisen v1_1 or higher
        df = pd.read_csv(p)
    else:
        df = _load_cytoplasmic_v1(taxon_id)

    return df


def compartment_itzhak_2016_determined(taxon_id=9606):
    """
    Compartment data from Itzhak et al.

    In contrast to the compartment_itzhak_2016_global_classifier function,
    this function only contains the ~4,500 gnes where location could be
    determined by the actual study

    Input:
        taxon_id    (safety check for usage on humans, as original data)

    Output:
        df          databframe, boolean columns individual compartments

    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'papers/itzhak_2016/itzhak2016_compartment_nombool_ncbi_gene.csv.gz')
    df = pd.read_csv(p)
    return df


def compartment_itzhak_2016_global_classifier(taxon_id=9660):
    """
    Compartment data from Itzhak et al.

    recommendation: use compartment_itzhak_2016_determined instead

    In contrast to the compartment_itzhak_2016_determined function,
    this function also contains around 3000 recoreds where location
    could not be determined, and two columns indicating abasence of
    determination of subcellular localization.

    Note that with a single exception there is no compartment allocated
    to any record which also carries one one of the tags for nondetermiantion

    Input:
        taxon_id    (safety check for usage on humans, as original data)

    Output:
        df          boolean columns individual compartments – including
                    two columns about determination status

    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'papers/itzhak_2016/' +
        'itzhak2016_global_classifier_nombool_ncbi_gene.csv.gz')
    df = pd.read_csv(p)

    def add_global(x):
        if x.startswith('Compartment'):
            x = 'Global ' + x
        return x

    df.columns = [add_global(x) for x in df.columns]
    return df


def compartment_itzhak_2016_global_scores(taxon_id=9606):
    """
    Prediction accuracy, and broad classification by nucleo-cyplasmic
    localization

    Input:
        taxon_id    int (safety check)

    Output:
        df      ordinal numbers

    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'papers/itzhak_2016/' +
        'itzhak2016_localization_stats_ordnum_ncbi_gene.csv.gz')
    df = pd.read_csv(p)

    def add_global(x):
        if not x.startswith('gene_ncbi'):
            x = 'Localization Itzhak ' + x
        return x

    df.columns = [add_global(x) for x in df.columns]
    return df


def compartment_thul_2017_main(taxon_id=9606):
    """
    Main localization as reported by Thul et al.
    (human protein atlas)

    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    df = _fetch_thul_2017(
        scope='main',
        confidence=[   # don't consider uncertain
            'approved',
            'supported',
            'validated',
        ])
    df = df.reset_index()

    def add_global(x):
        if not x.startswith('gene_ncbi'):
            x = 'Location Main Thul  ' + x
        return x

    return df


def compartment_thul_2017_observed(taxon_id=9606):
    """
    Any localization observed by Thul et al.
    (human protein atlas)

    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    df = _fetch_thul_2017(
        scope='observed',
        confidence=[   # don't consider uncertain
            'approved',
            'supported',
            'validated',
        ])
    df = df.reset_index()

    def add_global(x):
        if not x.startswith('gene_ncbi'):
            x = 'Location Observed Thul  ' + x
        return x

    return df


def essentiality_blomen_2015(taxon_id=9606):
    """
    Human gene essentiality as obtained from insertion into genes.
    Original source data: Blomen et al. 2015

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'papers/blomen_2015/blomen_2015_fitness_ncbi_gene.csv.gz')
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def essentiality_blomen_2015_sh(taxon_id=9606):
    """
    Human gene essentiality as obtained shRNA screen
    Original source data: Hart et al. 2015

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'hart_2015_hart2015_hct116_shRNA_ordnum_ordnum_gene_ncbi.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/hart_2015/{}'.format(n))

    df = pd.read_csv(p)
    df = df.dropna()  # to account for unknown design of sh library
    return df


def essentiality_hart_2015_crispr(taxon_id=9606):
    """
    Human gene essentiality as obtained from insertion into genes.
    Original source data: Hart et al. 2015

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    def load_file(code):
        n = 'hart_2015_hart2015_{}_ordnum_ordnum_gene_ncbi.csv.gz'.format(code)
        p = inout.get_path(
            'geisen',
            'papers/hart_2015/{}'.format(n))

        if os.path.exists(p):
            df = pd.read_csv(p)
        else:
            print(p)
            raise EnvironmentError('Did not find data for {}'.format(code))
        return df

    ds = ['a375_ko', 'dld1', 'gbm', 'hct116', 'hela', 'rpe1']

    dn = load_file(ds[0])

    for d in ds[1:]:
        dnn = load_file(d)
        dn = pd.merge(
            dn,
            dnn,
            left_on='gene_ncbi',
            right_on='gene_ncbi',
            how='outer')

    return dn


def essentiality_wang_2015_pvalue(taxon_id=9606):
    """
    Human gene essentiality as obtained from insertion into genes.
    Original source data: Wang et al. 2015

    p-values as given by authors

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'wang_2015_pvalue_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/wang_2015/{}'.format(n))

    df = pd.read_csv(p)
    return df


def essentiality_wang_2015_score(taxon_id=9606):
    """
    Human gene essentiality as obtained from insertion into genes.
    Original source data: Wang et al. 2015

    Confidence score as given by authors

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'wang_2015_cs_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/wang_2015/{}'.format(n))

    df = pd.read_csv(p)
    return df


def gdi(taxon_id=9606):
    """
    GDI mutagenesis intolerance
    source: GDI through genealacart

    Input:
        taxon_id    int (safety)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError(
            'Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'genealacart/genealacart_intolerance_gdi.gz')

    df = pd.read_csv(p)
    return df


def genbank_gene(taxon_id):
    """
    Features from genes; Source is one of the most complete genbank
    releases (manually selected one per species). Features have been
    extracted for a manually selected set of species (roughly
    corresponding to heavily studied species.)

    Input:
        taxon_id    int
    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'genbank/genbank_gene/genbank_gene_{}.csv'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def genbank_genomic_cds(taxon_id):
    """
    Features for predicted coding sequences; Source is one of the
    most complete genbank releases (manually selected one per species).
    Features have been extracted for a manually selected set of species
    (roughly corresponding to heavily studied species.)

    For genbank_genomic_cds the features correspond to individual
    nucleotides and individual codons.

    Note that genomic CDS may not be defined for some species within
    the original data source (genbank)

    Input:
        taxon_id    int
    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'genbank/genomic_cds/genbank_genomic_cds_{}.csv'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        print(p)
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def genbank_genomic_rna(taxon_id):
    """
    Features for genomic RNA; Source is one of the
    most complete genbank releases (manually selected one per species).
    Features have been extracted for a manually selected set of species
    (roughly corresponding to heavily studied species.)

    For genbank_genomic_rna the features correspond to
    individual nucleotides.

    Note that genomic RNA may not be defined for some species within
    the original data source (genbank)

    Input:
        taxon_id    int
    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'genbank/genomic_rna/genbank_genomic_rna_{}.csv'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def genbank_validated_rna(taxon_id):
    """
    Features for RNA with high support by genbank curators;
    Source is one of the most complete genbank releases
    (manually selected one per species). Features have been extracted
    for a manually selected set of species
    (roughly corresponding to heavily studied species.)

    For genbank_validated_rna the features correspond to
    individual nucleotides, codons, and cdon bias.

    Note that validated RNA may not be defined for some species within
    the original data source (genbank)

    Input:
        taxon_id    int
    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'genbank/validated_rna/genbank_validated_rna_{}.csv'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


# def homologenes_for_all_protein_coding(taxon_id):
#     """
#     Presence of homologs in homologene; Note that this
#     function will expand list of genes to all protein
#     coding genes in taxon_id

#     Input:
#         taxon_id    int
#     Output
#         df          boolean
#     """

#     df = _fetch_homologs(taxon_id)
#     g = standardizer.reference_genes(taxon_id, 'p')
#     df = df.loc[g, :]
#     df = df.fillna(False)
#     df = df.reset_index()

#     return df


def homologenes_only_for_genes_in_homologene(taxon_id):
    """
    Presence of homologs in homologene; Note that this
    function will expand list of genes to all protein
    coding genes in taxon_id

    Input:
        taxon_id    int
    Output
        df          boolean
    """

    df = _fetch_homologs(taxon_id)
    df = df.fillna(False)
    df = df.reset_index()

    return df


def interactions_rolland_2014(taxon_id=9606):
    """
    Presence of interaction with self, or other proteins
    Original source data: Rolland et al. 2014

    Stratified by any support and multiple support

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers, and boolean
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'rolland_counts_of_interactions.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/rolland_2014/{}'.format(n))

    df = pd.read_csv(p)
    df = df.rename(columns={'Unnamed: 0': 'gene_ncbi'})

    return df


def protein_abundance_itzhak_2015(taxon_id=9606):
    """
    Protein abundance in a human cell line
    Original source data: Itzhak et al. 2015

    Input:
        taxon_id    int (safety check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'papers/itzhak_2016/' +
        'itzhak2016_protein_abundance_ordnum_ncbi_gene.csv.gz')
    df = pd.read_csv(p)

    def add_label(x):
        if not x.startswith('gene_ncbi'):
            x = 'Protein Itzhak ' + x
        return x

    df.columns = [add_label(x) for x in df.columns]
    return df


def radar_swissprot(taxon_id):
    """
    RADAR (repetitiveness) properties of the most repetitive motive
    derived for proteins in the the swissprot
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'radar/radar_swissprot_{}_gene_id.csv.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def radar_trembl(taxon_id):
    """
    RADAR (repetitiveness) properties of the most repetitive motive
    derived for proteins in the TrEMBL
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'radar/radar_trembl_{}_gene_id.csv.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def rvis(taxon_id=9606):
    """
    RVIS mutagenesis intolerance
    source: RVIS through genealacart

    Input:
        taxon_id    int (safety)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError(
            'Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'geisen',
        'genealacart/genealacart_intolerance_rvis.gz')

    df = pd.read_csv(p)
    return df


def seg_swissprot(taxon_id):
    """
    SEG (unstructuredness) properties derived for proteins in the swissprot
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'seg/seg_swissprot_{}_gene_id.csv.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def seg_trembl(taxon_id):
    """
    SEG (unstructuredness) properties derived for proteins in the TrEBML
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'seg/seg_trembl_{}_gene_id.csv.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def signalp_swissprot(taxon_id):
    """
    SignalP (signal peptide) properties derived for proteins in the swissprot
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'signalp/signalp_swissprot_{}_gene_id.csv.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def signalp_trembl(taxon_id):
    """
    SignalP (signal peptide) properties derived for proteins in the swissprot
    database. These features have been computed locally for a selected
    set of taxa (basically: the most heavily studied ones)

    Input:
        taxon_id    int

    Output:
        df          ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'signalp/signalp_trembl_{}_gene_id.csv.gz'.format(
            taxon_id))
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        raise EnvironmentError(
            'Did not find data for taxon {}'.format(taxon_id))
    return df


def transcript_abundance_flybase_cells(taxon_id=7227):
    """
    Abundance of transcripts in cell lines
    Source: Flybase

    NaN if original value below 1 (0 on log scale) or below detection threshold

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers, and NaN
                    (NaN if below reliable measurement threshold)
    """

    if taxon_id != 7227:
        raise EnvironmentError(
            'Only supports taxon 7227, Drosophila melanogaster')
    n = 'flybase_expression_cells_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'flybase/{}'.format(n))

    df = pd.read_csv(p)
    return df


def thermal_stability_leuenberger_2017(taxon_id):
    """
    Thermal stability of proteins; supports e.coli,
    T. thermophilus, S. cerevisiae, Human HeLa cells

    Input:
        taxon_id    int
    Output:
        df          ordinal

    """

    sheet_names = {
        9606: 'Human HeLa Cells',
        511145: 'E. coli',  # haven't comparted protein IDs to strain
        300852: 'T. thermophilus',  # haven't comparted protein IDs to strain
        559292: 'S. cerevisiae'  # haven't comparted protein IDs to strain
    }

    p = inout.get_path(
        'publications',
        'leuenberger2017/aai7825_Leuenberger_Table-S3.xlsx')
    df = pd.read_excel(p, sheet_name=sheet_names[taxon_id])

    df = df[['Protein_ID', 'Tm Peptide', 'Tm Protein', 'T 90% Unfolded']]
    df = df.rename(columns={
        'Protein_ID': 'protein_uniprot',
        'Tm Peptide': 'tm_peptide',
        'Tm Protein': 'tm_protein',
        'T 90% Unfolded': 'tm_90p'
    })

    df = utils.split_text_to_multiple_rows(
        df, 'protein_uniprot', ';')
    df = mapper.uniprot_protein_2_gene_ncbi(
        df.dropna(), 'median').reset_index()

    return df


def transcript_abundance_flybase_stages(taxon_id=7227):
    """
    Abundance of transcripts within developmental stages
    Source: Flybase

    NaN if original value below 1 (0 on log scale) or below detection threshold

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers
                    and NaN (below reliable measurement threshold)
    """

    if taxon_id != 7227:
        raise EnvironmentError(
            'Only supports taxon 7227, Drosophila melanogaster')
    n = 'flybase_expression_stages_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'flybase/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcript_abundance_flybase_tissues(taxon_id=7227):
    """
    Abundance of transcripts within tissues
    Source: Flybase

    NaN if original value below 1 (0 on log scale) or below detection threshold

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers
                    and NaN (below reliable measurement threshold)
    """

    if taxon_id != 7227:
        raise EnvironmentError(
            'Only supports taxon 7227, Drosophila melanogaster')
    n = 'flybase_expression_tissues_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'flybase/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcript_abundance_flybase_treatments(taxon_id=7227):
    """
    Abundance of transcripts under treatments
    Source: Flybase

    NaN if original value below 1 (0 on log scale) or below detection threshold

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers
                    and NaN (below reliable measurement threshold)
    """

    if taxon_id != 7227:
        raise EnvironmentError(
            'Only supports taxon 7227, Drosophila melanogaster')
    n = 'flybase_expression_treatments_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'flybase/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcript_abundance_gerstein(taxon_id=6239):
    """
    Abundance of transcripts under mixed conditoins
    Source: modENCODE through gesrtein lab

    NaN if original value below 1 (0 on log scale) or below detection threshold

    Note that their meta-annotation is horrible, and misleading, and does
    not allow to separate treatments from tissues etc.
    While modEncode was helpful, and relayed request to Gerstein lab,
    the Gerstein lab never replied (neither to modEncode help, nor me),
    although this highl-level dataset is shown as part of modEncode


    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers,
                    and NaN (below reliable measurement threshold)
    """

    if taxon_id != 6239:
        raise EnvironmentError('Only supports taxon 6239, C. elegans')
    p = inout.get_path(
        'geisen',
        'gerstein/gerstein_expression__ncbi_gene.csv.gz')

    df = pd.read_csv(p)
    return df


def transcript_abundance_gex_mantalek_170222(mask):
    """
    Tissue specific gene expression
    Source: EBML-EBI Expression Atlas (https://www.ebi.ac.uk/gxa/)
            Selected datasets manually downloaded by M. Antalek on
            170222 using cutoff of 0

    Unlike many other gene expression datasets, the returned gene
    expression values are not log-transformed


    Input:
        mask    int or str; e.g: taxon_id, or taxon_id-condition
                (where condition is sample specific, e.g.: 10116-female)

    Output:
        df      ordinal numbers
    """

    p = inout.get_path(
        'geisen',
        'gxa/matt_antalek_170222/*-{}[-_]*_gene.csv.gz'.format(
            mask))

    g = glob.glob(p)

    if len(g) == 0:
        raise EnvironmentError(
            'Did not find any dataset matching the mask {}'.format(
                mask))

    df = pd.read_csv(g[0])

    if len(g) > 1:
        for gn in g:
            df_n = pd.read_csv(gn)

            ref_c = 'gene_ncbi'
            if sorted(df[ref_c].values) == sorted(df_n[ref_c].values):
                df = pd.merge(
                    df, df_n, left_on=ref_c, right_on=ref_c, how='outer')

        if any(df.isnull().sum()):
            raise EnvironmentError(
                'Some measurements are lacking in one of the datasets.')

    return df


def transcript_abundance_uhlen_2015_cells(taxon_id=9606):
    """
    Transcript abundance for ~45 cell lines.
    Source: Uhlen et al. 2015

    NaN if original value below 1 (0 on log scale) or below detection threshold

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers,
                    and NaN (below reliable measurement threshold)
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'uhlen_2015_cells_levels_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/uhlen_2015/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcript_abundance_uhlen_2015_tissues(taxon_id=9606):
    """
    Transcript abundance for ~120 tissue samples (include
    biological replicates).
    Source: Uhlen et al. 2015

    NaN if original value below 1 (0 on log scale) or below detection threshold

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers,
                    and NaN (below reliable measurement threshold)
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'uhlen_2015_tissue_levels_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/uhlen_2015/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcript_detection_uhlen_2015_cells(taxon_id=9606):
    """
    Fraction of cell lines with expression >= 1FPKM
    Source: Uhlen et al. 2015

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'uhlen_2015_detected_in_cells_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/uhlen_2015/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcript_detection_uhlen_2015_tissues(taxon_id=9606):
    """
    Fraction of tissue samples with expression >= 1FPKM
    Source: Uhlen et al. 2015

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'uhlen_2015_detected_in_tissuess_ncbi_gene.csv.gz'
    p = inout.get_path(
        'geisen',
        'papers/uhlen_2015/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcription_factors_genealacart_encode(taxon_id=9606):
    """
    Occurences of transcription factors in enhancers derived from Encode
    Source: Genealacart

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'genealacart_encode_tfs_by_gene.gz'
    p = inout.get_path(
        'geisen',
        'genealacart/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcription_factors_genealacart_promoters(taxon_id=9606):
    """
    Occurences of transcription factors in promoter regions
    Source: Genealacart

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordinal numbers
    """

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')
    n = 'genealacart_promoters_tfs_by_gene.gz'
    p = inout.get_path(
        'geisen',
        'genealacart/{}'.format(n))

    df = pd.read_csv(p)
    return df


def transcript_halflife_tani_2012_assuming_48h_for_stable(taxon_id=9606):
    """
    Transcript halflife in HeLa cells, as measured by Tani et al.
    Note that the orgiginal data is on RNA, and different RNA of
    same gene can get pooled (here: median). Thus this function uses
    an educated guess for > 24h; more specifically, 48h hours will be used
    so that, if one transcript is >24 and there is only one <24h,
    the presence of a long lived transcript will be shown in the
    aggregated readout of both RNA species of the same gene

    Input:
        taxon_id    int (saftey check)

    Output:
        df          ordnum (having given categorial value >24 
                        the auxiliary value 48) - see above
    """

    hours_to_assume_for_over_24 = 48   # manually selected by Thomas Stoeger

    if taxon_id != 9606:
        raise EnvironmentError('Only supports taxon 9606, Homo sapiens')

    p = inout.get_path(
        'publications',
        'tani2012/Tani_Supp_Tables_revised2.xls')

    df = pd.read_excel(p, sheet_name='Table S1', skiprows=3)
    f = df['t1/2 (h)'] != 'N.D.'
    df = df.loc[f, ['RepName', 't1/2 (h)']].rename(columns={
        'RepName': 'rna_ncbi',
        't1/2 (h)': 'rna_halflife_h'
    })

    df['rna_ncbi'] = df['rna_ncbi'].str.strip(',')
    df = utils.split_text_to_multiple_rows(df, 'rna_ncbi', ',')

    df['rna_halflife_h'] = df['rna_halflife_h'].replace(
        '>24', hours_to_assume_for_over_24).astype(float)

    df = df.groupby('rna_ncbi').agg(np.median)
    df = mapper.rna_ncbi_2_gene_ncbi(df, 'median')
    df = df.reset_index()

    return df


# GENE SET SPECIFIC FUNCTIONS (supporting more arguments)


def _fetch_homologs(taxon_id):
    """
    Fetches homologs, and creates occurence matrix where columns
    are different organisms

    Input:
        taxon_id        int
    Output:
        df              bool
    """

    hg = relations.homologene()
    h = hg[hg['taxon_ncbi'] == taxon_id]['homologene_group'].values
    hg = hg[hg['homologene_group'].isin(h)]
    cis_group = hg[hg['taxon_ncbi'] == taxon_id][[
        'gene_ncbi', 'homologene_group']]
    hg = hg[['homologene_group', 'taxon_ncbi']]
    df = pd.merge(hg, cis_group, how='outer')

    d = {x: meta.taxon_name(x).replace(' ', '_').replace(
        '-', '_').lower() for x in df['taxon_ncbi'].unique()}
    d = pd.Series(d).to_frame('conservation').reset_index().rename(
        columns={'index': 'taxon_ncbi'})
    df = pd.merge(df, d)

    df = df[df['taxon_ncbi'] != taxon_id][[
        'gene_ncbi', 'conservation']].drop_duplicates()

    df['annotation_id'] = df['conservation']
    df = df.rename(columns={'conservation': 'annotation_name'})

    hg = relations.homologene()
    h = hg[hg['taxon_ncbi'] == taxon_id]['homologene_group'].values
    hg = hg[hg['homologene_group'].isin(h)]
    cis_group = hg[hg['taxon_ncbi'] == taxon_id][[
        'gene_ncbi', 'homologene_group']]
    hg = hg[['homologene_group', 'taxon_ncbi']]
    df = pd.merge(hg, cis_group, how='outer')

    d = {x: meta.taxon_name(x).replace(' ', '_').replace(
        '-', '_').lower() for x in df['taxon_ncbi'].unique()}
    d = pd.Series(d).to_frame('conservation').reset_index().rename(
        columns={'index': 'taxon_ncbi'})
    df = pd.merge(df, d)

    df = df[df['taxon_ncbi'] != taxon_id][[
        'gene_ncbi', 'conservation']].drop_duplicates()
    df.loc[:, 'is_present'] = True
    df = df.pivot(index='gene_ncbi', columns='conservation',
                  values='is_present')

    df = df.fillna(False)

    df.columns = ['homologene_{}'.format(x) for x in df.columns]

    return df


def _fetch_thul_2017(scope, confidence):
    """
    Loads subcellular localization as obtained by Thul et al. 2017
    (human protein atlas). Also adds their prediction on secretion.
    Will substitute ensemble gene IDs of original data with
    ncbi gene IDs (Enterz IDs)

    Input:

    scope           scope of localization; either:
                        observed    any localization
                        main        main localization(s)
    confidence      any, or a list containing elements of
                        approved
                        supported
                        validated
                        uncertain

    Output:
        df          data frame with entrez gene in index,
                    and boolean occurence matrix for
                    localizations

    """

    p_in_dir = inout.get_path(
        'publications',
        'thul2017')

    df_secreted = pd.read_excel(
        os.path.join(
            p_in_dir,
            'aal3321_Thul_SM_table_S9.xlsx'
        ),
        sheet_name='HPA prediction')[
            ['Ensembl']].rename(columns={'Ensembl': 'gene_ensembl'})
    df_secreted.loc[:, 'df_secreted'] = True

    df_cellular_localization_full_table = pd.read_excel(
        os.path.join(
            p_in_dir,
            'aal3321_Thul_SM_table_S6.xlsx'
        ),
        sheet_name='Protein location results').rename(
        columns={'ENSG': 'gene_ensembl'})

    if confidence is not 'any':
        df_cellular_localization_full_table['Reliability'] = \
            df_cellular_localization_full_table['Reliability'].str.lower()

        f = df_cellular_localization_full_table['Reliability'].isin(confidence)
        df_cellular_localization_full_table = \
            df_cellular_localization_full_table.loc[f, :]

    c = sorted([
        'gene_ensembl', 'Nucleus', 'Nucleoplasm',
        'Nuclear bodies', 'Nuclear speckles', 'Nuclear membrane', 'Nucleoli',
        'Nucleoli (Fibrillar center)', 'Cytosol', 'Cytoplasmic bodies',
        'Rods and Rings', 'Lipid droplets', 'Aggresome', 'Mitochondria',
        'Microtubules', 'Microtubule ends', 'Microtubule organizing center',
        'Centrosome', 'Mitotic spindle', 'Cytokinetic bridge', 'Midbody',
        'Midbody ring', 'Intermediate filaments', 'Actin filaments',
        'Focal Adhesions', 'Endoplasmic reticulum', 'Golgi apparatus',
        'Vesicles', 'Plasma membrane', 'Cell Junctions'])

    df_cellular_localization = df_cellular_localization_full_table.loc[:, c]

    df_cellular_main = df_cellular_localization_full_table.loc[
        :,
        ['gene_ensembl', 'IF main protein location']
    ]

    df_cellular_main = utils.split_text_to_multiple_rows(
        df_cellular_main,
        'IF main protein location',
        r';')

    df_cellular_main.loc[:, 'present'] = True

    df_cellular_main = df_cellular_main.pivot(
        index='gene_ensembl',
        columns='IF main protein location',
        values='present'
    )

    df_cellular_main = df_cellular_main.reset_index().loc[:, c].fillna(False)

    def to_entrez(df):
        df = gene_mapper.gene_ensembl_2_gene_ncbi_unambiguously(df, 9606)
        return df

    df_secreted = to_entrez(df_secreted)
    df_cellular_localization = to_entrez(df_cellular_localization)
    df_cellular_main = to_entrez(df_cellular_main)
    df_secreted = df_secreted.rename(
        columns={'df_secreted': 'Secreted (prediction)'})

    if scope == 'observed':
        df = df_cellular_localization
    elif scope == 'main':
        df = df_cellular_main
    else:
        raise EnvironmentError(
            'Scope must be observed or main.')

    if df.index.value_counts().max() > 1:
        raise EnvironmentError(
            'Something wrong')

    add_secretion = True
    if add_secretion:
        df = pd.merge(
            df,
            df_secreted,
            left_index=True,
            right_index=True,
            how='left')
        df[
            'Secreted (prediction)'] = df[
            'Secreted (prediction)'].fillna(False)

    df = df.sort_index()

    return df
