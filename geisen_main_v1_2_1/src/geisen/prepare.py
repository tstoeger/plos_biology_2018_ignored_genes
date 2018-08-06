import glob
import os
import re
import sys

import numpy as np
import pandas as pd

from Bio import SeqIO

from geisen import mapper
from geisen import settings
from geisen import downloader

import geisen.inout as io

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Functions that can not be placed in shared repository;
# This includes unpublished work of others, or code
# whose usage is restrained because of legal permissions
s_dir = './../../../geisen_illegal_if_public/code/src/'
sys.path[0] = s_dir

# from Sophia Liu, Amaral Lab
import ENCCalculator as CodonBiasCalculator


'''
############### TIER 1 ###########################
biogrid
biosystems
genbank_gene_features
genbank_validated_rna_features
genbank_cds_rna_features
genbank_genomic_rna_features
gene_info
gene2ensembl
gene2go
gene2pubmed
generifs_basic
genetic_testing_registry
uniprot_databases
uniprot_id_mapper
taxdmp
'''


def biogrid():
    """
    Will download the biogrid protein interaction database
    - Separates heterologous interactions from within-species interaction
    - Removes non gene_ncbi gene names
    """

    # Initialize Absent Output
    out_p = 'biogrid'
    out_p = io.get_output_path(out_p)
    io.ensure_presence_of_directory(out_p)

    stored_locations = downloader.download_data_set(
        'biogrid', folder_contains_dots=True)

    if len(stored_locations) != 1:
        raise EnvironmentError(
            'biogrid download is ambiguous')
    else:
        p = stored_locations[0]

    c = [  # Everything except synonys and alternative names (symols)
        '#BioGRID Interaction ID', 'Entrez Gene Interactor A',
        'Entrez Gene Interactor B', 'BioGRID ID Interactor A',
        'BioGRID ID Interactor B', 'Experimental System',
        'Experimental System Type', 'Author', 'Pubmed ID',
        'Organism Interactor A', 'Organism Interactor B', 'Throughput',
        'Score',
        'Modification', 'Phenotypes', 'Qualifications', 'Tags',
        'Source Database']

    df = pd.read_table(p, low_memory=False, usecols=c)
    df = df.rename(
        columns={'#BioGRID Interaction ID': 'BioGRID Interaction ID'})

    is_heterologous = df.loc[
        :, 'Organism Interactor A'] != df.loc[:, 'Organism Interactor B']
    df_heterologous = df.loc[is_heterologous, :]
    o = os.path.join(out_p, 'biogrid_heterologous.csv.gz')
    df_heterologous.to_csv(o, compression='gzip', index=False)

    df_homologous = df.loc[~is_heterologous, :]
    df_homologous = df_homologous.drop('Organism Interactor B', axis=1)

    grouped = df_homologous.groupby('Organism Interactor A')

    for taxon, d in grouped:
        d = d.reset_index(drop=True)
        o = os.path.join(
            out_p, 'biogrid_homologous_taxon_{}.csv.gz'.format(taxon))
        d.to_csv(o, index=False, compression='gzip')


def biosystems():
    """
    Will download biosystems from ncbi, and format it for further usage.
    Note that biosystems exists in two versions on NCBI (with and without GO).
    Here the version without GO is used to avoid duplications with gene2go
    - adhere to science of science naming conventions
    - save individual data sources (e.g.: KEGG) as flat file
    - merge gene2biosystem table with biosystem metainfo
    - add column names from documentation (manual)
    - remove entries without annotated source (around 1%-2%)
    """

    # Initialize Absent Output
    out_p = 'ncbi/biosystems'
    out_p = io.get_output_path(out_p)
    io.ensure_presence_of_directory(out_p)

    number_of_files = io.check_number_of_files_in_directory(out_p, '*.gz')
    if number_of_files > 0:
        raise EnvironmentError(
            'At least one output file is already present.')

    # Download Dataset
    stored_locations = downloader.download_data_set('biosystems')

    if len(stored_locations) == 2:

        if 'biosystems_gene' in stored_locations[0]:
            p_gene = stored_locations[0]
        elif 'biosystems_gene' in stored_locations[1]:
            p_gene = stored_locations[1]
        else:
            raise ValueError('Could not biosystems_gene in download.')

        if 'bsid2info' in stored_locations[0]:
            p_info = stored_locations[0]
        elif 'bsid2info' in stored_locations[1]:
            p_info = stored_locations[1]
        else:
            raise ValueError('Could not biosystems_info in download.')
    else:
        raise ValueError(
            'Anticipated two datasets. Check settings file.')

    # Import datasets

    df_gene = pd.read_table(
        p_gene,
        names=[   # From biosystems documentation (readme)
            'bsid_ncbi',
            'gene_ncbi',
            'score'],
        encoding='ISO-8859-1'
    )

    df_info = pd.read_table(
        p_info,
        names=[   # From biosystems documentation (readme)
            'bsid_ncbi',
            'source database of biosystem',
            'source database accession',
            'name',
            'type of biosystem',
            'taxonomic scope of biosystem',
            'NCBI taxid',
            'description of biosystem'],
        encoding='ISO-8859-1'
    )

    # Create output
    df = pd.merge(
        df_gene,
        df_info,
        left_on='bsid_ncbi',
        right_on='bsid_ncbi',
        how='left')

    # Clean output
    c = ['bsid_ncbi', 'gene_ncbi', 'score', 'source database accession',
         'source database of biosystem', 'name', 'type of biosystem',
         'taxonomic scope of biosystem']
    f = df['source database of biosystem'].notnull()   # approx. 1%-2%
    df = df.loc[f, c]
    df = df.rename(columns={
        'score': 'bsid_2_gene_score',
        'source database accession': 'accession',
        'source database of biosystem': 'bsid_source',
        'name': 'biosystem_name',
        'type of biosystem': 'biosystem_type',
        'taxonomic scope of biosystem': 'biosystem_scope'
    })
    df = df.drop_duplicates(['bsid_ncbi', 'gene_ncbi'])

    # Save output
    grouped = df.groupby('bsid_source')
    for source, d in grouped:
        d = d.reset_index(drop=True)
        o = os.path.join(
            out_p, 'biosystem_{}.gz'.format(source.lower().replace(' ', '_')))
        d.to_csv(o, index=False, compression='gzip')

    # Patch output
    _biosystems_patch_v1_1()


def _biosystems_patch_v1_1():
    """
    Will separate biosystems by taxon
    """

    p_gi = io.get_output_path('ncbi/gene_info/taxon_gene.csv')
    df_gene_info = pd.read_csv(p_gi)
    df_gene_info.loc[:, 'taxon_ncbi'] = df_gene_info.loc[
        :, 'taxon_ncbi'].astype(int)

    def split_biosystem(code):

        p = io.get_output_path(
            'ncbi/biosystems/biosystem_{}.gz'.format(code))

        p_orig_folder, p_orig_filename = os.path.split(p)
        p_subfolder = p_orig_filename.replace('.gz', '')
        p_out_folder = os.path.join(p_orig_folder, p_subfolder)
        io.ensure_presence_of_directory(p_out_folder)

        df = pd.read_csv(p)

        # Isolate ID to term as separate lookup
        id_to_term = df.loc[
            :, ['accession', 'biosystem_name']].drop_duplicates()
        if id_to_term['accession'].value_counts().max() > 1:
            raise EnvironmentError(
                'Some duplicate entries in accession to biosystem_name')
        o = os.path.join(p_out_folder, '{}_accession_to_term.csv.gz'.format(
            p_subfolder))
        id_to_term.to_csv(o, index=False, compression='gzip')

        # Save GO by taxon
        df = df.drop('biosystem_name', axis=1)
        df = pd.merge(
            df,
            df_gene_info,
            left_on='gene_ncbi',
            right_on='gene_ncbi',
            how='left')
        grouped = df.groupby('taxon_ncbi')
        for taxon, d in grouped:
            d = d.reset_index(drop=True)
            o = os.path.join(
                p_out_folder, '{}_taxon_{}.csv.gz'.format(
                    p_subfolder, int(taxon)))
            d.to_csv(o, index=False, compression='gzip')

        os.remove(p)

    databases = [
        'biocyc', 'kegg', 'lipid_maps', 'pathway_interaction_database',
        'reactome', 'wikipathways']

    for d in databases:
        split_biosystem(d)


def genbank_gene_features(taxon_id):
    """
    Creates taxon-specific comma separated tables,
    listing basic properties of genes (e.g.: GC content, length)

    If a gene exists in multiple copies (within valid replicons)
    the median will be taken. This appears to happen if a
    given gene exists on mulitple sex chromosomes

    The format of genbank genome files needs manual inspection

    Input:
        taxon_id       int, ncbi_taxon ID (e.g.: 9606 for human)

    Output:
        /out/genbanke_gene_TAXONID.csv
    """

    print('Initiating download of genomic sequence')
    p_out_download = downloader.download_genome(taxon_id, 'genomic')

    print('Initiating isolation of standard replicons')
    p_replicions = downloader.isolate_valid_replicons(p_out_download, taxon_id)

    # Initialize extraction of features
    gl = os.path.join(p_replicions, '*.gbff')
    files = glob.glob(gl)

    agg_results = []

    for currf in files:

        fi = open(currf, 'r')
        record = SeqIO.read(fi, "gb")

        for feature in record.features:

            if feature.type == 'gene':

                res = dict()

                # Obtain external reference
                # Note on 17-02-11 it appears that all taxa
                # would support Entrez Gene ID as one of multiple
                # qualifiers; Thus look specifically for Entrez
                # Gene IDs; Prior versions would sometimes use
                # alternate IDs (this has now been removed on 17-02-11)

                gene_id, found = _extract_gene_id_from_qualifier(
                    feature.qualifiers['db_xref'])

                if not(found):
                    raise EnvironmentError(
                        'Anticipated Entrez Gene ID, but failed to find one'
                        'within the reference genome for at least one '
                        'gene entry')
                else:
                    res['qualifier'] = gene_id

                # Obtain Sequence
                seq = feature.extract(record.seq)

                # Extract features of sequence
                res.update(_count_ACGT_and_length(seq))

                agg_results.append(res)

    df = pd.DataFrame(agg_results)

    # some genes exist on multiple chromosomes (e.g. X and Y)
    grped = df.groupby('qualifier')
    median_per_qualifier = grped.agg(np.median)
    median_per_qualifier.columns = [
        'Genbank__gene: {}'.format(j) for j in median_per_qualifier.columns]
    median_per_qualifier.index.name = 'gene_ncbi'
    median_per_qualifier.sort_index(inplace=True)

    n = 'genbank_gene_' + str(taxon_id) + '.csv'
    n = os.path.join('genbank_gene', n)

    o = io.get_output_path(os.path.join('genbank', n))
    print('Exporing genbank features to {} .'.format(o))
    io.ensure_presence_of_directory(o)
    io.ensure_absence_of_file(o)
    median_per_qualifier.to_csv(o)


def genbank_validated_rna_features(taxon_id):
    """
    Creates taxon-specific comma separated tables,
    listing several properties of validated RNA transcripts.
    This includes: length of RNA and coding sequence,
    nucleotides in RNA and coding sequence. Metrics of
    codon bias.

    This analys will skip RNA, for which there is little support
    according to NIH, and RNA which is tagged as partial by NIH.

    Note that for some taxa RefSeq / NIH does not provide
    any (non-computational) description of RNAs. This
    function will thus not be suitable for those taxa.

    The present script will only process transcripts for which an
    ncbi (entrez) gene is defined within the refseq genbank file.
    While this somewhat opposes the original science of science v0.1
    it allows to extend to more species with less ambiguity
    since NIH's mapping between gene IDs and transript IDs appears
    somewhat limited to few heavily studied taxa, with the origin
    of the mapping not being unambiguously defined in a 1:1 relationship
    for other taxa (or requiring online APIs whose potential bias
    would be difficult to see or undestand)

    Input:
        taxon_id       int, ncbi_taxon ID (e.g.: 9606 for human)

    Output:
        /out/genbank/validated_rna/genbank_validated_rna_TAXONID.csv
    """

    taxa_without_validated_rna = [
        511145,     # Absent from server
        11676,      # Absent from server
        6239,       # Only partial RNAs
        284812,     # Only Provisional RNAs
        386585,     # Absent from server
        224308,     # Absent from server
        83332,      # Absent from server
    ]

    if taxon_id in taxa_without_validated_rna:
        print('Taxon {} does not support validated RNA'.format(
            taxon_id))
        return

    print('Initiating download of experimental RNA sequences')
    p_out_download = downloader.download_genome(taxon_id, 'rna')

    print('Decompress RNA Sequence')  # Required for Biopython
    p_rna = downloader.ungzip(p_out_download, return_path=True)

    print('Start to extract properties of the RNA')
    # Open RNA file, and extract features
    handle = open(p_rna, mode='r')

    agg_results = []  # Initialize collector
    for record in SeqIO.parse(handle, "gb"):
        proceed = True

        if record.annotations['comment'].startswith('VALIDATED'):
            proceed = True
        elif record.annotations['comment'].startswith('REVIEWED'):
            proceed = True
        else:
            proceed = False

        if proceed:
            if 'partial mRNA' in record.description:
                proceed = False
            else:
                proceed = True

        if proceed:
            for feature in record.features:
                if feature.type == 'CDS':

                    seq_full = str(record.seq)
                    seq_cds = str(feature.extract(record.seq))

                    proceed = _appears_to_be_canonical_cds(seq_cds)

                    if proceed:
                        if not _has_only_ACGT(seq_cds):
                            proceed = False

                    if proceed:
                        # Extract Ncbi (Entrez) gene ID
                        gene_ncbi, proceed = _extract_gene_id_from_qualifier(
                            feature.qualifiers['db_xref'])

                    if proceed:
                        props = dict()
                        props.update({'rna_ncbi.version': record.id})
                        props.update({'gene_ncbi': gene_ncbi})

                        x = _get_codon_usage(seq_cds)
                        x = _insert_name_in_key(x, 'cds')
                        props.update(x)

                        x = _count_ACGT_and_length(seq_cds)
                        x = _insert_name_in_key(x, 'cds')
                        props.update(x)

                        x = _count_ACGT_and_length(seq_full)
                        x = _insert_name_in_key(x, 'full')
                        props.update(x)

                        x = _get_codon_bias(seq_cds)
                        props.update(x)

                        agg_results.append(props)

    # Save results
    handle.close()

    if not(agg_results):
        raise ValueError('Did not find any verified RNA for taxon {}'.format(
            taxon_id))
    else:
        df = pd.DataFrame(agg_results)
        # df = df.set_index(keys='rna_ncbi.version', verify_integrity=True)

        df = df.drop('rna_ncbi.version', axis=1)
        grouped = df.groupby('gene_ncbi')
        df = grouped.agg(np.median)
        df.sort_index(inplace=True)

        df.columns = [
            'Genbank_validated_RNA: {}'.format(j) for j in df.columns]

        n = 'genbank_validated_rna_' + str(taxon_id) + '.csv'
        o = io.get_output_path(os.path.join('genbank', 'validated_rna', n))

        print('Exporing genbank features of verified RNA to {} .'.format(o))
        io.ensure_presence_of_directory(o)
        io.ensure_absence_of_file(o)
        df.to_csv(o)


def genbank_cds_rna_features(taxon_id):
    """
    Creates taxon-specific comma separated tables,
    listing several properties of predicted coding sequence (CDS)
    of RNA transcripts.
    This includes: length of RNA and nucleotides

    Note that for some taxa RefSeq / NIH does not link to
    Entrez (NCBI gene ID within the genomically predicted
    RNAs – here those are skipped – though manual curation
    might be possible)

    Will only process coding sequences that look like a canonical
    coding sequence (multple of three nucleotides, ATG start, and
    TAA or TGA or TAG stop) (thus no frameshift etc.)

    Input:
        taxon_id       int, ncbi_taxon ID (e.g.: 9606 for human)

    Output:
        /out/genbank/genomic_rna/genbank_genomic_rna_TAXONID.csv
    """

    # taxa_without_genomic_rna = [
    #     511145,  # no Entrez
    #     559292,  # no Entrez
    #     3702,  # no Entrez
    #     9031,  # no Entrez
    #     6239,  # no Entrez
    #     284812,  # no Entrez
    #     386585,  # no Entrez
    #     9986,  # no Entrez
    # ]

    taxa_without_genomic_rna = [
        # 511145,  # no Entrez
        # 559292,  # no Entrez
        # 3702,  # no Entrez
        # 9031,  # no Entrez
        # 6239,  # no Entrez
        # 284812,  # no Entrez
        # 386585,  # no Entrez
        # 9986,  # no Entrez
    ]

    if taxon_id in taxa_without_genomic_rna:
        print('Taxon {} does not support coding sequence RNA'.format(
            taxon_id))
        return

    print('Initiating download of experimental RNA sequences')
    p_out_download = downloader.download_genome(
        taxon_id, 'cds_from_genomic')

    print('Decompress predicted CDS RNA Sequence')  # Required for Biopython
    p_rna = downloader.ungzip(p_out_download, return_path=True)

    print('Start to extract properties of the RNA')
    # Open RNA file, and extract features
    handle = open(p_rna, mode='r')

    agg_results = []  # Initialize collector
    for record in SeqIO.parse(handle, "fasta"):

        props = dict()

        save = True

        rna_ncbi = re.findall('^.*\|([a-zA-Z_0-9.]*)', record.id)
        if len(rna_ncbi) == 0:
            save = False
            # print(record.id)
            # raise ValueError('Did not find RefSeq ID')
        elif len(rna_ncbi) > 1:
            save = False
            # print(record.id)
            # raise ValueError('Found multiple RefSeq ID. Had expected one.')
        else:
            props.update({'rna_ncbi.version': rna_ncbi[0]})

        entrez = re.findall('.*GeneID:([0-9]*).*', record.description)
        if len(entrez) == 0:
            save = False
            # raise ValueError('Did not find Entrez gene ID')
        elif len(entrez) > 1:
            save = False
            # raise ValueError(
            #     'Found multiple Entrez gene IDs. Had expected one.')
        else:
            props.update({'gene_ncbi': entrez[0]})

        if save:
            seq = str(record.seq)
            # Only process sequences that look like complete
            # coding sequences with canonical start and stop
            # and no frameshift
            save = _appears_to_be_canonical_cds(seq)

        if save:

            x = _get_codon_usage(seq)
            x = _insert_name_in_key(x, 'predicted_cds_rna')
            props.update(x)

            x = _count_ACGT_and_length(seq)
            x = _insert_name_in_key(x, 'predicted_cds_rna')
            props.update(x)

            # Introduced in Geisen v1_1
            x = _get_codon_bias(seq)
            x = _insert_name_in_key(x, 'predicted_cds_rna')
            props.update(x)

            agg_results.append(props)

    # Save results
    handle.close()

    df = pd.DataFrame(agg_results)
    # df = df.set_index(keys='rna_ncbi.version', verify_integrity=True)
    # df.columns = ['Genbank_genomic_cds: {}'.format(j) for j in df.columns]
    # df.sort_index(inplace=True)

    df = df.drop('rna_ncbi.version', axis=1)
    grouped = df.groupby('gene_ncbi')
    df = grouped.agg(np.median)
    df.sort_index(inplace=True)

    n = 'genbank_genomic_cds_' + str(taxon_id) + '.csv'
    o = io.get_output_path(os.path.join('genbank', 'genomic_cds', n))

    print('Exporting genbank features of CDS from genomic to {} .'.format(o))
    io.ensure_presence_of_directory(o)
    io.ensure_absence_of_file(o)
    df.to_csv(o)


def genbank_genomic_rna_features(taxon_id):
    """
    Creates taxon-specific comma separated tables,
    listing several properties of genomically predicted
    RNA transcripts.
    This includes: length of RNA and nucleotides

    Note that for some taxa RefSeq / NIH does not link to
    Entrez (NCBI gene ID within the genomically predicted
    RNAs – here those are skipped – though manual curation
    might be possible)

    Input:
        taxon_id       int, ncbi_taxon ID (e.g.: 9606 for human)

    Output:
        /out/genbank/genomic_rna/genbank_genomic_rna_TAXONID.csv
    """

    # taxa_without_genomic_rna = [
    #     511145,  # no Entrez
    #     559292,  # no Entrez
    #     9913,  # no Entrez
    #     11676,  # no file on server
    #     9031,  # no Entrez
    #     6239,  # no Entrez
    #     9823,  # no Entrez
    #     8355,  # no Entrez
    #     284812,  # no Entrez
    #     9615,  # no Entrez
    #     9986,  # no Entrez
    # ]

    taxa_without_genomic_rna = [
        11676,  # no file on server  (Human immunodeficiency virus 1)
    ]

    if taxon_id in taxa_without_genomic_rna:
        print('Taxon {} does not support genomic RNA'.format(
            taxon_id))
        return

    print('Initiating download of predicted genomic RNA sequences')
    p_out_download = downloader.download_genome(
        taxon_id, 'rna_from_genomic')

    print('Decompress genomic RNA Sequence')  # Required for Biopython
    p_rna = downloader.ungzip(p_out_download, return_path=True)

    print('Start to extract properties of the RNA')
    # Open RNA file, and extract features
    handle = open(p_rna, mode='r')

    agg_results = []  # Initialize collector
    for record in SeqIO.parse(handle, "fasta"):

        props = dict()

        save = True

        rna_ncbi = re.findall('^.*\|([a-zA-Z_0-9.]*)', record.id)
        if len(rna_ncbi) == 0:
            print(record.id)
            save = False
            # raise ValueError('Did not find RefSeq ID')
        elif len(rna_ncbi) > 1:
            print(record.id)
            save = False
            raise ValueError('Found multiple RefSeq ID. Had expected one.')
        else:
            props.update({'rna_ncbi.version': rna_ncbi[0]})

        entrez = re.findall('.*GeneID:([0-9]*).*', record.description)
        if len(entrez) == 0:
            save = False
            # raise ValueError('Did not find Entrez gene ID')
        elif len(entrez) > 1:
            save = False
            # raise ValueError(
            #     'Found multiple Entrez gene IDs. Had expected one.')
        else:
            props.update({'gene_ncbi': entrez[0]})

        if save:
            seq = str(record.seq)

            x = _count_ACGT_and_length(seq)
            x = _insert_name_in_key(x, 'predicted_genomic_RNA')
            props.update(x)

            agg_results.append(props)

    # Save results
    handle.close()

    df = pd.DataFrame(agg_results)
    # df = df.set_index(keys='rna_ncbi.version', verify_integrity=True)
    # df.columns = ['Genbank_genomic_RNA: {}'.format(j) for j in df.columns]
    # df.sort_index(inplace=True)

    df = df.drop('rna_ncbi.version', axis=1)
    grouped = df.groupby('gene_ncbi')
    df = grouped.agg(np.median)
    df.sort_index(inplace=True)

    n = 'genbank_genomic_rna_' + str(taxon_id) + '.csv'
    o = io.get_output_path(os.path.join('genbank', 'genomic_rna', n))

    print('Exporting genbank features of genomic RNA to {} .'.format(o))
    io.ensure_presence_of_directory(o)
    io.ensure_absence_of_file(o)
    df.to_csv(o)


def gene_info():
    """
    Will download gene_info from ncbi, and format it for further usage
    - adhere to science of science naming conventions
    - remove placeholders (no genes, but required for submision of GeneRIFs)
    - save individual taxa as flat files (compressed)
       ( to account for pandas not being able to save gene_info in
        HDF5 https://github.com/pandas-dev/pandas/issues/15318 )
    - save full table as gene_info_full.gz
    """

    # Initialize Absent Output
    out_p = 'ncbi/gene_info'
    out_p = io.get_output_path(out_p)
    io.ensure_presence_of_directory(out_p)

    number_of_files = io.check_number_of_files_in_directory(out_p, '*.gz')
    if number_of_files > 0:
        raise EnvironmentError(
            'At least one output file is already present.')

    # Download Dataset
    stored_locations = downloader.download_data_set('gene_info')

    if len(stored_locations) == 1:
        p = stored_locations[0]
    else:
        raise ValueError(
            'Anticipated a single dataset. Check settings file.')

    print('Importing gene_info file')
    df = pd.read_table(p, sep='\t', header=0)

    print('Cleaning gene_info')
    df = df.rename(
        columns={
            '#tax_id': 'taxon_ncbi',
            'GeneID': 'gene_ncbi',
            'Symbol': 'symbol_ncbi',
            'Symbol_from_nomenclature_authority': 'symbol_authority'})
    f = df['symbol_ncbi'] == 'NEWENTRY'  # Used for submitted GeneRIFs
    df = df.loc[~f, :]

    # Safety check for unique identity of genes
    # (note: that the following code paragraph could be adjusted
    # to enable an omission of ambiguous entries)
    orig_amount_of_genes = df.shape[0]
    df = df.drop_duplicates('gene_ncbi', keep=False)
    if df.shape[0] != orig_amount_of_genes:
        raise ValueError('Gene_info lists some genes several times')

    print('Exporting gene_info file')
    out_joined_file = os.path.join(out_p, 'gene_info_full.gz')
    df.to_csv(out_joined_file, index=False, compression='gzip')

    print('Exporting taxon_and_gene file')
    out_joined_file = os.path.join(out_p, 'taxon_gene.csv')
    df.loc[:, ['taxon_ncbi', 'gene_ncbi']].to_csv(out_joined_file, index=False)

    print('Exporting gene_info of individual taxa. Wait approx. 5min.')
    grouped = df.groupby('taxon_ncbi')

    for taxon, d in grouped:
        d = d.reset_index(drop=True)
        o = os.path.join(out_p, 'gene_info_taxon_{}.gz'.format(taxon))
        d.to_csv(o, index=False, compression='gzip')


def gene2ensembl():
    """
    Will download gene2ensembl, and format it for further usage
    - adhere to science of science naming conventions
    """

    # Initialize Absent Output
    out = 'ncbi/gene2ensembl.gz'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)

    # Download Dataset
    stored_locations = downloader.download_data_set('gene2ensembl')

    if len(stored_locations) == 1:
        p = stored_locations[0]
    else:
        raise ValueError(
            'Anticipated a single dataset. Check settings file.')

    df = pd.read_table(p, sep='\t', header=0)
    df = df.rename(columns={
        '#tax_id': 'taxon_ncbi',
        'GeneID': 'gene_ncbi',
        'Ensembl_gene_identifier': 'gene_ensembl',
        'RNA_nucleotide_accession.version': 'rna_ncbi.version',
        'Ensembl_rna_identifier': 'rna_ensembl',
        'protein_accession.version': 'protein_ncbi.version',
        'Ensembl_protein_identifier': 'protein_ensembl'})
    df = df.drop_duplicates()
    df.to_csv(out, index=False, compression='gzip')


def gene2go():
    """
    Will download gene2go, and format it for further usage
    - adhere to science of science naming conventions
    - add meta-tags:
        - Negating support: Support is negating
        - Any negating support: For given combination of GO and gene there
          is both positive support (note: can be possible according to
          GO documentation; e.g.: if in different tissues)
        - Temporary Evidence (max. 12m): Several GO annotation tags
          indicate that the support should should be temporary (see
          documentation; note that in contrast to common believe the tags
          of many annotation tags do not simply indicate, whether
          data was computationally acquired
        - Unmapped Evidence: Corresponds to places where either MedLine,
          or GO do not have a mapped evidence code (which can be
          a problem of Medline, or ancient GO entries)
    - prevent duplicates
    """

    # Initialize Absent Output
    out = 'ncbi/gene2go.gz'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)

    # Download Dataset
    stored_locations = downloader.download_data_set('gene2go')

    if len(stored_locations) == 1:
        p = stored_locations[0]
    else:
        raise ValueError(
            'Anticipated a single dataset. Check settings file.')

    print('Opening gene2go')
    df = pd.read_table(p, sep='\t', header=0)

    print('Extracing additional metadata')
    # Negating Evidence, stars with NOT
    list_of_unique_qualifiers = df['Qualifier'].unique()
    negative_qualifiers = [
        j for j in list_of_unique_qualifiers if j.startswith('NOT')]
    has_negating_qualifier = df['Qualifier'].isin(negative_qualifiers)
    has_negating_evidence = df['Evidence'].isin(['ND'])
    is_negating = (has_negating_qualifier) | (has_negating_evidence)
    df.loc[:, 'Negating support'] = is_negating

    # Approx. 1% of negated records also carry a positive record
    df.loc[:, 'tmp'] = df['GeneID'].astype(str) + '_' + df['GO_ID']
    u = df.loc[is_negating, 'tmp'].unique()
    has_a_negating = df['tmp'].isin(u)
    df['Any negating support'] = has_a_negating
    df = df.drop('tmp', axis=1)

    # Expiring Evidence
    expiring_evidence = ['IEA', 'RCA']      # Expiring
    df.loc[:, 'Temporary Evidence (max. 12m)'] = df[
        'Evidence'].isin(expiring_evidence)

    # Unmapped Evidence
    unmapped_evidence = ['-', 'NR']         # Unmapped
    df.loc[:, 'Unmapped Evidence'] = df['Evidence'].isin(unmapped_evidence)

    print('clean data')
    # Clean up
    df = df.rename(columns={
        '#tax_id': 'taxon_ncbi',
        'GeneID': 'gene_ncbi'})
    df = df.drop_duplicates([
        'taxon_ncbi', 'gene_ncbi', 'GO_ID', 'Evidence', 'Qualifier'])

    print('Exporting data')
    df.to_csv(out, index=False, compression='gzip')

    # Patch
    _gene2go_patch_v1_1()


def _gene2go_patch_v1_1():
    """
    Patch to geisen v1_1:
    - Splits gene2go to taxa, and isolates lookup of id and term
        (to allow faster loading)
    - removes deprecated gene2go from geisen_v1
    """

    # Initiate output folder
    out_p = 'ncbi/gen2ego'
    out_p = io.get_output_path(out_p)
    io.ensure_presence_of_directory(out_p)

    # Get gene2go from geisen_v1
    p_orig = 'ncbi/gene2go.gz'
    p_orig = io.get_output_path(p_orig)
    df = pd.read_csv(p_orig)

    # Isolate ID to term as separate lookup
    id_to_term = df.loc[:, ['GO_ID', 'GO_term']].drop_duplicates()
    if id_to_term['GO_ID'].value_counts().max() > 1:
        raise EnvironmentError(
            'Some duplicate entries in GO_ID to GO_term')
    o = os.path.join(out_p, 'go_id_to_term.csv.gz')
    id_to_term.to_csv(o, index=False, compression='gzip')

    # Save GO by taxon
    df = df.drop('GO_term', axis=1)
    grouped = df.groupby('taxon_ncbi')
    for taxon, d in grouped:
        d = d.reset_index(drop=True)
        o = os.path.join(
            out_p, 'gene2go_taxon_{}.csv.gz'.format(taxon))
        d.to_csv(o, index=False, compression='gzip')

    # Cleanup: remvoe file of geisen_v1
    os.remove(p_orig)


def gene2pubmed():
    """
    Will download gene2pubmed from ncbi, and format it for further usage
    - adhere to science of science naming conventions
    - save as hdf5, where every column could be used for filtering data
    """

    # Initialize Absent Output
    out = 'ncbi/gene2pubmed.h5'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)

    # Download Dataset
    stored_locations = downloader.download_data_set('gene2pubmed')

    if len(stored_locations) == 1:
        p = stored_locations[0]
    else:
        raise ValueError(
            'Anticipated a single dataset. Check settings file.')

    # Impose Science of Science convention
    df = pd.read_table(p, sep='\t', header=0)
    df = df.rename(
        columns={
            '#tax_id': 'taxon_ncbi',
            'GeneID': 'gene_ncbi',
            'PubMed_ID': 'pubmed_id'})
    df = df.drop_duplicates()

    # Export to HDF5 for selective loading
    df.to_hdf(
        out,
        'table',
        mode='w',
        append=True,
        data_columns=df.columns)


def generifs_basic():
    """
    Will download generifs_basic from ncbi, and format it for further usage
    - adhere to science of science naming conventions
    - remove all, but the youngest, entry, if duplicates
    """

    # Initialize Absent Output
    out = 'ncbi/generifs_basic.gz'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)
    print(out)

    # Download Dataset
    stored_locations = downloader.download_data_set('generifs_basic')

    # Safety Check
    if len(stored_locations) == 1:
        p = stored_locations[0]
    else:
        raise ValueError(
            'Anticipated a single dataset. Check settings file.')

    # Format
    df = pd.read_table(p, sep='\t', header=0, low_memory=False)
    df = df.rename(columns={
        '#Tax ID': 'taxon_ncbi',
        'Gene ID': 'gene_ncbi'
    })

    # Remove duplicate entries
    df = df.sort_values(by='last update timestamp', ascending=True)
    df = df.reset_index(drop=True)
    df = df.drop_duplicates(
        ['taxon_ncbi', 'gene_ncbi', 'PubMed ID (PMID) list', 'GeneRIF text'],
        keep='last')
    df = df.drop('last update timestamp', axis=1)
    df = df.sort_values(
        ['taxon_ncbi', 'gene_ncbi', 'PubMed ID (PMID) list', 'GeneRIF text'])
    df = df.reset_index(drop=True)

    # Export
    df.to_csv(out, index=False, compression='gzip')


def genetic_testing_registry():
    """
    Downloads datasets from genetic testing registry
    """

    df_settings = settings.retreive_general_download_settings()

    out_folder = io.get_output_path('nih/gtr')
    io.ensure_presence_of_directory(out_folder)

    p_in = df_settings.loc[
        'genetic_testing_registry_test_version',
        'location_on_server']
    _, basename = os.path.split(p_in)
    p_out = os.path.join(out_folder, basename)
    downloader.download(p_in, p_out)

    p_in = df_settings.loc[
        'genetic_testing_registry_test_condition',
        'location_on_server']
    _, basename = os.path.split(p_in)
    p_out = os.path.join(out_folder, basename)
    downloader.download(p_in, p_out)


def uniprot_databases():
    """
    Downloads uniprot databases. Note that this can be quite
    lengthy (taking up to 1.5h on a 40 MBps connetion)

    Create complete download of uniprot databases
    - Swissprot: manually curated
    - TREMBL: computationally predicted
    """

    # SWISSPROT: high confidence manually curated
    print('Initiating download of swissprot sequences')
    p_down = downloader.download_data_set('swissprot')

    if len(p_down) > 1:
        raise ValueError(
            'Download does not appear unambiguous'
            'Please only list a single swissport in the settings file.')
    else:
        p_down = p_down[0]

    print('Decompress swissprot Sequences')  # Required for Biopython
    du = downloader.ungzip(p_down, return_path=True)

    # TREMBL: computationally curated
    print(
        'Initiating download of trebml sequences.'
        'Note that this might require several hours.')
    p_down = downloader.download_data_set('trembl')

    if len(p_down) > 1:
        raise ValueError(
            'Download does not appear unambiguous'
            'Please only list a single trebml in the settings file.')
    else:
        p_down = p_down[0]

    print('Decompress trebml Sequences')  # Required for Biopython
    du = downloader.ungzip(p_down, return_path=True)
    del du


def uniprot_id_mapper():
    """
    Prepare Uniprot's ID mapper
    - download
    - dismisses most of the naming sources, leaving only Entrez
    - discard all entries which have a nan for taxon or gene or uniprot
    - will save mapper as hdf5 that is indexable by any column
    - enforce science of science nomenclature
    - columns are: taxon_ncbi, gene_ncbi, protein_uniprot
    """

    # Initialize Absent Output
    out_p = 'uniprot/uniprot_id_mapper.h5'
    out_p = io.get_output_path(out_p)
    io.ensure_presence_of_directory(out_p)
    io.ensure_absence_of_file(out_p)

    p_raw_download = downloader.download_data_set('uniprot_ids')

    if len(p_raw_download) > 1:
        raise ValueError(
            'Download does not appear unambiguous'
            'Please only list a single uniprot_ids in the settings file.')
    else:
        p_raw_download = p_raw_download[0]

    print('Loading downloaded Uniprot ID mapper')
    df = pd.read_table(
        p_raw_download,
        sep='\t',
        names=[   # from README file, largely redundant with NIH gene info
            'UniProtKB-AC', 'UniProtKB-ID', 'GeneID (EntrezGene)', 'RefSeq',
            'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc',
            'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL',
            'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO',
            'Additional PubMed'],
        usecols=['UniProtKB-AC', 'GeneID (EntrezGene)', 'NCBI-taxon'],
        low_memory=True)

    print('Tidying Uniprot ID mapper.')
    df = df.rename(columns={
        'UniProtKB-AC': 'protein_uniprot',
        'GeneID (EntrezGene)': 'gene_ncbi',
        'NCBI-taxon': 'taxon_ncbi'})

    df = df.dropna()

    f = df['gene_ncbi'].str.contains(';')
    df = df.loc[~(f == True), :]
    df.loc[:, 'gene_ncbi'] = df.loc[:, 'gene_ncbi'].astype(int)

    df = df.drop_duplicates()

    print('Saving Uniprot ID mapper to hdf5')
    df.to_hdf(
        out_p,
        'table',
        mode='w',
        append=True,
        data_columns=df.columns)


def taxdmp():
    """
    Will download taxdmp from ncbi, and format it for further usage
    Will only use the following files contained in taxdmp:
    - names: Names of taxa
    - nodes: Would allow reconstruction, contains level of hierarchy
    Processing of those files will
    - adhere to science of science naming conventions
    - add column names from documentation (manual)
    - save individual taxa as hdf5, searchable by taxon_ncbi

    """

    # Download Dataset
    stored_locations = downloader.download_data_set('taxdmp')

    # Safety Check
    if len(stored_locations) == 1:
        p = stored_locations[0]
    else:
        raise ValueError(
            'Anticipated a single dataset. Check settings file.')

    # unzip into same folder
    downloader.unzip(p)

    # Define Supprt function, for reading taxonomy files
    def _read_ncbi_taxonomy_file(filepath, column_names):
        """
        Reads an ncbi taxonomy files (note which have unusual line breaks)
            filepath        path to taxonomy file
            column_names    list containing names of columns:
                            must be manually inferred from NCBI's documentation
                            (readme.txt on their server)
        """

        if not os.path.exists(filepath):
            raise EnvironmentError(
                'Could not find ' + filepath + ' . Please check spelling.\n')
        else:
            df = pd.read_table(
                filepath, sep='\t\|\t',
                header=None,
                names=column_names,
                engine='python')
            df[df.columns[-1]] = df[
                df.columns[-1]].map(lambda x: x.replace("\t|", ""))
        return df

    # Import names.dmp
    p, _ = os.path.split(p)

    p_names = os.path.join(p, 'names.dmp')

    column_names = [     # from readme.txt:
        'tax_id',        # the id of node associated with this name
        'name_txt',      # name of taxon
        'unique name',   # the unique variant of this name if name not unique
        'name class']    # synonym, common name, ...

    df = _read_ncbi_taxonomy_file(p_names, column_names)
    df = df[df['name class'] == 'scientific name']

    if not all(df['tax_id'].value_counts() == 1):
        raise ValueError(
            'Scientific names are not unambiguous for individual taxon.' +
            'Please check if NIH changed format of files on their server')

    df = df.rename(columns={
        'tax_id': 'taxon_ncbi',
        'name_txt': 'taxon_name'})
    df = df.drop_duplicates('taxon_ncbi', keep=False)

    out = 'ncbi/taxon_names.h5'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)

    df.to_hdf(
        out,
        'table',
        mode='w',
        append=True,
        data_columns=['taxon_ncbi'])

    # Import nodes.dmp
    p_names = os.path.join(p, 'nodes.dmp')

    column_names = [        # from readme.txt
        'tax_id',           # node id in GenBank taxonomy database
        'parent tax_id',    # parent node id in GenBank taxonomy database
        'rank',             # superkingdom, kingdom, ...
        'embl code',        # locus-name prefix; not unique
        'division id',      # see division.dmp file
        'inherited div flag  (1 or 0)',  # inherits division from parent
        'genetic code id',               # see gencode.dmp file
        'inherited GC  flag  (1 or 0)',  # genetic code from parent
        'mitochondrial genetic code id',  # see gencode.dmp file
        'inherited MGC flag  (1 or 0)',   # mitochondrial gencode from parent
        'GenBank hidden flag (1 or 0)',   # suppressed in GenBank entry lineage
        'hidden subtree root flag (1 or 0)',  # no sequence data yet
        'comments']          # free-text comments and citations

    df = _read_ncbi_taxonomy_file(p_names, column_names)
    df = df.rename(columns={'tax_id': 'taxon_ncbi'})
    df = df.drop_duplicates('taxon_ncbi', keep=False)

    out = 'ncbi/taxon_nodes.h5'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)

    df.to_hdf(
        out,
        'table',
        mode='w',
        append=True,
        data_columns=['taxon_ncbi'])


'''
############### TIER 2 ###########################

The generation of these datasets requires that TIER 1 functions
have been executed first (and the respective TIER 1 datasets
have been created succesfully)

aminoacid_features
chromosomes
homologene
interpro
isolate_taxa_protein_fasta

'''


def chromosomes():
    """
    Extracts information about chromosomes

    Requirement:
        gene_info

    Output:
        ncbi/chromosomes.h5  with taxon_ncbi being index-able
    """

    out = 'ncbi/chromosomes.h5'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)

    p_in = io.get_output_path(
        os.path.join('ncbi', 'gene_info', 'gene_info_full.gz'))

    if not os.path.exists(p_in):
        raise ValueError(
            'Did not find gene_info_full.gz'
            'Please run gene_info before')

    df = pd.read_csv(
        p_in, usecols=['gene_ncbi', 'taxon_ncbi', 'chromosome'])

    df = df.drop_duplicates()

    df.to_hdf(
        out,
        'table',
        mode='w',
        append=True,
        data_columns=['taxon_ncbi'])

    # Patch
    _chromosomes_patch_v1_1()


def _chromosomes_patch_v1_1():
    """
    Separates information about chromosomes to per-taxon
    gnu zipped csv. The primary purpose is to reduce the overall file
    size of the project.
    """

    p = 'ncbi/chromosomes.h5'
    p = io.get_output_path(p)
    out_p = io.get_output_path('ncbi/chromosomes')
    io.ensure_presence_of_directory(out_p)

    df = pd.read_hdf(p, 'table')
    grouped = df.groupby('taxon_ncbi')

    for taxon, d in grouped:
        d = d.reset_index(drop=True)
        o = os.path.join(out_p, 'chromosomes_taxon_{}.csv.gz'.format(
            int(taxon)))
        d.to_csv(o, index=False, compression='gzip')

    os.remove(p)


def flybase_expression():
    """
    Flybase maintains a well maintained collection of gene expresssion
    within Drosophila melanogaster.

    Will separate individual datasets, such as modENCODE tissue sequences
    (note: not all gene expression data contained in flybase
    will be separated)
    """

    p_out = 'flybase'
    p_out = io.get_output_path(p_out)
    io.ensure_absence_of_file(p_out)
    io.ensure_presence_of_directory(p_out)

    p_down = downloader.download_data_set('flybase_expression')

    if len(p_down) > 1:
        raise ValueError(
            'Download does not appear unambiguous'
            'Please only list a single flybase_expression'
            'in the settings file.')
    else:
        p_down = p_down[0]

    df = pd.read_table(p_down, header=5)
    df = df.iloc[:-1:, :]  # remove download date

    unique_flybase_ids = df['FBgn#'].unique()

    datasets_to_consider = {
        'stages': 'modENCODE_mRNA-Seq_U',
        'tissues': 'modENCODE_mRNA-Seq_tissues',
        'treatments': 'modENCODE_mRNA-Seq_treatments',
        'cells': 'modENCODE_mRNA-Seq_cell.B',
    }

    for k, d in datasets_to_consider.items():
        f = df['Parent_library_name'] == d

        dff = df[f]

        agg = pd.DataFrame(index=unique_flybase_ids)

        unique_datasets = dff['RNASource_name'].unique()

        for t in unique_datasets:
            f = dff['RNASource_name'] == t
            e = dff[f].set_index('FBgn#', verify_integrity=True)
            agg[t] = e['RPKM_value']

        agg.index.name = 'gene_ensembl'  # true on: 170308
        agg.columns = ['flybase__{}: {}'.format(k, j) for j in agg.columns]

        # treat 0 as non defined / very low expression
        # see http://flybase.org/reports/FBrf0221009.html
        f = agg == 0
        agg[f] = np.nan
        agg = agg.apply(np.log10)

        def add_tag_to_column_name(df):
            prefix = 'log10_RNA_expression'
            df.columns = ['{}_{}'.format(prefix, j) for j in df.columns]
            return df

        agg = add_tag_to_column_name(agg)
        agg_entrez = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
            df=agg,
            taxon_id=7227)  # Drosophila melanogaster

        v = 'flybase_expression_{}'.format(k)  # label of datase
        _save_orig_and_ncbi_gene_mapped_tables(
            p_dir=p_out,
            filebase=v,
            df_orig=agg,
            df_ncbi=agg_entrez)


def gerstein_expression():
    """
    The Gerstein laboratory maintains a partially outdated (in its locus
    annotation) gene expression dataset from modENCODE. Sample
    annotation is in a terrible state, and, although modENCODE help
    would forward / CC them and direclty adress them, they were not
    able to give proper annotation for many samples (only their
    short abbreviation) (or simply ignored that request...)

    Anyways, this funciton will load the official modENCODE
    high level represenation for C. elegans
    """

    p_out = 'gerstein'
    p_out = io.get_output_path(p_out)
    io.ensure_absence_of_file(p_out)
    io.ensure_presence_of_directory(p_out)

    p_down = downloader.download_data_set('gerstein_expression')

    if len(p_down) > 1:
        raise ValueError(
            'Download does not appear unambiguous'
            'Please only list a single gerstein_expression'
            'in the settings file.')
    else:
        p_down = p_down[0]

    df = _load_gerstein_expression(p_down)
    df_entrez = mapper.locustag_2_gene_ncbi_unambiguously(
        df, taxon_id=6239)

    v = 'gerstein_expression_'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=df,
        df_ncbi=df_entrez)


def homologene():
    """
    Will download homologene from ncbi, and format it for further usage.
    Note that homologene maps can map the same gene to a different
    taxon, with MedLine usually showing the specific strain and
    homoogene the taxonomy ID of the species
    - adhere to science of science naming conventions
    - map to taxonomy ID as used in Medline
    - remove entries / taxa that can not be mapped to Medline (note:
        if amount of genes is below a given threshold, default 100)
    - remove most meta-columns of homologene
    """

    # Setting: some genes in homologene do map to strain that
    # has not been used in MedLine -> require certain amount of
    # genes in homlologene to avoid keeping these taxa
    minimal_amount_of_required_genes = 100

    # Initialize Absent Output
    out = 'ncbi/homologene.gz'
    out = io.get_output_path(out)
    io.ensure_absence_of_file(out)
    io.ensure_presence_of_directory(out)
    print(out)

    # Prepare gene_info, that will be required for processing
    p = os.path.join(io.get_output_path(   # here: use ouput of others as input
        'ncbi/gene_info'), 'taxon_gene.csv')
    df_gene_info = pd.read_csv(p, usecols=['taxon_ncbi', 'gene_ncbi'])
    df_gene_info = df_gene_info.rename(
        columns={'taxon_ncbi': 'taxon_gene_info'})

    # Download Dataset (homologene)
    stored_locations = downloader.download_data_set('homologene')

    # Safety Check
    if len(stored_locations) == 2:
        if stored_locations[0].endswith('homologene.data'):
            p = stored_locations[0]
        elif stored_locations[1].endswith('homologene.data'):
            p = stored_locations[1]
        else:
            raise ValueError('Could not find homoologene.data.')
    else:
        raise ValueError(
            'Anticipated two datasets. Check settings file.')

    df = pd.read_table(p, sep='\t', names=[
        'homologene_group', 'taxon_ncbi', 'gene_ncbi',
        'symbol_ncbi', 'prot_gi', 'protein_ncbi.version'])
    df = df.rename(columns={'taxon_ncbi': 'taxon_homologene'})
    df = df.drop_duplicates()

    d = pd.merge(
        df,
        df_gene_info,
        left_on='gene_ncbi',
        right_on='gene_ncbi',
        how='left')

    genes_in_taxa = d['taxon_gene_info'].value_counts()
    allowed_taxa = genes_in_taxa[
        genes_in_taxa >= minimal_amount_of_required_genes].index

    d = d[d['taxon_gene_info'].isin(allowed_taxa)]
    d['taxon_gene_info'] = d['taxon_gene_info'].astype(int)
    d = d.rename(columns={'taxon_gene_info': 'taxon_ncbi'})
    d = d.loc[:, ['homologene_group', 'taxon_ncbi', 'gene_ncbi']]
    d.to_csv(out, index=False, compression='gzip')


def interpro():
    """
    Places the interpro database in a format with
    which one can easily work with.

    Interpro is a large database that unites several databases
    on proteins. Note that individual databases can have
    entries which refer to the same interpro ID (and thus
    the same property of proteins)

    DBCODE   |    name
    ------------------------------------------------
    cd          Conserved Domain
    G3          CATH Superfamily
    MF          Hamap
    PD          ProDom
    PF          Pfam
    PI          Protein Information Resource
    PR          PRINTS
    PS          Prosite
    PT          Panther
    SF          Structure Function Linkage Database
    SM          SMART
    SS          SUPERFAMILY
    TI          TIGR

    Requirement:
    - uniprot_id_mapper must have been run first

    Output:
    - interpro/interpro_DBCODE.h5 files index-able by taxon
    - interpro/interpro_names.csv.gz    linking interpro id to name

    """

    # Checking for absent output
    out_p = 'interpro/interpro_names.csv.gz'
    out_p = io.get_output_path(out_p)
    io.ensure_presence_of_directory(out_p)
    io.ensure_absence_of_file(out_p)

    p_down = downloader.download_data_set('interpro')

    print('Initiate download of Interpro')
    if len(p_down) > 1:
        raise ValueError('Did not find an unambiguous download.')
    else:
        p_down = p_down[0]

    # Make datasets

    print('Start to extract names of interpro codes')

    df_n = pd.read_table(
        p_down,
        names=[
            'protein_uniprot', 'interpro_id', 'interpro_name',
            'source_id', 'undefined_A', 'undefined_B'],
        usecols=['interpro_id', 'interpro_name'])
    df_n = df_n.drop_duplicates()
    df_n.to_csv(out_p, index=False, compression='gzip')
    del df_n

    print('Start to extract databases within Interpro')

    df_i = pd.read_table(
        p_down,
        names=[
            'protein_uniprot', 'interpro_id', 'interpro_name',
            'source_id', 'undefined_A', 'undefined_B'],
        usecols=['protein_uniprot', 'interpro_id', 'source_id'])

    source_codes = [  # Manually curated on 170215 by TS, see docstring
        'cd', 'G3', 'MF', 'PD', 'PF', 'PI', 'PR',
        'PS', 'PT', 'SF', 'SM', 'SS', 'TI']

    df_u = _load_uniprot_mapper()

    # process individual source databases separately (to save RAM)
    for source_db in source_codes:
        f = df_i['source_id'].str.startswith(source_db)
        df_if = df_i[f]

        df = pd.merge(
            df_if,
            df_u,
            left_on='protein_uniprot',
            right_on='protein_uniprot',
            how='inner')

        df = df.drop_duplicates()

        p = io.get_output_path('interpro/interpro_{}.h5'.format(source_db))

        df.to_hdf(
            p,
            'table',
            mode='w',
            append=True,
            data_columns=['taxon_ncbi'])

    print('Finished exporting Interpro')


def isolate_taxa_protein_fasta():
    """
    Takes reference protein databases and extracts taxon specific
    FASTA files.

    - Sequentially processes swissport (experimentally verified),
      then trembl (computationally predicted).
    - Uses taxa defined in stettings file cfg/genome_links.txt
    - Creates joined records of all proteins of one taxon and
      recoreds separted into separte files (requirements of different
      bioinformatic tools)

    Output:
        virtualexchange/batches

    """

    # List of manually defiend taxa, that shall be subjected
    # to analyis of genomic features (170228: This corresponds
    # to the 20 most studied taxa)
    s = settings.retreive_genome_download_settings()
    reference_taxa = list(s.keys())

    # Uniprot databases that shall be use
    reference_dbs = ['swissprot', 'trembl']

    for db in reference_dbs:
        for taxon in reference_taxa:

            print('Start processing', taxon, 'of', db, 'database.')

            _isolate_taxon_protein_fasta_individual(taxon, db)
            _isolate_taxon_protein_fasta_joined(taxon, db)


"""
################ TIER 3 FUNCTIONS ################

aminoacid_features

note:
also run bash scripts in virtual machine for:

radar
seg
signalp


"""


def aminoacid_features(taxon_id, protein_db):

    """
    Extracts properties of aminoacids (indluding derived properites,
    such as gravy and isolectric point), and saves result as a comma
    separated file.

    Requirement:
        species specific fasta file containing all proteins

    Input:
        taxon_id       integer, ncbi_taxon_id (e.g.: 9606 for Homo sapiens)

    Output:
        /aminoacids/aminoacids_TAXONID.csv
    """

    # For input of batch jobs (run on dedicated virtual box)
    ext = os.path.join(
        'virtualexchange',
        'batches',
        protein_db,
        'all_in',
        str(taxon_id),
        'all.fasta')
    pa = io.get_internal_path(ext)

    species_of_interest = mapper.ncbi_taxon_2_uniprot_taxon(taxon_id)
    species_filter = 'OS=' + species_of_interest

    joined_per_protein = []

    for seq_record in SeqIO.parse((pa), 'fasta'):
        if species_filter in seq_record.description:

            p = _get_aminoacid_features(seq_record.seq)
            d = _extract_uniprot_from_uniprot_fasta_header(
                seq_record.description)

            d.update(p)
            joined_per_protein.append(d)

    if len(joined_per_protein) == 0:
        print('Did not find entries for', taxon_id, 'in', protein_db)
        return

    df = pd.DataFrame(joined_per_protein)
    df = df.set_index('protein_uniprot', verify_integrity=True)
    df.columns = ['Aminoacids_{}: {}'.format(
        protein_db, j) for j in df.columns]

    n = 'aminoacids_{}_{}_orig_id.gz'.format(protein_db, taxon_id)
    n = os.path.join('aminoacids', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df.to_csv(n, index=True, compression='gzip')

    # Create dataset aggregated by ncbi gene ID
    df_g = mapper.uniprot_protein_2_gene_ncbi(df, 'median')
    n = 'aminoacids_{}_{}_gene_id.gz'.format(protein_db, taxon_id)
    n = os.path.join('aminoacids', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df_g.to_csv(n, index=True, compression='gzip')


"""
################ TIER 4 FUNCTIONS ################


retreive_seg      (needs output generated by proteinalaysis @ virtualbox)
retreive_signalp  (needs output generated by proteinalaysis @ virtualbox)
retreive_radar    (needs output generated by proteinalaysis @ virtualbox)
"""


def retreive_radar(taxon_id, protein_db):
    """
    Obtains the best ranked RADAR prediction
    (intra-molecule similarity) from gene-specific
    computations of one taxon,
    and saves them as a comma-separated file

    Requirement:
        Predictions run through RADAR

    Input:
        taxon_id    integer, ncbi_taxon_id (e.g.: 9606 for Homo sapiens)

    Output:
        secondary/radar/radar_TAXONID.csv
        count_wrongs    integer with number of
        files that could not be processed

    """

    ext = os.path.join(
        'virtualexchange',
        'batches',
        protein_db,
        'radar',
        str(taxon_id),
        '*.radar')
    m = io.get_internal_path(ext)
    files = glob.glob(m)

    count_wrongs = 0
    agg_results = []

    for currf in files:

        fi = open(currf, 'r', encoding='ISO-8859-1')

        lines = fi.readlines()

        res = dict()

        (drive, path) = os.path.splitdrive(currf)
        (path, file) = os.path.split(path)
        res['protein_uniprot'] = os.path.splitext(file)[0]

        if len(lines) < 2:

            res['number_of_repeats'] = np.nan
            res['total_score'] = np.nan
            res['length'] = np.nan

            count_wrongs += 1

        elif 'No repeats' in lines[1]:
            res['number_of_repeats'] = 0
            res['total_score'] = 0
            res['length'] = 0

        elif '---------------------------' in lines[1]:

            # corresponds to position of top ranked
            # RADAR prediticion (hightest level)
            t = lines[3]
            m = '^ *([0-9]*)\| *([.0-9]*)\| *([0-9]*)\|.*\|.*\|.*\|.*$'

            ma = re.search(m, t)
            res['number_of_repeats'] = int(ma.group(1))
            res['total_score'] = float(ma.group(2))
            res['length'] = int(ma.group(3))

        else:

            res['number_of_repeats'] = np.nan
            res['total_score'] = np.nan
            res['length'] = np.nan

            count_wrongs += 1

        fi.close()

        agg_results.append(res)

    print('Taxon {} : Could not process {} RADAR entries'.format(
        taxon_id, count_wrongs))

    if len(agg_results) == 0:
        print('Did not find entries for', taxon_id, 'in', protein_db)
        return

    df = pd.DataFrame(agg_results)
    df = df.set_index(['protein_uniprot'])
    df.columns = ['Radar_{}_best: {}'.format(
        protein_db, j) for j in df.columns]
    df = df.dropna()       # delete entries with nan

    n = 'radar_{}_{}_orig_id.csv.gz'.format(protein_db, taxon_id)
    n = os.path.join('radar', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df.to_csv(n, index=True, compression='gzip')

    df_g = mapper.uniprot_protein_2_gene_ncbi(df, 'median')
    n = 'radar_{}_{}_gene_id.csv.gz'.format(protein_db, taxon_id)
    n = os.path.join('radar', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df_g.to_csv(n, index=True, compression='gzip')

    return


def retreive_seg(taxon_id, protein_db):
    """
    Obtains the SEG (Sequence Complexity) from
    gene-specific computations of one taxon,
    and saves them as a comma-separated file

    For mapping to ncbi_gene_id, the median will be used

    Requirement:
        Predictions run through SEG

    Input:
        taxon_id        int; ncbi taxonomy ID
        protein_db      str; used uniprot database
                            e.g. swissprot or trembl

    Output:
        seg/            gzipped csv with SEG
    """

    ext = os.path.join(
        'virtualexchange',
        'batches',
        protein_db,
        'seg',
        '{}.seg'.format(taxon_id))
    pa = io.get_internal_path(ext)

    agg = []
    for seq_record in SeqIO.parse((pa), 'fasta'):

        res, labels = _get_seg_metrics(seq_record)
        agg.append(res)

    if len(agg) == 0:
        print('Did not find entries for', taxon_id, 'in', protein_db)
        return

    df = pd.DataFrame(agg, columns=labels)
    df = df.set_index('protein_uniprot', verify_integrity=True)

    df.columns = ['SEG_{}: {}'.format(protein_db, j) for j in df.columns]

    n = 'seg_{}_{}_orig_id.csv.gz'.format(protein_db, taxon_id)
    n = os.path.join('seg', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df.to_csv(n, index=True, compression='gzip')

    df_g = mapper.uniprot_protein_2_gene_ncbi(df, 'median')
    n = 'seg_{}_{}_gene_id.csv.gz'.format(protein_db, taxon_id)
    n = os.path.join('seg', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df_g.to_csv(n, index=True, compression='gzip')


def retreive_signalp(taxon_id, protein_db):
    """
    Obtains the signalP (Signal Peptide, and multi-transmembrane)
    from gene-specific computations of one taxon,
    and saves them as a comma-separated file

    Requirement:
        Predictions run through signalP

    Input:
        taxon_id    integer, ncbi_taxon_id (e.g.: 9606 for Homo sapiens)

    Output:
        secondary/signalp/signalp_TAXONID.csv

    """

    ext = os.path.join(
        'virtualexchange',
        'batches',
        protein_db,
        'signalp',
        str(taxon_id),
        '*.signalp')
    pa = io.get_internal_path(ext)

    files = glob.glob(pa)

    agg_results = []

    for currf in files:
        fi = open(currf, 'r')
        row = fi.readlines()
        fi.close()

        header_lines = 2
        curr_results = row[header_lines:]

        for p in curr_results:

            s = re.sub('\s+', ' ', p).strip()
            e = s.split()

            m = '^.*\|(.*)\|.*$'
            ma = re.search(m, e[0])
            uniprot_id = ma.group(1)

            # For detailed explanation see Petersen et al., NMETH, 2011
            r = {
                'protein_uniprot': uniprot_id,
                'cmax': float(e[1]),
                'start_pos_of_mature_protein': int(e[2]),
                'cleaved': e[9] == 'Y',
                'has_four_transmembrane_residues': e[11] == 'SignalP-TM'}

            agg_results.append(r)

    if len(agg_results) == 0:
        print('Did not find entries for', taxon_id, 'in', protein_db)
        return

    df = pd.DataFrame(agg_results)
    df = df.set_index('protein_uniprot', verify_integrity=True)
    df.columns = ['SignalP_{}: {}'.format(protein_db, j) for j in df.columns]

    n = 'signalp_{}_{}_orig_id.csv.gz'.format(protein_db, taxon_id)
    n = os.path.join('signalp', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df.to_csv(n, index=True, compression='gzip')

    df_g = mapper.uniprot_protein_2_gene_ncbi(df, 'median')
    n = 'signalp_{}_{}_gene_id.csv.gz'.format(protein_db, taxon_id)
    n = os.path.join('signalp', n)
    n = io.get_output_path(n)
    io.ensure_presence_of_directory(n)
    df_g.to_csv(n, index=True, compression='gzip')


'''
############## SUPPORT FUNCTIONS #####################

_get_aminoacid_features
_appears_to_be_canonical_cds
_count_ACGT_and_length
_get_codon_bias
_get_codon_usage
_extract_gene_id_from_qualifier
_extract_uniprot_from_uniprot_fasta_header
_has_only_ACGT
_insert_name_in_key
_isolate_taxon_protein_fasta_individual
_isolate_taxon_protein_fasta_joined
_lookup_path_of_uniprot_databases

'''


def _get_aminoacid_features(seq):
    """
    Takes a protein amino acid sequence, and extracts properties of those
    aminoacids. It will ignore undefined aminoacids completely and remove
    them, essentially leaving slightly truncated protein

    Input:
        seq     str, one-letter code of amino acids

    Output:
        amino_acid_properties   dict, with diverse properties of amino acids,
                                such as counts, or derived biophysical
                                properties
    """

    # obtain string from sequene object
    seq = str(seq)

    # Remove amiguous amino acids
    if ('B' in seq) | ('X' in seq) | ('Z' in seq):
        has_undefined_amino_acid = True
        seq = re.sub('[BXZ]', '', seq)
    else:
        has_undefined_amino_acid = False

    length_of_protein = len(seq)

    x = ProteinAnalysis(seq)
    fr = x.get_amino_acids_percent()

    pp = dict()  # protein properties
    pp['basic'] = fr['R'] + fr['K'] + fr['H']
    pp['acidic'] = fr['D'] + fr['E']
    pp['charged'] = pp['acidic'] + pp['acidic']
    pp['polar'] = (
        fr['G'] + fr['S'] + fr['Y'] + fr['C'] +
        fr['N'] + fr['T'] + fr['Q'])
    pp['polar_uncharged'] = fr['S'] + fr['T'] + fr['N'] + fr['Q']
    pp['hydrophobic'] = (
        fr['A'] + fr['V'] + fr['I'] + fr['L'] +
        fr['M'] + fr['F'] + fr['Y'] + fr['W'])
    pp['aromatic'] = fr['F'] + fr['W'] + fr['Y']
    pp['helix_affine'] = (
        fr['V'] + fr['I'] + fr['Y'] +
        fr['F'] + fr['W'] + fr['L'])
    pp['turn_affine'] = fr['N'] + fr['P'] + fr['G'] + fr['S']
    pp['sheet_affine'] = fr['E'] + fr['M'] + fr['A'] + fr['L']
    amino_acid_properties = pp

    count_of_pyrrolysine = seq.count('O')
    count_of_selenocystein = seq.count('U')
    fr['O'] = count_of_pyrrolysine / length_of_protein
    fr['U'] = count_of_selenocystein / length_of_protein

    # some features can not be computed for pyrrolysine and selenocyteine
    # Manually restore or keep features (length, molecular weight)
    if (count_of_pyrrolysine) | (count_of_selenocystein):

        # Remove non-canonical amino acids
        seq = re.sub('[OU]', '', seq)
        x_non_extended_amino_acids = ProteinAnalysis(seq)

        # Masses from http://web.expasy.org/findmod/findmod_masses.html#AA
        average_mw_U = 150.0388
        average_mw_O = 237.3018

        mw = x_non_extended_amino_acids.molecular_weight()
        mw = (
            mw +
            (count_of_selenocystein * average_mw_U) +
            (count_of_pyrrolysine * average_mw_O))

        gravy_non_extended_amino_acids = x_non_extended_amino_acids.gravy()

    else:  # Only canonical amino acids

        mw = x.molecular_weight()
        gravy_non_extended_amino_acids = x.gravy()

    amino_acid_properties['molecular_weight'] = (
        mw)
    amino_acid_properties['gravy_ignoring_O_and_U'] = (
        gravy_non_extended_amino_acids)
    amino_acid_properties['isoelectric_point'] = (
        x.isoelectric_point())
    amino_acid_properties['has_undefined_amino_acid'] = (
        has_undefined_amino_acid)

    amino_acid_properties.update(
        fr)
    amino_acid_properties.update(
        {'amount_measured_amino_acids': length_of_protein})

    return amino_acid_properties


def _appears_to_be_canonical_cds(seq):
    """
    Checks if a given sequence starts with start codon,
    and ends with stop codon, and is consisting of an
    integer multiple of three, which suggests absence
    of frameshifts

    Input:
        seq     str, string of nucleotides

    Output:
        appears_canonical    bool;
    """

    appears_canonical = False

    if len(seq) % 3 == 0:

        if seq.startswith('ATG'):
            if _has_only_ACGT(seq):

                if seq.endswith('TAA'):
                    appears_canonical = True
                elif seq.endswith('TGA'):
                    appears_canonical = True
                elif seq.endswith('TAG'):
                    appears_canonical = True

    return appears_canonical


def _count_ACGT_and_length(seq):
    """
    - count fraction of Adenine, Cytosine, Guanine, and Thymidin
    – count fraction of Cytosin + Guanin (CG content)
    - count sum of Adenine, Cytosine, Guanine, and Thymidin
        (ignore undefined nucleotides)

    Input:
        seq     str, sequence

    Output:
        res     dictionary, containing
                    - A
                    - C
                    - G
                    - T
                    - CG
                    - SumACGT
    """

    A = 0
    C = 0
    G = 0
    T = 0

    for j in seq:

        if j == 'A':
            A += 1
        elif j == 'C':
            C += 1
        elif j == 'G':
            G += 1
        elif j == 'T':
            T += 1

    tot_ACGT = A + C + G + T

    res = dict()
    res['A'] = A / tot_ACGT
    res['C'] = C / tot_ACGT
    res['G'] = G / tot_ACGT
    res['T'] = T / tot_ACGT
    res['CG'] = res['C'] + res['G']
    res['SumACGT'] = tot_ACGT

    return res


def _extract_gene_id_from_qualifier(list_of_qualifiers):
    """
    Will parse a list of qualifiers and return ncbi gene ID.
    If Gene ID is ambiguous, or not found, nan will be returned

    Input:
        list_of_qualifiers  list of strings
                            (note: e.g in refseq genome format)

    Output:
        gene_id     ncbi (entrez) gene ID
        found_unambiguous   bool, reports whether a single
                            unambigous gene was found
    """

    gene_id = np.nan
    found_unambiguous = False
    found_genes = 0

    for entry in list_of_qualifiers:
        if entry.startswith('GeneID:'):
            if found_genes == 0:
                gene_id = int(entry[7:])
                found_genes += 1
                found_unambiguous = True
            if found_genes > 1:
                new_gene_id = int(entry[7:])
                if gene_id != new_gene_id:
                    found_unambiguous = False
                    gene_id = np.nan

    return gene_id, found_unambiguous


def _extract_uniprot_from_uniprot_fasta_header(description):
    """
    Extrects the uniprot ID from a uniprot header

    Input:
        description     str, uniprot fasta header

    Output:
        desc            dictionary with metainfo

    """
    m = '^.*\|(.*)\|.*PE=([0-9]*) SV=([0-9]*)$'
    ma = re.search(m, description)
    desc = dict()
    desc['protein_uniprot'] = (ma.group(1))

    return desc


def _get_codon_bias(seq):
    '''
    Obtains several metrics for codon usage bias. Needs
    code code from companion repository which can not be
    shared publically. More specifically, it uses unpublished
    code developped by Sophia Liu from Luis Amaral's lab

    Input:
        seq     str, coding sequence

    Output:
        codon_bias  dict, contains several metrics of codon bias

    '''

    o = CodonBiasCalculator.ENC_Calculator(seq)

    # Part of Science of Biology v0.1
    codon_bias = dict()
    codon_bias['codon_bias_wright1990'] = o.get_Nc_Wright()
    codon_bias['codon_bias_novembre2002'] = o.get_Nc_Novembre()
    codon_bias['codon_bias_sun2013'] = o.get_Nc_Sun()

    # Further metrics
    codon_bias['codon_bias_cdc'] = o.get_CDC()
    codon_bias['codon_bias_rcbs'] = o.get_RCBS()
    codon_bias['codon_bias_scuo'] = o.get_SCUO()

    return codon_bias


def _get_codon_usage(seq_cds):
    """
    Counts the usage of standard codons. Note that codons
    with non-standard nucleotides (e.g: N indicative for undefined)
    will be igonred.

    Input:
        seq_cds     str, coding sequence

    Output:
        dict_codons dict, contains usage of individual codons
                            (amount of codon / amount of all codons)

    """

    length_of_cds = len(seq_cds)
    if length_of_cds % 3 != 0:
        raise ValueError(
            'Length of seq_cds needs to be divisible by three.'
            'Please check seq_cds is indeed a complete coding sequence'
            'Note that codons consist of three bases.')

    codons = (seq_cds[n: n + 3] for n in range(0, length_of_cds, 3))

    # explicit checking for codons:
    # Thus ignore non-standard letters (e.g.: sequencing errors)

    dict_codons = {
        'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0,
        'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0,
        'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0,
        'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0,
        'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0,
        'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0,
        'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 'GGA': 0, 'GGC': 0,
        'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0,
        'TAA': 0, 'TAC': 0, 'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0,
        'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0,
        'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0}

    total_canonical_codons = 0
    for codon in codons:
        dict_codons[codon] += 1
        total_canonical_codons += 1

    for key, value in dict_codons.items():
        dict_codons[key] = value / total_canonical_codons

    return dict_codons


def _get_seg_metrics(seq_record):
    """
    Parses outputs of seg (program for low compexity extraction)

    Input:
        seq_record    Biopython Sequence Object

    Output:
        res     list; of extracted features
        labels  list; labels corresponding to res

    """

    header = seq_record.description
    m = '^.*\|(.*)\|.*$'
    ma = re.search(m, header)
    uniprot_id = ma.group(1)

    sequence = seq_record.seq
    lc = [j == 'x' for j in sequence]  # x indicates low complexity

    current_region = 0
    longest_region = 0
    number_regions = 0
    number_regions_larger_5 = 0
    number_regions_larger_10 = 0
    number_regions_larger_20 = 0
    number_regions_larger_40 = 0

    len_protein = len(lc)  # Keep original length prior modification
    lc.append(False)       # Initiate stop for algorithm

    is_in_region = False

    for p in lc:
        if p:
            current_region += 1
            is_in_region = True
        else:

            if is_in_region:

                number_regions += 1

                if current_region > 5:
                    number_regions_larger_5 += 1

                if current_region > 10:
                    number_regions_larger_10 += 1

                if current_region > 20:
                    number_regions_larger_20 += 1

                if current_region > 40:
                    number_regions_larger_40 += 1

                if current_region > longest_region:
                    longest_region = current_region

            current_region = 0
            is_in_region = False

    total_unstructured = sum(lc)
    fra_unstructured = total_unstructured / len_protein
    fra_longest = longest_region / len_protein

    res = [
        uniprot_id,
        total_unstructured, fra_unstructured, longest_region,
        fra_longest, number_regions, number_regions_larger_5,
        number_regions_larger_10, number_regions_larger_20,
        number_regions_larger_40]

    labels = [
        'protein_uniprot',
        'Total', 'Fraction_of_protein', 'Longest',
        'Fraction_of_protein_by_longest', 'Regions', 'Regions_larger5',
        'Regions_larger10', 'Regions_larger20', 'Regions_larger40']

    return res, labels


def _has_only_ACGT(seq):
    """
    Checks if a sequence only consists of known ACGT nucleotides

    Input:
        seq     str, string of nucleotides

    Output:
        has_only_ACGT   bool; True if sequence only consists of ACGT

    """
    search = re.compile(r'[^ACGT]').search

    has_only_ACGT = not bool(search(seq))

    return has_only_ACGT


def _insert_name_in_key(orig_dict, prefix):
    """
    Inserts a prefix and _ in front of the name of every key
    in a given dictionary

    Input:
        orig_dict   dict
        prefix      str, will be used as prefix for keys

    Output:
        renamed_dict    dict, where keys now start with prefix_

    """

    renamed_dict = dict(
        (prefix + '_' + key, value) for (key, value) in orig_dict.items())
    return renamed_dict


def _isolate_taxon_protein_fasta_individual(taxon_id, protein_db):
    """
    Extracts single protein sequences - belonging to the taxon defined by
    taxon_id - as FASTA for starting batch processing (e.g.: by RADAR)
    Wil ignore protein fragments.

    Creates output folder for programs that are run subsequentially
    from a virtual box:
        - radar         repetitiveness of proteins
        - signalp       presence of signal peptide (for secretion of proteins)

    Requirements:
        uniprot databases

    Input:
        taxon_id       integer, ncbi_taxon_id (e.g.: 9606 for Homo sapiens)
        protein_db     str, name of protein database (e.g.: uniprot)

    Output:
        virtualexchange/batches/protein/batch_in/TAXONID
    """

    # Initialize
    species_of_interest = mapper.ncbi_taxon_2_uniprot_taxon(taxon_id)
    species_filter = 'OS=' + species_of_interest

    fp = _lookup_path_of_uniprot_databases(protein_db)

    # For input of batch jobs (run on dedicated virtual box)
    ext = os.path.join(
        'virtualexchange', 'batches', protein_db, 'batch_in', str(taxon_id))
    out_path = io.get_internal_path(ext)
    io.ensure_presence_of_directory(os.path.join(out_path))

    # For collecting batch jobs (run on dedicated virtual box)
    ext = os.path.join(
        'virtualexchange', 'batches', protein_db, 'radar', str(taxon_id))
    coll_path = io.get_internal_path(ext)
    io.ensure_presence_of_directory(os.path.join(coll_path))

    # For collecting batch jobs (run on dedicated virtual box)
    ext = os.path.join(
        'virtualexchange', 'batches', protein_db, 'signalp', str(taxon_id))
    coll_path = io.get_internal_path(ext)
    io.ensure_presence_of_directory(os.path.join(coll_path))

    for seq_record in SeqIO.parse((fp), 'fasta'):
        if species_filter in seq_record.description:
            if _is_protein_complete(seq_record.description):
                uniprot_id = re.search(
                    '^.*\|(.*)\|.*$', seq_record.id).group(1)
                o = os.path.join(out_path, '{}.fasta'.format(uniprot_id))

                SeqIO.write(seq_record, o, 'fasta')


def _isolate_taxon_protein_fasta_joined(taxon_id, protein_db):
    """
    Extracts  protein sequences - belonging to the taxon defined by
    taxon_id - as FASTA (e.g.: for processing by SEG).
    Wil ignore protein fragments.

    Creates output folder for programs that are run subsequentially
    from a virtual box:
        - seg       for unstructured-ness of proteins


    Requirements:
        uniprot databases

    Input:
        taxon_id       integer, ncbi_taxon_id (e.g.: 9606 for Homo sapiens)
        protein_db     str, name of protein database (e.g.: uniprot)

    Output:
        virtualexchange/batches/protein/batch_in/TAXONID/all.fasta
    """

    # Initialize
    species_of_interest = mapper.ncbi_taxon_2_uniprot_taxon(taxon_id)
    species_filter = 'OS=' + species_of_interest

    fp = _lookup_path_of_uniprot_databases(protein_db)

    # For input of batch jobs (run on dedicated virtual box)
    ext = os.path.join(
        'virtualexchange',
        'batches',
        protein_db,
        'all_in',
        str(taxon_id),
        'all.fasta')
    out_path = io.get_internal_path(ext)
    io.ensure_presence_of_directory(os.path.join(out_path))

    # For collecting batch jobs (run on dedicated virtual box)
    ext = os.path.join(
        'virtualexchange', 'batches', protein_db, 'seg')
    coll_path = io.get_internal_path(ext)
    io.ensure_presence_of_directory(os.path.join(coll_path))

    output_handle = open(out_path, 'w')
    for seq_record in SeqIO.parse((fp), 'fasta'):
        if species_filter in seq_record.description:
            if _is_protein_complete(seq_record.description):
                SeqIO.write(seq_record, output_handle, 'fasta')

    output_handle.close()


def _is_protein_complete(description):
    """
    Returns wheter (Fragment) was not found in descripion
    """

    found = '(Fragment)' in description
    is_complete = not(found)
    return is_complete


def _load_gerstein_expression(p_expression_file):
    """
    Loads gerstein lab excel file with expression data into
    a pandas dataframe; removes extra columns that are not needed
    """

    df = pd.read_excel(p_expression_file)

    # unamibigous, and formally correct, nomenclature
    df = df.rename(columns={'Gene': 'LocusTag'})
    df['LocusTag'] = 'CELE_' + df['LocusTag'].copy()

    # Define columns that contain transcription (RNA) data
    # list of columns is manually curated, and excludes
    # epigenetic traces
    allowed_gene_expression = [
        'N2_EE_50-0',
        'N2_EE_50-30',
        'N2_EE_50-60',
        'N2_EE_50-90',
        'N2_EE_50-120',
        'N2_EE_50-150',
        'N2_EE_50-180',
        'N2_EE_50-210',
        'N2_EE_50-240',
        'N2_EE_50-300',
        'N2_EE_50-330',
        'N2_EE_50-360',
        'N2_EE_50-390',
        'N2_EE_50-420',
        'N2_EE_50-450',
        'N2_EE_50-480',
        'N2_EE_50-510',
        'N2_EE_50-540',
        'N2_EE_50-570',
        'N2_EE_50-600',
        'N2_EE_50-630',
        'N2_EE_50-660',
        'N2_EE_50-690',
        'N2_EE_50-720',
        'EmMalesHIM8',
        'EmMalesHIM8-2',
        'EmMalesHIM8_EmMalesHIM8-2',
        'N2_4cell_EE_RZ-56',
        'N2_E2-E8_sorted',
        'N2_EE_DSN-51',
        'EE',
        'N2_EE-2',
        'EE_N2_EE-2',
        'N2_EE_RZ-54',
        'LE',
        'N2_LE-1',
        'LE_N2_LE-1',
        'L1',
        'N2_L1-1',
        'L1_N2_L1-1',
        'LIN35',
        'N2_L2_RZ-53',
        'N2_L2_DSN-50',
        'L2',
        'N2_L2-4',
        'L2_N2_L2-4',
        'L3',
        'N2_L3-1',
        'L3_N2_L3-1',
        'DauerEntryDAF2',
        'DauerEntryDAF2-2',
        'DauerEntryDAF2-1-1',
        'DauerEntryDAF2-4-1',
        'DauerEntryDAF2_DauerEntryDAF2-2_DauerEntryDAF2-1-1_DauerEntryDAF2-4-1',
        'DauerDAF2',
        'DauerDAF2-2',
        'DauerDAF2-2-1',
        'DauerDAF2-5-1',
        'DauerDAF2_DauerDAF2-2_DauerDAF2-2-1_DauerDAF2-5-1',
        'DauerExitDAF2-2',
        'DauerExitDAF2-3-1',
        'DauerExitDAF2-6-1',
        'DauerExitDAF2-2_DauerExitDAF2-3-1_DauerExitDAF2-6-1',
        'L4',
        'L4b',
        'L4_L4b',
        'L4JK1107soma',
        'L4JK1107soma-2',
        'L4JK1107soma_L4JK1107soma-2',
        'L4MALE',
        'L4MALE5',
        'L4MALE_L4MALE5',
        'YA',
        'N2_Yad-1',
        'YA_N2_Yad-1',
        'N2_YA_RZ-1',
        'AdultSPE9',
        'N2_Ad_gonad-1-RZLI',
        'PharyngealMuscle',
        'DC-1-5',
        'DC-2-12',
        'OPDC-2-12',
        'EF-1-24',
        'PL-2-24',
        'Hsph',
        'HsphEcoliCntl',
        'SmacDb10',
        'SmacDb10EcoliCntl',
        'Harpo',
        'HarpoEcoliCntl',
        'DMM386-NSML_NSMR-nrn_L1',
        'DMM401_N2all_L1-DSN',
        'DMM402_N2all_L1-DSN',
        'DMM408_Amot_nrn_L2-DSN',
        'DMM414_Amot_nrn_L2-DSN',
        'DMM415_Amot_nrn_L2-DSN',
        'DMM408_414_415_Amot_nrn_L2-DSN',
        'DSN-Negative-Positive',
        'DMM387-NSML_NSMR-nrn_L1-V',
        'DMM383-all-nrn_L1-V',
        'DMM391-all-nrn_L1-V',
        'DMM383_391-all-nrn_L1-V',
        'DMM401-N2all_L1-V',
        'DMM402-N2all_L1-V',
        'DMM401_402-N2all_L1-V',
        'DMM389_NSM_L1',
        'DMM381_all-nrn_L1',
        'DMM383_all-nrn_L1',
        'DMM239_Z1Z4_Em',
        'DMM260_N2ref_EE',
        'AG1201'
    ]

    # Create DataFrame with transcript expression data
    columns_to_consider = ['LocusTag'] + allowed_gene_expression
    df = df.loc[:, columns_to_consider]

    df = df.drop_duplicates('LocusTag', keep=False)
    df = df.set_index('LocusTag')

    prefix = 'log10_RNA_expression__gerstein_expression_'
    df.columns = ['{}_{}'.format(prefix, j) for j in df.columns]

    threshold = 1
    f = df < threshold
    df[f] = np.nan
    df = df.apply(np.log10)

    return df


def _load_uniprot_mapper():
        """
        Loads uniprot mpper completely

        Output:
            df_uniprot_mapper   dataframe with mapper
        """

        # Checking for presence of required file (generated by
        # prepare.uniprot_id_mapper() , a tier 1 prepare function)
        p_uniprot_mapper = io.get_output_path('uniprot/uniprot_id_mapper.h5')
        if not os.path.exists(p_uniprot_mapper):
            raise EnvironmentError(
                '_load_uniprot_mapper() requires uniprot_id_mapper')

        # Load uniprot mapper
        df_uniprot_mapper = pd.read_hdf(p_uniprot_mapper, 'table')

        return df_uniprot_mapper


def _lookup_path_of_uniprot_databases(database_name):
    """
    Retreives path to uniprot databases, such as swiss-prot or trmbl

    Input:
        database_name   str; name of uniprot database
                            'swiss-prot'    manually curated
                            'trmbl'         computationally predicted
    Output:
        p_uniprot_database      str; path to selected uniprot database
    """

    dbs = {
        'swissprot': 'uniprot_sprot.fasta',
        'trembl': 'uniprot_trembl.fasta'
    }

    b = (
        'downloads/uniprot/pub/databases/uniprot/'
        'current_release/knowledgebase/complete')

    if database_name in dbs.keys():
        ext = os.path.join(b, dbs[database_name])
        p_uniprot_database = io.get_internal_path(ext)

    else:
        raise ValueError(
            'Did not find lookup for uniprot database'
            'Presently supported databases are swiss-prot and trmbl.')

    return p_uniprot_database


def _save_orig_and_ncbi_gene_mapped_tables(
        p_dir, filebase, df_orig, df_ncbi, prefix=None):

    """
    Saves two dataframes into the same folder, the first dataframe
    will have postfix _orig, whereas the other one will have postfix
    _ncbi_gene.

    Will export dataframes as gzipped csv, with index

    Input:
        p_dir       str, path to folder of export
        file_base   str, file base (without extensions or postfix)
        df_orig     dataframe, original dataframe
        df_ncbi     dataframe, mapped to ncbi_gene
        prefix      Insert prefix into data columns; by default None (ignore)

    Output:
        p_dir/      folder with exported dataframes

    """

    if prefix is not None:
        df_orig.columns = [prefix + y for y in df_orig.columns]
        df_ncbi.columns = [prefix + y for y in df_ncbi.columns]

    df_orig.to_csv(
        os.path.join(p_dir, '{}_orig.csv.gz'.format(filebase)),
        compression='gzip',
        index=True)

    df_ncbi.to_csv(
        os.path.join(p_dir, '{}_ncbi_gene.csv.gz'.format(filebase)),
        compression='gzip',
        index=True)
