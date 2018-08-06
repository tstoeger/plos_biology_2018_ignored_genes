import os
import re

import numpy as np
import pandas as pd

from geisen import mapper

from geisen.prepare import _save_orig_and_ncbi_gene_mapped_tables


import geisen.inout as io


def blomen_2015():
    """
    Extracts fitness phenotypes from Blomen et al., and saves them
    together with their NCBI gene ID. Will only retreive the insertions
    of the crispr cassettes, and will do so for KBM7 and HAP1 cells.
    """

    p_out = io.get_output_path('papers/blomen_2015')
    io.ensure_presence_of_directory(p_out)

    def _tidy_blomen(file_path, cellline):
        s = cellline + '_full_dataset'
        d = pd.read_excel(file_path, sheetname=s, header=1)
        d['tot.insertions'] = d['tot.sense'] + d['tot.anti']
        d['selected'] = d['selected'] == 'YES'
        d = d.drop('GENE_SYMBOL', axis=1)
        d = d.set_index('ENSEMBL_ID')
        c = 'Blomen2015__' + cellline
        d.columns = [c + ': {}'.format(j) for j in d.columns]
        return d

    fp_KBM7 = io.get_geisen_manual_data_path(
        'out/papers/Blomen2015/aac7557_SM_Table_S1.xlsx')
    cellline = 'KBM7'
    k = _tidy_blomen(fp_KBM7, cellline)

    fp_HAP1 = io.get_geisen_manual_data_path(
        'out/papers/Blomen2015/aac7557_SM_Table_S2.xlsx')
    cellline = 'HAP1'
    h = _tidy_blomen(fp_HAP1, cellline)

    blomen2015 = pd.concat(
        [k, h],
        join='outer',
        verify_integrity=True,
        axis=1)
    blomen2015.index.name = 'gene_ensembl'  # science of biology nomenclature

    # Selecet features which describe insertions, rather
    # than ratios
    # (note: in science of biology v.0.1 this was part of the predict module)
    c = [
        'Blomen2015__KBM7: tot.sense',
        'Blomen2015__KBM7: tot.anti',
        'Blomen2015__KBM7: p.val',
        'Blomen2015__KBM7: q.val',
        'Blomen2015__HAP1: tot.anti',
        'Blomen2015__HAP1: p.val',
        'Blomen2015__HAP1: q.val']
    blomen2015 = blomen2015.loc[:, c]

    v = 'blomen_2015_fitness_orig'
    blomen2015.to_csv(
        os.path.join(p_out, '{}.csv.gz'.format(v)),
        compression='gzip',
        index=True)

    blomen2015_entrez = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
        blomen2015, taxon_id=9606)

    v = 'blomen_2015_fitness_ncbi_gene'
    blomen2015_entrez.to_csv(
        os.path.join(p_out, '{}.csv.gz'.format(v)),
        compression='gzip',
        index=True)


def hart_2015():
    """
    Extracts fitness phenotypes from Hart et al., and saves them
    together with their NCBI gene ID.

    will isoalte individual datasets as separte fiels
    """

    p_out = io.get_output_path('papers/hart_2015')
    io.ensure_presence_of_directory(p_out)

    p_in = io.get_geisen_manual_data_path(
        'out/papers/hart2015/mmc3_TSDeletedThoseWithExcelToDateConversion.xlsx'
    )

    hart2015 = pd.read_excel(p_in)
    hart2015 = hart2015.rename(columns={'Gene': 'symbol_ambiguous'})
    hart2015 = hart2015.set_index('symbol_ambiguous', verify_integrity=True)
    hart2015.columns = ['Hart2015: {}'.format(j) for j in hart2015.columns]

    hart2015_entrez = mapper.symbol_2_gene_ncbi(
        hart2015,
        taxon_id=9606,  # Homo sapiens
        how='median')

    out_settings = {  # cell-line : column name
        'hart2015_hct116_ordnum': 'Hart2015: BF_hct116',
        'hart2015_hela_ordnum': 'Hart2015: BF_hela',
        'hart2015_gbm_ordnum': 'Hart2015: BF_gbm',
        'hart2015_rpe1_ordnum': 'Hart2015: BF_rpe1',
        'hart2015_dld1_ordnum': 'Hart2015: BF_dld1',
        'hart2015_a375_ko_ordnum': 'Hart2015: BF_a375_GeCKo',
        'hart2015_hct116_shRNA_ordnum': 'Hart2015: BF_hct116_shRNA'
    }

    for cellline, dataset in out_settings.items():

        v = 'hart_2015_{}_ordnum_orig'.format(cellline)
        h = hart2015.loc[:, [dataset]]
        h.to_csv(
            os.path.join(p_out, '{}.csv.gz'.format(v)),
            compression='gzip',
            index=True)

        v = 'hart_2015_{}_ordnum_gene_ncbi'.format(cellline)
        h = hart2015_entrez.loc[:, [dataset]]
        h.to_csv(
            os.path.join(p_out, '{}.csv.gz'.format(v)),
            compression='gzip',
            index=True)


def itzhak_2016():
    """
    Protein localization, and abundance, as measured for HeLa cells
    by Itzhak et al. 2016
    """

    p_out = io.get_output_path('papers/itzhak_2016')
    io.ensure_presence_of_directory(p_out)

    p = io.get_geisen_manual_data_path(
        'out/papers/itzhak2016/'
        'elife-16950-supp1-v3-download-hela-spatial-proteome.csv')

    df = pd.read_csv(p)

    r = {
        'Lead Gene name':
            'symbol_ambiguous',
        'Lead Protein ID':
            'protein_uniprot',
        'Non-cytosolic pool1 ':
            'Non-cytosolic pool',
        'Global classifier2':
            'Global classifier',
        'Sub compart-ment Prediction':
            'Subcompartment Prediction',
        ' Contribution to cell protein mass [ppm]':
            'Contribution to cell protein mass [ppm]'
    }

    c = [
        'symbol_ambiguous',
        'Prediction Confidence',
        'Subcompartment Prediction',
        'Lead Protein name',
        'Mol. weight [kDa]',
        'Sequence length (AA)',
        'Total MS/MS Count',
        'Organellar profiles in how many maps?']

    df = df.rename(columns=r)
    df = df.drop(c, axis=1)

    df['Cytosolic Pool'] = df['Cytosolic Pool'].map(
        lambda x: int(x.rstrip('%')))
    df['Non-cytosolic pool'] = df['Non-cytosolic pool'].map(
        lambda x: int(x.rstrip('%')))

    df['Estimated Copy number per cell'] = df[
        'Estimated Copy number per cell'].str.replace(
            ',', '').astype(int)

    df['Compartment Prediction'] = df['Compartment Prediction'].fillna(
        value='not determined')
    df = df.set_index('protein_uniprot', verify_integrity=True)

    pr = 'Itzhak2016_'

    v = 'itzhak2016_compartment_nombool'
    f = df['Compartment Prediction'].isin(['not determined', 'No Prediction'])
    y = _nominal_ser_2_boolean_df(
        df.loc[~f, 'Compartment Prediction'])
    d = mapper.uniprot_protein_2_gene_ncbi(df=y, how='any')
    _save_orig_and_ncbi_gene_mapped_tables(p_out, v, y, d, pr)

    v = 'itzhak2016_global_classifier_nombool'
    y = _nominal_ser_2_boolean_df(
        df.loc[:, 'Compartment Prediction'])
    d = mapper.uniprot_protein_2_gene_ncbi(df=y, how='any')
    _save_orig_and_ncbi_gene_mapped_tables(p_out, v, y, d, pr)

    v = 'itzhak2016_localization_cytoplasm'
    y = df.loc[:, ['Cytosolic Pool']]  # adds up to 100 with non-cytoplasmic
    d = mapper.uniprot_protein_2_gene_ncbi(df=y, how='median')
    _save_orig_and_ncbi_gene_mapped_tables(p_out, v, y, d, pr)

    v = 'itzhak2016_localization_stats_ordnum'
    y = df.loc[:, ['Prediction Score']]
    d = mapper.uniprot_protein_2_gene_ncbi(df=y, how='median')
    _save_orig_and_ncbi_gene_mapped_tables(p_out, v, y, d, pr)

    v = 'itzhak2016_protein_abundance_ordnum'
    y = df.loc[:, [
        'Estimated Copy number per cell',
        'Copy number Abundance Percentile',
        'Median cellular con-centration [nM]',
        'Contribution to cell protein mass [ppm]']]
    d = mapper.uniprot_protein_2_gene_ncbi(df=y, how='median')
    _save_orig_and_ncbi_gene_mapped_tables(p_out, v, y, d, pr)


def lek_2016():
    """
    ExAc database, as published by Lek et al. 2016

    Output:
        lek2016_aberration_ordnum       enrichemnt of aberrations
        lek2016_aniticipation_ordnum    anticipated background rates
    """

    p_out = io.get_output_path('papers/lek_2016')
    io.ensure_presence_of_directory(p_out)

    # high level representation (at transcript level)
    p = io.get_geisen_manual_data_path(
        'out/papers/lek2016/nature19057-SI Table 13.xlsx')
    # data sheet with information on all genes
    df = pd.read_excel(p, sheetname='Gene Constraint')

    # reformatting
    df = df.rename(
        columns={'transcript': 'rna_ensembl'})  # controlled vocabulary
    df['rna_ensembl'] = df['rna_ensembl'].replace(
        '\..*$', '', regex=True)  # ignore versions of transcripts

    v = 'lek2016_aberration_ordnum'
    df_aberration = df[[
        'rna_ensembl',
        'syn_z',
        'mis_z',
        'lof_z',
        'pLI',
        'pRec',
        'pNull']].set_index('rna_ensembl')
    per_gene_aberration = mapper.rna_ensembl_2_gene_ncbi(
        df_aberration, how='median')
    _save_orig_and_ncbi_gene_mapped_tables(
        p_out,
        filebase=v,
        df_orig=df_aberration,
        df_ncbi=per_gene_aberration)

    v = 'lek2016_anticipation_ordnum'
    df_anticipation = df[[
        'rna_ensembl',
        'exp_syn',
        'exp_mis',
        'exp_lof']].set_index('rna_ensembl')
    per_gene_anticipation = mapper.rna_ensembl_2_gene_ncbi(
        df_anticipation, how='median')
    _save_orig_and_ncbi_gene_mapped_tables(
        p_out,
        filebase=v,
        df_orig=df_anticipation,
        df_ncbi=per_gene_anticipation)


def rolland_2014():
    """
    Processes supplemental data of Rolland et al. 2014
    (binary interaction; three methods) to extract:
    - interactions with same gene or other genes
        (stratified by support level)
    - binary interactin table (note: of genes with at least one interaction)
    - list of genes, which were tested

    Requirement:
        papers/rolland2014/mmc3.xlsx

    Output:
        rolland_considered_genes
        rolland_counts_of_interactions
        rolland_table_binary_interactions

    """

    p_in = io.get_geisen_manual_data_path('out/papers/rolland2014/mmc3.xlsx')
    p_out = io.get_output_path('papers/rolland_2014')
    io.ensure_presence_of_directory(p_out)

    sheets_of_interest = ['2B', '2G']
    rolland = pd.read_excel(p_in, sheetname=sheets_of_interest)

    bait_table = rolland['2B']

    considered_entrez = []
    count_of_invalid_baits = 0

    # Considered Genes
    for row in bait_table.itertuples():
        t = row.Tsdummyheader  # Had manually inserted header
        ma = re.search('entrez_gene_id=(.*)\|', t)
        if ma:
            matched = ma.group(1)
            if matched == 'NA':
                count_of_invalid_baits += 1
            else:
                attach = int(matched)
                considered_entrez.append(attach)

    considered_entrez = list(set(considered_entrez))
    print(
        'Rolland2014: Ignored {} baits that do not map to a gene.'.format(
            count_of_invalid_baits))

    v = 'rolland_considered_genes'
    df = pd.DataFrame(data=list(considered_entrez), columns=[v])
    df.to_csv(
        os.path.join(p_out, '{}.csv.gz'.format(v)),
        compression='gzip',
        index=False)

    # Create table where each gene of a non-self interaction occurrs
    # once as _ida, and once as _idb; note that this was ignored
    # by accident in science of biology v0.1
    interaction_table = rolland['2G']
    c = ['entrez_gene_ida', 'entrez_gene_idb', 'screens_found']
    f = interaction_table['screens_found'] > 0
    df = interaction_table.loc[f, c]
    df_i = df.iloc[:, [1, 0, 2]].copy()
    df_j = pd.concat([df, df_i], axis=0, ignore_index=True)
    df_j = df_j.drop_duplicates()  # safety to avoid counting self twice

    v = 'rolland_table_binary_interactions'
    df.to_csv(
        os.path.join(p_out, '{}.csv.gz'.format(v)),
        compression='gzip',
        index=False)

    # Count occurences (note: code for readability rather than speed)
    df = pd.DataFrame(
        index=considered_entrez,
        columns=[
            'self_interaction_any_evidence',
            'self_interaction_multiple_evidence',
            'trans_interaction_any_evidence',
            'trans_interaction_multiple_evidence',
        ])

    df = df.fillna(False)  # Python internally treates False and 0 as same
    df = df.sort_index()

    for row in df_j.itertuples():
        ix, id_a, id_b, support = row

        if id_a == id_b:
            df.loc[id_a, 'self_interaction_any_evidence'] = True
        else:
            df.loc[id_a, 'trans_interaction_any_evidence'] += 1

            if support > 1:
                if id_a == id_b:
                    df.loc[id_a, 'self_interaction_multiple_evidence'] = True
                else:
                    df.loc[id_a, 'trans_interaction_multiple_evidence'] += 1

    v = 'trans_interaction_multiple_evidence'  # appears to never occur
    if not(any(df[v])):
        df = df.drop(v, axis=1)

    df.columns = ['Rolland2014: {}'.format(j) for j in df.columns]

    v = 'rolland_counts_of_interactions'
    df.index.name = 'gene_ncbi'
    df.to_csv(
        os.path.join(p_out, '{}.csv.gz'.format(v)),
        compression='gzip',
        index=True)


def rosenfeld_2013():
    """
    Patent data on human genes. Note that companies usually patent
    an n-mer sequence, and its variants, thus they do not really
    patent individual genes, but sequences that have some similarity
    to genes.
    """

    p_in = io.get_geisen_manual_data_path(
        'out/papers/rosenfeld2013/13073_2013_415_MOESM1_ESM.XLS')
    df = pd.read_excel(p_in, skiprows=3)
    df = df.drop_duplicates()
    df = df.rename(columns={
        'Patent': 'patent',
        'Matching Gene': 'symbol_ncbi'})
    df_entrez = mapper.symbol_2_gene_ncbi(df, 9606, 'substitute')

    p_out = io.get_output_path('papers/rosenfeld_2013')
    io.ensure_presence_of_directory(p_out)

    v = 'rosenfeld_2013_patents'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=df,
        df_ncbi=df_entrez)


def thul_2017():
    """
    Protein subcellular localization from human protein
    atlas
    """

    p_in = io.get_geisen_manual_data_path(
        'out/papers/thul2017/aal3321_Thul_SM_table_S6.xlsx')
    p_out = io.get_output_path('papers/uhlen_2015')
    io.ensure_presence_of_directory(p_out)

    df = pd.read_excel(p_in)

    col = [
        'ENSG', 'Nucleus', 'Nucleoplasm', 'Nuclear bodies',
        'Nuclear speckles', 'Nuclear membrane', 'Nucleoli',
        'Nucleoli (Fibrillar center)', 'Cytosol', 'Cytoplasmic bodies',
        'Rods and Rings', 'Lipid droplets', 'Aggresome', 'Mitochondria',
        'Microtubules', 'Microtubule ends', 'Microtubule organizing center',
        'Centrosome', 'Mitotic spindle', 'Cytokinetic bridge', 'Midbody',
        'Midbody ring', 'Intermediate filaments', 'Actin filaments',
        'Focal Adhesions', 'Endoplasmic reticulum', 'Golgi apparatus',
        'Vesicles', 'Plasma membrane', 'Cell Junctions', 'Reliability']

    df = df.loc[:, col]
    df = df.rename(columns={'ENSG': 'gene_ensembl'})
    df = df.set_index('gene_ensembl', verify_integrity=True)

    df_entrez = \
        mapper.gene_ensembl_2_gene_ncbi_unambiguously(
            df,
            taxon_id=9606)

    v = 'thul_2017_subcellular_localization'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=df,
        df_ncbi=df_entrez)


def uhlen_2015():
    """
    - RNA transcirpt data form human protein atlas.
    - log transform fpkm
    - Expession treshold is 1 fpkm (0 in log transform), as
        in original paper

    """

    p_in = io.get_geisen_manual_data_path(
        'out/papers/uhlen2015/1260419_Excel_TablesS1-S18.xlsx')
    p_out = io.get_output_path('papers/uhlen_2015')
    io.ensure_presence_of_directory(p_out)

    def get_single_sheet(name_of_sheet):
        df = pd.read_excel(p_in, sheetname=[name_of_sheet])
        df = df[name_of_sheet]
        return df

    df_cell_lines = get_single_sheet('S11. FPKM Cell-lines')
    df_tissues = get_single_sheet('S18. Full FPKM dataset, tissues')

    def tidy_and_index(df):
        df = df.drop('gene_name', axis=1)
        df = df.set_index(['enstid'])   # They use wrong name, as identifiers
        df.index.name = 'gene_ensembl'  # are actually genes (each occurs once)
        threshold_used_by_Uhlen_2015 = 1    # Take author's detection threshold
        default_for_not_detected = np.nan   # and ignore values below

        f = df < threshold_used_by_Uhlen_2015
        df[f] = default_for_not_detected

        return df

    def log10_fun(x):
        y = x.applymap(lambda x: np.log10(x))
        return y

    df_cell_lines = tidy_and_index(df_cell_lines)
    df_tissues = tidy_and_index(df_tissues)
    df_cell_lines_log10 = log10_fun(df_cell_lines)
    df_tissues_log10 = log10_fun(df_tissues)

    df_cell_lines_log10.columns = ['uhlen_2015_cells_log10fpkm: {}'.format(
        j) for j in df_cell_lines_log10.columns]
    df_tissues_log10.columns = ['uhlen_2015_cells_log10fpkm: {}'.format(
        j) for j in df_tissues_log10.columns]

    # From Science of Biology v.0.1 / Predict module
    uhlen2015_tissues_levels = df_tissues_log10
    uhlen2015_cells_levels = df_cell_lines_log10

    uhlen2015_cells_levels.columns = [j.replace(
        '.MEAN', '') for j in uhlen2015_cells_levels.columns]

    def get_detected_fraction(df):
        d = 1 - df.isnull().sum(axis=1) / df.shape[1]
        return d

    detected_in_cells = get_detected_fraction(
        uhlen2015_cells_levels).to_frame(
        'uhlen_2015_fraction_detection_cells')
    detected_in_tissues = get_detected_fraction(
        uhlen2015_tissues_levels).to_frame(
        'uhlen_2015_fraction_detection_tissues')

    detected_in_cells_entrez = \
        mapper.gene_ensembl_2_gene_ncbi_unambiguously(
            detected_in_cells, taxon_id=9606)

    detected_in_tissues_entrez = \
        mapper.gene_ensembl_2_gene_ncbi_unambiguously(
            detected_in_tissues, taxon_id=9606)

    # correct identity of cell line, also see:
    # http://www.proteinatlas.org/learn/cellines
    uhlen2015_cells_levels = uhlen2015_cells_levels.rename(
        columns={
            'uhlen_2015_cells_log10fpkm: km3':
            'uhlen_2015_cells_log10fpkm: reh'
        })

    uhlen2015_cells_levels_entrez = \
        mapper.gene_ensembl_2_gene_ncbi_unambiguously(
            uhlen2015_cells_levels,
            taxon_id=9606)  # science of biology v.0.1 did log again

    uhlen2015_tissues_levels_entrez = \
        mapper.gene_ensembl_2_gene_ncbi_unambiguously(
            uhlen2015_tissues_levels,
            taxon_id=9606)

    v = 'uhlen_2015_detected_in_cells'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=detected_in_cells,
        df_ncbi=detected_in_cells_entrez)

    v = 'uhlen_2015_detected_in_tissuess'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=detected_in_tissues,
        df_ncbi=detected_in_tissues_entrez)

    v = 'uhlen_2015_cells_levels'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=uhlen2015_cells_levels,
        df_ncbi=uhlen2015_cells_levels_entrez)

    v = 'uhlen_2015_tissue_levels'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=uhlen2015_tissues_levels,
        df_ncbi=uhlen2015_tissues_levels_entrez)


def wang_2015():
    """
    Wang et al. 2015 (loss of function mutation monitoring fitness)

    """

    p_in = io.get_geisen_manual_data_path(
        'out/papers/wang2015/aac7041_SM_Table_S3.xlsx')
    p_out = io.get_output_path('papers/wang_2015')
    io.ensure_presence_of_directory(p_out)

    df = pd.read_excel(p_in)
    df = df.drop('sgRNAs included', axis=1)
    df = df.rename(columns={'Gene': 'symbol_ambiguous'})
    df = df.set_index('symbol_ambiguous', verify_integrity=True)

    # Remove K562 CS cells, as 39 of the 63 cell specific hits, are artifact
    # of genome location (see publication)
    excl = ['K562 CS', 'K562 adjusted p-value']
    df = df.drop(excl, axis=1)

    df.columns = ['Wang2015: {}'.format(j) for j in df.columns]

    c = [
        'Wang2015: KBM7 CS',
        'Wang2015: Jiyoye CS',
        'Wang2015: Raji CS']
    wang_cs = df.loc[:, c]

    wang_cs_entrez = mapper.symbol_2_gene_ncbi(
        wang_cs,
        taxon_id=9606,  # Homo sapiens
        how='median')

    c = [
        'Wang2015: KBM7 adjusted p-value',
        'Wang2015: Jiyoye adjusted p-value',
        'Wang2015: Raji adjusted p-value']
    wang_pvalue = df.loc[:, c]

    wang_pvalue_entrez = mapper.symbol_2_gene_ncbi(
        wang_pvalue,
        taxon_id=9606,  # Homo sapiens
        how='median')

    v = 'wang_2015_cs'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=wang_cs,
        df_ncbi=wang_cs_entrez)

    v = 'wang_2015_pvalue'
    _save_orig_and_ncbi_gene_mapped_tables(
        p_dir=p_out,
        filebase=v,
        df_orig=wang_pvalue,
        df_ncbi=wang_pvalue_entrez)


def _nominal_ser_2_boolean_df(ser_nominal):
    """
    Convert PANDAS dataframe with nominal values to a
    PANDAS dataframe with boolean cateogries. This is needed
    for some models; This implmentation also works on cases where
    scikit's own implementation breaks. Obains the base of the
    labels of columns fro the .name attirbute of ser_nominal

    Input:
        ser_nominal      Series containing nominal values

    Output:
        df_is_in_categorty     Dataframe with multiple columns,
                                    where each column corresponds
                                    to one category

    """

    is_present = pd.DataFrame(index=ser_nominal.index)
    unique_categories = list(set(ser_nominal))

    for x in range(len(unique_categories)):
        c = unique_categories[x]
        b = ser_nominal == c
        is_present[c] = b

    is_present.columns = ['{}: {}'.format(
        ser_nominal.name, j) for j in is_present.columns]
    df_is_in_categorty = is_present

    return df_is_in_categorty
