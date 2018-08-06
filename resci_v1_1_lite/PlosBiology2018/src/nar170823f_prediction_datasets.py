from resci_formatting_language import *


from access_biology_data import annotation
from access_biology_data import properties as pr

import resci_non_standard_predictors as npr



def retreive_biophysics(ref_genes, taxon_id):

    d = dict()

    d['aminoacids'] = \
        swiss_complemented_with_trembl(
            pr.aminoacids_swissprot(taxon_id),
            pr.aminoacids_trembl(taxon_id))
    d['chromosomes'] = \
        replace_natsorted_first_column_by_integers(
            index_by_gene_ncbi(
                pr.chromosomes(taxon_id)))
    d['genbank_gene'] = \
        index_by_gene_ncbi(
            pr.genbank_gene(taxon_id))
    d['genbank_genomic_cds'] = \
        index_by_gene_ncbi(
            pr.genbank_genomic_cds(taxon_id))
    d['genbank_validated_rna'] = \
        index_by_gene_ncbi(
            pr.genbank_validated_rna(taxon_id))
    d['radar'] = \
        swiss_complemented_with_trembl(
            pr.radar_swissprot(taxon_id),
            pr.radar_trembl(taxon_id))
    d['seg'] = \
        swiss_complemented_with_trembl(
            pr.seg_swissprot(taxon_id),
            pr.seg_trembl(taxon_id))
    d['signalp'] = \
        swiss_complemented_with_trembl(
            pr.signalp_swissprot(taxon_id),
            pr.signalp_trembl(taxon_id))

    return d


def retreive_related_genes(ref_genes, taxon_id):

    d = dict()

    d['homologenes'] = \
        replace_nan_by_false(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.homologenes_only_for_genes_in_homologene(taxon_id)),
                ref_genes))
    d['interactors'] = \
        index_by_gene_ncbi(
        pr.interactions_rolland_2014(taxon_id))

    return d


def retreive_own_fame(ref_genes, taxon_id):

    d = dict()

    d['own_fame'] = \
        filter_by_reference_index(
            index_by_gene_ncbi(
                npr.own_rpo_log_fame(taxon_id)),
            ref_genes)

    return d


def retreive_related_homologenes(ref_genes, taxon_id):

    d = dict()

    d['homologenes'] = \
        replace_nan_by_false(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.homologenes_only_for_genes_in_homologene(taxon_id)),
                ref_genes))

    return d


def retreive_related_interactors(ref_genes, taxon_id):

    d = dict()

    d['interactors'] = \
        index_by_gene_ncbi(
        pr.interactions_rolland_2014(taxon_id))

    return d


def retreive_homologene_discoveries(ref_genes, taxon_id):

    d = dict()

    d['homologene_preceding_solo_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_discovery_in_homologenes(
                        taxon_id, 'first_solo_year')),
                ref_genes))

    d['homologene_preceding_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_discovery_in_homologenes(
                        taxon_id, 'first_year')),
                ref_genes))

    return d


def retreive_homologene_description(ref_genes, taxon_id):

    d = dict()

    d['homologene_description_preceding_solo_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_description_in_homologenes(
                        taxon_id, 'first_solo_year')),
                ref_genes))

    d['homologene_description_preceding_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_description_in_homologenes(
                        taxon_id, 'first_year')),
                ref_genes))

    return d


def retreive_homologene_literature(ref_genes, taxon_id):

    d = dict()

    d['homologene_papers'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_homologenes(
                        taxon_id, 'papers')),
                ref_genes)))

    d['homologene_attention'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_homologenes(
                        taxon_id, 'attention')),
                ref_genes)))

    return d


def retreive_literature_of_rolland_2014_interactors(ref_genes, taxon_id):

    d = dict()

    d['interactors_papers'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            filter_by_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_rolland_2014_interactors(
                        taxon_id, 'papers')),
                ref_genes)))

    d['interactors_attention'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            filter_by_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_rolland_2014_interactors(
                        taxon_id, 'attention')),
                ref_genes)))

    return d


def retreive_years_of_first_paper(ref_genes, taxon_id):

    d = dict()

    d['year_of_first_paper'] = \
        filter_by_reference_index(
        index_by_gene_ncbi(
            npr.year_of_description(
                taxon_id, 'first_year')),
        ref_genes)

    return d


def retreive_human_experiments(ref_genes, taxon_id):
    if taxon_id != 9606:
        raise EnvironmentError(
            'retreive_human_experiments only supports human.')

    d = dict()

    d['allelepool_lek_2016_aberration'] = \
        index_by_gene_ncbi(
            pr.allelepool_lek_2016_aberration(taxon_id))
    d['allelepool_lek_2016_anticipation'] = \
        index_by_gene_ncbi(
            pr.allelepool_lek_2016_anticipation(taxon_id))
    d['compartment_itzhak'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.compartment_itzhak_2016_determined(taxon_id)),
                ref_genes))
    d['cytoplasmic_itzhak'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.compartment_itzhak_2016_cytoplasmic(taxon_id)),
                ref_genes))
    d['essentiality_blomen'] = \
        index_by_gene_ncbi(
            pr.essentiality_blomen_2015(taxon_id))
    d['essentiality_hart_2015_crispr'] = \
        index_by_gene_ncbi(
            pr.essentiality_hart_2015_crispr(taxon_id))
    # d['essentiality_blomen_2015_sh'] = \    # biased set of short hairpins
    #     index_by_gene_ncbi(
    #         pr.essentiality_blomen_2015_sh(taxon_id))
    d['essentiality_wang_2015_score'] = \
        index_by_gene_ncbi(
            pr.essentiality_wang_2015_score(taxon_id))
    # d['gdi'] = \
    #     replace_nan_by_minus_one(
    #     enforce_reference_index(
    #         index_by_gene_ncbi(
    #             pr.gdi()),
    #         ref_genes))
    d['protein_abundance_itzhak_2015'] = \
        replace_nan_by_minus_one(   # indirectly contains addtional information
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.protein_abundance_itzhak_2015(taxon_id)),
                ref_genes))
    # d['rvis'] = \
    #     replace_nan_by_minus_one(
    #         enforce_reference_index(
    #             index_by_gene_ncbi(
    #                 pr.rvis()),
    #             ref_genes))
    d['thermal_stability_2017'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.thermal_stability_leuenberger_2017(taxon_id)),
                ref_genes))
    # d['transcript_abundance_gex_mantalek_170222'] = \  # human: from uhlen
    #     index_by_gene_ncbi(
    #         pr.transcript_abundance_gex_mantalek_170222(taxon_id))
    d['transcript_abundance_uhlen_2015_cells'] = \
        replace_nan_by_minus_one(
            index_by_gene_ncbi(
                pr.transcript_abundance_uhlen_2015_cells(taxon_id)))
    d['transcript_abundance_uhlen_2015_tissues'] = \
        replace_nan_by_minus_one(
            index_by_gene_ncbi(
                pr.transcript_abundance_uhlen_2015_tissues(taxon_id)))
    d['transcript_detection_uhlen_2015_cells'] = \
        index_by_gene_ncbi(
            pr.transcript_detection_uhlen_2015_cells(taxon_id))
    d['transcript_detection_uhlen_2015_tissues'] = \
        index_by_gene_ncbi(
            pr.transcript_detection_uhlen_2015_tissues(taxon_id))
    d['transcript_stability'] = \
        replace_nan_by_minus_one(    # indirectly contains information
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.transcript_halflife_tani_2012_assuming_48h_for_stable(taxon_id)),
                ref_genes))

    return d


def retreive_human_biased_experiments(ref_genes, taxon_id):
    if taxon_id != 9606:
        raise EnvironmentError(
            'retreive_human_experiments only supports human.')

    d = dict()
    d['compartment_hpa'] = \
        index_by_gene_ncbi(
            pr.compartment_thul_2017_observed(taxon_id))

    return d


def retreive_human_regulators(ref_genes, taxon_id):

    if taxon_id != 9606:
        raise EnvironmentError(
            'retreive_human_experiments only supports human.')

    d = dict()

    d['transcription_factors_genealacart_encode'] = \
        replace_nan_by_zero(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.transcription_factors_genealacart_encode(taxon_id)),
                ref_genes))
    d['transcription_factors_genealacart_promoters'] = \
        replace_nan_by_zero(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.transcription_factors_genealacart_promoters(taxon_id)),
                ref_genes))
    d['micro_rnas'] = \
        replace_nan_by_zero(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.mirdb(taxon_id)),
                ref_genes))

    d['amount_transcription_factors_encode'] = count_nonzero_in_row(
        d['transcription_factors_genealacart_encode'], 'amount_transcription_factors_encode')
    d['amount_transcription_factors_promoters'] = count_nonzero_in_row(
        d['transcription_factors_genealacart_promoters'], 'amount_transcription_factors_promoters')
    d['amount_micro_rnas'] = count_nonzero_in_row(
        d['micro_rnas'], 'amount_micro_rnas')

    d['micro_rnas'] = require_at_least_ten_non_zero_records(d['micro_rnas'])

    return d


def retreive_human_disease(ref_genes, taxon_id):

    if taxon_id != 9606:
        raise EnvironmentError(
            'retreive_human_experiments only supports human.')

    d = dict()

    d['unified_disease'] = \
        replace_nan_by_false(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.disease_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='unified_disease'),
                ref_genes)
    )
    d['orphanet_disease'] = \
        replace_nan_by_false(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.orphanet_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='orphanet_disease'),
                ref_genes)
    )
    d['human_phenotype'] = \
        replace_nan_by_false(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.human_phenotype_genealacart(
                        taxon_id=taxon_id, add_absenece=False)[[
                            'gene_ncbi',
                            'human_phenotype_genealacart: human_phenotype_name']],
                    col_name='human_phenotype_genealacart: human_phenotype_name'),
                ref_genes)
    )
    d['omim_disease'] = \
        replace_nan_by_false(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.omim_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='omim_disease'),
                ref_genes)
    )

    d['amount_unified_disease'] = count_nonzero_in_row(
        d['unified_disease'], 'amount_unified_diseases')
    d['amount_orphanet_diseases'] = count_nonzero_in_row(
        d['orphanet_disease'], 'amount_orphanet_diseases')
    d['amount_human_phenotypes'] = count_nonzero_in_row(
        d['human_phenotype'], 'amount_human_phenotypes')
    d['amount_omim_diseases'] = count_nonzero_in_row(
        d['omim_disease'], 'amount_omim_diseases')

    d['unified_disease'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['unified_disease']),
        'unfied_disease')
    d['orphanet_disease'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['orphanet_disease']),
        'orphanet')
    d['human_phenotype'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['human_phenotype']),
        'phenotype')
    d['omim_disease'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['omim_disease']),
        'omim')

    return d


def retreive_all_unified_disease(ref_genes, taxon_id):

    if taxon_id != 9606:
        raise EnvironmentError(
            'retreive_human_experiments only supports human.')

    d = dict()

    d['unified_disease'] = \
        replace_nan_by_false(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.disease_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='unified_disease'),
                ref_genes)
    )

    d['amount_unified_disease'] = count_nonzero_in_row(
        d['unified_disease'], 'amount_unified_diseases')

    return d


def retreive_all_omim_disease(ref_genes, taxon_id):

    if taxon_id != 9606:
        raise EnvironmentError(
            'retreive_human_experiments only supports human.')

    d = dict()

    d['omim_disease'] = \
        replace_nan_by_false(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.omim_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='omim_disease'),
                ref_genes)
    )

    d['amount_omim_diseases'] = count_nonzero_in_row(
        d['omim_disease'], 'amount_omim_diseases')

    return d


def retreive_domains(ref_genes, taxon_id):

    d = dict()

    d['interpro'] = \
        replace_nan_by_zero(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.interpro(
                        taxon_id=taxon_id)[['gene_ncbi', 'interpro_name']].drop_duplicates(),
                    col_name='interpro_name'),
                ref_genes
            ))

    d['amount_interpro'] = count_nonzero_in_row(
        d['interpro'], 'amount_interpro')

    d['interpro'] = require_at_least_twenty_non_zero_records(d['interpro'])

    return d


def retreive_single_sets(ref_genes, taxon_id):
    """
    retreives individual predictor sets;
    In contrast to joint sets there will be no reasonable substituiotn
    for absent values (e.g.: on rna stability: if no measurement
    since below expression value, this will not be substtiuted
    with a low expression value, as no expression dataset
    is incldued in same joint set of later application)

    """

    d = dict()

    d['aminoacids'] = \
        swiss_complemented_with_trembl(
            pr.aminoacids_swissprot(taxon_id),
            pr.aminoacids_trembl(taxon_id))
    d['chromosomes'] = \
        replace_natsorted_first_column_by_integers(
            index_by_gene_ncbi(
                pr.chromosomes(taxon_id)))
    d['genbank_gene'] = \
        index_by_gene_ncbi(
            pr.genbank_gene(taxon_id))
    d['genbank_genomic_cds'] = \
        index_by_gene_ncbi(
            pr.genbank_genomic_cds(taxon_id))
    d['genbank_validated_rna'] = \
        index_by_gene_ncbi(
            pr.genbank_validated_rna(taxon_id))
    d['radar'] = \
        swiss_complemented_with_trembl(
            pr.radar_swissprot(taxon_id),
            pr.radar_trembl(taxon_id))
    d['seg'] = \
        swiss_complemented_with_trembl(
            pr.seg_swissprot(taxon_id),
            pr.seg_trembl(taxon_id))
    d['signalp'] = \
        swiss_complemented_with_trembl(
            pr.signalp_swissprot(taxon_id),
            pr.signalp_trembl(taxon_id))

    # d['own_fame'] = \
    #     filter_by_reference_index(
    #         index_by_gene_ncbi(
    #             npr.own_rpo_log_fame(taxon_id)),
    #         ref_genes)

    d['homologenes'] = \
        replace_nan_by_false(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.homologenes_only_for_genes_in_homologene(taxon_id)),
                ref_genes))

    d['interactors'] = \
        index_by_gene_ncbi(
        pr.interactions_rolland_2014(taxon_id))

    d['homologene_preceding_solo_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_discovery_in_homologenes(
                        taxon_id, 'first_solo_year')),
                ref_genes))

    d['homologene_preceding_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_discovery_in_homologenes(
                        taxon_id, 'first_year')),
                ref_genes))

    d['homologene_description_preceding_solo_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_description_in_homologenes(
                        taxon_id, 'first_solo_year')),
                ref_genes))

    d['homologene_description_preceding_year'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.year_of_description_in_homologenes(
                        taxon_id, 'first_year')),
                ref_genes))

    d['homologene_papers'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_homologenes(
                        taxon_id, 'papers')),
                ref_genes)))

    d['homologene_attention'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            enforce_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_homologenes(
                        taxon_id, 'attention')),
                ref_genes)))

    d['interactors_papers'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            filter_by_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_rolland_2014_interactors(
                        taxon_id, 'papers')),
                ref_genes)))

    d['interactors_attention'] = \
        replace_nan_by_minus_two(
        replace_less_than_0p1_by_nan_and_log_transform(
            filter_by_reference_index(
                index_by_gene_ncbi(
                    npr.literature_of_rolland_2014_interactors(
                        taxon_id, 'attention')),
                ref_genes)))

    d['year_of_first_paper'] = \
        filter_by_reference_index(
        index_by_gene_ncbi(
            npr.year_of_description(
                taxon_id, 'first_year')),
        ref_genes)

    d['allelepool_lek_2016_aberration'] = \
        index_by_gene_ncbi(
            pr.allelepool_lek_2016_aberration(taxon_id))
    d['allelepool_lek_2016_anticipation'] = \
        index_by_gene_ncbi(
            pr.allelepool_lek_2016_anticipation(taxon_id))
    d['compartment_itzhak'] = \
        replace_nan_by_minus_one(
            filter_by_reference_index(
                index_by_gene_ncbi(
                    pr.compartment_itzhak_2016_determined(taxon_id)),
                ref_genes))
    d['cytoplasmic_itzhak'] = \
        replace_nan_by_minus_one(
            filter_by_reference_index(
                index_by_gene_ncbi(
                    pr.compartment_itzhak_2016_cytoplasmic(taxon_id)),
                ref_genes))
    d['essentiality_blomen'] = \
        index_by_gene_ncbi(
            pr.essentiality_blomen_2015(taxon_id))
    d['essentiality_hart_2015_crispr'] = \
        index_by_gene_ncbi(
            pr.essentiality_hart_2015_crispr(taxon_id))
    d['essentiality_blomen_2015_sh'] = \
        index_by_gene_ncbi(
            pr.essentiality_blomen_2015_sh(taxon_id))
    d['essentiality_wang_2015_score'] = \
        index_by_gene_ncbi(
            pr.essentiality_wang_2015_score(taxon_id))
    d['gdi'] = \
        replace_nan_by_minus_one(
        filter_by_reference_index(
            index_by_gene_ncbi(
                pr.gdi()),
            ref_genes))
    d['protein_abundance_itzhak_2015'] = \
        replace_nan_by_minus_one(   # indirectly contains addtional information
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.protein_abundance_itzhak_2015(taxon_id)),
                ref_genes))
    d['rvis'] = \
        replace_nan_by_minus_one(
            filter_by_reference_index(
                index_by_gene_ncbi(
                    pr.rvis()),
                ref_genes))
    d['thermal_stability_2017'] = \
        replace_nan_by_minus_one(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.thermal_stability_leuenberger_2017(taxon_id)),
                ref_genes))
    # d['transcript_abundance_gex_mantalek_170222'] = \  # human: from uhlen
    #     index_by_gene_ncbi(
    #         pr.transcript_abundance_gex_mantalek_170222(taxon_id))
    d['transcript_abundance_uhlen_2015_cells'] = \
        replace_nan_by_minus_one(
            index_by_gene_ncbi(
                pr.transcript_abundance_uhlen_2015_cells(taxon_id)))
    d['transcript_abundance_uhlen_2015_tissues'] = \
        replace_nan_by_minus_one(
            index_by_gene_ncbi(
                pr.transcript_abundance_uhlen_2015_tissues(taxon_id)))
    d['transcript_detection_uhlen_2015_cells'] = \
        index_by_gene_ncbi(
            pr.transcript_detection_uhlen_2015_cells(taxon_id))
    d['transcript_detection_uhlen_2015_tissues'] = \
        index_by_gene_ncbi(
            pr.transcript_detection_uhlen_2015_tissues(taxon_id))
    d['transcript_stability'] = \
        replace_nan_by_minus_one(    # indirectly contains information
            filter_by_reference_index(
                index_by_gene_ncbi(
                    pr.transcript_halflife_tani_2012_assuming_48h_for_stable(taxon_id)),
                ref_genes))

    d['compartment_hpa'] = \
        filter_by_reference_index(
        index_by_gene_ncbi(
            pr.compartment_thul_2017_observed(taxon_id)), ref_genes)

    d['transcription_factors_genealacart_encode'] = \
        replace_nan_by_zero(
            filter_for_genealacart(
                enforce_reference_index(
                    index_by_gene_ncbi(
                        pr.transcription_factors_genealacart_encode(taxon_id)),
                    ref_genes)))
    d['transcription_factors_genealacart_promoters'] = \
        replace_nan_by_zero(
        filter_for_genealacart(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.transcription_factors_genealacart_promoters(taxon_id)),
                ref_genes)))
    d['micro_rnas'] = \
        replace_nan_by_zero(
            enforce_reference_index(
                index_by_gene_ncbi(
                    pr.mirdb(taxon_id)),
                ref_genes))

    d['amount_transcription_factors_encode'] = count_nonzero_in_row(
        d['transcription_factors_genealacart_encode'], 'amount_transcription_factors_encode')
    d['amount_transcription_factors_promoters'] = count_nonzero_in_row(
        d['transcription_factors_genealacart_promoters'], 'amount_transcription_factors_promoters')
    d['amount_micro_rnas'] = count_nonzero_in_row(
        d['micro_rnas'], 'amount_micro_rnas')

    d['micro_rnas'] = require_at_least_ten_non_zero_records(d['micro_rnas'])

    d['unified_disease'] = \
        replace_nan_by_false(
        filter_for_genealacart(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.disease_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='unified_disease'),
                ref_genes))
    )
    d['orphanet_disease'] = \
        replace_nan_by_false(
        filter_for_genealacart(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.orphanet_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='orphanet_disease'),
                ref_genes))
    )
    d['human_phenotype'] = \
        replace_nan_by_false(
        filter_for_genealacart(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.human_phenotype_genealacart(
                        taxon_id=taxon_id, add_absenece=False)[[
                            'gene_ncbi',
                            'human_phenotype_genealacart: human_phenotype_name']],
                    col_name='human_phenotype_genealacart: human_phenotype_name'),
                ref_genes))
    )
    d['omim_disease'] = \
        replace_nan_by_false(
        filter_for_genealacart(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.omim_genealacart(
                        taxon_id=taxon_id, add_absenece=False),
                    col_name='omim_disease'),
                ref_genes)
        ))

    d['amount_unified_disease'] = count_nonzero_in_row(
        d['unified_disease'], 'amount_unified_diseases')
    d['amount_orphanet_diseases'] = count_nonzero_in_row(
        d['orphanet_disease'], 'amount_orphanet_diseases')
    d['amount_human_phenotypes'] = count_nonzero_in_row(
        d['human_phenotype'], 'amount_human_phenotypes')
    d['amount_omim_diseases'] = count_nonzero_in_row(
        d['omim_disease'], 'amount_omim_diseases')

    d['unified_disease'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['unified_disease']),
        'unfied_disease')
    d['orphanet_disease'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['orphanet_disease']),
        'orphanet')
    d['human_phenotype'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['human_phenotype']),
        'phenotype')
    d['omim_disease'] = prefix_column_names(
        require_at_least_ten_non_zero_records(
            d['omim_disease']),
        'omim')

    d['interpro'] = \
        replace_nan_by_zero(
            enforce_reference_index(
                pivot_gene_ncbi_to_predictor(
                    df=annotation.interpro(
                        taxon_id=taxon_id)[['gene_ncbi', 'interpro_name']].drop_duplicates(),
                    col_name='interpro_name'),
                ref_genes
            ))

    d['amount_interpro'] = count_nonzero_in_row(
        d['interpro'], 'amount_interpro')

    d['interpro'] = require_at_least_twenty_non_zero_records(d['interpro'])

    return d


# DEPRECATED


# The following block is tricky since imputation
# of absent values might be fair in one application,
# but wrong /misleading in other applications (where
# things could get shaped indirectly through solo-publication)


# IN ADDTION THE NAME OF THE KEY IS MISLEADING

# def retreive_own_description(ref_genes, taxon_id):

#     d = dict()

#     d['homologene_description_preceding_solo_year'] = \
#         filter_by_reference_index(
#         index_by_gene_ncbi(
#             npr.year_of_description(
#                 taxon_id, 'first_solo_year')),
#         ref_genes)

#     d['homologene_description_preceding_year'] = \
#         filter_by_reference_index(
#         index_by_gene_ncbi(
#             npr.year_of_description(
#                 taxon_id, 'first_year')),
#         ref_genes)

#     return d
