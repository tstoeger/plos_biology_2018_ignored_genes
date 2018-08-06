import pandas as pd


# from access_biology_data import meta
from access_science_shared import inout


def genome_rnai(taxon_id):
    """
    Will load genome_rnai for datasets mappable to a
    Pubmed ID

    Input:
        taxon_id  e.g.: 9606 for human

    Output:
        df     Genome RNAi table
    """

    if taxon_id == 9606:
        pheno_set = 'GenomeRNAi_v17_Homo_sapiens.txt'
    else:
        raise EnvironmentError(
            'This function presenlty only supports human.')

    # import genome RNAi dataset
    p = inout.get_path(
        'genome_rnai', pheno_set)

    agg = []

    is_valid_pubmed = False
    pubmed_id = 'not_set'
    with open(p, 'r') as rea:
        for line in rea:
            line = line.strip('\n')
            if line.startswith('#Pubmed ID='):
                pubmed_id = line[len('#Pubmed ID='):]
                is_valid_pubmed = len(pubmed_id) > 1
            elif line.startswith('#Screen Title'):
                pubmed_id = 'invalid'
            elif not line.startswith('#'):
                if not line.startswith('//'):
                    if is_valid_pubmed:
                        h = line.strip('\n').split('\t')
                        h = h + [pubmed_id]
                        agg.append(h)

    co = ['Stable ID', 'Entrez ID', 'Gene ID', 'Gene Symbol', 'Reagent ID',
          'Score', 'Phenotype', 'Conditions', 'Follow Up', 'Comment',
          'pubmed_id']
    df = pd.DataFrame(data=agg, columns=co)

    f = df['Phenotype'] == 'Inconclusive'
    df = df.loc[~f, :]

    df = df.rename(columns={
        'Entrez ID': 'gene_ncbi',
        'Gene ID': 'gene_id_ambiguous',
        'Gene Symbol': 'gene_symbol_ambiguous',
        'Stable ID': 'genome_rnai_id',
        'Phenotype': 'phenotype',
        'Reagent ID': 'reagent_id'
    })
    df = df.reset_index(drop=True)

    return df
