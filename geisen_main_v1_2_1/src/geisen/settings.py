import os
import re
import inspect
import pandas as pd


def get_data_path(requested='out'):
    mydir = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe())))
    data_dir = os.path.join(mydir, '..', '..', 'data')
    abs_path = os.path.abspath(data_dir)
    if requested == 'out':
        return os.path.join(abs_path, 'out')
    elif requested == 'internal':
        return os.path.join(abs_path, 'internal')


def retreive_general_download_settings():
    '''
    Loads settings file with individual download links,
    Excludes settings for genomes

    Output:
        df      DataFrame with setting file with 'name_of_dataset' as index
    '''

    file_path = _get_settings_file_path('download_links.txt')

    df = pd.read_table(
        file_path,
        sep='\t',
        header=0,
        index_col='name_of_dataset')

    return df


def retreive_genome_download_settings(taxon_of_interest=None):
    """
    Retrieves settings about downloading individual genomes. If
    taxon_of_interest is None, it will retreive information for
    all taxa in settings file, otherwise only for the specfied
    one.

    Input:
        taxon_of_interest   int; optional

    Output:
        genome_settings     dict, if taxon_of_interest is not defined:
                                for every taxon in settings file,
                                otherwise only for one

    The settigns file genome_links.txt has to be created
    manually since to ensure that the data is curated and
    not accidentally intrdocucing a mistake (for instance:
    that unmapped contigs are ot mistaken with microchromosomes
    of some species). Note that the settings file does
    not contain unmapped genes.

    The format of the settings file is:
    ## ... comment (e.g.: infos about curation)

    # xxxxx ... taxon ID (where x is taxonomy number)

    ftp://  ... folder of refseq genome built / annotation

    aaaa:bbbbbb  or
    aaaa:bbbbbb/ccccc where a: arbritrary lay name
                            b: refseq id
                            c: optional genbank id
    e.g.:

        # 9606

        ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10

        chromosome 1:NC_000001.11/CM000663.2
        chromosome 2:NC_000002.12/CM000664.2
        chromosome 3:NC_000003.12/CM000665.2

    """

    file_path = _get_settings_file_path('genome_links.txt')

    category = 'start'
    settings = dict()
    replicons = dict()
    info_of_taxon = dict()
    taxon = ''

    with open(file_path) as genome_links:
        for li in genome_links:

            li = li.strip()

            if li.startswith('##'):
                # Corresponds to comment in settings
                pass

            elif li.startswith('# '):

                if category == 'replicon':
                    info_of_taxon['replicons'] = replicons
                    settings[taxon] = info_of_taxon

                if category == 'link':
                    raise ValueError('settings file does not adhere.')

                taxon = int(li[2:])
                category = 'taxon'

                replicons = list()

            elif li.startswith('ftp://'):

                if category != 'taxon':
                    raise ValueError(
                        'Download link does not follow taxon definition.')

                info_of_taxon = dict()
                replicons = dict()

                info_of_taxon['link'] = li

                category = 'link'

            elif len(li) > 1:

                if category in ['link', 'replicon']:
                    pass
                else:
                    raise ValueError(
                        'Settings file in wrong format.'
                        'Replicons have to follow a link or a replicon.')

                # Extract lay name, and refseq name
                # (ignore genbank name same contig / replicon)
                lay_name = re.findall('(^.*):', li)[0]
                refseq_name = re.findall('.*:([A-Z_.0-9]*)', li)[0]

                replicons[lay_name] = refseq_name
                category = 'replicon'

    if category == 'replicon':
        info_of_taxon['replicons'] = replicons
        settings[taxon] = info_of_taxon
    else:
        raise ValueError('Settings file did not end with Replicon.')

    if taxon_of_interest is None:
        pass
    elif taxon_of_interest in settings.keys():
        settings = settings[taxon_of_interest]
    else:
        raise ValueError(
            'Did not find taxon {} in settings file {} '.format(
                taxon_of_interest, genome_links))

    return settings


def _get_settings_file_path(filebase):
    """
    Gets absolute path, which sits ./../../cfg relative
    to this file

    Input:
        filebase  str; name for file (including extension)
        file_path str; full path to the settings file

    """

    mydir = os.path.dirname(
        os.path.abspath(
            inspect.getfile(inspect.currentframe())))
    file_path = os.path.join(
        mydir, '..', '..', 'cfg', filebase)

    return file_path
